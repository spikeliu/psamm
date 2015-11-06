# This file is part of PSAMM.
#
# PSAMM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PSAMM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PSAMM.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

"""Flux coupling analysis

Described in [Burgard04]_.
"""

import math
import logging

from six import iteritems

from psamm.lpsolver import lp


logger = logging.getLogger(__name__)


class CouplingClass(object):
    """Enumeration of coupling types."""

    Inconsistent = object()
    """Reaction is flux inconsistent (v1 is always zero)."""

    Uncoupled = object()
    """Uncoupled reactions."""

    DirectionalForward = object()  # v1 -> v2
    """Directionally coupled from reaction 1 to 2."""

    DirectionalReverse = object()  # v2 -> v1
    """Directionally coupled from reaction 2 to 1."""

    Partial = object()  # v1 <-> v2
    """Partially coupled reactions."""

    Full = object()  # v1 <=> v2
    """Fully coupled reactions."""


class FluxCouplingProblem(object):
    """A specific flux coupling analysis applied to a metabolic model.

    Args:
        model: MetabolicModel to apply the flux coupling problem to
        bounded: Dictionary of reactions with minimum flux values
        solver: LP solver instance to use.
    """

    def __init__(self, model, bounded, solver):
        self._prob = solver.create_problem()

        # Define t variable
        self._prob.define('t', lower=0)
        t = self._prob.var('t')

        # Define flux variables
        for reaction_id in model.reactions:
            self._prob.define(('vbow', reaction_id))

            lower, upper = model.limits[reaction_id]
            if reaction_id in bounded:
                lower = bounded[reaction_id]
            flux_bow = self._prob.var(('vbow', reaction_id))
            self._prob.add_linear_constraints(
                flux_bow >= t * lower, flux_bow <= t * upper)

        # Define mass balance constraints
        massbalance_lhs = {compound: 0 for compound in model.compounds}
        for (compound, reaction_id), value in iteritems(model.matrix):
            flux_bow = self._prob.var(('vbow', reaction_id))
            massbalance_lhs[compound] += flux_bow * value
        for compound, lhs in iteritems(massbalance_lhs):
            self._prob.add_linear_constraints(lhs == 0)

        self._reaction_constr = None

    def solve(self, reaction_1, reaction_2, positive):
        """Return the flux coupling between two reactions.

        The flux coupling is returned as a tuple indicating the minimum and
        maximum value of the v1/v2 reaction flux ratio. A value of None as
        either the minimum or maximum indicates that optimization is
        infeasible in that direction because the v2 might be blocked. A value
        of +/- infinity (float) can be returned for either value.

        The v2 is constrained to either positive values or negative values
        depending on the ``positive`` parameter.
        """
        # Update objective for reaction_1
        self._prob.set_linear_objective(self._prob.var(('vbow', reaction_1)))

        sign = 1 if positive else -1

        # Update constraint for reaction_2
        if self._reaction_constr is not None:
            self._reaction_constr.delete()

        reaction_2_vbow = self._prob.var(('vbow', reaction_2))
        self._reaction_constr, = self._prob.add_linear_constraints(
            reaction_2_vbow == sign)

        results = []
        for i, sense in enumerate(
                (lp.ObjectiveSense.Minimize, lp.ObjectiveSense.Maximize)):
            result = self._prob.solve(sense)
            value = None
            if result:
                value = sign * result.get_value(('vbow', reaction_1))
            elif result.unbounded:
                if sense == lp.ObjectiveSense.Minimize:
                    value = sign * -float('inf')
                else:
                    value = sign * float('inf')

            results.append(value)

        if not positive:
            return tuple(reversed(results))

        return tuple(results)


def classify_coupling(coupling):
    """Return a constant indicating the type of coupling.

    Depending on the type of coupling, one of the constants from
    :class:`.CouplingClass` is returned.

    Args:
        coupling: Tuple of minimum and maximum flux ratio
    """
    lower, upper = coupling

    if math.isinf(lower) and math.isinf(upper):
        return CouplingClass.Uncoupled
    elif math.isinf(lower) or math.isinf(upper):
        return CouplingClass.DirectionalReverse
    elif lower == 0.0 and upper == 0.0:
        return CouplingClass.Inconsistent
    elif lower <= 0.0 and upper >= 0.0:
        return CouplingClass.DirectionalForward
    elif abs(lower - upper) < 1e-6:
        return CouplingClass.Full
    else:
        return CouplingClass.Partial


def reversed_coupling(coupling):
    """Determine the reverse coupling values for a reaction pair.

    Args:
        coupling: Tuple of minimum and maximum flux ratio
    """
    def inverse(v):
        if v is None:
            return 0.0
        elif v == 0.0:
            return float('inf')
        return 1.0 / v

    lower, upper = coupling
    return inverse(upper), inverse(lower)
