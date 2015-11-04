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
# Copyright 2014-2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

import logging

from ..command import (Command, SolverCommandMixin, CommandError,
                       FilePrefixAppendAction)
from .. import fluxanalysis, fluxcoupling
from ..util import MaybeRelative

logger = logging.getLogger(__name__)


class FluxCouplingCommand(SolverCommandMixin, Command):
    """Find flux coupled reactions in the model.

    This identifies any reaction pairs where the flux of one reaction
    constrains the flux of another reaction. The reactions can be coupled in
    three distinct ways depending on the ratio between the reaction fluxes.
    The reactions can be fully coupled (the ratio is static and non-zero);
    partially coupled (the ratio is bounded and non-zero); or directionally
    coupled (the ratio is non-zero).
    """

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--reaction', default=None, action=FilePrefixAppendAction,
            type=str, help='Only include the given reactions')
        parser.add_argument(
            '--threshold', default=MaybeRelative('100%'), type=MaybeRelative,
            help='Flux threshold on objective reaction')
        parser.add_argument(
            '--mode', default='couplings', choices=['couplings', 'groups'],
            help='Find flux couplings or coupled groups')
        super(FluxCouplingCommand, cls).init_parser(parser)

    def run(self):
        solver = self._get_solver()

        max_reaction = self._model.get_biomass_reaction()
        if max_reaction is None:
            raise CommandError('The biomass reaction was not specified')

        if self._args.reaction is None:
            reactions = sorted(self._mm.reactions)
        else:
            reactions = self._args.reaction

        threshold = self._args.threshold
        if threshold.relative:
            fba_fluxes = dict(fluxanalysis.flux_balance(
                self._mm, max_reaction, tfba=False, solver=solver))
            threshold.reference = fba_fluxes[max_reaction]

        logger.info('Setting objective threshold to {}'.format(threshold))

        logger.info('Identifying blocked reactions...')

        self._blocked = {}
        for reaction_id, (minimum, maximum) in fluxanalysis.flux_variability(
                self._mm, reactions, {max_reaction: float(threshold)},
                tfba=False, solver=solver):
            blocked_neg = minimum >= 0.0
            blocked_pos = maximum <= 0.0
            self._blocked[reaction_id] = blocked_neg, blocked_pos

        self._fcp = fluxcoupling.FluxCouplingProblem(
            self._mm, {max_reaction: float(threshold)}, solver)

        self._coupled = {}
        self._groups = []

        if self._args.mode == 'couplings':
            self._find_couplings(reactions)
        elif self._args.mode == 'groups':
            self._find_groups(reactions)
        else:
            raise CommandError('Invalid mode: {}'.format(self._args.mode))

        logger.info('Coupled groups:')
        for i, group in enumerate(
                g for g in self._groups if g is not None and len(g) > 1):
            logger.info('Group {}: {}'.format(i, ', '.join(sorted(group))))

    def _find_couplings(self, reactions):
        logger.info('Identifying coupled reactions...')

        for reaction1 in reactions:
            # Skip reactions that are completely blocked
            if self._blocked[reaction1] == (True, True):
                logger.info('Skipping reaction {} because it is'
                            ' completely blocked'.format(reaction1))
                continue

            for reaction2 in reactions:
                for positive in (True, False):
                    blocked = self._blocked[reaction2][positive]
                    if blocked:
                        continue

                    lower, upper, text = self._check_reactions(
                        reaction1, reaction2, positive)
                    print('{}\t{}\t{}\t{}\t{}\t{}'.format(
                        reaction1, reaction2, 'pos' if positive else 'neg',
                        lower, upper, text))

    def _find_groups(self, reactions):
        for i, reaction1 in enumerate(reactions):
            if reaction1 in self._coupled:
                continue

            for reaction2 in reactions[i+1:]:
                if (reaction2 in self._coupled and
                        (self._coupled[reaction2] ==
                         self._coupled.get(reaction1))):
                    continue

                self._check_reactions(reaction1, reaction2)

        for reaction in reactions:
            if reaction in self._coupled:
                print('{}\t{}'.format(reaction, self._coupled[reaction]))

    def _check_reactions(self, reaction1, reaction2, positive):
        logger.debug('Solving for {}, {} {}'.format(
            reaction1, 'positive' if positive else 'negative', reaction2))

        lower, upper = self._fcp.solve(reaction1, reaction2, positive)

        logger.debug('Result: {}, {}'.format(lower, upper))

        coupling = fluxcoupling.classify_coupling((lower, upper))
        if coupling in (fluxcoupling.CouplingClass.DirectionalForward,
                        fluxcoupling.CouplingClass.DirectionalReverse):
            text = 'Directional, v1 / v2 in [{}, {}]'.format(lower, upper)
            if (coupling == fluxcoupling.CouplingClass.DirectionalReverse and
                    not self._mm.is_reversible(reaction1) and lower == 0.0):
                text = None
        elif coupling == fluxcoupling.CouplingClass.Full:
            text = 'Full: v1 / v2 = {}'.format(lower)
            self._couple_reactions(reaction1, reaction2)
        elif coupling == fluxcoupling.CouplingClass.Partial:
            text = 'Partial: v1 / v2 in [{}; {}]'.format(lower, upper)
            self._couple_reactions(reaction1, reaction2)
        else:
            text = None

        return lower, upper, text

    def _couple_reactions(self, reaction1, reaction2):
        logger.debug('Couple {} and {}'.format(reaction1, reaction2))

        if reaction1 in self._coupled and reaction2 in self._coupled:
            if self._coupled[reaction1] == self._coupled[reaction2]:
                return
            logger.debug('Merge groups {}, {}'.format(
                self._coupled[reaction1], self._coupled[reaction2]))
            group_index = len(self._groups)
            group = (self._groups[self._coupled[reaction1]] |
                     self._groups[self._coupled[reaction2]])
            logger.debug('New group is {}: {}'.format(
                group_index, sorted(group)))
            self._groups.append(group)
            self._groups[self._coupled[reaction1]] = None
            self._groups[self._coupled[reaction2]] = None
            for reaction in group:
                self._coupled[reaction] = group_index
        elif reaction1 in self._coupled:
            group_index = self._coupled[reaction1]
            group = self._groups[group_index]
            logger.debug('Put {} into existing group {}: {}'.format(
                reaction2, self._coupled[reaction1], group))
            self._coupled[reaction2] = group_index
            group.add(reaction2)
        elif reaction2 in self._coupled:
            group_index = self._coupled[reaction2]
            group = self._groups[group_index]
            logger.debug('Put {} into existing group {}: {}'.format(
                reaction1, self._coupled[reaction2], group))
            self._coupled[reaction1] = group_index
            group.add(reaction1)
        else:
            group = set([reaction1, reaction2])
            group_index = len(self._groups)
            logger.debug('Creating new group {}'.format(group_index))
            self._coupled[reaction1] = group_index
            self._coupled[reaction2] = group_index
            self._groups.append(group)
