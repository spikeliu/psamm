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

from __future__ import print_function

import re
import json
import sqlite3
import logging

from ..command import Command, FilePrefixAppendAction
from ..reaction import Compound, Reaction

logger = logging.getLogger(__name__)


def parse_compound(s):
    """Parse a compound specification with optional compartment"""
    m = re.match(r'([^\[]+)\[(\w+)\]', s)
    if m is not None:
        return Compound(m.group(1), compartment=m.group(2))
    return Compound(s)


def filter_search_term(s):
    return re.sub(r'[^a-z0-9]+', '', s.lower())


class SearchCommand(Command):
    """Search for reactions and compounds in the model."""

    @classmethod
    def init_parser(cls, parser):
        """Initialize argument parser"""
        subparsers = parser.add_subparsers(title='Search domain')

        # Compound subcommand
        parser_compound = subparsers.add_parser(
            'compound', help='Search in compounds')
        parser_compound.set_defaults(which='compound')
        parser_compound.add_argument(
            '--id', '-i', dest='id', metavar='id', type=str,
            help='Compound ID')
        parser_compound.add_argument(
            '--name', '-n', dest='name', metavar='name', type=str,
            help='Name of compound')

        # Reaction subcommand
        parser_reaction = subparsers.add_parser(
            'reaction', help='Search in reactions')
        parser_reaction.set_defaults(which='reaction')
        parser_reaction.add_argument(
            '--id', '-i', dest='id', metavar='id', type=str,
            help='Reaction ID')
        parser_reaction.add_argument(
            '--compound', '-c', dest='compound', metavar='compound',
            action=FilePrefixAppendAction, type=str, default=[],
            help='Compound ID (optionally including compartment)')

    def run(self):
        """Run search command."""

        self._conn = sqlite3.connect(
            ':memory', detect_types=sqlite3.PARSE_DECLTYPES)
        self._conn.execute('PRAGMA foreign_keys = ON;')
        c = self._conn.cursor()

        c.execute('''
CREATE TABLE compound (
  id TEXT PRIMARY KEY,
  name TEXT,
  norm_name TEXT,
  formula TEXT,
  properties TEXT,
  filemark TEXT);''')
        c.execute('''
CREATE TABLE reaction (
  id TEXT PRIMARY KEY,
  name TEXT,
  norm_name TEXT,
  properties TEXT,
  filemark TEXT);''')
        c.execute('''
CREATE TABLE compound_appearance (
  compound TEXT REFERENCES compound (id),
  reaction TEXT REFERENCES reaction (id),
  compartment TEXT,
  PRIMARY KEY (compound, reaction));''')

        def default(o):
            if isinstance(o, Reaction):
                return str(o)
            raise TypeError

        for compound in self.model.parse_compounds():
            name = compound.properties.get('name')

            norm_name = None
            if name is not None:
                norm_name = filter_search_term(name)

            filemark = None
            if compound.filemark is not None:
                filemark = str(compound.filemark)

            c.execute('''
INSERT INTO compound (id, name, norm_name, formula, properties, filemark)
  VALUES (?,?,?,?,?,?);''',
                      (compound.id, name, norm_name,
                       compound.properties.get('formula'),
                       json.dumps(compound.properties),
                       filemark))

        for reaction in self.model.parse_reactions():
            name = reaction.properties.get('name')

            norm_name = None
            if name is not None:
                norm_name = filter_search_term(name)

            filemark = None
            if reaction.filemark is not None:
                filemark = str(reaction.filemark)

            c.execute('''
INSERT INTO reaction (id, name, norm_name, properties, filemark)
  VALUES (?,?,?,?,?);''',
                      (reaction.id, name, norm_name,
                       json.dumps(reaction.properties, default=default),
                       filemark))

            for compound , _ in reaction.equation.compounds:
                c.execute('''
INSERT OR IGNORE INTO compound_appearance
    (compound, reaction, compartment)
  VALUES (?,?,?);''',
                          (compound.name, reaction.id, compound.compartment))

        c.execute('CREATE INDEX compound_norm_name ON compound (norm_name);')
        c.execute('CREATE INDEX reaction_norm_name ON reaction (norm_name);')
        c.execute('CREATE INDEX compound_formula ON compound (formula);')

        self._conn.commit()

        which_command = self.args.which
        if which_command == 'compound':
            self._search_compound()
        elif which_command == 'reaction':
            self._search_reaction()

    def _search_compound(self):
        where_clause = []
        parameters = []

        if self.args.id is not None:
            where_clause.append('id = ?')
            parameters.append(self.args.id)

        if self.args.name is not None:
            where_clause.append('norm_name LIKE ?')
            parameters.append(filter_search_term(self.args.name) + '%')

        if len(where_clause) > 0:
            where_clause = 'WHERE ' + (' AND '.join(where_clause))
        else:
            where_clause = ''

        query = '''
SELECT id, properties, filemark FROM compound
  {} ORDER BY id;'''.format(where_clause)
        c = self._conn.cursor()
        c.execute(query, parameters)

        # Show results
        for row in c:
            compound_id, properties, filemark = row[:3]
            properties = json.loads(properties)
            prop_names = set(properties) - {'id'}
            print('id: {}'.format(compound_id))
            for prop in sorted(prop_names):
                print('{}: {}'.format(prop, properties[prop]))
            if filemark is not None:
                print('Defined in {}'.format(filemark))
            print()

    def _search_reaction(self):
        where_clause = []
        where_params = []
        having_clause = []
        having_params = []

        if self.args.id is not None:
            where_clause.append('id = ?')
            where_params.append(self.args.id)

        # Match all compounds
        match_compounds_clause = []
        if len(self.args.compound) > 0:
            having_clause.append('compound_count = ?')
            having_params.append(len(self.args.compound))

        for compound_spec in self.args.compound:
            compound_clause = []
            compound = parse_compound(compound_spec)
            compound_clause.append('compound_appearance.compound = ?')
            where_params.append(compound.name)
            if compound.compartment is not None:
                compound_clause.append('compound_appearance.compartment = ?')
                where_params.append(compound.compartment)
            match_compounds_clause.append(' AND '.join(compound_clause))

        if len(match_compounds_clause) > 0:
            where_clause.append(' OR '.join(
                '(' + c + ')' for c in match_compounds_clause))

        where_clause.append('compound_appearance.reaction = reaction.id')
        if len(where_clause) > 0:
            where_clause = 'WHERE ' + (' AND '.join(
                '(' + c + ')' for c in where_clause))
        else:
            where_clause = ''

        if len(having_clause) > 0:
            having_clause = 'HAVING ' + (' AND '.join(
                '(' + c + ')' for c in having_clause))
        else:
            having_clause = ''

        parameters = where_params + having_params
        query = '''
SELECT id, properties, filemark,
    COUNT(compound_appearance.compound) AS compound_count
  FROM reaction, compound_appearance
  {} GROUP BY id {} ORDER BY id;'''.format(where_clause, having_clause)
        c = self._conn.cursor()
        c.execute(query, parameters)

        # Show results
        for row in c:
            reaction_id, properties, filemark = row[:3]
            properties = json.loads(properties)
            prop_names = set(properties) - {'id'}

            print('id: {}'.format(reaction_id))
            #print('equation: {}'.format(
            #    reaction.equation))
            #translated_equation = reaction.equation.translated_compounds(
            #    lambda x: compound_name.get(x, x))
            #if reaction.equation != translated_equation:
            #    print('equation (compound names): {}'.format(
            #        translated_equation))
            for prop in sorted(prop_names):
                print('{}: {}'.format(prop, properties[prop]))
            if filemark is not None:
                print('Defined in {}'.format(filemark))
            print()
