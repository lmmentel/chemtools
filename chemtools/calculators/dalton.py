# -*- coding: utf-8 -*-

#The MIT License (MIT)
#
#Copyright (c) 2014 Lukasz Mentel
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

import os
from collections import Counter
from subprocess import Popen, call

from .calculator import Calculator, InputTemplate, parse_objective

from ..basisset import get_l


class Dalton(Calculator):

    'Wrapper for running Dalton program'

    def __init__(self, name='Dalton', **kwargs):

        self.name = name
        super(Dalton, self).__init__(**kwargs)

        self.daltonpath = os.path.dirname(self.executable)

    def parse(self, fname, objective, regularexp=None):
        '''
        Parse a value from the output file ``fname`` based on the
        ``objective``.

        If the value of the ``objective`` is ``regexp`` then the
        ``regularexp`` will be used to parse the file.
        '''

        regexps = {
            'hf total energy': r'^@\s+Final HF energy:\s*(\-?\d+\.\d+)',
            'cisd total energy': r'\d+\s*\d+\s*(\-?\d+\.\d+).*converged',
            'accomplished': r'End of Wave Function Section',
        }

        if objective == 'regexp':
            if regularexp is None:
                raise ValueError("<regularexp> needs to be specified for objective='regexp'")
            toparse = regularexp
        else:
            toparse = regexps.get(objective, None)
            if toparse is None:
                raise ValueError("Specified objective: '{0:s}' not supported".format(objective))

        return parse_objective(fname, toparse)

    def run(self, fname):
        '''
        Run a single job

        Args:
            fname : dict
                A dictionary with keys ``mol`` and ``dal`` and their
                respective file name strings as values

        Returns:
            out : str
                Name of the dalton output file
        '''

        dalbase = os.path.splitext(fname['dal'])[0]
        molbase = os.path.splitext(fname['mol'])[0]

        command = [self.executable] + self.runopts + [dalbase, molbase]
        call(command)

        return dalbase + '_' + molbase + '.out'

    def run_multiple(self, fnames):
        '''
        Spawn two single jobs as paralell processes
        '''

        procs = []
        outputs = []
        for fname in fnames:

            dalbase = os.path.splitext(fname['dal'])[0]
            molbase = os.path.splitext(fname['mol'])[0]
            outputs.append(dalbase + '_' + molbase + '.out')

            command = [self.executable] + self.runopts + [dalbase, molbase]
            process = Popen(command)
            procs.append(process)

        for proc in procs:
            proc.wait()

        return outputs

    def write_input(self, fname, template, basis, mol, core):
        '''
        Write dalton input files: ``fname.dal`` and ``system.mol``

        Args:
            fname : str
                Name of the input file ``.dal``
            template : dict
                Dictionary with templates for the ``dal`` and ``mol``
                with those strings as keys and actual templates as
                values
            basis : dict
                An instance of
                :py:class:`BasisSet <chemtools.basisset.BasisSet>`
                class or a dictionary of
                :py:class:`BasisSet <chemtools.basisset.BasisSet>`
                objects with element symbols as keys
            mol : :py:class:`chemtools.molecule.Molecule`
                Molecule object with the system geometry
            core : str
                Core definition
        '''

        # Dalton uses atomic units for xyz coordinats by default

        daltemplate = template['dal']
        moltemplate = template['mol']

        # loop over different elements (not atoms)
        atomtypes = Counter([a.symbol for a in mol.atoms])
        out = ''
        for symbol, count in atomtypes.items():
            atoms = [a for a in mol.atoms if a.symbol == symbol]
            atombasis = basis[symbol]
            atombasis.sort()
            # get max angular momentum + 1 and construct block string
            maxb = max([get_l(s) for s in atombasis.functions.keys()]) + 1
            block = str(maxb) + ' 1' * maxb
            out += 'Atoms={0:d} Charge={1:.1f} Block={2:s}\n'.format(count,
                                                                    float(atoms[0].atomic_number),
                                                                    block)
            for i, atom in enumerate(atoms, start=1):
                out += '{0:4s} {1:15.8f} {2:15.8f} {3:15.8f}\n'.format(atom.symbol+str(i),
                                                              atom.xyz[0], atom.xyz[1], atom.xyz[2])

            out += atombasis.to_dalton()

        molsubs = {'basis' : out}
        moltemp = InputTemplate(moltemplate)
        dalsubs = {'core' : core}
        daltemp = InputTemplate(daltemplate)

        with open(fname['mol'], 'w') as fmol:
            fmol.write(moltemp.substitute(molsubs))

        with open(fname['dal'], 'w') as fdal:
            fdal.write(daltemp.substitute(dalsubs))

    def __repr__(self):
        return "\n".join(["<Dalton(",
                          "\tname={},".format(self.name),
                          "\tdaltonpath={},".format(self.daltonpath),
                          "\texecutable={},".format(self.executable),
                          "\tscratch={},".format(self.scratch),
                          "\trunopts={},".format(str(self.runopts)),
                          ")>\n"])
