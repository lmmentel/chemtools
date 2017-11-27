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

from __future__ import print_function

from subprocess import Popen, call
import os

from .calculator import Calculator, InputTemplate, parse_objective


class Molpro(Calculator):
    '''
    Wrapper for the Molpro program.
    '''

    def __init__(self, name="Molpro", **kwargs):
        self.name = name
        self.inpext = '.inp'
        super(Molpro, self).__init__(**kwargs)

        self.molpropath = os.path.dirname(self.executable)

    def parse(self, fname, objective, regularexp=None):
        '''
        Parse a value from the output file ``fname`` based on the
        ``objective``.

        If the value of the ``objective`` is ``regexp`` then the
        ``regularexp`` will be used to parse the file.
        '''

        regexps = {
            'hf total energy'      : r'!RHF STATE \d+\.\d+ Energy\s+(\-?\d+\.\d+)',
            'mp2 total energy'     : r'!MP2 total energy\s+(\-?\d+\.\d+)',
            'ccsd total energy'    : r'!CCSD total energy\s+(\-?\d+\.\d+)',
            'ccsd(t) total energy' : r'!CCSD\(T\) total energy\s+(\-?\d+\.\d+)',
            'cisd total energy'    : r'!(?:RHF-R)?CISD\s+(?:total\s+)?energy\s+(\-?\d+\.\d+)',
            'fci total energy'     : r'!FCI STATE \d+\.\d+ Energy\s+(\-?\d+\.\d+)',
            'accomplished'         : r'\s*error',
        }

        if objective == 'regexp':
            toparse = regularexp
        else:
            toparse = regexps.get(objective, None)
            if toparse is None:
                raise ValueError("Specified objective: '{0:s}' not supported".format(objective))

        return parse_objective(fname, toparse)

    def run(self, inpfile):
        '''
        Run a single molpro job interactively - without submitting to the queue.
        '''

        if "-o" in self.runopts:
            outfile = self.runopts[self.runopts.index("-o") + 1]
        else:
            outfile = os.path.splitext(inpfile)[0] + ".out"

        command = [self.executable, inpfile] + self.runopts
        call(command)

        return outfile

    def run_multiple(self, inputs):
        '''
        Run a single molpro job interactively - without submitting to the queue.
        '''

        procs = []
        outputs = [os.path.splitext(inp)[0] + ".out" for inp in inputs]
        for inpfile, outfile in zip(inputs, outputs):
            opts = []
            opts.extend([self.executable, inpfile] + self.runopts)
            out = open(outfile, 'w')
            process = Popen(opts, stdout=out, stderr=out)
            out.close()
            procs.append(process)

        for p in procs: p.wait()

        return outputs

    def accomplished(self, fname):
        '''
        Return True if the job completed without errors
        '''

        # since the regexp is search for errors is None is found it is assumed
        # that the calcualtion is accomplished
        return self.parse(fname, 'accomplished') is None

    def __repr__(self):
        return "\n".join(["<Molpro(",
                          "\tname={},".format(self.name),
                          "\tmolpropath={},".format(self.molpropath),
                          "\texecutable={},".format(self.executable),
                          "\tscratch={},".format(self.scratch),
                          "\trunopts={},".format(str(self.runopts)),
                          ")>\n"])

    def write_input(self, fname=None, template=None, mol=None, basis=None, core=None):
        '''
        Write the molpro input to "fname" file based on the information from the
        keyword arguments.

        Args:
            mol : :py:class:`chemtools.molecule.Molecule`
                Molecule object instance
            basis : dict or :py:class:`BasisSet <chemtools.basisset.BasisSet>`
                An instance of :py:class:`BasisSet <chemtools.basisset.BasisSet>` class or a
                dictionary of :py:class:`BasisSet <chemtools.basisset.BasisSet>` objects with
                element symbols as keys
            core : list of ints
                Molpro core specification
            template : :py:class:`str`
                Template of the input file
            fname : :py:class:`str`
                Name of the input file to be used
        '''

        temp = InputTemplate(template)

        if isinstance(basis, dict):
            bs_str = "".join(x.to_molpro() for x in basis.values())
        else:
            bs_str = basis.to_molpro()

        if core is not None:
            core = "core,{0:s}\n".format(",".join([str(x) for x in core]))
        else:
            core = ''

        subs = {
            'geometry' : mol.molpro_rep(),
            'basis' : "basis={\n"+bs_str+"\n}\n",
            'core' : core,
        }

        with open(fname, 'w') as inp:
            inp.write(temp.substitute(subs))
