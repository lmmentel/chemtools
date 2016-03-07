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
from .code import Code
from subprocess import Popen, PIPE
import os
import re
import sys
from string import Template

class Psi4Template(Template):

    delimiter = '%'
    idpattern = r'[a-z][_a-z0-9]*'

    def get_keys(self):
        '''Parse the string and return a dict with possible keys to substitute.
        For most use case only the `named` fields are interesting'''

        keys = {}
        match = self.pattern.findall(self.template)
        for k, v in self.pattern.groupindex.items():
            keys[k] = [x[v-1] for x in match if x[v-1] != '']
        return keys

class Psi4(Code):
    '''
    Wrapper for the Psi4 program.
    '''

    def __init__(self, name="Psi4", **kwargs):
        self.name = name
        super(Psi4, self).__init__(**kwargs)

        self.psi4path = os.path.dirname(self.executable)

    def write_input(self, fname=None, template=None, mol=None, bs=None, core=""):
        '''
        Write the psi4 input to "fname" file based on the information from the
        keyword arguments.
        '''

        mi = Psi4Input(fname=fname, template=template)
        mi.write_input(mol=mol, bs=bs, core=core)

    def run(self, inpfile):
        '''
        Run a single Psi4 job interactively - without submitting to the queue.
        '''

        outfile = os.path.splitext(inpfile)[0] + ".out"
        errfile = os.path.splitext(inpfile)[0] + ".err"
        opts = []
        opts.extend([self.executable, inpfile] + self.runopts)

        process = Popen(opts, stdout=PIPE, stderr=PIPE)
        out, err = process.communicate()
        ferr = open(errfile, 'w')
        ferr.write(out)
        ferr.write("{0:s}\n{1:^80s}\n{0:s}\n".format("="*80, "Error messages:"))
        ferr.write(err)
        ferr.close()

        return outfile

    def run_multiple(self, inputs):
        '''
        Run a single Psi4 job interactively - without submitting to the queue.
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

    def parse(self, output, method, objective, regexp=None):
        '''
        Parse Psi4 output file to get the objective.
        '''

        parser = Psi4OutputParser(output)

        if objective == "total energy":
            if method == "hf":
                return parser.get_hf_total_energy()
            elif method == "cisd":
                return parser.get_cisd_total_energy()
        elif objective == "correlation energy":
                return parser.get_cisd_total_energy() - parser.get_hf_total_energy()
        elif objective == "core energy":
            if method == "cisd":
                return parser.get_cisd_total_energy()
        elif objective == "regexp":
            return parser.get_variable(regexp)
        else:
            raise ValueError("unknown objective in prase {0:s}".format(objective))

    def accomplished(self, outfile=None):
        '''
        Return True if Psi4 job finished without errors.
        '''

        if outfile is not None:
            parser = Psi4OutputParser(outfile)
        else:
            raise ValueError("oufile needs to be specified")
        return parser.accomplished()

    def __repr__(self):
        return "\n".join(["<Psi4(",
                        "\tname={},".format(self.name),
                        "\tpsi4path={},".format(self.psi4path),
                        "\texecutable={},".format(self.executable),
                        "\tscratch={},".format(self.scratch),
                        "\trunopts={},".format(str(self.runopts)),
                        ")>\n"])

class Psi4Input(object):
    '''
    Reading, parsing and writing of Psi4 input file.
    '''

    def __init__(self, fname=None, template=None):

        self.fname = fname
        self.template = template

    def write_input(self, mol=None, bs=None, core=None):
        '''
        Write the Psi4 input to "fname" file based on the information from the
        keyword arguments.

        Args:
          mol : chemtools.molecule.Molecule
            Molecule object instance
          bs : chemtools.basisset.BasisSet
            BasisSet class instance or list of those instances
          core : list of ints
            Psi4 core specification
        '''

        temp = Psi4Template(self.template)

        if isinstance(bs, list):
            bs_str = "".join(x.to_gaussian() for x in bs)
        else:
            bs_str = bs.to_gaussian()

        if core is not None:
            core = "core,{0:s}\n".format(",".join([str(x) for x in core]))
        else:
            core = ''

        subs = {
            'geometry' : mol.molpro_rep(),
            'basisset' : bs_str + "\n",
            'core' : core,
        }

        with open(self.fname, 'w') as inp:
            inp.write(temp.substitute(subs))

class Psi4OutputParser(object):

    '''Class for parsing molro output files'''

    def __init__(self, out=None):

        self.output = out
        self.outexists()

    def outexists(self):

        '''Check if the out file exists.'''

        if os.path.exists(self.output):
            return True
        else:
            sys.exit("Psi4 out file: {0:s} doesn't exist in {1:s}".format(
                     self.output, os.getcwd()))

    def get_hf_total_energy(self):

        '''Return the total HF energy.'''

        with open(self.output, 'r') as out:
            data = out.read()

        hfre = re.compile(r'@RHF Final Energy:\s+(?P<energy>\-?\d+\.\d+)', flags=re.M)
        match = hfre.search(data)
        if match:
            return float(match.group("energy"))

    def get_ccsd_total_energy(self):

        '''Return the total CCSD energy.'''

        with open(self.output, 'r') as out:
            data = out.read()

        ccre = re.compile(r'!CCSD total energy\s+(?P<energy>\-?\d+\.\d+)', flags=re.M)
        match = ccre.search(data)
        if match:
            return float(match.group("energy"))

    def get_ccsdt_total_energy(self):

        '''Return the total CCSD(T) energy.'''

        with open(self.output, 'r') as out:
            data = out.read()

        ccre = re.compile(r'!CCSD\(T\) total energy\s+(?P<energy>\-?\d+\.\d+)', flags=re.M)
        match = ccre.search(data)
        if match:
            return float(match.group("energy"))

    def get_cisd_total_energy(self):

        '''Return the total CISD energy.'''

        with open(self.output, 'r') as out:
            data = out.read()

        cire = re.compile(r'\s+\* CISD total energy:\s*(?P<energy>\-?\d+\.\d+)', flags=re.M)
        match = cire.search(data)
        if match:
            return float(match.group("energy"))

    def get_fci_total_energy(self):

        '''Return the total HF energy.'''

        with open(self.output, 'r') as out:
            data = out.read()

        fcire = re.compile(r'!FCI STATE \d+\.\d+ Energy\s+(?P<energy>\-?\d+\.\d+)', flags=re.M)
        match = fcire.search(data)
        if match:
            return float(match.group("energy"))

    def get_variable(self, rawstring):
        with open(self.output, 'r') as out:
            data = out.read()

        genre = re.compile(rawstring, flags=re.M)
        match = genre.search(data)
        if match:
            return float(match.group(1))

    def accomplished(self):

        '''Check if the job terminated succesfully.'''

        with open(self.output, 'r') as out:
            data = out.read()

        successre = re.compile(r'\*\*\* PSI4 exiting successfully.', flags=re.I)

        match = successre.search(data)
        if match:
            return True
        else:
            return False
