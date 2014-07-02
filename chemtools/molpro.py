import basisset as bas
from code import Code
from subprocess import Popen, PIPE
import os
import re

class Molpro(Code):
    '''
    Generic class holding a Molpro object.
    '''

    def __init__(self, **kwargs):
        super(Molpro, self).__init__(**kwargs)

    def write_input(self, inpfile=None, core=None, bs=None, inpdata=None, mol=None):

        inpdata = re.sub('geometry', mol.molpro_rep(), inpdata, flags=re.I)
        inpdata = re.sub('basis', bas.write_molpro_basis(bs), inpdata, flags=re.I)
        if core:
            inpdata = re.sub("core","core,{0:s}\n".format(",".join([str(x) for x in core])), inpdata, flags=re.I)
        else:
            inpdata = re.sub("core","", inpdata, flags=re.I)

        with open(inpfile, 'w') as inp:
            inp.write(inpdata)

    def run(self, inpfile):
        '''
        Run a single molpro job interactively - without submitting to the queue.
        '''

        if "-o" in self.runopts:
            outfile = self.runopts[self.runopts.index("-o") + 1]
        else:
            outfile = os.path.splitext(inpfile)[0] + ".out"
        errfile = os.path.splitext(outfile)[0] + ".err"
        opts = []
        opts.extend([self.execpath, inpfile] + self.runopts)

        process = Popen(opts, stdout=PIPE, stderr=PIPE)
        out, err = process.communicate()
        ferr = open(errfile, 'w')
        ferr.write(out)
        ferr.write("{0:s}\n{1:^80s}\n{0:s}\n".format("="*80, "Error messages:"))
        ferr.write(err)
        ferr.close()

        self.outfile = outfile

    def run_multiple(self, inputs):
        '''
        Run a single molpro job interactively - without submitting to the queue.
        '''

        procs = []
        outputs = [os.path.splitext(inp)[0] + ".out" for inp in inputs]
        for inpfile, outfile in zip(inputs, outputs):
            opts = []
            opts.extend([self.execpath, inpfile] + self.runopts)
            out = open(outfile, 'w')
            process = Popen(opts, stdout=out, stderr=out)
            out.close()
            procs.append(process)

        for p in procs: p.wait()

        return outputs

    def parse(self, method, objective, regexp=None):
        '''
        Parser molpro output file to get the objective.
        '''

        parser = MolproOutputParser(self.outfile)

        if objective == "total energy":
            if method == "hf":
                return parser.get_hf_total_energy()
            elif method == "cisd":
                return parser.get_cisd_total_energy()
        elif objective == "correlation energy":
                return parser.get_cisd_total_energy() - parser.get_hf_total_energy()
        elif objective == "regexp":
            return parser.get_variable(regexp)
        else:
            sys.exit("<parse>: unknown objective {0:s}".format(objective))

    def parse_tote(self, outfile):
        parser = MolproOutputParser(outfile)
        return parser.get_cisd_total_energy()

    def isok(self, outfile=None):
        '''
        Check if molpro job finished without errors.
        '''

        if outfile:
            parser = MolproOutputParser(outfile)
        else:
            parser = MolproOutputParser(self.outfile)
        return parser.terminatedOK()

class MolproOutputParser(object):

    '''Class for parsing molro output files'''

    def __init__(self, out=None):

        self.output = out
        self.outexists()

    def outexists(self):

        '''Check if the out file exists.'''

        if os.path.exists(self.output):
            return True
        else:
            sys.exit("Molpro out file: {0:s} doesn't exist in {1:s}".format(
                     self.output, os.getcwd()))

    def get_hf_total_energy(self):

        '''Return the total HF energy.'''

        with open(self.output, 'r') as out:
            data = out.read()

        hfre = re.compile(r'!RHF STATE \d+\.\d+ Energy\s+(?P<energy>\-?\d+\.\d+)', flags=re.M)
        match = hfre.search(data)
        if match:
            return float(match.group("energy"))

    def get_cisd_total_energy(self):

        '''Return the total CISD energy.'''

        with open(self.output, 'r') as out:
            data = out.read()

        cire = re.compile(r'!CISD total energy\s+(?P<energy>\-?\d+\.\d+)', flags=re.M)
        match = cire.search(data)
        if match:
            return float(match.group("energy"))

    def get_variable(self, rawstring):
        with open(self.output, 'r') as out:
            data = out.read()

        genre = re.compile(rawstring, flags=re.M)
        match = genre.search(data)
        if match:
            return float(match.group(1))

    def terminatedOK(self):

        '''Check if the job terminated succesfully.'''

        with open(self.output, 'r') as out:
            data = out.read()

        errorre = re.compile(r'\s*error', flags=re.I)

        match = errorre.search(data)
        if match:
            return False
        else:
            return True

