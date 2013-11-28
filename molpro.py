import re
import os

class OutputParser(object):

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

    def terminatedOK(self):

        '''Check if the job terminated succesfully.'''

        with open(self.output, 'r') as out:
            data = out.read()

        errorre = re.compile(r'\?\s*error', flags=re.I)

        match = errorre.search(data)
        if match:
            return False
        else:
            return True

