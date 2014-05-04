
from code import Code
from gamessus import GamessParser
from subprocess import Popen

import os

class Dmft(Code):

    def __init__(self, logfile):

        '''
        Initialize the dmft class object, logfile is the gamess log file name.
        '''

        super(Dmft, self).__init__(**kwargs)
        self.logfile = logfile

    def write_input(self, functional=None, a1=None, b1=None, a2=None, b2=None, printlevel=None, analyze=None):

        '''Write dmft input based on information in the gamess log file.'''

        inputdata = self.parse_gamess()

        if functional:
            inputdata["functional"]  = functional
        if a1:
            inputdata["a1"]          = a1
        if b1:
            inputdata["b1"]          = b1
        if a2:
            inputdata["a2"]          = a2
        if b2:
            inputdata["b2"]          = b2
        if printlevel:
            inputdata["print_level"] = printlevel
        if analyze:
            inputdata["analyze"]     = ".true."

        filebase        = "dmft_fun{0:<s}".format(str(inputdata["functional"]))
        self.inputfile  = filebase + ".inp"
        self.outputfile = filebase + ".out"
        self.nosfile    = filebase + ".nos"

        inp = open(self.inputfile, 'w')
        inp.write("&input\n")
        for key in sorted(inputdata.keys()):
            inp.write("\t{:s}={:s}\n".format(key, str(inputdata[key])))
        inp.write("/\n")
        inp.close()

    def parse_gamess(self):

        '''Parse gamess-us log file to get the neccessary data to write dmft input
            and set defaults.'''

        gp = GamessParser(self.logfile)

        inputdata = {
            "a1"             : 0.0,
            "b1"             : 0.0,
            "a2"             : 0.0,
            "b2"             : 0.0,
            "analyze"        : ".false.",
            "dictfile"       : "'{}'".format(gp.dictionary),
            "exportNOs"      : ".false.",
            "functional"     : 2,
            "loadNOs"        : 1,
            "loadIntegrals"  : 0,
            "nbasis"         : gp.get_number_of_mos(),
            "nuclear_charge" : gp.get_electrons()+gp.get_charge(),
            "print_level"    : 1,
            "restart"        : 0,
            "title"          : "'{}'".format(os.path.splitext(self.logfile)[0]),
            "total_charge"   : gp.get_charge(),
            "twointfile"     : "'{}'".format(gp.twoemofile),
        }
        return inputdata


    def run(self):

        '''Run a single dmft job.'''

        out = open(self.outputfile, 'w')
        process = Popen([self.execpath, self.inputfile], stdout=out, stderr=out)
        process.wait()
        out.close()
