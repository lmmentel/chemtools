'''
Module for handling Gamess-US related jobs,:
    Gamess          : running and submitting jobs, writing inputs,
    GamessInpParser : parsing the input file,
    GamessLogParser : parsing the output file,
    GamessDatParser : parsing data from the gamess PUNCH (*.dat) file
    GamessReader    : reading gamess bianry files.
'''

from __future__ import print_function

from chemtools.code import Code
from subprocess import Popen
from collections import OrderedDict
import numpy as np
import os
import re
import sys

# fortran modules nedded for GamessReader Class
try:
    import dictfile as df
    import twoe
except:
    pass

class Gamess(Code):

    '''Container object for Gamess-us jobs.'''

    def __init__(self, name="GamessUS", version="00", **kwargs):
        super(Gamess, self).__init__(**kwargs)
        self.name = name
        self.gamesspath = os.path.dirname(self.executable)
        self.version = version

        if os.path.isfile(os.path.join(self.gamesspath, 'ddikick.x')):
            self.ddikick = os.path.join(self.gamesspath, 'ddikick.x')

        if os.path.isfile(os.path.join(self.gamesspath, 'rungms')):
            self.rungms = os.path.join(self.gamesspath, 'rungms')
        else:
            raise IOError('Could not find "rungms" under {0:s}'.format(self.gamesspath))

    @property
    def version(self):
        """Return the version."""
        return self._version

    @version.setter
    def version(self, value):
        """Check if the version exist and set it if it does."""

        versions = []
        for item in os.listdir(self.gamesspath):
            match = re.match(r'[a-z]+\.([0-9]+)\.x', item)
            if match:
                versions.append(match.group(1))

        if value in versions:
            self._version = value
        else:
            raise IOError('GamessUS version {0:s} not found in {1:s}'.format(value, self.gamesspath))

    def remove_dat(self, inpfile):
        '''
        Remove the gamess dat file if it exists in the scratch directory.
        '''

        datfile = os.path.splitext(inpfile)[0] + ".dat"
        if os.path.exists(os.path.join(self.scratch, datfile)):
            os.remove(os.path.join(self.scratch, datfile))


    def run(self, inpfile, logfile=None, nproc=1, remove_dat=True):
        '''
        Run a single gamess job interactively - without submitting to the
        queue.
        '''

        if remove_dat:
            self.remove_dat(inpfile)

        if logfile is None:
            logfile = os.path.splitext(inpfile)[0] + ".log"

        out = open(logfile, 'w')
        process = Popen([self.executable, inpfile, self.version, str(nproc)], stdout=out, stderr=out)
        process.wait()
        out.close()
        return logfile

    def run_multiple(self, inputs):
        pass

    def accomplished(self, outfile):
        '''
        Return True if Gamess(US) job finished without errors.
        '''

        parser = GamessLogParser(outfile)
        return parser.accomplished()

    def parse(self):
        pass

    def write_input(self, inpfile=None, core=None, bs=None, inpdata=None, mol=None):
        '''
        Write a file containing gamess input.
        '''

        if isinstance(bs, list):
            basstr = "".join(x.write_gamess() for x in bs)
        else:
            basstr = bs.write_gamess()
        inpdata = re.sub('geometry', mol.gamess_rep(), inpdata, flags=re.I)
        inpdata = re.sub('basis', basstr, inpdata, flags=re.I)
        if core:
            inpdata = re.sub("core","{0:s}\n".format(str(core)), inpdata, flags=re.I)
        else:
            inpdata = re.sub("core","", inpdata, flags=re.I)

        with open(inpfile, 'w') as inp:
            inp.write(inpdata)

class GamessInput(object):
    '''
    A class for parsing and writing gamess-us input files.
    '''

    def __init__(self, fname=None, inpdata=None):
        '''
        Initialize the class.
        '''

        self.fname = fname
        if inpdata is not None:
            self.inpdata = inpdata
        else:
            self.inpdata = {}
        self.end = " $end\n"
        # not nested groups of input blocks (not parsed into a dict of dicts)
        self._notnested = ["$data", "$vec", "$ecp"]

    def parse(self):
        '''
        Parse gamess input file into a dictionary of dictionaries, where the
        highest level entries are gamess namelist fileds and that contain
        dictionaries of options. All key are converted to lowercase.
        '''

        with open(self.fname, 'r') as finp:
            contents = finp.read()
        return self.parse_from_string(contents)

    def parse_from_string(self, inpstr):
        '''
        Parse gamess input file into a dictionary of dictionaries, where the
        highest level entries are gamess namelist fileds and that contain
        dictionaries of options. All key are converted to lowercase. For
        example if the following input was parsed:
        '''

        pat = re.compile(r'(?P<block>\$[a-zA-Z]{3,6})\s+(?P<entries>.*?)\$END', flags=re.DOTALL)

        iterator = pat.finditer(inpstr)
        for match in iterator:
            if match.group("block").lower() not in self._notnested:
                self.inpdata[match.group("block").lower()] = {}
                fields = [s.strip() for s in match.group("entries").split("\n")]
                for field in fields:
                    if not field.startswith("!"):
                        for line in field.split():
                            key, value = line.split("=")
                            self.inpdata[match.group("block").lower()][key.lower()] = value
            elif match.group("block").lower() == "$data":
                self.inpdata["$data"] = self.parse_data(match.group("entries"))
            elif match.group("block").lower() in ["$vec", "$ecp"]:
                self.inpdata[match.group("block").lower()] = match.group("entries")
        return True

    def parse_data(self, datastr, parse_basis=False):
        '''
        Parse $DATA block specified in the gamess input file. The parser
        assumes that in $DATA input the atom positions are specified using
        cartesian coordiantes and returns a list of dictionaries with parsed
        information about each specified atom.

        Args:
            datastr (str)
                string with the contents of the $data block

        Returns:
            datadict (dict)
        '''

        block = re.compile(r'(?P<label>[a-zA-Z]{1,2}[0-9]{0,2})\s*'
                  +r'(?P<atomic>\d+\.\d+)'
                  +'(?P<xyz>(\s+\-?\d+\.\d+){3})\s*'
                  +'(?P<basis>.*?)\n\s*\n', flags=re.S)

        datadict = dict()

        datadict["title"] = datastr.split('\n')[0].strip()
        datadict["group"] = datastr.split('\n')[1].strip()
        datadict["atoms"] = list()

        itfound = block.finditer(datastr)
        for m in itfound:
            datadict["atoms"].append({'label'  : m.group('label'),
             'atomic' : m.group('atomic'),
             'xyz'    : tuple(float(x) for x in m.group('xyz').split()),
             'basis'  : m.group('basis'),})
        return datadict

    def write_input(self, inpdict, fname=None, mol=None, bs=None):
        '''
        Write a gamess input file under the name <inpfile> based on the
        information fstored in the dictionary <inpdict>.

        Args:
            inpfile (str): name of the file to be written,
            inpdict (dict): dictionary with the gamess input specification.
        '''

        if not isinstance(inpdict, dict):
            raise TypeError("expected a dictionary but got {0:s}".format(type(inpdict)))

        inpstr = ""

        # write nested namelist groups
        for key, value in sorted(inpdict.items()):
            if key not in  self._notnested:
                inpstr += " {0:<s}\n".format(key)
                for kkey, vvalue in sorted(value.items()):
                    inpstr += "    {k:s}={v:s}\n".format(k=kkey, v=str(vvalue))
                inpstr += self.end
        #write $data card
        inpstr += self.write_data(inpdict["$data"], mol, bs)

        with open(self.fname, "w") as finp:
            finp.write(inpstr)

    def write_data(self, datadict, mol=None, bs=None):
        '''
        Args:
        =====
        bs: (dict)
            dictionary of BasisSet objects
        '''

        data = ""
        data += " {0:s}\n".format("$data")
        data += "{0:s}\n".format(datadict["title"])
        data += "{0:s}\n\n".format(datadict["group"])
        if mol is not None:
            for atom in mol.unique():
                data += atom.gamess_rep()
                if bs is not None:
                    data += bs[atom.symbol].write_gamess()
        data += self.end
        return data

    def set_gamess_input(self, dinp, mol, bs, code, core):

        if "$contrl" in dinp.keys():
            dinp["$contrl"]["icharg"] = mol.charge
            dinp["$contrl"]["mult"]   = mol.multiplicity
            if mol.multiplicity == 1:
                dinp["$contrl"]["scftyp"] = "rhf"
            elif mol.multiplicity > 1:
                dinp["$contrl"]["scftyp"] = "rohf"
            if code["method"].lower() == "hf":
                dinp["$contrl"] = {key:value for key, value in dinp["$contrl"].items() if key != "cityp"}
        else:
            sys.exit("no $contrl group in the gamess input string")
        if "$cidet" in dinp.keys():
            dinp["$cidet"]["nact"] = bs.get_no_functions(bs) - core
            dinp["$cidet"]["ncore"] = core
            dinp["$cidet"]["nels"] = mol.electrons - core*2
        if "$ormas" not in dinp.keys():
            dinp["$ormas"] = {}
        if code["method"] == "cisd":
            # sum of singly and doubly occupied  orbitals
            noccupied = mol.multiplicity - 1 + (mol.electrons-(mol.multiplicity - 1))/2
            dinp["$ormas"]["nspace"] = 2
            dinp["$ormas"]["mstart(1)"] = ','.join([str(core+1),str(noccupied+1)])
            dinp["$ormas"]["mine(1)"] = ",".join([str((mol.electrons-core*2)-2), str(0)])
            dinp["$ormas"]["maxe(1)"] = ",".join([str(mol.electrons-core*2), str(2)])
        elif code["method"] == "mrcisd":
            dinp["$ormas"]["nspace"] = 3
            dinp["$ormas"]["mstart(1)"] = ','.join([str(core+1),str(mol.electrons-core+1)])
            dinp["$ormas"]["mine(1)"] = ",".join([str(mol.electrons-2), str(0)])
            dinp["$ormas"]["maxe(1)"] = ",".join([str(mol.electrons), str(2)])

        return dinp

    def write_with_vec(self, inpfile, vecstr):
        '''
        Write new gamess input based on exisiting input and natural orbitals
        from $VEC section in PUNCH file.


        Write a new input file for gamess based on previously prased and/or
        constructed dictionary containing input specification and append
        orbitals from a previous run as starting orbitals. The starting
        orbitals should be stored in an ASCII PUNCH file whose name is given
        under "datfile" variable.

        Args:
            inpfile (string)
                name of the input file to be created,
            inpdict (dictionary)
                dict containing the inpu specification, either parsed from
                using "GamessInpParser.parse" or constructed by hand,
            datfile (string)
                name of the PUNCH file that stores the orbitals genereated in a
                previous run,
        '''

        self.inpdata['$contrl']['scftyp'] = 'none'
        self.inpdata['$guess'] = {}
        self.inpdata['$guess']['guess'] = 'moread'
        self.inpdata['$guess']['norb'] = str(get_naos_nmos(vecstr)[1])

        inp = open(inpfile, "a")
        inp.write("\n $vec\n")
        inp.write(vecstr)
        inp.write(" $end\n")
        inp.close()

    def print_inpdict(self, inpdict):
        '''
        Neat print of the dictionary with parsed input
        '''

        for key, value in sorted(inpdict.items()):
            if isinstance(value, dict):
                print(key)
                for kkey, vvalue in sorted(value.items()):
                    print("\t{0:<10s} : {1:}".format(kkey, vvalue))
            else:
                print(key, '\n', value)

    def parse_gamess_basis(self, basis_str):
        '''
        Parse the basis set into a list of dictionaries from a string in
        gamess format.
        '''

        pat = re.compile(r'^\s*(?P<shell>[SPDFGHIspdfghi])\s*(?P<nf>[1-9]+)')

        bslines = basis_str.split("\n")

        functions = OrderedDict()

        for i, line in enumerate(bslines):
            match = pat.search(line)
            if match:
                shell, nf = match.group("shell").lower(), match.group("nf")
                exps, indxs, coeffs = self.parse_function(bslines[i+1:i+int(nf)+1])
                if shell in functions.keys():
                    lasti = len(functions[shell]["exponents"])
                    functions[shell]["exponents"].extend(exps)
                    functions[shell]["contractedfs"].append({"indices" : [lasti + i - 1 for i in indxs],
                                                            "coefficients" : coeffs})
                else:
                    functions[shell] = OrderedDict()
                    functions[shell]["exponents"] = exps
                    functions[shell]["contractedfs"] = list()
                    functions[shell]["contractedfs"].append({"indices" : [i - 1 for i in indxs],
                                                            "coefficients" : coeffs})
        return functions

    @staticmethod
    def parse_function(los):
        '''
        Parse a basis set function information from list of strings into
        three lists containg:
            exponents
            indices
            coefficients

        Remeber that python doesn't recognise the '1.0d-3' format where 'd' or 'D' is used
        to the regex subsitution has to take care of that.
        '''

        real = re.compile(r'[dD]')

        indxs = [int(item.split()[0]) for item in los]
        exps = [float(real.sub('E', item.split()[1])) for item in los]
        coeffs = [float(real.sub('E', item.split()[2])) for item in los]

        return (exps, indxs, coeffs)

class GamessLogParser(object):
    '''
    Methods for parsing gamess-us log file.
    '''

    def __init__(self, log):
        self.logfile = log

    @property
    def logfile(self):
        return self._logfile

    @logfile.setter
    def logfile(self, value):
        '''
        Check if the logfile exists first.
        '''
        if os.path.exists(value):
            self._logfile = value
        else:
            raise ValueError("File: {} does not exist".format(value))

    def logexists(self):
        '''
        Check if the log file exists.
        '''

        if os.path.exists(self.logfile):
            return True
        else:
            sys.exit("Gamess log file: {0:s} doesn't exist in {1:s}".format(
                     self.logfile, os.getcwd()))

    def parse(self, regex):
        '''
        Parsing function based on querries in the form of regular expressions that
        return the match object.
        '''

        cpatt = re.compile(regex)
        with open(self.logfile, 'r') as log:
            content = log.read()
        return cpatt.search(content)

    def accomplished(self):

        '''Check if a job teminated normally.'''

        regex = r'TERMINATED NORMALLY'
        match = self.parse(regex)
        if match:
            return True
        else:
            return False

    def get_version(self):
        '''
        Get the version of the GAMESS(US) package.
        '''

        regex = r'.*\s+GAMESS VERSION =\s*(.*?)\s*\*'
        match = self.parse(regex)
        if match:
            return match.group(1)
        else:
            return None

    def get_charge(self):

        '''Get total charge.'''

        regex = r'CHARGE OF MOLECULE\s+=\s*(?P<charge>\d+)'
        match = self.parse(regex)
        if match:
            return int(match.group("charge"))

    def get_electrons(self):

        '''Get number of electrons.'''

        regex = r'NUMBER OF ELECTRONS\s+=\s*(?P<nele>\d+)'
        match = self.parse(regex)
        if match:
            return int(match.group("nele"))

    def get_homo(self):

        '''Get the orbital index of homo orbital (indexing starts from zero).'''

        if int(self.get_electrons()) % 2 == 0:
            return int(self.get_electrons())/2 - 1
        else:
            sys.exit("open shell handling not implemented")

    def get_number_of_atoms(self):

        '''Get total number of atoms from gamess log file.'''

        regex = r'TOTAL NUMBER OF ATOMS\s+=\s*(?P<nat>\d+)'
        match = self.parse(regex)
        if match:
            return int(match.group("nat"))

    def get_number_of_aos(self):

        '''Get the number of primitive cartesian gaussian basis functions from
        Gamess log file'''

        regex = r'NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS =\s*(?P<nao>\d+)'
        match = self.parse(regex)
        if match:
            return int(match.group("nao"))

    def get_number_of_core_mos(self):

        '''Get the number of core molecular orbitals from the Gamess log file'''

        regex = r'NUMBER OF CORE MOLECULAR ORBITALS\s*=\s*(?P<ncore>\d+)'
        match = self.parse(regex)
        if match:
            return int(match.group("ncore"))

    def get_number_of_mos(self):

        '''Get the number of molecular orbitals from Gammess log file.'''

        ispher_patt    = r'ISPHER=\s*(?P<ispher>\-?\d{1}).*'
        var_space_patt = r'.*VARIATION SPACE IS\s*(?P<nmo>\d+).*'
        c_ispher       = re.compile(ispher_patt)
        c_var_space    = re.compile(var_space_patt)

        with open(self.logfile, 'r') as log:
            lines = log.read()

        match = c_ispher.search(lines)
        if match:
            ispher = int(match.group("ispher"))

        if ispher == -1:
            return self.get_number_of_aos()
        elif ispher == 1:
            match = c_var_space.search(lines)
            if match:
                n_mo = int(match.group("nmo"))
                return n_mo
        else:
            sys.exit("wrong ispher found: {0:d}".format(ispher))

    def get_linear_deps(self):

        '''Get number of linearly dependent combinations dropped.'''

        regex = r'NUMBER OF LINEARLY DEPENDENT MOS DROPPED=\s*(?P<lindep>\d+)'
        match = self.parse(regex)
        if match:
            return int(match.group('lindep'))

    def get_scf_type(self):

        '''Get the information on SCFTYP used in the gamess job.'''

        regex = r'.*SCFTYP=(?P<scftyp>[A-Z]+)'
        match = self.parse(regex)
        if match:
            return match.group("scftyp")

    def get_cc_type(self):

        '''Get the information on CCTYP used in the gamess job.'''

        regex = r'.*CCTYP =(?P<cctyp>[A-Z\(\)\-]+)'
        match = self.parse(regex)
        if match:
            return match.group("cctyp")
        else:
            return None

    def get_ci_type(self):

        '''Get the information on CITYP used in the gamess job.'''

        regex = r'.*CITYP =(?P<cityp>[A-Z]+)'
        match = self.parse(regex)
        if match:
            return match.group("cityp")

    def get_mplevel(self):

        '''Get the information on MPLEVL used in the gamess job.'''

        regex = r'.*MPLEVL=\s*(?P<mplevl>\d+)'
        match = self.parse(regex)
        if match:
            return int(match.group("mplevl"))

    def get_hf_total_energy(self):

        '''Return the total HF energy.'''

        regex = r'FINAL R[O]*HF ENERGY IS\s+(?P<energy>\-?\d+\.\d+)'
        match = self.parse(regex)
        if match:
            return float(match.group("energy"))

    def get_cc_total_energy(self):
        '''
        Independently of CCTYP value return the "HIGHEST LEVEL RESULT" total energy.
        '''

        regex = r'COUPLED-CLUSTER ENERGY E(.*?)\s*=\s*(?P<energy>\-?\d+\.\d+)'
        match = self.parse(regex)
        if match:
            return float(match.group("energy"))

    def get_energy_components(self, method):
        '''
        Read the summary of the energies printed in the gamess log file at the
        property section corresponding to a particular "method".
        '''

        if method.lower() in ["hf", "scf", "hfscf"]:
            if self.get_scf_type() != "NONE":
                header = 'PROPERTY VALUES FOR THE {0:<5s} SELF-CONSISTENT FIELD WAVEFUNCTION'.format(self.get_scf_type())
            else:
                sys.exit("No HF calculation was performed, check the log file once again.")
        elif method.lower() in ["ci"]:
            if self.get_ci_type().lower() in ["guga", "ormas", "fsoci"]:
                header  = '{0:<5s} CI PROPERTIES'.format(self.get_ci_type())
        else:
            sys.exit("Wrong method in <get_energy_components>: {0:s}".format(method))

        with open(self.logfile, 'r') as log:
            data = log.readlines()

        return self.parse_pairs(self.slice_after(data, header, 22))

    def get_lz_values(self):
        '''
        Get values of Lz operator for every orbital and translate them into
        labels of orbitals.
        '''

        lz_patt = re.compile(r'\s*MO\s*(?P<index>\d*)\s*\(\s*(?P<shell>\d+)\s*\)\s*HAS LZ\(WEIGHT\)=\s*(?P<lz>-?\d*\.\d+)')
        with open(self.logfile, 'r') as log:
            lz_string = self.slice_between(log.read(), 'LZ VALUE ANALYSIS FOR THE MOS', 'EIGENVECTORS')

        res = list()
        for line in lz_string.split('\n'):
            match = lz_patt.search(line)
            if match:
                index = int(match.group('index')) - 1
                shell = int(match.group('shell'))
                lz   = int(abs(float(match.group('lz'))))
                res.append({"index" : index, "shell" : shell, "lz" : lz})
        return res

    def get_ao_labels(self, orbitals=None):
        '''
        Retrieve the information about atomic basis set
        '''

        with open(self.logfile, 'r') as log:
            orb_string = self.slice_between(log.read(), 'EIGENVECTORS', 'END OF RHF CALCULATION')

        ao_patt = re.compile(r'\s*(?P<index>\d+)\s{2}(?P<symbol>[a-zA-Z]{1,2})\s*(?P<center>\d+)\s*(?P<component>[A-Z]{1,})')
        res = list()
        for line in orb_string.split('\n'):
            match = ao_patt.search(line)
            if match:
                index = int(match.group('index')) - 1
                symbol = match.group('symbol')
                center = int(match.group('center'))
                component = match.group('component')
                res.append({"index" : index, "symbol" : symbol, "center" : center, "component" : component})
        return res

    @staticmethod
    def parse_pairs(los, sep="="):
        '''
        Parse a given list of strings "los" into a dictionary based on
        separation by "sep" character and return the dictionary.
        '''
        out = []
        for line in los:
            if sep in line:
                (name, value) = line.split(sep)
                out.append((name.strip(), float(value)))
        return dict(out)

    @staticmethod
    def slice_after(seq, item, num):
        '''
        Return "num" elements of a sequence "seq" present after the item "item".
        '''
        it = iter(seq)
        for element in it:
            if item in element:
                return [next(it) for i in range(num)]

    @staticmethod
    def slice_between(string, start, end):
        istart = string.index(start)
        iend = string.index(end)
        return string[istart+len(start):iend]

class GamessFortranReader(object):
    '''
    Class for holding method for reading gamess binary files:
        $JOB.F08 : two electron integrals over AO's,
        $JOB.F09 : two electron integrals over MO's,
        $JOB.F10 : the dictionary file with 1-e integrals, orbitals etc.,
        $JOB.F15 : GUGA and ORMAS two-electron reduced density matrix,

    TODO:
        CI coefficients, and CI hamiltonian matrix elements.
    '''
    # add an option to choose which reader should be used for dictionary file
    # and sequential files, the options are:
    # - wrapped fortran code that need to be compiled and installed or
    # - native reader written in python using DictionaryFile and SequentialFile
    # classes
    def __init__(self, log):
        self.logfile    = log
        i = self.logfile.index("log")
        self.filebase   = self.logfile[:i-1]
        self.datfile    = self.filebase + ".dat"
        self.twoeaofile = self.filebase + ".F08"
        self.twoemofile = self.filebase + ".F09"
        self.dictionary = self.filebase + ".F10"
        self.rdm2file   = self.filebase + ".F15"
        self.gp         = GamessLogParser(log=self.logfile)

    def get_onee_size(self, aos=True):
        '''
        Get the size of the vector holding upper (or lower) triangle
        of a square matrix of size naos or nmos.
        '''
        if aos:
            n = self.gp.get_number_of_aos()
        else:
            n = self.gp.get_number_of_mos()
        return n*(n+1)/2

    def get_twoe_size(self):
        '''
        Get the size of the 1d vector holding upper (or lower) triangle
        of a supermatrix of size nmos (2RDM and two-electrons integrals).
        '''
        n = self.get_onee_size(aos=False)
        return n*(n+1)/2

    def read_rdm2(self, filename=None, nmo=None):

        '''Read the 2rdm from the gamess-us file'''

        # initialize numpy array to zeros
        rdm2 = np.zeros(self.get_twoe_size(), dtype=float)

        # use gamess module to read the integrals from the file -filename-
        if filename:
            if os.path.exists(filename):
                print("Reading {}".format(filename))
                twoe.integrals.readinao(rdm2, filename)
                return rdm2
            else:
                sys.exit("File '{0:s}' doesn't exist, exiting...".format(filename))
        elif os.path.exists(self.rdm2file):
            print("Reading {}".format(self.rdm2file))
            twoe.integrals.readinao(rdm2, self.rdm2file)
            return rdm2
        else:
            sys.exit("File '{0:s}' doesn't exist, exiting...".format(self.rdm2file))

    def read_twoemo(self, filename=None, nmo=None):

        '''Read the two electron integrals from the gamess-us file'''

        ints = np.zeros(self.get_twoe_size(), dtype=float)

    def read_twoemo_fortran(self, filename=None, nmo=None):

        '''Read the two electron integrals from the gamess-us file'''

        # initialize numpy array to zeros
        ints = np.zeros(self.get_twoe_size(), dtype=float)
        # use gamess module to read the integrals from the file -filename-
        if filename:
            if os.path.exists(filename):
                twoe.integrals.readinmo(ints, filename)
                return ints
            else:
                sys.exit("File '{0:s}' doesn't exist, exiting...".format(filename))
        elif os.path.exists(self.twoemofile):
            twoe.integrals.readinmo(ints, self.twoemofile)
            return ints
        else:
            sys.exit("File '{0:s}' doesn't exist, exiting...".format(self.twoemofile))


    def read_H(self):
        '''
        Read the bare nucleus hamiltonian integrals form section 11 of gamess-us
        dictionary file.
        '''
        vec = np.zeros(self.get_onee_size(), dtype=float)
        df.readit(self.dictionary, vec, 11)
        return vec

    def read_S(self):
        '''
        Read the overlap integrals form section 12 of gamess-us dictionary
        file.
        '''
        vec = np.zeros(self.get_onee_size(), dtype=float)
        df.readit(self.dictionary, vec, 12)
        return vec

    def read_T(self):
        '''
        Read the kinetic energy integrals form section 13 of gamess-us
        dictionary file.
        '''
        vec = np.zeros(self.get_onee_size(), dtype=float)
        df.readit(self.dictionary, vec, 13)
        return vec


    def read_occupations(self):
        '''
        Get the natural orbital occupation numbers from section 21 of the
        gamess-us dictionary file.
        '''
        vec = np.zeros(self.gp.get_number_of_mos(), dtype=float)
        df.readit(self.dictionary, vec, 21)
        return vec

    def read_mos(self):
        '''
        Read the Hartree-Fock MO's from the section 15 of the gamess dictionary
        file.
        '''

        mat = np.zeros(self.gp.get_number_of_aos()*self.gp.get_number_of_mos(), dtype=float)
        df.readit(self.dictionary, mat, 15)
        mat = mat.reshape((self.gp.get_number_of_aos(), self.gp.get_number_of_mos()), order='F')
        return mat

    def read_orbital_energies(self):
        '''
        Read orbital energies (HF) from the section 17 of the gamess dictionary
        file.
        '''

        vec = np.zeros(self.gp.get_number_of_mos(), dtype=float)
        df.readit(self.dictionary, vec, 17)
        return vec

    def read_nos(self):
        '''
        Read Natural Orbitals form section 19 of the gamess dictionary file.
        '''

        mat = np.zeros(self.gp.get_number_of_aos()*self.gp.get_number_of_mos(), dtype=float)
        df.readit(self.dictionary, mat, 19)
        mat  = mat.reshape((self.gp.get_number_of_aos(), self.gp.get_number_of_mos()), order='F')
        return mat

    @staticmethod
    def factor(i, j, k, l):
        '''
        Based on the orbitals indices return the factor that takes into account
        the index permutational symmetry.
        '''
        if i == j and k == l and i == k:
            fijkl = 1.0
        elif i == j and k == l:
            fijkl = 2.0
        elif (i == k and j == l) or (i == j and i == k) or (j == k and j == l) or (i == j or k == l):
            fijkl = 4.0
        else:
            fijkl = 8.0
        return fijkl

    # this function should be moved somewhere else


    def print_twoe(self, twoe, nbf):
        '''Print the two-electron integrals.'''
        ij=0
        for i in xrange(nbf):
            for j in xrange(i+1):
                ij += 1
                kl = 0
                for k in xrange(nbf):
                    for l in xrange(k+1):
                        kl += 1
                        if ij >= kl:
                            if abs(twoe[self.ijkl(i,j,k,l)]) > 1.0e-10:
                                print("{0:3d}{1:3d}{2:3d}{3:3d} {4:25.14f}".format(
                                    i, j, k, l, twoe[self.ijkl(i,j,k,l)]))

class GamessDatParser(object):

    def __init__(self, datfile):
        self.datfile = datfile

    @property
    def datfile(self):
        return self._datfile

    @datfile.setter
    def datfile(self, value):
        if os.path.exists(value):
            self._datfile = value
        else:
            raise ValueError("File: {} does not exist".format(value))

    def get_occupations(self):
        '''
        Parse the occupation numbers from the ascii PUNCH file (*.dat).
        '''

        with open(self.datfile, 'r') as dat:
            data = dat.read()

        no_patt = re.compile(r'\$OCC(NO)?(?P<occ>.*?)\$END', flags=re.DOTALL)
        match = no_patt.search(data)
        nooc = list()
        if match:
            for line in match.group("occ").split('\n'):
                nooc.extend([float(x) for x in line.split()])
            return np.asarray(nooc)
        else:
            raise ValueError('No section with occupation numbers found.')

    def parse_data(self):
        '''
        Parse $DATA block from the dat file.
        '''

        with open(self.datfile, 'r') as dat:
            datastr = self.find_between(dat.read(), '$DATA', '$END')
        gip = GamessInput()
        return gip.parse_data(datastr.lstrip(' \n'))

    def get_orbitals(self, method):
        '''
        Parse the natural orbitals from the ascii PUNCH file (*.dat).
        Args:
            method (str)
                acceptable values are:
                    for scf orbitals: "scf", "hf", "rhf", "rohf", "uhf", "gvb"
                    for mcscf orbitals: "mcscfmos", "mcscfnos"
                    for ci orbitals: "ci", "aldet", "fsoci", "guga", "genci", "ormas""
        Returns:
            orbs (numpy.array)
                2D numpy array of shape (naos, nmos) containing ortbials
                expansion coefficients
        '''

        with open(self.datfile, 'r') as dat:
            data = dat.read()

        if method.lower() in ['scf', 'hf', 'rhf', 'rohf', 'uhf', 'gvb']:
            try:
                orbi = data.index('ORBITALS ---')
            except:
                raise ValueError('SCF orbitals header not found, check dat file')
        elif method.lower() in ['ci', 'aldet', 'fsoci', 'guga', 'genci', 'ormas']:
            try:
                orbi = data.index('NO-S')
            except:
                raise ValueError('CI orbitals header not found, check dat file')
        elif method.lower() == 'mcscfmos':
            try:
                orbi = data.index('OPTIMIZED MCSCF MO-S')
            except:
                raise ValueError('MCSCF orbitals header not found, check dat file')
        elif method.lower() == 'mcscfnos':
            try:
                orbi = data.index('NATURAL ORBITALS OF MCSCF')
            except:
                raise ValueError('MCSCF orbitals header not found, check dat file')
        else:
            raise ValueError("Don't know what to do with: '{0:s}'".format(method))
        vecstr = self.find_between(data[orbi:], '$VEC', '$END').strip(" \n\t\r")
        if vecstr == "":
            raise ValueError("No $VEC section found for method: '{0:s}'".format(method))
        # add a single whitespace to account for the removed space in strip
        return self.parse_orbitals(" "+vecstr)

    def parse_orbitals(self, vecstr, clength=15):
        '''
        Parse dat orbitals into numpy array of the shape (naos, nmos)

        Parse orbitals given a string obtained from the $VEC section of the
        PUNCH file (*.dat) into a numpy array of the shape (naos, nmos) where
        "naos" and "nmos" are the number of atomic orbitals and number of
        molecular orbitals respectively. The shape is deduced from the last
        line of the string.

        Args:
            vecstr (str)
                string with the contents fo the gamess $VEC block
            clength (int)
                total length of the coefficient string as stored by the gamess
                format, by default gamess stores orbital coefficients in
                'e15.8' fortran format so the total length is 15.

        Returns:
            orbs (numpy.array)
                2D numpy array of shape (naos, nmos) containing ortbials
                expansion coefficients
        '''

        naos, nmos, nlines = get_naos_nmos(vecstr)
        orblines = vecstr.split('\n')
        orbs = np.zeros((naos, nmos), dtype=float)
        counter = -1
        for i in range(0, nmos):
            for j in range(0, nlines):
                counter += 1
                nitems = int(len(orblines[counter][5:]))/clength
                orbs[5*j:5*(j+1), i] = [float(orblines[counter][5+15*n:5+15*(n+1)]) for n in range(nitems)]
        return orbs


    @staticmethod
    def find_between(string, first, last):
        '''
        return a slice of the `string` between phrases `first` and `last`.
        '''
        try:
            start = string.index(first) + len(first)
            end = string.index(last, start)
            return string[start:end]
        except ValueError:
            return ""

def get_naos_nmos(vecstr, clength=15):
    '''
    Get the number of AO's and MO's from a string in $VEC block.

    Args:
        vecstr (str)
            string with the contents fo the gamess $VEC block
        clength (int)
            total length of the coefficient string as stored by the gamess
            format, by default gamess stores orbital coefficients in
            'e15.8' fortran format so the total length is 15.

    Returns:
        naos (int)
            number of atomic orbitals
        naos (int)
            number of moelcular orbitals
        nlines (int)
            number of lines per molecular orbital
    '''

    veclines = len(vecstr.split('\n'))
    lineit = iter(vecstr.split('\n'))
    nlines = 0
    while lineit.next()[:2].strip() == '1':
        nlines += 1
    if nlines == 0:
        raise ValueError("'nlines' cannot be zero, check vecstr in 'get_naos_nmos'")
    naos = 5*(nlines - 1) + len(vecstr.split('\n')[nlines-1][5:])/clength
    nmos = veclines/nlines
    return naos, nmos, nlines

def take(seq, num):
    '''
    Iterate over a sequence "seq" "num" times and return the list of the
    elements iterated over.
    '''
    return [next(seq) for _ in range(num)]


#"""
#BinaryFile: A class for accessing data to/from large binary files
#=================================================================
#
#The data is meant to be read/write sequentially from/to a binary file.
#One can request to read a piece of data with a specific type and shape
#from it.  Also, it supports the notion of Fortran and C ordered data,
#so that the returned data is always well-behaved (C-contiguous and
#aligned).
#
#This is also seeking capable.
#
#:Author:   Francesc Alted
#:Contact:  faltet@pytables.org
#:Created:  2010-03-18
#:Acknowledgment: Funding for the development of this code is provided
#through the Norwegian Research Council VAUUAV project #184724, 2010
#
#"""

class BinaryFile(object):

    """Class representing a binary file (C or Fortran ordered)."""

    def __init__(self, filename, mode="r", order="fortran"):
        """Open the `filename` for write/read binary files.

        The `mode` can be 'r', 'w' or 'a' for reading (default),
        writing or appending.  The file will be created if it doesn't
        exist when opened for writing or appending; it will be
        truncated when opened for writing.  Add a '+' to the mode to
        allow simultaneous reading and writing.

        `order` specifies whether the file is is written in 'fortran'
        or 'c' order.
        """
        self.mode = mode + "b"
        self.file = open(filename, mode=self.mode, buffering=1)
        """The file handler."""
        if order not in ['fortran', 'c']:
            raise ValueError("order should be either 'fortran' or 'c'.")
        self.order = order
        """The order for file ('c' or 'fortran')."""


    def read(self, dtype, shape=(1,)):
        """Read an array of `dtype` and `shape` from current position.

        `shape` must be any tuple made of integers or even () for scalars.

        The current position will be updated to point to the end of
        read data.
        """
        if not hasattr(dtype, "itemsize"):
            dtype = np.dtype(dtype)
        if type(shape) is int:
            shape = (shape,)
        if type(shape) is not tuple:
            raise ValueError("shape must be a tuple")
        length = dtype.itemsize
        if shape is not ():
            length *= np.array(shape).prod()

        # Correct the shape in case dtype is multi-dimensional
        shape = shape + dtype.shape

        if shape in (1, (1,)):
            order = "c"
        else:
            order = self.order

        # Read the data from file
        data = self.file.read(length)
        if len(data) < length:
            raise EOFError("Asking for more data than available in file.")
        # Convert string into a regular array
        data = np.ndarray(shape=shape, buffer=data, dtype=dtype.base)

        # If original data file is in fortran order, reverse the shape first
        if order == "fortran":
            shape = [i for i in shape[::-1]]
        data = data.reshape(shape)

        if shape == ():
            # Retrieve the scalar out of the 0-dim array
            data = data[()]
        else:
            # If original data file is in fortran order, do a transpose.
            # As the shape was reversed previously, we get the original
            # shape again.
            if order == "fortran":
                data = data.transpose().copy()
            # Do an additional copy just in case the array is not
            # wel-behaved (i.e., it is not aligned or not contiguous).
            elif not data.flags.behaved:
                data = data.copy()
        return data


    def write(self, arr):
        """Write an `arr` to current position.

        The current position will be updated to point to the end of
        written data.
        """
        # Transpose data if case we need to
        if (self.order == "fortran") != (arr.flags.fortran):
            arr = arr.transpose().copy()
        # Write the data to file
        self.file.write(arr.data)


    def seek(self, offset, whence=0):
        """Move to new file position.

        Argument offset is a byte count.  Optional argument whence
        defaults to 0 (offset from start of file, offset should be >=
        0); other values are 1 (move relative to current position,
        positive or negative), and 2 (move relative to end of file,
        usually negative, although many platforms allow seeking beyond
        the end of a file).  If the file is opened in text mode, only
        offsets returned by tell() are legal.  Use of other offsets
        causes undefined behavior.
        """
        self.file.seek(offset, whence)


    def tell(self):
        "Returns current file position, an integer (may be a long integer)."
        return self.file.tell()

class SequentialFile(BinaryFile):

    def __init__(self, filename):
        '''
        Initialize the class with the superclass method.
        '''
        super(SequentialFile, self).__init__(filename)
        log = os.path.splitext(filename)[0] + ".log"
        glp = GamessLogParser(log=log)
        self.nao = glp.get_number_of_aos()
        self.nmo = glp.get_number_of_mos()
        self.core = glp.get_number_of_core_mos()
        if self.nmo < 255:
            self.large_labels = False
        else:
            self.large_labels = True

    def ijkl(self, i, j, k, l):
        '''
        Based on the four orbital indices i,j,k,l return the address
        in the 1d vector.
        '''
        ij = max(i, j)*(max(i, j) - 1)/2 + min(i, j)
        kl = max(k, l)*(max(k, l) - 1)/2 + min(k, l)
        return max(ij, kl)*(max(ij, kl) - 1)/2 + min(ij, kl) - 1

    def readseq(self, buff_size=15000, int_size=8, mos=False, skip_first=False):
        '''
        Read FORTRAN sequential unformatted file with two-electron quantities:
            * two electron integrals over AO's: .F08 file
            * two electron integrals over MO's: .F09 file
            * elements of the two particle density matrix: .F15 file

        Input
        =====
            filename (str)
                name of the file to read,
            buffSize (int)
                size of the buffer holding values to be read, in gamessus it is
                stored under "NINTMX" variable and in Linux version is equal to
                15000 which is the default value,
            large_labels (bool)
                a flag indicating if large labels should were used, if largest
                label "i" (of the MO) is i<255 then large_labels should be False
                (case "LABSIZ=1" in gamessus), otherwise set to True (case
                "LABSIZ=2" in gamess(us),
            skip_first (bool)
                skips the first record of the file is set to True,

        Output
        ======
            numpy 1D array holding the values
        '''

        if int_size == 4:
            if self.large_labels:
                indexBuffSize = 2*buff_size
            else:
                indexBuffSize = buff_size
        elif int_size == 8:
            if self.large_labels:
                indexBuffSize = buff_size
            else:
                indexBuffSize = (buff_size + 1)/2
        else:
            raise ValueError

        if mos:
            if skip_first:
                nmo = self.nmo-self.core
                n1 = nmo*(nmo+1)/2
                self.seek(8+8*n1)
                #hmo = self.read('f8', shape=(nmo*(nmo+1)/2))
                #pos = self.tell()
                #self.seek(pos+4)

            nt = self.nmo*(self.nmo+1)/2
            ints = np.zeros(nt*(nt+1)/2, dtype=float, order='F')
        else:
            self.seek(0)
            nt = self.nao*(self.nao+1)/2
            ints = np.zeros(nt*(nt+1)/2, dtype=float, order='F')


        int_type = np.dtype('i'+str(int_size))
        index_buffer = np.zeros(indexBuffSize, dtype=int_type, order='F')
        value_buffer = np.zeros(buff_size, dtype=float, order='F')

        pos = self.tell()
        self.seek(pos+4)
        length = 1
        while length > 0:
            length = self.read(int_type)
            if length > buff_size:
                raise ValueError('the read record length: {0:10d} greater that the buffer size {1:10d}'.format(int(length), buff_size))
            index_buffer = self.read(int_type, shape=(indexBuffSize, ))
            value_buffer = self.read('f8', shape=(buff_size, ))

            for m in range(1, abs(length)+1):
                if int_size == 4:
                    if self.large_labels:
                        label1 = int(index_buffer[2*m-1])
                        label2 = int(index_buffer[2*m])
                        i = label1 >> 16
                        j = label1 & 65535
                        k = label2 >> 16
                        l = label2 & 65535
                    else:
                        label = int(index_buffer[m-1])
                        i = label >> 24
                        j = label >> 16 & 255
                        k = label >>  8 & 255
                        l = label       & 255
                elif int_size == 8:
                    if self.large_labels:
                        label = int(index_buffer[m-1])
                        i = label >> 48
                        j = label >> 32 & 65535
                        k = label >> 16 & 65535
                        j = label       & 65535
                    else:
                        if m % 2 == 0:
                            label = int(index_buffer[m/2-1])
                            i = label >> 24 & 255
                            j = label >> 16 & 255
                            k = label >>  8 & 255
                            l = label       & 255
                        else:
                            label = int(index_buffer[m/2])
                            i = label >> 56 & 255
                            j = label >> 48 & 255
                            k = label >> 40 & 255
                            l = label >> 32 & 255
                ints[self.ijkl(i, j, k, l)] = value_buffer[m-1]
            pos = self.tell()
            self.seek(pos+4)
        return ints

class GamessReader(object):
    '''
    Class for holding method for reading gamess binary files:
        $JOB.F08 : two electron integrals over AO's,
        $JOB.F09 : two electron integrals over MO's,
        $JOB.F10 : the dictionary file with 1-e integrals, orbitals etc.,
        $JOB.F15 : GUGA and ORMAS two-electron reduced density matrix,

    TODO:
        CI coefficients, and CI hamiltonian matrix elements.
    '''
    # add an option to choose which reader should be used for dictionary file
    # and sequential files, the options are:
    # - wrapped fortran code that need to be compiled and installed or
    # - native reader written in python using DictionaryFile and SequentialFile
    # classes
    def __init__(self, log):
        self.logfile    = log
        i = self.logfile.index("log")
        self.filebase   = self.logfile[:i-1]
        self.datfile    = self.filebase + ".dat"
        self.twoeaofile = self.filebase + ".F08"
        self.twoemofile = self.filebase + ".F09"
        self.dictionary = self.filebase + ".F10"
        self.rdm2file   = self.filebase + ".F15"
        self.gp         = GamessLogParser(log=self.logfile)

    def get_onee_size(self, aos=True):
        '''
        Get the size of the vector holding upper (or lower) triangle
        of a square matrix of size naos or nmos.
        '''
        if aos:
            n = self.gp.get_number_of_aos()
        else:
            n = self.gp.get_number_of_mos()
        return n*(n+1)/2

    def get_twoe_size(self):
        '''
        Get the size of the 1d vector holding upper (or lower) triangle
        of a supermatrix of size nmos (2RDM and two-electrons integrals).
        '''
        n = self.get_onee_size(aos=False)
        return n*(n+1)/2

    def read_rdm2(self, filename=None, nmo=None):

        '''Read the 2rdm from the gamess-us file'''

        # initialize numpy array to zeros
        rdm2 = np.zeros(self.get_twoe_size(), dtype=float)


    def read_twoeao(self, filename=None, nmo=None):

        '''Read the two electron integrals from the gamess-us file'''

        ints = np.zeros(self.get_twoe_size(), dtype=float)


from collections import namedtuple

rec = namedtuple('record', ['name', 'dtype'])
records = {
    1 : rec("atomic coordinates", "f8"),
    2 : rec("enrgys", "f8"),
    3 : rec("gradient vector", "f8"),
    4 : rec("hessian matrix", "f8"),
    5 : rec("not used", ""),
    6 : rec("not used", ""),
    7 : rec("ptr", "f8"),
    8 : rec("dtr", "f8"),
    9 : rec("ftr", "f8"),
   10 : rec("gtr", "f8"),
   11 : rec("bare nucleus", "f8"),
   12 : rec("overlap", "f8"),
   13 : rec("kinetic energy", "f8"),
   14 : rec("alpha fock matrix", "f8"),
   15 : rec("alpha orbitals", "f8"),
   16 : rec("alpha density matrix", "f8"),
   17 : rec("alpha energies or occupation numbers", "f8"),
   18 : rec("beta fock matrix", "f8"),
   19 : rec("beta orbitals", "f8"),
   20 : rec("beta density matrix", "f8"),
   21 : rec("beta energies or occupation numbers", "f8"),
   22 : rec("error function interpolation table", "f8"),
   23 : rec("old alpha fock matrix", "f8"),
   24 : rec("older alpha fock matrix", "f8"),
   25 : rec("oldest alpha fock matrix", "f8"),
   26 : rec("old beta fock matrix", "f8"),
   27 : rec("older beta fock matrix", "f8"),
   28 : rec("odest beta fock matrix", "f8"),
   29 : rec("vib 0 gradient in FORCE", "f8"),
   30 : rec("vib 0 alpha orbitals in FORCE", "f8"),
   31 : rec("Vib 0 beta  orbitals in FORCE", "f8"),
   32 : rec("Vib 0 alpha density matrix in FORCE", "f8"),
   33 : rec("Vib 0 beta  density matrix in FORCE", "f8"),
   34 : rec("dipole derivative tensor in FORCE", "f8"),
   35 : rec("frozen core Fock operator", "f8"),
   36 : rec("RHF/UHF/ROHF Lagrangian", "f8"),
   37 : rec("floating point part of common block /OPTGRD/", "f8"),
   38 : rec("integer part of common block /OPTGRD/", "i8"),
   39 : rec("ZMAT of input internal coords", "f8"),
   40 : rec("IZMAT of input internal coords", "i8"),
   41 : rec("B matrix of redundant internal coords", "f8"),
   42 : rec("pristine core Fock matrix in MO basis (see 87)", "f8"),
   43 : rec("Force constant matrix in internal coordinates", "f8"),
   44 : rec("SALC transformation", "f8"),
   45 : rec("symmetry adapted Q matrix", "f8"),
   46 : rec("S matrix for symmetry coordinates", "f8"),
   47 : rec("ZMAT for symmetry internal coords", "f8"),
   48 : rec("IZMAT for symmetry internal coords", "i8"),
   49 : rec("B matrix", "f8"),
   50 : rec("B inverse matrix", "f8"),
   95 : rec("x dipole integrals in AO basis", "f8"),
   96 : rec("y dipole integrals in AO basis", "f8"),
   97 : rec("z dipole integrals in AO basis", "f8"),
  251 : rec("static polarizability tensor alpha", "f8"),
  252 : rec("X dipole integrals in MO basis", "f8"),
  253 : rec("Y dipole integrals in MO basis", "f8"),
  254 : rec("Z dipole integrals in MO basis", "f8"),
  255 : rec("alpha MO symmetry labels", "S8"),
  256 : rec("beta MO symmetry labels", "S8"),
  379 : rec("Lz integrals", "f8"),
}

class DictionaryFile(BinaryFile):
    '''
    Wrapper for reading GAMESS(US) dictionary file (*.F10).
    '''

    def __init__(self, filename, irecln=4090, int_size=8):
        """
        Args:

            irecln: int, is the record length that is used by GAMESS(US) when
                writing the dicitonary file. It is a platform dependent
                variable that is called IRECLN in the GAMESS(US) code. In
                GAMESS(US) it is set by the function NRASIZ(UNIT) in iolib.src
                and for dictionary file (unit=10) and UNX it is equal to 4090,
                for all other files and UNX it is equal to 2048. If you use
                GAMESS(US) on a different platform check the NRASIZ(UNIT)
                function for the proper value and supply it when instantiating
                the class.
            int_size: int, is the integer size (in bytes) that the GAMESS(US)
                was compiled with.
        """
        super(DictionaryFile, self).__init__(filename)

        self.irecln = irecln
        self.int_size = int_size
        # read the first record with the information about the
        # structure of the dictionary file
        self.irecst = self.read(np.dtype('i'+str(self.int_size)))
        self.ioda   = self.read(np.dtype('i'+str(self.int_size)), shape=(950,))
        self.ifilen = self.read(np.dtype('i'+str(self.int_size)), shape=(950,))
        self.iss    = self.read(np.dtype('i'+str(self.int_size)))
        self.ipk    = self.read(np.dtype('i'+str(self.int_size)))

    def read_record(self, nrec, dtype=None):
        '''
        Read a logical record 'rec' from the dictionary file and return a numpy
        array of type defined in the 'records' list, and size defined through
        'self.ifilen' array.
        '''

        if self.ioda[nrec-1] < 0:
            raise IOError("Record {0} was not previously written, IODA[{0}]={1}".format(nrec, self.ioda[nrec-1]))

        self.seek(8*self.irecln*(int(self.ioda[nrec-1])-1))
        if dtype is not None:
            return self.read(dtype, shape=(self.ifilen[nrec-1],))
        else:
            return self.read(records[nrec].dtype, shape=(self.ifilen[nrec-1],))

def tri2full(vector, shape):
    '''
    Convert a triagonal matrix whose elements are stored in the `vector` into a 
    rectangular matrix of the shape given by `shape` tuple.
    '''

    nrow, ncol = shape
    matrix = np.zeros((nrow, ncol), dtype=float, order='F')

    ij = -1
    for i in range(nrow):
        for j in range(i + 1):
            ij += 1
            matrix[i, j] = matrix[j, i] = vector[ij]
    return matrix
