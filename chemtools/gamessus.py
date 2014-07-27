
'''
Module for handling Gamess-US related jobs,:
    Gamess          : running and submitting jobs, writing inputs,
    GamessInpParser : parsing the input file,
    GamessLogParser : parsing the output file,
    GamessDatParser : parsing data from the gamess PUNCH (*.dat) file
    GamessReader    : reading gamess bianry files.
'''

from chemtools.code import Code
from subprocess import Popen
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

    def __init__(self, **kwargs):
        super(Gamess, self).__init__(**kwargs)

        _versions = []

        for item in os.listdir(self.execpath):
            match = re.match(r'[a-z]+\.([0-9]+)\.x', item)
            if match:
                _versions.append(match.group(1))

        if "version" in kwargs.keys():
            if kwargs["version"] in _versions:
                self.version = kwargs["version"]
            else:
                sys.exit('gamess version {0:s} not found in {1:s}'.format(kwargs["version"], self.execpath))
        else:
            self.version = _versions[0]

        if os.path.isfile(os.path.join(self.execpath, 'ddikick.x')):
            self.ddikick = os.path.join(self.execpath, 'ddikick.x')

        if os.path.isfile(os.path.join(self.execpath, 'rungms')):
            self.rungms = os.path.join(self.execpath, 'rungms')
        else:
            sys.exit('Could not find "rungms" under {0:s}'.format(self.execpath))

    def remove_dat(inpfile):
        '''
        Remove the gamess dat file if it exists in the scratch directory.
        '''

        datfile = os.path.splitext(inpfile)[0] + ".dat"
        if remove_dat:
            if os.path.exists(os.path.join(self.scratch, datfile)):
                os.remove(os.path.join(self.scratch, datfile))


    def run(self, inpfile, logfile, nproc=1, remove_dat=True):
        '''
        Run a single gamess job interactively - without submitting to the
        queue.
        '''

        if remove_dat:
            self.remove_dat(inpfile)

        out = open(logfile, 'w')
        process = Popen([self.rungms, inpfile, self.version, str(nproc)], stdout=out, stderr=out)
        process.wait()
        out.close()

    def run_multiple(self, inputs):
        pass

    def submit(self):

        '''Submit a single gamess job into the queue.'''

        pass

    def accomplished(self, outfile=None):
        '''
        Return True if Gamess(US) job finished without errors.
        '''

        if outfile is not None:
            parser = GamessLogParser(outfile)
        else:
            parser = GamessLogParser(self.logfile)
        return parser.accomplished()

    def parse(self):
        pass

    def write_input(self):
        pass

class GamessInpParser(object):
    '''
    A class for parsing and writing gamess-us input files.
    '''

    def __init__(self):
        '''
        Initialize the class.
        '''

        self.end = " $end"
        # not nested groups of input blocks (not parsed into a dict of dicts)
        self._notnested = ["$data", "$vec", "$ecp"]

    def parse(self, inpfile):
        '''
        Parse gamess input file into a dictionary of dictionaries, where the
        highest level entries are gamess namelist fileds and that contain
        dictionaries of options. All key are converted to lowercase.
        '''

        with open(inpfile, 'r') as finp:
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

        dinput = {}

        iterator = pat.finditer(inpstr)
        for match in iterator:
            if match.group("block").lower() not in self._notnested:
                dinput[match.group("block").lower()] = {}
                fields = [s.strip() for s in match.group("entries").split("\n")]
                for field in fields:
                    if not field.startswith("!"):
                        for line in field.split():
                            key, value = line.split("=")
                            dinput[match.group("block").lower()][key.lower()] = value
            elif match.group("block").lower() == "$data":
                dinput["$data"] = self.parse_data(match.group("entries"))
            elif match.group("block").lower() in ["$vec", "$ecp"]:
                dinput[match.group("block").lower()] = match.group("entries")
        return dinput

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


    def write_input(self, inpfile, inpdict):
        '''
        Write a gamess input file under the name <inpfile> based on the
        information fstored in the dictionary <inpdict>.

        Args:
            inpfile (str): name of the file to be written,
            inpdict (dict): dictionary with the gamess input specification.
        '''

        if not isinstance(inpdict, dict):
            raise TypeError("expected a dictionary but got {0:s}".format(type(inpdict)))
        inp = open(inpfile, 'w')

        # write nested namelist groups
        for key, value in sorted(inpdict.items()):
            if key not in  self._notnested:
                inp.write(" {0:<s}\n".format(key))
                for kkey, vvalue in sorted(value.items()):
                    inp.write("    {k:s}={v:s}\n".format(k=kkey, v=str(vvalue)))
                inp.write(self.end + '\n')
        #write $data card
        inp.write(" $data\n")
        inp.write(inpdict['$data'])
        inp.write(self.end + '\n')
        inp.close()

    def set_gamess_input(dinp, mol, bs, code, core):

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
            dinp["$cidet"]["nact"]  = bas.get_no_functions(bs) - core
            dinp["$cidet"]["ncore"] = core
            dinp["$cidet"]["nels"]  = mol.electrons - core*2
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

    def write_with_vec(self, inpfile, inpdict, datfile):
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

        gdp = GamessDatParser(datfile)
        vecstr = gdp.get_nos()

        inpdict['$contrl']['scftyp'] = 'none'
        inpdict['$guess'] = {}
        inpdict['$guess']['guess'] = 'moread'
        inpdict['$guess']['norb'] = str(gdp.get_naos_nmos(vecstr)[1])

        self.write(inpfile, inpdict)

        inp = open(inpfile, "a")
        inp.write("\n $vec\n")
        inp.write(gdp.get_nos())
        inp.write(" $end\n")
        inp.close()

    def print_inpdict(self, inpdict):
        '''
        Neat print of the dictionary with parsed input
        '''

        for key, value in sorted(inpdict.items()):
            if isinstance(value, dict):
                print key
                for kkey, vvalue in sorted(value.items()):
                    print "\t{0:<10s} : {1:}".format(kkey, vvalue)
            else:
                print key, '\n', value

    def parse_gamess_basis(self, basisset):
        '''
        Parse the basis set into a list of dictionaries from a string in
        gamess format.
        '''

        pat = re.compile(r'^\s*(?P<shell>[SPDFGHIspdfghi])\s*(?P<nf>[1-9]+)')

        bs = []
        bslines = rawbs.split("\n")

        for i, line in enumerate(bslines):
            match = pat.search(line)
            if match:
                bs.extend(self.parse_function(match.group("shell"),
                                        atomic,
                                        bslines[i+1:i+int(match.group("nf"))+1]))
        return bs

    def parse_function(self, shell, atomic, los):

        '''Parse basis info from list of strings'''

        exps   = [float(item.split()[1]) for item in los]
        coeffs = [float(item.split()[2]) for item in los]

        bs =[{"shell"  : _shells.index(shell.upper()),
            "typ"    : "",
            "exps"   : exps,
            "coeffs" : coeffs,
            "atomic" : atomic,
            }]
        return bs

class GamessLogParser(object):

    '''methods for parsing gamess-us log file'''

    def __init__(self, log):
        self.logfile = log

    @property
    def logfile(self):
        return self._logfile

    @logfile.setter
    def logfile(self, value):
        if os.path.exists(value):
            self._logfile = value
        else:
            raise ValueError("File: {} does not exist".format(value))

    def logexists(self):

        '''Check if the log file exists.'''

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

        regex = r'NUMBER OF LINEARLY DEPENDENT MOS DROPPED=\s*(\d+)'
        match = self.parse(regex)
        if match:
            return match.group(1)

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

    @staticmethod
    def parse_pairs(los, sep="="):
        '''
        Parse a given list of strings "los" into a dictionary based on separation
        by "sep" character and return the dictionary.
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

class GamessReader(object):

    '''Class for holding method for reading gamess binary files:
        $JOB.F08 : two electron integrals over AO's,
        $JOB.F09 : two electron integrals over MO's,
        $JOB.F10 : the dictionary file with one electron integrals, orbitals etc.,
        $JOB.F15 : GUGA and ORMAS two-electron reduced density matrix,

        TODO:
        CI coefficients, and CI hamiltonian matrix elements.'''

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

    def factor(i,j,k,l):
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

    def ijkl(self, i,j,k,l):
        '''
        Based on the four orbital indices i,j,k,l return the address
        in the 1d vector.
        '''
        ij = max(i, j)*(max(i, j) + 1)/2 + min(i, j)
        kl = max(k, l)*(max(k, l) + 1)/2 + min(k, l)
        return max(ij, kl)*(max(ij, kl) + 1)/2 + min(ij, kl)

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
                                print "{0:3d}{1:3d}{2:3d}{3:3d} {4:25.14f}".format(
                                    i, j, k, l, twoe[self.ijkl(i,j,k,l)])

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
        gip = GamessInpParser()
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

        naos, nmos, nlines = self.get_naos_nmos(vecstr)
        orblines = vecstr.split('\n')
        orbs = np.zeros((naos, nmos), dtype=float)
        counter = -1
        for i in range(0, nmos):
            for j in range(0, nlines):
                counter += 1
                nitems = int(len(orblines[counter][5:]))/clength
                orbs[5*j:5*(j+1), i] = [float(orblines[counter][5+15*n:5+15*(n+1)]) for n in range(nitems)]
        return orbs

    def get_naos_nmos(self, vecstr, clength=15):
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

    @staticmethod
    def find_between(s, first, last):
        try:
            start = s.index(first) + len(first)
            end = s.index(last, start)
            return s[start:end]
        except ValueError:
            return ""

def take(seq, num):
    '''
    Iterate over a sequence "seq" "num" times and return the list of the
    elements iterated over.
    '''
    return [next(seq) for i in range(num)]


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
            raise ValueError, "order should be either 'fortran' or 'c'."
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
            raise ValueError, "shape must be a tuple"
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
            raise EOFError, "Asking for more data than available in file."
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
        super(SequrntialFile, self).__init__(filename)
        log = os.path.splitext(filename)[0] + ".log"
        glp = GamessLogParser(log=log)
        self.nmo = glp.get_number_of_mos()
        if self.nmo < 255:
            self.large_labels = 1
        else:
            self.large_labels = 2


    def readseq(self, buff_size=15000, int_size=8, mos=False):
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
                size of the buffer holding values to be read, in gamess(us) it is
                stored under "NINTMX" variable and in Linux version is equal to
                15000 which is the default value,
            large_labels (bool)
                a flag indicating if large labels should were used, if largest
                label "i" (of the MO) is i<255 then large_labels should be False
                (case "LABSIZ=1" in gamess(us)), otherwise set to True (case
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
            self.seek(4)
            hmo = self.read('f8', shape=(self.nmo*(self.nmo+1)/2))
            pos = self.tell()
            self.seek(pos+4)

        int_type = np.dtype('i'+str(int_size))
        index_buffer = np.zeros(indexBuffSize, dtype=int_type, order='F')
        value_buffer = np.zeros(buff_size, dtype=float, order='F')

        length = 1
        while length > 0:
            pos = self.tell()
            print pos
            self.seek(pos+4)
            length = self.read(int_type)
            if length > buff_size:
                raise ValueError('the read record length: {0:10d} greater that the buffer size {1:10d}'.format(length, buff_size))
            print length
            index_buffer = self.read(int_type, shape=(indexBuffSize))
            value_buffer = self.read('f8', shape=(buff_size))

            pos = self.tell()
            self.seek(pos+4)

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
                #print(i,j,k,l)

        #print values
        return hmo

from collections import namedtuple

rec = namedtuple('record', ['name', 'dtype'])
records = [
    rec("atomic coordinates", "f8"),
    rec("enrgys", "f8"),
    rec("gradient vector", "f8"),
    rec("hessian matrix", "f8"),
    rec("not used", ""),
    rec("not used", ""),
    rec("ptr", "f8"),
    rec("dtr", "f8"),
    rec("ftr", "f8"),
    rec("gtr", "f8"),
    rec("bare nucleus", "f8"),
    rec("overlap", "f8"),
    rec("kinetic energy", "f8"),
    rec("alpha fock matrix", "f8"),
    rec("alpha orbitals", "f8"),
    rec("alpha density matrix", "f8"),
    rec("alpha energies or occupation numbers", "f8"),
    rec("beta fock matrix", "f8"),
    rec("beta orbitals", "f8"),
    rec("beta density matrix", "f8"),
    rec("beta energies or occupation numbers", "f8"),
    rec("error function interpolation table", "f8"),
    rec("old alpha fock matrix", "f8"),
    rec("older alpha fock matrix", "f8"),
    rec("oldest alpha fock matrix", "f8"),
    rec("old beta fock matrix", "f8"),
    rec("older beta fock matrix", "f8"),
    rec("odest beta fock matrix", "f8"),
    rec("vib 0 gradient in FORCE", "f8"),
    rec("vib 0 alpha orbitals in FORCE", "f8"),
    rec("Vib 0 beta  orbitals in FORCE", "f8"),
    rec("Vib 0 alpha density matrix in FORCE", "f8"),
    rec("Vib 0 beta  density matrix in FORCE", "f8"),
    rec("dipole derivative tensor in FORCE", "f8"),
    rec("frozen core Fock operator", "f8"),
    rec("RHF/UHF/ROHF Lagrangian", "f8"),
    rec("floating point part of common block /OPTGRD/", "f8"),
    rec("integer part of common block /OPTGRD/", "i8"),
    rec("ZMAT of input internal coords", "f8"),
]


    #int 40. IZMAT of input internal coords
    #    41. B matrix of redundant internal coords
    #    42. pristine core Fock matrix in MO basis (see 87)
    #    43. Force constant matrix in internal coordinates.
    #    44. SALC transformation
    #    45. symmetry adapted Q matrix
    #    46. S matrix for symmetry coordinates
    #    47. ZMAT for symmetry internal coords
    #int 48. IZMAT for symmetry internal coords
    #    49. B matrix
    #    50. B inverse matrix
    #    51. overlap matrix in Lowdin basis,
    #        temp Fock matrix storage for ROHF
    #    52. genuine MOPAC overlap matrix
    #    53. MOPAC repulsion integrals
    #    54. exchange integrals for screening
    #    55. orbital gradient during SOSCF MCSCF
    #    56. orbital displacement during SOSCF MCSCF
    #    57. orbital hessian during SOSCF MCSCF
    #    58. reserved for Pradipta
    #    59. Coulomb integrals in Ruedenberg localizations
    #    60. exchange integrals in Ruedenberg localizations
    #    61. temp MO storage for GVB and ROHF-MP2
    #    62. temp density for GVB
    #    63. dS/dx matrix for hessians
    #    64. dS/dy matrix for hessians
    #    65. dS/dz matrix for hessians
    #    66. derivative hamiltonian for OS-TCSCF hessians
    #    67. partially formed EG and EH for hessians
    #    68. MCSCF first order density in MO basis
    #    69. alpha Lowdin populations
    #    70. beta Lowdin populations
    #    71. alpha orbitals during localization
    #    72. beta orbitals during localization
    #    73. alpha localization transformation
    #    74. beta localization transformation
    #    75. fitted EFP interfragment repulsion values
    #    76. model core potential information
    #    77. model core potential information
    #    78. "Erep derivative" matrix associated with F-a terms
    #    79. "Erep derivative" matrix associated with S-a terms
    #    80. EFP 1-e Fock matrix including induced dipole terms
    #    81. interfragment dispersion values
    #    82. MO-based Fock matrix without any EFP contributions
    #    83. LMO centroids of charge
    #    84. d/dx dipole velocity integrals
    #    85. d/dy dipole velocity integrals
    #    86. d/dz dipole velocity integrals
    #    87. unmodified h matrix during SCRF or EFP, AO basis
    #    88. PCM solvent operator contribution to Fock
    #    89. EFP multipole contribution to one e- Fock matrix
    #    90. ECP coefficients
    #int 91. ECP labels
    #    92. ECP coefficients
    #int 93. ECP labels
    #    94. bare nucleus Hamiltonian during FFIELD runs
    #    95. x dipole integrals, in AO basis
    #    96. y dipole integrals, in AO basis
    #    97. z dipole integrals, in AO basis
    #    98. former coords for Schlegel geometry search
    #    99. former gradients for Schlegel geometry search
    #   100. dispersion contribution to EFP gradient
    #
    #     records 101-248 are used for NLO properties
    #
    #101. U'x(0)       149. U''xx(-2w;w,w)   200. UM''xx(-w;w,0)
    #102.   y          150.    xy            201.    xy
    #103.   z          151.    xz            202.    xz
    #104. G'x(0)       152.    yy            203.    yz
    #105.   y          153.    yz            204.    yy
    #106.   z          154.    zz            205.    yz
    #107. U'x(w)       155. G''xx(-2w;w,w)   206.    zx
    #108.   y          156.    xy            207.    zy
    #109.   z          157.    xz            208.    zz
    #110. G'x(w)       158.    yy            209. U''xx(0;w,-w)
    #111.   y          159.    yz            210.    xy
    #112.   z          160.    zz            211.    xz
    #113. U'x(2w)      161. e''xx(-2w;w,w)   212.    yz
    #114.   y          162.    xy            213.    yy
    #115.   z          163.    xz            214.    yz
    #116. G'x(2w)      164.    yy            215.    zx
    #117.   y          165.    yz            216.    zy
    #118.   z          166.    zz            217.    zz
    #119. U'x(3w)      167. UM''xx(-2w;w,w)  218. G''xx(0;w,-w)
    #120.   y          168.     xy           219.    xy
    #121.   z          169.     xz           220.    xz
    #122. G'x(3w)      170.     yy           221.    yz
    #123.   y          171.     yz           222.    yy
    #124.   z          172.     zz           223.    yz
    #125. U''xx(0)     173. U''xx(-w;w,0)    224.    zx
    #126.    xy        174.    xy            225.    zy
    #127.    xz        175.    xz            226.    zz
    #128.    yy        176.    yz            227. e''xx(0;w,-w)
    #129.    yz        177.    yy            228.    xy
    #130.    zz        178.    yz            229.    xz
    #131. G''xx(0)     179.    zx            230.    yz
    #132.    xy        180.    zy            231.    yy
    #133.    xz        181.    zz            232.    yz
    #134.    yy        182. G''xx(-w;w,0)    233.    zx
    #135.    yz        183.    xy            234.    zy
    #136.    zz        184.    xz            235.    zz
    #137. e''xx(0)     185.    yz            236. UM''xx(0;w,-w) 
    #138.    xy        186.    yy            237.     xy
    #139.    xz        187.    yz            238.     xz
    #140.    yy        188.    zx            239.     yz
    #141.    yz        189.    zy            240.     yy
    #142.    zz        190.    zz            241.     yz
    #143. UM''xx(0)    191. e''xx(-w;w,0)    242.     zx
    #144.     xy       192.    xy            243.     zy
    #145.     xz       193.    xz            244.     zz
    #146.     yy       194.    yz
    #147.     yz       195.    yy
    #148.     zz       196.    yz
    #                  197.    zx
    #                  198.    zy
    #                  199.    zz
    #
    #    245. old NLO Fock matrix
    #    246. older NLO Fock matrix
    #    247. oldest NLO Fock matrix
    #    249. polarizability derivative tensor for Raman
    #    250. transition density matrix in AO basis
    #    251. static polarizability tensor alpha
    #    252. X dipole integrals in MO basis
    #    253. Y dipole integrals in MO basis
    #    254. Z dipole integrals in MO basis
    #    255. alpha MO symmetry labels
    #    256. beta MO symmetry labels
    #    257. not used
    #    258. Vnn gradient during MCSCF hessian
    #    259. core Hamiltonian from der.ints in MCSCF hessian
    #260-261. reserved for Dan
    #    262. MO symmetry integers during determinant CI
    #    263. PCM nuclei/induced nuclear Charge operator
    #    264. PCM electron/induced nuclear Charge operator
    #    265. pristine alpha guess (MOREAD or Huckel+INSORB)
    #    266. EFP/PCM IFR sphere information
    #    267. fragment LMO expansions, for EFP Pauli
    #    268. fragment Fock operators, for EFP Pauli
    #    269. fragment CMO expansions, for EFP charge transfer
    #    270. not used
    #    271. orbital density matrix in divide and conquer
    #int 272. subsystem data during divide and conquer
    #    273. old alpha Fock matrix for D&C Anderson-like DIIS
    #    274. old  beta Fock matrix for D&C Anderson-like DIIS
    #    275. not used
    #    276. Vib 0 Q matrix    in FORCE
    #    277. Vib 0 h integrals in FORCE
    #    278. Vib 0 S integrals in FORCE
    #    279. Vib 0 T integrals in FORCE
    #    280. Zero field LMOs during numerical polarizability
    #    281. Alpha zero field dens. during num. polarizability
    #    282. Beta zero field dens. during num. polarizability
    #    283. zero field Fock matrix. during num. polarizability
    #    284. Fock eigenvalues for multireference PT
    #    285. density matrix or Fock matrix over LMOs
    #    286. oriented localized molecular orbitals
    #    287. density matrix of oriented LMOs
    #290-299. not used
    #    301. Pocc during MP2 (RHF or ZAPT) or CIS grad
    #    302. Pvir during MP2 gradient (UMP2= 411-429)
    #    303. Wai during MP2 gradient
    #    304. Lagrangian Lai during MP2 gradient
    #    305. Wocc during MP2 gradient
    #    306. Wvir during MP2 gradient
    #    307. P(MP2/CIS)-P(RHF) during MP2 or CIS gradient
    #    308. SCF density during MP2 or CIS gradient
    #    309. energy weighted density in MP2 or CIS gradient
    #    311. Supermolecule h during Morokuma
    #    312. Supermolecule S during Morokuma
    #    313. Monomer 1 orbitals during Morokuma
    #    314. Monomer 2 orbitals during Morokuma
    #    315. combined monomer orbitals during Morokuma
    #    316. RHF density in CI grad; nonorthog. MOs in SCF-MI
    #    317. unzeroed Fock matrix when MOs are frozen
    #    318. MOREAD orbitals when MOs are frozen
    #    319. bare Hamiltonian without EFP contribution
    #    320. MCSCF active orbital density
    #    321. MCSCF DIIS error matrix
    #    322. MCSCF orbital rotation indices
    #    323. Hamiltonian matrix during QUAD MCSCF
    #    324. MO symmetry labels during MCSCF
    #    325. final uncanonicalized MCSCF orbitals
    #326-329. not used
    #    330. CEL matrix during PCM
    #    331. VEF matrix during PCM
    #    332. QEFF matrix during PCM
    #    333. ELD matrix during PCM
    #    334. PVE tesselation info during PCM
    #335-339. not used
    #    340. DFT alpha Fock matrix
    #    341. DFT beta Fock matrix
    #    342. DFT screening integrals
    #    343. DFT: V aux basis only
    #    344. DFT density gradient d/dx integrals
    #    345. DFT density gradient d/dy integrals
    #    346. DFT density gradient d/dz integrals
    #    347. DFT M[D] alpha density resolution in aux basis
    #    348. DFT M[D] beta density resolution in aux basis
    #    349. DFT orbital description
    #    350. overlap of true and auxiliary DFT basis
    #    351. previous iteration DFT alpha density
    #    352. previous iteration DFT beta density
    #    353. DFT screening matrix (true and aux basis)
    #    354. DFT screening integrals (aux basis only)
    #    355. h in MO basis during DDI integral transformation
    #    356. alpha symmetry MO irrep numbers if UHF/ROHF
    #    357. beta  symmetry MO irrep numbers if UHF/ROHF
    #358-369. not used
    #    370. left transformation for pVp
    #    371. right transformation for pVp
    #    370. basis A (large component) during NESC
    #    371. basis B (small component) during NESC
    #    372. difference basis set A-B1 during NESC
    #    373. basis N (rel. normalized large component)
    #    374. basis B1 (small component) during NESC
    #    375. charges of non-relativistic atoms in NESC
    #    376. common nuclear charges for all NESC basis
    #    377. common coordinates for all NESC basis
    #    378. common exponent values for all NESC basis
    #    372. left transformation for V  during RESC
    #    373. right transformation for V during RESC
    #    374. 2T, T is kinetic energy integrals during RESC
    #    375. pVp integrals during RESC
    #    376. V integrals during RESC
    #    377. Sd, overlap eigenvalues during RESC
    #    378. V, overlap eigenvectors during RESC
    #    379. Lz integrals
    #    380. reserved for Ly integrals.
    #    381. reserved for Lx integrals.
    #    382. X, AO orthogonalisation matrix during RESC
    #    383. Td, eigenvalues of 2T during RESC
    #    384. U, eigenvectors of kinetic energy during RESC
    #    385. exponents and contraction for the original basis
    #int 386. shell integer arrays for the original basis
    #    387. exponents and contraction for uncontracted basis
    #int 388. shell integer arrays for the uncontracted basis
    #    389. Transformation to contracted basis
    #    390. S integrals in the internally uncontracted basis
    #    391. charges of non-relativistic atoms in RESC
    #    392. copy of one e- integrals in MO basis in SO-MCQDPT
    #    393. Density average over all $MCQD groups in SO-MCQDPT
    #    394. overlap integrals in 128 bit precision
    #    395. kinetic ints in 128 bit precision, for relativity
    #396-400. not used
    #    401. dynamic polarizability tensors
    #    402. GVB Lagrangian
    #    403. MCSCF Lagrangian
    #    404. GUGA CI Lagrangian (see 308 for CIS)
    #    405. not used
    #    406. MEX search state 1 alpha orbitals
    #    407. MEX search state 1 beta orbitals
    #    408. MEX search state 2 alpha orbitals
    #    409. MEX search state 2 beta orbitals
    #    410. not used
    #    411. alpha Pocc during UMP2 gradient (see 301-309)
    #    412. alpha Pvir during UMP2 gradient
    #    413. alpha Wai during UMP2 gradient
    #    414. alpha Lagrangian Lai during UMP2 gradient
    #    415. alpha Wocc during UMP2 gradient
    #    416. alpha Wvir during UMP2 gradient
    #    417. alpha P(MP2/CIS)-P(RHF) during UMP2/USFTDDFT grad
    #    418. alpha SCF density during UMP2/USFTDDFT gradient
    #    419. alpha energy wghted density in UMP2/USFTDDFT grad
    #    420. not used
    #421-429. same as 411-419, for beta orbitals
    #    430. not used
    #440-469. reserved for NEO
    #    470. QUAMBO expansion matrix
    #    471. excitation vectors for FMO-TDDFT
    #    472. X+Y in MO basis during TD-DFT gradient
    #    473. X-Y in MO basis during TD-DFT gradient
    #    474. X+Y in AO basis during TD-DFT gradient
    #    475. X-Y in AO basis during TD-DFT gradient
    #    476. excited state density during TD-DFT gradient
    #    477. energy-weighted density in AO basis for TD-DFT
    #478-489. not used
    #    490. transition Lagrangian right hand side during NACME
    #    491. gradients vectors during NACME
    #    492. NACME vectors during NACME
    #    493. difference gradient in conical intersection search
    #    494. derivative coupling vector in CI search
    #    495. mean energy gradient in CI search
    #    496. unused
    #    497. temp storage of gradient of 1st state in CI search
    #498-500. not used
    #    501. A2 cavity data in COSMO
    #    502. A3 cavity data in COSMO
    #    503. AMTSAV cavity data in COSMO
    #504-510. not used
    #    511. effective polarizability in LRD
    #    512. C6 coefficients in LRD
    #    513. C8 coefficients in LRD
    #    514. C10 coefficients in LRD
    #    515. atomic pair LRD energy
    #    520. Malmqvist factorized orb transformation (wrt 325)
    #    521. SVD localized orthogonal orbitals
    #    522. SVD localized nonorthogonal orbitals
    #    523. initial-to-SVD LMO nonorthogonal transformatio
    #    524. SVD LMO nonorthogonal-to-orthogonal transformation
    #    525. initial-to-SVD LMO orthog transformation (wrt 15)
    #    526. 1st order density for orthogonal SVD localized MOs
    #    527. collective orbital reordering for Malmqvist
    #    528. atom-to-orbital assignment for SVD orbitals
    #    529. Malmqvist re-ordered set of SVD LMOs
    #    530. oriented SVD density in the order of record 527
    #    531. oriented or SVD atom-to-orbital assignment for CT
    #    532. block zapped 'standard Fock operator' in AO basis
    #    533. overlap of stored atom's MBS with current basis
    #534-540. not used
    #    541. pristine MCSCF orbs during diabatization
    #    542. reference geometry orbs during diabatization
    #    543. PT2 state rotation during diabatization
    #    544. PT2 state energies during diabatization
    #    545. PT2's CAS-CI largest CI coefs, in diabatization
    #    532-950. not used

class DictionaryFile(BinaryFile):
    '''
    Handler for teh gamess(us) dictionary file.
    '''

    def __init__(self, filename, irecl=4090):
        self.irecl = irecl
        super(DictionaryFile, self).__init__(filename)

        # read the first record with the information about the
        # structure of the dictionary file
        self.irecst = self.read(np.dtype('i8'))
        self.ioda   = self.read(np.dtype('i8'), shape=(950,))
        self.ifilen = self.read(np.dtype('i8'), shape=(950,))
        self.iss    = self.read(np.dtype('i8'))
        self.ipk    = self.read(np.dtype('i8'))

    def read_record(self, nrec, dtype=None):
        '''
        Read a logical record 'rec' from the dictionary file and return a numpy
        array of type defined in the 'records' list, and size defined through
        'self.ifilen' array.
        '''

        if self.ioda[nrec-1] < 0:
            raise ValueError("Record {0} was not previously written, IODA[{0}]={1}".format(nrec, self.ioda[nrec-1]))

        self.seek(8*4090*(int(self.ioda[nrec-1])-1))
        if dtype is not None:
            return self.read(dtype, shape=(self.ifilen[nrec-1],))
        else:
            return self.read(records[nrec-1].dtype, shape=(self.ifilen[nrec-1],))
