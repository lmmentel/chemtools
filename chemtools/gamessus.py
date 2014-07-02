
'''
Module for handling Gamess-US related jobs,:
    Gamess          : running and submitting jobs, writing inputs,
    GamessInpParser : parsing the input file,
    GamessLogParser : parsing the output file,
    GamessDatParser : parsing data from the gamess PUNCH (*.dat) file
    GamessReader    : reading gamess bianry files.
'''

from code import Code
from molecule import Molecule
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

    def run(self, inpfile, nproc=1, remove_dat=True):
        '''
        Run a single gamess job interactively - without submitting to the
        queue.
        '''

        basfile = os.path.splitext(inpfile)[0]
        datfile = basfile + ".dat"
        logfile = basfile + ".log"
        if remove_dat:
            if os.path.exists(os.path.join(self.scratch, datfile)):
                os.remove(os.path.join(self.scratch, datfile))

        out = open(logfile, 'w')
        process = Popen([self.rungms, inpfile, self.version, str(nproc)], stdout=out, stderr=out)
        process.wait()
        out.close()

    def run_multiple(self, inputs):
        pass

    def submit(self):

        '''Submit a single gamess job into the queue.'''

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
        # not nested group of input (not parsed into a dict of dicts)
        self._notnested = ["$data", "$vec"]

    def parse_inp(inpfile):
        '''parse the input file'''

        with open(inpfile, 'r') as finp:
            contents = finp.read()
        return self.parse(contents)

    def parse(self, inpstr):
        '''
        Parse gamess input file into a dictionary of dictionaries, where the
        highest level entries are gamess namelist fileds and that contain
        dictionaries of options. All key are converted to lowercase. For
        example if the following input was parsed:

        >>> from gamessus import GamessInpParser as GIP
        >>> gip = GIP()
        >>> gip.parse(""" $CONTRL scftyp=rhf units=bohr\
            runtyp=energy   cityp=ormas $END\
            $SYSTEM TIMLIM=1000 mwords=500 $END""")
        {"$contrl" : {"scftyp" : "rhf",
                    "units"  : "bohr",
                    "runtyp" : "energy",
                    "cityp"  : "ormas"},
        "$system" : {"timlim" : "1000",
                    "mwords" : "500"},}
        '''

        pat = re.compile(r'(?P<block>\$[a-zA-Z]{3,6})\s+(?P<entries>.*?)\$END', flags=re.S)

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
                dinput["$data"] = match.group("entries")
        return dinput

    def parse_data(self, inpdict):
        '''
        Parse $DATA block specified in the gamess input file. The parser
        assumes that in $DATA input the atom positions are specified using
        cartesian coordiantes and returns a list of dictionaries with parsed
        information about each specified atom.

        Args:
            inpdict (dict)
                dictionary with parsed input contents

        Returns:
            atoms (list of dicts)
                list of dictionaries holding parsed data about each atom
                specified in the $data block
            title (str)
                title as specified in the $data block
            group (str)
                point group symmetry as specified in the $data group
        '''

        block = re.compile(r'(?P<label>[a-zA-Z]{1,2}[0-9]{0,2})\s*'
                  +r'(?P<atomic>\d+\.\d+)'
                  +'(?P<xyz>(\s+\-?\d+\.\d+){3})\s*'
                  +'(?P<basis>.*?)\n\s*\n', flags=re.S)

        title = inpdict['$data'].split('\n')[0]
        group = inpdict['$data'].split('\n')[1]

        atoms = []
        itfound = block.finditer(inpdict['$data'])
        for m in itfound:
            atoms.append({'label'  : m.group('label'),
             'atomic' : m.group('atomic'),
             'xyz'    : m.group('xyz'),
             'basis'  : m.group('basis'),})
        return atoms, title, group


    def write(self, inpfile, inpdict):
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

class GamessLogParser(object):

    '''methods for parsing gamess-us log file'''

    def __init__(self, log):
        self.logfile = log
        self.logexists()

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

    def terminatedOK(self):

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

        regex = r'.*CCTYP =(?P<cctyp>[A-Z]+)'
        match = self.parse(regex)
        if match:
            return match.group("cctyp")

    def get_ci_type(self):

        '''Get the information on CITYP used in the gamess job.'''

        regex = r'.*CITYP =(?P<cityp>[A-Z]+)'
        match = self.parse(regex)
        if match:
            return match.group("cityp")

    def get_hf_total_energy(self):

        '''Return the total HF energy.'''

        with open(self.logfile, 'r') as log:
            data = log.read()

        regex = r'FINAL R[O]*HF ENERGY IS\s+(?P<energy>\-?\d+\.\d+)'
        match = self.parse(regex)
        if match:
            return float(match.group("energy"))

    def get_ormas_total_energy(self):

        '''Return the total ORMAS CI energy.'''

        with open(self.logfile, 'r') as log:
            data = log.read()

        orcire = re.compile(r'ORMAS CI PROPERTIES.*TOTAL ENERGY \=\s*(?P<energy>\-?\d+\.\d+).*END OF', flags=re.S)
        match = orcire.search(data)
        if match:
            return float(match.group("energy"))

    def get_ci_total_energy(self):

        '''Return the total CI energy.'''

        with open(self.logfile, 'r') as log:
            data = log.read()
        energies = re.findall(r'^\s+TOTAL ENERGY =\s*(\-?\d+\.\d+)', data,
                              re.MULTILINE)
        return float(energies[-1])

    def get_ccsd_total_energy(self):

        '''Return total CCSD energy.'''

        with open(self.logfile, 'r') as log:
            data = log.read()
        ccsdt_re = r'CCSD\s{1,4}ENERGY:\s*(\-?\d+\.\d+)'
        compre = re.compile(ccsdt_re)
        match = compre.search(data, re.MULTILINE)
        if match:
            return float(match.group(1))

    def get_ccsdt_total_energy(self):

        '''Return total CCSD(T) energy.'''

        with open(self.logfile, 'r') as log:
            data = log.read()
        ccsdt_re = r'CCSD\(T\) ENERGY:\s*(\-?\d+\.\d+)'
        compre = re.compile(ccsdt_re)
        match = compre.search(data, re.MULTILINE)
        if match:
            return float(match.group(1))

    def get_ci_ee_energy(self):
        with open(self.logfile, 'r') as log:
            data = log.read()
        energies = re.findall(r'^\s+TWO ELECTRON ENERGY =\s*(\-?\d+\.\d+)', data, re.MULTILINE)
        return float(energies[-1])

    def get_ci_oe_energy(self):
        with open(self.logfile, 'r') as log:
            data = log.read()
        energies = re.findall(r'^\s+ONE ELECTRON ENERGY =\s*(\-?\d+\.\d+)', data, re.MULTILINE)
        return float(energies[-1])

    def get_ci_nucrep_energy(self):
        with open(self.logfile, 'r') as log:
            data = log.read()
        energies = re.findall(r'^\s+NUCLEAR REPULSION ENERGY =\s*(\-?\d+\.\d+)', data, re.MULTILINE)
        return float(energies[-1])

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

        return parse_pairs(slice_after(data, header, 22))

class GamessReader(object):

    '''Class for holding method for reading gamess binary files:
        $JOB.F08 : two electron integrals over AO's,
        $JOB.F09 : two electron integrals over MO's,
        $JOB.F10 : the dictionary file with one electron integrals, orbitals etc.,
        $JOB.F15 : GUGA and ORMAS two-electron reduced density matrix,

        TODO:
        CI coefficients, and CI hamiltonian amtrix elements.'''

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

    def get_occupations(self):
        '''
        Parse the occupation numbers from the ascii PUNCH file (*.dat).
        '''

        with open(self.datfile, 'r') as dat:
            data = dat.read()

        no_patt = re.compile(r'\$OCC(.*?)\$END', flags=re.DOTALL)
        match = no_patt.search(data)
        nooc = []
        if match:
            for line in match.group(1).split('\n'):
                nooc.extend([float(x) for x in line.split()])
            return np.asarray(nooc)
        else:
            sys.exit('No section with occupation numbers found.')

    def get_nos(self):
        '''
        Parse the natural orbitals from the ascii PUNCH file (*.dat).
        '''

        with open(self.datfile, 'r') as dat:
            data = dat.read()

        no_patt = re.compile(r'NO.*\$VEC\s+\n(.*?)\$END', flags=re.DOTALL)
        match = no_patt.search(data)
        if match:
            return match.group(1)
        else:
            sys.exit('No section with natural orbitals found.')

    def parse_nos(self):
        '''
        Parse the orbitals read from the $JOB.dat file in the ASCII format into
        an 2 dimensional array.
        '''

        return self.parse_orbitals(self.get_nos())

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

        naos, nmosi, nlines = self.get_naos_nmos(vecstr)

        orblines = orbstr.split('\n')[:-1]

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
        Get the number of AO's and MO's from the $VEC string.

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
                number of lines permolecular orbital
        '''

        veclines = len(vecstr.split('\n')) - 1

        vecit = iter(vecstr.split('\n'))
        nlines = 0
        while vecit.next()[:2] == ' 1':
            nlines += 1

        naos = 5*(nlines - 1) + len(vecstr.split('\n')[nlines-1][5:])/clength
        nmos = veclines/nlines
        return naos, nmos, nlines


def take(seq, num):
    '''
    Iterate over a sequence "seq" "num" times and return the list of the
    elements iterated over.
    '''
    return [next(seq) for i in range(num)]

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

def slice_after(seq, item, num):
    '''
    Return "num" elements of a sequence "seq" present after the item "item".
    '''
    it = iter(seq)
    for element in it:
        if item in element:
            return [next(it) for i in range(num)]

