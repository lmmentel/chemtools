'''
Module for handling Gamess-US related jobs,:
    Gamess          : running and submitting jobs, writing inputs,
    GamessInpParser : parsing the input file,
    GamessLogParser : parsing the output file,
    GamessDatParser : parsing data from the gamess PUNCH (*.dat) file
'''

from __future__ import print_function

from chemtools.code import Code
from subprocess import Popen
from collections import OrderedDict
import numpy as np
import os
import re
import sys


class Gamess(Code):

    '''Container object for Gamess-us jobs.'''

    def __init__(self, name="GamessUS", version="00", **kwargs):
        super(Gamess, self).__init__(**kwargs)
        self.name = name
        self.rungms = self.executable
        self.gmspath = os.path.dirname(self.executable)
        self.version = version

        if os.path.isfile(os.path.join(self.gmspath, 'ddikick.x')):
            self.ddikick = os.path.join(self.gmspath, 'ddikick.x')

    @property
    def version(self):
        """Return the version."""
        return self._version

    @version.setter
    def version(self, value):
        """Check if the version exist and set it if it does."""

        versions = []
        for item in os.listdir(self.gmspath):
            match = re.match(r'[a-z]+\.([0-9]+)\.x', item)
            if match:
                versions.append(match.group(1))

        if value in versions:
            self._version = value
        else:
            raise IOError('GamessUS version {0:s} not found in {1:s}'.format(value, self.gmspath))

    def remove_dat(self, inpfile):
        '''
        Remove the gamess dat file if it exists in the scratch directory.
        '''

        datfile = os.path.splitext(inpfile)[0] + ".dat"
        if os.path.exists(os.path.join(self.scratch, datfile)):
            os.remove(os.path.join(self.scratch, datfile))


    def run(self, inpfile, logfile=None):
        '''
        Run a single gamess job interactively - without submitting to the
        queue.
        '''

        if self.runopts["remove_dat"]:
            self.remove_dat(inpfile)

        if logfile is None:
            logfile = os.path.splitext(inpfile)[0] + ".log"

        out = open(logfile, 'w')
        process = Popen([self.executable, inpfile, self.version, str(self.runopts["nproc"])], stdout=out, stderr=out)
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

    def parse(self, output, method, objective, regexp=None):
        '''
        Parser molpro output file to get the objective.
        '''

        parser = GamessLogParser(output)

        if objective == "total energy":
            if method == "hf":
                return parser.get_hf_total_energy()
            elif method == "cisd":
                energies = self.get_energy_components(method)
                return energies["TOTAL ENERGY"]
        elif objective == "correlation energy":
                energies = self.get_energy_components(method)
                return energies["TOTAL ENERGY"] - parser.get_hf_total_energy()
        elif objective == "core energy":
            if method == "cisd":
                energies = self.get_energy_components(method)
                return energies["TOTAL ENERGY"]
        elif objective == "regexp":
            return parser.get_variable(regexp)
        else:
            raise ValueError("unknown objective in prase {0:s}".format(objective))

    def write_input(self, fname=None, template=None, core=None, bs=None, mol=None):
        '''
        Write a file containing gamess input.
        '''

        gi = GamessInput(fname=fname, template=template)
        gi.write_input(mol=mol, bs=bs, core=core)

    def __repr__(self):
        return "\n".join(["<Gamess(",
                        "\tname={},".format(self.name),
                        "\tgmspath={},".format(self.gmspath),
                        "\trungms={},".format(self.rungms),
                        "\tversion={},".format(self.version),
                        "\tscratch={},".format(self.scratch),
                        "\trunopts={},".format(str(self.runopts)),
                        ")>"])

class GamessInput(object):
    '''
    A class for parsing and writing gamess-us input files.
    '''

    def __init__(self, fname=None, template=None):
        '''
        Initialize the class.
        '''

        self.fname = fname
        self.template = template
        self.end = " $end\n"
        # not nested groups of input blocks (not parsed into a dict of dicts)
        self._notnested = ["$data", "$vec", "$ecp"]

    @property
    def template(self):
        return self._template

    @template.setter
    def template(self, value):
        if isinstance(value, dict):
            self._template = value
        elif value is None:
            self._template = {}
        else:
            raise TypeError("expected a dictionary but got {0:s}".format(type(value)))

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

    def write_input(self, mol=None, bs=None, core=None):
        '''
        Write a gamess input file under the name <inpfile> based on the
        information fstored in the dictionary <inpdict>.

        Args
        ----
        mol : Molecule
            Molecule object
        bs : BasisSet
            BasisSet object or a list of BasisSet objects
        '''

        inpstr = ""

        # write nested namelist groups
        for key, value in sorted(self.template.items()):
            if key not in  self._notnested:
                inpstr += " {0:<s}\n".format(key)
                for kkey, vvalue in sorted(value.items()):
                    inpstr += "    {k:s}={v:s}\n".format(k=kkey, v=str(vvalue))
                inpstr += self.end
        #write $data card
        inpstr += self.write_data(mol=mol, bsl=bs)

        with open(self.fname, "w") as finp:
            finp.write(inpstr)

    def write_data(self, mol=None, bsl=None):
        '''
        Return the $DATA part of the input based on the information in the $data
        dict.

        Args
        ----
        mol : Molecule
            Molecule object
        bs : BasisSet
            BasisSet object or a list of BasisSet objects
        '''

        # converrt the list of BasisSet object into a dict wilt element symbols
        # as keys

        if isinstance(bsl, list):
            bsd = {b.element : b for b in bsl}
        else:
            bsd = {bsl.element : bsl}

        data = ""
        data += " {0:s}\n".format("$data")
        data += "{0:s}\n".format(self.template["$data"]["title"])
        data += "{0:s}\n\n".format(self.template["$data"]["group"])
        if mol is not None:
            for atom in mol.unique():
                data += atom.gamess_rep()
                data += bsd[atom.symbol].write_gamess()
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

    def get_variable(self, rawstring):
        with open(selflogfile, 'r') as out:
            data = out.read()

        genre = re.compile(rawstring, flags=re.M)
        match = genre.search(data)
        if match:
            return float(match.group(1))

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
