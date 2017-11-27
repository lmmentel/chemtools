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

'''
Module for handling Gamess-US related jobs:
- Gamess          : running and submitting jobs, writing inputs,
- GamessInpParser : parsing the input file,
- GamessLogParser : parsing the output file,
- GamessDatParser : parsing data from the gamess PUNCH (.dat) file
'''

from __future__ import print_function

import argparse
import itertools
import math
import os
import re
import sys
from copy import copy
from subprocess import call
import numpy as np

from .calculator import Calculator, InputTemplate
from ..parsetools import slicebetween, sliceafter, parsepairs, getlines


class GamessUS(Calculator):

    '''Container object for Gamess-us jobs.'''

    def __init__(self, name="GamessUS", version="00", runopts=None, **kwargs):
        super(GamessUS, self).__init__(**kwargs)
        self.name = name
        self.inpext = '.inp'
        self.rungms = self.executable
        self.gmspath = os.path.dirname(self.executable)
        self.version = version
        self.runopts = runopts

        if os.path.isfile(os.path.join(self.gmspath, 'ddikick.x')):
            self.ddikick = os.path.join(self.gmspath, 'ddikick.x')

    @property
    def version(self):
        'Return the version.'
        return self._version

    @version.setter
    def version(self, value):
        'Check if the version exist and set it if it does.'

        versions = []
        for item in os.listdir(self.gmspath):
            match = re.match(r'[a-z]+\.([0-9]+)\.x', item)
            if match:
                versions.append(match.group(1))

        if value in versions:
            self._version = value
        else:
            raise IOError('GamessUS version {0:s} not found in {1:s}'.format(value, self.gmspath))

    @property
    def runopts(self):
        'Return the runopts'
        return self._runopts

    @runopts.setter
    def runopts(self, value):
        'Set the runopts dict'

        if value is None:
            self._runopts = ['1']
        else:
            self._runopts = value

    def get_command(self, inpfile):
        'Return the command to execute'

        return [self.executable, os.path.splitext(inpfile)[0], self.version] +\
            self.runopts

    def remove_dat(self, inpfile):
        '''
        Remove the gamess dat file if it exists in the scratch directory.
        '''

        datfile = os.path.splitext(inpfile)[0] + ".dat"
        if os.path.exists(os.path.join(self.scratch, datfile)):
            os.remove(os.path.join(self.scratch, datfile))

    def run(self, inpfile, logfile=None, remove_dat=True):
        '''
        Run a single gamess job interactively - without submitting to the
        queue.
        '''

        if remove_dat:
            self.remove_dat(inpfile)

        if logfile is None:
            logfile = os.path.splitext(inpfile)[0] + ".log"

        with open(logfile, 'w') as fobj:
            call(self.get_command(inpfile), stdout=fobj, stderr=fobj)

        return logfile

    def run_multiple(self, inputs):
        pass

    def accomplished(self, outfile):
        '''
        Return True if Gamess(US) job finished without errors.
        '''

        parser = GamessLogParser(outfile)
        return parser.accomplished()

    def parse(self, output, objective, regexp=None):
        '''
        Parser GAMESS(US) output file to get the objective.
        '''

        parser = GamessLogParser(output)

        if objective == "hf total energy":
            return parser.get_hf_total_energy()
        elif objective == "cisd total energy":
            energies = parser.get_energy_components('ci')
            return energies["TOTAL ENERGY"]
        elif objective == "correlation energy":
            energies = parser.get_energy_components('ci')
            return energies["TOTAL ENERGY"] - parser.get_hf_total_energy()
        elif objective == "regexp":
            return parser.get_variable(regexp)
        else:
            raise ValueError("unknown objective in prase {0:s}".format(objective))

    def write_input(self, fname, template=None, mol=None, basis=None,
                    core=None):
        '''
        Write the molpro input to "fname" file based on the information from
        the keyword arguments.

        Args:
            mol : :py:class:`chemtools.molecule.Molecule`
                Molecule object instance
            basis : dict or :py:class:`BasisSet <chemtools.basisset.BasisSet>`
                An instance of :py:class:`BasisSet <chemtools.basisset.BasisSet>`
                class or a dictionary of :py:class:`BasisSet <chemtools.basisset.BasisSet>`
                objects with element symbols as keys
            core : list of ints
                Molpro core specification
            template : :py:class:`str`
                Template of the input file
            fname : :py:class:`str`
                Name of the input file to be used
        '''

        out = ' $data\nbasis optimization\n{0:s}'.format(mol.symmetry)
        if mol.symmetry.lower() == 'c1':
            out += '\n'
        else:
            out += '\n\n'

        for i, atom in enumerate(mol.atoms):
            if i in mol.unique_labels:
                atombasis = basis[atom.symbol]
                out += atom.gamess_rep()
                out += atombasis.to_gamessus()

        out += ' $end'

        subs = {'basis': out, 'core': core}
        temp = InputTemplate(template)

        with open(fname, 'w') as finp:
            finp.write(temp.substitute(subs))

    def __repr__(self):
        return "\n".join([
            "<Gamess(",
            "\tname={},".format(self.name),
            "\tgmspath={},".format(self.gmspath),
            "\trungms={},".format(self.rungms),
            "\tversion={},".format(self.version),
            "\tscratch={},".format(self.scratch),
            "\trunopts={},".format(str(self.runopts)),
            ")>\n"])


class GamessInput(object):
    '''
    A class for parsing and writing gamess-us input files.
    '''

    def __init__(self, fname=None, parsed=None):
        '''
        Initialize the class.
        '''

        self.fname = fname
        self.end = " $end\n"
        self.parsed = parsed
        # not nested groups of input blocks (not parsed into a dict of dicts)
        self._notnested = ["$data", "$vec", "$ecp"]

    @property
    def parsed(self):
        'Return the `parsed` attribute'
        return self._parsed

    @parsed.setter
    def parsed(self, value):
        if isinstance(value, dict):
            self._parsed = value
        elif value is None:
            self._parsed = {}
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

        pat = re.compile(r'(?P<block>\$[a-zA-Z]{3,6})\s+(?P<entries>.*?)\$END',
                         flags=re.DOTALL | re.IGNORECASE)

        iterator = pat.finditer(inpstr)
        for match in iterator:
            if match.group("block").lower() not in self._notnested:
                self.parsed[match.group("block").lower()] = {}
                fields = [s.strip() for s in match.group("entries").split("\n")]
                for field in fields:
                    if not field.startswith("!"):
                        for line in field.split():
                            key, value = line.split("=")
                            self.parsed[match.group("block").lower()][key.lower()] = value
            elif match.group("block").lower() == "$data":
                self.parsed["$data"] = self.parse_data(match.group("entries"))
            elif match.group("block").lower() in ["$vec", "$ecp"]:
                self.parsed[match.group("block").lower()] = match.group("entries")
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

        block = re.compile(r'(?P<label>[a-zA-Z]{1,2}[0-9]{0,2})\s*' +
                           r'(?P<atomic>\d+\.\d+)' +
                           r'(?P<xyz>(\s+\-?\d+\.\d+){3})\s*' +
                           r'(?P<basis>.*?)\n\s*\n', flags=re.S)

        datadict = dict()

        datadict["title"] = datastr.split('\n')[0].strip()
        datadict["group"] = datastr.split('\n')[1].strip()
        datadict["atoms"] = list()

        itfound = block.finditer(datastr)
        for match in itfound:
            datadict["atoms"].append({
                'label': match.group('label'),
                'atomic': float(match.group('atomic')),
                'xyz': tuple(float(x) for x in match.group('xyz').split()),
                'basis': match.group('basis'),
            })

        return datadict

    def parsed2str(self):
        '''
        Return a string with formatted string containing the GAMESS(US) input
        based on the previously parsed data.
        '''

        inpstr = ""
        # write nested namelist groups
        for key, value in sorted(self.parsed.items()):
            if key not in self._notnested:
                inpstr += " {0:<s}\n".format(key)
                for kkey, vvalue in sorted(value.items()):
                    inpstr += "    {k:s}={v:s}\n".format(k=kkey, v=str(vvalue))
                inpstr += self.end
        # write $data card
        inpstr += self.data2str()
        return inpstr

    def data2str(self, header=True):
        '''
        Return the $DATA card of the input as a formatted string based on the
        previously parsed data.a
        '''

        data = " {0:s}\n".format("$data")
        if header:
            data += "{0:s}\n".format(self.parsed["$data"]["title"])
            data += "{0:s}\n".format(self.parsed["$data"]["group"])
            if self.parsed['$data']['group'].lower() != 'c1':
                data += '\n'
        for atom in self.parsed['$data']['atoms']:
            data += '{0:5s}{1:5.1f}{2:12.5f}{3:12.5f}{4:12.5f}\n'.format(
                atom['label'], atom['atomic'], atom['xyz'][0], atom['xyz'][1],
                atom['xyz'][2])
            data += atom['basis'] + '\n\n'
        data += self.end
        return data

    def write_input(self, filename):
        '''
        Write a gamess input file under the name <inpfile> based on the
        information fstored in the dictionary <inpdict>.

        Args:
          filename : str
            Name of the input file
        '''

        with open(filename, "w") as finp:
            finp.write(self.parsed2str())

    def write_with_vec(self, filename, vecstr, keepscftyp=False):
        '''
        Write new gamess input based on exisiting input and natural orbitals
        from $VEC section in PUNCH file.


        Write a new input file for gamess based on previously prased and/or
        constructed dictionary containing input specification and append
        orbitals from a previous run as starting orbitals. The starting
        orbitals should be stored in an ASCII PUNCH file whose name is given
        under "datfile" variable.

        Args:
            filename : str
                Name of the input file to be created
            vecstr : str
                String with the vectors to be written
        '''

        parsed = copy(self.parsed)

        if not keepscftyp:
            self.parsed['$contrl']['scftyp'] = 'none'
        self.parsed['$guess'] = {}
        self.parsed['$guess']['guess'] = 'moread'
        self.parsed['$guess']['norb'] = str(get_naos_nmos(vecstr)[1])

        with open(filename, "w") as finp:

            finp.write(self.parsed2str())

            finp.write("\n $vec\n")
            finp.write(vecstr)
            finp.write("\n $end\n")
            finp.close()

        self.parsed = parsed

    def set_gamess_input(self, dinp, mol, bs, code, core):

        if "$contrl" in dinp.keys():
            dinp["$contrl"]["icharg"] = mol.charge
            dinp["$contrl"]["mult"] = mol.multiplicity
            if mol.multiplicity == 1:
                dinp["$contrl"]["scftyp"] = "rhf"
            elif mol.multiplicity > 1:
                dinp["$contrl"]["scftyp"] = "rohf"
            if code["method"].lower() == "hf":
                dinp["$contrl"] = {key: value for key, value in dinp["$contrl"].items() if key != "cityp"}
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

    def print_parsed(self):
        '''
        Neat print of the dictionary with parsed input
        '''

        for key, value in sorted(self.parsed.items()):
            if isinstance(value, dict):
                print(key)
                for kkey, vvalue in sorted(value.items()):
                    print("\t{0:<10s} : {1:}".format(kkey, vvalue))
            else:
                print(key, '\n', value)


class GamessLogParser(object):
    '''
    Methods for parsing gamess-us log file.
    '''

    def __init__(self, log):
        self.logfile = log

    @property
    def logfile(self):
        'Return the value of the `logfile` attribute'
        return self._logfile

    @logfile.setter
    def logfile(self, value):
        'Check if the logfile exists first.'
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
        Parsing function based on querries in the form of regular expressions
        that return the match object.
        '''

        cpatt = re.compile(regex)
        with open(self.logfile, 'r') as log:
            content = log.read()
        return cpatt.search(content)

    def accomplished(self):
        '''Check if a job teminated normally.'''

        regex = r'TERMINATED NORMALLY'
        match = self.parse(regex)
        return match is not None

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
        '''
        Get the orbital index of homo orbital (indexing starts from zero).
        '''

        if int(self.get_electrons()) % 2 == 0:
            return int(self.get_electrons()) / 2 - 1
        else:
            sys.exit("open shell handling not implemented")

    def get_number_of_atoms(self):
        '''Get total number of atoms from gamess log file.'''

        regex = r'TOTAL NUMBER OF ATOMS\s+=\s*(?P<nat>\d+)'
        match = self.parse(regex)
        if match:
            return int(match.group("nat"))

    def get_number_of_aos(self):
        '''
        Get the number of primitive cartesian gaussian basis functions from
        Gamess log file
        '''

        regex = r'NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS =\s*(?P<nao>\d+)'
        match = self.parse(regex)
        if match:
            return int(match.group("nao"))

    def get_number_of_core_mos(self):
        '''
        Get the number of core molecular orbitals from the Gamess log file
        '''

        regex = r'NUMBER OF CORE MOLECULAR ORBITALS\s*=\s*(?P<ncore>\d+)'
        match = self.parse(regex)
        if match:
            return int(match.group("ncore"))

    def get_number_of_mos(self):
        '''
        Get the number of molecular orbitals from Gammess log file.
        '''

        ispher_patt = r'ISPHER=\s*(?P<ispher>\-?\d{1}).*'
        var_space_patt = r'.*VARIATION SPACE IS\s*(?P<nmo>\d+).*'
        c_ispher = re.compile(ispher_patt)
        c_var_space = re.compile(var_space_patt)

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
        '''
        Get number of linearly dependent combinations dropped.
        '''

        regex = r'NUMBER OF LINEARLY DEPENDENT MOS DROPPED=\s*(?P<lindep>\d+)'
        match = self.parse(regex)
        if match:
            return int(match.group('lindep'))

    def get_scf_type(self):
        '''
        Get the information on SCFTYP used in the gamess job.
        '''

        regex = r'.*SCFTYP=(?P<scftyp>[A-Z]+)'
        match = self.parse(regex)
        if match:
            return match.group("scftyp")

    def get_cc_type(self):
        '''
        Get the information on CCTYP used in the gamess job.
        '''

        regex = r'.*CCTYP =(?P<cctyp>[A-Z\(\)\-]+)'
        match = self.parse(regex)
        if match:
            return match.group("cctyp")
        else:
            return None

    def get_ci_type(self):
        '''
        Get the information on CITYP used in the gamess job.
        '''

        regex = r'.*CITYP =(?P<cityp>[A-Z]+)'
        match = self.parse(regex)
        if match:
            return match.group("cityp")

    def get_mplevel(self):
        '''
        Get the information on MPLEVL used in the gamess job.
        '''

        regex = r'.*MPLEVL=\s*(?P<mplevl>\d+)'
        match = self.parse(regex)
        if match:
            return int(match.group("mplevl"))

    def get_hf_total_energy(self):
        '''
        Return the total HF energy.
        '''

        regex = r'FINAL R[O]*HF ENERGY IS\s+(?P<energy>\-?\d+\.\d+)'
        match = self.parse(regex)
        if match:
            return float(match.group("energy"))

    def get_cc_total_energy(self):
        '''
        Independently of CCTYP value return the "HIGHEST LEVEL RESULT" total
        energy.
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
                header = '{0:<5s} CI PROPERTIES'.format(self.get_ci_type())
        else:
            sys.exit("Wrong method in <get_energy_components>: {0:s}".format(method))

        with open(self.logfile, 'r') as log:
            data = log.readlines()

        return parsepairs(sliceafter(data, header, 22))

    def get_lz_values(self):
        '''
        Get values of Lz operator for every orbital and translate them into
        labels of orbitals.
        '''

        mos = r'\s*MO\s*(?P<index>\d*)\s*\(\s*(?P<shell>\d+)\s*\)'
        weights = r'\s*HAS LZ\(WEIGHT\)=\s*(?P<lz>-?\d*\.\d+)'
        lz_patt = re.compile(mos + weights)

        locstr = self.get_loc_strings('lz values')
        lines = getlines(self.logfile, locstr)

        res = list()
        for line in lines:
            match = lz_patt.search(line)
            if match:
                index = int(match.group('index')) - 1
                shell = int(match.group('shell'))
                lz = int(abs(float(match.group('lz'))))
                res.append({"index": index, "shell": shell, "lz": lz})
        return res

    def get_ao_labels(self, orbs='hf orbs'):
        '''
        Retrieve the information about atomic basis set
        '''

        indsym = r'\s*(?P<index>\d+)\s{2}(?P<symbol>[a-zA-Z]{1,2})'
        cencom = r'\s*(?P<center>\d+)\s*(?P<component>[A-Z]{1,})'
        ao_patt = re.compile(indsym + cencom)

        istart = 3
        if orbs == 'local orbs':
            istart = 1
        locstr = self.get_loc_strings(orbs)
        lines = getlines(self.logfile, locstr)
        lines = lines[istart:istart + self.get_number_of_aos()]

        res = list()
        for line in lines:
            match = ao_patt.search(line)
            if match:
                index = int(match.group('index')) - 1
                symbol = match.group('symbol')
                center = int(match.group('center'))
                component = match.group('component')
                res.append({
                    "index": index,
                    "symbol": symbol,
                    "center": center,
                    "component": component})
        return res

    def get_variable(self, rawstring):
        'wrapper around a regex search method'

        with open(self.logfile, 'r') as out:
            data = out.read()

        genre = re.compile(rawstring, flags=re.M)
        match = genre.search(data)
        if match:
            return float(match.group(1))

    def get_loc_strings(self, name):

        locators = {
            'csfs'      : [('DETERMINANT CONTRIBUTION TO CSF', 2),
                           (' ...... END OF -DRT- GENERATION ......', -5)],
            'ci coeffs' : [('      CSF      COEF    OCCUPANCY (IGNORING CORE)', 2),
                           (' ...... END OF CI-MATRIX DIAGONALIZATION ......', 0)],
            'initial orbs' : [('INITIAL GUESS ORBITALS', 3),
                              (' ...... END OF INITIAL ORBITAL SELECTION ......', 0)],
            'hf orbs'   : [('EIGENVECTORS', 3),
                           (' ...... END OF {} CALCULATION ......'.format(self.get_scf_type()), 0)],
            'ci orbs'   : [('NATURAL ORBITALS IN ATOMIC ORBITAL', 3),
                           (' ...... END OF DENSITY MATRIX CALCULATION ......', 0)],
            'local orbs': [('LOCALIZED ORBITALS', 3),
                           (' LOCALIZATION ...', 0)],
            'lz values' : [('LZ VALUE ANALYSIS FOR', 2),
                           ('EIGENVECTORS', -2)],
        }

        if name in locators.keys():
            return locators[name]
        else:
            raise ValueError('wrong "name", should be one of: {}'.format(', '.join(locators.keys())))

    def get_orbital_labels(self, orbs, minlength=5):
        '''Return the labels of the specified orbitals

        Args:
          orbs : str
            String key for extracting locators
          minlength : int
            Minimal length of a line to be kept

        Returns:
          out : list
        '''

        locstr = self.get_loc_strings(orbs)
        lines = getlines(self.logfile, locstr)
        lines = [l for l in lines if len(l) > minlength]

        nao = self.get_number_of_aos()

        inds = list(itertools.chain(*[l.split() for l in lines[::nao + 3]]))
        evs = list(itertools.chain(*[l.split() for l in lines[1::nao + 3]]))
        syms = list(itertools.chain(*[l.split() for l in lines[2::nao + 3]]))

        return zip(inds, evs, syms)

    def get_ci_coeffs(self):
        '''
        Parse CI coefficients from a list of lines containing the output of
        GAMESS(US) calculation

        Returns:
          coeffs : dict
            Dictionary with parsed coeficients and their indexes in the
            internal GAMESS(US) ordering.
        '''

        locstr = self.get_loc_strings('ci coeffs')
        lines = getlines(self.logfile, locstr)

        patt = re.compile(r'\s*(\d+)\s*(-?\d*\.\d+)')
        index = []
        coeff = []
        for line in lines:
            match = patt.search(line)
            if match:
                index.append(match.group(1))
                coeff.append(match.group(2))

        return dict(zip([int(x) for x in index], [float(x) for x in coeff]))

    def get_csfs(self, withcoeffs=True):
        '''
        Parse CSF information from a list of lines from GAMESS(US) output.
        '''

        locstr = self.get_loc_strings('csfs')
        lines = getlines(self.logfile, locstr)

        ne = self.get_electrons() - 2 * self.get_number_of_core_mos()

        csfre = r'^\s*(CSF\s*(?P<csfno>\d+)\:)?'
        coeffre = r'\s*C\(\s*(?P<detno>\d+)\)=\s*(?P<coeff>-?\d*\.\d+)'
        orbsre = r'\s+\:\s{2}(?P<orbs>.*)'
        patt = re.compile(csfre + coeffre + orbsre)

        csfstr = "".join(lines).lstrip()
        csfs = [s for s in re.split(r'CASE VECTOR =\s*\d+\s*', csfstr) if len(s) > 0]

        out = []
        for csfstr in csfs:
            endvec = csfstr.find('FOR MS=')
            vector = csfstr[:endvec].strip()
            csf = {'vector': vector}
            csf['dets'] = []
            csfstart = csfstr.find('CSF')
            for i, line in enumerate(csfstr[csfstart:].split('\n')):
                match = patt.search(line)
                if match:
                    gd = match.groupdict()
                    # in the first line of the csf specification there is the CSF number
                    # (index) that should be stored
                    if i == 0:
                        csf['no'] = int(gd['csfno'])
                        spatial = [abs(int(x)) for x in detsplit(gd['orbs'], ne)]
                        isort = np.argsort(np.array(spatial))
                        csf['spatial'] = [spatial[i] for i in isort]
                    spin = det_to_spin_coupling(gd['orbs'], ne)
                    csf['dets'].append((int(gd['detno']), float(gd['coeff']), "".join([spin[i] for i in isort])))
            out.append(csf)

        if withcoeffs:
            coeffs = self.get_ci_coeffs()
            merged = [dict(csf.items() + [('coeff', coeffs.get(csf['no'], None))]) for csf in out]
            return merged
        else:
            return out


def detsplit(det, norb):
    '''
    Split a string representation of the spatial part of the determinant `det`
    into `n` parts.

    Args:
        det : str
            String representation of the determinant
        norb : int
            Number of orbitals in the determinant

    Returns:
        out : list
            List of determinant components
    '''

    out = det.split()
    if len(out) != norb:
        nchar = int(math.floor(len(det) / float(norb)))
        out = [det[nchar * i: nchar * i + nchar] for i in range(norb)]
    return out


def det_to_spin_coupling(det, norb):
    '''
    Convert a string representaiton of a determinant to a string with spin
    coupling.

    It is assumed here that a minus sign means beta spin. `a` will denote
    :math:`\alpha` and `b` will denote :math:`\beta`

    Args:
      det : str
        String representation of a determinant, e.g. `-3  4  2 -1`
      norb : int
        Number of orbitals that should be present in `det` equal to the number
        of active electrons

    Returns:
      sc : str
        String representation of a spin coupling
    '''

    sc = ''

    detspl = detsplit(det, norb)

    for oi in detspl:
        if '-' in oi:
            sc += 'b'
        else:
            try:
                sc += 'a'
            except:
                sc += '*'
    return sc


class GamessDatParser(object):
    'Parser for the GAMESS(US) dat (.F10) file'

    def __init__(self, datfile):
        self.datfile = datfile

    @property
    def datfile(self):
        'Return the value of the `datfile` attribute'
        return self._datfile

    @datfile.setter
    def datfile(self, value):
        if os.path.exists(value):
            self._datfile = value
        else:
            raise ValueError("File: {} does not exist".format(value))

    def get_occupations(self):
        '''
        Parse the occupation numbers from the ascii PUNCH file (.dat).
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
            datastr = slicebetween(dat.read(), '$DATA', '$END')
        gip = GamessInput()
        return gip.parse_data(datastr.lstrip(' \n'))

    def get_vec_string(self, method):
        '''
        Parse the natural orbitals from the ascii PUNCH file (.dat).

        Args:
          method : str
            acceptable values are:
                - for scf orbitals: "scf", "hf", "rhf", "rohf", "uhf", "gvb"
                - for mcscf orbitals: "mcscfmos", "mcscfnos"
                - for ci orbitals: "ci", "aldet", "fsoci", "guga", "genci",
                  "ormas""
                - for localized orbitals: `local`

        Returns:
            orbs : numpy.array
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
        elif method.lower() == 'local':
            try:
                orbi = data.index('LOCALIZED ORBITALS')
            except:
                raise ValueError('LOCALIZED orbitals header not found, check dat file')
        else:
            raise ValueError("Don't know what to do with: '{0:s}'".format(method))
        vecstr = slicebetween(data[orbi:], '$VEC', '$END').strip(" \n\t\r")
        if vecstr == "":
            raise ValueError("No $VEC section found for method: '{0:s}'".format(method))
        else:
            # add a single whitespace to account for the removed space in strip
            return ' ' + vecstr

    def get_orbitals(self, method):
        '''
        Parse the natural orbitals from the ascii PUNCH file (.dat).

        Args:
          method : str
            acceptable values are:
                - for scf orbitals: "scf", "hf", "rhf", "rohf", "uhf", "gvb"
                - for mcscf orbitals: "mcscfmos", "mcscfnos"
                - for ci orbitals: "ci", "aldet", "fsoci", "guga", "genci",
                  "ormas""
                - for localized orbitals: `local`

        Returns:
            orbs : numpy.array
              2D numpy array of shape (naos, nmos) containing ortbials
              expansion coefficients
        '''

        vecstr = self.get_vec_string(method)
        return self.parse_orbitals(vecstr)

    @staticmethod
    def parse_orbitals(vecstr, clength=15):
        '''
        Parse dat orbitals into numpy array of the shape (naos, nmos)

        Parse orbitals given a string obtained from the $VEC section of the
        PUNCH file (.dat) into a numpy array of the shape (naos, nmos) where
        "naos" and "nmos" are the number of atomic orbitals and number of
        molecular orbitals respectively. The shape is deduced from the last
        line of the string.

        Args:
          vecstr : str
            string with the contents fo the gamess $VEC block
          clength : int
            total length of the coefficient string as stored by the gamess
            format, by default gamess stores orbital coefficients in
            `e15.8` fortran format so the total length is 15.

        Returns:
          orbs : numpy.array
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
                nitems = int(len(orblines[counter][5:])) // clength
                orbs[5 * j: 5 * (j + 1), i] = \
                    [float(orblines[counter][5 + 15 * n: 5 + 15 * (n + 1)])
                     for n in range(nitems)]
        return orbs


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

    veclines = vecstr.split('\n')
    noveclines = len(veclines)
    lineit = iter(veclines)
    nlines = 0

    while next(lineit)[:2].strip() == '1':
        nlines += 1

    if nlines == 0:
        raise ValueError("'nlines' cannot be zero, check vecstr in 'get_naos_nmos'")

    naos = 5 * (nlines - 1) + len(vecstr.split('\n')[nlines - 1][5:]) // clength
    nmos = noveclines // nlines
    return naos, nmos, nlines


def to_gamess_vec(coeffs, cfmt='15.8e'):
    '''
    Given a 2D array return a string with the elements converted to the
    GAMESS(US) $VEC format

    Args:
        coeffs : numpy.array
            A two-dimensional array of type `float`

        cmft : str
            Default format for the elements of `coeffs`

    Returns:
        vec : str
            A string representation of the `coeffs` using GAMESS(US) $VEC
            format
    '''

    if coeffs.dtype != np.float64:
        raise ValueError('array should be of type np.flaot64, got: {}'.format(coeffs.dtype))

    nao, nmo = coeffs.shape

    if nao % 5 == 0:
        nlines = nao / 5
    else:
        nlines = 1 + nao / 5

    vec = ''
    for i in range(nmo):
        if i + 1 >= 100:
            ilab = i + 1 % 100
        else:
            ilab = i + 1
        for l in range(nlines):
            if l + 1 > 1000:
                llab = l + 1 % 1000
            else:
                llab = l + 1
            if l < nlines - 1:
                coeffrow = "".join(["{0:{1:}}".format(coeffs[j + l * 5, i], cfmt) for j in range(5)])
            else:
                coeffrow = "".join(["{0:{1:}}".format(coeffs[j + l * 5, i], cfmt) for j in range(0, 5 - ((l + 1) * 5 - nao))])
            vec += "{0:2d}{1:3d}{2:s}\n".format(ilab, llab, coeffrow)
    return vec


def writeorbinp():
    '''
    A script for writing the GAMESS(US) input with guess orbitals from a
    previous run.
    '''

    parser = argparse.ArgumentParser()

    parser.add_argument('inpfile')
    parser.add_argument('datfile')
    parser.add_argument('orbitals', default='ci')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()

    gip = GamessInput(args.inpfile)
    gip.parse()

    gdp = GamessDatParser(args.datfile)
    vectr = gdp.get_vec_string(args.orbitals)

    if args.output is None:
        args.output = os.path.splitext(args.inpfile)[0] + '_{}.inp'.format(args.orbitals)
    gip.write_with_vec(args.output, vectr)
    print('Input written to: {}'.format(args.output))
