#!/usr/bin/env python
from basis import Basis
from molecule import Molecule
from subprocess import Popen
import numpy as np
import os
import re
import sys


class Gamess(object):

    '''Container object for Gamess-us jobs.'''

    def __init__(self, molecule=None, basis=None, inputname=None,
                 workdir=os.getcwd(), template=None):
        self.molecule = molecule
        self.basis    = basis
        self.workdir  = workdir
        self.scratch  = "/home/lmentel/scratch"
        self.template = template
        if inputname:
            self.filebase   = os.path.splitext(inputname)[0]
        else:
            self.filebase   = "_".join([self.molecule.name, self.basis.name])
        self.inputfile  = self.filebase + ".inp"
        self.outputfile = self.filebase + ".log"
        self.datfile    = self.filebase + ".dat"
        self.twoeaofile = self.filebase + ".F08"
        self.twoemofile = self.filebase + ".F09"
        self.dictionary = self.filebase + ".F10"
        self.rdm2file   = self.filebase + ".F15"

    def run(self, executable, version="00", nproc="1", remove_dat=True):

        '''Run a single gamess job interactively - without submitting to the
        queue.'''

        # should add checking for the existance of the .dat file before running

        if remove_dat:
            if os.path.exists(os.path.join(self.scratch, self.datfile)):
                os.remove(os.path.join(self.scratch, self.datfile))

        out = open(self.outputfile, 'w')
        process = Popen([executable, self.inputfile, version, nproc], stdout=out, stderr=out)
        process.wait()
        out.close()

    def run_startno(self, executable, version="00", nproc="1"):

        '''Run gamess calculation starting ffrom the orbitals in the dat file.'''

        # should add checking for the existance of the .dat file before running

        startnoinput = self.write_startno_input()
        startnooutput = os.path.splitext(startnoinput)[0] + ".log"

        out = open(startnooutput, 'w')
        process = Popen([executable, startnoinput, version, nproc], stdout=out, stderr=out)
        process.wait()
        out.close()

    def submit(self):

        '''Submit a single gamess job into the queue.'''

        pass

    def norbitals(self, spherical=True):

        '''Get the total number of atomic orbitals (spherical of cartesian)
           from the basis set information.'''

        nao = 0
        for atom in self.molecule.atoms:
            nao += self.basis.get_number_of_aos(atom.atomic, spherical)
        return nao

    def get_atoms_and_basis(self):

        '''Print all the atoms with their corresponding basis set into a format
        of Gamess-us $DATA card and return it as a string "atandbas".'''

        atandbas = ""
        atoms = [self.molecule.atoms[i] for i in self.molecule.unique]
        for atom in atoms:
            atandbas = atandbas + "{s:<4s}{a:>5.2f}{x:>15.6f}{y:>15.6f}{z:>15.6f}\n".format(s=atom.symbol,
                        a=float(atom.atomic), x=atom.xyz.x, y=atom.xyz.y, z=atom.xyz.z)
            atandbas = atandbas + self.basis.get_basis(atom.atomic) + '\n'
        atandbas = atandbas + ' $END'
        return atandbas

    def write_input(self, spherical=False, core=0):

        '''Write gamess input file based on information on molecule and basis set.'''

        with open(os.path.join(self.workdir,self.template), 'r') as tmp:
            contents = tmp.read()
        contents = re.sub(r'icharg=(\d+)?', r'icharg={0:s}'.format(str(self.molecule.charge)), contents, flags=re.I)
# write also ispher variable
        contents = re.sub(r'nact=(\d+)?', r'nact={0:<d}'.format(self.norbitals(spherical)-core), contents, flags=re.I)
        contents = re.sub(r'nels=(\d+)?', r'nels={0:<d}'.format(self.molecule.nele()-core*2), contents, flags=re.I)
        contents = re.sub(r'ncore=(\d+)?', r'ncore={0:<d}'.format(core), contents, flags=re.I)
        contents = re.sub("atomsandbasis", self.get_atoms_and_basis(), contents, flags=re.I)
        inp = open(self.inputfile, 'w')
        inp.write(contents)
        inp.close()

    def write_startno_input(self):

        '''Write gamess input file with starting orbitals from a previous run,
            stored in the self.datname file.'''

        gp = GamessParser(self.outputfile)

        with open(self.inputfile, 'r') as inp:
            inpcontent = inp.read()

        inpcontent = re.sub(r'scftyp=[A-Za-z]*', r'scftyp=none', inpcontent, flags=re.I)

        startnofile = self.filebase + "_NO.inp"
        newinp = open(startnofile, "w")

        newinp.write(inpcontent)
        newinp.write(" $GUESS guess=moread norb={:<d} $END\n\n".format(gp.get_number_of_mos()))
        newinp.write(" $VEC\n")

        dat    = open(self.datfile, 'r')

        citype = gp.get_ci_type()

        if citype == "GUGA":
            datheader = 'GUGA-CI'
        elif citype == "ORMAS":
            datheader = "- - - NO-S OF CI STATE"

        line = dat.readline()
        while not datheader in line and line != "":
            line = dat.readline()

        while not "$VEC" in line and line != "":
            line = dat.readline()

        while not "$END" in line and line != "":
            line = dat.readline()
            newinp.write(line)

        dat.close()
        newinp.close()
        return startnofile

    def __repr__(self):
        return "Gamess-US job object:\n\tMolecule  : {0:s}\n\tBasis     : {1:s}\n\tInput file: {2:s}".format(self.molecule.name, self.basis.name, self.inputfile)


class GamessParser(object):

    '''Object holding tools for parsing gmaess-us log file.'''

    def __init__(self, log=None, inp=None):
        if log:
            self.logfile    = log
            i = self.logfile.index("log")
            self.filebase   = self.logfile[:i-1]
            self.inputfile  = self.filebase + ".inp"
        elif inp:
            self.inputfile  = inp
            self.filebase   = os.path.splitext(self.inputfile)[0]
            self.logfile    = self.filebase + ".log"
        self.logexists()
        self.datfile    = self.filebase + ".dat"
        self.twoeaofile = self.filebase + ".F08"
        self.twoemofile = self.filebase + ".F09"
        self.dictionary = self.filebase + ".F10"
        self.rdm2file   = self.filebase + ".F15"

    def logexists(self):

        '''Check if the log file exists.'''

        if os.path.exists(self.logfile):
            return True
        else:
            sys.exit("Gamess log file: {0:s} doesn't exist in {1:s}".format(
                     self.logfile, os.getcwd()))

    def terminatedOK(self):

        '''Check if a job teminated normally.'''

        patt = r'TERMINATED NORMALLY'
        compatt = re.compile(patt)
        with open(self.logfile, 'r') as log:
            lines = log.read()
        match = compatt.search(lines)
        if match:
            return True
        else:
            return False


    def get_charge(self):

        '''Get total charge.'''

        patt = r'CHARGE OF MOLECULE\s+=\s*(?P<charge>\d+)'
        compatt = re.compile(patt)
        with open(self.logfile, 'r') as log:
            lines = log.read()
        match = compatt.search(lines)
        if match:
            return int(match.group("charge"))

    def get_electrons(self):

        '''Get number of electrons.'''

        patt = r'NUMBER OF ELECTRONS\s+=\s*(?P<nele>\d+)'
        compatt = re.compile(patt)
        with open(self.logfile, 'r') as f:
            lines = f.read()
        match = compatt.search(lines)
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

        with open(self.logfile, 'r') as log:
            data = log.read()

        nat_re = r'TOTAL NUMBER OF ATOMS\s+=\s*(?P<nat>\d+)'
        compre   = re.compile(nat_re)
        match    = compre.search(data)
        if match:
            return int(match.group("nat"))

    def get_number_of_aos(self):

        '''Get the number of primitive cartesian gaussian basis functions from
        Gamess log file'''

        patt    = r'NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS =\s*(?P<nao>\d+)'
        compatt = re.compile(patt)
        with open(self.logfile, 'r') as f:
            lines = f.read()
        match = compatt.search(lines)
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

        with open(self.logfile, 'r') as log:
            lines = log.read()

        linre = re.compile(r'NUMBER OF LINEARLY DEPENDENT MOS DROPPED=\s*(\d+)')

        match = linre.search(lines)
        if match:
            return match.group(1)

    def get_scf_type(self):

        '''Get the information on SCFTYP used in the gamess job.'''

        with open(self.logfile, 'r') as log:
            data = log.read()

        scftyp_re = r'.*SCFTYP=(?P<scftyp>[A-Z]+)'
        compre   = re.compile(scftyp_re)
        match    = compre.search(data)
        if match:
            return match.group("scftyp")

    def get_cc_type(self):

        '''Get the information on CCTYP used in the gamess job.'''

        with open(self.logfile, 'r') as log:
            data = log.read()

        cctyp_re = r'.*CCTYP =(?P<cctyp>[A-Z]+)'
        compre   = re.compile(cctyp_re)
        match    = compre.search(data)
        if match:
            return match.group("cctyp")

    def get_ci_type(self):

        '''Get the information on CITYP used in the gamess job.'''

        with open(self.logfile, 'r') as log:
            data = log.read()

        cityp_re = r'.*CITYP =(?P<cityp>[A-Z]+)'
        compre   = re.compile(cityp_re)
        match    = compre.search(data)
        if match:
            return match.group("cityp")

    def get_hf_total_energy(self):

        '''Return the total HF energy.'''

        with open(self.logfile, 'r') as log:
            data = log.read()

        hfre = re.compile(r'FINAL RHF ENERGY IS\s+(?P<energy>\-?\d+\.\d+)', flags=re.M)
        match = hfre.search(data)
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

    def get_onee_size(self):
        '''Get the size of the vector holding upper (or lower) triangle
           of a square matrix of size nmos.'''
        n = self.get_number_of_mos()
        return n*(n+1)/2

    def get_twoe_size(self):
        '''Get the size of the 1d vector holding upper (or lower) triangle
           of a supermatrix of size nmos (2RDM and two-electrons integrals) .'''
        n = self.get_onee_size()
        return n*(n+1)/2
