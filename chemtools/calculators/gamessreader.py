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
    GamessReader    : reading gamess binary files.
'''

import os
import sys
import numpy as np
from collections import namedtuple

from .gamessus import GamessLogParser

# add an option to choose which reader should be used for dictionary file
# and sequential files, the options are:
# - wrapped fortran code that need to be compiled and installed or
# - native reader written in python using DictionaryFile and SequentialFile
# classes

# fortran modules nedded for GamessReader Class
# from gamessusfortranext import twoe


class GamessFortranReader(object):
    '''
    Class for holding method for reading gamess binary files:
        $JOB.F08 : two electron integrals over AO's,
        $JOB.F09 : two electron integrals over MO's,
        $JOB.F10 : the dictionary file with 1-e integrals, orbitals etc.,
        $JOB.F15 : GUGA and ORMAS two-electron reduced density matrix,
    '''

    def __init__(self, log):
        self.logfile = log
        i = self.logfile.index("log")
        self.filebase = self.logfile[: i - 1]
        self.twoeaofile = self.filebase + ".F08"
        self.twoemofile = self.filebase + ".F09"
        self.rdm2file = self.filebase + ".F15"
        self.gp = GamessLogParser(log=self.logfile)

    def get_onee_size(self, aos=True):
        '''
        Get the size of the vector holding upper (or lower) triangle
        of a square matrix of size naos or nmos.
        '''
        if aos:
            n = self.gp.get_number_of_aos()
        else:
            n = self.gp.get_number_of_mos()
        return n * (n + 1) / 2

    def get_twoe_size(self, aos=False):
        '''
        Get the size of the 1d vector holding upper (or lower) triangle
        of a supermatrix of size nmos (2RDM and two-electrons integrals).
        '''
        n = self.get_onee_size(aos)
        return n * (n + 1) / 2

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
            twoe.integrals.readinao(rdm2, self.rdm2file)
            return rdm2
        else:
            sys.exit("File '{0:s}' doesn't exist, exiting...".format(self.rdm2file))

    def read_twoeao(self, filename=None, nmo=None):

        '''Read the two electron integrals from the gamess-us file'''

        # initialize numpy array to zeros
        ints = np.zeros(self.get_twoe_size(aos=True), dtype=float)
        # use gamess module to read the integrals from the file -filename-
        if filename:
            if os.path.exists(filename):
                twoe.integrals.readinao(ints, filename)
                return ints
            else:
                sys.exit("File '{0:s}' doesn't exist, exiting...".format(filename))
        elif os.path.exists(self.twoeaofile):
            twoe.integrals.readinao(ints, self.twoeaofile)
            return ints
        else:
            raise OSError("File '{0:s}' doesn't exist, exiting...".format(self.twoeaofile))

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


def ijkl(i, j, k, l):
    '''
    Based on the four orbital indices i, j, k, l return the address
    in the 1d vector.
    '''

    ij = max(i, j) * (max(i, j) + 1) / 2 + min(i, j)
    kl = max(k, l) * (max(k, l) + 1) / 2 + min(k, l)

    return max(ij, kl) * (max(ij, kl) + 1) / 2 + min(ij, kl)


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


def print_twoe(twoe, nbf):
    '''Print the two-electron integrals.'''

    ij = 0
    for i in xrange(nbf):
        for j in xrange(i + 1):
            ij += 1
            kl = 0
            for k in xrange(nbf):
                for l in xrange(k + 1):
                    kl += 1
                    if ij >= kl:
                        if abs(twoe[ijkl(i, j, k, l)]) > 1.0e-10:
                            print("{0:3d}{1:3d}{2:3d}{3:3d} {4:25.14f}".format(
                                i, j, k, l, twoe[ijkl(i, j, k, l)]))


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

        length = int(length)

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

    def __init__(self, filename, logfile=None):
        '''
        Initialize the class with the superclass method.
        '''
        super(SequentialFile, self).__init__(filename)
        if logfile is None:
            logfile = os.path.splitext(filename)[0] + ".log"
        glp = GamessLogParser(log=logfile)
        self.nao = glp.get_number_of_aos()
        self.nmo = glp.get_number_of_mos()
        self.core = glp.get_number_of_core_mos()
        if self.nao < 255:
            self.large_labels = False
        else:
            self.large_labels = True

    def ijkl(self, i, j, k, l):
        '''
        Based on the four orbital indices i,j,k,l return the address
        in the 1d vector.
        '''

        ij = max(i, j) * (max(i, j) - 1) / 2 + min(i, j)
        kl = max(k, l) * (max(k, l) - 1) / 2 + min(k, l)
        return max(ij, kl) * (max(ij, kl) - 1) / 2 + min(ij, kl) - 1

    def get_index_buffsize(self, buff_size, int_size):
        ''' Return the index buffer size for reading 2-electron integrals'''

        if int_size == 4:
            if self.large_labels:
                return 2 * buff_size
            else:
                return buff_size
        elif int_size == 8:
            if self.large_labels:
                return buff_size
            else:
                return (buff_size + 1) / 2
        else:
            raise ValueError('wrong "int_size": {}'.format(int_size))

    def readseq(self, buff_size=15000, int_size=8, mos=False,
                skip_first=False):
        '''
        Read FORTRAN sequential unformatted file with two-electron
        quantities:

        - two electron integrals over AO's: .F08 file
        - two electron integrals over MO's: .F09 file
        - elements of the two particle density matrix: .F15 file

        Args:
            buff_size : int
                size of the buffer holding values to be read, in
                gamessus it is stored under ``NINTMX`` variable and in
                Linux version is equal to 15000 which is the default
                value
            large_labels : bool
                a flag indicating if large labels should were used, if
                largest label ``i`` (of the MO) is ``i<255`` then
                ``large_labels`` should be ``False``
                (case ``LABSIZ=1`` in gamessus), otherwise set to ``True``
                (case ``LABSIZ=2`` in gamess(us),
            skip_first : bool
                skips the first record of the file is set to True,

        Returns:
            numpy 1D array holding the values
        '''

        indexBuffSize = self.get_index_buffsize(buff_size, int_size)

        self.seek(0)
        if mos:
            nmo = self.nmo - self.core
            nt = self.nmo * (self.nmo+1)/2
            if skip_first:
                #n1 = nmo*(nmo+1)/2
                #self.seek(8+8*n1)
                self.seek(self.tell() + 4)
                Hoe = self.read('f8', shape=(nt, ))
                # advance four trailing bytes
                self.seek(self.tell() + 4)

        else:
            nt = self.nao * (self.nao + 1) / 2

        ints = np.zeros(nt * (nt + 1) / 2, dtype=float, order='F')

        int_type = np.dtype('i' + str(int_size))
        index_buffer = np.zeros(indexBuffSize, dtype=int_type, order='F')
        value_buffer = np.zeros(buff_size, dtype=float, order='F')

        length = 1
        while length > 0:
            self.seek(self.tell() + 4)
            length = self.read(int_type)
            if length > buff_size:
                raise ValueError('the read record length: {0:10d} is greater that the buffer size {1:10d}'.format(int(length), buff_size))

            index_buffer = self.read(int_type, shape=(indexBuffSize, ))
            value_buffer = self.read('f8', shape=(buff_size, ))

            for m in range(1, abs(length) + 1):
                if int_size == 4:
                    if self.large_labels:
                        label1 = int(index_buffer[2 * m - 1])
                        label2 = int(index_buffer[2 * m])
                        i = label1 >> 16
                        j = label1 & 65535
                        k = label2 >> 16
                        l = label2 & 65535
                    else:
                        label = int(index_buffer[m - 1])
                        i = label >> 24
                        j = label >> 16 & 255
                        k = label >>  8 & 255
                        l = label       & 255
                elif int_size == 8:
                    if self.large_labels:
                        label = int(index_buffer[m - 1])
                        i = label >> 48
                        j = label >> 32 & 65535
                        k = label >> 16 & 65535
                        l = label       & 65535
                    else:
                        if m % 2 == 0:
                            label = int(index_buffer[m / 2 - 1])
                            i = label >> 24 & 255
                            j = label >> 16 & 255
                            k = label >>  8 & 255
                            l = label       & 255
                        else:
                            label = int(index_buffer[m / 2])
                            i = label >> 56 & 255
                            j = label >> 48 & 255
                            k = label >> 40 & 255
                            l = label >> 32 & 255
                ints[self.ijkl(i, j, k, l)] = value_buffer[m - 1]
            self.seek(self.tell() + 4)
        return ints

    def read_ci_coeffs(self):
        '''
        Read CI coefficients from file NFT12 and return them as a numpy array of floats.
        '''

        self.seek(4)
        nstates = self.read('i8')
        nconfs = self.read('i8')

        title = self.read('a80', shape=(2,))

        # advance 4 bytes for the end of the record and additional
        # 4 for the start of the new record
        self.seek(self.tell() + 8)
        return self.read('f8', shape=(nconfs * nstates,))


class GamessReader(object):
    '''
    Class for holding method for reading gamess binary files:
        - $JOB.F08 : two electron integrals over AO's,
        - $JOB.F09 : two electron integrals over MO's,
        - $JOB.F15 : GUGA and ORMAS two-electron reduced density matrix,

    TODO:
        CI coefficients, and CI hamiltonian matrix elements.
    '''
    # add an option to choose which reader should be used for dictionary file
    # and sequential files, the options are:
    # - wrapped fortran code that need to be compiled and installed or
    # - native reader written in python using DictionaryFile and SequentialFile
    # classes

    def __init__(self, log):
        self.logfile = log
        i = self.logfile.index("log")
        self.filebase = self.logfile[:i-1]
        self.datfile = self.filebase + ".dat"
        self.twoeaofile = self.filebase + ".F08"
        self.twoemofile = self.filebase + ".F09"
        self.dictionary = self.filebase + ".F10"
        self.rdm2file = self.filebase + ".F15"
        self.gp = GamessLogParser(log=self.logfile)

    def get_onee_size(self, aos=True):
        '''
        Get the size of the vector holding upper (or lower) triangle
        of a square matrix of size naos or nmos.
        '''
        if aos:
            n = self.gp.get_number_of_aos()
        else:
            n = self.gp.get_number_of_mos()
        return n * (n + 1) / 2

    def get_twoe_size(self):
        '''
        Get the size of the 1d vector holding upper (or lower) triangle
        of a supermatrix of size nmos (2RDM and two-electrons integrals).
        '''
        n = self.get_onee_size(aos=False)
        return n * (n + 1) / 2

    def read_rdm2(self, filename=None, nmo=None):

        '''Read the 2rdm from the gamess-us file'''

        # initialize numpy array to zeros
        rdm2 = np.zeros(self.get_twoe_size(), dtype=float)

    def read_twoeao(self, filename=None, nmo=None):

        '''Read the two electron integrals from the gamess-us file'''

        ints = np.zeros(self.get_twoe_size(), dtype=float)


rec = namedtuple('record', ['name', 'dtype'])
records = {
    1: rec("atomic coordinates", "f8"),
    2: rec("enrgys", "f8"),
    3: rec("gradient vector", "f8"),
    4: rec("hessian matrix", "f8"),
    5: rec("not used", ""),
    6: rec("not used", ""),
    7: rec("ptr", "f8"),
    8: rec("dtr", "f8"),
    9: rec("ftr", "f8"),
   10: rec("gtr", "f8"),
   11: rec("bare nucleus", "f8"),
   12: rec("overlap", "f8"),
   13: rec("kinetic energy", "f8"),
   14: rec("alpha fock matrix", "f8"),
   15: rec("alpha orbitals", "f8"),
   16: rec("alpha density matrix", "f8"),
   17: rec("alpha energies or occupation numbers", "f8"),
   18: rec("beta fock matrix", "f8"),
   19: rec("beta orbitals", "f8"),
   20: rec("beta density matrix", "f8"),
   21: rec("beta energies or occupation numbers", "f8"),
   22: rec("error function interpolation table", "f8"),
   23: rec("old alpha fock matrix", "f8"),
   24: rec("older alpha fock matrix", "f8"),
   25: rec("oldest alpha fock matrix", "f8"),
   26: rec("old beta fock matrix", "f8"),
   27: rec("older beta fock matrix", "f8"),
   28: rec("odest beta fock matrix", "f8"),
   29: rec("vib 0 gradient in FORCE", "f8"),
   30: rec("vib 0 alpha orbitals in FORCE", "f8"),
   31: rec("Vib 0 beta  orbitals in FORCE", "f8"),
   32: rec("Vib 0 alpha density matrix in FORCE", "f8"),
   33: rec("Vib 0 beta  density matrix in FORCE", "f8"),
   34: rec("dipole derivative tensor in FORCE", "f8"),
   35: rec("frozen core Fock operator", "f8"),
   36: rec("RHF/UHF/ROHF Lagrangian", "f8"),
   37: rec("floating point part of common block /OPTGRD/", "f8"),
   38: rec("integer part of common block /OPTGRD/", "i8"),
   39: rec("ZMAT of input internal coords", "f8"),
   40: rec("IZMAT of input internal coords", "i8"),
   41: rec("B matrix of redundant internal coords", "f8"),
   42: rec("pristine core Fock matrix in MO basis (see 87)", "f8"),
   43: rec("Force constant matrix in internal coordinates", "f8"),
   44: rec("SALC transformation", "f8"),
   45: rec("symmetry adapted Q matrix", "f8"),
   46: rec("S matrix for symmetry coordinates", "f8"),
   47: rec("ZMAT for symmetry internal coords", "f8"),
   48: rec("IZMAT for symmetry internal coords", "i8"),
   49: rec("B matrix", "f8"),
   50: rec("B inverse matrix", "f8"),
   69: rec("alpha Lowdin populations", "f8"),
   70: rec("beta Lowdin populations", "f8"),
   71: rec("alpha orbitals during localization", "f8"),
   72: rec("betha orbitals during localization", "f8"),
   73: rec("alpha localization transformation", "f8"),
   74: rec("beta localization transformation", "f8"),
   95: rec("x dipole integrals in AO basis", "f8"),
   96: rec("y dipole integrals in AO basis", "f8"),
   97: rec("z dipole integrals in AO basis", "f8"),
  251: rec("static polarizability tensor alpha", "f8"),
  252: rec("X dipole integrals in MO basis", "f8"),
  253: rec("Y dipole integrals in MO basis", "f8"),
  254: rec("Z dipole integrals in MO basis", "f8"),
  255: rec("alpha MO symmetry labels", "S8"),
  256: rec("beta MO symmetry labels", "S8"),
  286: rec("oriented localized molecular orbitals", "f8"),
  379: rec("Lz integrals", "f8"),
}


class DictionaryFile(BinaryFile):
    '''
    Wrapper for reading GAMESS(US) dictionary file (\*.F10).
    '''

    def __init__(self, filename, irecln=4090, int_size=8):
        """

        Args:
          irecln: int
            Rrecord length that is used by GAMESS(US) when
            writing the dicitonary file. It is a platform dependent
            variable that is called IRECLN in the GAMESS(US) code. In
            GAMESS(US) it is set by the function NRASIZ(UNIT) in iolib.src
            and for dictionary file (unit=10) and UNX it is equal to 4090,
            for all other files and UNX it is equal to 2048. If you use
            GAMESS(US) on a different platform check the NRASIZ(UNIT)
            function for the proper value and supply it when instantiating
            the class.
          int_size: int
            Integer size (in bytes) that the GAMESS(US) was compiled with.
        """
        super(DictionaryFile, self).__init__(filename)

        self.irecln = irecln
        self.int_size = int_size
        # read the first record with the information about the
        # structure of the dictionary file
        self.irecst = self.read(np.dtype('i' + str(self.int_size)))
        self.ioda = self.read(np.dtype('i' + str(self.int_size)), shape=(950,))
        self.ifilen = self.read(np.dtype('i' + str(self.int_size)), shape=(950,))
        self.iss = self.read(np.dtype('i' + str(self.int_size)))
        self.ipk = self.read(np.dtype('i' + str(self.int_size)))

    def read_record(self, nrec, dtype=None):
        '''
        Read a logical record 'rec' from the dictionary file and return a numpy
        array of type defined in the 'records' list, and size defined through
        'self.ifilen' array.
        '''

        if self.ioda[nrec-1] < 0:
            raise IOError("Record {0} was not previously written, IODA[{0}]={1}".format(nrec, self.ioda[nrec-1]))

        self.seek(8 * self.irecln * (int(self.ioda[nrec - 1]) - 1))
        if dtype is not None:
            return self.read(dtype, shape=(self.ifilen[nrec - 1],))
        else:
            return self.read(records[nrec].dtype, shape=(self.ifilen[nrec - 1],))


def tri2full(vector, sym=True):
    '''
    Convert a triagonal matrix whose elements are stored in the `vector` into a
    rectangular matrix of the shape given by `shape` tuple.
    '''

    # get the shape of the symmetric matrix from solving n^2 + n - 2s = 0 equation
    # where n is the number of rows/columns and s is the size of 1D vector

    n = int((np.sqrt(8.0 * vector.size + 1) - 1.0) / 2.0)
    matrix = np.zeros((n, n), dtype=float, order='F')

    ij = -1
    for i in range(n):
        for j in range(i + 1):
            ij += 1
            if sym:
                matrix[i, j] = matrix[j, i] = vector[ij]
            else:
                matrix[i, j] = vector[ij]
                matrix[j, i] = -vector[ij]

    return matrix
