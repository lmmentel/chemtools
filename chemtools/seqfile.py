"""
BinaryFile: A class for accessing data to/from large binary files
=================================================================

The data is meant to be read/write sequentially from/to a binary file.
One can request to read a piece of data with a specific type and shape
from it.  Also, it supports the notion of Fortran and C ordered data,
so that the returned data is always well-behaved (C-contiguous and
aligned).

This is also seeking capable.

:Author:   Francesc Alted
:Contact:  faltet@pytables.org
:Created:  2010-03-18
:Acknowledgment: Funding for the development of this code is provided
through the Norwegian Research Council VAUUAV project #184724, 2010

"""

import docopt
import numpy as np

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


def readseq(filename, buffSize=15000, intSize=8, largeLabels=False, skipFirst=False):
    '''

    Input
    =====
        filename (str)
            name of the file to read,
        buffSize (int)
            size of the buffer holding values to be read, in gamess(us) it is
            stored under "NINTMX" variable and in Linux version is equal to
            15000 which is the default value,
        largeLabels (bool)
            a flag indicating if large labels should were used, if largest
            label "i" (of the MO) is i<255 then largeLabels should be False
            (case "LABSIZ=1" in gamess(us)), otherwise set to True (case
            "LABSIZ=2" in gamess(us),
        skipFirst (bool)
            skips the first record of the file is set to True,

    Output
    ======
        numpy 1D array holding the values
    '''

    if intSize == 4:
        if largeLabels:
            indexBuffSize = 2*buffSize
        else:
            indexBuffSize = buffSize
    elif intSize == 8:
        if largeLabels:
            indexBuffSize = buffSize
        else:
            indexBuffSize = (buffSize + 1)/2
    else:
        raise ValueError

    intType = np.dtype('i'+str(intSize))
    indexBuffer = np.zeros(indexBuffSize, dtype=intType, order='F')
    valueBuffer = np.zeros(buffSize, dtype=float, order='F')

    seqf = BinaryFile(filename)

    if skipFirst:
        seqf.seek(92)  # this works for now but should be tested on other platforms and for other values of intsize
    n = seqf.read(intType)
    print n
    print seqf.tell()
    indexBuffer = seqf.read(intType, shape=(indexBuffSize))
    print indexBuffer[:31]
    print seqf.tell()
    valueBuffer = seqf.read('f8', shape=(buffSize))
    print valueBuffer[:32]

    for i in range(abs(n)):
        if intSize == 4:
            if largeLabels:
                pass
            else:
                pass
        elif intSize == 8:
            if largeLabels:
                pass
            else:
                pass


def main():
    '''
    read gamess(us) fortran sequential unformatted file

    Usage:
        seqfile.py <file>
    '''

    args = docopt.docopt(main.__doc__, help=True)

    readseq(args['<file>'], skipFirst=True)

if __name__ == "__main__":
    main()
