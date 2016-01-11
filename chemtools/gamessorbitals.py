'Orbitals class'

import numpy as np
import pandas as pd

from chemtools.gamessus import GamessLogParser
from chemtools.gamessreader import DictionaryFile, tri2full

class Orbitals(pd.DataFrame):
    '''
    A convenience class for handling GAMESS(US) orbitals.
    '''

    def __init__(self, *args, **kwargs):

        super(Orbitals, self).__init__(*args, **kwargs)

    @classmethod
    def from_files(cls, name=None, logfile=None, dictfile=None):
        '''
        Initialize the `Orbitals` instance based on orbital information parsed from the `logfile`
        and read from the `dictfile`.

        Args:
            name : str
                One of `hf` or `ci`
            logfile : str
                Name of the GAMESS(US) log file
            dictfile : str
                Name of the GAMESS(US) dictionary file .F10
        '''

        if name == 'hf':
            evec_record = 15
            eval_record = 17
            syml_record = 255
        elif name == 'ci':
            evec_record = 19
            eval_record = 21
            syml_record = 256
        else:
            raise ValueError('name should be one either "hf" or "ci"')

        # parse the number of aos and mos
        glp = GamessLogParser(logfile)
        nao = glp.get_number_of_aos()
        nmo = glp.get_number_of_mos()

        # read the relevant record from the dictfile
        df = DictionaryFile(dictfile)
        # read the orbitals
        mosv = df.read_record(evec_record)
        mos = mosv[:nao*nmo].reshape((nao, nmo), order='F')
        # read the eigenvectors
        ev = df.read_record(eval_record)
        ev = ev[:nmo]
        # read symmetry labels
        symlab = df.read_record(syml_record)
        symlab = symlab[:nmo]

        data = dict([
                ('symlabels', symlab),
                ('eigenvals', ev),
                ('gindex', range(1, nmo + 1)),
        ])

        dataframe = cls(data=data)

        # remove whitespace from symmetry labels
        dataframe['symlabels'] = dataframe['symlabels'].str.strip()

        dataframe.nao = nao
        dataframe.nmo = nmo
        dataframe.name = name
        dataframe.logfile = logfile
        dataframe.dictfile = dictfile
        dataframe.coeffs = pd.DataFrame(data=mos)

        return dataframe

    def assign_lz_values(self, decimals=6, tolv=1.0e-2):
        '''
        Determine the eigenvalues of the Lz operator for each nondegenerate or
        a combination of degenerate orbitals and assign them to `lzvals` column
        in the DataFrame

        The Lz integrals over atomic orbitals are read from the dictionary file
        record No. 379.

        Args:
          decimals : int
            Number of decimals to keep when comparing float eigenvalues
          tolv : float
            Threshold for keeping wights of the eiegenvalues (squared eigenvector components)

        WARNING:
            currently this only works if the eigenvalues on input are sorted
        '''

        lzmap = {0: 'sigma', 1: 'pi', 2: 'delta', 3: 'phi'}

        df = DictionaryFile(self.dictfile)
        # read Lz integrals from the dictionary file
        lzvec = df.read_record(379)
        # expand triangular to square matrix
        lzints = tri2full(lzvec, sym=False) 

        unq, indices = check_duplicates(self.eigenvals, decimals=decimals)

        out = []

        mos = self.coeffs.values

        for IG, (i, j) in enumerate(indices, start=1):
            # transform the integrals to the MO basis
            lzmo = np.dot(mos[:, i:j].T, np.dot(lzints, mos[:, i:j]))

            # tranform to complex type
            lzcplx = np.zeros_like(lzmo, dtype=np.complex128)
            lzcplx.imag = -lzmo

            # diagonalize the lzcplx matrix
            w, v = np.linalg.eig(lzcplx)

            wij2 = np.multiply(np.abs(v), np.abs(v))
            #x, y = np.where(wij2 > tolv)
            idxs = [np.nonzero(row > tolv)[0] for row in wij2]
            for row, comp in enumerate(idxs):
                out.append([(IG, w.real[m], wij2[row, m]*100.0) for m in comp])

        evs = [np.unique(np.round(np.abs(np.array([x[1] for x in row])), decimals=1)) for row in out]
        lzvals = [int(x[0]) if x.size == 1 else np.NaN for x in evs]

        self['lzvals'] = [int(x[0]) if x.size == 1 else np.NaN for x in evs]
        self['lzlabels'] = self['lzvals'].map(lzmap)
        return out

    def lzmos(self, ETOLLZ=1.0e-6, TOLV=1.0e-2):
        '''
        A python rewrite of the GAMESS(US) LZMOS subroutine (from symorb.src) for analyzing the Lz composition of the
        orbtials

        Args:
            ETOLLZ : float
            TOLV : float

        '''

        lzmap = {0: 'sigma', 1: 'pi', 2: 'delta', 3: 'phi'}

        df = DictionaryFile(self.dictfile)
        # read Lz integrals from the dictionary file
        lzvec = df.read_record(379)
        # expand triangular to square matrix
        lzints = tri2full(lzvec, sym=False)

        mos = self.coeffs.values
        ev = self.eigenvals.values
        out = []

        L = 0
        ILAST = 0
        IG = 0
        N = 1
        for I in range(1, self.nmo + 1):
            EI = 0
            if I < self.nmo: EI = ev[I]

            if (abs(ev[I - 1] - EI) > ETOLLZ or I > self.nmo):

                IG += 1
                # transform the integrals to the MO basis
                lzmo = np.dot(mos[:, ILAST:ILAST+N].T, np.dot(lzints, mos[:, ILAST:ILAST+N]))
                # tranform to complex type
                lzcplx = np.zeros_like(lzmo, dtype=np.complex128)
                lzcplx.imag = -lzmo

                # diagonalize the lzcplx matrix
                w, v = np.linalg.eig(lzcplx)

                for j in range(v.shape[0]):
                    m = -1
                    temp = []
                    for k in range(v.shape[1]):
                        wij = np.abs(v[j, k])
                        if wij*wij > TOLV:
                            m += 1
                            temp.append(wij*wij*100.0)
                    out.append([(IG, w.real[mi], temp[mi]) for mi in range(m+1)])

                ILAST = I
                N = 1
            else:
                N = N + 1

        evs = [np.unique(np.round(np.abs(np.array([x[1] for x in row])), decimals=1)) for row in out]
        lzvals = [int(x[0]) if x.size == 1 else np.NaN for x in evs]
        # should check here if the length of the evs is the same as the number of rows of the dataframe

        self['lzvals'] = [int(x[0]) if x.size == 1 else np.NaN for x in evs]
        self['lzlabels'] = self['lzvals'].map(lzmap)
        self['ig'] = [x[0][0] for x in out]
        return None


def check_duplicates(a, decimals=6):
    '''
    This funciton assumes that the array `a` is sorted

    http://stackoverflow.com/questions/25264798/checking-for-and-indexing-non-unique-duplicate-values-in-a-numpy-array
    http://stackoverflow.com/questions/5426908/find-unique-elements-of-floating-point-array-in-numpy-with-comparison-using-a-d
    '''

    unq, unq_idx, unq_cnt = np.unique(np.round(a, decimals=decimals), return_inverse=True, return_counts=True)

    indices = np.zeros((len(unq), 2), dtype=np.int)
    for ig in range(len(unq)):
        idx = np.nonzero(unq_idx == ig)[0]
        indices[ig, :] = [np.min(idx), np.max(idx) + 1]

    return unq, indices
