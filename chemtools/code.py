'''
An abstract class that should be subclassed when addding a new code interface.

All the methods should be implemented for the basisopt module to work with the
new code.
'''

from abc import ABCMeta, abstractmethod
import os

class Code():
    __metaclass__ = ABCMeta
    '''
    Abstract class that should be subclassed when adding a new code interace.
    '''

    def __init__(self, executable=None, runopts=None, scratch=None):

        self.executable = executable
        self.runopts = runopts
        self.scratch = scratch

    @property
    def executable(self):
        return self._executable

    @executable.setter
    def executable(self, value):
        if os.path.exists(value):
            if os.path.isdir(value):
                    raise OSError('{0} is a directory not executable'.format(value))
            elif os.path.isfile(value):
                if os.access(value, os.X_OK):
                    self._executable = value
                else:
                    raise OSError('file: {0} is not executable'.format(value))
            else:
                raise OSError('unknown type of the executable: {0}'.format(value))
        else:
            raise ValueError("executable: '{0}' does not exists".format(value))

    @abstractmethod
    def parse():
        pass

    @abstractmethod
    def run():
        '''
        run s single job
        '''
        pass

    @abstractmethod
    def run_multiple():
        pass

    @abstractmethod
    def accomplished():
        '''
        return True if the job completed without errors
        '''
        pass

    @abstractmethod
    def write_input():
        '''
        write the input file in the format of the code used
        '''
        pass
