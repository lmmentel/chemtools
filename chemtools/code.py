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

    def __init__(self, name=None, executable=None, execpath=None, runopts=None, scratch=None):
        self.name = name
        self.execpath = execpath
        self.executable = executable
        self.runopts = runopts
        self.scratch = scratch

    @property
    def execpath(self):
        return self._execpath

    @execpath.setter
    def execpath(self, value):
        if os.path.exists(value):
            self._execpath = value
        else:
            raise ValueError("execpath: '{}' does not exists".format(value))

    @property
    def executable(self):
        return self._executable

    @executable.setter
    def executable(self, value):
        if os.path.exists(os.path.join(self.execpath, value)):
            self._executable = value
        else:
            raise ValueError("executable: '{0}' does not exists at {1}".format(value, self.execpath))

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
