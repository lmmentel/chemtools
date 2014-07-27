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

    def __init__(self, name=None, execpath=None, runopts=None, scratch=None, **kwargs):
        self.name = name
        self.execpath = execpath
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
            raise ValueError("File: '{}' does not exists".format(value))

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
