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
                    raise OSError('{0} is a directory, not executable'.format(value))
            elif os.path.isfile(value):
                if os.access(value, os.X_OK):
                    self._executable = value
                else:
                    raise OSError('file: {0} is not executable'.format(value))
            else:
                raise OSError('unknown type of the executable: {0}'.format(value))
        else:
            raise ValueError("executable: '{0}' does not exists".format(value))

    @property
    def scratch(self):
        return self._scratch

    @scratch.setter
    def scratch(self, path):

        if path is None:
            self._scratch = os.path.join(os.getenv('HOME'), 'scratch')
        elif os.path.exists(path):
            self._scratch = path
        else:
            raise ValueError("Scratch directory doesn't exist: {}".format(path))

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
