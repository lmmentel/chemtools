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
An abstract class that should be subclassed when addding a new code interface.

All the methods should be implemented for the basisopt module to work with the
new code.
'''

import os
import re

from string import Template

from abc import ABCMeta, abstractmethod


class Calculator():
    '''
    Abstract class that should be subclassed when adding a new code interace.
    '''
    __metaclass__ = ABCMeta

    def __init__(self, exevar=None, executable=None, runopts=None,
                 scratch=None):

        self.exevar = exevar

        if executable is None and exevar is None:
            self.exe_found = False
        elif executable is None and exevar is not None:
            self.executable = os.getenv(self.exevar)
            self.exe_found = True
        elif executable is not None and exevar is None:
            self.executable = executable
            self.exe_found = True
        elif executable is not None and exevar is not None:
            exefromvar = os.getenv(self.exevar)
            if os.path.abspath(exefromvar) == os.path.abspath(executable):
                self.executable = exefromvar
                self.exe_found = True
            else:
                raise ValueError('conflicting <executable> and <exevar> options specified: '\
                                 '{0:s} != {1:s}'.format(executable, exefromvar))

        self.runopts = runopts
        self.scratch = scratch

    @property
    def executable(self):
        'Set the executable'
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
        'Return scratch'
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
    def parse(self, fname, objective, regexp=None):
        '''
        Parse the value of the ``objective`` from the file ``fname``
        '''
        pass

    @abstractmethod
    def run(self, *args):
        '''
        run a single job
        '''

        pass

    @abstractmethod
    def run_multiple(self, *args):
        '''
        Run multiple jobs
        '''
        pass

    def accomplished(self, fname):
        '''
        Return True if the job completed without errors
        '''

        return self.parse(fname, 'accomplished') is True

    @abstractmethod
    def write_input(self):
        '''
        write the input file in the format of the code used
        '''
        pass


def parse_objective(fname, regexp, reflags=re.I):
    '''
    Wrapper for the parse methods

    Args:
        fname : :py:class:`str`
            Name of the file to be parsed
        regexp : :py:class:`str`
            Regular expression as raw string
        reflags :
            Flags for the regular expression compilation
    '''

    if not os.path.exists(fname):
        raise ValueError("Output file: {0:s} doesn't exist in {1:s}".format(fname, os.getcwd()))

    rec = re.compile(regexp, flags=reflags)
    with open(fname, 'r') as fobj:
        for line in fobj:
            match = rec.search(line)
            if match:
                # check if a group was captured or not
                if match.lastindex is not None and match.lastindex > 0:
                    return float(match.group(1))
                else:
                    return True
        else:
            return None


class InputTemplate(Template):

    'Modified string.Template class to be used for input rendering'

    delimiter = '%'
    idpattern = r'[a-z][_a-z0-9]*'

    def get_keys(self):
        '''
        Parse the string and return a dict with possible keys to substitute.
        For most use case only the `named` fields are interesting.
        '''

        keys = {}
        match = self.pattern.findall(self.template)
        for key, val in self.pattern.groupindex.items():
            keys[key] = [x[val - 1] for x in match if x[val - 1] != '']
        return keys


