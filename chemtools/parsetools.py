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
Module with convenience functions for parsing files
'''

from __future__ import print_function

import itertools
from collections import defaultdict


def contains(string, query):
    'Check if `string` contains `query`'
    return string.find(query) > -1


def locatelinenos(filename, tolocate):
    '''
    Given a file and a list of strings return a dict with string as keys
    and line numbers in which they appear a values.

    Args:
      filename : str
        Name of the file
      tolocate : list of tuples
        List of tuples with strings to find (queries) as first elements and
        integer offset values as second

    Returns:
      out : dict
        Dictionary whose keys are indices corresponding to item in input list
        and values are lists of line numbers in which those string appear

    TODO:
        - add option to ignore the case of the strings to search
    '''

    out = defaultdict(list)
    for lineno, line in enumerate(open(filename, 'r')):
        for idx, (query, offset) in enumerate(tolocate):
            if contains(line, query):
                out[idx].append(lineno + offset)
    return out


def getlines(filename, tolocate):
    '''
    Return the lines from the files based on `tolocate`

    Args:
      filename : str
        Name of the file
      tolocate : list of tuples
        List of tuples with strings to find (queries) as first elements and
        integer offset values as second

    Return:
    '''

    located = locatelinenos(filename, tolocate)

    if len(tolocate) == len(located):
        for k, v in located.items():
            if len(v) > 1:
                raise ValueError('multiple lines found for "{0}": {1}'.format(
                    tolocate[k][0], ', '.join([str(x) for x in v])))

        startlno = min(list(itertools.chain(*located.values())))
        endlno = max(list(itertools.chain(*located.values())))
        return getchunk(filename, startlno, endlno)
    else:
        # TODO: this needs to be corrected to be more informative
        raise ValueError('len(tolocate) != len(located): {0} != {1}'.format(
            len(tolocate), len(located)))


def getchunk(filename, startlno, endlno):
    '''
    Get a list of lines from a file between specified line numbers `startlno`
    and `endlno`.

    Args:
      filename : str
        Name of the file to process
      startlno : int
        Number of the first line to obtain
      endlno : int
        Number of the last line to obtain

    Returns:
      lines : list
        A list of lines from the file `filename` between line numbers\
        `startlno` and `endlno`
    '''

    fobj = open(filename, 'r')
    fileiter = iter(fobj)

    for _ in range(startlno):
        next(fileiter)

    return [next(fileiter) for _ in range(endlno - startlno)]


def take(seq, num):
    '''
    Iterate over a sequence `seq` `num` times and return the list of the
    elements iterated over.
    '''
    return [next(seq) for _ in range(num)]


def parsepairs(los, sep="="):
    '''
    Parse a given list of strings "los" into a dictionary based on
    separation by "sep" character and return the dictionary.
    '''

    out = []
    for line in los:
        if sep in line:
            (name, value) = line.split(sep)
            out.append((name.strip(), float(value)))
    return dict(out)


def sliceafter(seq, item, num):
    '''
    Return "num" elements of a sequence "seq" present after the item "item".
    '''

    it = iter(seq)
    for element in it:
        if item in element:
            return [next(it) for _ in range(num)]


def slicebetween(string, start, end):
    '''
    Return a slice of the `string` between phrases `start` and `end`.
    '''

    istart = string.index(start)
    iend = string[istart:].index(end)
    return string[istart + len(start):istart + iend]
