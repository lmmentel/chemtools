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

from collections import defaultdict

def contains(string, query):
    'Check if `string` contains `query`'
    return string.find(query) > -1

def locate_lines(filename, strings):
    '''
    Given a file and a list of strings return a dict with string as keys
    and line numbers in which they appear a values.

    Args:
      filename : str
        Name of the file
      string : dict
        Dictionary with labels as keys and strings to find (queries) as values

    Returns:
      out : dict
        Dictionary whose keys are passed string and values are lists of line
        numbers in which those string appear

    TODO:
        - add option to ignore the case of the strings to search
        - add option to specify offset
    '''

    out = defaultdict(list)
    for lineno, line in enumerate(open(filename, 'r')):
        for label, query in strings.items():
            if contains(line, query):
                out[label].append(lineno)
    return out

def get_chunk(filename, start, end):
    '''
    Get a list of lines from a file between specified line numbers `start` and
    `end`.

    Args:
      filename : str
        Name of the file to process
      start : int
        Number of the first line to obtain
      end : int
        Number of the last line to obtain

    Returns:
      lines : list
        A list of lines from the file `filename` between line numbers `start`
        and `end`
    '''

    fobj = open(filename, 'r')
    fileiter = iter(fobj)

    for _ in range(start):
        next(fileiter)

    return [next(fileiter) for _ in range(end - start)]

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
    iend = string.index(end)
    return string[istart+len(start):iend]

