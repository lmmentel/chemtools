# -*- coding: utf-8 -*-

from argparse import ArgumentParser
import os

from chemtools.basisset import BasisSet


def bsprint():
    '''
    CLI script to read a basis set from a file in internal format and print
    in a specified format
    '''

    parser = ArgumentParser()
    parser.add_argument("file",
                        help="file name with a pickled BasisSet object")
    parser.add_argument("-f", "--format",
                        choices=["cfour", "dalton", "gamessus", "gaussian",
                                 "json", "molpro", "nwchem"],
                        default="default")
    args = parser.parse_args()

    bs = BasisSet.from_pickle(args.file)

    if args.format == 'default':
        print(str(bs))
    else:
        writer = "to_" + args.format
        method = getattr(bs, writer)
        print(method())


def bsconvert():
    '''
    CLI script to convert between different basis set formats
    '''

    parser = ArgumentParser(description='Convert basis set between formats of different programs')
    parser.add_argument("filename",
                        help="file name with a basis set, default='pickle'")
    parser.add_argument("-if",
                        "--inputformat",
                        choices=["gamessus", "gaussian", "json", "molpro",
                                 "pickle"],
                        help="Basis set input format",
                        default="pickle")
    parser.add_argument("-of",
                        "--outputformat",
                        choices=["cfour", "dalton", "gamessus", "gaussian",
                                 "json", "molpro", "nwchem", "pickle"],
                        help="Basis set output format",
                        default="molpro")
    args = parser.parse_args()

    name = os.path.splitext(args.filename)[0]

    if args.inputformat == "pickle":
        bsets = BasisSet.from_pickle(args.filename)
    elif args.inputformat == "json":
        bsets = BasisSet.from_json(args.filename)
    else:
        bsets = BasisSet.from_file(fname=args.filename, fmt=args.inputformat,
                                   name=name)

    if args.outputformat == "pickle":
        if isinstance(bsets, dict):
            for elem, bset in bsets.items():
                bset.to_pickle(name + '_' + elem + '.pkl')
        elif isinstance(bsets, BasisSet):
            bsets.to_pickle(name + '_' + bsets.element + '.pkl')
    else:
        writer = "to_" + args.outputformat

        if isinstance(bsets, dict):
            for elem, bset in bsets.items():
                method = getattr(bset, writer)
                print(method())
        elif isinstance(bsets, BasisSet):
            method = getattr(bsets, writer)
            print(method())
        else:
            raise ValueError('Something went wrong')
