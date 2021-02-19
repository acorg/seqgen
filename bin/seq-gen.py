#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
from os import path
from json.decoder import JSONDecodeError
from seqgen import Sequences


def pngFile(filepath):
    _, file_extension = path.splitext(filepath)
    if file_extension != '.png':
        parser.error(f'File {filepath!r} needs to have a .png extension.')
    return filepath


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=('Create genetic sequences according to a '
                 'JSON specification file and write them to stdout.'))

parser.add_argument(
    '--specification', metavar='FILENAME', default=sys.stdin, type=open,
    required=True,
    help=('The name of the JSON sequence specification file. Standard input '
          'will be read if no file name is given.'))

parser.add_argument(
    '--defaultIdPrefix', metavar='PREFIX', default=Sequences.DEFAULT_ID_PREFIX,
    help=('The default prefix that sequence ids should have (for those that '
          'are not named individually in the specification file) in the '
          'resulting FASTA. Numbers will be appended to this value.'))

parser.add_argument(
    '--defaultLength', metavar='N', default=Sequences.DEFAULT_LENGTH, type=int,
    help=('The default length that sequences should have (for those that do '
          'not have their length given in the specification file) in the '
          'resulting FASTA.'))

parser.add_argument(
    '--makeTree', action='store_true', default=False,
    help=("Make a rooted tree based on the given specification and display it "
          "in a .png file. Providing a filename with a .png extension via "
          "'--treeFile' is necessary. Recombinant sequences will not be used "
          "as they cannot be properly represented in a Newick tree."))

parser.add_argument(
    '--treeFile', type=pngFile,
    help=("The path to a .png file where the tree implied by the "
          "specification file will be stored. Option '--makeTree' must be "
          "set."))

parser.add_argument(
    '--noSequences', action='store_false', dest='makeSequences', default=True,
    help=("Don't create sequences. You will want to set this option when you "
          "only need to make a Newick tree (see --makeTree and -treeFile "
          "options)."))

args = parser.parse_args()

if bool(args.makeTree) != bool(args.treeFile):
    parser.error('Either both "--makeTree" and "--treeFile" or neither must '
                 'be specified.')

try:
    sequences = Sequences(args.specification,
                          defaultLength=args.defaultLength,
                          defaultIdPrefix=args.defaultIdPrefix)
except JSONDecodeError:
    print('Could not parse your specification JSON. Stacktrace:',
          file=sys.stderr)
    raise
else:
    if args.makeSequences:
        for sequence in sequences:
            print(sequence.toString('fasta'), end='')
    if args.makeTree:
        sequences.buildNewickTree(args.treeFile)
