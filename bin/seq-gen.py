#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
from json.decoder import JSONDecodeError
from seqgen import Sequences

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=('Create genetic sequences according to a '
                 'JSON specification file and write them to stdout.'))

parser.add_argument(
    '--specification', metavar='FILENAME', default=sys.stdin, type=open,
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

args = parser.parse_args()

try:
    sequences = Sequences(args.specification,
                          defaultLength=args.defaultLength,
                          defaultIdPrefix=args.defaultIdPrefix)
except JSONDecodeError:
    print('Could not parse your specification JSON. Stacktrace:',
          file=sys.stderr)
    raise
else:
    for sequence in sequences:
        print(sequence.toString('fasta'), end='')
