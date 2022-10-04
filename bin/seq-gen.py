#!/usr/bin/env python

import sys
import argparse
from json.decoder import JSONDecodeError
from seqgen import Sequences

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=(
        "Create genetic sequences according to a "
        "JSON specification file and write them to stdout."
    ),
)

parser.add_argument(
    "--specification",
    metavar="FILENAME",
    default=sys.stdin,
    type=open,
    help=(
        "The name of the JSON sequence specification file. Standard input "
        "will be read if no file name is given."
    ),
)

parser.add_argument(
    "--defaultIdPrefix",
    metavar="PREFIX",
    default=Sequences.DEFAULT_ID_PREFIX,
    help=(
        "The default prefix that sequence ids should have (for those that "
        "are not named individually in the specification file) in the "
        "resulting FASTA. Numbers will be appended to this value."
    ),
)

parser.add_argument(
    "--format",
    choices=("fasta", "fastq"),
    default="fasta",
    help=(
        f"Set the default output format. The output format can be set (via the "
        f"specification file) for each set of reads, if desired. This option just "
        f"sets the default. If the format is 'fastq', the quality for each "
        f"nucleotide will be a constant, according to the value given to --quality "
        f"(or {Sequences.DEFAULT_QUALITY} if --quality is not used."
    ),
)

parser.add_argument(
    "--quality",
    metavar="N",
    help=(
        f"The quality value to use. This will result in FASTQ output. The value will "
        f"be converted to a single character, according to the Phred scale. The "
        f"numeric value you give will be added to the value for '!' to get the "
        f"character that will be used for all quality scores. So use 0 for the lowest "
        f"quality or, e.g., 30 for a reasonably high quality. If --fastq is used but "
        f"--quality is not, a value of {Sequences.DEFAULT_QUALITY} will be used."
    ),
)

parser.add_argument(
    "--defaultLength",
    metavar="N",
    default=Sequences.DEFAULT_LENGTH,
    type=int,
    help=(
        "The default length that sequences should have (for those that do "
        "not have their length given in the specification file) in the "
        "resulting FASTA."
    ),
)

args = parser.parse_args()

try:
    sequences = Sequences(
        args.specification,
        defaultLength=args.defaultLength,
        defaultIdPrefix=args.defaultIdPrefix,
        defaultQuality=args.quality,
        _format=args.format,
    )
except JSONDecodeError:
    print("Could not parse your specification JSON. Stacktrace:", file=sys.stderr)
    raise
else:
    sequences.write()
