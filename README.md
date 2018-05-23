This repo contains a Python script, [bin/seq-gen.py](bin/seq-gen.py), to
generate genetic sequences based on directives given in a JSON file.
Python (2.7, 3.5, 3.6 are all known to work).

I wrote this because I wanted to be able to easily generate alignments of
FASTA files with known properties, and to then examine the output of
phylogenetic inference or recombination analysis programs run on the
sequences.

## Installation

```sh
$ pip install seqgen
```

You can also get the source
[from PyPI](https://pypi.org/project/seqgen/) or
[on Github](https://github.com/acorg/seqgen).

## Usage

In summary: put a specification for a set of sequences into a
[JSON](http://json.org/)-formatted file. Pass this to `bin/seq-gen.py`,
either using the `--specification` option or on standard input.

Run with `--help` or `-h` to see all options:

```sh
$ seq-gen.py --help
usage: seq-gen.py [-h] [--specification FILENAME] [--defaultIdPrefix PREFIX]
                  [--defaultLength N]

Create genetic sequences according to a JSON specification file and write them
to stdout.

optional arguments:
  -h, --help            show this help message and exit
  --specification FILENAME
                        The name of the JSON sequence specification file.
                        Standard input will be read if no file name is given.
                        (default: stdin)
  --defaultIdPrefix PREFIX
                        The default prefix that sequence ids should have (for
                        those that are not named individually in the
                        specification file) in the resulting FASTA. Numbers
                        will be appended to this value. (default: seq-id-)
  --defaultLength N     The default length that sequences should have (for
                        those that do not have their length given in the
                        specification file) in the resulting FASTA. (default:
                        100)
```

## Sequence specification

Your JSON specifies what sequences you want created.

As an example, the JSON shown below specifies the following:

* Create a random sequence `A` of length 100 nucelotides.
* Make 19 (approximately) 1% mutants from `A`.
* Make a new sequence `B` approximately 15% different from `A`.
* Make 19 (approximately) 1% mutants from `B`.
* Make a recombinant sequence whose first 30 nucleotides are from `A`
  and last 70 nucleotides from `B`.

```json
{
    "variables": {
        "length": 100
    },
    "sequences": [
        {
            "id": "seq-A",
            "length": "%(length)d"
        },
        {
            "from id": "seq-A",
            "id prefix": "seq-A-mutant-",
            "length": "%(length)d",
            "count": 19,
            "mutation rate": 0.01
        },
        {
            "id": "seq-B",
            "from id": "A",
            "length": "%(length)d",
            "mutation rate": 0.15
        },
        {
            "from id": "seq-B",
            "id prefix": "seq-B-mutant-",
            "length": "%(length)d",
            "count": 19,
            "mutation rate": 0.01
        },
        {
            "id": "recombinant",
            "sections": [
                {
                    "from id": "seq-A",
                    "start": 1,
                    "length": 30
                },
                {
                    "from id": "seq-B",
                    "start": 31,
                    "length": 70
                }
            ]
        }
    ]
}
```

Note: in JSON parlance, "object" is what Python programmers call a
"dict". In the text below I am describing the JSON specification, so I'm
using "object". If you're more comfortable thinking of the JSON as having
"dicts", go right ahead.

<a id="convenience"></a>
The `variables` section is optional. For convenience, if no `variables`
section is given, the JSON may simply be a `list` of sequence objects.

There are many more small examples of usage in
[test/testSequences.py](test/testSequences.py).

The sequence specification must be a list of objects.  The full list of
specification keys you can put in a object in the JSON is as follows:

* `alphabet`: A string of characters from which sequences will be drawn
  (see `random nt` and `random aa` below).
* `count`: The number of sequences to generate from this object.
* `description`: The sequence description. This will be appended to the
  FASTA id (separated by a space). Note that if a sequence has a
  description and you want to refer to it using `from id` (see below), you
  must include its id and description, separated by a single space.
* `id`: The FASTA id to give this sequence.
* `id prefix`: The prefix of the FASTA ids to give a set of sequences. A
  count will be appended to this prefix. This is useful when you specify a
  `count` value.
* `from id`: The sequence should be based on another (already named)
  sequence in the JSON file. The value given should either be the exact
  `id` of another sequence or else be the `id prefix` of another sequence
  followed by a number. E.g., if you specify `"count": 10` and `"id
  prefix": "my-id"` for a set of ten sequences, and want to refer to the
  3rd of them when specifying another sequence, you would use `"from id":
  "my-id-3"`. The default id prefixes for sequences that are not given an
  id explicity is `seq-id-`.
* `length`: The sequence length.
* `mutation rate`: A mutation rate to apply to the sequence.
* `random aa`: The sequence should be made of random amino acids.
* `random nt`: The sequence should be made of random nucleotides. This is
  the default.
* `ratchet`: If `true` and a count and mutation rate are given, sequences
  after the first will be mutants of the previous sequence. This can be
  used to build a series of sequences that are successive mutants of one
  another. For the name derivation, see
  [Muller's Ratchet](https://en.wikipedia.org/wiki/Muller's_ratchet).
* `sections`: Gives a list of sequence sections used to build up another
  sequence (see example above). Each section is a JSON object that may
  contain the keys present in a regular sequence object, excluding `count`,
  `description`, `id`, `id prefix`, `name`, and `sections`. The main idea
  here is to allow you to easily specify sequences that are recombinants.
* `sequence`: Give the exact sequence of nucleotides or amino acids.
* `sequence file`: Specify a FASTA file to get the sequence from. Currently
  only the first sequence in the file is used.
* `skip`: If `true` the sequence will not be output. This is useful either
  for temporarily omitting a sequence or for just giving a sequence (e.g.,
  one read from a file) an id so it can be used in the construction of
  another sequence.

All specification keys are optional. A completely empty specification
object will get you a sequence of the default length, with a default id,
composed of random nucleotides. Hence:

```sh
$ echo '[{},{}]' | seq-gen.py
>seq-id-1
ACTCGTGCTATAGGGCGAATATCGCAAAATGCTCACATACCCAATAGCTTAGGAATAGTTCCTGTCGGGGCGCTCGTTGATTTAAGTCAATGAGCATCCT
>seq-id-2
CTTGAGATGATTCGGCAACGTTAGCCGATAGATCATGGAAGGAATACGGCTAAAATATTCAGGTAATTAATGGATACGTCCTAGATAAGTAGAATCGAAT
```

(Note that this example takes advantage of the convenience <a
href="#convenience">mentioned above</a>).

Although the code will complain about unknown keys, it does not detect
cases where you specify a sequence in two different ways. You'll have to
play around and/or read the code in
[seqgen/sequences.py](seqgen/sequences.py) to see the order in which the
various sequence specification keys are acted on.

## Development

To run the tests:

```sh
$ make check
```

or if you have [Twisted](https://twistedmatrix.com/trac/) installed, you
can use its `trial` test runner, via

```sh
$ make tcheck
```

You can also use

```sh
$ tox
```

to run tests for various versions of Python.
