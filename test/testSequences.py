from unittest import TestCase
from six.moves import builtins
from six import assertRaisesRegex, PY3, StringIO
from seqgen.sequences import Sequences
from dark.aa import AA_LETTERS

try:
    from unittest.mock import mock_open, patch
except ImportError:
    from mock import mock_open, patch

open_ = ('builtins' if PY3 else '__builtin__') + '.open'


class TestSequences(TestCase):
    """
    Test the Sequences class.
    """
    @patch(open_, new_callable=mock_open)
    def testNonExistentSpecificationFile(self, mock):
        """
        Passing a specification filename that does not exist must raise a
        FileNotFoundError (PY3) or IOError (PY2).
        """
        errorClass = builtins.FileNotFoundError if PY3 else IOError
        mock.side_effect = errorClass('abc')
        error = '^abc$'
        assertRaisesRegex(self, errorClass, error, Sequences, spec='filename')

    @patch(open_, new_callable=mock_open, read_data="not JSON")
    def testNotJSON(self, mock):
        """
        If the specification is not valid JSON, a ValueError must be raised.
        """
        if PY3:
            error = '^Expecting value: line 1 column 1 \(char 0\)$'
        else:
            error = '^No JSON object could be decoded$'
        assertRaisesRegex(self, ValueError, error, Sequences,
                          spec='filename')

    @patch(open_, new_callable=mock_open, read_data='{"xxx": 33}')
    def testNoSequencesKey(self, mock):
        """
        If the JSON specification has no 'sequences' key, a ValueError
        must be raised.
        """
        error = "^The specification JSON must have a 'sequences' key\.$"
        assertRaisesRegex(self, ValueError, error, Sequences,
                          spec='filename')

    @patch(open_, new_callable=mock_open,
           read_data='[{"name": "a"}, {"name": "a"}]')
    def testDuplicatedName(self, mock):
        """
        If a duplicate sequence name is present in the JSON, a ValueError
        must be raised.
        """
        s = Sequences(spec='filename')
        error = "^Name 'a' is duplicated in the JSON specification\.$"
        assertRaisesRegex(self, ValueError, error, list, s)

    def testNoSequences(self):
        """
        If no sequences are specified, none should be returned.
        """
        s = Sequences(StringIO('[]'))
        self.assertEqual([], list(s))

    def testOneSequenceNameOnly(self):
        """
        If only one sequence is specified, and only by name, one sequence
        should be created, it should have the default length, the expected
        id, and it should be entirely composed of nucleotides.
        """
        s = Sequences(StringIO('[{"name": "a"}]'))
        (read,) = list(s)
        self.assertEqual(Sequences.DEFAULT_ID_PREFIX + '1', read.id)
        self.assertEqual(Sequences.DEFAULT_LENGTH, len(read.sequence))
        self.assertEqual(set(), set(read.sequence) - set('ACGT'))

    def testOneSequenceIdOnly(self):
        """
        If only one sequence is specified, and only by id, one sequence
        should be created, it should have the default length, the expected
        id, and it should be entirely composed of nucleotides.
        """
        s = Sequences(StringIO('[{"id": "the-id"}]'))
        (read,) = list(s)
        self.assertEqual('the-id', read.id)
        self.assertEqual(Sequences.DEFAULT_LENGTH, len(read.sequence))
        self.assertEqual(set(), set(read.sequence) - set('ACGT'))

    def testOneSequenceSequenceOnly(self):
        """
        If only one sequence is specified, and only by its sequence, one
        sequence should be created, and it should have the specified sequence.
        """
        s = Sequences(StringIO('[{"sequence": "ACCG"}]'))
        (read,) = list(s)
        self.assertEqual('ACCG', read.sequence)

    # The following doesn't work under Python 3. Grrr.
    @patch(open_, new_callable=mock_open, read_data='>id1\nACCT\n')
    def xxx_testOneSequenceSequenceFileOnly(self, mock):
        """
        If only one sequence is specified, and only by its sequence filename,
        one sequence should be read and created, and it should have the
        specified sequence.
        """
        s = Sequences(StringIO('[{"sequence file": "xxx.fasta"}]'))
        (read,) = list(s)
        self.assertEqual('ACCG', read.sequence)

    @patch(open_, new_callable=mock_open)
    def xxx_testOneSequenceSequenceFileOnlyUnknownFile(self, mock):
        """
        If only one sequence is specified, and only by its sequence filename,
        but the file is unknown, ValueError must be raised.
        """
        s = Sequences(StringIO('[{"sequence file": "xxx.fasta"}]'))
        errorClass = builtins.FileNotFoundError if PY3 else IOError
        mock.side_effect = errorClass('abc')
        error = '^abc$'
        assertRaisesRegex(self, errorClass, error, list, s)

    def testOneSequenceIdAndCountGreaterThanOne(self):
        """
        If only one sequence is specified with an id, a ValueError must be
        raised if its count is greater than one.
        """
        s = Sequences(StringIO('''{
            "sequences": [
                {
                    "id": "the-id",
                    "count": 6
                }
            ]
        }'''))
        error = ("^Sequence with id 'the-id' has a count of 6\. If you want "
                 "to specify one sequence with an id, the count must be 1\. "
                 "To specify multiple sequences with an id prefix, use "
                 "'id prefix'\.$")
        assertRaisesRegex(self, ValueError, error, list, s)

    def testUnknownSpecificationKey(self):
        """
        If an unknown key is given in a sequence specification, a ValueError
        must be raised.
        """
        error = "^Sequence specification 1 contains an unknown key: dog\.$"
        assertRaisesRegex(self, ValueError, error, Sequences, StringIO('''{
            "sequences": [
                {
                    "dog": "xxx"
                }
            ]
        }'''))

    def testUnknownSectionKey(self):
        """
        If an unknown key is given in a sequence specification section, a
        ValueError must be raised.
        """
        for key in 'id', 'id prefix', 'name':
            error = ("^Section 1 of sequence specification 1 contains an "
                     "unknown key: %s\.$" % key)
            assertRaisesRegex(self, ValueError, error, Sequences, StringIO('''{
                "sequences": [
                    {
                        "sections": [
                            {
                                "%s": "xxx"
                            }
                        ]
                    }
                ]
            }''' % key))

    def testOneSequenceIdOnlyDefaultLength(self):
        """
        If only one sequence is specified, and only by id, one sequence
        should be created, and it should have the length passed in
        defaultLength.
        """
        s = Sequences(StringIO('''{
            "sequences": [
                {
                    "id": "the-id"
                }
            ]
        }'''), defaultLength=500)
        (read,) = list(s)
        self.assertEqual('the-id', read.id)
        self.assertEqual(500, len(read.sequence))

    def testOneSequenceIdOnlyDefaultIdPrefix(self):
        """
        If only one sequence is specified, and only by length, one sequence
        should be created, it should have the length passed, and its id taken
        from passed defaultIdPrefix.
        """
        s = Sequences(StringIO('''{
            "globals": {
                "id prefix": "the-prefix."
            },
            "sequences": [
                {
                    "length": 5
                }
            ]
        }'''), defaultIdPrefix='the-prefix.')
        (read,) = list(s)
        self.assertEqual('the-prefix.1', read.id)
        self.assertEqual(5, len(read.sequence))

    def testOneSequenceAAOnly(self):
        """
        If only one sequence is specified, and only by indicating that it
        should be amino acids, one sequence should be created, it should have
        the default length, the expected id, and it should be entirely
        composed of nucleotides.
        """
        s = Sequences(StringIO('[{"name": "a", "random aa": true}]'))
        (read,) = list(s)
        self.assertEqual(Sequences.DEFAULT_ID_PREFIX + '1', read.id)
        self.assertEqual(Sequences.DEFAULT_LENGTH, len(read.sequence))
        self.assertEqual(set(), set(read.sequence) - set(AA_LETTERS))

    def testTwoSequencesSecondFromName(self):
        """
        If only one sequence is specified by name and a second refers to it
        by name, the second sequence should be the same as the first.
        """
        s = Sequences(StringIO('''[
            {
                "name": "a"
            },
            {
                "from name": "a"
            }
        ]'''))
        (read1, read2) = list(s)
        self.assertEqual(read1.sequence, read2.sequence)

    def testOneSequenceByLength(self):
        """
        If only one sequence is specified, and only by giving its length,
        only one sequence should be created, and it should have the given
        length.
        """
        s = Sequences(StringIO('[{"length": 55}]'))
        (read,) = list(s)
        self.assertEqual(55, len(read.sequence))

    def testOneSequenceRandomNTs(self):
        """
        A sequence must be able to be composed of random NTs.
        """
        s = Sequences(StringIO('''{
            "sequences": [
                {
                    "random nt": true
                }
            ]
        }'''))
        (read,) = list(s)
        self.assertEqual(set(), set(read.sequence) - set('ACGT'))

    def testOneSequenceIdPrefix(self):
        """
        A sequence must be able to be given just using an id prefix.
        """
        s = Sequences(StringIO('''{
            "sequences": [
                {
                    "id prefix": "xxx-"
                }
            ]
        }'''))
        (read,) = list(s)
        self.assertEqual('xxx-1', read.id)

    def testTwoSequences(self):
        """
        If two sequences are requested (only by giving a count) they should
        have the expected lengths and ids.
        """
        s = Sequences(StringIO('[{"count": 2}]'))
        (read1, read2) = list(s)
        self.assertEqual(Sequences.DEFAULT_LENGTH, len(read1.sequence))
        self.assertEqual(Sequences.DEFAULT_ID_PREFIX + '1', read1.id)
        self.assertEqual(Sequences.DEFAULT_LENGTH, len(read2.sequence))
        self.assertEqual(Sequences.DEFAULT_ID_PREFIX + '2', read2.id)

    def testOneSequenceLengthIsAVariable(self):
        """
        If only one sequence is specified, and only by giving its length,
        (as a variable) one sequence should be created, and it should have the
        given length.
        """
        s = Sequences(StringIO('''{
            "variables": {
                "len": 200
            },
            "sequences": [
                {
                    "length": "%(len)d"
                }
            ]
        }'''))
        (read,) = list(s)
        self.assertEqual(200, len(read.sequence))

    def testOneSectionWithLength(self):
        """
        A sequence must be able to be built up from sections, with just one
        section given by length.
        """
        s = Sequences(StringIO('''{
            "sequences": [
                {
                    "sections": [
                        {
                            "length": 40
                        }
                    ]
                }
            ]
        }'''))
        (read,) = list(s)
        self.assertEqual(40, len(read.sequence))

    def testOneSectionWithSequence(self):
        """
        A sequence must be able to be built up from sections, with just one
        section given by a sequence.
        """
        s = Sequences(StringIO('''{
            "sequences": [
                {
                    "sections": [
                        {
                            "sequence": "ACTT"
                        }
                    ]
                }
            ]
        }'''))
        (read,) = list(s)
        self.assertEqual('ACTT', read.sequence)

    def testOneSectionRandomNTs(self):
        """
        A sequence must be able to be built up from sections, with just one
        section of random NTs.
        """
        s = Sequences(StringIO('''{
            "sequences": [
                {
                    "sections": [
                        {
                            "random nt": true
                        }
                    ]
                }
            ]
        }'''))
        (read,) = list(s)
        self.assertEqual(set(), set(read.sequence) - set('ACGT'))

    def testOneSectionRandomAAs(self):
        """
        A sequence must be able to be built up from sections, with just one
        section of random AAs.
        """
        s = Sequences(StringIO('''{
            "sequences": [
                {
                    "sections": [
                        {
                            "random aa": true
                        }
                    ]
                }
            ]
        }'''))
        (read,) = list(s)
        self.assertEqual(set(), set(read.sequence) - set(AA_LETTERS))

    def testTwoSectionsWithLengths(self):
        """
        A sequence must be able to be built up from sections, with two
        sections given by length.
        """
        s = Sequences(StringIO('''{
            "sequences": [
                {
                    "sections": [
                        {
                            "length": 40
                        },
                        {
                            "length": 10
                        }
                    ]
                }
            ]
        }'''))
        (read,) = list(s)
        self.assertEqual(50, len(read.sequence))

    def testSectionWithNameReference(self):
        """
        A sequence must be able to be built up from sections, with two
        sections given by length.
        """
        s = Sequences(StringIO('''{
            "sequences": [
                {
                    "name": "xxx",
                    "sequence": "ACCGT"
                },
                {
                    "sections": [
                        {
                            "from name": "xxx"
                        },
                        {
                            "length": 10
                        }
                    ]
                }
            ]
        }'''))
        (read1, read2) = list(s)
        self.assertEqual(15, len(read2.sequence))
        self.assertTrue(read2.sequence.startswith('ACCGT'))

    def testSectionWithUnknownNameReference(self):
        """
        If a sequence is built up from sections and a referred to sequence
        does not exist, a ValueError must be raised.
        """
        s = Sequences(StringIO('''{
            "sequences": [
                {
                    "sections": [
                        {
                            "from name": "xxx"
                        }
                    ]
                }
            ]
        }'''))
        error = ("^Sequence section refers to name 'xxx' of "
                 "non-existent other sequence\.$")
        assertRaisesRegex(self, ValueError, error, list, s)

    def testSectionWithNameReferenceTooShort(self):
        """
        If a sequence is built up from sections and a referred to sequence
        is too short for the desired length, a ValueError must be raised.
        """
        s = Sequences(StringIO('''{
            "sequences": [
                {
                    "name": "xxx",
                    "sequence": "ACCGT"
                },
                {
                    "sections": [
                        {
                            "from name": "xxx",
                            "length": 10
                        }
                    ]
                }
            ]
        }'''))
        error = ("^Sequence specification refers to sequence name 'xxx', "
                 "starting at index 1 with length 10, but 'xxx' is not long "
                 "enough to support that\.$")
        assertRaisesRegex(self, ValueError, error, list, s)

    def testNamedRecombinant(self):
        """
        It must be possible to build up and name a recombinant.
        """
        s = Sequences(StringIO('''{
            "sequences": [
                {
                    "name": "xxx",
                    "sequence": "ACCA"
                },
                {
                    "name": "yyy",
                    "sequence": "GGTT"
                },
                {
                    "id": "recombinant",
                    "sections": [
                        {
                            "from name": "xxx",
                            "start": 1,
                            "length": 3
                        },
                        {
                            "from name": "yyy",
                            "start": 2,
                            "length": 2
                        }
                    ]
                }
            ]
        }'''))
        (read1, read2, read3) = list(s)
        self.assertEqual('recombinant', read3.id)
        self.assertEqual('ACCGT', read3.sequence)

    def testRecombinantFromFullOtherSequences(self):
        """
        It must be possible to build up a recombinant that is composed of two
        other sequences by only giving the names of the other sequences.
        """
        s = Sequences(StringIO('''{
            "sequences": [
                {
                    "name": "xxx",
                    "sequence": "ACCA"
                },
                {
                    "name": "yyy",
                    "sequence": "GGTT"
                },
                {
                    "id": "recombinant",
                    "sections": [
                        {
                            "from name": "xxx"
                        },
                        {
                            "from name": "yyy"
                        }
                    ]
                }
            ]
        }'''))
        (read1, read2, read3) = list(s)
        self.assertEqual('ACCAGGTT', read3.sequence)

    def testOneSequenceSequenceMutated(self):
        """
        A sequence should be be able to be mutated.
        """
        sequence = 'A' * 100
        s = Sequences(StringIO('''[{
            "sequence": "%s",
            "mutation rate": 1.0
        }]
        ''' % sequence))
        (read,) = list(s)
        # All bases should have been changed, due to a 1.0 mutation rate.
        diffs = sum((a != b) for (a, b) in zip(sequence, read.sequence))
        self.assertEqual(len(sequence), len(read.sequence))
        self.assertEqual(diffs, len(read.sequence))