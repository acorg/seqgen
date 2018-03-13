from six import string_types, PY3
from six.moves import builtins
from json import load
from random import choice, uniform
from dark.aa import AA_LETTERS
from dark.fasta import FastaReads
from dark.reads import Read


class Sequences(object):
    """
    Create genetic sequences from a JSON specification.

    @param spec: A C{str} filename or an open file pointer to read the
        specification from.
    @raise json.decoder.JSONDecodeError: If the specification JSON cannot
        be read.
    @raise ValueError: If the specification JSON is an object but does not
        have a 'sequences' key.
    """
    NT = list('ACGT')
    AA = list(AA_LETTERS)
    DEFAULT_LENGTH = 100
    DEFAULT_ID_PREFIX = 'seq-id-'
    LEGAL_SPEC_KEYS = {
        'count',
        'description',
        'id',
        'id prefix',
        'from name',
        'length',
        'mutation rate',
        'name',
        'random aa',
        'random nt',
        'sections',
        'sequence',
        'sequence file',
    }
    LEGAL_SPEC_SECTION_KEYS = {
        'from name',
        'length',
        'mutation rate',
        'random aa',
        'random nt',
        'start',
        'sequence',
        'sequence file',
    }

    def __init__(self, spec, defaultLength=None, defaultIdPrefix=None):
        self._defaultLength = defaultLength or self.DEFAULT_LENGTH
        self._defaultIdPrefix = defaultIdPrefix or self.DEFAULT_ID_PREFIX
        self._readSpecification(spec)

    def _readSpecification(self, spec):
        """
        Read the specification in C{spec}.

        @param spec: A C{str} filename or an open file pointer to read the
            specification from.
        @raise KeyError: if the specification JSON is an object and does not
            have a 'sequences' key.
        """
        if isinstance(spec, string_types):
            with open(spec) as fp:
                j = load(fp)
        else:
            j = load(spec)

        if isinstance(j, list):
            vars_, sequenceSpecs = {}, j
        else:
            try:
                vars_, sequenceSpecs = j.get('variables', {}), j['sequences']
            except KeyError:
                raise ValueError("The specification JSON must have a "
                                 "'sequences' key.")

        self._vars = vars_
        self._sequenceSpecs = list(map(self._expandSpec, sequenceSpecs))
        self._checkKeys()
        self._names = {}

    def _checkKeys(self):
        """
        Check that all specification dicts only contain legal keys.

        @param sequenceSpec: A C{dict} with information about the sequences
            to be produced.
        @raise ValueError: If an unknown key is found.
        """
        for specCount, spec in enumerate(self._sequenceSpecs, start=1):
            unexpected = set(spec) - self.LEGAL_SPEC_KEYS
            if unexpected:
                raise ValueError(
                    'Sequence specification %d contains %sunknown key%s: %s.' %
                    (specCount, 'an ' if len(unexpected) == 1 else '',
                     '' if len(unexpected) == 1 else 's',
                     ', '.join(sorted(unexpected))))
            try:
                sections = spec['sections']
            except KeyError:
                pass
            else:
                for sectionCount, section in enumerate(sections, start=1):
                    unexpected = set(section) - self.LEGAL_SPEC_SECTION_KEYS
                    if unexpected:
                        raise ValueError(
                            'Section %d of sequence specification %d contains '
                            '%sunknown key%s: %s.' %
                            (sectionCount, specCount,
                             'an ' if len(unexpected) == 1 else '',
                             '' if len(unexpected) == 1 else 's',
                             ', '.join(sorted(unexpected))))

    def _expandSpec(self, sequenceSpec):
        """
        Recursively expand all string values in a sequence specification.

        @param sequenceSpec: A C{dict} with information about the sequences
            to be produced.
        @return: A C{dict} with all string values expanded.
        """
        new = {}
        for k, v in sequenceSpec.items():
            if isinstance(v, string_types):
                value = v % self._vars
                # If a substitution was done and the converted string is
                # all digits, convert to int. Or if it can be converted to
                # a float do that.
                if value != v:
                    if all(str.isdigit(x if PY3 else str(x)) for x in value):
                        value = int(value)
                    else:
                        try:
                            value = float(value)
                        except ValueError:
                            pass

            elif isinstance(v, dict):
                value = self._expandSpec(v)
            else:
                value = v
            new[k] = value
        return new

    def _specToRead(self, spec):
        """
        Get a sequence from a specification.

        @param spec: A C{dict} with keys/values specifying a sequence.
        @raise ValueError: If the section spec refers to a non-existent other
            sequence, or to part of another sequence but the requested part
            exceeds the bounds of the other sequence. Or if the C{spec} does
            not have a 'length' key when no other sequence is being referred
            to.
        @return: A C{dark.Read} instance.
        """
        nt = True
        length = spec.get('length', self._defaultLength)

        if 'from name' in spec:
            name = spec['from name']
            try:
                namedRead = self._names[name]
            except KeyError:
                raise ValueError("Sequence section refers to name '%s' of "
                                 "non-existent other sequence." % name)
            else:
                # The start offset in the spec is 1-based. Convert to 0-based.
                index = int(spec.get('start', 1)) - 1
                # Use the given length (if any) else the length of the
                # named read.
                length = spec.get('length', len(namedRead))
                sequence = namedRead.sequence[index:index + length]

                if len(sequence) != length:
                    raise ValueError(
                        "Sequence specification refers to sequence name '%s', "
                        "starting at index %d with length %d, but '%s' "
                        "is not long enough to support that." %
                        (name, index + 1, length, name))

                read = Read(None, sequence)

        elif 'sequence' in spec:
            read = Read(None, spec['sequence'])

        elif 'sequence file' in spec:
            noFileClass = builtins.FileNotFoundError if PY3 else IOError
            reads = iter(FastaReads(spec['sequence file']))
            try:
                read = next(reads)
            except StopIteration:
                raise ValueError("Sequence file '%s' is empty." %
                                 spec['sequence file'])
            except noFileClass:
                raise ValueError("Sequence file '%s' could not be read." %
                                 spec['sequence file'])

        elif spec.get('random aa'):
            read = Read(None, ''.join(choice(self.AA) for _ in range(length)))
            nt = False

        else:
            # Assume random nucleotides are wanted.
            read = Read(None, ''.join(choice(self.NT) for _ in range(length)))

        try:
            rate = spec['mutation rate']
        except KeyError:
            pass
        else:
            read.sequence = self._mutate(read.sequence, rate, nt)

        return read

    def _mutate(self, sequence, rate, nt):
        """
        Mutate a sequence at a certain rate.

        @param sequence: A C{str} nucleotide or amino acid sequence.
        @param rate: A C{float} mutation rate.
        @param nt: If C{True} the sequence is nucleotides, else amino acids.
        @return: The mutatated C{str} sequence.
        """
        result = []
        if nt:
            possibles = set(self.NT)
        else:
            possibles = set(self.AA)
        for current in sequence:
            if uniform(0.0, 1.0) < rate:
                result.append(choice(list(possibles - {current})))
            else:
                result.append(current)

        return ''.join(result)

    def _readsForSpec(self, spec, start):
        """
        Yield reads for a given specification.

        @param sequenceSpec: A C{dict} with information about the sequences
            to be produced.
        @param start: The C{int} starting count for sequence ids (that do not
            supply their own id).
        """
        nSequences = spec.get('count', 1)
        for count in range(nSequences):
            id_ = None
            if 'sections' in spec:
                sequence = ''
                for section in spec['sections']:
                    sequence += self._specToRead(section).sequence
            else:
                read = self._specToRead(spec)
                sequence = read.sequence
                id_ = read.id

            if id_ is None:
                try:
                    id_ = spec['id']
                except KeyError:
                    id_ = '%s%d' % (
                        spec.get('id prefix', self._defaultIdPrefix),
                        start + count)
                else:
                    # If an id is given, the number of sequences requested must
                    # be one.
                    if nSequences != 1:
                        raise ValueError(
                            "Sequence with id '%s' has a count of %d. If you "
                            "want to specify one sequence with an id, the "
                            "count must be 1. To specify multiple sequences "
                            "with an id prefix, use 'id prefix'." %
                            (id_, nSequences))

            try:
                id_ = id_ + ' ' + spec['description']
            except KeyError:
                pass

            read = Read(id_, sequence)

            # Keep a reference to this result if it is named.
            if count == 0 and 'name' in spec:
                name = spec['name']
                if name in self._names:
                    raise ValueError("Name '%s' is duplicated in the JSON "
                                     "specification." % name)
                else:
                    self._names[name] = read

            yield read

    def __iter__(self):
        startCount = 1
        for sequenceSpec in self._sequenceSpecs:
            thisCount = 0
            for read in self._readsForSpec(sequenceSpec, startCount):
                yield read
                thisCount += 1
            startCount += thisCount
