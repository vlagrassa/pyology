from array import array

from ._transform import _translate
from ._constants import base_to_bin, bin_to_dna, bin_to_rna

class Seq(array):

    # Class variables: Information about the DNA bases
    from ._constants import bases


    ### Initialization ###

    def __new__(cls, seq=None):
        return array.__new__(cls, 'B')

    def __init__(self, seq=None, rna=False):

        from ._reader import SeqReader

        self._is_rna = rna

        self.translate = self._instance_translate

        if seq is not None:
            if isinstance(seq, self.__class__):
                self.extend(seq)

            elif isinstance(seq, array):
                self.extend(seq)

            elif isinstance(seq, str):
                super().extend(bytearray(seq, 'utf8'))

            elif isinstance(seq, SeqReader):
                for subseq in seq.as_DNA():
                    self.extend(subseq)

            else:
                try:
                    self.fromlist(seq)
                except:
                    try:
                        self.append(seq)
                    except:
                        raise TypeError(f"Couldn't construct Seq from {seq}.")


    ### Sequence Alignment ###
    from .align import align as align
    align = staticmethod(align)


    ### Translating between DNA and RNA ###
    from ._transform import dna_to_rna, rna_to_dna
    dna_to_rna = staticmethod(dna_to_rna)
    rna_to_dna = staticmethod(rna_to_dna)

    from ._transform import base_count


    def scan(self, window_size):
        for i in range(0, len(self), window_size):
            yield self[i : i + window_size]



    ### String Functions ###

    def __repr__(self):
        if len(self) > 10:
            return f"<{self[:10].tobytes().decode('utf8')}...>"
        else:
            return f"<{self}>"

    def __str__(self):
        return super().tobytes().decode("utf8")


    ### Concatenation ###
    
    # Concatenation (+ operator)
    def __add__(self, other):

        # Concat two Seqs
        if isinstance(other, self.__class__):
            new_seq = Seq()
            new_seq.extend(self)
            new_seq.extend(other)
            return new_seq

        # Allow concatenation with a string
        elif isinstance(other, str):
            new_seq = Seq()
            new_seq.extend(self)
            new_seq.extend(bytearray(other, 'utf8'))
            return new_seq

        else:
            return None

    # Concatenation (+ operator), where Sequence is the second element
    def __radd__(self, other):

        # Concat two Seqs
        if isinstance(other, self.__class__):
            new_seq = Seq()
            new_seq.extend(other)
            new_seq.extend(self)
            return new_seq

        # Allow concatenation with a string
        elif isinstance(other, str):
            new_seq = Seq(other)
            new_seq.extend(self)
            return new_seq

        else:
            return None

    # Concatenation, in-place (+= operator)
    def __iadd__(self, other):

        # Concat two Seqs
        if isinstance(other, self.__class__):
            self.extend(other)
            return self

        # Allow concatenation with a string
        elif isinstance(other, str):
            self.extend(bytearray(other, 'utf8'))
            return self

        else:
            return None


    # Static class method for concatenating two of more sequences
    @staticmethod
    def concat(*seqs):

        if len(seqs) == 1:
            seqs = seqs[0]

        full_seq = Seq()
        for seq in seqs:
            full_seq += seq
        return full_seq


    ### Repetition ###

    # Repetition (* operator)
    def __mul__(self, other):
        seq = Seq()
        for i in range(other):
            seq.extend(self)
        return seq
        # return Seq(str(self) * other)


    ### Logical Operations - Matching Sequences Pairwise ###


    # Compute the combination of both sequences
    # If two bases match, keep them, otherwise insert a gap
    def __and__(self, other):

        from ._reader import SeqReader

        if isinstance(other, SeqReader):
            return other.__rand__(self)

        if not isinstance(other, Seq):
            other = Seq(other)

        return Seq([
            # x if x == y else ord('-')
            bin_to_dna[base_to_bin[x] & base_to_bin[y]]
            for x, y in zip(self, other)
        ])

    def __rand__(self, other):
        return Seq(other) & self

    def __iand__(self, other):
        self = self & other
        return self


    # Compute the sequence that represents both arguments
    def __or__(self, other):

        from ._reader import SeqReader

        if isinstance(other, SeqReader):
            return other.__ror__(self)

        return Seq([
            bin_to_dna[base_to_bin[x] | base_to_bin[y]]
            for x, y in zip(self, other)
        ])

    def __ror__(self, other):
        return Seq(other) | self

    def __ior__(self, other):
        self = self | other
        return self


    # Keep matching bases, overwriting gaps; replace mismatches with gaps
    # Like and, except for gaps
    def __xor__(self, other):

        from ._reader import SeqReader

        if isinstance(other, SeqReader):
            return other.__rxor__(self)

        gap = ord('-')

        return Seq([
            x if x == y or y == gap else (y if x == gap else gap)
            for x, y in zip(self, other)
        ])

    def __rxor__(self, other):
        return Seq(other) ^ self

    def __ixor__(self, other):
        self = self ^ other
        return self


    # For each position, remove possibility of the base in `other` existing there
    def __sub__(self, other):
        return Seq([
            bin_to_dna[base_to_bin[x] & ~base_to_bin[y]]
            for x, y in zip(self, other)
        ])

    def __rsub__(self, other):
        return Seq(other) - self

    def __isub__(self, other):
        self = self - other
        return self



    ### Complements ###

    # Compute the complement of a sequence
    def __invert__(self):
        return Seq([
            bin_to_dna[((base_to_bin[b] >> 2) & 0b0011) | ((base_to_bin[b] << 2) & 0b1100)]
            for b in self[::-1]
        ])

    # Check whether two sequences are complements of one another
    def __mod__(self, other):
        if isinstance(other, self.__class__):
            return self == ~other
        elif isinstance(other, str):
            return ~self == other
        else:
            raise TypeError(f"unsupported operand type(s) for %: '{self.__class__.__name__}' and '{other.__class__.__name__}'")



    ### Bitshifting -- Add gaps to either end ###

    def __lshift__(self, n):
        return self + ('-' * n)

    def __ilshift__(self, n):
        self += ('-' * n)
        return self

    def __rshift__(self, n):
        return ('-' * n) + self

    def __irshift__(self, n):
        self.reverse()
        self += ('-' * n)
        self.reverse()
        return self



    # Perform some functions by just converting to string

    def __hash__(self):
        return str(self).__hash__()

    def __contains__(self, other):
        return str(self).__contains__(str(other))



    ### Replace ###

    def replace(self, source, target):
        return Seq([ target if b == source else b for b in self ])

    def transform(self, lookup):
        return Seq([ lookup.get(b, b) for b in self ])



    ### Translation ###

    def to_DNA(self):
        return Seq.rna_to_dna(self)

    def to_RNA(self):
        return Seq.dna_to_rna(self)

    def as_DNA(self, complement=False, offset=0, window_size=None):

        # Get the default window size, if none is provided
        if window_size is None:
            window_size = len(self)

        yield from self.scan(window_size=window_size)

    def as_RNA(self, complement=False, offset=0, window_size=None):

        # Get the default window size, if none is provided
        if window_size is None:
            window_size = len(self)

        yield from self.to_RNA().scan(window_size=window_size)

    @staticmethod
    def translate(seq, **kwargs):
        if kwargs.get("as_generator", False):
            return _translate(seq, **kwargs)
        else:
            return Seq(_translate(seq.to_RNA(), **kwargs))

    def _instance_translate(self, **kwargs):
        kwargs["as_generator"] = False
        return Seq( _translate(self.to_RNA(), **kwargs) )


    ### Indexing ###

    def __getitem__(self, elem):
        return Seq(super().__getitem__(elem))

    def get_raw(self, elem):
        return super().__getitem__(elem)

    # def __iter__(self):
    #     # print( [ chr(b) for b in iter(super(array, self)) ] )

    #     print(type(self))
    #     print(super(array, self))

    #     for b in super(array, self):
    #         print(b)
    #     # print( [ self ] )
    #     # print("a")
    #     return iter("abc")
    #     # return super().__iter__()

