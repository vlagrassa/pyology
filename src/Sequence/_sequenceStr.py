
class SeqStr:

    # Class variables: Information about the DNA bases
    bases     = "ACGT"
    bases_inv = {  "A": "T",  "C": "G",  "G": "C",  "T": "A"  }
    
    def __init__(self, seq=None):
        if seq is None: self._seq = ""
        else:           self._seq = seq
    
    def __invert__(self):
        return SeqStr("".join([ self.__class__.bases_inv[b] for b in self._seq[::-1] ]))

    def __len__(self):
        return len(self._seq)


    ### String Functions ###

    def __repr__(self):
        if len(self._seq) > 10:
            return f"<{self._seq[:10]}...>"
        else:
            return f"<{self._seq}>"

    def __str__(self):
        return self._seq


    ### Comparisons ###

    # Equality (== operator)
    def __eq__(self, other):
        if isinstance(other, self.__class__): return self._seq == other._seq
        else:                                 return False

    # Inequality (!= operator)
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __lt__(self, other):
        return self._seq < other._seq


    ### Concatenation ###
    
    # Concatenation (+ operator)
    def __add__(self, other):
        if   isinstance(other, self.__class__): return SeqStr(self._seq + other._seq)
        elif isinstance(other, str):            return SeqStr(self._seq + other)
        else:                                   return None

    # Concatenation (+ operator), where Sequence is the second element
    def __radd__(self, other):
        if   isinstance(other, self.__class__): return SeqStr(other._seq + self._seq)
        elif isinstance(other, str):            return SeqStr(other + self._seq)
        else:                                   return None

    # Static class method for concatenating two of more sequences
    @staticmethod
    def concat(*seqs):

        if len(seqs) == 1:
            seqs = seqs[0]

        full_seq = SeqStr()
        for seq in seqs:
            full_seq += seq
        return full_seq


    ### Repetition ###

    # Repetition (* operator)
    def __mul__(self, other):
        return SeqStr(self._seq * other)


    ### Matching with Other Sequences ###

    # Matching operator
    def __and__(self, other):
        return SeqStr(''.join([
            x if x == y else '\u00b7'
            for x, y in zip(self._seq, other._seq)
        ]))

    # Modulus - check whether two Seqs are complements of one another
    def __mod__(self, other):
        if isinstance(other, self.__class__):
            return self == ~other
        else:
            raise TypeError(f"unsupported operand type(s) for %: '{self.__class__.__name__}' and '{other.__class__.__name__}'")

    def __hash__(self):
        return self._seq.__hash__()
    
    def __getitem__(self, elem):
        return SeqStr(self._seq.__getitem__(elem))
