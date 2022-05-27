from ._constants import bases_inv, bases_inv_rna, codon_table, protein_codes



def invert_gen(source_gen, is_RNA=False):
    """
    Invert a generator.
    """
    # from ._sequence import Seq

    # # Get the desired lookup table
    # if is_RNA:    lookup = bases_inv_rna
    # else:         lookup = bases_inv

    # Return transformed generator
    # return ( seq.transform( lookup ) for seq in source_gen )
    return ( ~seq[::-1] for seq in source_gen )



### DNA <---> RNA ###


def rna_to_dna(rna_generator, as_string=False):
    if as_string:
        return ( chunk.replace("U", "T") for chunk in rna_generator )
    else:
        return ( chunk.replace(ord("U"), ord("T")) for chunk in rna_generator )


def dna_to_rna(seq, as_generator=False):

    if isinstance(seq, str):
        t = "T"
        u = "U"

    else:
        t = ord("T")
        u = ord("U")

    if as_generator:
        return ( chunk.replace(t, u) for chunk in seq )
    else:
        return seq.replace(t, u)



### RNA ----> Amino Acids ###

def _translate(seq, as_generator=False, from_DNA=False, frame=1):

    # If coming from DNA, insert intermediate transformer to RNA
    if from_DNA:
        seq = dna_to_rna(seq, as_generator=as_generator)

    if as_generator:

        # Take RNA sequence in groups of 3 and map to codons, ignoring if no map exists
        return (
            codon_table.get(str(codon), None) for codon in seq if len(codon) == 3
        )

    else:

        limit = ( (len(seq)) // 3) * 3

        return ''.join([
            codon_table.get(str(seq[i:i+3]), "?") for i in range(0, limit, 3)

            # codon_table.get(str( seq.subseq(i, i+3) ), None) for i in range(0, len(seq), 3)
            # str(seq[i:i+3]) for i in range(0, len(seq), 3)
            # str(seq[i:i+3]) + "-" for i in range(0, len(seq), 3)
        ])



def base_count(self, keys=protein_codes):

    keys = set([ ord(key) for key in keys ])

    base_count = { key: 0 for key in keys }

    for base in self:
        if base in keys:
            base_count[base] += 1

    base_count = { chr(key): val for key, val in base_count.items() }

    return base_count
