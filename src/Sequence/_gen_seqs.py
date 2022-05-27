def gen_palindromes(n):
    """
    Generator for sequences that are their own complement.
    """
    assert n % 2 == 0, \
        "Palindromes must have even length. (What would we do with the center base?)"

    # Base case: Generate an empty sequence
    if n == 0: yield Seq()

    # Step case: Append all combinations of base & inverted base around inner recursion
    else:
        yield from (
            base + inner + ~base
            for base  in Bases
            for inner in Seq.gen_palindromes(n-2)
        )


def gen_unique_seqs(n):
    """
    Generator for sequences that are not their own complement.
    """
    # Consider all permutations of the bases, yielding only those less than their own inverse
    #   - Since (<) is not commutative, only one of the two can ever be generated
    #   - Since a palindromic sequence is equal to its own inverse, it can't be
    #     less than it, meaning no palindromic sequences can be generated
    return (
        seq for seq in permutations(Bases, n, Seq()) if seq < ~seq
    )