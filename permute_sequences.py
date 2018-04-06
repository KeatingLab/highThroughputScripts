import argparse
import random
import os

amino_acids = "ACEDGFIHKMLNQPSRTWVY"
bases = "AGCT"

def detect_language(template):
    for char in template:
        if char not in amino_acids and char in bases:
            print("Inferring DNA.")
            return bases
        elif char not in bases and char in amino_acids:
            print("Inferring protein.")
            return amino_acids
        elif len(char) > 1:
            try:
                # This is a set of allowed characters, so check each character individually
                return detect_language(char)
            except AssertionError:
                pass
    assert False, "Can't detect whether to use amino acids or nucleotides"

def all_permutations_recursive(template, base, language):
    if len(template) == 0:
        yield base
    elif template[0] == '*':
        for i, char in enumerate(language):
            if len(base) == 0:
                print("{} / {} first positions".format(i, len(language)))
            yield from all_permutations_recursive(template[1:], base + char, language)
    elif len(template[0]) > 1:
        # This item in the template contains a list of possible values
        for i, char in enumerate(template[0]):
            yield from all_permutations_recursive(template[1:], base + char, language)
    else:
        yield from all_permutations_recursive(template[1:], base + template[0], language)

def all_permutations(template, alphabet=None):
    """Computes and yields all permutations of the given DNA/protein template in a memory-efficient manner.

    This function automatically detects whether to use the DNA or the protein alphabet depending on the content of the template, unless you specify one or the other using the 'alphabet' parameter. The template supports specifying one possible residue at a position, a set of possible residues, or any residue.

    Parameters:
    template -- An iterable that specifies the possible values at each position. Each index in the iterable should be one of the following: (1) a character in the DNA/protein alphabet, (2) an iterable of characters, or (3) the '*' character denoting a wildcard.
    alphabet -- 'aa', 'dna', a different string of characters, or None. If it is not None, the alphabet for permutation will be set to the according language, or to the custom alphabet you specify. If None, the alphabet will be determined automatically.

    Returns:
    A generator that permutes all possible residues according to the given template.

    Examples:
    >>> list(all_permutations("A*", alphabet="dna"))
    ["AA", "AG", "AC", "AT"]

    >>> list(all_permutations(["A", "CG", "*"], alphabet="dna"))
    ["ACA", "ACG", "ACC", "ACT", "AGA", "AGG", "AGC", "AGT"]
    """
    if alphabet == 'aa':
        language = amino_acids
    elif alphabet == 'dna':
        language = bases
    elif alphabet is not None:
        language = alphabet
    else:
        language = detect_language(template)

    yield from all_permutations_recursive(template, "", language)

if __name__ == '__main__':
    pass
