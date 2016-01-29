# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE
    What is the right content for a header comment? 
@author: YOUR NAME HERE 
    Me? I'm Apurva Raman
    apurvaraman.github.com

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'

    checks if G maps to C
    >>> get_complement('G')
    'C'

    checks if T maps to A
    >>> get_complement('T')
    'A'

    checks if a non-nucleotide raises an error
    >>> get_complement('Q')
    Traceback (most recent call last):
    ...
    ValueError: A non-nucleotide character was entered.

    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    else:
        raise ValueError("A non-nucleotide character was entered.")
    #FINISHED: Implemented and tested
    pass


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    
    checks that strings with non-nucleotide characters in them raise an error
    >>> get_reverse_complement("ATGCCCGQWERTY")
    Traceback (most recent call last):
    ...
    ValueError: A non-nucleotide character was entered.
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """

    build_complement = '' #initializes an empty string that the complements will be added to
    for i in dna:
        build_complement += get_complement(i)
    reverse_complement_dna = build_complement[::-1]
    return reverse_complement_dna

    # FINISHED: IMPLEMENTED AND TESTED
    pass


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    checks that the function returns the whole string if there is no STOP codon
    >>> rest_of_ORF('ATGAGA')
    'ATGAGA'
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    #The STOP codons are TAA, TAG, and TGA
    for i in range(len(dna) / 3): #iterating through every frame
        if dna[3*i:3*i+3] == 'TAG' or dna[3*i : 3*i+3] == 'TAA' or dna[3*i : 3*i+3] == 'TGA': #checking if one of the frames is the STOP codon
            return dna[:3*i]
    return dna
    # FINISHED: IMPLEMENTED AND TESTED
    pass


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    checks that if there is no START codon in dna, the result will be an empty list
    >>> find_all_ORFs_oneframe("CATCATCATCATCAT")
    []
    checks that the resulting list starts with the START codon
    >>> find_all_ORFs_oneframe("TTTCTTATTATGGTTTCTCAT")
    ['ATGGTTTCTCAT']
    checks that the resulting list ends with the STOP codon
    >>> find_all_ORFs_oneframe("TTTCTTATTATGGTTTCTTAACAT")
    ['ATGGTTTCT']
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    ORFs = []
    index = 0
    while index < len(dna):
        if dna[index: index + 3] == 'ATG':
            new_ORF = len(rest_of_ORF(dna[index:]))
            ORFs.append(rest_of_ORF(dna[index:]))
            index += new_ORF
        else:
            index += 3
    return ORFs
    # FINISHED: IMPLEMENTED AND TESTED
    pass


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    checks that if there is no START codon in dna, the result will be an empty list
    >>> find_all_ORFs('CATCATCATCATCAT')
    []
    checks that the resulting list ends with the STOP codon
    >>> find_all_ORFs('ATGATGTAA')
    ['ATGATG']
    >>> find_all_ORFs('CTTTCTTATTATGGTTTCTCAT')
    checks that an "out of frame" START codon does not result in an empty list
    ['ATGGTTTCTCAT']
    >>> find_all_ORFs('ATGCATGAATGTAG')
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    ORFs_allframes = []
    ORFs_allframes.extend(find_all_ORFs_oneframe(dna)) #first frame
    ORFs_allframes.extend(find_all_ORFs_oneframe(dna[1:])) #second frame
    ORFs_allframes.extend(find_all_ORFs_oneframe(dna[2:])) #third frame
    return ORFs_allframes
    # FINISHED: IMPLEMENTED AND TESTED
    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    checks that the result ignores nested ORFs and shows strands from both ORFs from the START to the STOP codon
    >>> find_all_ORFs_both_strands('TTAATGATGATGCATCATTAA')
    ['ATGATGATGCATCAT', 'ATGATGCATCATCAT']
    >>> find_all_ORFs_both_strands('ATGCGAATGTAGCATCAAA')
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    ORFs_both_strands = []
    ORFs_both_strands.extend(find_all_ORFs(dna)) #first strand
    ORFs_both_strands.extend(find_all_ORFs(get_reverse_complement(dna))) #second strand
    return ORFs_both_strands

    # TODO: implement this
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(find_all_ORFs_both_strands, globals(), verbose = True)
