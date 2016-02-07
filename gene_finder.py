# -*- coding: utf-8 -*-
"""
This code analyzes a DNA sequence and outputs snippets of DNA that are likely to be protein-coding genes.

@author: Apurva Raman
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

    >>> get_complement('')
    Traceback (most recent call last):
    ...
    ValueError: No nucleotide was entered.

    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == '':
        raise ValueError("No nucleotide was entered.")
    else:
        raise ValueError("A non-nucleotide character was entered.")

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

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'

    checks that entering an empty dna string raises a ValueError
    >>> longest_ORF("")
    Traceback (most recent call last):
    ...
    ValueError: max() arg is an empty sequence

    """
    longest_ORF_string = max(find_all_ORFs_both_strands(dna))
    return longest_ORF_string

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    shuffle_longest_ORF = []
    for i in range(num_trials):
        shuffled_dna = shuffle_string(dna) #shuffle the input dna using the shuffle__string helper function.
        longest_length = len(longest_ORF(shuffled_dna)) 
        shuffle_longest_ORF.append(longest_length)
        longest_length_all_runs = max(shuffle_longest_ORF)
    return longest_length_all_runs

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
        no doctests added: 
        just handling these example strand should be sufficient to demonstrate the function.
        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    codon = '' 
    amino_acid_sequence = ''
    for x in range(len(dna) / 3):
        codon = dna[(3 * x):(3 * x + 3)]
        if len(codon) == 3:
            amino_acid = aa_table[codon]
        amino_acid_sequence += amino_acid
        x += 3
    return amino_acid_sequence


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    randomness_threshold = longest_ORF_noncoding(dna,1500)
    all_ORFs = find_all_ORFs_both_strands(dna)
    non_random_ORFs = [i for i in all_ORFs if len(i) > randomness_threshold]
    protein = [coding_strand_to_AA(j) for j in non_random_ORFs]
    return protein


from load import load_seq
dna = load_seq("./data/X73525.fa")
print gene_finder(dna)


#if __name__ == "__main__":
    #import doctest
    #doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose = True)
