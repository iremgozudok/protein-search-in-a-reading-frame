from structures import *
from collections import Counter


# DNA validation
def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq



def translate_seq(seq, init_pos=0):
    return [DNA_Codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)]


def reverse_complement(seq):
    """Swapping adenine with thymine and guanine with cytosine. Reversing newly generated string"""
    return ''.join([DNA_ReverseComplement[nuc] for nuc in seq][::-1])


def gen_reading_frames(seq):
    """Generate the six reading frames of a DNA sequence, including the reverse complement"""
    frames = []
    frames.append(translate_seq(seq, 0))
    frames.append(translate_seq(seq, 1))
    frames.append(translate_seq(seq, 2))
    frames.append(translate_seq(reverse_complement(seq), 0))
    frames.append(translate_seq(reverse_complement(seq), 1))
    frames.append(translate_seq(reverse_complement(seq), 2))
    return frames


def proteins_from_rf(aa_seq):
    """Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""
    current_prot =[]
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            # Stop accumulating amino acids if _ - Stop was found
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            # Start accumulating amino acids if M - Start was found
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins