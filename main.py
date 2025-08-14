from structures import *
from DNAToolkit import *
import random

randDNASeq = ''.join([random.choice(Nucleotides)
                      for nuc in range(50)])

DNAStr = validateSeq(randDNASeq)

print(f"\nSequence: {DNAStr}\n")
print(f"Sequence Lenght: {len(DNAStr)}\n")