from structures import *
from DNAToolkit import *

print(f"All prots in 6 open reading frames:")
for prot in all_proteins_from_orfs(NM_001313993_2, 0, 0, True):
    print(f"{prot}")