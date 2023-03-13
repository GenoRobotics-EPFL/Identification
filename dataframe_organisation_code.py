#%%

import pandas as pd
genes = pd.DataFrame( columns= ["Family", "species", "sequence"])
genes.loc[len(genes)] = ['Homo', 'Homo Sapiens', 'ATGCGTAGCTAGCT']

# %%
def extract_species(description):
    seperated = description.split(" ")
    return seperated[1]+ " " + seperated[2]

def extract_family(description):
    seperated = description.split(" ")
    return seperated[1]

# %%
from Bio import SeqIO
gene_path = "psbA-trnH_sequence.fasta"
seq_objects = SeqIO.parse(gene_path, "fasta")
i=0
for record in seq_objects:
    genes.loc[len(genes)] = [extract_family(record.description), 
                             extract_species(record.description),
                             record.seq]
    i += 1
    if i == 3:
        break
print(genes)


# %%
