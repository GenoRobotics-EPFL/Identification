#%%

#tasks:
#find a way to identify the family from the description
#attribute of the Bio.seq.seq object.

#ask how to organize the databases:
# - organize csv files by family -> code isnt provided in utils,
# one file of hundreds of csv files for each family
# - one csv file with all genes of all families
# and selection of family happens inside code
# 


# adapt build_database function to incorporate options 
# depending on gene and write comments for annick.

import pandas as pd
from Bio import SeqIO
import os.path as ospath



def build_database(option):
    """
        returns a dataframe with the genes (columns: Family, Species, Sequence)
        is dependent on the conservation of the "Identification file" organization.

    """

    genes = pd.DataFrame( columns= ["Family", "species", "sequence"])
    directory = ospath.abspath(os.getcwd())
    directory += '\Database'
    for file in os.listdir(directory):
        gene_path = ospath.join(directory, file)
        seq_objects = SeqIO.parse(gene_path, "fasta")
        i=0
        for record in seq_objects:
            genes.loc[len(genes)] = [extract_family(record.description,option), 
                                    extract_species(record.description, option),
                                    record.seq]
            i+=1
            if i==3:
                break
    return genes


def extract_species(description):
    seperated = description.split(" ")
    return seperated[1]+ " " + seperated[2]

def extract_family(description):
    seperated = description.split(" ")
    return seperated[1]


print(build_database())



# %%
