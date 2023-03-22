#tasks:
#find a way to identify the family from the description
#attribute of the Bio.seq.seq object.

#ask how to organize the databases:
# - organize csv files by family -> code isnt provided in utils,
# one file of hundreds of csv files for each family
# - one csv file with all genes of all families
# and selection of family happens inside code


import pandas as pd
from Bio import SeqIO
import os.path as ospath
import os



def build_database():
    """
    Organizes the GenoRobotics species database into families 
    Note
    ----
    The databases should be in a separate folder (in the working directory) named "Database"
    Parameters
    ----------
        db_name : string, name of the data base to organize 
                options -> "matk", "psba_trnh", "rbcl", "its"
    Returns
    -------
    Dataframe with the sorted genes (columns: Family, Species, Sequence)
    """

    genes = pd.DataFrame(columns= ["Gene", "Family", "species", "sequence"])
    directory = ospath.abspath(os.getcwd() + "/Database") # path to the database 

    for file in os.listdir(directory):
        print(file)
        gene_path = ospath.join(directory, file)
        seq_objects = SeqIO.parse(gene_path, "fasta")

        i=0
        for record in seq_objects:
            # .loc(genes) = append (create an entry at end)
            genes.loc[len(genes)] = [extract_gene(file), 
                                    extract_order(record.description, file), 
                                    extract_species(record.description, file),
                                    record.seq]
            i+=1 
            if i==3:
                break
    return genes 
    # right now, all dataframes end up in one variable -> genes (I added gene type to figure out which is what)

def extract_gene(db_name):
    """
    Extracts the species name from the corresponding fasta file entry
    Parameters
    ----------
        Description = fasta formatted description of the gene sequence
        db_name = name of gene, whose database we are extracting its species 
    Returns
    -------
    Species name, string
    """
    if db_name == "sequences_matK_800-1550.fasta":
        return "matk"
    elif db_name == "psbA-trnH_sequence.fasta":
        return "psba-trnh"
    elif db_name == "rcbL_sequence.fasta":
        return "rbcl"
    elif db_name == "its":
        return "its"
    else:
        raise ValueError("No valid file name")

def extract_species(description, db_name):
    """
    Extracts the species name from the corresponding fasta file entry
    Parameters
    ----------
        Description = fasta formatted description of the gene sequence
        db_name = name of gene, whose database we are extracting its species 
    Returns
    -------
    Species name, string
    """
    separated = description.split(" ")
    if db_name == "sequences_matK_800-1550.fasta":
        return separated[1] + ' ' + separated[2]
    elif db_name == "psbA-trnH_sequence.fasta":
        return separated[1] + ' ' + separated[2]
    elif db_name == "rcbL_sequence.fasta":
        return separated[1] + ' ' + separated[2]
    elif db_name == "its":
        return
    else:
        raise ValueError("No valid file name")

def extract_order(description, db_name):
    """
    Extracts the order name from the corresponding fasta file entry
    Parameters
    ----------
        Description = fasta formatted description of the gene sequence
        db_name = name of gene, whose database we are extracting its order
    Returns
    -------
    Order name, string
    """
    separated = description.split(" ")
    if db_name == "sequences_matK_800-1550.fasta":
        return separated[1] 
    elif db_name == "psbA-trnH_sequence.fasta":
        return separated[1] 
    elif db_name == "rcbL_sequence.fasta":
        return separated[1] 
    elif db_name == "its":
        return
    else:
        raise ValueError("No valid file name")

print(build_database())