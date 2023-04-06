#tasks:
#find a way to identify the family from the description
#attribute of the Bio.seq.seq object.

#ask how to organize the databases:
# - organize csv files by family -> code isnt provided in utils,
# one file of hundreds of csv files for each family
# - one csv file with all genes of all families
# and selection of family happens inside code

Entrez.email = 'eva.frossard@epfl.ch'

import pandas as pd
from Bio import SeqIO
import os.path as ospath
import os
from Bio import Entrez

hierarchy = pd.read_csv(os.getcwd() + "/hierarchy_report.csv")

def build_database():
    """
    Creates a new folder that organizes the Genorobotics database by family.
    The folder is called Database_by_family.
    This folder contains:
    -a csv with the gene entries that were not treated properly and need manual intervention
    - sub-folders containing the name of the family, each containing
    four csv files for each of the barcoding genes: matK, psba_trnh, rbcl, its.

    Note
    ----
    The initial databases should be in a separate folder (in the working directory) named "Database"

    Returns
    -------
    None
    """

    initial_dir = ospath.abspath(os.getcwd() + "/Database") # path to the initial database
    os.mkdir("Database_by_family") #creates a folder for saving the new_database
    creating_dir = ospath.abspath(os.getcwd() + "/Database_by_family") # and saves its dir

    #creates dataframe to deal with undefined gene entries
    manual_db = pd.DataFrame(columns = ["gene", "description", "sequence"])
    manual_db.to_csv(creating_dir+"//manual_treatment")

    for file in os.listdir(initial_dir):

        gene_path = ospath.join(initial_dir, file) #path to any of the four gene files
        seq_objects = SeqIO.parse(gene_path, "fasta") #creates iterable SeqRecord object
        gene = extract_gene(file)

        for record in seq_objects:

            #elements for the dataframe,
            gene_id = extract_id(record.description)
            (extraction_code, genus, species) = extract_binomial(record.description)
            sequence = record.seq

            if extraction_code == 0:
                new_row = pd.DataFrame({"gene":[gene], "description":[record.description], "sequence": [sequence]})
                new_row.to_csv(ospath.join(creating_dir + "//manual_treatment"), mode = 'a', header = False)
            else:
                new_row = pd.DataFrame({"gene_id":[gene_id], "genus": [genus], "species":[species], "sequence": [sequence]})
                family = get_families(get_taxids([genus + species]))
            
                if family not in os.listdir(creating_dir):
                    #the family's folder has not yet been created
                    os.mkdir(ospath.join(creating_dir,family))
                    family_dir = ospath.join(creating_dir,family)
                    new_row.to_csv(family_dir + "\\" +gene)

            
                else:
                    #the family's folder has already been created
                    family_dir = ospath.join(creating_dir,family)
                    
                    if gene not in os.listdir(family_dir):
                        #the gene's file has not yet been created
                        new_row.to_csv(ospath.join(family_dir, gene))
                    else:
                        #the gene's file has already been created
                        new_row.to_csv(ospath.join(family_dir, gene), mode = 'a', header=False) #mode a to append instead of overwriting the file (mode w)
            
            

def extract_gene(db_name):
    """
    Extracts the gene name from the corresponding fasta file entry
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
    separated = separate(description) #list of words
    (code, genus_name) = treat_word(separated[1]) #first element of seperated is the id, second is the genus

    if code != 2:
         return (code,genus_name)
    
    while code == 2:
         separated.pop(1)
         (code, genus_name) = treat_word(separated[1])
        
    return (code, genus_name)

def extract_binomial(description):
    """
    Extracts the order name from the corresponding fasta file entry. also returns an 
    int code indicating whether the extraction was successful according to the following convention:
    
    0 - contains a . , probably an abbreviation for the genus/species in question, manual treatment needed
    1 - probably a valid bionomial nomenclature
    

    Note
    ----------
    If the integer code is 1, PLEASE DISCARD THE GENUS AND THE SPECIES, they will probably give unpredictable values

    Parameters
    ----------
        Description = fasta formatted description of the gene sequence

    Returns
    -------
        a tuple containing:
        - the integer code (index: 0)
        - the genus (index: 1)
        - the species (index: 2 )
    """
    separated = separate(description) #list of words

    (genus_code, genus_name) = treat_word(separated[1]) #first element of seperated is the id, second is the genus
    
    while genus_code == 2:
        #elimnating keywords from the list and skipping to next word
        separated.pop(1) 
        (genus_code, genus_name) = treat_word(separated[1])
    
    separated.pop(1) #popping the genus to deal with the species now

    (species_code, species_name) = treat_word(separated[1]) #first element of seperated is the id, second is now the species
    
    while species_code == 2:
        #elimnating keywords from the list and skipping to next word
        separated.pop(1) 
        (species_code, species_name) = treat_word(separated[1])
    
    new_code = species_code*genus_code #if there is an error in extracting any of the two, indicate an error
    return(new_code, genus_name, species_name)
    
    

    




def separate(description):
     """
    creates a list of the different words of the description seperated by spaces
    ----------
        Description = fasta formatted description of the gene sequence
    Returns
    -------
    list in question
    """
     return description.split(" ")

def treat_word(word):
     """
    used to treat a word (genus or species) to integrate it in the database. also returns an 
    int code indicating whether the treatment was successful according to the following convention:
    0 - contains a . , probably an abbreviation for the genus/species in question, manual treatment needed
    1 - not a keyword, does not contain a . , probably a valid genus/species name
    2 - keyword (unverified, sp., var., uncultured), probably skip to the next word

    Parameters
    ----------
        word = word to be treated as a string

    Returns
    -------
        a tuple containing:
        - the integer code (index: 0)
        - the word after treatement (index: 1)
    """

     code = 1
     new_word = word.lower()

     strip_list =["[", "]", "(", ")", "'"]
     for symbol in strip_list:
         new_word = new_word.replace(symbol,"")

     skipping_keywords = ["sp.", "var." "uncultured", "unverified"]

     if "." in new_word:
         code = 0
     for keyword in skipping_keywords:
         if keyword in new_word:
             code = 2
     return (code, new_word)
    

def extract_id(description):
     separated = separate(description)
     id = separated[0]
     id = id.replace(">", "")
     return id

"""
def get_family(genus):
     for i in range(0,len(hierarchy)):
         if genus in hierarchy.genuses[i]:
             return hierarchy.Family[i]
     return "undefined"
"""

def get_taxids(species_list):
    taxids = []
    for species_name in species_list:
        handle = Entrez.esearch(db='taxonomy', term=species_name)
        record = Entrez.read(handle)
        if not len(record['IdList']) == 0:
            taxid = record['IdList'][0]
            taxids.append(taxid)
    return taxids

def get_families(taxids):
    handle = Entrez.efetch(db="taxonomy", id=taxids, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    family_names = []
    for record in records:
        try:
            family_name = record['LineageEx'][17]['ScientificName']
            family_names.append(family_name)
        except:
            family_names.append(None)
    return family_names




build_database()