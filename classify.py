#tasks:
#find a way to identify the family from the description
#attribute of the Bio.seq.seq object.

#ask how to organize the databases:
# - organize csv files by family -> code isnt provided in utils,
# one file of hundreds of csv files for each family
# - one csv file with all genes of all families
# and selection of family happens inside code

genus_to_family = {}

import pandas as pd
from Bio import SeqIO
import os.path as ospath
import os
import time
from Bio import Entrez

Entrez.email = 'eva.frossard@epfl.ch'

def build_database():
    """
    Creates a new folder that organizes the Genorobotics database by family.
    The folder is called Database_by_family.
    This folder contains:
    - a csv with the gene entries that were not treated properly and need manual intervention
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
        #this for loop has 4 iterations, one for each of the four gene files
        print("begin: " + file)
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
                #print(genus+ " " + species)
                family = get_family(genus + " " + species)
                
            
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
        
        print("end: " + file)    
            

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
    if db_name == "sequences_matK_800-1550.fasta" or db_name == "matk.fasta":
        return "matk.csv"
    elif db_name == "psbA-trnH_sequence.fasta" or db_name == "psba-trnh.fasta":
        return "psba-trnh.csv"
    elif db_name == "rcbL_sequence.fasta" or db_name == "rbcl.fasta":
        return "rbcl.csv"
    elif db_name == "its" or db_name == "its.fasta":
        return "its.csv"
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

     skipping_keywords = ["sp.", "var." "uncultured", "unverified", "cf.", "aff."]

     if "." in new_word:
         code = 0
     for keyword in skipping_keywords:
         if keyword in new_word:
             code = 2
     return (code, new_word)
    

def extract_id(description):
    """
    Extracts the unique NCBI id from the description of a gene sequence entry on GenBank.

    Parameters
    ----------
        Description = fasta formatted description of the gene sequence

    Returns
    -------
        the ID
    """
    separated = separate(description)
    id = separated[0]
    id = id.replace(">", "")
    return id

def get_taxids(binomial_list):
    """
    Accesses the NCBI taxonomy database through the Entrez library and retrieves the corresponding taxonomy ids
    from a list of species name. The names are given using the binomial nomenclature (genus + species).
    If a name is not found in the taxonomy, the taxid is set to -1

    Parameters
    ----------
        binomial_list: a list with the species' names

    Returns
    -------
        a list of the corresponding taxids
    """
    taxids = []
    for binomial in binomial_list:
        handle = Entrez.esearch(db='taxonomy', term=binomial)
        record = Entrez.read(handle)
        if not len(record['IdList']) == 0:
            taxid = record['IdList'][0]
            taxids.append(taxid)
        else:
            taxid = -1
            taxids.append(taxid)

    return taxids

def get_families(taxids):
    """
    Accesses the NCBI taxonomy database through the Entrez library and retrieves the corresponding families
    from a list of taxids.
    If a name is not found in the taxonomy, the family's name is replaced by "Not Found"

    Parameters
    ----------
        taxids: a list of taxids

    Returns
    -------
        a list of the corresponding family names
    """
    handle = Entrez.efetch(db="taxonomy", id=taxids, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    family_names = []
    for record in records:
        try:
            family_name = treat_lineage_ex(record['LineageEx'])
            #print(type(record['LineageEx']))
            family_names.append(family_name)
        except:
            family_names.append("Not Found")
    return family_names


def treat_lineage_ex(lineage):
    """
    Extracts the family from a Bio.Entrez.Parser.ListElement object corresponding to the search result from
    looking through the NCBI taxonomy database. 

    Parameters
    ----------
    lineage = 'Bio.Entrez.Parser.ListElement' object, is the 'LineageEx' element of the search result.


    Returns
    -------
    the name of the family

    """
    family = "Not Found"
    for level in lineage:
        if level['Rank'] == "family":
            family = level['ScientificName']
    return family

def get_family(binomial):
    """
    get the family of a species from its binomial nomenclature by searching through the NCBI taxonomy db.
    Returns "Not Found" if the family was not found in the db. This can happen in three cases:

    1- the binomial passed as argument is not valid, might indicate a problem in extracting the binomial
    from the fasta description, (get_taxids(binomial[0] returns -1)
    2- a taxid for the binomial was found, but no corresponding desciption of lineage.
    3- a description of lineage was found but no family member. 

    It is still unclear Whether cases 2 and 3 are possible since the taxonomy db should have a description lineage
    for each taxid it can provide.

    Parameters
    ----------
    binomial = binomial nomenclature of species (genus + species)


    Returns
    -------
    the name of the family
    
    """
    if binomial[0] in genus_to_family.keys():
        return genus_to_family[binomial[0]]
    else:
        taxid = get_taxids([binomial])[0]
        family = get_families([taxid])[0]
        genus_to_family[binomial[0]] = family
        return family

    


start=time.time()
build_database()
end = time.time()
print(f'Temps d\'exécution : ', end-start, 's')