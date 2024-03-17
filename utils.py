from Bio import SeqIO, Entrez
import pandas as pd
import os.path as ospath

Entrez.email = 'ghassan.abboud@epfl.ch'

#-------------------------------
#extract information from ncbi
#-------------------------------

def build_search_term(gene_name, min_len=0, max_len=-1):
    """ generate a search query for NCBI's nucleotide database (GenBank). Bio.Entrez does not provide filters such as sequence length and gene name
    so a search query keywords must be used instead.
    Used to make database queries in extract_database and download_database
    
    Parameters
    ----------
        gene_name: (str) name of the gene of interest
        min_len: (int, optional) lower bound of the sequence length filter, by default 0
        max_len: (int, optional) upper bound of the sequence length filter, by default -1 means no limit
    
    Returns
    ----------
        term: (str) search query
    """
    term = f"{gene_name}[Gene Name]"
    if max_len == -1:
        term += f" NOT 0:{min_len}[Sequence Length]"
    else:
        term += f" AND {min_len}:{max_len}[Sequence Length]"
    return term


def extract_database(gene_name, seqlen_start=0, seqlen_stop=-1):
    '''Request data on GenBank from NCBI (with Entrez module), with filtering options

    Parameters
    ----------
        gene_name: (str) name of the gene to extract from GenBank
        seqlen_start: (int, optional) lower bound of the sequence length filter, by default 0
        seqlen_stop: (int, optional) upper bound of the sequence length filter, by default -1 means no limit

    Returns
    ----------
        seq_record: (str) raw data of all sequences
    '''
    handle = Entrez.esearch(db= "nucleotide", term= build_search_term(gene_name, seqlen_start, seqlen_stop), retmax = 10000)
    record = Entrez.read(handle)
    matching_requests = record["IdList"]
    request_counter = 0
    while len(record["IdList"]) == 10000:
        request_counter += 1
        handle = Entrez.esearch(db= "nucleotide", term= build_search_term(gene_name, seqlen_start, seqlen_stop)
                                , retstart = request_counter*10000,  retmax = 10000)
        record = Entrez.read(handle)
        matching_requests += record["IdList"]
    seq_handle = Entrez.efetch(db="nucleotide", id=matching_requests, retmode = "fasta", rettype = "fasta")
    seq_record = seq_handle.read()
    return seq_record


def download_database(gene_name, seqlen_start=0, seqlen_stop=-1):
    '''Download data to a .fasta file, with filtering options

    Args:
        gene_name: (str) name of the gene to extract from GenBank
        seqlen_start: (int) lower bound of the sequence length filter
        seqlen_stop: (int) upper bound of the sequence length filter

    Returns:
        data_path: (str) path to downloaded data .fasta file
    '''
    data_path = f"{gene_name}[{seqlen_start},{seqlen_stop}].fasta"
    with open(data_path, mode="w") as file:
        file.write(extract_database(gene_name, seqlen_start, seqlen_stop))
    return data_path


def parse_database(path):
    #Create a dictionnary from the database's file
    #Return a dict : 
    #Description of the sequence are the keys
    #DNA sequence are the values

    sequence = {}
    if path.endswith(".fasta"):
        for record in SeqIO.parse(path, "fasta"):
            sequence[record.description] = record.seq
    elif path.endswith(".csv"):
        db = pd.read_csv(path)
        for i in range(len(db.index)):
            species = db.loc[i]["genus"] + " " + db.loc[i]["species"]
            sequence[species] = db.loc[i]["sequence"]
    else:
        raise FileNotFoundError("Invalid file extension, must be .fasta or .csv")
    return sequence


def parse_sequence(path):
    #Put the sequence to align in a list

    sequence = []
    for record in SeqIO.parse(path, "fasta"):
        sequence = record.seq
    return sequence


def replace_nucleotide(seq):
    #Replace all the uncertain nucleotide by '-'

    table =  ["A", "T", "C", "G"]
    seq_clean = ""
    for i in range (len(seq)): 
        if seq[i] in table:
            seq_clean += seq[i]
        elif seq[i] not in table:
            seq_clean += '-'
    return seq_clean

def sequence_uppercase(seq):
    #Put all the nucleotide in uppercase

    seq_up = seq.upper()
    return seq_up


def import_data(database, sequence, option = '1'):
    #Compute the different function to preprocess all the sequence

    if option == '1':
        # Option 1 = Comparison of a DNA sequence with whole the database
        db =  parse_database(database)
        seq = parse_sequence(sequence)
        T = sequence_uppercase(seq)
        T = replace_nucleotide(T)
        
    if option == '2': 
        # Option 2 = Comparison of two DNA sequence
        db =  parse_sequence(database)
        seq = parse_sequence(sequence)
        T = replace_nucleotide(seq)
    return db, T
