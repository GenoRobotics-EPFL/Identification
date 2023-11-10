from Bio import Entrez
from Bio import SeqIO
import pandas as pd

Entrez.email = 'ghassan.abboud@epfl.ch'



def build_search_term(gene_name, min_len = 0, max_len= -1):
    term = f"{gene_name}[Gene Name]"
    if max_len == -1:
        for i in range(min_len):
            term += f" NOT {i}[Sequence Length]"
    else:
        term += f" AND ({min_len}[Sequence Length]"
        for i in range(min_len +1, max_len+1):
            term += f" OR {i}[Sequence Length]"
        term += ")"
    return term

def extract_database(gene_name, seqlen_start = 0, seqlen_stop = -1):
    '''Request data on GenBank from NCBI (with Entrez module), with filtering options
    
    Args:
        gene_name: (str) name of the gene to extract from GenBank
        seqlen_start: (int) lower bound of the sequence length filter
        seqlen_stop: (int) upper bound of the sequence length filter
    
    Returns:
        seq_record: (str) raw data of all sequences'''
    handle = Entrez.esearch(db= "nucleotide", term= build_search_term(gene_name, seqlen_start, seqlen_stop))
    record = Entrez.read(handle)
    matching_requests = record["IdList"]
    seq_handle = Entrez.efetch(db="nucleotide", id=matching_requests, retmode = "fasta", rettype = "fasta")
    seq_record = seq_handle.read()
    return seq_record

def download_database(gene_name, seqlen_start = 0, seqlen_stop = -1):
    '''Download data to a .fasta file, with filtering options
    
    Args:
        gene_name: (str) name of the gene to extract from GenBank
        seqlen_start: (int) lower bound of the sequence length filter
        seqlen_stop: (int) upper bound of the sequence length filter
    
    Returns:
        data_path: (str) path to downloaded data .fasta file'''
    data_path = f"{gene_name}[{seqlen_start},{seqlen_stop}].fasta"
    with open(f"{gene_name}[{seqlen_start},{seqlen_stop}].fasta", mode = "w") as file:
        file.write(extract_database(gene_name,seqlen_start, seqlen_stop))
    return data_path

def parse_data(path):
    '''Creates a dataframe (pandas) from raw data
    
    Args:
        path: (str) path to raw data file
        
    Returns:
        df: pandas dataframe of sequences'''
    # Parse the FASTA file into SeqRecord object
    sequences = SeqIO.parse(path, "fasta")

    sequence_data = []
    r = []
    # Extract from each record a dictionnary of features
    for record in sequences:
        r.append(record)
        sequence_data.append({'ID': record.id,
                              'Sequence': str(record.seq),
                              'Description': record.description})

    # Transform the list of dictionnaries to a pandas DataFrame
    df = pd.DataFrame(sequence_data)
    return df


def calculate_avg_length(db):
    #return avg_length
    a=3
def calculate_unverified_percentage(db):
    #lol
    a=5