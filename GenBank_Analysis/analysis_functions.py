from Bio import Entrez
from Bio import SeqIO

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
    handle = Entrez.esearch(db= "nucleotide", term= build_search_term(gene_name, seqlen_start, seqlen_stop))
    record = Entrez.read(handle)
    matching_requests = record["IdList"]
    seq_handle = Entrez.efetch(db="nucleotide", id=matching_requests, retmode = "fasta", rettype = "fasta")
    seq_record = seq_handle.read()
    return seq_record

def download_database(gene_name, seqlen_start = 0, seqlen_stop = -1):
        with open(f"{gene_name}[{seqlen_start},{seqlen_stop}].fasta", mode = "w") as file:
            file.write(extract_database(gene_name,seqlen_start, seqlen_stop))


def calculate_avg_length(db):
    #return avg_length
    a=3
def calculate_unverified_percentage(db):
    #lol
    a=5