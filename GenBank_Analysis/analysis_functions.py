from Bio import Entrez
from Bio import SeqIO
import matplotlib.pyplot as plt
import os.path as ospath
import os
import numpy as np

Entrez.email = 'ghassan.abboud@epfl.ch'


def calculate_unverified_percentage(path_to_fasta):
    seq_objects = SeqIO.parse(path_to_fasta, "fasta")  # converts to a seq object
    gene_amount = 0
    unverified_number = 0
    for genes in seq_objects:
        gene_amount += 1
        if "UNVERIFIED" in genes.description:
            unverified_number += 1

    unverified_percentage = unverified_number / gene_amount * 100

    return unverified_percentage


def show_length_graph(path_to_fasta):
    # takes the path to the fasta file as a parameter
    seq_objects = SeqIO.parse(path_to_fasta, "fasta")  # converts to a seq object
    sizes_and_recurrence = {}
    for genes in seq_objects:  # creating a dict with key:size and value:recurrence amount
        if len(genes) in sizes_and_recurrence:
            sizes_and_recurrence[len(genes)] += 1
        else:
            sizes_and_recurrence[len(genes)] = 1
    sizes_and_recurrence = sorted([(k, v) for k, v in sizes_and_recurrence.items()])  # sorting the data
    sizes = [k for k, v in sizes_and_recurrence]  # sizes
    recurrences = [v for k, v in sizes_and_recurrence]  # recurrences

    fig, ax = plt.subplots()
    ax.bar(sizes, recurrences)
    ax.set_title('Distribution of Lengths of Genes')
    plt.xlabel("Lengths of Genes")
    plt.ylabel("Amount of Genes")
    plt.show()

#show_length_graph("/Users/ilgazarslan/Documents/Genorobotics/Identification/Database/matk.fasta")


# does the same thing as length grapher but not efficient at all
def length_grapher_histogram(path_to_fasta):
    gene_lengths_amount = 0
    gene_lengths = []
    seq_objects = SeqIO.parse(path_to_fasta, "fasta")
    for genes in seq_objects:
        gene_lengths_amount += 1
        gene_lengths.append(len(genes))

    plt.hist(gene_lengths, bins=gene_lengths_amount, color='skyblue', edgecolor='black')
    plt.xlabel('Gene Length')
    plt.ylabel('Frequency')
    plt.title('Gene Length Distribution')
    plt.show()


def calculate_avg_length(path_to_fasta):
    seq_objects = SeqIO.parse(path_to_fasta, "fasta")
    number_of_genes = 0
    total_length = 0
    for genes in seq_objects:
        number_of_genes += 1
        total_length += len(genes)
    average = total_length / number_of_genes
    return average

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
