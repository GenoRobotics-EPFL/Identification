from Bio import Entrez
from Bio import SeqIO, AlignIO
import matplotlib.pyplot as plt
import os.path as ospath
import os
import numpy as np
import pandas as pd
from difflib import SequenceMatcher


#-------------------------------
#parse data
#-------------------------------


def parse_data(path):
    '''Creates a dataframe (pandas) from raw data

    Args:
        path: (str) path to raw data file

    Returns:
        df: pandas dataframe of sequences'''
    # Parse the FASTA file into SeqRecord object
    sequences = SeqIO.parse(path, "fasta")

    sequence_data = []
    # Extract from each record a dictionnary of features
    for record in sequences:
        sequence_data.append({'ID': record.id,
                              'Sequence': str(record.seq),
                              'Description': record.description.lower()})

    # Transform the list of dictionnaries to a pandas DataFrame
    df = pd.DataFrame(sequence_data)
    return df


#-------------------------------
#calculate database characteristics
#-------------------------------


def calculate_unverified_percentage(path_to_fasta):
    """calculate the percentage of sequences in a fasta file that are labeled 'unverified' in their
    description, this indicates that they were not checked through the translation to protein.
    
    Parameters
    ----------
    path_to_fasta: (str) path to the fasta file
    Returns
    ----------
    unverified_percentage: (float) percentage of unverfied sequences
    """
    seq_objects = SeqIO.parse(path_to_fasta, "fasta")  # converts to a seq object
    gene_amount = 0
    unverified_number = 0
    for gene in seq_objects:
        gene_amount += 1
        if "unverified" in gene.description.lower():
            unverified_number += 1

    unverified_percentage = unverified_number / gene_amount * 100

    return unverified_percentage


def show_length_graph(path_to_fasta):
    # takes the path to the fasta file as a parameter
    seq_objects = SeqIO.parse(path_to_fasta, "fasta")  # converts to a seq object
    sizes_and_recurrence = {}
    for gene in seq_objects:  # creating a dict with key:size and value:recurrence amount
        if len(gene.seq) in sizes_and_recurrence:
            sizes_and_recurrence[len(gene.seq)] += 1
        else:
            sizes_and_recurrence[len(gene.seq)] = 1

    sizes_and_recurrence = sorted([(k, v) for k, v in sizes_and_recurrence.items()])  # sorting the data
    sizes = [k for k, v in sizes_and_recurrence]  # sizes
    recurrences = [v for k, v in sizes_and_recurrence]  # recurrences

    file_name = os.path.basename(path_to_fasta)  # name of the file

    fig, ax = plt.subplots()
    ax.bar(sizes, recurrences)
    ax.set_title(f'Distribution of Lengths of Genes - {file_name}')
    plt.xlabel("Length of Genes")
    plt.ylabel("number of GenBank entries")
    plt.show()


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
    """
    calculate average length of the sequences in a fasta file
    Parameters
    ----------
        path_to_fasta: (str) path of the fasta file
    Returns
    ----------
        average: (float) average length in bp
    """
    seq_objects = SeqIO.parse(path_to_fasta, "fasta")
    number_of_genes = 0
    total_length = 0
    for genes in seq_objects:
        number_of_genes += 1
        total_length += len(genes)
    average = total_length / number_of_genes
    return average


def nucleotide_ratio(path):
    """ calculate the ratio of the four nucleotides for each sequence of a given DataFrame and add the ratios as
    additional columns of the Dataframe
     
    Parameters
    ----------
        df: gene DataFrame
    Returns
    ----------
        Nucleotide ration of all sequences """
    df = parse_data(path)
    df['A'] = df.Sequence.apply(lambda x: x.count('A') / len(x) * 100)
    df['T'] = df.Sequence.apply(lambda x: x.count('T') / len(x) * 100)
    df['G'] = df.Sequence.apply(lambda x: x.count('G') / len(x) * 100)
    df['C'] = df.Sequence.apply(lambda x: x.count('C') / len(x) * 100)
    ratio = df[['A', 'T', 'G', 'C']]
    return ratio


# a general function which goes through every gene file in a directory :
def database_file_analysis(directory, method_12345):
    ''' creates a list containing every path to the file in the directory,
        then goes through every gene with the indicated function as the second parameter:

        1 : calculate_unverified_percentage
        2 : show_length_graph
        3 : calculate_avg_length
        4 : nucleotide_ratio
        5 : ????? the last consensus function which is not finished yet
                                                                    '''

    file_paths = []
    file_names = []

    # Walk through the directory tree using os.walk
    for root, dirs, files in os.walk(directory):
        for file in files:
            # Construct the full path of the file
            file_path = os.path.join(root, file)
            file_paths.append(file_path)
            file_names.append(file)

    information = []
    n = len(file_paths)

    if method_12345 == 1:
        print("The percentage of UNVERIFIED sequences are :")
        for i in range(n):
            information.append(calculate_unverified_percentage(file_paths[i]))
            print("\t", file_names[i], " : ", information[i])

    elif method_12345 == 2:
        for i in range(n):
            show_length_graph(file_paths[i])

    elif method_12345 == 3:
        print("The average length of genes are :")
        for i in range(n):
            information.append(calculate_avg_length(file_paths[i]))
            print("\t", file_names[i], " : ", information[i])

    elif method_12345 == 4:
        print("The percentage of the nucleotides in each sequence are :")
        for i in range(n):
            print("\t", file_names[i], " : \n", nucleotide_ratio(file_paths[i]))


def get_consensus_sequence(fasta_path):
    # find a efficient way of PROPERLY making consensus
    # chatgpt is useless in this , obviously gotta use biopython + other libraries
    print("idk")


# a way of calculating the similarity between sequences for the consensus similarity function:
# Jaccard distance --> the measure unit (the similarity between two sequences)
def calculate_jaccard_distance(sequence1, sequence2):
    matcher = SequenceMatcher(None, sequence1, sequence2)
    similarity = matcher.ratio()
    return 1 - similarity  # Convert similarity to distance


# function which calculates the jaccard distance for every sequence:
def calculate_distances_to_consensus(fasta_path, consensus_sequence):
    sequences = [str(seq.seq) for seq in SeqIO.parse(fasta_path, "fasta")]
    distances = [calculate_jaccard_distance(seq, consensus_sequence) for seq in sequences]
    return distances


# a function to plot these distances:
def plot_scatter(distances, output_filename="scatter_plot.png"):
    plt.figure(figsize=(10, 6))
    plt.scatter(range(1, len(distances) + 1), distances, color='blue', alpha=0.7)
    plt.title('Jaccard Distances to Consensus Sequence')
    plt.xlabel('Sequence Index')
    plt.ylabel('Jaccard Distance')
    # plt.savefig(output_filename)
    plt.show()

# distances = calculate_distances_to_consensus(fasta_path, consensus_sequence)
# plot_scatter(distances)


database_file_analysis(os.getcwd()+ "\\Database", 3)
