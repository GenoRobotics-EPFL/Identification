import os
import shutil
from collections import Counter
from Bio import SeqIO
import multiprocessing
from mash_distance import *
from scipy.cluster.hierarchy import linkage, fcluster
import subprocess
from Bio.Seq import Seq

##file management

def create_folder_and_delete_content(folder_path):
    create_folder(folder_path)
    delete_all_clusters(folder_path)

def create_folder(folder_path):
    """Create a folder if it doesn't exist."""
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"Folder created: {folder_path}")

def delete_all_clusters(output_folder):
    file_list = os.listdir(output_folder)

    """Iterate through the files in the output folder and delete all those that may conflict with the files we generate"""
    for file_name in file_list:
        item_path = os.path.join(output_folder, file_name)
        if file_name.startswith('cluster'):
            if os.path.isfile(item_path):
                os.remove(item_path)
                print(f"Deleted file: {item_path}")
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
                print(f"Deleted folder: {item_path}")
        elif file_name.startswith('group'):
            if os.path.isfile(item_path):
                os.remove(item_path)
                print(f"Deleted file: {item_path}")
        elif file_name.startswith('primalscheme'):
            if os.path.isdir(item_path):
                shutil.rmtree(item_path)
                print(f"Deleted folder: {item_path}")
        elif file_name.startswith('iteration'):
            if os.path.isdir(item_path):
                shutil.rmtree(item_path)
                print(f"Deleted folder: {item_path}")


##clusters

def generate_clusters(sequences, kmers_nucleotides, distance_within_clusters, output_folder, nb_clusters_generated):
    """Generate clusters based on the given sequences and k-mers nucleotides.

    Args:
        sequences (list): list of sequences remaining.
        kmers_nucleotides (list): list k-mers nucleotides for sequences remaining.
        distance_within_clusters (float): The distance threshold for clustering.
        output_folder (str): The output folder path.
        nb_clusters_generated (int): The number of clusters generated.

    """

    dist_matrix = compute_mash_distance_matrix(kmers_nucleotides)
    linkage_matrix = linkage(dist_matrix, 'centroid')
    

    clusters = fcluster(linkage_matrix, t=distance_within_clusters, criterion='distance')

    save_largest_clusters(sequences, clusters, output_folder, nb_clusters_generated)



def save_largest_clusters(sequences, clusters, output_folder, nb_clusters_generated=3):
    """
    arrange the sequences according to the 3 largest groups in the clusters file
    Args:
        sequences (array): sequences of the dataset, unaligned
        clusters (array): array assigning each sequence of the dataset to a cluster number
        output_folder (string): path of output folder
        nb_clusters_generated (integer): number of clusters to generate (default 3)
    """
    cluster_counts = Counter(clusters)
    for cluster_index, (number, _) in enumerate(cluster_counts.most_common(nb_clusters_generated)):
        filter_and_save_cluster(sequences, clusters, number, cluster_index, output_folder)


    
def filter_and_save_cluster(sequences, clusters, cluster_number, cluster_index, output_folder):
    """
        filter the sequences according to a cluster number and save it in a new fasta file
    Args:
        sequences (array): sequences of the dataset, unaligned
        clusters (array): array assigning each sequence of the dataset to a cluster number
        cluster_number (integer): number of the current cluster
        cluster_index (integer): index of the cluster output file (only output file)
    """
    filtered_sequences = [record for record, cluster in zip(sequences, clusters) if cluster == cluster_number]
    print(f"cluster {cluster_index} has  {len(filtered_sequences)} sequences \n")
    SeqIO.write(
        filtered_sequences,
        f"{output_folder}/cluster{cluster_index}.fasta",
        "fasta",
    )




##stats

def get_mean_and_deviation(filename_in):
    """Calculate the mean and standard deviation of sequence lengths from the given file.

    """
    alignements = list(SeqIO.parse(filename_in, "fasta"))
    num_sequences = len(alignements)

    lengths = list(map(len, alignements))
    length_mean = sum(lengths) / len(lengths) 
    standard_deviation = sum(abs(i - length_mean) for i in lengths) / len(lengths) 

    print(f"Length stats for cluster of path {filename_in} :\n" + 
        f"Mean : {length_mean}\n" + f"Standard deviation : {standard_deviation}")

def get_number_pairs_found(pair_dictionnary):
    return sum(len(pairs) for amplicon_range, pairs in pair_dictionnary.items())

def percentage_round(number):
    return round(number*100, 2)

##alignement

def run_clustal_command(filename_in, filename_out):
    """Run the Clustal Omega command to perform multiple sequence alignment.

    """
    clustalo_command = [
    "clustalo",
    "-infile", filename_in,
    "--iterations", "1",
    "-align",
    "--threads", str(multiprocessing.cpu_count()),
    "-outfile", filename_out,
    "--outfmt", "fasta"
    ]

    if not os.path.exists(filename_out):
        subprocess.run(clustalo_command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)



##regex/coverage

def produce_regex_with_mismatches(primer, mismatches_allowed):
    return f"({primer}){{s<={mismatches_allowed}}}"

def generate_regex(primerL, primerR, allowed_mismatches=6):
    primer_reversed = str(Seq(primerR.upper()).reverse_complement())
    primer_forward = primerL.upper()
    return produce_regex_with_mismatches(primer_forward, allowed_mismatches), produce_regex_with_mismatches(primer_reversed, allowed_mismatches)

