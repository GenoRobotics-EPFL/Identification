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


##clusters

def generate_clusters(records, k_length, distance_within_clusters, output_folder, NUMBER_CLUSTERS):
    kmers_nucleotides = retrieve_kmers(records, k_length)


    dist_matrix = compute_mash_distance_matrix(kmers_nucleotides)
    linkage_matrix = linkage(dist_matrix, 'centroid')
    

    clusters = fcluster(linkage_matrix, t=distance_within_clusters, criterion='distance')

    save_largest_clusters(records, clusters, output_folder, NUMBER_CLUSTERS)



def save_largest_clusters(records, clusters, output_folder, number_clusters=3):
    """_summary_
    arrange the records according to the 3 largest groups in the clusters file
    Args:
        records (array): sequences of the dataset, unaligned
        clusters (array): array assigning each sequence of the dataset to a cluster number
        output_folder (string): path of output folder
        number_clusters (integer): number of clusters to generate (default 3)
    """
    cluster_counts = Counter(clusters)
    cluster_index = 1
    for number, _ in cluster_counts.most_common(number_clusters):
        filter_and_save_cluster(records, clusters, number, cluster_index, output_folder)

        cluster_index += 1


    
def filter_and_save_cluster(records, clusters, cluster_number, cluster_index, output_folder):
    """_summary_
        filter the records according to a cluster number and dave it in a new fasta file
    Args:
        records (array): sequences of the dataset, unaligned
        clusters (array): array assigning each sequence of the dataset to a cluster number
        cluster_number (integer): number of the current cluster
        cluster_index (integer): index of the cluster output file (only output file)
    """
    filtered_records = [record for record, cluster in zip(records, clusters) if cluster == cluster_number]
    print("cluster {} has  {} records".format(cluster_number, len(filtered_records)))
    SeqIO.write(filtered_records, "{}/cluster{}.fasta".format(output_folder, cluster_index), "fasta") ##save into file




##stats

def get_mean_and_deviation(filename_in):
    alignements = list(SeqIO.parse(filename_in, "fasta"))
    num_sequences = len(alignements)

    lengths = list(map(len, alignements))
    length_mean = sum(lengths) / len(lengths) 
    standard_deviation = sum(abs(i - length_mean) for i in lengths) / len(lengths) 

    print(f"Length stats for cluster of path {filename_in} :\n" + 
        f"Mean : {length_mean}\n" + f"Standard deviation : {standard_deviation}")

def get_number_pairs_found(pair_dictionnary):
    #pairs of primers found
    number_pairs = 0
    for amplicon_range, pairs in pairs_dictionnary.items():
        number_pairs += len(pairs)
    print(f"We have found {number_pairs} pairs of primers")

##alignement

def run_clustal_command(filename_in, filename_out):

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
        subprocess.run(clustalo_command, check=True, stdout=subprocess.DEVNULL)



##regex/coverage

def produce_regex_with_mismatches(primer, mismatches_allowed):
    return f"({primer}){{s<={mismatches_allowed}}}"

def generate_regex(primerL, primerR, allowed_mismatches):
    primer_reversed = Seq(primerR.upper()).reverse_complement()
    primer_forward = primerL.upper()
    return produce_regex_with_mismatches(primer_forward, allowed_mismatches), produce_regex_with_mismatches(primer_reversed, allowed_mismatches)

