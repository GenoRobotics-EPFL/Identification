

import os
from collections import defaultdict
from file_utils import *
from mash_distance import *
from primer_extraction import *
import concurrent.futures
import argparse
from Bio.Seq import Seq
import csv

def write_csv(result_pairs, args):
    """Write the result pairs to a CSV file.

    Args:
        result_pairs (dict): A dictionary containing the result pairs of primers
        args: The arguments object.
    """
    delete_all_clusters(args.output_folder) ##delete all clusters generate by primalschem and alignement to leave only the output csv file
    csv_file_path = os.path.join(args.output_folder, "output.csv")

    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)

        ##write csv header
        csv_writer.writerow(['Amplicon Size Range', 'Forward Primer Seq', 'Reverse Primer Seq', 'Coverage [%]'])
        for amplicon_range, value in result_pairs.items():
            for primer_tuple, primer_coverage in value:
                write_pair_of_primers(amplicon_range, primer_tuple, primer_coverage, csv_writer)

    print(f"Data has been written to {csv_file_path}")

def write_pair_of_primers(amplicon_range, primer_tuple, primer_coverage, csv_writer):
    """Write a pair of primers to a CSV file.

    Args:
        amplicon_range (str): The range of the amplicon.
        primer_tuple (tuple): A tuple containing the forward and reverse primer sequences.
        coverage (float): The coverage percentage of the primer pair against the dataset.
        csv_destination_file (str): The path to the output CSV file.

    """
    forward_primer_seq = primer_tuple[0]['seq']
    reverse_primer_seq = primer_tuple[1]['seq']
    primer_coverage = percentage_round(primer_coverage)
    csv_writer.writerow([amplicon_range, forward_primer_seq, reverse_primer_seq, primer_coverage])

def main(args):
    """Run the main process for primer design.

    Args:
        args: The arguments object (output_folder, fasta_path, k_length, distance_allowed_within_clusters, umber_clusters generated, number_pairs returned, min and max amplicon range)

    """
    ##prepare output folder
    create_folder_and_delete_content(args.output_folder)


    indices_covered = set() ##indices of sequences covered by found primers
    result_pairs = defaultdict(list) ##result pairs have shape {amplicon size range : [((primer_forward, primer_reverse), coverage)]}

    ##create dataset and kmers
    all_sequences = list(SeqIO.parse(args.filename_fasta, "fasta"))
    nb_total_sequences = len(all_sequences)
    kmers_nucleotides = retrieve_kmers(all_sequences, args.k_length)


    for iteration_steps in range(1, args.number_pairs + 1):
        ##create folder for iteration
        print(f"\n\n Step {iteration_steps}\n")
        iteration_folder = os.path.join(args.output_folder, f"iteration{iteration_steps}")
        create_folder_and_delete_content(iteration_folder)

        ##filter sequences and kmers if already covered by primers
        filtered_sequences = [all_sequences[i] for i in range(nb_total_sequences) if i not in indices_covered]
        kmers_filtered = [kmers_nucleotides[i] for i in range(nb_total_sequences) if i not in indices_covered]

        ##generate clusters according to kmers and mash distance
        generate_clusters(filtered_sequences, kmers_filtered, args.distance_within_clusters, iteration_folder, args.number_clusters)

        ##align each cluster using clustal
        align_all_clusters(iteration_folder, args.number_clusters)

        ##new indices covered, best pair
        new_indices, best_pair = generate_and_test_primers(iteration_folder, args.number_clusters, indices_covered, result_pairs, all_sequences)
        
        for k, v in best_pair.items():
            result_pairs[k].append(v)
            
        indices_covered.update(new_indices)
        print(f"Total covered : {len(indices_covered)} out of {nb_total_sequences} sequences, total coverage : {percentage_round(len(indices_covered)/nb_total_sequences)}%")
    write_csv(result_pairs, args)

def generate_and_test_primers(iteration_folder, number_clusters, indices_covered, found_primers, all_sequences):
    new_pairs = defaultdict(list)
    with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.map(lambda j: generate_primer_for_cluster(j, iteration_folder, new_pairs), range(number_clusters))

    return get_primer_coverage_against_dataset([r.seq for r in all_sequences], new_pairs, indices_covered, found_primers)

def align_all_clusters(iteration_folder, number_clusters):
    for i in range(args.number_clusters):
            print(f"Aligning cluster {i}\n")
            run_clustal_command(
                f"{iteration_folder}/cluster{i}.fasta",
                f"{iteration_folder}/cluster{i}.aln",
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Your script description')
    parser.add_argument('-o', '--output-folder', type=str, default="output", help='Output folder path')
    parser.add_argument('-f', '--filename-fasta', type=str, required=True, help='path to fasta file')
    parser.add_argument('-k', '--k-length', type=int, default=15, help='Length of k-mers')
    parser.add_argument('-d', '--distance-within-clusters', type=float, default=0.1, help='Distance within clusters')
    parser.add_argument('-n', '--number-clusters', type=int, default=5, help='Number of clusters')
    parser.add_argument('-p', '--number-pairs', type=int, default=3, help='Number of pairs returned')
    parser.add_argument('-min_range', '--min_range', type=int, default=3, help='Minimum amplicon size for primer search')
    parser.add_argument('-max_range', '--max_range', type=int, default=3, help='Maximum amplicon size for primer search')

    args = parser.parse_args()
    main(args)