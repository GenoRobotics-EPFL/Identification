

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
    delete_all_clusters(args.output_folder)
    csv_file_path = f"{args.output_folder}/output.csv"
    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)

        # Write the header
        csv_writer.writerow(['Range', 'Forward Primer Seq', 'Reverse Primer Seq', 'Coverage [%]'])

        # Write the data
        for key, value in result_pairs.items():
            for pair in value:
                forward_primer_seq = pair[0][0]['seq']
                reverse_primer_seq = pair[0][1]['seq']
                coverage = round(pair[1]*100, 2)
                csv_writer.writerow([key, forward_primer_seq, reverse_primer_seq, coverage])

    print(f"Data has been written to {csv_file_path}")

def main(args):
    create_folder(args.output_folder)
    delete_all_clusters(args.output_folder)
    indices_covered = set()
    result_pairs = defaultdict(list)
    alignements = list(SeqIO.parse(args.filename_fasta, "fasta"))
    nb_total_records = len(alignements)
    kmers_nucleotides = retrieve_kmers(alignements, args.k_length)
    for iteration_steps in range(1, args.number_pairs + 1):

        print(f"\n\nStep {iteration_steps}\n")
        output_folder = f"{args.output_folder}/iteration{iteration_steps}"
        create_folder(output_folder)
        records = [alignements[i] for i in range(nb_total_records) if i not in indices_covered]
        kmers_filtered = [kmers_nucleotides[i] for i in range(nb_total_records) if i not in indices_covered]

        generate_clusters(records, kmers_filtered, args.distance_within_clusters, output_folder, args.number_clusters)


        for i in range(1, args.number_clusters + 1):
            print(f"Aligning cluster {i}\n")
            run_clustal_command("{}/cluster{}.fasta".format(output_folder, i),
                                "{}/cluster{}.aln".format(output_folder, i))

        pairs_dictionary = defaultdict(list)

        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.map(lambda j: generate_primer_for_cluster(j, output_folder, pairs_dictionary),
                        range(1, args.number_clusters + 1))
        new_indices, best_pair = test_primer_coverage_against_dataset([r.seq for r in alignements], pairs_dictionary, indices_covered, result_pairs)
        for k, v in best_pair.items():
            result_pairs[k].append(v)
        indices_covered.update(new_indices)
        print(f"Total covered : {len(indices_covered)} out of {nb_total_records} records, total coverage : {round(100*len(indices_covered)/nb_total_records, 2)}%")
    write_csv(result_pairs, args)


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