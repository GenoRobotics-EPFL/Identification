
from Bio.Seq import Seq
import os
from collections import defaultdict
from file_utils import *
from mash_distance import *
from primer_extraction import *
import concurrent.futures
import argparse


def main(args):
    create_folder(args.output_folder)
    delete_all_clusters(args.output_folder)

    records = list(SeqIO.parse(args.filename_fasta, "fasta"))

    generate_clusters(records, args.k_length, args.distance_within_clusters, args.output_folder, args.number_clusters)

    for i in range(1, args.number_clusters + 1):
        print(f"\nAligning cluster {i}")
        run_clustal_command("{}/cluster{}.fasta".format(args.output_folder, i),
                            "{}/cluster{}.aln".format(args.output_folder, i))

    pairs_dictionary = defaultdict(list)

    with concurrent.futures.ThreadPoolExecutor() as executor:
        executor.map(lambda j: generate_primer_for_cluster(j, args.output_folder, pairs_dictionary),
                     range(1, args.number_clusters + 1))

    test_primer_coverage_against_dataset([r.seq for r in records], pairs_dictionary)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Your script description')
    parser.add_argument('-o', '--output-folder', type=str, required=True, help='Output folder path')
    parser.add_argument('-f', '--filename-fasta', type=str, required=True, help='path to fasta file')
    parser.add_argument('-k', '--k-length', type=int, default=15, help='Length of k-mers')
    parser.add_argument('-d', '--distance-within-clusters', type=float, default=0.5, help='Distance within clusters')
    parser.add_argument('-n', '--number-clusters', type=int, default=3, help='Number of clusters')

    args = parser.parse_args()
    main(args)