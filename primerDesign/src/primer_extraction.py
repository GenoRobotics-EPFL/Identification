import os
import subprocess
import pandas as pd
from file_utils import *
import regex
import primer3
from collections import defaultdict
from cross_dimer import has_dimer


"""
returns True if the primer proposed has not conflict with the pairs already found
"""
def check_bindings(primer_forward, primer_reverse, pairs_found):
    """
    Check if the proposed primer has conflicts with the pairs already found.

    Args:
        primer_forward (str): The forward primer sequence.
        primer_reverse (str): The reverse primer sequence.
        pairs_found (dict): Dictionary of pairs already found.

    Returns:
        bool: True if the primer has no conflicts, False otherwise.
    """
    return not any(
        has_dimer(primer_forward, pair[0].seq) or
        has_dimer(primer_forward, pair[1].seq) or
        has_dimer(primer_reverse, pair[0].seq) or
        has_dimer(primer_reverse, pair[1].seq)
        for pairs in pairs_found.values() for pair, _ in pairs
    ) and not has_dimer(primer_forward, primer_reverse)

def generate_primer_for_cluster(cluster_index, output_folder, pairs_dictionnary, min_range=300, max_range=500):
    """
    Generate primers for a cluster.

    Args:
        cluster_index (int): The index of the cluster.
        output_folder (str): The output folder path.
        pairs_dictionnary (dict): Dictionary to store the primer pairs.
        min_range (int): The minimum range for amplicon.
        max_range (int): The maximum range for amplicon.

    """

    cluster_aln = f"{output_folder}/cluster{cluster_index}.aln"
    # Design primers
    for j in range(min_range, max_range, 80): ##amplicon between 300 and 1000 seems a reasonable range
        command = [
            "primalscheme", "multiplex",
            "-a", str(j),
            "-a", str(j + 80),
            "-f", cluster_aln,
            "-o", f"{output_folder}/primalscheme_results{cluster_index}"
        ]
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # Meant to be used on a Linux terminal

        primer_results_path = f"{output_folder}/primalscheme_results{cluster_index}/scheme.primer.tsv"
        if os.path.exists(primer_results_path):
            with open(primer_results_path) as f:
                data = pd.read_csv(f, sep='\t')
                pair_key = f"{j}-{j + 80}"
                pairs_dictionnary[pair_key].extend(
                    (data.iloc[l], data.iloc[l + 1]) for l in range(0, len(data), 2)
                )

def get_primer_coverage_against_dataset(all_sequences, new_primer_pairs, indices_covered, pairs_found):
    """
    Calculate the coverage of primer pairs against the dataset of sequences.

    Args:
        all_sequences (list): List of all sequences in dataset.
        new_primer_pairs (dict): Dictionary of new primer pairs.
        indices_covered (list): List of indices already covered by primers.
        pairs_found (dict): Dictionary of pairs already found.

    Returns:
        tuple: A tuple containing the best covered indices and the best pair.

    """

    nucleotides = [seq for i, seq in enumerate(all_sequences) if i not in indices_covered]
    new_primer_pairs = {k: [pair for pair in v if check_bindings(pair[0].seq, pair[1].seq, pairs_found)]
                    for k, v in new_primer_pairs.items()}

    allowed_mismatches = 6
    best_covered_indices = []
    max_coverage = {}
    best_pair = {}

    if not new_primer_pairs:
        print("\nno additional primers found \n")
        return best_covered_indices, {}

    for amplicon_range, pairs in new_primer_pairs.items():
        max_coverage[amplicon_range] = 0.0
        for pair in pairs:
            primerR = pair[1].seq
            primerL = pair[0].seq
            left_regex, right_regex = generate_regex(primerL, primerR, allowed_mismatches)

            covered_indices = [j for j, seq in enumerate(all_sequences) if j not in indices_covered and
                               regex.findall(left_regex, str(seq)) and regex.findall(right_regex, str(seq))]

            coverage = len(covered_indices) / len(all_sequences)
            if coverage > max_coverage[amplicon_range]:
                max_coverage[amplicon_range] = coverage
                best_pair[amplicon_range] = (pair, coverage)
                if coverage == max(max_coverage.values()):
                    best_covered_indices = covered_indices

    print(f"maximal number covered {len(best_covered_indices)} out of {len(nucleotides)} records")
    return best_covered_indices, best_pair
