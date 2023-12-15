import os
import subprocess
import pandas as pd
from file_utils import *
import regex
import primer3
from collections import defaultdict

MINIMAL_MELTING_TEMPERATURE = 50

def assert_no_self_binding(primer_sequence):

    # Run primer3 and check hairpin formation
    result = primer3.calcHairpin(primer_sequence)
    return result.tm < MINIMAL_MELTING_TEMPERATURE

def assert_no_cross_binding(primer_forward, primer_reverse):

    # Run primer3 and check hairpin formation
    result = primer3.calcHeterodimer(primer_forward, primer_reverse)
    return result.tm < MINIMAL_MELTING_TEMPERATURE

def assert_pair_primer_bindings(primer_forward, primer_reverse):

    reverse_complement = str(Seq(primer_reverse).reverse_complement())

    return assert_no_self_binding(primer_forward) and assert_no_self_binding(reverse_complement) and assert_no_cross_binding(primer_forward, reverse_complement)

def check_bindings(primer_forward, primer_reverse, pairs_found):
    if (pairs_found):
        for tuple in pairs_found.values():
            for forward, reverse in tuple:
                forward = forward.seq
                reverse = reverse.seq
                if not (assert_pair_primer_bindings(primer_forward, reverse) and assert_pair_primer_bindings(forward, primer_reverse)
                    and assert_no_cross_binding(primer_forward, forward) and assert_no_cross_binding(reverse, primer_reverse)):
                    return False
    return assert_pair_primer_bindings(primer_forward, primer_reverse)

def generate_primer_for_cluster(cluster_index, output_folder, pairs_dictionnary, min_range=300, max_range=500):
    cluster_aln = f"{output_folder}/cluster{cluster_index}.aln"
    
    # Design primers
    for j in range(min_range, max_range, 80): ##amplicon between 300 and 1000 seems a reasonable range
        command = [
            "primalscheme", "multiplex",
            "-a", str(j),
            "-a", str(j + 80),
            "-f", f"{cluster_aln}",
            "-o", f"{output_folder}/primalscheme_results{cluster_index}"
        ]
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # Meant to be used on a Linux terminal
        
        primer_results_path = f"{output_folder}/primalscheme_results{cluster_index}/primer.primer.tsv"
        if os.path.exists(primer_results_path):
            data = pd.read_csv(primer_results_path, sep='\t')
            
            # Create pairs of rows
            pair_key = "{}-{}".format(j, j + 80)
            pairs_dictionnary[pair_key].extend(
                (data.iloc[l], data.iloc[l + 1]) for l in range(0, len(data), 2)
            )


def test_primer_coverage_against_dataset(alignements, primer_pairs, indices_covered, pairs_found):

    nucleotides = [alignements[i] for i in range(len(alignements)) if i not in indices_covered] ##filter recods for already found primers
    primer_pairs = {k : [pair for pair in v if check_bindings(pair[0].seq, pair[1].seq, pairs_found)] for k, v in primer_pairs.items()} ##eliminate pairs who fail self/cross binding

    allowed_mismatches = 6
    best_covered_indices = []
    max_coverage = {}
    best_pair = {"": 0}
    if not primer_pairs:
        print("\nno additional primers found \n")
        return best_covered_indices, {}
    for amplicon_range, pairs in primer_pairs.items():
        max_coverage[amplicon_range] = 0.0
        if (pairs):
            best_pair_for_amplicon_size = {amplicon_range : pairs[0]}
            for i in range(len(pairs)):
                covered_indices = []
                primerR = pairs[i][1].seq
                primerL = pairs[i][0].seq
    
                
                left_regex, right_regex = generate_regex(primerL, primerR, allowed_mismatches)
                
                count = 0
                for j in range(len(alignements)):
                    if (j not in indices_covered):
                        left_match = regex.findall(left_regex, str(alignements[j]))
                        right_match = regex.findall(right_regex, str(alignements[j]))
                    
                        if (left_match and right_match):
                            count += 1
                            covered_indices.append(j)
                
                coverage = count/len(nucleotides)
                if (max_coverage[amplicon_range] < coverage):
                    max_coverage[amplicon_range] = coverage
                    if (max(max_coverage.values()) == coverage):
                        best_covered_indices = covered_indices
                        best_pair = {amplicon_range: pairs[i]}
                    best_pair_for_amplicon_size = {amplicon_range : pairs[i]}
            print(f"pair {list(best_pair_for_amplicon_size.values())[0][0].seq} and {list(best_pair_for_amplicon_size.values())[0][1].seq} for range {amplicon_range} and coverage {round(max_coverage[amplicon_range]*100, 2)}%")
    print(f"maximal number covered {len(best_covered_indices)} out of {len(nucleotides)} records")
    return best_covered_indices, best_pair