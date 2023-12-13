import os
import subprocess
import pandas as pd
from file_utils import *
import regex

def generate_primer_for_cluster(cluster_index, output_folder, pairs_dictionnary):
    cluster_aln = f"{output_folder}/cluster{cluster_index}.aln"
    
    # Design primers
    for j in range(300, 1000, 80): ##amplicon between 300 and 1000 seems a reasonable range
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


def test_primer_coverage_against_dataset(nucleotides, primer_pairs):
    allowed_mismatches = 6
    for amplicon_range, pairs in primer_pairs.items():
        max_coverage = 0.0
        if (pairs):
            best_pair = pairs[0]
            for i in range(len(pairs)):
                primerR = pairs[i][1].seq
                primerL = pairs[i][0].seq
    
                
                left_regex, right_regex = generate_regex(primerL, primerR, allowed_mismatches)
                
                count = 0
                for nucleotide in nucleotides:
                    left_match = regex.findall(left_regex, str(nucleotide))
                    right_match = regex.findall(right_regex, str(nucleotide))
                
                    if (left_match and right_match):
                        count += 1
                coverage = count/len(nucleotides)
                if (max_coverage < coverage):
                    max_coverage = coverage
                    best_pair = pairs[i]
            print(f"pair {best_pair[0].seq} and {best_pair[1].seq} for range {amplicon_range} and coverage {round(max_coverage*100, 2)}%")