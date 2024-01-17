#!/usr/bin/env python3

##inputs : first argument : cluster_to_plot, second argument : ouptput_file.png
#small util to plot a consensus starting from multiple sequences aligned  in a clustal file (all sequences must have same length)
#you will see on X the indices of the bp, on Y in blue the percentage the consensus obtains (at 100% it means that all bp for that index are the same)
#added to that you will see in orange the percentage of 'holes' in the multiple sequence alignment (represented by the '-' symbol)

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import Counter

def plot_consensus(primer_index, cluster_filename_input, output_figure):
    # Read sequences from the input file (may vary the encoding)
    sequences = list(SeqIO.parse(cluster_filename_input, "fasta"))

    # Check if all sequences have the same length
    sequence_length = len(sequences[0])
    if any(len(seq) != sequence_length for seq in sequences):
        raise ValueError("All sequences must have the same length.")

    print(f"For primer {primer_index}, the amplicon size is {sequence_length}")
    # Convert sequences into a NumPy array
    alignements = np.array([list(str(seq.seq)) for seq in sequences]).T

    # Initialize variables for plotting
    x = range(sequence_length)
    count_a, count_t, count_c, count_g = [], [], [], []
    plt.figure(figsize=(20, 8))
    plt.title(f"Cluster {cluster_filename_input}")
    plt.xlabel("Index of base")
    plt.ylabel("Basis percentage")

    basis = ['A', 'T', 'C', 'G']
    # Iterate over each position in the sequences
    for i in range(sequence_length):
        counts = Counter(alignements[i])
        count = counts['A']
        count_a.append(100 * count / (len(alignements[i])))
        count = counts['T']
        count_t.append(100 * count / (len(alignements[i])))
        count = counts['C']
        count_c.append(100 * count / (len(alignements[i])))
        count = counts['G']
        count_g.append(100 * count / (len(alignements[i])))


    plt.bar(x, count_a, label="A")
    plt.bar(x, count_t, bottom=count_a, label="T")
    plt.bar(x, count_c, bottom=count_t, label="C")
    plt.bar(x, count_g, bottom=count_c, label="G")

    plt.legend()
    plt.xticks(x)

    # Save the plot to the specified output figure path
    plt.savefig(output_figure)

def treat_arguments():
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python script.py <cluster_filename_input.fasta> <output_figure.png>")
        sys.exit(1)

    # Get input and output paths from command line arguments
    cluster_filename_input = sys.argv[1]
    output_figure = sys.argv[2]

    # Check if the input file exists
    if not os.path.isfile(cluster_filename_input):
        print(f"Error: Input file '{cluster_filename_input}' not found.")
        sys.exit(1)

    # Check if the output path ends with '.png'
    if not output_figure.lower().endswith('.png'):
        print("Error: Output figure path must end with '.png'.")
        sys.exit(1)

    return cluster_filename_input, output_figure

def main():
    # Process command line arguments
    cluster_filename_input, output_figure = treat_arguments()

    # Generate and display the consensus plot
    plot_consensus(cluster_filename_input, output_figure)

if __name__ == "__main__":
    main()
