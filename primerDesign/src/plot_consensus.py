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

def plot_consensus(cluster_filename_input, output_figure):
    # Read sequences from the input file (may vary the encoding)
    sequences = list(SeqIO.parse(cluster_filename_input, "clustal"))

    # Check if all sequences have the same length
    sequence_length = len(sequences[0])
    if any(len(seq) != sequence_length for seq in sequences):
        raise ValueError("All sequences must have the same length.")

    # Convert sequences into a NumPy array
    alignements = np.array([list(str(seq.seq)) for seq in sequences]).T

    # Initialize variables for plotting
    x = []
    consensus_y = []
    holes_y = []

    # Iterate over each position in the sequences
    for i in range(sequence_length):
        counts = Counter(alignements[i])
        most_common, count = counts.most_common(2)[0]

        # Skip positions where all sequences have a gap ("-")
        if (most_common == "-" and count == len(alignements[i])):
            count = 0
            continue
        elif (most_common == "-"):
            most_common, count = counts.most_common(2)[1]

        # Populate the plotting variables
        x.append(i)
        consensus_y.append(100 * count / (len(alignements[i])))
        holes_y.append(100 * counts["-"] / (len(alignements[i])))

    # Create the consensus plot
    plt.figure(figsize=(20, 8))
    plt.title("Cluster {}".format(cluster_filename_input))
    plt.xlabel("Index of base")
    plt.ylabel("Consensus percentage")
    plt.bar(x, consensus_y, label="Consensus")
    plt.bar(x, holes_y, bottom=consensus_y, label="Unknown")
    plt.legend()
    plt.xticks(x)

    # Save the plot to the specified output figure path
    plt.savefig(output_figure)
    plt.show()

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
        print(f"Error: Output figure path must end with '.png'.")
        sys.exit(1)

    return cluster_filename_input, output_figure

def main():
    # Process command line arguments
    cluster_filename_input, output_figure = treat_arguments()

    # Generate and display the consensus plot
    plot_consensus(cluster_filename_input, output_figure)

if __name__ == "__main__":
    main()
