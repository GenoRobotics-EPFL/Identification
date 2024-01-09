import subprocess
import time
import pandas as pd

# Specify the number of times you want to run the command
num_runs = 5

gene_to_test = "ITS2"

gene_path = f"examples/{gene_to_test}_test/{gene_to_test}_test.fasta"

number_clusters = [1, 3, 5, 7]
kmer_lengths = [4, 8, 15, 30, 60]
distance_within_clusters = [0.001, 0.01, 0.1, 0.3, 0.5]
mash_lengths = [50, 100, 200, 400]

def read_coverage(csv_path):
    df = pd.read_csv(csv_path)
    coverage_sum = df['Coverage [%]'].sum()
    print(f"Coverage : {round(coverage_sum, 2)}%\n")

def test_clusters_nb():
    print("testing number clusters")
    for nb_cluster in number_clusters:
        print(f"Testing with {nb_cluster} clusters")

        # Use subprocess to run the command
        start_time = time.time()
        command = f"time python src/primerDesigner.py -f {gene_path} -o stats/output/{gene_to_test}/cluster_{nb_cluster} -n {nb_cluster}"
        subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        end_time = time.time()

        # Calculate and print the elapsed time for each run
        elapsed_time = end_time - start_time
        read_coverage(f"stats/output/{gene_to_test}/cluster_{nb_cluster}/output.csv")
        print(f"Elapsed Time to test {nb_cluster} clusters: {elapsed_time:.2f} seconds\n")

def test_kmers_length():
    print("testing kmer_length")
    for kmer_length in kmer_lengths:
        print(f"Testing with {kmer_length} length")

        # Use subprocess to run the command
        start_time = time.time()
        command = f"time python src/primerDesigner.py -f {gene_path} -o stats/output/{gene_to_test}/kmer_length_{kmer_length} -k {kmer_length}"
        subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        end_time = time.time()

        # Calculate and print the elapsed time for each run
        elapsed_time = end_time - start_time
        read_coverage(f"stats/output/{gene_to_test}/kmer_length_{kmer_length}/output.csv")
        print(f"Elapsed Time to test kmer_length {kmer_length}: {elapsed_time:.2f} seconds\n")

def test_cluster_dist():
    print("testing cluster_dist")
    for cluster_dist in distance_within_clusters:
        print(f"Testing with {cluster_dist} cluster distance")

        # Use subprocess to run the command
        start_time = time.time()
        command = f"time python src/primerDesigner.py -f {gene_path} -o stats/output/{gene_to_test}/cluster_dist_{cluster_dist} -d {cluster_dist}"
        subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        end_time = time.time()

        # Calculate and print the elapsed time for each run
        elapsed_time = end_time - start_time
        read_coverage(f"stats/output/{gene_to_test}/cluster_dist_{cluster_dist}/output.csv")
        print(f"Elapsed Time to test cluster_dist {cluster_dist}: {elapsed_time:.2f} seconds\n")


def test_mash_length():
    print("testing mash_length")
    for mash_length in mash_lengths:
        print(f"Testing with {mash_length} mash length")

        # Use subprocess to run the command
        start_time = time.time()
        command = f"time python src/primerDesigner.py -f {gene_path} -o stats/output/{gene_to_test}/mash_length_{mash_length} -m {mash_length}"
        subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        end_time = time.time()

        # Calculate and print the elapsed time for each run
        elapsed_time = end_time - start_time
        read_coverage(f"stats/output/{gene_to_test}/mash_length_{mash_length}/output.csv")
        print(f"Elapsed Time to test mash_length {mash_length}: {elapsed_time:.2f} seconds\n")


test_clusters_nb()
print("\n\n")
test_kmers_length()
print("\n\n")
test_cluster_dist()
print("\n\n")
test_mash_length()

print("All runs completed.")
