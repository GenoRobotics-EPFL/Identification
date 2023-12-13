from Bio import SeqIO
import multiprocessing
from itertools import chain
from functools import partial
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

def extract_kmers_from_sequence(sequence, k_length):
    sequence = [*(''.join(sequence)).strip('-')]
    kmers = {tuple(sequence[i:i + k_length]) for i in range(min(len(sequence) - k_length + 1, 200))}
    kmers |= {tuple(sequence[i:i + k_length]) for i in range(max(0, len(sequence) - 200), len(sequence) - k_length + 1)}
    return kmers

def retrieve_kmers(records, k_length):
    nucleotides = list(map(lambda x : extract_kmers_from_sequence(x, k_length), [record.seq for record in records]))
    return nucleotides


def mash_distance(x, y):
    intersection = len(x.intersection(y))
    union = len(x.union(y))
    jaccard_similarity = intersection / union
    return jaccard_similarity

def calculate_mash_for_index(args):
    index, nucleotides, nb_elements = args
    array_output = []
    for j in range(index + 1, nb_elements):
        similarity = 1.0 - mash_distance(nucleotides[index], nucleotides[j])
        array_output.append(similarity)
    return array_output


def compute_mash_distance_matrix(nucleotides):
    num_elements = len(nucleotides)
    args_list = [(i, nucleotides, num_elements) for i in range(num_elements)]

    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        dist_matrix = list(tqdm(executor.map(calculate_mash_for_index, args_list, chunksize=10), total=num_elements))

    return list(chain.from_iterable(dist_matrix))


