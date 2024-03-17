# Basic alignment and identification

In this early version of the project, BioPython's global alingmnent function (pairwise2) as well as a handmade Needleman-Wuntsch algorithm are used to find the best matches for each gene. This was later replaced by BLASTn, which dramatically increased the speed of alignment.

## Algorithm Needleman-Wunsch
The Needleman-Wunsch algorithm is an algorithm that performs a maximum global alignment of two DNA sequences. It is commonly used in bioinformatics to align protein or nucleotide sequences. The Needleman-Wunsch algorithm is an example of dynamic programming, it guarantees to find the maximum score alignment. To determine the maximum score alignment, a two-dimensional array, or matrix, is used. There is one row for each character in sequence A, and one column for each character in sequence B. So, if we align sequences of size n and m, the execution time of the algorithm is O(nm), and the memory space used is O(nm) too.

The algorithm is executed in three steps: 

1. Calculation of the similarity matrix, matrix S. For each cell of the matrix we assign a score Z if the nucleotide of sequence A is equal to that of sequence B and a score Y if there is a substitution

2. Then we calculate the optimal matrix, matrix M. For each cell of M(i,j), we take the maximum between 3 cases: 
M(i,j) = max (M(i-1, j-1)+S(i,j), M(i-1, j)+g, M(i, j-1) + g)
with g being the weight assigned to the gap.

<img width = "200" src = ../Images/NW1.png>

3. We look at the matrix M from the last cell to the first. For each cell we go to the best score between the 3 before. So for the cell M(i, j), we look at M(i, j-1), M(i-1, j) and M(i-1, j-1) and we keep only the best score.

<img width = "300" src = ../Images/NW2.png>

Finally this gives us an optimal alignment of the two sequences with gaps and shifts.


## Implementation
Library: Only three library are used: Numpy, Pandas and Biopython. From Biopython we only use SeqIO to parse the fasta files and pairwise2 to align two DNA sequences. 

utils.py : All the function use to preprocess the file : 
- parse_sequence : Parse the database in a dictionnary, Sequence information are the keys and DNA sequence are the values.
- parse_database : Parse the sequence in a list. 
- replace_nucleotide : Verify if the sequence provide is completely correct, if not change the nucleotide by a '-'. This choice is a personal choice to make the computation easier and reduce the runtime. Could be change to improve the accuracy. 
- nucleotide_uppercase : Put all the sequence in uppercase

algo.py : All the function about the algorithm
- global_align : Needleman-Wunsch algorithm, compare two sequence to find the alignment and the matching score. Could be a bit long. Work with score_align() to compute the alignment score
- bioalign : Needleman-Wunsch algorithm coding in the library biopython, compare the sequence with all the database to find the sequence with the best alignment and the best matching score. It's the fastest method. 
- align : Compute the alignment and the score using the different function above

Identification_GR.py : main function to run the programm.



## Utilisation
To use this project you have to clone the main repository and use the following command line in the terminal.

**\python Identification_GR.py [OPTION]**

OPTION = 3 possible choices of algorithm : 
- [1] : Sequence identification, identifies a sequence within the database, the 10 best matching scores will return.
- [2] : Sequence alignment, returns the best alignment and the matching score (in percent). 
Uses global_algo, a Needleman-Wunsch Algorithm coded internally.
- [3] : Sequence alignment, returns the best alignment and the matching score (in percent). 
Uses a Needleman-Wunsch Algorithm with the librart Biopython.

**\python Identification_GR.py [option = 1] [sequence_path] [gene_name] [--family_name]**
- SEQUENCE_PATH = The path of your sequence (must be a fasta file, with extension ".fasta")
- GENE_NAME = Name of the gene, choices: "matk", "rbcl", "psbA-trnh", "its"
- --FAMILY_NAME = Optional, name of the Family taxon (if known) 

**\python Identification_GR.py [option = 2 or 3] [template_path] [consensus_path]**
- TEMPLATE_PATH = The path of your template sequence (must be a fasta file)
- CONSESUS_PATH = The path of your consensus sequence (must be a fasta file)



## Result
This is an example result using a test file.

Result for matK with option 1: 

<img width="665" src="https://github.com/GenoRobotics-EPFL/Identification/assets/102163457/6c3aff39-1c7e-4c44-945a-56038e5d0aea">