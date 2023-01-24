# DNA identification 

## Goal of the project
### GenoRobotics
GenoRobotics is an interdisciplinary project associating engineers, scientists and students towards the development of a portable device enabling field DNA analysis. The main goal of the association is to develop a miniaturized tool to automatically process samples into the field and accelerate biodiversity identification

### Bioinformatic
The bioinformaitc team of the association was split in three different project : 
- Database project : Create our own database to store all the information collected and not depend on an external agent.
- Consensus sequence : Find the consensus sequence of the sequence extracted from the plant.
- Identification project : Compare this consensus sequence with the database to find wich plante it is. 

### Identification Project
The goal of this project was to find a way to compare a DNA sequence with a large database as fastest as possible. 
The first step was about the database. We had to find a way to download the ncbi database to use it locally on our computer.
The second step was to find an algorithm efficient to compare the DNA sequence with several tens of thousands sequence. We have to compare the 4 gene extracted from the plant (1, 2, 3, 4) with microneedle with this database to find which plant it is,  or,  if this plant is not in the listed, the closest plant to it. To do that we use the Needleman-Wunsch algorithm. 


## Download the useful database
To improve the efficiency of the algorithm you have to download only the part of the database of your interest. With the microneedle technologie, we extract 4 genes (1, 2, 3, 4) so it's important to download the database of these genes. 
To download the best database possible please follow this few rules :
- Go on the website: https://www.ncbi.nlm.nih.gov/nuccore
- Enter the gene yu need (MatK, )
- Select the length of the sequence of interest with [Sequence length] on the left (allow to select less sequence and only the sequence of interest)
- Download it to have the information and the DNA sequence 
Click and send to (corner top right) > Complete Record > File > Format = Fasta > Sort by : (as you want it's not important)

#### ADD ncbi1 and NCBI2

## Algorithme Needleman-Wunsch
The Needleman-Wunsch algorithm is an algorithm that performs a maximum global alignment of two DNA sequences. It is commonly used in bioinformatics to align protein or nucleotide sequences. The Needleman-Wunsch algorithm is an example of dynamic programming, it guarantees to find the maximum score alignment. To determine the maximum score alignment, a two-dimensional array, or matrix, is used. There is one row for each character in sequence A, and one column for each character in sequence B. So, if we align sequences of size n and m, the execution time of the algorithm is O(nm), and the memory space used is O(nm) too.

The algorithm is executed in three steps: 

1. Calculation of the similarity matrix, matrix S. For each cell of the matrix we assign a score \alpha if the nucleotide of sequence A is equal to that of sequence B and a score \beta if there is a substitution

2. Then we calculate the optimal matrix, matrix O. For each cell of O(i,j), we take the maximum between 3 cases: 
O(i,j) = max (O(i-1, j-1)+S(i,j), O(i-1, j)+g, O(i, j-1) + g)
with g being the weight assigned to the gap. 
##### ADD NW1

3. We look at the matrix O from the last cell to the first. For each cell we go to the best score between the 3 before. So for the cell O(i, j), we look at O(i, j-1), O(i-1, j) and O(i-1, j-1) and we keep only the best score. 
#### ADD NW2

Finally this gives us an optimal alignment of the two sequences with gaps and shifts.

Translated with www.DeepL.com/Translator (free version)

## Implementation
Library: Only two library are used: Numpy and Biopython. From Biopython we only use SeqIO to parse the fasta files and pairwise2 to align two DNA sequence. 

utils.py : All the function use to preprocess the file : 
- parse_sequence : Parse the database in a dictionnary, Sequence information are the keys and DNA sequence are the values.
- parse_database : Parse the sequence in a list. 
- replace_nucleotide : Verify if the sequence provide is completely correct, if not change the nucleotide by a '-'. 

algo.py : All the function about the algorithm
- global_align : Needleman-Wunsch algorithm, compare two sequence to find the alignment and the matching score. Could be a bit long. 
- bioalign : Needleman-Wunsch algorithm coding in the library biopython, compare the sequence with all the database to find the sequence with the best alignment and the best matching score. It's the fastest method. 

Blastn_GR.py : main function to run the programm.


## Utilisation
To use this project you have to clone the main repository and use the following command line in the terminal.
\python Blastn_GR.py [DATABASE_PATH] [SEQUENCE_PATH] [ALGORITHM]
###### Rajouter des print screen des résultats pour les 3 options
DATABASE_PATH = The path of your database (must be a fasta file)
SEQUENCE_PATH = The path of your sequence (must be a fasta file)
ALGORTHM = 3 possible choice of algorithm : 
- [1] : To compare your sequence with all the database, only the best matching score will return
- [2] : To align your sequence with another sequence, return the best alignment and the matching score (in percent). 
Use global_algo, a Needleman-Wunsch Algorithm coding by myself 
- [3] : To align your sequence with another sequence, return the best alignment in a different way and the matching score. 
Use a Needleman-Wunsch Algorithm using the librart Biopython.

## Possible Amelioration


...