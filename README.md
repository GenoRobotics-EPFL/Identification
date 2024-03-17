# DNA identification 

### GenoRobotics

GenoRobotics is an interdisciplinary project associating engineers, scientists and students towards the development of a portable device enabling field DNA analysis. The main goal of the association is to develop a miniaturized tool to automatically process samples into the field and accelerate biodiversity identification

### Bioinformatics

The bioinformatics team of the association was split in three different project : 
- Database project : Create our own database to store all the information collected and not depend on an external agent.
- Consensus sequence : Find the consensus sequence of the sequence extracted from nanopore sequencing of a plant.
- Identification project : Compare this consensus sequence with a database to find wich plant it is. 

### Identification Project

The goal of this project was to find a way to compare a DNA sequence with a large database as fast as possible. 
The first step was about the database. We had to find a way to download the ncbi database to use it locally on our computer. Next, we have decided to create a second database organized into the family taxon to improve runtime. 
The second step was to find an algorithm efficient to compare the DNA sequence with several tens of thousands sequence. We have to compare the 4 gene extracted from the plant (MatK, rbcL, psbA-trnH, ITS) with nanopore sequencing with this database to find which plant it is,  or,  if this plant is not in the database, the closest plant to it. To do that we use the Needleman-Wunsch algorithm. 


### Download the useful database

To improve the efficiency of the algorithm you have to download only the part of the database of your interest. With the nanopore sequencing using MinION technology, we extract 4 genes (MatK, rbcL, psbA-trnH, ITS) so it's important to download the database of these genes. 
To download the best database possible please follow this few rules :
- Go on the website: https://www.ncbi.nlm.nih.gov/nuccore
- Enter the gene you need (MatK, rbcL, psbA-trnH, ITS)
- Select the length of the sequence of interest with [Sequence length] on the left (allow to select less sequence and only the sequence of interest)
Example of length : 
    - MatK : 750 to 1500 -> ~110k Sequence, ~90Mo
    - rbcL : 600 to 1000 -> ~90k Sequence, ~80Mo
    - psbA-trnH : 400 to 800 -> ~60k Sequence, ~40Mo
    - its: 1000 to 35000
- Download it to have the information and the DNA sequence 
Click and send to (corner top right) > Complete Record > File > Format = Fasta > Sort by : (as you want it's not important)
<img width = "660" src = Images/ncbi.png>
<img width = "250" src = Images/ncbi2.png>

## Sub-folders

### Basic Identification with BioPython

In this early version of the Identification code, Biopython's `pairwise2` algorithm is used to align the query sequence with all the sequences in the database. A handmade version of the Needleman-Wuntsch algorithm, while much slower, provides a basis for understanding how this algorithm works. A command lined execution is provided and explained in [the subfolder's README](basic_identification/README.md)

### reorganization_by_family

In this mini-project, we try to reorganize the databases extracted from GenBank according to the species' family. The rationale is that a lot of time can be saved if the query sequence is only aligned with reference sequences of plants from the same family. This can be used if the family of a sample can be determined from morphological indicators but the exact species remains ambiguous. Each folder represents a family and the files it contains are in csv format. This formal is not compatible with the use of a BLASTn as an alignment tool. Moreover, the efficiency and speed of BLASTn also puts into question the usefulness of using this sorted version of the database. More details on the sorting of the database in [the subfolder's README](reorganization_by_family/README.md)

### GenBank Analysis

In this unfinished mini-project, we try to establish a routine to analyse a database of sequences extracted from GenBank. It determines average sequence length, nucleotide percentage and similarity between sequences (unfinished part). This could be useful to determine the optimal range of sequences to extract a database that has all necessary information without being too large.





