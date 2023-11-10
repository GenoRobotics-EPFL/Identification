# BLASTn How To :
## Part 1: Installing BLASTn
### Download BLAST:

Visit the NCBI BLAST download page: [BLAST Download](https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/).
Choose the appropriate version for your operating system (Windows, macOS, or Linux), and follow the instructions provided for your specific OS

### Configure BLAST:

(Windows) To be able to call BLASTn directly from your terminal, you need to add several user environment variables.

- A modified path environment variable to indicate the location of installed blast+ programs, with "C:\users\Milou\Desktop\blast-2.10.0+\bin\;" or whatever the absolute path to your blast directory is, prepended to its existing value
- A new BLASTDB environment variable as pointer to database location, with "C:\users\Milou\Desktop\blast-2.10.0+\db\", or whatever the absolute path to your blast db directory is, as its value
- A new BLASTDB_LMDB_MAP_SIZE, with 1000000 as its value (needed to optimize makeblastdb operation when creating new database files)

## Part 2: Setting Up a Local Database
### Download Gene Sequences: (by @AlexMlld, the same as for the previous algorithm)

To improve the efficiency of the algorithm you have to download only the part of the database of your interest. With the nanopore sequencing using MinION technology, we extract 4 genes (MatK, rbcL, psbA-trnH, ITS) so it's important to download the database of these genes. To download the best database possible please follow this few rules :

- Go on the [NCBI website](https://www.ncbi.nlm.nih.gov/nuccore)
- Enter the gene you need (MatK, rbcL, psbA-trnH, ITS)
- Select the length of the sequence of interest with [Sequence length] on the left (allow to select less sequence and only the sequence of interest) Example of length :
- MatK : 750 to 1500 -> ~110k Sequence, ~90Mo
- rbcL : 600 to 1000 -> ~90k Sequence, ~80Mo
- psbA-trnH : 400 to 800 -> ~60k Sequence, ~40Mo
- Download it to have the information and the DNA sequence Click and send to (corner top right) > Complete Record > File > Format = Fasta > Sort by : (as you want it's not important)
- Put all those fasta files in the db directory you previously added to your environment variables.

### Create a Local Database:

Use the makeblastdb command-line tool included with BLAST : 

makeblastdb -in `<db name.fasta>` -dbtype nucl -parse_seqids -out `<output name>`

You can check it was correctly installed by asking infos about the resulting db : 

 blastdbcmd -db `<db name>` -info

## Part 3: Using BLASTn with Your Local Database
### Prepare the Query Sequence:

Go in the root folder of blast, and put here a fasta file you want to query in your db. We included several fasta files so you can test by yourself. 

You can also extract a fasta from your database with the following command : 

blastdbcmd -db `<db name>` -entry `<gene name>` > `<output file name.fasta>`

### Run BLASTn:

Use the blastn command-line tool : 

blastn -query `<input file.fasta>` -db `<db name>` -out `<output name.txt>` -max_target_seqs 5

Check your .txt result to interpret results ! 

# Full tutorial and more resources

https://www.ncbi.nlm.nih.gov/books/NBK52637/ 
