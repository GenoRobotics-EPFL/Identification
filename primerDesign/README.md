# Automatic Primer Designer

## Goal of the project
The goal of this project was to permit, given any dataset of genome sequences, to generate pairs of primers that have a coverage as high as possible over the whole dataset.

## File Structure - Algorithm
The algorithm we use has the following shape :
1. first we take some dataset of sequences (we want a dataset not larger than 10'000 sequences, because the algorithm is $O(n^{2})$. It performs better if the sequences are close to each other.
2. first we compute the mash distance : reference here : https://mash.readthedocs.io/en/latest/distances.html, using hashed kmers
3. Then we use those distances to cluster the sequences into groups, using the scipy fcluster library
4. Each cluster is then aligned using clustal omega
5. For each aligned cluster we run the primalscheme program, which returns some candidate primers, that we add up in a list
6. we test each pair of primers on the whole dataset, and retain the one that covers the most elements.
7. We repeat the algorithm on the remaining datasets while ensuring that no cross or self-binding occurs between the primers chosen.

## Usage - Installation
After having cloned the project, you will need to install python3, pip and clustal omega (clustalo) via your favorite package installer.

Then install the required libraries using this command 

```console
sudo pip install bio biopython tqdm mmh3 scipy pandas regex primer3-py
```

To use primalscheme and be able to configure it properly, you will to clone [this source github repository](https://github.com/aresti/primalscheme), navigate to its root folder.

Under the path *primalschem/src/primalscheme/config.py* you will find the configuration file, where you can tune the different parameters for RPA. Here are the main differences on my side :

- primer size : 30-36, optimal 33
- %gc : 30-70, optimal 55
- maximal alignement ap percent : 0.1
- primer maximum mismatches : 6

Then you will be ready to install primalscheme by executing the following commands :

```console
sudo pip install flit
flit install --pth-file
```

## Usage - Command line arguments
You can use this command to run the program.
```console
  python src/primerDesigner.py -f filename-fasta
```
You must provide filename-fasta, th path to your dataset of sequences in __fasta__ format.

And you could provide some additional arguments if needed :
* ```-o output_folder``` Set the output folder (default is './output'). Please note that __all data from previous runs in output folder will be deleted__
* ```-k kmer_length``` Set the length of the kmers used in mash distance (study the mash distance algorithm before tuning that argument). Default value is 15.
* ```-d distance-within_cluster``` Set the distance allowed within cluster for the _fcluster_ function. This value must be between 0 and 1, and the smaller it is the smaller the cluster will be, reducing computation time without harming primer quality that much. Default value is 0.1.
* ```-n number-cluster``` Set the number of groups retrieved at point 3 of the algorithm. The default value is 3.
* ```-p numer_pairs``` Set the number of primer pairs finally returned by the algorithm. The default value is 3.
* ```-min_range min_range``` and ```-max_range max_range``` Set the minimal and maximal amplicon size. Please note that this is an iterative process the larger that range will be, the longer the program will run. Default values are (min_range=300, max_range=1000)

## Citing

I used the following external programs : 
- [Clustal Omega](https://europepmc.org/article/MED/35412617) for multiple sequence alignements
- [primalscheme](https://www.nature.com/articles/nprot.2017.066) for primer generation

