from Bio import SeqIO
import pandas as pd
import os.path as ospath


def parse_database(path):
    #Create a dictionnary from the database's file
    #Return a dict : 
    #Description of the sequence are the keys
    #DNA sequence are the values

    sequence = {}
    if path.endswith(".fasta"):
        for record in SeqIO.parse(path, "fasta"):
            sequence[record.description] = record.seq
    elif path.endswith(".csv"):
        db = pd.read_csv(path)
        for i in range(len(db.index)):
            species = db.loc[i]["genus"] + " " + db.loc[i]["species"]
            sequence[species] = db.loc[i]["sequence"]
    else:
        raise FileNotFoundError("Invalid file extension, must be .fasta or .csv")
    return sequence


def parse_sequence(path):
    #Put the sequence to align in a list

    sequence = []
    for record in SeqIO.parse(path, "fasta"):
        sequence = record.seq
    return sequence


def replace_nucleotide(seq):
    #Replace all the uncertain nucleotide by '-'

    table =  ["A", "T", "C", "G"]
    seq_clean = ""
    for i in range (len(seq)): 
        if seq[i] in table:
            seq_clean += seq[i]
        elif seq[i] not in table:
            seq_clean += '-'
    return seq_clean

def sequence_uppercase(seq):
    #Put all the nucleotide in uppercase

    seq_up = seq.upper()
    return seq_up


def import_data(database, sequence, option = '1'):
    #Compute the different function to preprocess all the sequence

    if option == '1':
        # Option 1 = Comparison of a DNA sequence with whole the database
        db =  parse_database(database)
        seq = parse_sequence(sequence)
        T = sequence_uppercase(seq)
        T = replace_nucleotide(T)
        
    if option == '2': 
        # Option 2 = Comparison of two DNA sequence
        db =  parse_sequence(database)
        seq = parse_sequence(sequence)
        T = replace_nucleotide(seq)
    return db, T
