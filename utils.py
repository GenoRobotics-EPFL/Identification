from Bio import SeqIO


def parse_database(path):
    #Create a dictionnary from the database's file
    #Description of the sequence are the keys
    #DNA sequence are the values
    sequence = {}
    for record in SeqIO.parse(path, "fasta"):
        sequence[record.description] = record.seq
    return sequence


def parse_sequence(path):
    #Put the sequence to align in a list
    sequence = []
    for record in SeqIO.parse(path, "fasta"):
        sequence = record.seq
    return sequence


def replace_nucleotide(seq):
    #Replace all the false nucleotide by '-'
    table =  ["A", "T", "C", "G"]
    seq_clean = ""
    for i in range (len(seq)): 
        if seq[i] in table:
            seq_clean += seq[i]
        elif seq[i] not in table:
            seq_clean += '-'
    return seq_clean


def import_data(database, sequence, option = '1'):
    #Compute the different function to preprocess all the sequence
    if option == '1':
        # Option 1 = Comparison of a DNA sequence with whole the database
        db =  parse_database(database)
        seq = parse_sequence(sequence)
        T = replace_nucleotide(seq)
    if option == '2': 
        # Option 2 = Comparison of two DNA sequence
        db =  parse_sequence(database)
        seq = parse_sequence(sequence)
        T = replace_nucleotide(seq)
    return db, T
