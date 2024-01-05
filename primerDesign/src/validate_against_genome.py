from file_utils import generate_regex
from Bio import SeqIO, Seq
import csv
import regex

def validate_against_chloroplast(primers, ref_chloroplast_path):
    ref_seq = list(SeqIO.parse(ref_chloroplast_path, "fasta"))[0].seq
    for forward, reverse in primers:
        regex_forward, regex_reverse = generate_regex(forward, reverse)
        if (regex.findall(regex_forward, str(ref_seq)) and
            regex.findall(regex_reverse, str(ref_seq))):
            return True
    return False
    
def extract_primer_from_csv(csv_path):
    primer_tuples = []
    with open(csv_path, 'r', newline='') as csv_file:
        csv_reader = csv.reader(csv_file)
        next(csv_reader)

        for row in csv_reader:
            forward_primer_seq = row[1]
            reverse_primer_seq = row[2]

            primer_tuples.append((forward_primer_seq, reverse_primer_seq))

    return primer_tuples

def validate_against_genome(ITS_primer_path, matK_primer_path, 
                            rbcL_primer_path, psbA_trnH_primer_path,
                            ref_chloroplast_path, ref_ITS_path):
    
    if (validate_against_chloroplast(extract_primer_from_csv(psbA_trnH_primer_path), ref_chloroplast_path)):
        print("psbA-trnH validated !")
    else :
        print("psbA-trnH not validated..")

    if (validate_against_chloroplast(extract_primer_from_csv(matK_primer_path), ref_chloroplast_path)):
        print("matK validated !")
    else :
        print("matK not validated..")

    if (validate_against_chloroplast(extract_primer_from_csv(rbcL_primer_path), ref_chloroplast_path)):
        print("rbcL validated !")
    else :
        print("rbcL not validated..")
    
    if (validate_against_chloroplast(extract_primer_from_csv(ITS_primer_path), ref_chloroplast_path)):
        print("ITS validated !")
    else :
        print("ITS not validated..")
    return 0

validate_against_genome("output/ITS.csv", "output/matK.csv", "output/rbcL.csv", "output/psbA-trnH.csv", "examples/solanum_lycopersicum/chloroplast.fasta", "")