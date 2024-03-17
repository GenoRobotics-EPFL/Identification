import sys
sys.path.append("..")
from utils import *
from algo import *
import time
import os.path as ospath
import os
import argparse



def run_identification(sequence_path, gene_name, family_name):
    """
    Function that runs the identification search (algorithm option "1")
    Parameters
    ----------
        sequence_path: path to the sequence to identify (".fasta")
        gene_name: name of the sequenced gene ("matk", "rbcl", "trnh-psba", "its")
        family_name: name of the family taxon 
            can be "None"
    
    Returns
    -------
    10 best alignments with scores and corresponding taxonomy 
    """
    # OPTION 1: comparison of a sequence with all the database

    # Managing errors
    if not ospath.exists(sequence_path):
        raise TypeError('The sequence path does not exist')

    else:
        start = time.time()
            # In case no family name was provided
        if family_name == None:
            print(os.getcwd())
            database_path =  ospath.dirname(os.getcwd()) + "\\Database\\" + gene_name + ".fasta"
            db, seq =  import_data(database_path, sequence_path, option = '1')

            # In case a family name is provided
        else: 
            family_db = ospath.dirname(os.getcwd()) + "\\reorganization_by_family\\Database_by_family\\" + family_name + "\\" + gene_name + ".csv"
            if ospath.exists(family_db):
                db, seq =  import_data(family_db, sequence_path, option = '1')

        ranking = align(db, seq, "bioalign")
        end = time.time()
        ranking = ranking.head(10)

        print(ranking)
        # create a function that returns a dictionary with all taxonomy of our species result
        # ranking_taxonomy = get_taxonomy(ranking["Species"])
        id_result = [(ranking["Species"], ranking["Alignment score"])]
        print(f'Execution time : ', end-start, 's')

        return id_result
        


def run_alignment(option, template_path, consensus_path):
    """
    Function that runs alignment (algorithm option "2" or "3")
    Parameters
    ---------
        option: algorithm choice 


    """

    if not ospath.exists(template_path):
        raise TypeError('The template sequence path does not exist')
    if not ospath.exists(consensus_path):
        raise TypeError('The consensus sequence path does not exist')
    
    # OPTION 2: Comparison between 2 sequence using global_align
    elif option == "2":
        start = time.time()
        db, seq =  import_data(template_path, consensus_path, option = '2') 
        align(db, seq, "global_align")
        end = time.time()
        print(f'Execution time : ', end-start, 's')

    # OPTION 3: Comparison between 2 sequence using biopython
    elif option == "3":
        start = time.time()
        db, seq =  import_data(template_path, consensus_path, option = '2') 
        align(db, seq,"bioalign", option ="A")
        end = time.time()
        print(f'Execution time : ', end-start, 's')



if __name__ == "__main__":
    #Run all the programm
    #Run the options (see the Readme file for more information on it)
    #function time : allow to see the runtime of each different algorithm 
    #print the result: 
        #   DataFrame with the species and the score for option 1 
        #   Alignment and score for option 2 and 3

    parser = argparse.ArgumentParser(epilog="See '<option> --help' to read about a specific sub-command.")
    subparser = parser.add_subparsers(dest = "option", help = "Sub-commands")

    search = subparser.add_parser("1", help = "Sequence Identification")
    search.add_argument("sequence_path", help = "Path to the consensus sequence")
    search.add_argument("gene_name", help = "Name of sequenced gene", choices = ["matk", "rbcl", "psbA-trnh", "its"])
    search.add_argument("--family_name", help = "If known, family name of sequenced species")

    compare2 = subparser.add_parser("2", help = "Sequence comparison with NW algo")
    compare2.add_argument("template_path", help = "Path to the template sequence")
    compare2.add_argument("consensus_path", help = "Path to the consensus sequence")

    compare3 = subparser.add_parser("3", help = "Sequence comparison with Bioalign")
    compare3.add_argument("template_path", help = "Path to the template sequence")
    compare3.add_argument("consensus_path", help = "Path to the consensus sequence")
    
    args = parser.parse_args()

    if args.option == "1":
        run_identification(args.sequence_path, args.gene_name, args.family_name)
    elif args.option == "2" or "3":
        run_alignment(args.option, args.template_path, args.consensus_path)
