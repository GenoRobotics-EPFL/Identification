from utils import *
from algo import *
import pandas as pd
import time
import os.path as ospath
import os
import argparse

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

    # OPTION 1: comparison of a sequence with all the database
    if args.option == "1":

        # Managing errors
        if not ospath.exists(args.sequence_path):
            raise TypeError('The sequence path does not exist')

        else:
            start = time.time()
                # In case no family name was provided
            if args.family_name == None:
                database_path =  ospath.abspath(os.getcwd() + "/Database/" + args.gene_name)
                db, seq =  import_data(database_path, args.sequence_path, option = '1')

                # In case a family name is provided
            else: 
                family_db = ospath.abspath(os.getcwd() + "/Database_by_family/" + args.family_name + "/" + args.gene_name + ".csv")
                if ospath.exists(family_db):
                    db, seq =  import_data(family_db, args.sequence_path, option = '1')

            result = align(db,seq,"bioalign")
            end = time.time()
            print(result.head(10))
            print(f'Execution time : ', end-start, 's')
        

    elif args.option == "2" or args.option =="3":

        if not ospath.exists(args.model_path):
            raise TypeError('The template sequence path does not exist')
        if not ospath.exists(args.consensus_path):
            raise TypeError('The consensus sequence path does not exist')
        
        # OPTION 2: Comparison between 2 sequence using global_align
        elif args.option == "2":
            start = time.time()
            db, seq =  import_data(args.template_path, args.consensus_path, option = '2') 
            align(db, seq, "global_align")
            end = time.time()
            print(f'Execution time : ', end-start, 's')

        # OPTION 3: Comparison between 2 sequence using biopython
        elif args.option == "3":
            start = time.time()
            db, seq =  import_data(args.template_path, args.consensus_path, option = '2') 
            align(db, seq,"bioalign", option ="A")
            end = time.time()
            print(f'Execution time : ', end-start, 's')
