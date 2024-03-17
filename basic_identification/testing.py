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
            database_path =  ospath.abspath(os.getcwd() + "/Database/" + gene_name + ".fasta")
            print(database_path)
            db, seq =  import_data(database_path, sequence_path, option = '1')

            # In case a family name is provided
        else: 
            family_db = ospath.abspath(os.path.dirname(os.getcwd()) + "/Database_by_family/" + family_name + "/" + gene_name + ".csv")
            family_db = r"C:/Users/ghass/projects/genorobotics/Identification/reorganization_by_family/Database_by_family/" + family_name + "/" +gene_name +".csv"
            print(family_db)
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