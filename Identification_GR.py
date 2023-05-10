from utils import *
from algo import *
import pandas as pd
import time
import sys
import os.path as ospath
import os

if __name__ == "__main__":
    #Run all the programm
    #Run the options (see the Readme file for more information on it)
    #function time : allow to see the runtime of each different algorithm 
    #print the result: 
        #   DataFrame with the species and the score for option 1 
        #   Alignment and score for option 2 and 3

    # Check error of the command line
    if len(sys.argv) != 4 and len(sys.argv) != 5:
        raise TypeError('Error! Verify the number of arguments, must be equal to 4 or 5')

    # OPTION 1: comparison of a sequence with all the database
    elif sys.argv[1] == '1':
        print(len(sys.argv))
        # Managing errors
        if not ospath.exists(sys.argv[2]):
            raise TypeError('The sequence path does not exist')
        elif sys.argv[3] != "matk" and sys.argv[3] != "psba-trnh" and sys.argv[3] != "rbcl" and sys.argv[3] != "its":
            raise TypeError("The gene name must be one of the following: 'matk', 'psba-trnh', 'rbcl', 'its'")
        else:
            start = time.time()
            # In case no family name was provided
            if len(sys.argv) == 4:
                database_path =  ospath.abspath(os.getcwd() + "/Database/" + sys.argv[3]  + ".fasta")
                db, seq =  import_data(database_path, sys.argv[2], option = '1')
                

            # In case a family name is provided
            else: 
                family_db = ospath.abspath(os.getcwd() + "/Database_by_family/" + sys.argv[4] + "/" + sys.argv[3]+ ".csv")
                if not ospath.exists(sys.argv[2]):
                    raise TypeError('The family specific csv file path does not exist')
                if ospath.exists(family_db):
                    db, seq =  import_data(family_db, sys.argv[2] , option = '1')
            
            print("the test sequence is:")
            print(seq)
            result = align(db,seq,"bioalign")
            end = time.time()
            print(result.head(10))
            print(f'Temps d\'exécution : ', end-start, 's')
        

    # OPTION 2: Comparison between 2 sequence using global_align
    elif sys.argv[1] == '2':
        if not ospath.exists(sys.argv[2]):
            raise TypeError('The sequence path does not exist')
        if not ospath.exists(sys.argv[3]):
            raise TypeError('The template sequence path does not exist')
        else:
            start=time.time()
            db, seq =  import_data(sys.argv[2], sys.argv[3], option = '2') 
            align(db, seq, "global_align")
            end = time.time()
            print(f'Temps d\'exécution : ', end-start, 's')
    
    # OPTION 3: Comparison between 2 sequence using biopython
    elif sys.argv[1] == '3':
        if not ospath.exists(sys.argv[2]):
            raise TypeError('The sequence path does not exist')
        if not ospath.exists(sys.argv[3]):
            raise TypeError('The template sequence path does not exist')
        else:
            start=time.time()
            db, seq =  import_data(sys.argv[2], sys.argv[3], option = '2') 
            align(db, seq,"bioalign", option ="A")
            end = time.time()
            print(f'Temps d\'exécution : ', end-start, 's')

    else:
        raise TypeError("[option] must be equal to '1', '2' or '3'")
