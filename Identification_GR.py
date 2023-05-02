from utils import *
from algo import *
import pandas as pd
import time
import sys
import os.path as ospath
import os

if __name__ == "__main__":
    #Run all the programm
    #Run the 3 options (see the Readme file for more information on it)
    #function time : allow to see the runtime of each different algorithm 
    #print the result: 
        #   DataFrame with the species and the score for option 1 
        #   Alignment and score for option 2 and 3

    #Check error of the command line
    if len(sys.argv) != 4 or len(sys.argv) != 5:
        print("Error !")
        print("Please verify the number of argument")
        print("Number of argument must be egual to 3 or 4:  [Sequence file path], [option (1, 2 or 3)], [Gene Name], [OPTIONAL: Family Name]" )
    
    #Comeison of a sequence with all the database
    elif sys.argv[2] == '1':
        start = time.time()
        if len(sys.argv) == 4:
            database_path =  ospath.abspath(os.getcwd() + "/Database/" + sys.argv[3])
            db, seq =  import_data(database_path, sys.argv[1], option = '1')
            ##code de alex
        else: 
            # boucler sur les noms des familles pour verifier si argv[4] est bien un nom de famille
            # sinon raise ValueError
            database_path =  ospath.abspath(os.getcwd() + "/Database/" + sys.argv[4] + "/" + sys.argv[3])
            db, seq =  import_data(database_path, sys.argv[1], option = '1')

        result = align(db,seq,"bioalign")
        end = time.time()
        print(result.head(10))
        print(f'Temps d\'exécution : ', end-start, 's')
    
    #Comparison between 2 sequence using global_align
    elif sys.argv[3] == '2':
        start=time.time()
        db, seq =  import_data(sys.argv[1], sys.argv[2], option = '2') 
        align(db, seq, "global_align")
        end = time.time()
        print(f'Temps d\'exécution : ', end-start, 's')
    
    #Comparison between 2 sequence using biopython
    elif sys.argv[3] == '3':
        start=time.time()
        db, seq =  import_data(sys.argv[1], sys.argv[2], option = '2') 
        align(db, seq,"bioalign", option ="A")
        end = time.time()
        print(f'Temps d\'exécution : ', end-start, 's')

