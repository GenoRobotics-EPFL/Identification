from utils import *
from algo import *
import pandas as pd
import time
import sys

if __name__ == "__main__":
    #Run all the programm
    #Run the 3 options (see the Readme file for more information on it)
    #function time : allow to see the runtime of each different algorithm 
    #print the result: DataFrame with the species and the score for option 1 and alignment and score for option 2 and 3

    #Check error of the command line
    if len(sys.argv) != 4:
        print("Error !")
        print("Please verify the number of argument")
        print("Number of argument must be egual to 3 : [Database file path], [Sequence file path], [option (1, 2 or 3)]")
    
    #Comparison of a sequence with all the database
    elif sys.argv[3] == '1':
        start = time.time()
        db, seq =  import_data(sys.argv[1], sys.argv[2], option = '1')
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

