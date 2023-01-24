from utils import *
from algo import *
import sys

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Error !")
        print("Please verify the number of argument")
        print("Number of argument must be egual to 3 : [Database file path], [Sequence file path], [option (1, 2 or 3)]")
    elif sys.argv[3] == '1':
        db, seq =  import_data(sys.argv[1], sys.argv[2], option = '1')
        align(db,seq,"bioalign")
    elif sys.argv[3] == '2':
        db, seq =  import_data(sys.argv[1], sys.argv[2], option = '2') 
        align(db, seq, "global_align")
    elif sys.argv[3] == '3':
        db, seq =  import_data(sys.argv[1], sys.argv[2], option = '2') 
        align(db, seq,"bioalign", option ="A")