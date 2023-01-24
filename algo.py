from Bio import pairwise2
import numpy as np

def bio_align(S, T, opt = 'S'): 
    # Needleman-Wunsch algorithm using the Biopython library
    # Option 'S' : Return only the score of the alignment
    # Option 'A' : Return the score and the alignment
    # S = Sequence 1
    # T = Sequence 2

    if opt == 'S': 
        score = pairwise2.align.globalms(S, T, 2, -1, -0.5, -0.1, score_only=True)
        return score
    elif opt == 'A': 
        alignment = pairwise2.align.globalms(S, T, 2, -1, -0.5, -0.1)
        return pairwise2.format_alignment(*alignment[0])
    return 0

def global_align(S, T, match = 1, mismatch = 1, gap = 0):
    # Needleman-Wunsch algorithm coding from scrath (less performant)
    # S = Sequence 1
    # T = Sequence 2

    len_S = len(S)
    len_T = len(T)

    # Init map
    map = np.zeros((len_S + 1, len_T + 1))  # Init map with all value at 0 
    map[:,0] = np.linspace(0,   len_S * gap,   len_S + 1) # Init first line
    map[0,:] = np.linspace(0, -len_T * gap, len_T + 1) # Init first column
    
    # Pointers to trace through an optimal aligment.
    Opti_map = np.zeros((len_S + 1, len_T + 1))  # Init Opti_map with all value at 0 
    Opti_map[:,0] = 3 # Init first line
    Opti_map[0,:] = 4 # Init first column

    # Temporary scores.
    tmp = np.zeros(3)
    for i in range  (len_S):
        for j in range(len_T):
            if S[i] == T[j]: 
                tmp[0] = map[i,j] + match
            else:
                tmp[0] = map[i,j] - mismatch
            tmp[1] = map[i,j+1] - gap
            tmp[2] = map[i+1,j] - gap
            tmax = np.max(tmp)
            map[i+1,j+1] = tmax
            if tmp[0] == tmax:
                Opti_map[i+1,j+1] += 2
            if tmp[1] == tmax:
                Opti_map[i+1,j+1] += 3
            if tmp[2] == tmax:
                Opti_map[i+1,j+1] += 4
        
    # Trace through an optimal alignment.
    i = len_S
    j = len_T
    rS = []
    rT = []
    while i > 0 or j > 0:
        if Opti_map[i,j] in [2, 5, 6, 9]:
            rS.append(S[i-1])
            rT.append(T[j-1])
            i -= 1
            j -= 1
        elif Opti_map[i,j] in [3, 5, 7, 9]:
            rS.append(S[i-1])
            rT.append('-')
            i -= 1
        elif Opti_map[i,j] in [4, 6, 7, 9]:
            rS.append('-')
            rT.append(T[j-1])
            j -= 1
    
    # Reverse the strings.
    rS = ''.join(rS)[::-1]
    rT = ''.join(rT)[::-1]
    score = score_align(rS, rT)
    return '\n'.join([rS, rT]), score

def score_align(S, T):
    #Compute the score of the alignment for global_align
    
    score = 0
    size = 0
    for i in range (len(S)): 
        if ((S[i] == '-' or T[i] == '-') and (S[i-1] != '-' and T[i-1] != '-')):
            size += 1
        if S[i] == T[i]: 
            size += 1
            score += 1
    accuracy = score/size*100
    print("Alignment accuracy : ", accuracy, "%")
    return accuracy

def align(db, seq, algo, option = 'S'):
    #Compute the alignment and the score using the different function above
    #algo 'biopython': use bio_align()
            #Option 'S' : Return only the best score alignment.
            #Option 'A' : Return the alignement with the score. 
    #algo 'global_algng': use global_align()
            #Return the alignment and the score.
    
    if algo == "bioalign": 
        if option == 'S':
            score = 0
            species = ""
            for key in db.keys(): 
                S = db[key]
                tmp = bio_align(S, seq)
                if tmp > score: 
                    score = tmp
                    species = key
            score = score/(2*len(db[species]))*100
            print("Alignment score : ", score, "%")
            print("Species : ", species)

        elif option == 'A':
            result = bio_align(db, seq, opt='A')
            print(result)

    elif algo == "global_align":
        result = global_align(db, seq)
        print(result)
