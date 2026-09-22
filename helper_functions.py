def global_alignment(seq1, seq2, scoring_function):
    #define variables
    d = 8
    n = len(seq1)
    m = len(seq2)
    #initialise pointer array
    pointer = [[None] * (m + 1) for _ in range(n + 1)]

    #initialise scoring matrix to 0
    a = [[0] * (m + 1) for _ in range(n + 1)] 
    #initialise left column with gap penalty
    for i in range(1, n+1): 
        a[i][0] = -i*8
        pointer[i][0] = (i-1,0)
    #initialise top row with gap penalty
    for j in range(1, m+1): 
        a[0][j] = -j*8
        pointer[0][j] = (0,j-1)

    #perform alignment
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = a[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
            gap_x = a[i-1][j] - d
            gap_y = a[i][j-1] - d 
            a[i][j] = max(match, gap_x, gap_y)
            if a[i][j] == match:
                pointer[i][j] = (i-1,j-1)
            elif a[i][j] == gap_x:
                pointer[i][j] = (i-1,j)
            else:
                pointer[i][j] = (i,j-1)
    
    #traceback to construct alignment
    i = n
    j = m
    k = 0
    identity = 0
    seq1a = []
    seq2a = []
    while i > 0 or j > 0:
        ip, jp = pointer[i][j]
        if ip == i:
            #gap in y
            seq1a.append("-")
            seq2a.append(seq2[j-1])
        elif jp == j:
            #gap in x
            seq1a.append(seq1[i-1])
            seq2a.append("-")
        else: 
            #match
            if seq1[i-1] == seq2[j-1]: identity += 1
            seq1a.append(seq1[i-1])
            seq2a.append(seq2[j-1])

        i, j = ip, jp
        k += 1

    id_score = 100*identity/k 

    seq1a = "".join(reversed(seq1a))
    seq2a = "".join(reversed(seq2a))

    return seq1a, seq2a, id_score

def local_alignment(seq1, seq2, scoring_function):
    #define variables
    d = 8
    n = len(seq1)
    m = len(seq2)
    #initialise pointer array
    pointer = [[None] * (m + 1) for _ in range(n + 1)]

    #initialise scoring matrix to 0 (does left column and top row automatically)
    a = [[0] * (m + 1) for _ in range(n + 1)]

    #perform alignment
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = a[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
            gap_x = a[i-1][j] - d
            gap_y = a[i][j-1] - d 
            a[i][j] = max(0, match, gap_x, gap_y)
            if a[i][j] == 0:
                pointer[i][j] = None
            elif a[i][j] == match:
                pointer[i][j] = (i-1,j-1)
            elif a[i][j] == gap_x:
                pointer[i][j] = (i-1,j)
            else:
                pointer[i][j] = (i,j-1)

    #traceback to construct alignment
    #find the maximum score in the scoring matrix
    max_score = 0
    max_i = 0
    max_j = 0
    for i in range(n+1):
        for j in range(m+1):
            if a[i][j] > max_score:
                max_score = a[i][j]
                max_i = i
                max_j = j

    #start traceback from the cell with the maximum score
    i = max_i
    j = max_j
    k = 0
    identity = 0
    seq1a = []
    seq2a = []
    while i > 0 and j > 0 and a[i][j] != 0:
        ip, jp = pointer[i][j]
        if ip == i:
            #gap in y
            seq1a.append("-")
            seq2a.append(seq2[j-1])
        elif jp == j:
            #gap in x
            seq1a.append(seq1[i-1])
            seq2a.append("-")
        else: 
            #match
            if seq1[i-1] == seq2[j-1]: identity += 1
            seq1a.append(seq1[i-1])
            seq2a.append(seq2[j-1])

        i, j = ip, jp
        k += 1

    id_score = 100*identity/k

    seq1a = "".join(reversed(seq1a))
    seq2a = "".join(reversed(seq2a))

    return seq1a, seq2a, id_score

#load the BLOSUM62 substitution matrix outside of scoring_function for efficiency
from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load("BLOSUM62")
def scoring_function_simple(aa_i,aa_j):
    return (blosum62[aa_i][aa_j])