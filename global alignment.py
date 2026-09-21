#seq 1 = #read from the sars_cov_2.fa file and extract characters 21562 -> 25384 from the sequence
#seq 2 = #use SeqIO.read to read the 'S' gene or product field from the 20 human coronaviruses
from Bio.Align import substitution_matrices

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """

    #define gap penalty
    d = 8
    #initialise pointer array
    pointer = [[None] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]

    #initialise scoring matrix to 0
    a = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)] 
    #initialise left column with gap penalty
    for i in range(1, len(seq1)+1): 
        a[i][0] = -i*8
        pointer[i][0] = (i-1,0)
    #initialise top row with gap penalty
    for j in range(1, len(seq2)+1): 
        a[0][j] = -j*8
        pointer[0][j] = (0,j-1)

    #perform alignment
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
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
    i = len(seq1) 
    j = len(seq2)
    k = 0
    identity = 0
    seq1a = ""
    seq2a = ""
    while i > 0 or j > 0:
        ip, jp = pointer[i][j]
        if ip == i:
            #gap in y
            seq1a += "-"
            seq2a += seq2[j-1]
        elif jp == j:
            #gap in x
            seq1a += seq1[i-1]
            seq2a += "-"
        else: 
            #match
            if seq1[i-1] == seq2[j-1]: identity += 1
            seq1a += seq1[i-1]
            seq2a += seq2[j-1]

        i, j = ip, jp
        k += 1

    id_score = 100*identity/k 

    seq1a = seq1a[::-1]
    seq2a = seq2a[::-1]

    return seq1a, seq2a, id_score

def scoring_function(aa_i,aa_j):
    blosum62 = substitution_matrices.load("BLOSUM62")
    score = blosum62[aa_i][aa_j]
    return (score)

seq1 = "AEMGDGPGILGS"
seq2 = "AEMVLIGDGILPGAV"
seq1a, seq2a, id_score = global_alignment(seq1,seq2,scoring_function)
print(f"{seq1a}\n{seq2a}\n{id_score}")
