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

    a = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)] #initialise scoring matrix to 0
    for i in range(1, len(seq1)+1): #initialise left column with gap penalty
        a[i][0] = -i*8
    for j in range(1, len(seq2)+1): #initialise top row with gap penalty
        a[0][j] = -j*8
        
  
    print(a)

    #raise NotImplementedError()

def scoring_function(aa_i,aa_j):
    blosum62 = substitution_matrices.load("BLOSUM62")
    score = blosum62[aa_i][aa_j]
    return (score)

seq1 = "AEMODOPOULOS"
seq2 = "KLPSRTMNE"
global_alignment(seq1,seq2,scoring_function)
