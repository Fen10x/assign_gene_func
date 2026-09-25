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

from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load("BLOSUM62") #load in the BLOSUM62 substitution matrix
def scoring_function(aa_i,aa_j):
    return (blosum62[aa_i][aa_j])

#scrape protein sequence from sars_cov_2.fa
from Bio import SeqIO
record = SeqIO.read("data/sars_cov_2.fa", "fasta")
dna = record.seq[421:667].reverse_complement() #get the reverse complement of the sequence
sars_cov_2 = dna.translate()

#scrape protein sequence from NCBI accession codes
accession_codes = {
    # 6 known human coronaviruses
    "Human-SARS": "NC_004718",
    "Bat-CoV RaTG13": "MN996532",
    "Pangolin-CoV MP789": "MT121216",
}

proteins = []
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
Entrez.email = "z5308203@ad.unsw.edu.au"
for name, accession in accession_codes.items():
    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="gb",
        retmode="text"
    )

    protein = SeqIO.read(handle, "genbank")
    handle.close()

    for feature in protein.features:
        if feature.type != "CDS":
            continue

        protein = SeqRecord(feature.extract(record.seq).translate())
        protein.id = name
        protein.description = feature.qualifiers.get("product", ["unknown"])[0]
        proteins.append(protein)

for protein in proteins:
    a, b, c = local_alignment(sars_cov_2, protein.seq, scoring_function)
    print(f"{protein.id}: {protein.description}, {c}, length: {len(a)}")
