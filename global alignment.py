#seq 1 = #read from the sars_cov_2.fa file and extract characters 21562 -> 25384 from the sequence
#seq 2 = #use SeqIO.read to read the 'S' gene or product field from the 20 human coronaviruses
from Bio import SeqIO
from Bio.Align import substitution_matrices
from Bio import Entrez

def global_alignment(seq1, seq2, scoring_function):
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

#scrape protein sequence from sars_cov_2.fa
record = SeqIO.read("data/sars_cov_2.fa", "fasta")
dna = record.seq[21561:25384]
sars_cov_2 = dna.translate()

#scrape protein sequence from 
accession_codes = {
    # 6 known human coronaviruses
    "Human-SARS": "NC_004718",
    "Human-MERS": "NC_019843",
    "Human-HCoV-OC43": "NC_006213",
    "Human-HCoV-229E": "NC_002645",
    "Human-HCoV-NL63": "NC_005831",
    "Human-HCoV-HKU1": "NC_006577",
    
    # Bat
    "Bat-CoV MOP1": "EU420138",
    "Bat-CoV HKU8": "NC_010438",
    "Bat-CoV HKU2": "NC_009988",
    "Bat-CoV HKU5": "NC_009020",
    "Bat-CoV RaTG13": "MN996532",
    "Bat-CoV-ENT": "NC_003045",
    
    # Other animals
    "Hedgehog-CoV 2012-174/GER/2012": "NC_039207",
    "Pangolin-CoV MP789": "MT121216",
    "Rabbit-CoV HKU14": "NC_017083",
    "Duck-CoV isolate DK/GD/27/2014": "NC_048214",
    "Feline infectious peritonitis virus": "NC_002306",  # cat
    "Giraffe-CoV US/OH3/2003": "EF424623",
    "Murine-CoV MHV/BHKR_lab/USA/icA59_L94P/2012": "KF268338",  # mouse
    "Equine-CoV Obihiro12-2": "LC061274",  # horse
}
spike_proteins = []
Entrez.email = "z5308203@ad.unsw.edu.au"
for name, accession in accession_codes.items():
    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="gb",
        retmode="text"
    )

    record = SeqIO.read(handle, "genbank")
    handle.close()

    for feature in record.features:
        if feature.type != "CDS":
            continue

        gene = feature.qualifiers.get("gene", [""])[0]
        product = feature.qualifiers.get("product", [""])[0]

        if gene == "S" or "spike protein" in product.lower():
            record.seq = feature.extract(record.seq).translate()
            record.description = name
            spike_proteins.append(record)
            #print(f"{name}: found Spike ({len(record.seq)} bp)")
            break

for spike_protein in spike_proteins:
    seq1a, seq2a, id_score = global_alignment(sars_cov_2, spike_protein, scoring_function)
    print(f"{spike_protein.description}: {id_score}")

#seq1a, seq2a, id_score = global_alignment(sars_cov_2,seq2,scoring_function)