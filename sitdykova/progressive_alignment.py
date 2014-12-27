import sys
import readers
import pairwise_alignment
import pgma
import neighbor_joining

if len(sys.argv) < 4:
    print("Usage: python progressive_alignment.py <seqs.fasta> <score_matrix> <upgma|wpgma|nj>")
    exit()
else:
    fasta_filename = sys.argv[1]
    matrix_filename = sys.argv[2]
    tree_type = sys.argv[3]
    print("Read fasta...")
    names, seqs = readers.read_fasta(fasta_filename)
    print("Read score matrix...")
    score_matrix = readers.read_matrix(matrix_filename)
    print("Align sequences...")
    if tree_type == "wpgma":
        names, seqs = pgma.pgma(names, seqs, score_matrix, 'w')
    elif tree_type == "upgma":
        names, seqs = pgma.pgma(names, seqs, score_matrix, 'u')
    elif tree_type == "nj":
        names, seqs = neighbor_joining.neigbor_joining(names, seqs, score_matrix)
    else:
        print("Error: Unkonwn option. Choose upgma, wpgma or nj.")
        exit()

    out_filename = fasta_filename.split('.')[0] + "_aligned.fasta"

    with open(out_filename, 'w') as out:
        for i in range(0, len(names)):
            out.write(names[i] + "\n")
            out.write(seqs[i] + "\n")

    print("Alignment saved if file " + out_filename)
