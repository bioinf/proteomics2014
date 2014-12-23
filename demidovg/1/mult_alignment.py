import substitution
import proteins
import math
import sys
from alignment import *
from clusterization import *


def return_distance(matrix_of_distances, key1, key2):
    if (key1 in matrix_of_distances and key2 in matrix_of_distances[key1]):
        return matrix_of_distances[key1][key2]
    elif (key2 in matrix_of_distances and key1 in matrix_of_distances[key2]):
        return matrix_of_distances[key2][key1]
    else:
        return None

def main():
    filename = sys.argv[1]
    protfile = sys.argv[2]
    method = int(sys.argv[3])
    #filename = "BLOSUM62.txt"
    #protfile = "test.fasta"

    # consts
    gap_open = -8
    gap_extention = -1
    gap_costs = (gap_open, gap_extention)

    blosum = substitution.Substitution_matrix(filename)
    prots = proteins.Sequences(protfile)

    alphabet = list(blosum.matrix.iterkeys())
    alphabet.append("-")
    alphabet.remove("*")
    matrix_of_distances = defaultdict(dict)
    weights_own_alignments = {}
    for key1 in prots.seqs.iterkeys():
        prot1 = prots.seqs[key1]
        weights_own_alignments[key1] = align_seq_to_seq(prot1, prot1, blosum, gap_costs)[1]
    for key1 in prots.seqs.iterkeys():
        for key2 in prots.seqs.iterkeys():
            prot1 = prots.seqs[key1]
            prot2 = prots.seqs[key2]
            if key1 == key2:
                matrix_of_distances[key1][key2] = 0
                matrix_of_distances[key2][key1] = 0
            elif return_distance(matrix_of_distances, key1, key2) == None:
                arr = align_seq_to_seq(prot1, prot2, blosum, gap_costs)
                x = arr[1]
                matrix_of_distances[key1][key2] = 1 - (1.0 * align_seq_to_seq(prot1, prot2, blosum, gap_costs)[1] /
                                                         math.sqrt(weights_own_alignments[key1] * weights_own_alignments[key2]))
                matrix_of_distances[key2][key1] = matrix_of_distances[key1][key2]
            else:
                continue
    clusterize = Clusterizator(matrix_of_distances, prots)
    if method == 0:
        clusterize.nj()
    elif method == 1:
        clusterize.upgma()
    else:
        clusterize.wpgma()


main()