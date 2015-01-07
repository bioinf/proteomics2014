#!/usr/bin/python

import sys
import logging
import ScoreMatrix
import Clustering

def read_input(filename):
    names, seqs = [], []

    with open(filename, 'r') as input:
        name, seq = "", ""
        for line in input.readlines():
            if line[0] == '>':
                names.append(name)
                seqs.append(seq)
                name = line.strip()[1:]
                seq = ""
            else:
                seq += line.strip()

    names.append(name)
    seqs.append(seq)
    return names[1:], seqs[1:]

def write_ouput(filename, names, seqs):
    with open(filename, 'w') as out:
        for i in range(len(names)):
            out.write("> " + names[i] + "\n")
            out.write(seqs[i].replace('*', '-') + "\n")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    if (len(sys.argv) < 5):
        sys.stderr.write("USAGE: main.py <matrix.txt> <in.fasta> <NJ | WPGMA | UPGMA> <out.fasta>\n")
        sys.exit(1)

    logging.info("START DO JOB")
    logging.info("Read matrix")
    blosum = ScoreMatrix.Blosum(sys.argv[1])

    logging.info("Read input file")
    names, seqs = read_input(sys.argv[2])

    methods = {"NJ": ("Neighbor joining", Clustering.neigbor_joining),
               "WPGMA": ("WPGMA", Clustering.wpgma),
               "UPGMA": ("UPGMA", Clustering.upgma)}

    logging.info("Determine method for build guide tree")
    if sys.argv[3] in methods:
        name, cluster_func = methods[sys.argv[3]]
        logging.info("Build guide tree with " + name)
    else:
        logging.error("Unknown method " + sys.argv[3])
        sys.exit(1)

    logging.info("Do job")
    result_names, result_seqs = cluster_func(names, seqs, blosum)

    logging.info("Write result alignment")
    write_ouput(sys.argv[4], result_names, result_seqs)
    logging.info("FINISH JOB")
