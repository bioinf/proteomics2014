#!/usr/bin/env python3

import sys
import logging

from readers import read_fasta, read_score_matrix
import clustering


METHODS = {
    'nj': clustering.neigbor_joining_cluster,
    'upgma': clustering.upgma_cluster,
    'wpgma': clustering.wpgma_cluster
}
logging.basicConfig(level=logging.INFO)


def usage():
    return ('Usage:\n'
            '  main.py <fasta> <score_matrix> <clustering type:({})>\n'.format('|'.join(METHODS.keys())))


def main():
    if len(sys.argv) != 4:
        logging.error('Wrong count of arguments')
        print(usage())
        sys.exit()

    _, fasta_filename, score_matrix_filename, method = sys.argv
    method = method.lower()

    if method not in METHODS:
        logging.error('Wrong method for cluster')
        print(usage())
        sys.exit()

    logging.info('Load fasta file')
    sequences = read_fasta(fasta_filename)
    logging.info('Finish load fasta file')

    logging.info('Load score matrix')
    score_matrix = read_score_matrix(score_matrix_filename)
    logging.info('Finish load score matrix')

    logging.info('Start progressive alignment with method ' + method)
    result_sequences = METHODS[method](sequences, score_matrix)
    logging.info('Finish progressive alignment')

    for sequence in result_sequences:
        print('>' + sequence.name)
        print(sequence.seq)

    logging.info('Done.')


if __name__ == '__main__':
    main()