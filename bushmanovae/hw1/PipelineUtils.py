__author__ = 'letovesnoi'

import argparse
import sys

def getArgs():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description = "ASSESSMENT OF THE PROGRESSIVE ALIGNMENT %(prog)s\r\n\r\nUsage:\r\n"
                                                   "python %(prog)s --fasta FASTA --matrix MATRIX --type TYPE --out OUT\r\n",
                                     conflict_handler='resolve', prog=sys.argv[0])

    groupInputData = parser.add_argument_group('input data')
    groupInputData.add_argument('--fasta', help='FASTA file with sequences', type=str, required=True)
    groupInputData.add_argument('--matrix', help='TXT alignment matrix from NCBI FTP', type=str, required=True)
    groupInputData.add_argument('--type', help='type of clustering', type=str, required=True, choices=['NJ', 'UPGMA', 'WPGMA'])

    groupBasic = parser.add_argument_group('basic options')
    groupBasic.add_argument('--out', help='output file with alignment', type=str, required=True)

    args = parser.parse_args()

    return args