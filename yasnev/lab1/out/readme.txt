--------------------------
DESCRIPTION
--------------------------

Progressive alignment

Input:
FASTA file with aminoacids, scoring matrix, clustering method

Output:
FASTA file with MSA

--------------------------
USAGE
--------------------------

Usage: java -jar lab1.jar <in.fasta> <matrix.txt> method=nj|upgma|wpgma
    <in.fasta>             - input file with aminoacid sequences in FASTA format
    <matrix.txt>           - scoring matrix, e.g. BLOSUM62
    method                 - clustering method (nj - neighbor-joining)

Options:
    -h, --help, ?           Show this help
    -o <file>               Output file
