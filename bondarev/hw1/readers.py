class Sequence:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def __getitem__(self, pos):
        return self.seq[pos]

    def __len__(self):
        return len(self.seq)


class ScoreMatrix(object):
    def __init__(self, matrix):
        self.matrix = matrix
        self.symbols = list(matrix.keys())

    def __getitem__(self, key):
        i, j = (x.upper() for x in key)
        return self.matrix[i][j]


def read_fasta(filename):
    names, seqs = [], []
    seq = []
    with open(filename, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                names.append(line[1:-1])  # remove > and \n
                if seq:
                    seqs.append(''.join(seq))
                seq = []
            else:
                seq.append(line.strip())
    if seq:
        seqs.append(''.join(seq))
    return (Sequence(name, seq) for name, seq in zip(names, seqs))


def read_score_matrix(filename):
    matrix = {}

    with open(filename, 'r') as matrix_file:
        line = matrix_file.readline()

        # skip comments
        while line.startswith('#'):
            line = matrix_file.readline()

        column_symbols = tuple(symbol.upper() if symbol != '*' else '-' for symbol in line.split())

        line = matrix_file.readline()
        while line:
            values = line.split()
            row_symbol = values[0].upper() if values[0] != '*' else '-'
            matrix[row_symbol] = {column_symbol: int(values[i + 1])
                                  for i, column_symbol in enumerate(column_symbols)}
            line = matrix_file.readline()
    return ScoreMatrix(matrix)
