def read_fasta(filename):
    names, seqs = [], []
    with open(filename, 'r') as f:
        name = "trash"
        seq = "trash"
        for line in f.readlines():
            if line[0] == '>':
                names.append(name)
                seqs.append(seq)
                name = line.strip()
                seq = ""
            else:
                seq += line.strip()
    names.append(name)
    seqs.append(seq)
    return names[1:], seqs[1:]

def read_matrix(filename):
    matix = {}
    with open(filename, 'r') as file:
        line = file.readline()
        while line[0] =='#':
            line = file.readline()
        amino_acids = line.split()
        line = file.readline()
        while line:
            blocks = line.split()
            matix[blocks[0]] = {}
            for i in range(0, len(amino_acids)):
                matix[blocks[0]][amino_acids[i]] = int(blocks[i + 1])
            line = file.readline()
    return matix
