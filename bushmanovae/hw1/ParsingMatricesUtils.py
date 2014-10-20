__author__ = 'letovesnoi'

def parseMatrix(inputFile):
    matrix = {}
    with open(inputFile, 'r') as fin:
        line = fin.readline().strip()
        while '#' in line:
            line = fin.readline().strip()
        aminoAcids =line.split()
        line = fin.readline().strip()
        num = 0
        while line != '':
            line = line.split(' ')
            for symbols in line:
                if symbols == '':
                    line.remove(symbols)
            line = line[1:]
            for iAcid in range(len(aminoAcids)):
                acid = aminoAcids[iAcid]
                matrix[acid, aminoAcids[num]] = float(line[iAcid])
            num += 1
            line = fin.readline().strip()
    return matrix