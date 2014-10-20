__author__ = 'letovesnoi'

def readFASTA(fileName):
    seq = []
    with open(fileName, 'r') as fin:
        temp = fin.readline().strip()
        while temp != '':
            temp = fin.readline().strip()
            seq.append('')
            while '>' not in temp:
                seq[-1] += temp
                if temp == '':
                    break
                temp = fin.readline().strip()
            seq[-1] = seq[-1].upper()
    return seq