import pairwise_alignment

class Node:
    def __init__(self, i, name,  seq):
        self.parent = i
        self.rank = 1
        self.count = 1
        self.index = i
        self.seq_names = [name]
        self.seqs = [seq]

def find(nodes, i):
    if nodes[i].parent != i:
        nodes[i].parent = find(nodes, nodes[i].parent)
    return nodes[i].parent

def join(nodes, i, j, score_matrix):
    x = nodes[find(nodes, i)]
    y = nodes[find(nodes, j)]
    if x.index == y.index:
        return -1, -1
    if x.rank == y.rank:
        x.rank += 1

    if x.rank < y.rank:
        x.parent = y.index
        y.count += x.count
        y.seq_names += x.seq_names
        y.seqs = pairwise_alignment.profile_alignment(y.seqs, x.seqs, score_matrix)
        res1 = y.index
        res2 = x.index
    else:
        y.parent = x.index
        x.count += y.count
        x.seq_names += y.seq_names
        x.seqs = pairwise_alignment.profile_alignment(x.seqs, y.seqs, score_matrix)
        res1 = x.index
        res2 = y.index

    nodes[x.index] = x
    nodes[y.index] = y
    return res1, res2 #first will present new cluster
