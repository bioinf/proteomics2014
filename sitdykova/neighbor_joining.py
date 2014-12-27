import node
import pairwise_alignment

def sum_scores(scores_list, indices):
    res = 0
    for k, v in indices.items():
        res += scores_list[k]
    return res

def neigbor_joining(seq_names, seqs, score_matrix):
    n = len(seqs)
    nodes = [node.Node(i, seq_names[i], seqs[i]) for i in range(0, n)]
    count = [1 for _ in range(0, n)]
    scores = [[0 for _ in range(0, n)] for _ in range(0, n)]
    minus_infinity = -10000000000
    for i in range(0, n):
        for j in range(i + 1, n):
            score = pairwise_alignment.seq_alignment(seqs[i], seqs[j], score_matrix)
            scores[i][j] = score
            scores[j][i] = score
    for i in range(0, n):
        scores[i][i] = 0

    clusters = {}
    for i in range(0, n):
        clusters[i] = 1

    for lolo in range(0, n - 1):
        argmax1 = -1
        argmax2 = -1
        max = minus_infinity
        r = len(clusters)
        for i in range(0, n):
            if i in clusters:
                for j in range(i + 1, n):
                    if j in clusters:
                        score = scores[i][j] * (r - 2) - sum_scores(scores[i], clusters) - sum_scores(scores[j], clusters)
                        if score > max:
                            max = score
                            argmax1 = i
                            argmax2 = j

        i1, i2 = node.join(nodes, argmax1, argmax2, score_matrix) # first will present new cluster

        #recalc_scores
        for k in range(0, n):
            if k != i1 and k!= i2:
                new_score = float(scores[i1][k] + scores[i2][k] - scores[i1][i2]) / 2
                scores[i1][k] = new_score
                scores[k][i1] = new_score

        del clusters[i2]

        count[i1] += count[i2]

    root_index = node.find(nodes, 0)
    return nodes[root_index].seq_names, nodes[root_index].seqs
