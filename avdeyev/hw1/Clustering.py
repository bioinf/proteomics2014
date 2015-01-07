import sys
import logging
from functools import partial

import alignment

logger = logging.getLogger()


class Cluster(object):
    def __init__(self, i, name, seq):
        self.ind = i
        self.names = [name]
        self.seqs = [seq]
        self.depth = 1
        self.parent = i


def pgma_modify_score(scores, i, j, indices):
    return scores[i][j]


def upgma_update_score(scores, i1, i2, k, count):
    if (i1 != k):
        return float(scores[i1][k] * count[i1] + scores[i2][k] * count[i2]) / (count[i1] + count[i2])
    else:
        return scores[i1][k]


def wpgma_update_score(scores, i1, i2, k, count):
    if (i1 != k):
        return float(scores[i1][k] + scores[i2][k]) / 2
    else:
        return scores[i1][k]


def neigbor_joining_modify_score(scores, i, j, indices):
    def my_sum(list):
        res = 0
        for i in indices:
            res += list[i]
        return res

    return scores[i][j] * (len(indices) - 2) - my_sum(scores[i]) - my_sum(scores[j])


def neigbor_joining_update_score(scores, i1, i2, k, count):
    if k != i1 and k != i2:
        return float(scores[i1][k] + scores[i2][k] - scores[i1][i2]) / 2
    else:
        return scores[i1][k]


def template_clustering(modify_score, update_score, names, seqs, score_matrix):
    n = len(seqs)
    clusters = [Cluster(i, names[i], seqs[i]) for i in range(n)]
    indices, count = set([i for i in range(n)]), [1 for _ in range(n)]
    scores = [[0 for _ in range(n)] for _ in range(n)]

    def find(ind):
        if clusters[ind].parent != ind:
            clusters[ind].parent = find(clusters[ind].parent)
        return clusters[ind].parent

    logger.info("Init score matrix")
    for i in range(n):
        for j in range(i + 1, n):
            scores[i][j] = scores[j][i] = alignment.sequence_alignment(seqs[i], seqs[j], score_matrix)
        scores[i][i] = -sys.maxsize


    for iter in range(0, n - 1):
        max_i, max_j, max = -1, -1, -sys.maxsize
        for i in range(n):
            if i not in indices:
                continue

            for j in range(i + 1, n):
                if j not in indices:
                    continue

                score = modify_score(scores, i, j, indices)
                if score > max:
                    max_i, max_j, max = i, j, score

        x1, x2 = clusters[find(max_i)], clusters[find(max_j)]
        logger.info("Start to join next two clusters " + str(x1.ind) + " " + str(x2.ind))
        new_i, new_j = merge(clusters, x1, x2, score_matrix)

        logger.info("Update score matrix")
        for k in range(n):
            scores[new_i][k] = scores[k][new_i] = update_score(scores, new_i, new_j, k, count)
        indices.remove(new_j)
        count[new_i] += count[new_j]

    return clusters[find(0)].names, clusters[find(0)].seqs

upgma = partial(template_clustering, pgma_modify_score, upgma_update_score)
wpgma = partial(template_clustering, pgma_modify_score, wpgma_update_score)
neigbor_joining = partial(template_clustering, neigbor_joining_modify_score, neigbor_joining_update_score)


def merge(clusters, x, y, score_matrix):
    if x.ind == y.ind:
        return -1, -1

    if x.depth == y.depth:
        x.depth += 1

    logger.info("Run profile alignment")
    if x.depth < y.depth:
        x.parent = y.ind
        y.names += x.names
        y.seqs = alignment.profile_alignment(y.seqs, x.seqs, score_matrix)
        res1, res2 = y.ind, x.ind
    else:
        y.parent = x.ind
        x.names += y.names
        x.seqs = alignment.profile_alignment(x.seqs, y.seqs, score_matrix)
        res1, res2 = x.ind, y.ind
    logger.info("Finish profile alignment")

    clusters[x.ind], clusters[y.ind] = x, y
    return res1, res2
