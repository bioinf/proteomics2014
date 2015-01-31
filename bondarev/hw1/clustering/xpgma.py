def get_middle_score(scores, i, j, _):
    return scores[i][j]


def upgma_calc_score(scores, u, v, w, cluster_sizes):
    if u == w:
        return scores[u][w]
    return ((scores[u][w] * cluster_sizes[u] +
             scores[v][w] * cluster_sizes[v]) /
            (cluster_sizes[u] + cluster_sizes[v]))


def wpgma_calc_score(scores, i1, i2, k, _):
    if i1 == k:
        return scores[i1][k]
    return (scores[i1][k] + scores[i2][k]) / 2
