def get_middle_score(cluster_set, score_matrix, i, j):
    return (score_matrix[i][j] * (len(cluster_set.existing) - 2)
            - sum(score_matrix[i][k] for k in cluster_set.existing)
            - sum(score_matrix[j][k] for k in cluster_set.existing))


def calc_score(score_matrix, i1, i2, k, _):
    if k == i1 or k == i2:
        return score_matrix[i1][k]
    return (score_matrix[i1][k] + score_matrix[i2][k] - score_matrix[i1][i2]) / 2
