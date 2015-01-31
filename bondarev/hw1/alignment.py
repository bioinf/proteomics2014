def get_global_alignment_score(seq1, seq2, score_matrix):
    m, n = len(seq1), len(seq2)
    d = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(n):
        d[i + 1][0] = score_matrix[seq2[i], '-'] * (i + 1)
    for j in range(m):
        d[0][j + 1] = score_matrix[seq1[j], '-'] * (j + 1)

    for i in range(n):
        for j in range(m):
            horizontal = d[i][j + 1] + score_matrix[seq1[j], '-']
            vertical = d[i + 1][j] + score_matrix[seq2[i], '-']
            diagonal = d[i][j] + score_matrix[seq1[j], seq2[i]]
            d[i + 1][j + 1] = max(horizontal, vertical, diagonal)

    return d[-1][-1]
