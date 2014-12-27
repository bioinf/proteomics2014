def seq_alignment(s1, s2, score_matrix):
    m = len(s1);
    n = len(s2);
    d = [[0 for _ in range(0, m + 1)] for _ in range(0, n + 1)]
    d[0][0] = score_matrix['*']['*']
    for i in range(1, n + 1):
        d[i][0] = score_matrix[s2[i - 1]]['*']
    for j in range(1, m + 1):
        d[0][j] = score_matrix[s1[j - 1]]['*']
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d[i][j] = max(d[i - 1][j] + score_matrix[s1[j - 1]]['*'], d[i][j - 1] + score_matrix[s2[i - 1]]['*'], d[i - 1][j - 1] + score_matrix[s1[j - 1]][s2[i - 1]])
    s1_aligned = ""
    s2_aligned = ""
    i = n
    j = m
    while i > 0 and j > 0:
        if d[i][j] == (d[i - 1][j - 1] + score_matrix[s1[j - 1]][s2[i - 1]]):
            s1_aligned = s1[j - 1] + s1_aligned
            s2_aligned = s2[i - 1] + s2_aligned
            i -= 1
            j -= 1
        elif d[i][j] == d[i - 1][j] + score_matrix[s1[j - 1]]['*']:
            s1_aligned = '*' + s1_aligned
            s2_aligned = s2[i - 1] + s2_aligned
            i -= 1
        else:
            s1_aligned = s1[j - 1] + s1_aligned
            s2_aligned = '*' + s2_aligned
            j -= 1
    if i == 0:
        while j > 0:
            s1_aligned = s1[j - 1] + s1_aligned
            s2_aligned = '*' + s2_aligned
            j -= 1
    if j == 0:
        while i > 0:
            s1_aligned = '*' + s1_aligned
            s2_aligned = s2[i - 1] + s2_aligned
            i -= 1
    return d[n][m]

def make_profile_by_alignment(alignment, score_matrix):
    profile = {}
    n = len(alignment[0])
    m = len(alignment)
    for acid in list(score_matrix.keys()):
        profile[acid] = []
        for i in range(0, n):
            count = 0
            for j in range(0, m):
                if alignment[j][i] == acid:
                    count += 1
            profile[acid].append(float(count)/m)
    return profile

def profile_alignment(alignment1, alignment2, score_matrix):
    profile1 = make_profile_by_alignment(alignment1, score_matrix)
    profile2 = make_profile_by_alignment(alignment2, score_matrix)
    m = len(alignment1[0]);
    n = len(alignment2[0]);
    d = [[0 for _ in range(0, m + 1)] for _ in range(0, n + 1)]
    d[0][0] = score_matrix['*']['*']
    for i in range(1, n + 1):
        sum = 0
        for acid in list(score_matrix.keys()):
            sum += profile2[acid][i - 1] * score_matrix[acid]['*']
        d[i][0] = sum
    for j in range(1, m + 1):
        sum = 0
        for acid in list(score_matrix.keys()):
            sum += profile1[acid][j - 1] * score_matrix[acid]['*']
        d[0][j] = sum
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            gap_in_profile2_score = d[i][j - 1]
            for acid in list(score_matrix.keys()):
                gap_in_profile2_score += profile1[acid][j - 1] * score_matrix[acid]['*']

            gap_in_profile1_score = d[i - 1][j]
            for acid in list(score_matrix.keys()):
                gap_in_profile1_score += profile2[acid][i - 1] * score_matrix[acid]['*']

            match_score = d[i - 1][j - 1]
            for acid1 in list(score_matrix.keys()):
                for acid2 in list(score_matrix.keys()):
                    match_score += profile1[acid1][j - 1] * profile2[acid2][i - 1] * score_matrix[acid1][acid2]

            d[i][j] = max(gap_in_profile1_score, gap_in_profile2_score, match_score)


    alignment1_aligned = ["" for _ in range(0, len(alignment1))]
    alignment2_aligned = ["" for _ in range(0, len(alignment2))]
    i = n
    j = m
    while i > 0 and j > 0:
        gap_in_profile1_score = d[i - 1][j]
        for acid in list(score_matrix.keys()):
            gap_in_profile1_score += profile2[acid][i - 1] * score_matrix[acid]['*']

        match_score = d[i - 1][j - 1]
        for acid1 in list(score_matrix.keys()):
            for acid2 in list(score_matrix.keys()):
                match_score += profile1[acid1][j - 1] * profile2[acid2][i - 1] * score_matrix[acid1][acid2]

        if d[i][j] == match_score:
            for k in range(0, len(alignment1)):
                alignment1_aligned[k] = alignment1[k][j - 1] + alignment1_aligned[k]
            for k in range(0, len(alignment2)):
                alignment2_aligned[k] = alignment2[k][i - 1] + alignment2_aligned[k]
            i -= 1
            j -= 1

        elif d[i][j] == gap_in_profile1_score:
            for k in range(0, len(alignment1)):
                alignment1_aligned[k] = '*' + alignment1_aligned[k]
            for k in range(0, len(alignment2)):
                alignment2_aligned[k] = alignment2[k][i - 1] + alignment2_aligned[k]
            i -= 1
        else:
            for k in range(0, len(alignment1)):
                alignment1_aligned[k] = alignment1[k][j - 1] + alignment1_aligned[k]
            for k in range(0, len(alignment2)):
                alignment2_aligned[k] = '*' + alignment2_aligned[k]
            j -= 1
    if i == 0:
        while j > 0:
            for k in range(0, len(alignment1)):
                alignment1_aligned[k] = alignment1[k][j - 1] + alignment1_aligned[k]
            for k in range(0, len(alignment2)):
                alignment2_aligned[k] = '*' + alignment2_aligned[k]
            j -= 1
    if j == 0:
        while i > 0:
            for k in range(0, len(alignment1)):
                alignment1_aligned[k] = '*' + alignment1_aligned[k]
            for k in range(0, len(alignment2)):
                alignment2_aligned[k] = alignment2[k][i - 1] + alignment2_aligned[k]
            i -= 1

    alignment = alignment1_aligned + alignment2_aligned
    return alignment