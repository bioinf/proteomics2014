
class Profile(object):
    def __init__(self, alignment, score_matrix):
        self.profile = {}
        for ac in list(score_matrix.keys()):
            self.profile[ac] = []
            for i in range(len(alignment[0])):
                count = 0
                for j in range(len(alignment)):
                    if alignment[j][i] == ac:
                        count += 1
                self.profile[ac].append(float(count) / len(alignment))

    def get_gap_score(self, ind, score_matrix):
        res = 0
        for a in list(score_matrix.keys()):
            res += self.profile[a][ind - 1] * score_matrix.get_elem(a, '*')
        return res

    def get_match_score(self, ind_i, profile, ind_j, score_matrix):
        res = 0
        for a1 in list(score_matrix.keys()):
            for a2 in list(score_matrix.keys()):
                res += self.profile[a1][ind_i - 1] * profile.profile[a2][ind_j - 1] * score_matrix.get_elem(a1, a2)
        return res


def profile_alignment(s1, s2, score_matrix):
    profile1, profile2 = Profile(s1, score_matrix), Profile(s2, score_matrix)
    matrix = _build_matrix(len(s1[0]), len(s2[0]), profile1, profile2, score_matrix)
    return _get_answer(s1, s2, profile1, profile2, matrix, score_matrix)


def sequence_alignment(s1, s2, score_matrix):
    m, n = len(s1), len(s2)
    dyn = [[0 for _ in range(0, m + 1)] for _ in range(0, n + 1)]

    for i in range(1, n + 1):
        dyn[i][0] = score_matrix.get_elem(s2[i - 1], '*')

    for j in range(1, m + 1):
        dyn[0][j] = score_matrix.get_elem(s1[j - 1], '*')

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            gap_in_i = dyn[i - 1][j] + score_matrix.get_elem(s1[j - 1], '*')
            gap_in_j = dyn[i][j - 1] + score_matrix.get_elem(s2[i - 1], '*')
            match_score = dyn[i - 1][j - 1] + score_matrix.get_elem(s1[j - 1], s2[i - 1])
            dyn[i][j] = max(gap_in_i, gap_in_j, match_score)

    return dyn[n][m]


def _build_matrix(m, n, profile1, profile2, score_matrix):
    dyn = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    for i in range(1, n + 1):
        dyn[i][0] = profile2.get_gap_score(i, score_matrix)

    for j in range(1, m + 1):
        dyn[0][j] = profile1.get_gap_score(j, score_matrix)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            gap_in_i = dyn[i - 1][j] + profile2.get_gap_score(i, score_matrix)
            gap_in_j = dyn[i][j - 1] + profile1.get_gap_score(j, score_matrix)
            match_score = dyn[i - 1][j - 1] + profile1.get_match_score(j, profile2, i, score_matrix)
            dyn[i][j] = max(gap_in_i, gap_in_j, match_score)

    return dyn


def _get_answer(s1, s2, profile1, profile2, dyn, score_matrix):
    alignment1, alignment2 = ["" for _ in range(len(s1))], ["" for _ in range(len(s2))]
    i, j = len(s2[0]), len(s1[0])

    def update_alignment(alignment, l, ind, s):
        for k in range(l):
            alignment[k] = s[k][ind - 1] + alignment[k]

    def update_by_sym(alignment, l, sym):
        for k in range(l):
            alignment[k] = sym + alignment[k]

    while i > 0 and j > 0:
        gap_in_profile1_score = dyn[i - 1][j] + profile2.get_gap_score(i, score_matrix)
        match_score = dyn[i - 1][j - 1] + profile1.get_match_score(j, profile2, i, score_matrix)

        if dyn[i][j] == match_score:
            update_alignment(alignment1, len(s1), j, s1)
            update_alignment(alignment2, len(s2), i, s2)
            i -= 1
            j -= 1
        elif dyn[i][j] == gap_in_profile1_score:
            update_by_sym(alignment1, len(s1), '*')
            update_alignment(alignment2, len(s2), i, s2)
            i -= 1
        else:
            update_alignment(alignment1, len(s1), j, s1)
            update_by_sym(alignment2, len(s2), '*')
            j -= 1

    while j > 0:
        update_alignment(alignment1, len(s1), j, s1)
        update_by_sym(alignment2, len(s2), '*')
        j -= 1

    while i > 0:
        update_by_sym(alignment1, len(s1), '*')
        update_alignment(alignment2, len(s2), i, s2)
        i -= 1

    return (alignment1 + alignment2)