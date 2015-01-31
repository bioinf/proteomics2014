from readers import Sequence


class Profile(object):
    def __init__(self, sequences, alphabet):
        self.length = len(sequences[0])  # all sequences in cluster have same length (by merge procedure)
        self.sequences = sequences
        self.alphabet = alphabet
        self.profile = [{} for _ in range(self.length)]  # table: column is position, row is symbol
        for pos in range(self.length):
            for symbol in alphabet:
                count = sum(sequence[pos].upper() == symbol for sequence in sequences)
                self.profile[pos][symbol] = count / len(sequences)

    def __getitem__(self, pos):
        return self.profile[pos]


def align(seqs1, seqs2, score_matrix):
        profile1 = Profile(seqs1, score_matrix.symbols)
        profile2 = Profile(seqs2, score_matrix.symbols)
        matrix = _get_dynamic_matrix(profile1, profile2, score_matrix)
        return _get_merged_sequences(profile1, profile2, matrix, score_matrix)


def _get_dynamic_matrix(profile1, profile2, score_matrix):
    cols, rows = profile1.length, profile2.length
    d = [[0] * (cols + 1) for _ in range(rows + 1)]

    for row in range(rows):
        d[row + 1][0] = d[row][0] + _get_gap_score(profile2, row, score_matrix)
    for col in range(cols):
        d[0][col + 1] = d[0][col] + _get_gap_score(profile1, col, score_matrix)

    for row in range(rows):
        for col in range(cols):
            horizontal = d[row][col + 1] + _get_gap_score(profile2, row, score_matrix)
            vertical = d[row + 1][col] + _get_gap_score(profile1, col, score_matrix)
            diagonal = d[row][col] + _get_match_score(profile1, col, profile2, row, score_matrix)
            d[row + 1][col + 1] = max(horizontal, vertical, diagonal)

    return d


def _get_match_score(profile1, pos1, profile2, pos2, score_matrix):
    return sum(profile1[pos1][symbol1] * profile2[pos2][symbol2] * score_matrix[symbol1, symbol2]
               for symbol1 in score_matrix.symbols
               for symbol2 in score_matrix.symbols)


def _get_gap_score(profile, pos, score_matrix):
    return sum(profile[pos][symbol] * score_matrix[symbol, '-'] for symbol in profile.alphabet)


def _get_merged_sequences(profile1, profile2, d, score_matrix):
    cnt1, cnt2 = len(profile1.sequences), len(profile2.sequences)
    aligned_seqs1 = [[] for _ in range(cnt1)]
    aligned_seqs2 = [[] for _ in range(cnt2)]
    row, col = profile2.length, profile1.length

    while row > 0 and col > 0:
        vertical_shift_score = d[row - 1][col] + _get_gap_score(profile2, row - 1, score_matrix)
        diagonal_shift_score = d[row - 1][col - 1] + _get_match_score(profile1, col - 1,
                                                                      profile2, row - 1, score_matrix)

        if d[row][col] == diagonal_shift_score:
            _add_column(aligned_seqs1, col, profile1.sequences)
            _add_column(aligned_seqs2, row, profile2.sequences)
            row -= 1
            col -= 1
        elif d[row][col] == vertical_shift_score:
            _add_gap(aligned_seqs1)
            _add_column(aligned_seqs2, row, profile2.sequences)
            row -= 1
        else:  # horizontal gap
            _add_column(aligned_seqs1, col, profile1.sequences)
            _add_gap(aligned_seqs2)
            col -= 1

    while col > 0:
        _add_column(aligned_seqs1, col, profile1.sequences)
        _add_gap(aligned_seqs2)
        col -= 1

    while row > 0:
        _add_gap(aligned_seqs1)
        _add_column(aligned_seqs2, row, profile2.sequences)
        row -= 1

    result_sequences = list(_make_sequences(aligned_seqs1, profile1))
    result_sequences.extend(_make_sequences(aligned_seqs2, profile2))
    return result_sequences


def _add_column(aligned_seqs, pos, reference_seqs):
    for aligned_seq, reference_seq in zip(aligned_seqs, reference_seqs):
        aligned_seq.append(reference_seq[pos - 1])


def _add_gap(aligned_seqs):
    for aligned_seq in aligned_seqs:
        aligned_seq.append('-')


def _make_sequences(aligned_seqs, profile):
    return (Sequence(profile_seq.name, ''.join(reversed(result_seq)))
            for profile_seq, result_seq in zip(profile.sequences, aligned_seqs))
