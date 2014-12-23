__author__ = 'german'

from collections import defaultdict
import align_table

def create_profile(sequences, alphabet):
    string_flag = True
    profile = defaultdict(list)
    len_of_align = 0
    if isinstance(sequences, str):
        len_of_align = len(sequences)
    else:
        len_of_align = len(sequences[0])
        string_flag = False

    for letter in alphabet:
        profile[letter] = [0 for i in range(len_of_align)]
    for i in range(len_of_align):
        for letter in alphabet:
            if string_flag:
                if sequences[i] == letter:
                    profile[letter][i] += 1
            else:
                for sequence in sequences:
                    if sequence[i] == letter:
                        profile[letter][i] += 1
    return profile

def find_backward_alignment(seq1, seq2, backward_matrix):
    cell = backward_matrix[len(seq1)][len(seq2)]

    i = len(seq1)
    j = len(seq2)
    first_string = ""
    second_string = ""

    while ( (i > 0)) and ( (j > 0)):
        new_i = cell[0]
        new_j = cell[1]
        first_letter = None
        second_letter = None
        if i == new_i:
            first_letter = "-"
        else:
            first_letter = seq1[new_i]
        if j == new_j:
            second_letter = "-"
        else:
            second_letter = seq2[new_j]
        cell = backward_matrix[new_i][new_j]
        first_string += first_letter
        second_string += second_letter
        i = new_i
        j = new_j
    while i > 0:
        first_string += seq1[i - 1]
        second_string += "-"
        i -= 1
    while j > 0:
        first_string += "-"
        second_string += seq2[j - 1]
        j -= 1
    first_string = first_string[::-1]
    second_string = second_string[::-1]
    return [first_string, second_string]

def find_backward_alignment_profiles(pack_of_seqs1, pack_of_seqs2, backward_matrix):
    seq_1 = False
    seq_2 = False
    if isinstance(pack_of_seqs2, str):
        pack_of_seqs2 = [pack_of_seqs2]
    if isinstance(pack_of_seqs1, str):
        pack_of_seqs1 = [pack_of_seqs1]
    first_height = len(pack_of_seqs1)

    second_height = len(pack_of_seqs2)

    i = len(pack_of_seqs1[0])
    j = len(pack_of_seqs2[0])

    cell = backward_matrix[i][j]

    first_string = ""
    second_string = ""


    while len(cell) > 2:
        if (cell[2] == "top"):
            new_i = i - 1
            new_j = j
            first_letter = "L"
            second_letter = "G"
        elif cell[2] == "left":
            new_i = i
            new_j = j - 1
            first_letter = "G"
            second_letter = "L"
        elif cell[2] == "diag":
            new_i = i - 1
            new_j = j - 1
            first_letter = "L"
            second_letter = "L"

        first_string += first_letter
        second_string += second_letter
        cell = backward_matrix[new_i][new_j]
        i = new_i
        j = new_j

    new_pack_of_seqs1 = ["" for i in range((first_height))]
    new_pack_of_seqs2 = ["" for i in range((second_height))]
    counter = 0

    for j in range(len(first_string) - 1, -1, -1):
            if first_string[j] == "G":
                for i in range(first_height):
                    new_pack_of_seqs1[i] += "-"
            else:
                for i in range(first_height):
                    new_pack_of_seqs1[i] += pack_of_seqs1[i][counter]
                counter += 1

    counter = 0
    for j in range(len(second_string) - 1, -1 , -1):
            if second_string[j] == "G":
                for i in range(second_height):
                    new_pack_of_seqs2[i] += "-"
            else:
                for i in range(second_height):
                    new_pack_of_seqs2[i] += pack_of_seqs2[i][counter]
                counter += 1
    return new_pack_of_seqs1 + new_pack_of_seqs2


def align_profile_to_sequence(profile1, sequence, blosum, gap_costs):
    alphabet = list(profile1.iterkeys())
    profile2 = create_profile(sequence, alphabet)
    return align_profile_to_profile(profile1, profile2, blosum, gap_costs)


def align_seq_to_seq(seq1, seq2, blosum, gap_costs):
    d = gap_costs[0]
    table_for_algo = align_table.align_table(seq1, seq2, 1)
    M = table_for_algo.table


    for i in xrange(1, len(seq1) + 1):
        M[i][0] = d * (i )
        table_for_algo.set_pred_top(i, 0)
    for j in xrange(1, len(seq2) + 1):
        M[0][j] = d * (j )
        table_for_algo.set_pred_left(0, j)

    for i in xrange(1, len(seq1) + 1):
        for j in xrange(1, len(seq2) + 1):
            M[i][j] = max(M[i - 1][j] + d, M[i][j - 1] + d, M[i - 1][j - 1] + blosum.wt(seq1[i - 1], seq2[j - 1]))
            if M[i][j] == M[i - 1][j] + d:
                table_for_algo.set_pred_top(i, j)
            elif M[i][j] == M[i][j - 1] + d:
                table_for_algo.set_pred_left(i, j)
            else:
                table_for_algo.set_pred_diag(i, j)
    return table_for_algo.predacessors, M[len(seq1)][len(seq2)]

def psp(profile1, profile2, i, j, blosum, gap_cost, first_height, second_height):
    alphabet = list(blosum.matrix.iterkeys())
    alphabet.pop(0)
    alphabet.append('-')
    psp = 0
    if profile1 == "-":
        return second_height * gap_cost
    if profile2 == "-":
        return first_height * gap_cost

    for x in alphabet:
        for y in alphabet:
            weight = gap_cost
            if (x != "-" and y != "_"):
                weight = blosum.wt(x, y)
            psp += (profile1[x][i] * profile2[y][j] * weight)
    return psp



def align_profile_to_profile(profile1, profile2, blosum, gap_costs):
    d = gap_costs[0]

    first_height = 0
    second_height = 0
    alphabet = list(blosum.matrix.iterkeys())
    alphabet.append("-")
    alphabet.remove("*")
    for x in alphabet:
        first_height += profile1[x][0]
        second_height += profile2[x][0]
    first_len = len(profile1['A'])
    second_len = len(profile2['A'])

    table_for_algo = align_table.align_table(first_len, second_len, 0)
    M = table_for_algo.table

    for i in xrange(1, first_len + 1):
        M[i][0] = i * psp(profile1, "-", i, 0, blosum, d, first_height, second_height)
        table_for_algo.set_pred_top(i , 0)
    for j in xrange(1, second_len + 1):
        M[0][j] = j * psp("-", profile2, 0, j, blosum, d, first_height, second_height)
        table_for_algo.set_pred_left(0, j)

    for i in xrange(1, first_len + 1):
        for j in xrange(1, second_len + 1):
            M[i][j] = max(
                M[i-1][j-1] + psp(profile1, profile2, i - 1, j - 1, blosum, d, first_height, second_height),
                M[i - 1][j] + psp(profile1, "-", i - 1 , j, blosum, d, first_height, second_height),
                M[i][j - 1] + psp("-", profile2, i, j - 1, blosum, d, first_height, second_height)
            )
            if M[i][j] == M[i][j - 1] + psp("-", profile2, i, j - 1, blosum, d, first_height, second_height):
                table_for_algo.set_pred_left(i, j)
            elif M[i][j] == M[i - 1][j] + psp(profile1, "-", i - 1, j, blosum, d, first_height, second_height):
                table_for_algo.set_pred_top(i, j)
            else:
                table_for_algo.set_pred_diag(i, j)

    return table_for_algo.predacessors, M[first_len][second_len]
