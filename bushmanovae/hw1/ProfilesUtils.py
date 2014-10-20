__author__ = 'letovesnoi'

# get tmpProfile for two sequences:
def getProfileSequences(alignments, acids):
    profile = {}
    for acid in acids:
        profile[acid] = []
        for i in range(len(alignments[0])):
            profile[acid].append(0.0)

    for i in range(len(alignments[0])):
        for alignment in alignments:
            acid = alignment[i]
            profile[acid][i] += 1
        for acid in acids:
            profile[acid][i] /= len(alignments)

    return profile

# get tmpProfile for one sequence:
def getProfileSequence(seq, acids):
    profile = {}
    for acid in acids:
        profile[acid] = []
        for i in range(len(seq)):
            profile[acid].append(0.0)

    for i in range(len(seq)):
        acid = seq[i]
        profile[acid][i] += 1

    return profile

# get tmpProfile for two profiles:
def getGlueProfile(acids, tmpProfile, profile1, profile2, number1, number2, back, i, j):
    if i == -1 or j == -1:
        if i == -1:
            while j > -1:
                for acid in acids:
                    tmpProfile[acid].append(profile2[acid][j] * number2)
                tmpProfile['*'][-1] += number1
                for acid in acids:
                    tmpProfile[acid][-1] /= (number1 + number2)
                j -= 1
        if j == -1:
            while i > -1:
                for acid in acids:
                    tmpProfile[acid].append(profile1[acid][i] * number1)
                tmpProfile['*'][-1] += number2
                for acid in acids:
                    tmpProfile[acid][-1] /= (number1 + number2)
                i -= 1
        return tmpProfile

    if back[i, j] == (i, j - 1):
        tmpProfile = getGlueProfile(acids, tmpProfile, profile1, profile2, number1, number2, back, i, j - 1)
        for acid in acids:
            tmpProfile[acid].append(profile2[acid][j] * number2)
        tmpProfile['*'][-1] += number1
        for acid in acids:
            tmpProfile[acid][-1] /= (number1 + number2)
    elif back[i, j] == (i - 1, j):
        tmpProfile = getGlueProfile(acids, tmpProfile, profile1, profile2, number1, number2, back, i - 1, j)
        for acid in acids:
            tmpProfile[acid].append(profile1[acid][i] * number1)
        tmpProfile['*'][-1] += number2
        for acid in acids:
            tmpProfile[acid][-1] /= (number1 + number2)
    else:
        tmpProfile = getGlueProfile(acids, tmpProfile, profile1, profile2, number1, number2, back, i - 1, j - 1)
        for acid in acids:
            tmpProfile[acid].append((profile1[acid][i] * number1 + profile2[acid][j] * number2) / (number1 + number2))

    return tmpProfile