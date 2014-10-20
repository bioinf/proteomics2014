__author__ = 'letovesnoi'

# get alignment for two sequences:
def tracebackSequences(str1, str2, back, v, w, i, j):
    if i == -1 or j == -1:
        if i == -1:
            while j > -1:
                str1 += '*'
                str2 += w[j]
                j -= 1
        if j == -1:
            while i > -1:
                str1 += v[i]
                str2 += '*'
                i -= 1
        return str1, str2

    if back[i, j] == (i, j - 1):
        str1, str2 = tracebackSequences(str1, str2, back, v, w, i, j - 1)
        str1 += '*'
        str2 += w[j]
    elif back[i, j] == (i - 1, j):
        str1, str2 = tracebackSequences(str1, str2, back, v, w, i - 1, j)
        str1 += v[i]
        str2 += '*'
    else:
        str1, str2 = tracebackSequences(str1, str2, back, v, w, i - 1, j - 1)
        str1 += v[i]
        str2 += w[j]
    return str1, str2

# get score for two sequences:
def getDistanceSequences(v, w, matrix):
    s = {}
    back = {}
    s[-1, -1] = 0
    for i in range(0, len(v)):
        s[i, -1] = s[i - 1, -1] + matrix[v[i], '*']
    for j in range(0, len(w)):
        s[-1, j] = s[-1, j - 1] + matrix['*', w[j]]
    for j in range(0, len(w)):
        for i in range(0, len(v)):
            s[i, j] = max(s[i - 1, j] + matrix[v[i], '*'],
                          s[i, j - 1] + matrix['*', w[j]],
                          s[i - 1, j - 1] + matrix[v[i], w[j]])
            if s[i, j] == s[i - 1, j - 1] + matrix[v[i], w[j]]:
                back[i, j] = (i - 1, j - 1)
            elif s[i, j] == s[i, j - 1] + matrix['*', w[j]]:
                back[i, j] = (i, j - 1)
            elif s[i, j] == s[i - 1, j] + matrix[v[i], '*']:
                back[i, j] = (i - 1, j)
    #str1, str2 = tracebackSequences(str1, str2, back, v, w, len(v) - 1, len(w) - 1)
    return s[len(v) - 1, len(w) - 1]

# get score and alignment two sequences:
def getAlignmentSequences(str1, str2, v, w, matrix):
    s = {}
    back = {}
    s[-1, -1] = 0
    for i in range(0, len(v)):
        s[i, -1] = s[i - 1, -1] + matrix[v[i], '*']
    for j in range(0, len(w)):
        s[-1, j] = s[-1, j - 1] + matrix['*', w[j]]
    for j in range(0, len(w)):
        for i in range(0, len(v)):
            s[i, j] = max(s[i - 1, j] + matrix[v[i], '*'],
                          s[i, j - 1] + matrix['*', w[j]],
                          s[i - 1, j - 1] + matrix[v[i], w[j]])
            if s[i, j] == s[i - 1, j - 1] + matrix[v[i], w[j]]:
                back[i, j] = (i - 1, j - 1)
            elif s[i, j] == s[i, j - 1] + matrix['*', w[j]]:
                back[i, j] = (i, j - 1)
            elif s[i, j] == s[i - 1, j] + matrix[v[i], '*']:
                back[i, j] = (i - 1, j)
    str1, str2 = tracebackSequences(str1, str2, back, v, w, len(v) - 1, len(w) - 1)
    return s[len(v) - 1, len(w) - 1], [str1, str2]

# get mismatch score:
def getScoreIJProfile(matrix, profile1, profile2, i1, i2):
    score = 0
    for acid1 in profile1:
        for acid2 in profile2:
            score += profile1[acid1][i1] * profile2[acid2][i2] * matrix[acid1, acid2]
    return score

# get gap score:
def getScoreGAPProfile(matrix, profile, iProfile, gapProfile):
    score = 0
    for acid in profile:
        score += profile[acid][iProfile] * matrix[acid, '*']
    return score

def getTracebackProfiles(matrix, profile1, profile2):
    traceback = {}
    s = {}
    s[-1, -1] = 0
    for i in range(0, len(profile1['*'])):
        s[i, -1] = s[i - 1, -1] + getScoreGAPProfile(matrix, profile1, i, profile2)
    for j in range(0, len(profile2['*'])):
        s[-1, j] = s[-1, j - 1] + getScoreGAPProfile(matrix, profile2, j, profile1)
    for i in range(0, len(profile1['*'])):
        for j in range(0, len(profile2['*'])):
            valueI1J = s[i - 1, j] + getScoreGAPProfile(matrix, profile1, i, profile2)
            valueIJ1 = s[i, j - 1] + getScoreGAPProfile(matrix, profile2, j, profile1)
            valueI1J1 = s[i - 1, j - 1] + getScoreIJProfile(matrix, profile1, profile2, i, j)
            s[i, j] = max(valueI1J, valueIJ1, valueI1J1)
            if s[i, j] == valueI1J1:
                traceback[i, j] = (i - 1, j - 1)
            elif s[i, j] == valueIJ1:
                traceback[i, j] = (i, j - 1)
            elif s[i, j] == valueI1J:
                traceback[i, j] = (i - 1, j)
    score = s[len(profile1['*']) - 1, len(profile2['*']) - 1]
    return traceback

# get alignment for two profiles:
def getGlueAlignment(back, tmpAlignment, alignment1, alignment2, i, j):
    if i == -1 or j == -1:
        if i == -1:
            while j > -1:
                for iAlignment in range(len(alignment1)):
                    tmpAlignment[iAlignment] += '*'
                for iAlignment in range(len(alignment2)):
                    tmpAlignment[len(alignment1) + iAlignment] += alignment2[iAlignment][j]
                j -= 1
        if j == -1:
            while i > -1:
                for iAlignment in range(len(alignment1)):
                    tmpAlignment[iAlignment] += alignment1[iAlignment][i]
                for iAlignment in range(len(alignment2)):
                    tmpAlignment[len(alignment1) + iAlignment] += '*'
                i -= 1
        return tmpAlignment

    if back[i, j] == (i, j - 1):
        tmpAlignment = getGlueAlignment(back, tmpAlignment, alignment1, alignment2, i, j - 1)
        for iAlignment in range(len(alignment1)):
            tmpAlignment[iAlignment] += '*'
        for iAlignment in range(len(alignment2)):
            tmpAlignment[len(alignment1) + iAlignment] += alignment2[iAlignment][j]
    elif back[i, j] == (i - 1, j):
        tmpAlignment = getGlueAlignment(back, tmpAlignment, alignment1, alignment2, i - 1, j)
        for iAlignment in range(len(alignment1)):
            tmpAlignment[iAlignment] += alignment1[iAlignment][i]
        for iAlignment in range(len(alignment2)):
            tmpAlignment[len(alignment1) + iAlignment] += '*'
    else:
        tmpAlignment = getGlueAlignment(back, tmpAlignment, alignment1, alignment2, i - 1, j - 1)
        for iAlignment in range(len(alignment1)):
            tmpAlignment[iAlignment] += alignment1[iAlignment][i]
        for iAlignment in range(len(alignment2)):
            tmpAlignment[len(alignment1) + iAlignment] += alignment2[iAlignment][j]

    return tmpAlignment


