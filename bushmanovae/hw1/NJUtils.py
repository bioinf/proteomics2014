__author__ = 'letovesnoi'

import Clustering

def getMaxSimilarityForNJ(d, seq):
    Q = {}
    n = len(seq)
    for (i, j) in d:
        sumI = 0
        sumJ = 0
        for k in seq:
            if i != k:
                sumI += d[i, k]
            if j != k:
                sumJ += d[j, k]
        Q[i, j] = (n - 2) * d[i, j] - sumI - sumJ
    maxSimilarity = Clustering.getMax(Q)
    return maxSimilarity

def getDeltaDistForNJ(d, f, g, u, seq):
    delta = {}
    n = len(seq)
    sumF = 0
    sumG = 0
    for k in seq:
        if f != k:
            sumF += d[f, k]
        if g != k:
            sumG += d[g, k]

    delta[u, f] = 0.5 * d[f, g]
    if n != 2:
        delta[u, f] += 1.0 / (2 * (n - 2)) * (sumF - sumG)
    delta[u, g] = d[f, g] - delta[u, f]
    return delta

def getDistUForNJ(d, newD, f, g, u, seq):
    for k in seq:
        newD[u, k] = 0.5 * (d[f, k] + d[g, k] - d[f, g])
        newD[k, u] = newD[u, k]
    return newD