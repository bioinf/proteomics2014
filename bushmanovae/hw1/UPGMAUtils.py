__author__ = 'letovesnoi'

import Clustering

def getMaxSimilarityForUPGMA(d):
    maxSimilarity = Clustering.getMax(d)
    return maxSimilarity

def getDeltaDistForUPGMA(d, f, g, u):
    delta = {}
    delta[u, f] = 0.5 * d[f, g]
    delta[u, g] = d[f, g] - delta[u, f]
    return delta

def getDistUForUPGMA(d, newD, f, g, u, seq, numF, numG):
    for k in seq:
        newD[u, k] = (numF * d[f, k] + numG * d[g, k]) * 1.0 / (numF + numG)
        newD[k, u] = newD[u, k]
    return newD