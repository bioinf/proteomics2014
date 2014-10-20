__author__ = 'letovesnoi'

import Clustering

def getMaxSimilarityForWPGMA(d):
    maxSimilarity = Clustering.getMax(d)
    return maxSimilarity

def getDeltaDistForWPGMA(d, f, g, u):
    delta = {}
    delta[u, f] = 0.5 * d[f, g]
    delta[u, g] = d[f, g] - delta[u, f]
    return delta

def getDistUForWPGMA(d, newD, f, g, u, seq):
    for k in seq:
        newD[u, k] = 0.5 * (d[f, k] + d[g, k])
        newD[k, u] = newD[u, k]
    return newD