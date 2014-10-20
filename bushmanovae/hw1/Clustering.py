__author__ = 'letovesnoi'

import AlignmentsUtils
import NJUtils
import UPGMAUtils
import WPGMAUtils

def constructEDM(seq, matrix):
    EDM = {}
    for i in range(len(seq)):
        for j in range(len(seq)):
            if i != j:
                EDM[seq[i], seq[j]] = AlignmentsUtils.getDistanceSequences(seq[i], seq[j], matrix)
    return EDM

def getMax(EDM):
    max = {}
    max['value'] = EDM[EDM.keys()[0]]
    max['key'] = EDM.keys()[0]
    for pair in EDM:
        if max['value'] < EDM[pair]:
            max['value'] = EDM[pair]
            max['key'] = pair
    return max

def updateEDMSTI(d, seq, tree, num, clusteringType):
    newD = {}

    if clusteringType == 'NJ':
        maxSimilarity = NJUtils.getMaxSimilarityForNJ(d, seq)
    elif clusteringType == 'UPGMA':
        maxSimilarity = UPGMAUtils.getMaxSimilarityForUPGMA(d)
    elif clusteringType == 'WPGMA':
        maxSimilarity = WPGMAUtils.getMaxSimilarityForWPGMA(d)

    f = maxSimilarity['key'][0]
    iF = seq.index(f)
    numF = num[iF]

    g = maxSimilarity['key'][1]
    if f == g:
        iG = iF + 1 + seq[iF + 1:].index(g)
    else:
        iG = seq.index(g)
    numG = num[iG]

    u = '(' + f + ':' + str(d[maxSimilarity['key']]) + ':' + g + ')'

    seq.remove(f)
    seq.remove(g)

    if clusteringType == 'NJ':
        delta = NJUtils.getDeltaDistForNJ(d, f, g, u, seq)
    elif clusteringType == 'UPGMA':
        delta = UPGMAUtils.getDeltaDistForUPGMA(d, f, g, u)
    elif clusteringType == 'WPGMA':
        delta = WPGMAUtils.getDeltaDistForWPGMA(d, f, g, u)

    tree.append({'1': u, '2': f, 'dist': delta[u, f]})
    tree.append({'1': u, '2': g, 'dist': delta[u, g]})

    for i in range(len(seq)):
        for j in range(len(seq)):
            if i != j:
                newD[seq[i], seq[j]] = d[seq[i], seq[j]]

    if clusteringType == 'NJ':
        newD = NJUtils.getDistUForNJ(d, newD, f, g, u, seq)
    elif clusteringType == 'UPGMA':
        newD = UPGMAUtils.getDistUForUPGMA(d, newD, f, g, u, seq, numF, numG)
    elif clusteringType == 'WPGMA':
        newD = WPGMAUtils.getDistUForWPGMA(d, newD, f, g, u, seq)

    seq.append(u)
    num.append(numF + numG + 1)

    return newD, seq, num, tree, iF, iG