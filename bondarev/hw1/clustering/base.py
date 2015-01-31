import logging

from alignment import get_global_alignment_score
import profile_alignment

logger = logging.getLogger()


def template_clustering(middle_score_fn, calc_score_fn, seqs, score_matrix):
    """
    :param middle_score_fn: cluster score if it is different from sequence score
     arguments: cluster_set, pair_scores and i, j: pair of cluster id for calc cluster score

    :param calc_score_fn: pair score function for cluster (u, v) with cluster w.
     arguments: score_matrix, u_id, v_id, w, _

    :param seqs: collection for Sequence
    :param score_matrix: ScoreMatrix object
    :return: collection of progressive aligned Sequences
    """
    seqs = list(seqs)
    logger.info('Pair align stage')
    pair_scores = _get_pair_align_score_matrix(seqs, score_matrix)
    logger.info('Finish pair align')

    logger.info('Clustering stage')
    logger.info('Make cluster set')
    cluster_set = ClusterSet(seqs)
    logger.info('Join clusters')
    for _ in range(len(seqs) - 1):   # count of operations for join all clusters
        logger.info('Find max score clusters')
        clust1, clust2 = cluster_set.get_max_score_clusters_pair(middle_score_fn, pair_scores)

        logger.info('Start join clusters {} and {}'.format(clust1.id, clust2.id))
        main_cluster_id, child_cluster_id = cluster_set.merge(clust1, clust2, score_matrix)

        logger.info('Update score matrix')
        for other in range(len(seqs)):
            score_i2k = calc_score_fn(pair_scores, main_cluster_id, child_cluster_id, other, cluster_set.sizes)
            pair_scores[main_cluster_id][other] = pair_scores[other][main_cluster_id] = score_i2k
        cluster_set.existing.remove(child_cluster_id)
        cluster_set.sizes[main_cluster_id] += cluster_set.sizes[child_cluster_id]

    return cluster_set.get_cluster_root(0).seqs


class Cluster:
    def __init__(self, id_, seq):
        self.id = id_
        self.seqs = [seq]
        self.depth = 1
        self.parent = id_

    def is_root(self):
        return self.id == self.parent


class ClusterSet:
    def __init__(self, seqs):
        n = len(seqs)
        self.clusters = [Cluster(i, seqs[i]) for i in range(n)]
        self.sizes = [1] * n
        self.existing = set(range(n))

    def __getitem__(self, i):
        return self.clusters[i]

    def get_cluster_root(self, id_):
        if self[id_].parent != id_:
            self[id_].parent = self.get_cluster_root(self[id_].parent).id
        return self.clusters[self[id_].parent]

    def merge(self, clust1, clust2, score_matrix):
        """
        :return: main cluster id and child cluster id.
            Second cluster included into main.
        """
        main, child = clust1, clust2
        if main.id == child.id:
            logger.warning('Equals cluster id {}'.format(child.id))
            return -1, -1

        if main.depth == child.depth:
            main.depth += 1

        if child.depth > main.depth:
            main, child = child, main

        # now main.depth > child.depth => child join into main
        child.parent = main.id
        logger.info('Run profile alignment')
        main.seqs = profile_alignment.align(main.seqs, child.seqs, score_matrix)
        logger.info('Finish profile alignment')

        self.clusters[main.id] = main
        self.clusters[child.id] = child
        return main.id, child.id

    def get_max_score_clusters_pair(self, middle_score_fn, pair_scores):
        n = len(self.clusters)
        max_i, max_j, max_value = -1, -1, float('-Inf')
        for i in range(n):
            if i not in self.existing:
                continue

            for j in range(i + 1, n):
                if j not in self.existing:
                    continue

                cluster_pair_score = middle_score_fn(self, pair_scores, i, j)
                if cluster_pair_score > max_value:
                    max_i, max_j, max_value = i, j, cluster_pair_score

        return self.get_cluster_root(max_i), self.get_cluster_root(max_j)


def _get_pair_align_score_matrix(seqs, score_matrix):
    n = len(seqs)
    scores = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            scores[i][j] = scores[j][i] = get_global_alignment_score(seqs[i].seq, seqs[j].seq, score_matrix)
        scores[i][i] = float('-Inf')
    return scores


