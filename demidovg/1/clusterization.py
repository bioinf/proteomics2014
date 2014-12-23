__author__ = 'german'
from collections import defaultdict
from alignment import *
import substitution
import copy

filename = "BLOSUM62.txt"

    # consts
gap_open = -8
gap_extention = -1
gap_costs = (gap_open, gap_extention)

blosum = substitution.Substitution_matrix(filename)
alphabet = list(blosum.matrix.iterkeys())
alphabet.append("-")
alphabet.remove("*")

class Cluster:
    block_of_sequences = None
    profile = None
    num_of_objects = 0
    name = None
    def __init__(self, name, sequence):
        self.block_of_sequences = sequence
        self.name = name
        self.num_of_objects += 1

    def merge(self, cluster1):
        self.num_of_objects += cluster1.num_of_objects
        self.name = (self.name, cluster1.name)
        new_block = []
        self.block_of_sequences = self.alignment(cluster1)
        print self.name

        self.profile = create_profile(self.block_of_sequences, alphabet)
        return self

    def alignment(self, cluster1):
        if self.profile == None and cluster1.profile == None:
            table = align_seq_to_seq(self.block_of_sequences, cluster1.block_of_sequences, blosum, gap_costs)[0]
            return find_backward_alignment(self.block_of_sequences, cluster1.block_of_sequences, table)
        elif self.profile != None and cluster1.profile == None:
            table = align_profile_to_sequence(self.profile, cluster1.block_of_sequences, blosum, gap_costs)[0]
            return find_backward_alignment_profiles(self.block_of_sequences, cluster1.block_of_sequences, table)
        elif self.profile == None and cluster1.profile != None:
            table = align_profile_to_sequence(cluster1.profile, self.block_of_sequences, blosum, gap_costs)[0]
            return find_backward_alignment_profiles(cluster1.block_of_sequences, self.block_of_sequences, table)
        else:
            table = align_profile_to_profile(self.profile, cluster1.profile, blosum, gap_costs)[0]
            return find_backward_alignment_profiles(self.block_of_sequences, cluster1.block_of_sequences, table)


class Clusterizator:
    matrix_of_distances = None
    info_clusters = {}
    min_dist = 1000000000
    def __init__(self, dist_matr, sequences):
        self.matrix_of_distances = dist_matr
        self.info_clusters = {}
        for key in self.matrix_of_distances.iterkeys():
            self.info_clusters[key] = Cluster(key, sequences.seqs[key])

    def find_min(self, matr_of_dist):
        min_dist = self.min_dist
        pair = None
        alphabet = list(matr_of_dist.iterkeys())
        for i in range(0, len(alphabet)):
            for j in range(i + 1, len(alphabet)):
                if alphabet[j] in matr_of_dist[alphabet[i]]:
                    if matr_of_dist[alphabet[i]][alphabet[j]] < min_dist:
                        pair = (alphabet[i], alphabet[j])
                        min_dist = matr_of_dist[alphabet[i]][alphabet[j]]
        return pair

    def q_matrix_calc(self):
        Q = defaultdict(dict)
        n = len(self.matrix_of_distances)
        min_i = None
        min_j = None
        min_dist = 10000

        for i in self.matrix_of_distances.iterkeys():
            for j in self.matrix_of_distances.iterkeys():
                minus = 0
                for k in self.matrix_of_distances.iterkeys():
                    minus += self.matrix_of_distances[i][k]
                    minus += self.matrix_of_distances[j][k]
                Q[i][j] = ((n - 2) * self.matrix_of_distances[i][j] - minus)
                if (Q[i][j] < min_dist and not i == j):
                    min_i = i
                    min_j = j
                    min_dist = Q[i][j]
        self.matrix_Q = Q

        return min_i, min_j

    def nj(self):
        min_i, min_j = self.q_matrix_calc()
        pair_to_merge = (min_i, min_j)
        merged_cluster = self.info_clusters[pair_to_merge[0]].merge(self.info_clusters[pair_to_merge[1]])
        self.info_clusters[pair_to_merge] = merged_cluster

        small_keys = copy.deepcopy(self.matrix_Q)

        while (len(small_keys)) > 2:
            for key in list(self.matrix_Q.iterkeys()):
                if key in small_keys.iterkeys():
                    if key not in pair_to_merge:
                        self.matrix_Q[merged_cluster.name][key] =  (self.matrix_Q[pair_to_merge[0]][key] +
                        self.matrix_Q[pair_to_merge[1]][key] - self.matrix_Q[pair_to_merge[0]][pair_to_merge[1]]) / 2
                        self.matrix_Q[key][merged_cluster.name] = self.matrix_Q[merged_cluster.name][key]

                        self.matrix_Q[key].pop((pair_to_merge[1]))
                        self.matrix_Q[key].pop((pair_to_merge[0]))
            self.matrix_Q.pop(pair_to_merge[0])
            self.matrix_Q.pop(pair_to_merge[1])

            small_keys = copy.deepcopy(self.matrix_Q)

            self.matrix_Q[merged_cluster.name][merged_cluster.name] = 0


            self.info_clusters[pair_to_merge] = merged_cluster
            self.info_clusters.pop(pair_to_merge[0])
            self.info_clusters.pop(pair_to_merge[1])

            pair_to_merge = Clusterizator.find_min(self, self.matrix_Q)
            merged_cluster = self.info_clusters[pair_to_merge[0]].merge(self.info_clusters[pair_to_merge[1]])
            """self.info_clusters[pair_to_merge] = merged_cluster
            self.info_clusters.pop(pair_to_merge[0])
            self.info_clusters.pop(pair_to_merge[1])"""

        lst = list(self.info_clusters.itervalues())
        for elem in lst[0].block_of_sequences:
            print elem


    def upgma(self):
        pair_to_merge = Clusterizator.find_min(self, self.matrix_of_distances)
        merged_cluster = self.info_clusters[pair_to_merge[0]].merge(self.info_clusters[pair_to_merge[1]])
        self.info_clusters[pair_to_merge] = merged_cluster

        small_keys = copy.deepcopy(self.matrix_of_distances)

        while (len(small_keys)) > 2:
            for key in list(self.matrix_of_distances.iterkeys()):
                if key in small_keys.iterkeys():
                    if key not in pair_to_merge:
                        self.matrix_of_distances[merged_cluster.name][key] =  (self.matrix_of_distances[pair_to_merge[0]][key] +
                        self.matrix_of_distances[pair_to_merge[1]][key]) / 2
                        self.matrix_of_distances[key][merged_cluster.name] = self.matrix_of_distances[merged_cluster.name][key]

                        self.matrix_of_distances[key].pop((pair_to_merge[1]))
                        self.matrix_of_distances[key].pop((pair_to_merge[0]))
            self.matrix_of_distances.pop(pair_to_merge[0])
            self.matrix_of_distances.pop(pair_to_merge[1])

            small_keys = copy.deepcopy(self.matrix_of_distances)

            self.matrix_of_distances[merged_cluster.name][merged_cluster.name] = 0


            self.info_clusters[pair_to_merge] = merged_cluster
            self.info_clusters.pop(pair_to_merge[0])
            self.info_clusters.pop(pair_to_merge[1])

            pair_to_merge = Clusterizator.find_min(self, self.matrix_of_distances)
            merged_cluster = self.info_clusters[pair_to_merge[0]].merge(self.info_clusters[pair_to_merge[1]])
            """self.info_clusters[pair_to_merge] = merged_cluster
            self.info_clusters.pop(pair_to_merge[0])
            self.info_clusters.pop(pair_to_merge[1])"""

        lst = list(self.info_clusters.itervalues())
        for elem in lst[0].block_of_sequences:
            print elem

    def wpgma(self):
        pair_to_merge = Clusterizator.find_min(self, self.matrix_of_distances)
        merged_cluster = self.info_clusters[pair_to_merge[0]].merge(self.info_clusters[pair_to_merge[1]])
        self.info_clusters[pair_to_merge] = merged_cluster

        small_keys = copy.deepcopy(self.matrix_of_distances)

        while (len(small_keys)) > 2:
            for key in list(self.matrix_of_distances.iterkeys()):
                if key in small_keys.iterkeys():
                    if key not in pair_to_merge:

                        u = self.info_clusters[pair_to_merge[0]].num_of_objects
                        v = self.info_clusters[pair_to_merge[1]].num_of_objects
                        self.matrix_of_distances[merged_cluster.name][key] = ((u * self.matrix_of_distances[pair_to_merge[0]][key] +
                        v * self.matrix_of_distances[pair_to_merge[1]][key]) / (u + v))
                        self.matrix_of_distances[key][merged_cluster.name] = self.matrix_of_distances[merged_cluster.name][key]

                        self.matrix_of_distances[key].pop((pair_to_merge[1]))
                        self.matrix_of_distances[key].pop((pair_to_merge[0]))
            self.matrix_of_distances.pop(pair_to_merge[0])
            self.matrix_of_distances.pop(pair_to_merge[1])

            small_keys = copy.deepcopy(self.matrix_of_distances)

            self.matrix_of_distances[merged_cluster.name][merged_cluster.name] = 0


            self.info_clusters[pair_to_merge] = merged_cluster
            self.info_clusters.pop(pair_to_merge[0])
            self.info_clusters.pop(pair_to_merge[1])

            pair_to_merge = Clusterizator.find_min(self, self.matrix_of_distances)
            merged_cluster = self.info_clusters[pair_to_merge[0]].merge(self.info_clusters[pair_to_merge[1]])
            """self.info_clusters[pair_to_merge] = merged_cluster
            self.info_clusters.pop(pair_to_merge[0])
            self.info_clusters.pop(pair_to_merge[1])"""

        lst = list(self.info_clusters.itervalues())
        for elem in lst[0].block_of_sequences:
            print elem
