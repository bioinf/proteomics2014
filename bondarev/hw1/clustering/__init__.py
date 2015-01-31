from functools import partial

from clustering.base import template_clustering
import clustering.neighbor_joining as nj
import clustering.xpgma

upgma_cluster = partial(template_clustering, xpgma.get_middle_score, xpgma.upgma_calc_score)
wpgma_cluster = partial(template_clustering, xpgma.get_middle_score, xpgma.wpgma_calc_score)
neigbor_joining_cluster = partial(template_clustering, nj.get_middle_score, nj.calc_score)
