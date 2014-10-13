#include "progressive_alignment.hpp"

std::shared_ptr<DistanceMap> ProgressiveAlignment::build_distances() const {
    std::shared_ptr <DistanceMap> distances_ptr(new DistanceMap());
    for (size_t i = 0; i != sequences_.size(); ++i) {
        for (size_t j = i + 1; j != sequences_.size(); ++j) {
            PSSM s1(sequences_[i].sequence, i);
            PSSM s2(sequences_[j].sequence, j);
            AlignmentMatrix aln_matrix;
            build_alignment_matrix(s1, s2, matrix_, aln_matrix);
            (*distances_ptr)[std::make_pair(i, j)] = get_alignment_score(aln_matrix);
        }
    }
    return distances_ptr;
}

void ProgressiveAlignment::build_guide_tree() {
    auto distances_ptr = build_distances();
    std::shared_ptr <Clusterization> clusterer;
    switch(type_) {
    case NJ_t:
        clusterer = std::shared_ptr <Clusterization> (new NeighbourJoining(*distances_ptr, sequences_.size()));
        break;
    case UPGMA_t:
        clusterer = std::shared_ptr <Clusterization> (new UPGMA(*distances_ptr, sequences_.size()));
        break;
    case WPGMA_t:
        clusterer = std::shared_ptr <Clusterization> (new WPGMA(*distances_ptr, sequences_.size()));
        break;
    }

    guide_tree_root_ = clusterer->Clusterize();
}

std::shared_ptr <PSSM> ProgressiveAlignment::align_node(NodePtr node) {
    if (node->child1_) {
        auto pssm1_ptr = align_node(node->child1_);
        auto pssm2_ptr = align_node(node->child2_);
        AlignmentMatrix aln_matrix;
        build_alignment_matrix(*pssm1_ptr, *pssm2_ptr, matrix_, aln_matrix);
        pssm1_ptr->merge_and_update_sequences(*pssm2_ptr, aln_matrix, sequences_);
        return pssm1_ptr;
    } else {
        return std::shared_ptr <PSSM>(new PSSM(sequences_[node->seq_ind_].sequence, node->seq_ind_));
    }
}

void ProgressiveAlignment::align() {
    build_guide_tree();
    align_node(guide_tree_root_);
}
