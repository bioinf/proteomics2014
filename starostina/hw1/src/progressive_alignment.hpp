#include "alignment.hpp"
#include "clusterization_methods.hpp"
#include "sequence.hpp"

enum ClusterizationType {
    NJ_t, UPGMA_t, WPGMA_t
};

class ProgressiveAlignment {
public:
    ProgressiveAlignment(ClusterizationType type, Matrix const &m, std::vector <Sequence> &sequences) :
        sequences_(sequences), type_(type), matrix_(m), guide_tree_root_(nullptr) { }

    std::shared_ptr <DistanceMap> build_distances() const;
    void build_guide_tree();
    std::shared_ptr<PSSM> align_node(NodePtr node);
    void align();

private:
    std::vector <Sequence> &sequences_;
    ClusterizationType type_;
    Matrix const &matrix_;
    NodePtr guide_tree_root_;
};
