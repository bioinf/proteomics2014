#pragma once

#include <memory>
#include <vector>
#include <set>
#include <map>

class Node;
typedef std::shared_ptr <Node> NodePtr;

struct Node {
    Node(NodePtr child1, NodePtr child2) :
        child1_(child1), child2_(child2) { }

    Node(size_t seq_ind) :
        child1_(nullptr), child2_(nullptr), seq_ind_(seq_ind) { }

    NodePtr child1_;
    NodePtr child2_;
    size_t seq_ind_;
};

typedef std::map <std::pair <size_t, size_t>, double> DistanceMap;
typedef std::map <std::pair <NodePtr, NodePtr>, double> ScoreMap;

class Clusterization {
public:
    Clusterization(DistanceMap const &distances, size_t clusters_num);
    virtual ~Clusterization() { }

    virtual std::pair <NodePtr, NodePtr> GetClosestClusters() = 0;
    virtual void MergeAndUpdateDistances(NodePtr i1, NodePtr i2) = 0;

    NodePtr Clusterize();

protected:
    void CreateDistances(DistanceMap const &distances);

    size_t clusters_num_;
    std::set <NodePtr> nodes_;
    ScoreMap distances_;
};

class NeighbourJoining : public Clusterization {
public:
    NeighbourJoining(DistanceMap const &distances, size_t clusters_num);

    virtual std::pair <NodePtr, NodePtr> GetClosestClusters();
    virtual void MergeAndUpdateDistances(NodePtr i, NodePtr j);

private:
    void BuildQ();

    ScoreMap Q_;
};

class UPGMA : public Clusterization {
public:
    UPGMA(DistanceMap const &distances, size_t clusters_num);

    virtual std::pair <NodePtr, NodePtr> GetClosestClusters();
    virtual void MergeAndUpdateDistances(NodePtr i, NodePtr j);

private:
    void InitSizes();

    std::map <NodePtr, size_t> sizes_;
};

class WPGMA : public Clusterization {
public:
    WPGMA(DistanceMap const &distances, size_t clusters_num) :
        Clusterization(distances, clusters_num) { }

    virtual std::pair <NodePtr, NodePtr> GetClosestClusters();
    virtual void MergeAndUpdateDistances(NodePtr i, NodePtr j);
};
