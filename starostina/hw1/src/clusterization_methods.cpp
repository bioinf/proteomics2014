#include <numeric>
#include <algorithm>

#include "clusterization_methods.hpp"

Clusterization::Clusterization(const DistanceMap &distances, size_t clusters_num) :
        clusters_num_(clusters_num) {
    CreateDistances(distances);
}

NodePtr Clusterization::Clusterize() {
    while (nodes_.size() != 1) {
        std::pair <NodePtr, NodePtr> closest = GetClosestClusters();
        MergeAndUpdateDistances(closest.first, closest.second);
    }
    return *(nodes_.begin());
}

void Clusterization::CreateDistances(const DistanceMap &distances)
{
    std::vector <NodePtr> nodes;
    for (size_t i = 0; i != clusters_num_; ++i) {
        nodes.push_back(NodePtr(new Node(i)));
    }

    for (auto it = distances.begin(); it != distances.end(); ++it) {
        size_t i = it->first.first;
        size_t j = it->first.second;
        distances_[std::make_pair(nodes[i], nodes[j])] = it->second;
    }
    nodes_.insert(nodes.begin(), nodes.end());
}

NeighbourJoining::NeighbourJoining(const DistanceMap &distances, size_t clusters_num) :
        Clusterization(distances, clusters_num) {
    BuildQ();
}

std::pair<NodePtr, NodePtr> NeighbourJoining::GetClosestClusters() {
    auto closest_it = Q_.begin();
    for (auto it = Q_.begin(); it != Q_.end(); ++it) {
        if (it->second > closest_it->second) {
            closest_it = it;
        }
    }
    return closest_it->first;
}

void NeighbourJoining::MergeAndUpdateDistances(NodePtr i, NodePtr j) {
    NodePtr new_node(new Node(i, j));
    nodes_.erase(i);
    nodes_.erase(j);
    auto dist_ij_it = distances_.find(std::make_pair(i, j));
    if (dist_ij_it == distances_.end()) {
        dist_ij_it = distances_.find(std::make_pair(j, i));
    }
    for (auto node_it = nodes_.begin(); node_it != nodes_.end(); ++node_it) {
        auto dist_i_it = distances_.find(std::make_pair(i, *node_it));
        if (dist_i_it == distances_.end()) {
            dist_i_it = distances_.find(std::make_pair(*node_it, i));
        }
        auto dist_j_it = distances_.find(std::make_pair(j, *node_it));
        if (dist_j_it == distances_.end()) {
            dist_j_it = distances_.find(std::make_pair(*node_it, j));
        }
        distances_.erase(dist_i_it);
        distances_.erase(dist_j_it);
        distances_[std::make_pair(*node_it, new_node)] = (dist_i_it->second + dist_j_it->second - dist_ij_it->second) / 2;
    }
    distances_.erase(dist_ij_it);
    nodes_.insert(new_node);
    BuildQ();
}

void NeighbourJoining::BuildQ() {
    Q_.clear();
    for (auto it = distances_.begin(); it != distances_.end(); ++it) {
        NodePtr i = it->first.first;
        NodePtr j = it->first.second;
        double sum_i = 0;
        double sum_j = 0;
        for (auto node_it = nodes_.begin(); node_it != nodes_.end(); ++node_it) {
            if (*node_it == i || *node_it == j) {
                continue;
            }
            auto dist_i_it = distances_.find(std::make_pair(i, *node_it));
            if (dist_i_it == distances_.end()) {
                dist_i_it = distances_.find(std::make_pair(*node_it, i));
            }
            auto dist_j_it = distances_.find(std::make_pair(j, *node_it));
            if (dist_j_it == distances_.end()) {
                dist_j_it = distances_.find(std::make_pair(*node_it, j));
            }
            sum_i += dist_i_it->second;
            sum_j += dist_j_it->second;
        }
        Q_[std::make_pair(i, j)] = (clusters_num_ - 2) * it->second - sum_i - sum_j;
    }
}

UPGMA::UPGMA(const DistanceMap &distances, size_t clusters_num) :
        Clusterization(distances, clusters_num) {
    InitSizes();
}

std::pair<NodePtr, NodePtr> UPGMA::GetClosestClusters() {
    auto closest_it = distances_.begin();
    for (auto it = distances_.begin(); it != distances_.end(); ++it) {
        if (it->second > closest_it->second) {
            closest_it = it;
        }
    }
    return closest_it->first;
}

void UPGMA::MergeAndUpdateDistances(NodePtr i, NodePtr j) {
    NodePtr new_node(new Node(i, j));
    size_t size_i = sizes_[i];
    size_t size_j = sizes_[j];
    nodes_.erase(i);
    nodes_.erase(j);
    sizes_.erase(i);
    sizes_.erase(j);
    for (auto node_it = nodes_.begin(); node_it != nodes_.end(); ++node_it) {
        auto dist_i_it = distances_.find(std::make_pair(i, *node_it));
        if (dist_i_it == distances_.end()) {
            dist_i_it = distances_.find(std::make_pair(*node_it, i));
        }
        auto dist_j_it = distances_.find(std::make_pair(j, *node_it));
        if (dist_j_it == distances_.end()) {
            dist_j_it = distances_.find(std::make_pair(*node_it, j));
        }
        distances_[std::make_pair(*node_it, new_node)] = (size_i * (dist_i_it->second) + size_j * (dist_j_it->second)) / (size_i + size_j);
        distances_.erase(dist_i_it);
        distances_.erase(dist_j_it);
    }

    auto dist_ij_it = distances_.find(std::make_pair(i, j));
    if (dist_ij_it == distances_.end()) {
        dist_ij_it = distances_.find(std::make_pair(j, i));
    }
    distances_.erase(dist_ij_it);
    nodes_.insert(new_node);
    sizes_[new_node] = size_i + size_j;
}

void UPGMA::InitSizes() {
    for (auto node_it = nodes_.begin(); node_it != nodes_.end(); ++node_it) {
        sizes_[*node_it] = 1;
    }
}

std::pair<NodePtr, NodePtr> WPGMA::GetClosestClusters() {
    auto closest_it = distances_.begin();
    for (auto it = distances_.begin(); it != distances_.end(); ++it) {
        if (it->second > closest_it->second) {
            closest_it = it;
        }
    }
    return closest_it->first;
}

void WPGMA::MergeAndUpdateDistances(NodePtr i, NodePtr j) {
    NodePtr new_node(new Node(i, j));
    nodes_.erase(i);
    nodes_.erase(j);
    for (auto node_it = nodes_.begin(); node_it != nodes_.end(); ++node_it) {
        auto dist_i_it = distances_.find(std::make_pair(i, *node_it));
        if (dist_i_it == distances_.end()) {
            dist_i_it = distances_.find(std::make_pair(*node_it, i));
        }
        auto dist_j_it = distances_.find(std::make_pair(j, *node_it));
        if (dist_j_it == distances_.end()) {
            dist_j_it = distances_.find(std::make_pair(*node_it, j));
        }
        distances_[std::make_pair(*node_it, new_node)] = ((dist_i_it->second) + (dist_j_it->second)) / 2;
        distances_.erase(dist_i_it);
        distances_.erase(dist_j_it);
    }

    auto dist_ij_it = distances_.find(std::make_pair(i, j));
    if (dist_ij_it == distances_.end()) {
        dist_ij_it = distances_.find(std::make_pair(j, i));
    }
    distances_.erase(dist_ij_it);
    nodes_.insert(new_node);
}
