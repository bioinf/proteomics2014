#ifndef Cluster_h__
#define Cluster_h__

#include <string>

class Cluster
{
public:
	size_t left;
	size_t right;
	double distance;
	size_t cluster_size;
	std::string cluster_id;
	std::string sequence;

	Cluster(size_t left, size_t right, double distance, size_t cluster_size, std::string const& cluster_id):
		left(left), right(right), distance(distance), cluster_size(cluster_size), cluster_id(cluster_id)
	{}

	Cluster(std::string const& cluster_id, std::string const& sequence):
		left(0), right(0), distance(0), cluster_size(1), cluster_id(cluster_id), sequence(sequence)
	{}

	Cluster():
		left(0), right(0), distance(0), cluster_size(0)
	{}
};

#endif // Cluster_h__