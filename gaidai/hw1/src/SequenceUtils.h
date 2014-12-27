#ifndef SequenceUtils_h__
#define SequenceUtils_h__

#include "Cluster.h"

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include <list>
#include <iterator>
#include <cmath>

typedef std::vector<std::map<char, int>> Profile;

// first - multiple alignment (by columns), second - profile (by columns)
typedef std::pair < std::vector<std::string>,  Profile> AlignmentProfile;

class SequenceUtils
{
public:
	SequenceUtils();
	~SequenceUtils();

	enum ClusterizationMethod {UPGMA, WPGMA, NJ};

	//sequences are aligned by dynamic programming and then edit distance is returned
	static int GetEditScore(std::string const& seq1, std::string const& seq2, 
							std::map<char, std::map<char, int>> const& substitution_matrix);
	static std::vector<std::string> MultipleAlignment(std::map<std::string, std::string> const& sequences,
													  std::map<char, std::map<char, int>> const& substitution_matrix, 
													  ClusterizationMethod method);

private:
	//updates distance in pairwise alignments relatively to clusters with specified number (new cluster) with UPGMA method
	static void RecalculateDistanceUPGMA(std::vector<bool> const& is_clustered,
										 std::vector<Cluster> const& clusters,
										 size_t cluster1,
										 size_t cluster2,
										 std::multimap<double, std::pair<size_t, size_t>, std::greater<double>>& pairwise_alignments,
										 std::vector <std::vector<double> >& pairwise_alignments_matrix);
	static void RecalculateDistanceNJ(std::vector<bool> const& is_clustered,
									  std::vector<Cluster> const& clusters,
									  size_t cluster1,
									  size_t cluster2,
									  std::vector <std::vector<double> >& pairwise_alignments_matrix);
	static double& GetMatrixValue(size_t i, size_t j, std::vector< std::vector<double> >& matrix);
	static std::pair<size_t, std::vector<Cluster>> ConstructGuideTree(std::vector<std::pair<std::string, std::string>> const& sequences,
									   std::map<char, std::map<char, int>> const& substitution_matrix,
									   ClusterizationMethod method);
	static std::pair<size_t, size_t> GetNJClosestAlignment(std::vector<std::vector<double>>& pairwise_alignment_matrix,
														   std::vector<bool> const& is_clustered);
	static std::pair<size_t, size_t> GetClosestAlignment(std::vector<std::vector<double>>& pairwise_alignment_matrix,
														 std::multimap<double, std::pair<size_t, size_t>, std::greater<double>>& pairwise_alignments,
														 ClusterizationMethod method,
														 std::vector<bool> const& is_clustered);
	static void RecalculateDistanceWPGMA(std::vector<bool> const& is_clustered,
										 std::vector<Cluster> const& clusters,
										 size_t cluster1,
										 size_t cluster2,
										 std::multimap<double, std::pair<size_t, size_t>, std::greater<double>>& pairwise_alignments,
										 std::vector <std::vector<double> >& pairwise_alignments_matrix);
	static void UpdateAlignmentsOrder(std::multimap<double, std::pair<size_t, size_t>, std::greater<double>>& pairwise_alignments, 
									  double key,
									  std::pair<size_t, size_t> value, double new_key);
	static AlignmentProfile ProfileAlignment(AlignmentProfile& left_profile, AlignmentProfile& right_profile, 
											 std::map<char, std::map<char, int>> const& substitution_matrix);
	static AlignmentProfile GetClusterProfile(size_t cluster, std::vector<Cluster>& guide_tree, 
											  std::map<char, std::map<char, int>> const& substitution_matrix);
	static int AlignProfileColumns(std::map<char, int> const& column1, std::map<char, int> const& column2,
									  std::map<char, std::map<char, int>> const& substitution_matrix);
	static std::map<char, int> MergeProfileColumns(std::map<char, int>& column1, std::map<char, int>& column2);
	static bool IsAllClustered(std::vector<bool> const& is_clustered);
	static void SetInfinity(std::vector <std::vector<double> >& pairwise_alignments_matrix, int cluster);
};

#endif // SequenceUtils_h__