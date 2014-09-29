#include "SequenceUtils.h"

using namespace std;

SequenceUtils::SequenceUtils()
{}

SequenceUtils::~SequenceUtils()
{}

int SequenceUtils::GetEditScore(std::string const& seq1, std::string const& seq2, 
								std::map<char, std::map<char, int>> const& substitution_matrix)
{
	string sequence1 = " " + seq1;
	string sequence2 = " " + seq2;

	vector< vector<int> > dynamic;
	dynamic.resize(sequence2.size());

	dynamic[0].resize(sequence1.size());
	dynamic[0][0] = 0;
	for (size_t i = 1; i < dynamic[0].size(); ++i)
	{
		dynamic[0][i] = substitution_matrix.at('*').at(sequence1[i]) + dynamic[0][i - 1];
	}

	for (size_t i = 1; i < dynamic.size(); ++i)
	{
		dynamic[i].resize(sequence1.size());
		dynamic[i][0] = substitution_matrix.at(sequence2[i]).at('*') + dynamic[i - 1][0];
		for (size_t j = 1; j < dynamic[i].size(); ++j)
		{
			dynamic[i][j] = max(substitution_matrix.at(sequence2[i]).at(sequence1[j]) + dynamic[i - 1][j - 1],
								max(substitution_matrix.at(sequence2[i]).at('*') + dynamic[i - 1][j],
								substitution_matrix.at('*').at(sequence1[j]) + dynamic[i][j - 1]));
		}
	}

	return dynamic.back().back();
}

std::vector<std::string> SequenceUtils::MultipleAlignment(std::map<std::string, std::string> const& sequences,
														  std::map<char, std::map<char, int>> const& substitution_matrix, 
														  ClusterizationMethod method)
{
	vector<pair<string, string>> ordered_sequences(sequences.begin(), sequences.end());
	// first - main cluster, second - vector of clusters itself
	pair<size_t, vector<Cluster>> res = ConstructGuideTree(ordered_sequences, substitution_matrix, method);

	return GetClusterProfile(res.first, res.second, substitution_matrix).first;
}

// cluster1, cluster2 - clusters that we are going to merge
void SequenceUtils::RecalculateDistanceNJ(std::vector<bool> const& is_clustered,
										  std::vector<Cluster> const& clusters,
										  size_t cluster1,
										  size_t cluster2,
										  std::vector <std::vector<double> >& pairwise_alignments_matrix)
{
	for (size_t i = 0; i < is_clustered.size(); ++i)
	{
		//need to be updated
		if (!is_clustered[i] && i != cluster1)
		{
			//distance from cluster 1 to updating cluster
			//we need to make sure that higher index goes first because we've allocated memory only to lower triangle of matrix
			double distance_1i = GetMatrixValue(cluster1, i, pairwise_alignments_matrix);
			double distance_2i = GetMatrixValue(cluster2, i, pairwise_alignments_matrix);
			double distance_12 = GetMatrixValue(cluster1, cluster2, pairwise_alignments_matrix);

			double new_distance = (distance_1i + distance_2i - distance_12) / 2;
			GetMatrixValue(cluster1, i, pairwise_alignments_matrix) = new_distance;
		}
	}
}

void SequenceUtils::UpdateAlignmentsOrder(std::multimap<double, std::pair<size_t, size_t>, std::greater<double>>& pairwise_alignments, double key,
						   std::pair<size_t, size_t> value, double new_key)
{
	auto range = pairwise_alignments.equal_range(key);
	size_t min_cluster = min(value.first, value.second);
	size_t max_cluster = max(value.first, value.second);
	for (auto it = range.first; it != range.second; ++it)
	{
		//make sure that we found required pair
		if (it->second.first == min_cluster && it->second.second == max_cluster)
		{
			pairwise_alignments.erase(it);
			pairwise_alignments.insert(make_pair(new_key, make_pair(min_cluster, max_cluster)));
			//there is unique value pairs in map so we can break loop
			break;
		}
	}
}

void SequenceUtils::RecalculateDistanceUPGMA(std::vector<bool> const& is_clustered,
											 std::vector<Cluster> const& clusters,
											 size_t cluster1,
											 size_t cluster2,
											 std::multimap<double, std::pair<size_t, size_t>, std::greater<double>>& pairwise_alignments,
											 std::vector <std::vector<double> >& pairwise_alignments_matrix)
{
	for (size_t i = 0; i < is_clustered.size(); ++i)
	{
		//need to be updated
		if (!is_clustered[i] && i != cluster1)
		{
			//distance from cluster 1 to updating cluster
			//we need to make sure that higher index goes first because we've allocated memory only to lower triangle of matrix
			double distance_1u = GetMatrixValue(cluster1, i, pairwise_alignments_matrix);
			double distance_2u = GetMatrixValue(cluster2, i, pairwise_alignments_matrix);
			size_t cluster1_size = clusters[cluster1].cluster_size;
			size_t cluster2_size = clusters[cluster2].cluster_size;
			double new_distance = (cluster1_size * distance_1u + cluster2_size * distance_2u) / (cluster1_size + cluster2_size);
			GetMatrixValue(cluster1, i, pairwise_alignments_matrix) = new_distance;

			//now we need to update corresponding key in multimap
			UpdateAlignmentsOrder(pairwise_alignments, distance_1u, make_pair(cluster1, cluster2), new_distance);
		}
	}
}

double& SequenceUtils::GetMatrixValue(size_t i, size_t j, vector< vector<double> >& matrix)
{
	return matrix[max(i, j)][min(i, j)];
}

std::pair<size_t, size_t> SequenceUtils::GetClosestAlignment(std::vector<std::vector<double>>& pairwise_alignment_matrix,
															 std::multimap<double, std::pair<size_t, size_t>, std::greater<double>>& pairwise_alignments,
															 ClusterizationMethod method)
{
	if (method == NJ)
	{
		return GetNJClosestAlignment(pairwise_alignment_matrix);
	}
	else
	{
		return pairwise_alignments.begin()->second;
	}
}

// calculates q-matrix and return alignment with min value
std::pair<size_t, size_t> SequenceUtils::GetNJClosestAlignment(std::vector<std::vector<double>>& pairwise_alignment_matrix)
{
	vector<vector<double>> q_matrix(pairwise_alignment_matrix);
	pair<size_t, size_t> closest_alignment;
	double min_distance = INFINITY;
	for (int i = 0; i < pairwise_alignment_matrix.size(); ++i)
	{
		double substrahend_i = 0;
		for (int k = 0; k < pairwise_alignment_matrix.size(); ++k)
		{
			if (k != i)
			{
				substrahend_i -= GetMatrixValue(k, i, pairwise_alignment_matrix);
			}
		}

		for (int j = i + 1; j < pairwise_alignment_matrix.size(); ++j)
		{
			double substrahend_j = 0;
			for (int k = 0; k < pairwise_alignment_matrix.size(); ++k)
			{
				if (k != j)
				{
					substrahend_j -= GetMatrixValue(k, j, pairwise_alignment_matrix);
				}
			}
			q_matrix[j][i] = pairwise_alignment_matrix[j][i] * (pairwise_alignment_matrix.size() - 2) - substrahend_i - substrahend_j;

			if (q_matrix[j][i] < min_distance)
			{
				min_distance = q_matrix[j][i];
				closest_alignment = make_pair(j, i);
			}
		}
	}

	return closest_alignment;
}

std::pair<size_t, std::vector<Cluster>> SequenceUtils::ConstructGuideTree(std::vector<std::pair<std::string, std::string>> const& sequences,
										   std::map<char, std::map<char, int>> const& substitution_matrix,
										   ClusterizationMethod method)
{
	//will be used for determine whether ith cluster is inside of other cluster or not
	//cluster leader will be a cluster with lowest number
	vector<bool> is_clustered(sequences.size());

	//matrix is symmetric so we will store only lower triangle
	vector< vector<double> > pairwise_alignments_matrix;
	pairwise_alignments_matrix.resize(sequences.size());
	for (size_t i = 1; i < pairwise_alignments_matrix.size(); ++i)
	{
		pairwise_alignments_matrix[i].resize(i);
	}

	//let's find all pairwise alignments and add them to search tree
	//key - edit distance, value - pair of cluster numbers with such distance
	multimap<double, pair<size_t, size_t>, greater<double> > pairwise_alignments;
	for (size_t i = 0; i < sequences.size(); ++i)
	{
		for (size_t j = i + 1; j < sequences.size(); ++j)
		{
			int edit_distance = GetEditScore(sequences[i].second, sequences[j].second, substitution_matrix);
			// we don't need it in NJ
			if (method != NJ)
			{
				pairwise_alignments.insert(make_pair(edit_distance, make_pair(i, j)));
			}
			pairwise_alignments_matrix[j][i] = edit_distance;
		}
	}

	//initial clustering
	vector<Cluster> clusters(sequences.size());
	for (size_t i = 0; i < clusters.size(); ++i)
	{
		clusters[i] = Cluster(sequences[i].first, sequences[i].second);
	}

	size_t last_cluster;
	while (pairwise_alignments.size() > 0)
	{
		pair<size_t, size_t> closest_alignment = GetClosestAlignment(pairwise_alignments_matrix, pairwise_alignments, method); 

		while (is_clustered[closest_alignment.first] || is_clustered[closest_alignment.second])
		{
			pairwise_alignments.erase(pairwise_alignments.begin());

			if (pairwise_alignments.size() == 0)
			{
				break;
			}

			closest_alignment = GetClosestAlignment(pairwise_alignments_matrix, pairwise_alignments, method);
		}

		if (pairwise_alignments.size() == 0)
		{
			break;
		}
		last_cluster = closest_alignment.first;

		//second number always greater than first so it's declared as clustered
		is_clustered[closest_alignment.second] = true;

		if (method == NJ)
		{
			RecalculateDistanceNJ(is_clustered,
								  clusters,
								  closest_alignment.first,
								  closest_alignment.second,
								  pairwise_alignments_matrix);
		}
		else if (method == WPGMA)
		{
			RecalculateDistanceWPGMA(is_clustered,
									 clusters,
									 closest_alignment.first,
									 closest_alignment.second,
									 pairwise_alignments,
									 pairwise_alignments_matrix);
		}
		else if (method == UPGMA)
		{
			RecalculateDistanceUPGMA(is_clustered,
									 clusters,
									 closest_alignment.first,
									 closest_alignment.second,
									 pairwise_alignments,
									 pairwise_alignments_matrix);
		}
		else
		{
			throw "Unsupported method";
		}
		

		clusters.push_back(Cluster());
		//place new cluster instead of leading, we need to keep reference to old, so save it at the end of vector
		swap(clusters.back(), clusters[closest_alignment.first]);
		//elements in vector could be moved only from first sequences.size() - 1 positions so we can safetly reference to them after moving
		//non-leading element also can't be moved after clustering
		clusters[closest_alignment.first].left = clusters.size() - 1;
		clusters[closest_alignment.first].right = closest_alignment.second;
		clusters[closest_alignment.first].cluster_size = clusters.back().cluster_size + clusters[closest_alignment.second].cluster_size;
		clusters[closest_alignment.first].distance = pairwise_alignments.begin()->first;
	}

	return make_pair(last_cluster, clusters);
}

void SequenceUtils::RecalculateDistanceWPGMA(std::vector<bool> const& is_clustered,
											 std::vector<Cluster> const& clusters,
											 size_t cluster1,
											 size_t cluster2,
											 std::multimap<double, std::pair<size_t, size_t>, std::greater<double>>& pairwise_alignments,
											 std::vector <std::vector<double> >& pairwise_alignments_matrix)
{
	for (size_t i = 0; i < is_clustered.size(); ++i)
	{
		//need to be updated
		if (!is_clustered[i] && i != cluster1)
		{
			//distance from cluster 1 to updating cluster
			//we need to make sure that higher index goes first because we've allocated memory only to lower triangle of matrix
			double distance_1u = GetMatrixValue(cluster1, i, pairwise_alignments_matrix);
			double distance_2u = GetMatrixValue(cluster2, i, pairwise_alignments_matrix);
			double new_distance = (distance_1u + distance_2u) / 2;
			GetMatrixValue(cluster1, i, pairwise_alignments_matrix) = new_distance;

			//now we need to update corresponding key in multimap
			UpdateAlignmentsOrder(pairwise_alignments, distance_1u, make_pair(cluster1, cluster2), new_distance);
		}
	}
}

int SequenceUtils::AlignProfileColumns(std::map<char, int> const& column1, std::map<char, int> const& column2,
										  std::map<char, std::map<char, int>> const& substitution_matrix)
{
	int result_score = 0;
	for (auto pair1 : column1)
	{
		for (auto pair2 : column2)
		{
			result_score += substitution_matrix.at(pair1.first).at(pair2.first) * pair1.second * pair2.second;
		}
		
	}

	return result_score;
}

std::map<char, int> SequenceUtils::MergeProfileColumns(std::map<char, int>& column1, std::map<char, int>& column2)
{
	for (auto pair : column2)
	{
		column1[pair.first] += pair.second;
	}

	return column1;
}

AlignmentProfile SequenceUtils::ProfileAlignment(AlignmentProfile& left_alignment_profile, AlignmentProfile& right_alignment_profile,
												 std::map<char, std::map<char, int>> const& substitution_matrix)
{
	vector<string> left_alignment = move(left_alignment_profile.first);
	Profile left_profile = move(left_alignment_profile.second);
	vector<string> right_alignment = move(right_alignment_profile.first);
	Profile right_profile = move(right_alignment_profile.second);

	vector< vector<int> > dynamic;
	map<char, int> left_gap_column;
	int left_alignment_size = left_alignment[0].size();
	int right_alignment_size = right_alignment[0].size();
	left_gap_column['*'] = left_alignment_size;
	map<char, int> right_gap_column;
	right_gap_column['*'] = right_alignment_size;

	enum Direction {left, up, diag};
	vector<vector<Direction>> backtrack;
	
	dynamic.resize(left_alignment.size() + 1);
	backtrack.resize(left_alignment.size() + 1);

	dynamic[0].resize(right_alignment.size() + 1);
	backtrack[0].resize(right_alignment.size() + 1);
	dynamic[0][0] = 0;
	for (size_t i = 1; i < dynamic[0].size(); ++i)
	{
		dynamic[0][i] = dynamic[0][i - 1] + AlignProfileColumns(left_gap_column, right_profile[i - 1], substitution_matrix);
		backtrack[0][i] = left;
	}

	for (size_t i = 1; i < dynamic.size(); ++i)
	{
		dynamic[i].resize(right_alignment.size() + 1);
		backtrack[i].resize(right_alignment.size() + 1);
		dynamic[i][0] = dynamic[i - 1][0] + AlignProfileColumns(left_profile[i - 1], right_gap_column, substitution_matrix);
		backtrack[i][0] = up;

		for (size_t j = 1; j < dynamic[i].size(); ++j)
		{
			int diag_value = AlignProfileColumns(left_profile[i - 1], right_profile[j - 1], substitution_matrix) + dynamic[i - 1][j - 1];
			int upper_value = AlignProfileColumns(left_profile[i - 1], right_gap_column, substitution_matrix) + dynamic[i - 1][j];
			int left_value = AlignProfileColumns(left_gap_column, right_profile[j - 1], substitution_matrix) + dynamic[i][j - 1];

			vector<int> values = {diag_value, upper_value, left_value};
			vector<Direction> directions = {diag, up, left};
			
			auto it = max_element(values.begin(), values.end());
			dynamic[i][j] = *it;
			backtrack[i][j] = directions[it - values.begin()];
		}
	}

	vector<string> alignment;
	Profile profile;
	int i = dynamic.size() - 1;
	int j = dynamic[0].size() - 1;
	while (i != 0 || j != 0)
	{
		if (backtrack[i][j] == diag)
		{
			alignment.push_back(left_alignment[i - 1] + right_alignment[j - 1]);
			profile.push_back(MergeProfileColumns(left_profile[i - 1], right_profile[j - 1]));
			--i;
			--j;
		}
		else if (backtrack[i][j] == up)
		{
			alignment.push_back(left_alignment[i - 1] + string(right_alignment_size, '-'));
			profile.push_back(MergeProfileColumns(left_profile[i - 1], right_gap_column));
			--i;
		}
		else if (backtrack[i][j] == left)
		{
			alignment.push_back(string(left_alignment_size, '-') + right_alignment[j - 1]);
			profile.push_back(MergeProfileColumns(left_gap_column, right_profile[j - 1]));
			--j;
		}
		else
		{
			throw "wrong direction";
		}
	}

	return make_pair(vector<string>(make_move_iterator(alignment.rbegin()), make_move_iterator(alignment.rend())),
					 Profile(make_move_iterator(profile.rbegin()), make_move_iterator(profile.rend())));
}

AlignmentProfile SequenceUtils::GetClusterProfile(size_t cluster, std::vector<Cluster>& guide_tree,
												  std::map<char, std::map<char, int>> const& substitution_matrix)
{
	if (guide_tree[cluster].cluster_size == 1)
	{
		// create one-string alignment and profile
		vector<string> alignment;
		Profile profile;
		for (char sym : guide_tree[cluster].sequence)
		{
			string next;
			next.push_back(sym);
			alignment.push_back(next);

			map<char, int> next_map;
			next_map[sym] = 1;
			profile.push_back(next_map);
		}

		return make_pair(alignment, profile);
	}
	
	AlignmentProfile left_profile = GetClusterProfile(guide_tree[cluster].left, guide_tree, substitution_matrix);
	AlignmentProfile right_profile = GetClusterProfile(guide_tree[cluster].right, guide_tree, substitution_matrix);
	return ProfileAlignment(left_profile, right_profile, substitution_matrix);
}