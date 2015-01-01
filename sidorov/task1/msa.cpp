#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <vector>
#include <set>
#include <stack>
#include <string>
#include <utility>
#include <algorithm>
#include <math.h>
#include <unistd.h>
#include <cstdlib>
#include <locale>
#include <queue>

#include "Tree.h"

using std::cout;
using std::endl;
using std::setprecision;
using std::map;
using std::set;
using std::vector;
using std::stack;
using std::string;
using std::stringstream;
using std::to_string;
using std::max;
using std::min;
using std::max_element;
using std::pair;
using std::make_pair;
using std::ifstream;
using std::ofstream;
using std::locale;
using std::toupper;

typedef string SeqId;
typedef vector<string> Alignment;

const size_t AA_number = 24; // 20 AA + 3 ( = |{B, Z, X}|) + 1 gap symbol
map<char, size_t> aa_indices = {{'A', 0}, {'R', 1}, {'N', 2}, {'D', 3}, 
                                {'C', 4}, {'Q', 5}, {'E', 6}, {'G', 7},
                                {'H', 8}, {'I', 9}, {'L', 10}, {'K', 11},
                                {'M', 12}, {'F', 13}, {'P', 14}, {'S', 15},
                                {'T', 16}, {'W', 17}, {'Y', 18}, {'V', 19},
                                {'B', 20}, {'Z', 21}, {'X', 22}, {'-', 23}};
const int max_column_count = 70;

void print_error(string error_message) {
    cout << "Error: " << error_message << endl;
}

int get_score(char aa_1, char aa_2, vector<vector<int>> & alignment_matrix) {
    size_t index_1 = aa_indices[aa_1];
    size_t index_2 = aa_indices[aa_2];
    return alignment_matrix[index_1][index_2];
}

int get_alignment_score(Alignment & alignment_1, Alignment & alignment_2, 
                        int index_1, int index_2, 
                        vector<vector<int>> & alignment_matrix) {
    size_t str_count_1 = alignment_1.size();
    size_t str_count_2 = alignment_2.size();
    char aa_1;
    char aa_2;
    int column_score = 0;
    for (size_t i = 0; i != str_count_1; ++i) {
        aa_1 = (index_1 < 0) ? '-' : alignment_1[i][index_1];
        for (size_t j = 0; j != str_count_2; ++j) {
            aa_2 = (index_2 < 0) ? '-' : alignment_2[j][index_2];
            column_score += get_score(aa_1, aa_2, alignment_matrix);
        }
    }
    return column_score;
}

void calc_alignment_table(Alignment & alignment_1, Alignment & alignment_2, 
                          vector<vector<int>> & alignment_matrix, 
                          vector<vector<pair<int, int>>> & table) {
    size_t al_1_length = alignment_1[0].length();
    size_t al_2_length = alignment_2[0].length();
    for (size_t i = 0; i <= al_1_length; ++i) {
        for (size_t j = 0; j <= al_2_length; ++j) {
            if (i == 0) {
                if (j == 0) { // fill out (0, 0)
                    table[i][j] = make_pair(get_alignment_score(alignment_1, 
                                                                alignment_2, 
                                                                -1, -1, 
                                                                alignment_matrix), 
                                            0);
                } else { // fill out line 0
                    int score = get_alignment_score(alignment_1, alignment_2, 
                                                    -1, j - 1, alignment_matrix);
                    table[i][j] = (j == 1) ? make_pair(score, 0) : 
                                             make_pair(table[i][j - 1].first + 
                                                       score, 0);

                }
            } else {
                if (j == 0) { // fill out column 0
                    int score = get_alignment_score(alignment_1, alignment_2, 
                                                    i - 1, -1, alignment_matrix);
                    table[i][j] = (i == 1) ? make_pair(score, 2) : 
                                             make_pair(table[i - 1][j].first + 
                                                       score, 2);
                } else { // fill out (i, j), i != 0, j != 0
                    int h_score = get_alignment_score(alignment_1, alignment_2, 
                                                      -1, j - 1, alignment_matrix);
                    int d_score = get_alignment_score(alignment_1, alignment_2, 
                                                      i - 1, j - 1, 
                                                      alignment_matrix);
                    int v_score = get_alignment_score(alignment_1, alignment_2, 
                                                      i - 1, -1, alignment_matrix);
                    int max_score = max({table[i][j - 1].first + h_score, 
                                         table[i - 1][j - 1].first + d_score, 
                                         table[i - 1][j].first + v_score});

                    int parent = 2;
                    if (max_score == table[i - 1][j - 1].first + d_score) {
                        parent = 1;
                    } else if (max_score == table[i][j - 1].first + h_score) {
                        parent = 0;
                    }
                    table[i][j] = make_pair(max_score, parent);
                }
            }
        }
    }
}

void restore_alignment(Alignment & alignment_1, Alignment & alignment_2, 
                       vector<vector<pair<int, int>>> & table, 
                       Alignment & result) {
    size_t index_1 = alignment_1[0].length() - 1;
    size_t index_2 = alignment_2[0].length() - 1;
    size_t column_index = alignment_2[0].length();
    size_t alignment_1_size = alignment_1.size();
    size_t alignment_2_size = alignment_2.size();

    for (vector<vector<pair<int, int>>>::reverse_iterator rit = table.rbegin();
         rit != table.rend(); ++rit) {

        if (rit == table.rend() - 1 && column_index == 0) {
            break;
        }

        size_t parent = (*rit)[column_index].second;
        while (parent == 0 && column_index != 0) {
            for (size_t i = 0; i != alignment_1_size; ++i) {
                result[i].push_back('-');
            }
            for (size_t i = 0; i != alignment_2_size; ++i) {
                result[i + alignment_1_size].push_back(alignment_2[i][index_2]);
            }
            --column_index;
            --index_2;
            parent = (*rit)[column_index].second;
        }

        switch (parent) {
            case 1: for (size_t i = 0; i != alignment_1_size; ++i) {
                        result[i].push_back(alignment_1[i][index_1]);
                    }
                    for (size_t i = 0; i != alignment_2_size; ++i) {
                        result[i + alignment_1_size].push_back(alignment_2[i][index_2]);
                    }
                    --column_index;
                    --index_1;
                    --index_2;
                    break;
            case 2: for (size_t i = 0; i != alignment_1_size; ++i) {
                        result[i].push_back(alignment_1[i][index_1]);
                    }
                    for (size_t i = 0; i != alignment_2_size; ++i) {
                        result[i + alignment_1_size].push_back('-');
                    }
                    --index_1;
                    break;
        }
    }
    size_t result_size = result.size();
    for (size_t i = 0; i != result_size; ++i) {
        reverse(result[i].begin(), result[i].end());
    }
}

void align(Alignment & alignment_1, Alignment & alignment_2, 
           Alignment & result, vector<vector<int>> & alignment_matrix) {
    size_t prof_1_length = alignment_1[0].length();
    size_t prof_2_length = alignment_2[0].length();
    vector<vector<pair<int, int>>> table(prof_1_length + 1, 
                                     vector<pair<int, int>>(prof_2_length + 1));
    calc_alignment_table(alignment_1, alignment_2, alignment_matrix, table);
    restore_alignment(alignment_1, alignment_2, table, result);
}

void traverse_guide_tree(Tree * node, vector<string> & str, 
                   vector<vector<int>> & alignment_matrix,
                   stack<vector<string>> & alignment_stack) {
    if (node->begin() != node->end()) {
        for (Tree::const_iterator it = node->begin(); it != node->end(); ++it) {
            traverse_guide_tree(*it, str, alignment_matrix, alignment_stack);
            if (it != node->begin()) {
                Alignment alignment_2 = alignment_stack.top();
                alignment_stack.pop();
                Alignment alignment_1 = alignment_stack.top();
                alignment_stack.pop();
                Alignment result = vector<string>(alignment_1.size() + 
                                                  alignment_2.size());
                align(alignment_1, alignment_2, result, alignment_matrix);
                alignment_stack.push(result);
            }
        }
    } else {
        alignment_stack.push(Alignment({str[node->get_num()]}));
    }
}

void msa(Tree * tree, vector<string> & str, 
         vector<vector<int>> & alignment_matrix, string output_fasta_file) {
    
    stack<Alignment> alignment_stack;
    traverse_guide_tree(tree, str, alignment_matrix, alignment_stack);
    Alignment final_alignment = alignment_stack.top();

    ofstream outfile;
    outfile.open(output_fasta_file);
    if (!outfile.good()) {
        print_error("Error on opening output file. Exit.");
    }

    int index = 0;
    int final_alignment_length = final_alignment[0].length();
    while (index < final_alignment_length) {
        for (Alignment::iterator it = final_alignment.begin(); 
             it != final_alignment.end(); ++it) {
            outfile << it->substr(index, max_column_count) << endl;
        }
        index += max_column_count;
        if (index < final_alignment_length) {
            outfile << endl;
        }
    }
    outfile.close();
}

double calc_cluster_dist(NodeNum current_id, NodeNum new_id, NodeNum id_1, 
                         NodeNum id_2, map<pair<NodeNum, NodeNum>, double> & dist, 
                         map<NodeNum, size_t> & counts, string cl_type) {

    size_t count_1 = counts[id_1];
    size_t count_2 = counts[id_2];
    NodeNum index_11 = min(id_1, current_id);
    NodeNum index_12 = max(id_1, current_id);
    NodeNum index_21 = min(id_2, current_id);
    NodeNum index_22 = max(id_2, current_id);
    double dist_current_1 = dist[make_pair(index_11, index_12)];
    double dist_current_2 = dist[make_pair(index_21, index_22)];
    if (cl_type == "wpgma") {
        return (count_1 * dist_current_1 + count_2 * dist_current_2) / 
               (count_1 + count_2);
    } else { // UPGMA
        return (dist_current_1 + dist_current_2) / 2;
    }
}

void build_pgma_tree(map<pair<NodeNum, NodeNum>, double> & dist, 
                     map<pair<NodeNum, NodeNum>, double> & clusters, 
                     map<NodeNum, Tree *> & subtrees, 
                     map<NodeNum, size_t> & counts, 
                     string cl_type, Tree * & tree) {

    size_t node_counter = counts.size();
    if (subtrees.size() == 1) {
        tree = subtrees.begin()->second;
        return;
    }

    while (subtrees.size() != 1) {
        // BLOSUM is assumed, so, distance is maximized
        double max_distance = dist.begin()->second;

        // find max distance and respective clusters (nodes)
        NodeNum id_1_tmp = (dist.begin()->first).first; //0;
        NodeNum id_2_tmp = (dist.begin()->first).second; //0;
        NodeNum id_1 = id_1_tmp;
        NodeNum id_2 = id_2_tmp;
        for (auto it = dist.begin(); it != dist.end(); ++it) {
            id_1_tmp = (*it).first.first;
            id_2_tmp = (*it).first.second;
            if ((max_distance < (*it).second)) {
              max_distance = (*it).second;
              id_1 = id_1_tmp;
              id_2 = id_2_tmp;
            }
        }

        NodeNum new_id = node_counter;

        subtrees[new_id] = new Tree(new_id);
        subtrees[new_id]->add_child(subtrees[id_1]);
        subtrees[new_id]->add_child(subtrees[id_2]);
        subtrees.erase(id_1);
        subtrees.erase(id_2);

        // calc distances for new cluster
        NodeNum current_id;
        for (map<NodeNum, Tree *>::iterator it = subtrees.begin(); 
             it != subtrees.end(); ++it) {
            current_id = it->first;
            if (current_id != id_1 && current_id != id_2 && current_id != new_id) {
                double distance = calc_cluster_dist(current_id, new_id, id_1, 
                                                    id_2, dist, counts, cl_type);
                dist[make_pair(current_id, new_id)] = distance;
            }
        }

        // remove distances for clusters id_1 and id_2
        vector<map<pair<NodeNum, NodeNum>, double>::iterator> it_to_delete;
        for (auto it = dist.begin(); it != dist.end(); ++it) {
            if ((*it).first.first == id_1 || (*it).first.first == id_2 ||
                (*it).first.second == id_1 || (*it).first.second == id_2) {
                
                it_to_delete.push_back(it);
            }
        }
        for (auto it = it_to_delete.begin(); it != it_to_delete.end(); ++it) {
            dist.erase(*it);
        }

        // correct node counts in clusters
        counts[new_id] = counts[id_1] + counts[id_2];
        counts.erase(id_1);
        counts.erase(id_2);

        // a new cluster was added
        ++node_counter;
    }
    tree = subtrees.begin()->second;
}

double sum_dists(NodeNum node, map<pair<NodeNum, NodeNum>, double> & dist) {
    double result = 0;
    for (auto node_pair : dist) {
        if (node_pair.first.first == node || node_pair.first.second == node) {
            result += node_pair.second;
        }
    }
    return result;
}

void remove_dists(NodeNum node, map<pair<NodeNum, NodeNum>, double> & dist) {
    for (auto it = dist.begin(); it != dist.end(); ) {
        if (it->first.first == node || it->first.second == node) {
            dist.erase(it++);
        } else {
            ++it;
        }
    }
}

void remove_from_Q(NodeNum node, map<pair<NodeNum, NodeNum>, double> & Q) {
    map<pair<NodeNum, NodeNum>, double>::iterator it = Q.begin();
    while (it != Q.end()) {
        if ((it->first).first == node || (it->first).second == node) {
           map<pair<NodeNum, NodeNum>, double>::iterator toErase = it;
           ++it;
           Q.erase(toErase);
        } else {
           ++it;
        }
    }
}

void nj_merge_nodes(NodeNum num_1, NodeNum num_2, Tree * & tree, 
                    size_t node_count) {
    Tree * node_1 = tree->find(num_1); // find can return nullptr
    Tree * node_2 = tree->find(num_2); // find can return nullptr
    Tree * node_1_parent = node_1->get_parent();
    Tree * node_2_parent = node_2->get_parent();
    node_1_parent->remove_child(node_1);
    node_2_parent->remove_child(node_2);
    Tree * new_node = new Tree(node_count);
    tree->add_child(new_node);
    new_node->add_child(node_1);
    new_node->add_child(node_2);
    node_1->set_parent(new_node);
    node_2->set_parent(new_node);
}

struct CompareSecond {
    bool operator()(const pair<pair<NodeNum, NodeNum>, double> & left, 
                    const pair<pair<NodeNum, NodeNum>, double> & right) const {
        return left.second < right.second;
    }
};

double get_dist(map<pair<NodeNum, NodeNum>, double> & dist, 
                NodeNum i, NodeNum j) {
    return (dist[make_pair(i, j)] != 0) ? dist[make_pair(i, j)] : 
                                          dist[make_pair(j, i)];
}

void build_nj_tree(map<pair<NodeNum, NodeNum>, double> & dist, 
                   map<pair<NodeNum, NodeNum>, double> & Q, int seq_count, 
                   Tree * & tree) {

    // construct initial tree and node list
    set<NodeNum> actual_nodes;
    for (NodeNum i = 0; i < seq_count; ++i) {
        Tree * child = new Tree(i);
        tree->add_child(child);
        actual_nodes.insert(i);
    }

    size_t actual_node_count = seq_count;
    size_t node_count = seq_count + 1;
    // build final tree from initial one
    for (int iter = 0; iter < seq_count - 3; ++iter) {
        // calc Q matrix from dist matrix
        for (auto dist_pair : dist) {
            NodeNum node_1 = dist_pair.first.first;
            NodeNum node_2 = dist_pair.first.second;
            Q[make_pair(node_1, node_2)] = (actual_node_count - 2) * 
                get_dist(dist, node_1, node_2) - sum_dists(node_1, dist) - 
                sum_dists(node_2, dist);
        }

        // find first MAX element in Q (BLOSUM assumed, so maximize)
        auto max_it = max_element(Q.begin(), Q.end());
        NodeNum max_i = max_it->first.first;
        NodeNum max_j = max_it->first.second;

        // merge found nodes and calc distances from the new node
        nj_merge_nodes(max_i, max_j, tree, node_count);

        NodeNum new_node = node_count;
        actual_nodes.insert(new_node);
        actual_nodes.erase(max_i);
        actual_nodes.erase(max_j);
        ++node_count;
        for (auto node : actual_nodes) {
            if (node != new_node) {
                dist[make_pair(new_node, node)] = 0.5 * 
                    (get_dist(dist, max_i, node) + get_dist(dist, max_j, node)
                    - get_dist(dist, max_i, max_j));
            }
        }

        // remove max_i and max_j nodes from dist matrix and from the list
        remove_dists(max_i, dist);
        remove_dists(max_j, dist);
        remove_from_Q(max_i, Q);
        remove_from_Q(max_j, Q);
        
        --actual_node_count;
    }
}

int calc_edit_dist(size_t index_1, size_t index_2, vector<string> & str,
                      vector<vector<int>> & alignment_matrix) {
    string seq_1 = str[index_1];
    string seq_2 = str[index_2];
    size_t seq_1_length = seq_1.length();
    size_t seq_2_length = seq_2.length();
    vector<vector<int>> table(seq_1_length + 1, 
                              vector<int>(seq_2_length + 1));
    for (size_t i = 0; i <= seq_1_length; ++i) {
        for (size_t j = 0; j <= seq_2_length; ++j) {
            if (i == 0) {
                if (j == 0) { // fill out (0, 0)
                    table[i][j] = 0;
                } else { // fill out line 0
                    int score = get_score('-', seq_2[j - 1], alignment_matrix);
                    table[i][j] = (j == 1) ? score : table[i][j - 1] + score; 
                }
            } else {
                if (j == 0) { // fill out column 0
                    int score = get_score(seq_1[i - 1], '-', alignment_matrix);
                    table[i][j] = (i == 1) ? score : table[i - 1][j] + score;
                } else { // fill out (i, j), i != 0, j != 0
                    int h_score = get_score(seq_1[i - 1], '-', alignment_matrix);
                    int d_score = get_score(seq_1[i - 1], seq_2[j - 1], 
                                            alignment_matrix);
                    int v_score = get_score('-', seq_2[j - 1], alignment_matrix);
                    table[i][j] = max({table[i][j - 1] + h_score, 
                                       table[i - 1][j - 1] + d_score, 
                                       table[i - 1][j] + v_score});
                }
            }
        }
    }
    return table[seq_1_length][seq_2_length];
}

void calc_distance_matrix(vector<string> & str, 
                          vector<vector<int>> & alignment_matrix,
                          map<pair<NodeNum, NodeNum>, double> & dist) {
    size_t str_count = str.size();
    for (size_t i = 0; i < str_count; ++i) {
        for (size_t j = i + 1; j < str_count; ++j) {
            dist[make_pair(i, j)] = calc_edit_dist(i, j, str, alignment_matrix);
        }
    }
}

bool read_input(string input_fasta_file, vector<string> & str, 
                vector<string> & seq_ids, string matrix_file, 
                vector<vector<int>> & alignment_matrix) {
    // read FASTA file
    ifstream fasta_infile;
    locale loc;
    fasta_infile.open(input_fasta_file);
    if (!fasta_infile.good()) {
        print_error("Error on opening FASTA input file (maybe, no such file). \
Exit.");
        return false;
    }
    string tmp;
    string seq;
    getline(fasta_infile, tmp); // read id
    seq_ids.push_back(tmp);
    while (!fasta_infile.eof()) {
        getline(fasta_infile, tmp); // read first line of the i-th DNA sequence
        while (!tmp.empty() && tmp[0] != '>') {
            seq += tmp;
            getline(fasta_infile, tmp);
        }
        if (seq != "") {
            size_t seq_length = seq.length();
            for (size_t index = 0; index < seq_length; ++index) {
                seq[index] = toupper(seq[index], loc);
            }
            str.push_back(seq);
        }
        seq = "";
        if (!tmp.empty()) {
            seq_ids.push_back(tmp);
        }
    }
    for (auto& id : seq_ids) {
        id = id.substr(1);
    }

    // read alignment matrix
    ifstream matrix_infile;
    matrix_infile.open(matrix_file);
    if (!matrix_infile.good()) {
        print_error("Error on opening alignment matrix input file (maybe, no \
such file). Exit.");
        return false;
    }
    alignment_matrix.resize(AA_number);
    size_t line_counter = 0;
    if (!matrix_infile.eof()) {
        getline(matrix_infile, tmp);
        while (tmp[0] == '#' || tmp[0] == ' ') { 
            // it meets NCBI format for alignment scoring matrices                                                        
            getline(matrix_infile, tmp); // read lines while they begin with #
                                         // or with ' ' (the head of the table)
        }
        double score;
        string aa_abbr; // AA abbreviation each matrix line begins with 
        for (size_t j = 0; j < AA_number; ++j) {
            stringstream tmp_stream(tmp);
            tmp_stream >> aa_abbr; // just is read, isn't used
            for (size_t i = 0; i < AA_number; ++i) {
                tmp_stream >> score;
                alignment_matrix[line_counter].push_back(score);
            }
            if (j < AA_number - 1) {
                getline(matrix_infile, tmp);
                ++line_counter;
            }
        }
    }

    return true;
}

void print_synopsis() {
	cout << "Usage: msa -f <input_fasta_file> -m <matrix_file> -c \
<cl_type> -o <output_fasta_file>" << endl;
    cout << "clusterization type: nj | upgma | wpgma" << endl;
}

int main(int argc, char ** argv) {
    string input_fasta_file;
    string output_fasta_file;
    string matrix_file;
    string cl_type;

	int opt;
	int f_flag = -1;
	int m_flag = -1;
	int c_flag = -1;
	int o_flag = -1;
	while ((opt = getopt(argc, argv, "f:m:c:o:")) != -1){
		switch (opt){
			case 'f':
				if (f_flag == 1) {
					print_synopsis();
					return EXIT_FAILURE;
				}
				f_flag = 1;
				input_fasta_file = optarg;
				break;
			case 'm':
				if (m_flag == 1) {
					print_synopsis();
					return EXIT_FAILURE;
				}
				m_flag = 1;
				matrix_file = optarg;
				break;
			case 'c':
				if (c_flag == 1) {
					print_synopsis();
					return EXIT_FAILURE;
				}
				c_flag = 1;
				cl_type = optarg;
				break;
			case 'o':
				if (o_flag == 1) {
					print_synopsis();
					return EXIT_FAILURE;
				}
				o_flag = 1;
				output_fasta_file = optarg;
				break;
			case '?':
				print_synopsis();
				return EXIT_FAILURE;
        };
	};

	if (f_flag == -1 || m_flag == -1 || c_flag == -1 || o_flag == -1) {
		print_synopsis();
		return EXIT_FAILURE;
	}

    if (cl_type != "nj" && cl_type != "upgma" && cl_type != "wpgma") {
        print_synopsis();
		return EXIT_FAILURE;
    }
    
    vector<string> str;
    vector<string> seq_ids;
    vector<vector<int>> alignment_matrix;

    if (read_input(input_fasta_file, str, seq_ids, matrix_file, 
        alignment_matrix)) {

        size_t seq_count = seq_ids.size();
        if (seq_count <= 1) {
            print_error("There must be 2 or more sequences in FASTA file. Exit.");
        }

        map<pair<NodeNum, NodeNum>, double> dist;
        calc_distance_matrix(str, alignment_matrix, dist);

        Tree * tree = nullptr;
        if (cl_type == "nj") { // NJ
            map<pair<NodeNum, NodeNum>, double> Q;
            NodeNum root = seq_count;
            tree = new Tree(root);
            build_nj_tree(dist, Q, seq_count, tree);
        } else { // PGMA
            map<pair<NodeNum, NodeNum>, double> clusters;
            map<NodeNum, Tree *> subtrees;
            map<NodeNum, size_t> counts;
            tree = new Tree();
            for (size_t i = 0; i < seq_count; ++i) {
                subtrees[i] = new Tree(i);
                counts[i] = 1;
            }
            build_pgma_tree(dist, clusters, subtrees, counts, cl_type, tree);
        }
        msa(tree, str, alignment_matrix, output_fasta_file);
    }
    return 0;
}

