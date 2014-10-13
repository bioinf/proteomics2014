#include <cassert>

#include "alignment.hpp"

PSSM::PSSM(std::string const &seq, size_t seq_ind) {
    add_sequence(seq, seq_ind);
}

size_t PSSM::size() const {
    return data_.size();
}

double PSSM::get_freq(size_t pos, char c) const {
    if (data_[pos].find(c) == data_[pos].end()) {
        return 0;
    }
    return (double)(data_[pos].at(c)) / seq_inds_.size();
}

const std::map<char, size_t> &PSSM::operator[](size_t pos) const {
    return data_[pos];
}

double PSSM::get_weight(size_t pos, char a, Matrix const & m) const {
    double w = 0;
    for (auto it = data_[pos].begin(); it != data_[pos].end(); ++it) {
        w += m.GetScore(a, it->first) * (double)it->second / seq_inds_.size();
    }
    return w;
}

double PSSM::get_weight(size_t pos, PSSM const &profile, size_t profile_pos, Matrix const &m) const {
    double w = 0;
    for (auto it = profile[profile_pos].begin(); it != profile[profile_pos].end(); ++it) {
        w += get_weight(pos, it->first, m) * it->second / profile.seq_inds_.size();
    }
    return w;
}

void PSSM::merge_and_update_sequences(PSSM const &p, AlignmentMatrix const &matrix, std::vector <Sequence> & seqs) {
    size_t i = data_.size();
    size_t j = p.data_.size();
    while (i > 0 || j > 0) {
        switch(matrix[i][j].dir) {
        case DIAG:
            for (auto it = p.data_[j - 1].begin(); it != p.data_[j - 1].begin(); ++it) {
                if (data_[i - 1].find(it->first) == data_[i - 1].end()) {
                    data_[i - 1][it->first] = 0;
                }
                data_[i - 1][it->first] += it->second;
            }
            --i;
            --j;
            break;
        case VERT:
            data_[i - 1]['-'] += p.seq_inds_.size();
            for (size_t k = 0; k != p.seq_inds_.size(); ++k) {
                size_t curr_ind = p.seq_inds_[k];
                seqs[curr_ind].sequence.insert(j, "-");
            }
            --i;
            break;
        case HOR:
            data_.insert(data_.begin() + i, p.data_[j - 1]);
            data_[i]['-'] += seq_inds_.size();
            for (size_t k = 0; k != seq_inds_.size(); ++k) {
                size_t curr_ind = seq_inds_[k];
                seqs[curr_ind].sequence.insert(i, 1, '-');
            }
            --j;
            break;
        }
    }
    seq_inds_.insert(seq_inds_.end(), p.seq_inds_.begin(), p.seq_inds_.end());
}

void PSSM::add_sequence(std::string const &seq, size_t seq_ind) {
    if (data_.empty()) {
        data_.assign(seq.size(), std::map <char, long unsigned int>());
    }
    assert(data_.size() == seq.size());
    seq_inds_.push_back(seq_ind);
    for (size_t i = 0; i != seq.size(); ++i) {
        char c = seq[i];
        if (data_[i].find(c) == data_[i].end()) {
            data_[i][c] = 0;
        }
        ++data_[i][c];
    }
}

void build_alignment_matrix(const PSSM &s1, const PSSM &s2, Matrix const &m, AlignmentMatrix &matrix) {
    matrix.assign(s1.size() + 1, std::vector <AlnCell>());
    for (size_t i = 0; i != s1.size() + 1; ++i) {
        matrix[i].assign(s2.size() + 1, AlnCell());
    }

    for (size_t i = 0; i != s1.size(); ++i) {
        matrix[i + 1][0].dir = VERT;
        matrix[i + 1][0].score = matrix[i][0].score + s1.get_weight(i, '-', m);
    }
    for (size_t i = 0; i != s2.size(); ++i) {
        matrix[0][i + 1].dir = HOR;
        matrix[0][i + 1].score = matrix[0][i].score + s2.get_weight(i, '-', m);
    }

    for (size_t i = 0; i != s1.size(); ++i) {
        for (size_t j = 0; j != s2.size(); ++j) {
            double diag_score = matrix[i][j].score + s1.get_weight(i, s2, j, m);
            double vert_score = matrix[i][j + 1].score + s1.get_weight(i, '-', m);
            double hor_score = matrix[i + 1][j].score + s2.get_weight(j, '-', m);
            if (diag_score > vert_score && diag_score > hor_score) {
                matrix[i + 1][j + 1].score = diag_score;
                matrix[i + 1][j + 1].dir = DIAG;
            } else if (vert_score > hor_score) {
                matrix[i + 1][j + 1].score = vert_score;
                matrix[i + 1][j + 1].dir = VERT;
            } else {
                matrix[i + 1][j + 1].score = hor_score;
                matrix[i + 1][j + 1].dir = HOR;
            }
        }
    }
}


double get_alignment_score(const AlignmentMatrix &matrix) {
    return matrix[matrix.size() - 1][matrix[0].size() - 1].score;
}
