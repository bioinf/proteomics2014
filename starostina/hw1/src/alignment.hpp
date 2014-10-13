#pragma once

#include <map>
#include <vector>

#include "matrix.hpp"
#include "sequence.hpp"

enum Dir {
    DIAG, VERT, HOR
};

struct AlnCell {
    AlnCell() : score(0), dir(DIAG) { }

    double score;
    Dir dir;
};

typedef std::vector <std::vector <AlnCell> > AlignmentMatrix;

class PSSM {
public:
    PSSM() { }
    PSSM(std::string const &seq, size_t seq_ind);

    size_t size() const;
    double get_freq(size_t pos, char c) const;
    std::map <char, size_t> const & operator[](size_t pos) const;
    double get_weight(size_t pos, char a, Matrix const &m) const;
    double get_weight(size_t pos, PSSM const &profile, size_t profile_pos, Matrix const &m) const;

    void add_sequence(std::string const &seq, size_t seq_ind);
    void merge_and_update_sequences(PSSM const &p, AlignmentMatrix const &matrix, std::vector <Sequence> &seqs);

private:
    std::vector <std::map <char, size_t> > data_;
    std::vector <size_t> seq_inds_;
};

void build_alignment_matrix(PSSM const & s1, PSSM const & s2, Matrix const &m, AlignmentMatrix &matrix);
double get_alignment_score(AlignmentMatrix const &matrix);
