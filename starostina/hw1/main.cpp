#include <iostream>
#include <string>
#include <fstream>

#include "src/matrix.hpp"
#include "src/fasta_utils.hpp"
#include "src/clusterization_methods.hpp"
#include "src/alignment.hpp"
#include "src/progressive_alignment.hpp"

int main(int argc, char ** argv) {
    if (argc == 1) {
        std::cout << "./main type matrix.txt in.fasta out.fasta" << std::endl;
        std::cout << "  - type can be 'NJ', 'UPGMA' or 'WPGMA'" << std::endl;
        return 0;
    }

    std::string type_str = argv[1];
    std::string matrix_file = argv[2];
    std::string in_file = argv[3];
    std::string out_file = argv[4];

    ClusterizationType type = NJ_t;
    if (type_str.compare("UPGMA") == 0) {
        type = UPGMA_t;
    } else if (type_str.compare("WPGMA") == 0) {
        type = WPGMA_t;
    }

    std::ifstream fmatrix(matrix_file);
    Matrix matrix;
    matrix.Read(fmatrix);
    fmatrix.close();

    std::ifstream fsequences(in_file);
    FastaReader reader;
    std::vector <Sequence> sequences;
    reader.Read(fsequences, sequences);
    fsequences.close();

    ProgressiveAlignment aligner(type, matrix, sequences);
    aligner.align();

    std::ofstream fout(out_file);
    FastaWriter writer;
    writer.Write(fout, sequences);

    return 0;
}
