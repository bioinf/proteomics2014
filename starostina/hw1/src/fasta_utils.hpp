#pragma once

#include <fstream>
#include <vector>

#include "sequence.hpp"

struct FastaReader {
    static void Read(std::ifstream &fin, std::vector <Sequence> &sequences);
};

struct FastaWriter {
    static void Write(std::ofstream &fout, std::vector <Sequence> const &sequences);
};
