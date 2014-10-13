#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <boost/algorithm/string.hpp>

class Matrix {
public:
    Matrix() { }

    void Read(std::ifstream &fin);
    char GetScore(char a, char b) const;

private:
    std::map <std::pair <char, char>, int> matrix_;
};
