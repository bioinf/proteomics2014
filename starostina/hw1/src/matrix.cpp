#include "matrix.hpp"

void Matrix::Read(std::ifstream &fin) {
    char buffer[1024];
    bool first = true;
    std::vector <char> letters;
    while (!fin.eof()) {
        std::vector <std::string> fields;
        fin.getline(buffer, 1024);
        if (buffer[0] == '#' || strlen(buffer) == 0) {
            continue;
        }
        if (first) {
            boost::split(fields, buffer, boost::is_any_of(" \t"));
            for (auto it = fields.begin(); it != fields.end(); ++it) {
                if (*it == "") {
                    continue;
                }
                if (*it == "*") {
                    letters.push_back('-');
                } else {
                    letters.push_back((*it)[0]);
                }
            }
            first = false;
        } else {
            boost::split(fields, buffer, boost::is_any_of(" \t"));
            char from = fields[0][0];
            if (from == '*') {
                from = '-';
            }
            size_t letter_i = 0;
            for (size_t i = 1; i != fields.size(); ++i) {
                if (fields[i] == "") {
                    continue;
                }
                int curr = std::atoi(fields[i].c_str());
                matrix_[std::make_pair(from, letters[letter_i])] = curr;
                ++letter_i;
            }
        }
    }
}


char Matrix::GetScore(char a, char b) const
{
    return matrix_.at(std::make_pair(a, b));
}
