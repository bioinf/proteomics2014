#include "fasta_utils.hpp"

void FastaReader::Read(std::ifstream &fin, std::vector<Sequence> &sequences)
{
    std::string tmp;
    std::string sequence;
    std::string id;
    while(!fin.eof()) {
        std::getline(fin, tmp);
        if (tmp[0] == '>') {
            if (!sequence.empty()) {
                sequences.push_back(Sequence(id, sequence));
            }
            id = tmp.substr(1);
            sequence.clear();
        } else {
            sequence.append(tmp);
        }
    }
    sequences.push_back(Sequence(id, sequence));
}

void FastaWriter::Write(std::ofstream &fout, const std::vector<Sequence> &sequences)
{
    for (auto it = sequences.begin(); it != sequences.end(); ++it) {
        fout << '<' << it->id << std::endl;
        fout << it->sequence << std::endl;
    }
}
