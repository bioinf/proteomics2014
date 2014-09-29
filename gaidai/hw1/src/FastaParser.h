#pragma once

#include <fstream>
#include <string>
#include <map>
#include <iostream>
#include <algorithm>

class Fasta_Parser
{
public:
	explicit Fasta_Parser(char const* filename);
	~Fasta_Parser();

    void Get_Next_Sequence_Fasta(std::pair<std::string, std::string>& sequence);
    void Get_All_Fasta_Sequences_Pairs(std::map<std::string, std::string>& result);

private:
	std::ifstream fasta;
};