#include "FastaParser.h"
#include "CommandLineParser.h"
#include "StringUtils.h"
#include "MatrixParser.h"
#include "SequenceUtils.h"
#include <algorithm>
#include <iostream>

using namespace std;

void PrintAlignment(vector<string> const& alignment)
{
	for (int i = 0; i < alignment[0].size(); ++i)
	{
		for (string col : alignment)
		{
			cout << col[i];
		}
		cout << endl;
	}
}

int main(int argc, char** argv)
{
	map<string, string> args = CommandLineParser::ParseNamedArgs(argc, argv);
	args["-method"] = StringUtils::ToLower(args["-method"]);
	string method_string = args["-method"];
	SequenceUtils::ClusterizationMethod method = SequenceUtils::UPGMA;
	if (method_string == "wpgma")
	{
		method = SequenceUtils::WPGMA;
	}
	if (method_string == "nj")
	{
		method = SequenceUtils::NJ;
	}

	Fasta_Parser fasta_parser(args["-file"].c_str());
	map<string, string> proteins;
	fasta_parser.Get_All_Fasta_Sequences_Pairs(proteins);

	map<char, map<char, int>> matrix = MatrixParser::ParseMatrix(args["-matrix"].c_str());
	
	vector<string> alignment = SequenceUtils::MultipleAlignment(proteins, matrix, method);
	PrintAlignment(alignment);
}