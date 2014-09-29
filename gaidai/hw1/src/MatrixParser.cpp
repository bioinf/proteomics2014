#include "MatrixParser.h"
#include "StringUtils.h"

#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

using namespace std;

MatrixParser::MatrixParser()
{}

MatrixParser::~MatrixParser()
{}

void RemoveEmpty(vector<string>& vec)
{
	vec.erase(remove_if(vec.begin(), vec.end(), mem_fn(&string::empty)), vec.end());
}

std::map<char, std::map<char, int>> MatrixParser::ParseMatrix(char const* filename)
{
	ifstream matrix_file(filename);
	string line;
	while (getline(matrix_file, line) && (line[0] == '#' || line.empty()))
	{
		//skip all comments and empty
	}

	vector<string> string_tokens = StringUtils::SplitAt(line, " ");
	
	// remove empty entries
	RemoveEmpty(string_tokens);

	// convert to vector<char>
	vector<char> aa_order;
	transform(string_tokens.begin(), string_tokens.end(), back_inserter(aa_order), [](string const& str)
	{
		return ::tolower(str[0]);
	});

	map<char, map<char, int>> result;
	for (char amino_acid : aa_order)
	{
		getline(matrix_file, line);
		string_tokens = StringUtils::SplitAt(line, " ");
		RemoveEmpty(string_tokens);

		vector<int> scores;
		transform(string_tokens.begin() + 1, string_tokens.end(), back_inserter(scores), [](string const& str)
		{
			return atoi(str.c_str());
		});
		for (int i = 0; i < aa_order.size(); ++i)
		{
			result[amino_acid][aa_order[i]] = scores[i];
		}
	}

	return result;
}