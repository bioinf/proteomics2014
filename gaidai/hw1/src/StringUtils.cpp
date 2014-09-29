#include "StringUtils.h"

#include <algorithm>

using namespace std;

StringUtils::StringUtils()
{}

StringUtils::~StringUtils()
{}

std::string StringUtils::ToLower(std::string const& arg)
{
	string answer;
	transform(arg.begin(), arg.end(), back_inserter(answer), ::tolower);
	return answer;
}

std::vector<std::string> StringUtils::SplitAt(std::string const& to_split, std::string const& delimiter)
{
	vector<string> result;
	int new_pos = 0;
	int old_pos = new_pos;

	while ((new_pos = to_split.find(delimiter, old_pos)) != string::npos)
	{
		result.push_back(to_split.substr(old_pos, new_pos - old_pos));
		old_pos = new_pos + delimiter.size();
	}
	result.push_back(to_split.substr(old_pos));

	return result;
}

std::vector<std::string> StringUtils::SplitAt(std::string const& to_split, char const* delimiter)
{
	return SplitAt(to_split, string(delimiter));
}