#include "CommandLineParser.h"

using namespace std;

CommandLineParser::CommandLineParser()
{}

CommandLineParser::~CommandLineParser()
{}

std::map<std::string, std::string> CommandLineParser::ParseNamedArgs(int argc, char** argv)
{
	map<string, string> result;
	string key = "";

	int args_counter = 0;
	for (int i = 1; i < argc; ++i)
	{
		string current_arg(argv[i]);
		if (!key.empty())
		{
			result.emplace(move(key), move(current_arg));
		}
		else if (current_arg[0] == '-')
		{
			key = move(current_arg);
		}
		else
		{
			result.emplace(to_string(args_counter++), move(current_arg));
		}
	}

	return move(result);
}