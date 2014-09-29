#ifndef CommandLineParser_h__
#define CommandLineParser_h__

#include <map>
#include <string>

class CommandLineParser
{
public:
	CommandLineParser();
	~CommandLineParser();

	// unnamed args could be obtained by their number
	static std::map<std::string, std::string> ParseNamedArgs(int argc, char** argv);
};

#endif // CommandLineParser_h__