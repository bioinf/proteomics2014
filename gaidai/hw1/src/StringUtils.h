#ifndef StringUtils_h__
#define StringUtils_h__

#include <string>
#include <vector>

class StringUtils
{
public:
	StringUtils();
	~StringUtils();

	static std::string ToLower(std::string const& arg);
	static std::vector<std::string> SplitAt(std::string const& to_split, std::string const& delimiter);
	static std::vector<std::string> SplitAt(std::string const& to_split, char const* delimiter);
};

#endif // StringUtils_h__
