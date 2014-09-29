#ifndef MatrixParser_h__
#define MatrixParser_h__

#include <map>

class MatrixParser
{
public:
	MatrixParser();
	~MatrixParser();

	static std::map<char, std::map<char, int>> ParseMatrix(char const* filename);
};

#endif // MatrixParser_h__