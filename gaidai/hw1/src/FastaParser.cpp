#include "FastaParser.h"

using std::ifstream;
using std::pair;
using std::string;
using std::map;
using std::cerr;
using std::endl;

Fasta_Parser::Fasta_Parser(char const* filename)
{
	fasta.open(filename);

    if (!fasta.good())
    {
      throw "Bad file";
    }
}

Fasta_Parser::~Fasta_Parser()
{
	fasta.close();
}

// sequence (out) - next pair (header, content) in fasta stream (in lower case)
// return empty pair if eof reached
void Fasta_Parser::Get_Next_Sequence_Fasta(pair<string, string>& sequence)
{
  sequence.first.clear();
  sequence.second.clear();

  string line;
  getline(fasta, line);
  sequence.first = line;

  while (getline(fasta, line))
  {
    sequence.second += line;
    char beginning = fasta.peek();
    if (beginning == '>')
    {
      break;
    }
  }

  transform(sequence.first.begin(), sequence.first.end(), sequence.first.begin(), ::tolower);
  transform(sequence.second.begin(), sequence.second.end(), sequence.second.begin(), ::tolower);
}

// result (out) - map (header, content) of all fasta records in file (in lower case)
void Fasta_Parser::Get_All_Fasta_Sequences_Pairs(map<string, string>& result)
{
  result.clear();

  while (fasta)
  {
    pair<string, string> next_sequence;
    Get_Next_Sequence_Fasta(next_sequence);
    result[next_sequence.first] = next_sequence.second;
  }
}