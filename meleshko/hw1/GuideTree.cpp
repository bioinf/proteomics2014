// GuideTree.cpp : Defines the entry point for the console application.
//
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include "fasta.h"
#include "cluster.h"
typedef std::unordered_map<char, std::unordered_map<char, int>> ScoringMatrix;


enum
{
  Left, Up, UpLeft
};

template<class T>
class Node
{
public:
  Node(T* _t)
  {
    t = _t;
  }
private:
  T* t;
  T* left;
  T* right;
  std::string name;
};

class Profile
{
public:
  Profile(std::string str)
  {
    profile.resize(str.length(), std::unordered_map<char, double>());
    count.resize(str.length(), std::unordered_map<char, int>());
    picture.resize(str.length());
    for (size_t i = 0; i < str.length(); ++i)
    {
      picture[i] = str.substr(i, 1);
      profile[i][str[i]] = 1;
      count[i][str[i]] = 1;
    }
    numberOflines = 1;
    score = 1;
    penalty = 0;
  }

  Profile(Profile a, Profile b, ScoringMatrix& matrix)
  {
    //Construct profile of alignment of two other profiles
    std::vector<std::vector<double>> matrixOfAlignment(a.getSize() + 1, std::vector<double>(b.getSize() + 1, -1000));
    std::vector<std::vector<int>> path(a.getSize() + 1, std::vector<int>(b.getSize() + 1));
    numberOflines = a.numberOflines + b.numberOflines;
    double newScore = 0;
    for (size_t i = 0; i < matrixOfAlignment.size(); ++i)
    {
      matrixOfAlignment[i][0] = (int)i * a.getNumberOflines() * b.getNumberOflines() * matrix['*']['A'];
      path[i][0] = Up;
    }
    for (size_t i = 0; i < matrixOfAlignment[0].size(); ++i)
    {
      matrixOfAlignment[0][i] = (int)i * a.getNumberOflines() * b.getNumberOflines() * matrix['*']['A'];
      path[0][i] = Left;
    }
    for (size_t i = 1; i < matrixOfAlignment.size(); ++i)
    {
      for (int j = 1; j < matrixOfAlignment[0].size(); ++j)
      {

        int leftScore = 0;
        for (auto d : b.getCount()[j - 1])
        {
          leftScore += d.second * matrix['*'][d.first];
        }

        if (matrixOfAlignment[i][j] < matrixOfAlignment[i][j - 1] + leftScore)
        {
          path[i][j] = Left;
          matrixOfAlignment[i][j] = matrixOfAlignment[i][j - 1] + leftScore;
        }

        int upScore = 0;
        for (auto d : a.getCount()[i - 1])
        {
          upScore += d.second * matrix['*'][d.first];
        }

        if (matrixOfAlignment[i][j] < matrixOfAlignment[i - 1][j] + upScore)
        {
          path[i][j] = Up;
          matrixOfAlignment[i][j] = matrixOfAlignment[i - 1][j] + upScore;
        }

        int upLeftScore = 0;
        

        for (auto d : a.getCount()[i - 1])
        {
          for (auto d2 : b.getCount()[j - 1])
          {
            upLeftScore += d.second * d2.second * matrix[d.first][d2.first];
          }
        }
        if (matrixOfAlignment[i][j] < matrixOfAlignment[i - 1][j - 1] + upLeftScore)
        {
          path[i][j] = UpLeft;
          matrixOfAlignment[i][j] = matrixOfAlignment[i - 1][j - 1] + upLeftScore;
        }
      }
    }
    score = matrixOfAlignment[matrixOfAlignment.size() - 1][matrixOfAlignment[0].size() - 1];

    //backtracking
    int YPos = matrixOfAlignment[0].size() - 1;
    int XPos = matrixOfAlignment.size() - 1;
    penalty = matrixOfAlignment[XPos][YPos];

    while (XPos != 0 || YPos != 0)
    {
      if (path[XPos][YPos] == Left)
      {

        YPos--;


        std::string startPicture = "";

        for (int i = 0; i < a.numberOflines; ++i)
        {
          startPicture.push_back('*');
        }
        picture.push_back(startPicture);
        picture[picture.size() - 1].append(b.getPicture()[YPos]);

        std::unordered_map<char, int> newCount = b.getCount()[YPos];
        newCount['*'] += a.numberOflines;
        std::unordered_map<char, double> newProfile;
        for (auto it : newCount)
          newProfile[it.first] = (double)it.second / (a.numberOflines + b.numberOflines);

        profile.push_back(newProfile);
        count.push_back(newCount);
        newScore += b.getProfile()[YPos]['*'];
        continue;
      }
      if (path[XPos][YPos] == Up)
      {

        XPos--;

        std::string startPicture = "";
        for (int i = 0; i < b.numberOflines; ++i)
        {
          startPicture.push_back('*');
        }
        picture.push_back(a.getPicture()[XPos]);
        picture[picture.size() - 1].append(startPicture);
        
        std::unordered_map<char, int> newCount = a.getCount()[XPos];
        newCount['*'] += b.numberOflines;
        std::unordered_map<char, double> newProfile;
        for (auto it : newCount)
          newProfile[it.first] = (double)it.second / (a.numberOflines + b.numberOflines);

        profile.push_back(newProfile);
        count.push_back(newCount);

        newScore += a.getProfile()[XPos]['*'];
        continue;
      }
      if (path[XPos][YPos] == UpLeft)
      {

        XPos--;
        YPos--;
        picture.push_back(a.getPicture()[XPos]);
        picture[picture.size() - 1].append(b.getPicture()[YPos]);

        std::unordered_map<char, int> newCount;
        for (auto it : a.getCount()[XPos])
        {
          newScore += a.getProfile()[XPos][it.first] * b.getProfile()[YPos][it.first];
          newCount[it.first] = it.second + b.getCount()[YPos][it.first];
        }
        for (auto it : b.getCount()[YPos])
        {
          if (newCount.find(it.first) == newCount.end())
            newCount[it.first] = it.second;
        }


        std::unordered_map<char, double> newProfile;
        for (auto it : newCount)
        {
          newProfile[it.first] = (double)it.second / (double)numberOflines;
        }
        profile.push_back(newProfile);
        count.push_back(newCount);
        continue;
      }
    }
    std::reverse(picture.begin(), picture.end());
    std::reverse(count.begin(), count.end());
  }

  size_t getSize()
  {
    return profile.size();
  }

  int getScore()
  {
    return score;
  }

  double getPenalty()
  {
    return penalty;
  }

  std::vector<std::string> getPicture()
  {
    return picture;
  }

  std::vector<std::unordered_map<char, double>>& getProfile()
  {
    return profile;
  }

  std::vector<std::unordered_map<char, int>>& getCount()
  {
    return count;
  }

  int getNumberOflines()
  {
    return numberOflines;
  }

  void print(std::ofstream& out)
  {
    for (int j = 0; j < picture[0].size(); ++j)
    {
      for (int i = 0; i < picture.size(); ++i)
      {
        out << picture[i][j];
      }
      out << std::endl;
    }
  }
private:
  std::vector<std::unordered_map<char, double>> profile;
  std::vector<std::unordered_map<char, int>> count;
  std::vector<std::string> picture;
  int numberOflines;
  int score;
  double penalty;
};

void usage()
{
  std::cerr << "Usage:" << std::endl;
  std::cerr << "GuideTree sequences.fasta ScoringMatrix.txt algoType" << std::endl;
  std::cerr << "algoType = NJ | WPGMA | UPGMA" << std::endl;
}

std::vector<char> &split(const std::string &s, char delim, std::vector<char> &elems)
{
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim))
  {
    if (item[0] != '\0')
    {
      elems.push_back(item[0]);
    }
  }
  return elems;
}

std::vector<char> split(const std::string &s, char delim)
{
  std::vector<char> elems;
  split(s, delim, elems);
  return elems;
}

void ParseNCBIMatrix(ScoringMatrix& matrix, std::string filename)
{
  std::ifstream matrixInput(filename);
  std::string str;
  do
  {
    std::getline(matrixInput, str);
  } 
  while (str[0] == '#');
  std::vector<char> supmas = split(str, ' ');
  for (int i = 0; i < supmas.size(); ++i)
  {
    char ch;
    matrixInput >> ch;
    for (int j = 0; j < supmas.size(); ++j)
    {
      int z;
      matrixInput >> z;
      matrix[supmas[i]][supmas[j]] = z;
    }
  }
}

struct alInfo
{
  int firstProfile;
  int secondProfile;
  int score;
  alInfo(int _firstProfile, int _secondProfile, int _score)
    : firstProfile(_firstProfile), secondProfile(_secondProfile), score(_score)
  {  }

};

struct cmpAlInfo {
  bool operator()(const alInfo& a, const alInfo& b) const {
    return a.score < b.score;
  }
};

void UPGMA(std::vector<FastaRecord>& sequences, ScoringMatrix& matrix)
{
  std::vector<Profile> profiles;
  std::vector<cluster> clusters;
  std::unordered_set<int> indexSet;
  std::unordered_set<int> deleted;
  for (int i = 0; i < sequences.size(); ++i)
  { 
    profiles.push_back(Profile(sequences[i].sequence_));
    clusters.push_back(cluster(i));
    indexSet.insert(i);
  }
  
  std::unordered_map<int, std::unordered_map<int, int>> d;

  std::priority_queue<alInfo, std::vector<alInfo>, cmpAlInfo> queue;

  for (int i = 0; i < sequences.size(); ++i)
  {
    for (int j = i + 1; j < sequences.size(); ++j)
    {
      int temp_score = Profile(profiles[i], profiles[j], matrix).getScore();
      d[i][j] = temp_score;
      d[j][i] = temp_score;
      queue.push(alInfo(i, j, temp_score));
    }
  }

  int numberOfRound = 0;
  while (numberOfRound != sequences.size() - 1)
  {
    alInfo info = queue.top();
    queue.pop();
    if (deleted.find(info.firstProfile) != deleted.end() || deleted.find(info.secondProfile) != deleted.end())
    {
      continue;
    }
    else
    {
      indexSet.erase(info.firstProfile);
      indexSet.erase(info.secondProfile);

      deleted.insert(info.firstProfile);
      deleted.insert(info.secondProfile);
    }
    profiles.push_back(Profile(profiles[info.firstProfile], profiles[info.secondProfile], matrix));
    d[info.firstProfile][info.secondProfile] = profiles[profiles.size() - 1].getScore();
    d[info.secondProfile][info.firstProfile] = profiles[profiles.size() - 1].getScore();

    clusters.push_back(clusters[info.firstProfile]);
    clusters[clusters.size() - 1].merge(clusters[info.secondProfile]);
    
    for (auto d3 : indexSet)
    {
      int sumOfScore = 0;
      for (auto j : clusters[clusters.size() - 1].getElems())
      {
          for (auto i : clusters[d3].getElems())
          {
            sumOfScore += d[i][j];
          }
          sumOfScore /= 2;
          d[d3][j] = sumOfScore;
          d[j][d3] = sumOfScore;
      }
      queue.push(alInfo(d3, clusters.size() - 1, sumOfScore));
    }

    indexSet.insert(clusters.size() - 1);
    numberOfRound++;
  }

  std::ofstream out("output.txt");
  profiles[profiles.size() - 1].print(out);

}

void WPGMA(std::vector<FastaRecord>& sequences, ScoringMatrix& matrix)
{
  std::vector<Profile> profiles;
  std::vector<cluster> clusters;
  std::unordered_set<int> indexSet;
  std::unordered_set<int> deleted;
  for (int i = 0; i < sequences.size(); ++i)
  {
    profiles.push_back(Profile(sequences[i].sequence_));
    clusters.push_back(cluster(i));
    indexSet.insert(i);
  }

  std::unordered_map<int, std::unordered_map<int, int>> d;

  std::priority_queue<alInfo, std::vector<alInfo>, cmpAlInfo> queue;

  for (int i = 0; i < sequences.size(); ++i)
  {
    for (int j = i + 1; j < sequences.size(); ++j)
    {
      int temp_score = Profile(profiles[i], profiles[j], matrix).getScore();
      d[i][j] = temp_score;
      d[j][i] = temp_score;
      queue.push(alInfo(i, j, temp_score));
    }
  }

  int numberOfRound = 0;
  while (numberOfRound != sequences.size() - 1)
  {
    alInfo info = queue.top();
    queue.pop();
    if (deleted.find(info.firstProfile) != deleted.end() || deleted.find(info.secondProfile) != deleted.end())
    {
      continue;
    }
    else
    {
      indexSet.erase(info.firstProfile);
      indexSet.erase(info.secondProfile);

      deleted.insert(info.firstProfile);
      deleted.insert(info.secondProfile);
    }
    profiles.push_back(Profile(profiles[info.firstProfile], profiles[info.secondProfile], matrix));
    d[info.firstProfile][info.secondProfile] = profiles[profiles.size() - 1].getScore();
    d[info.secondProfile][info.firstProfile] = profiles[profiles.size() - 1].getScore();

    clusters.push_back(clusters[info.firstProfile]);
    clusters[clusters.size() - 1].merge(clusters[info.secondProfile]);

    for (auto d3 : indexSet)
    {
      int sumOfScore = 0;
      for (auto j : clusters[clusters.size() - 1].getElems())
      {
        for (auto i : clusters[d3].getElems())
        {
          sumOfScore += d[i][j];
        }
        sumOfScore /= (clusters[info.firstProfile].getSize() * clusters[info.secondProfile].getSize());
        d[d3][j] = sumOfScore;
        d[j][d3] = sumOfScore;
      }
      queue.push(alInfo(d3, clusters.size() - 1, sumOfScore));

    }

    indexSet.insert(clusters.size() - 1);
    numberOfRound++;
  }

  std::ofstream out("output.txt");
  profiles[profiles.size() - 1].print(out);

}



void ComputeQ(std::unordered_map<int, std::unordered_map<int, int>>& d, std::unordered_map<int, std::unordered_map<int, int>>& Q, std::unordered_set<int>& indexSet, std::pair<int, int>& minIndexes)
{
  int minQ = -10000000;
  for (auto i : indexSet)
  {
    for (auto j : indexSet)
    {
      if (i >= j)
        continue;
      Q[i][j] = (d.size() - 1) * d[i][j];
      for (auto k : indexSet)
      {
        if (k != i)
        {
          Q[i][j] += d[i][k];
        }

        if (k != j)
        {
          Q[i][j] += d[j][k];
        }
      }
      Q[j][i] = Q[i][j];
      if (Q[i][j] > minQ)
      {
        minQ = Q[i][j];
        minIndexes.first = i;
        minIndexes.second = j;
      }
    }
  }

}

void NJ(std::vector<FastaRecord>& sequences, ScoringMatrix& matrix)
{
  std::vector<Profile> profiles;
  for (int i = 0; i < sequences.size(); ++i)
  {
    profiles.push_back(Profile(sequences[i].sequence_));
  }
  std::set<int> deleted;
  
  std::unordered_map<int, std::unordered_map<int, int>> d;
  std::unordered_map<int, std::unordered_map<int, int>> Q;


  int numberOfRounds = sequences.size() - 1;
  for (int i = 0; i < sequences.size(); ++i)
  {
    for (int j = i + 1; j < sequences.size(); ++j)
    {
      int temp_score = Profile(profiles[i], profiles[j], matrix).getScore();
      d[i][j] = temp_score;
      d[j][i] = temp_score;
    }
  }
  std::unordered_set<int> indexSet;
  for (int i = 0; i < sequences.size(); ++i)
  {
    indexSet.insert(i);
  }

  while (numberOfRounds != 0)
  {
    std::pair<int, int> minIndexes;
    ComputeQ(d, Q, indexSet, minIndexes);
    indexSet.erase(minIndexes.first);
    indexSet.erase(minIndexes.second);
    for (auto it : indexSet)
    {
      d[profiles.size()][it] = 0.5 * (d[minIndexes.first][it] + d[minIndexes.second][it] - d[minIndexes.first][minIndexes.second]);
      d[it][profiles.size()] = d[profiles.size()][it];
    }
    indexSet.insert(profiles.size());
    profiles.push_back(Profile(profiles[minIndexes.first], profiles[minIndexes.second], matrix));
    numberOfRounds--;
  }

  std::ofstream out("output.txt");
  profiles[profiles.size() - 1].print(out);
}

int main(int argc, char* argv[])
{
  if (argc != 4)
  {
    usage();
    return 1;
  }  
  FastaReader fReader(argv[1]);
  std::vector<FastaRecord> sequences;
  fReader.GetSequences(sequences);
  ScoringMatrix matrix;
  ParseNCBIMatrix(matrix, argv[2]);
  std::string type = argv[3];
  
  if (type == "NJ")
  {
    NJ(sequences, matrix);
  }
  else
  if (type == "WPGMA")
  {
    WPGMA(sequences, matrix);
  }
  else
  if (type == "UPGMA")
  {
    UPGMA(sequences, matrix);
  }
  else
  {
    std::cerr << "Wrong alignment type";
    return 3;
  }
	return 0;
}

