// proteomics2.cpp : Defines the entry point for the console application.
//
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <iomanip>
#include <random>

struct Point
{
public:
  double x;
  double y;
  double z;
  Point(double _x, double _y, double _z)
    : x(_x), y(_y), z(_z)
  { }
  Point()
    : x(0), y(0), z(0)
  {}
};


struct Atom
{
  Point coord;
  std::string prefix;
  std::string suffix;
  std::string type;
  Atom(std::string& atomPDB)
  {
    prefix = atomPDB.substr(0, 30);
    suffix = atomPDB.substr(54, 26);
    coord = Point(atof(atomPDB.substr(30, 7).c_str()), atof(atomPDB.substr(38, 7).c_str()), atof(atomPDB.substr(36, 7).c_str()));
    type = atomPDB.substr(76, 2);
    if (isspace(type[1]))
      type.pop_back();
  }
};

double distance(Point& a, Point& b)
{
  return sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2) + std::pow(a.z - b.z, 2));
}

double distance(Atom& a, Atom& b)
{
  return distance(a.coord, b .coord);
}

double distPointAxis(Point ax1, Point ax2, Point current)
{
  double xR = ax1.x - current.x;
  double yR = ax1.y - current.y;
  double zR = ax1.z - current.z;

  double crossRDirX = yR * (ax2.z - ax1.z) - zR * (ax2.y - ax1.y);
  double crossRDirY = zR * (ax2.x - ax1.x) - xR * (ax2.z - ax1.z);
  double crossRDirZ = xR * (ax2.y - ax1.y) - yR * (ax2.x - ax1.x);

  return sqrt(pow(crossRDirX, 2) + pow(crossRDirY, 2) + pow(crossRDirZ, 2)) / sqrt(pow(ax2.x - ax1.x, 2) + pow(ax2.y - ax1.y, 2) + pow(ax2.z - ax1.z, 2));
}

double crossProduct(Point ax1, Point ax2, Point ax3, Point ax4)
{
  return (ax2.x - ax1.x)*(ax4.x - ax3.x) + (ax2.y - ax1.y)*(ax4.y - ax3.y) + (ax2.y - ax1.y)*(ax4.y - ax3.y);
}

void normalize(Point& ax1, Point& ax2)
{
  double dist = distance(ax1, ax2);
  ax2.x = ax1.x + (ax2.x - ax1.x) / dist;
  ax2.y = ax1.y + (ax2.y - ax1.y) / dist;
  ax2.z = ax1.z + (ax2.z - ax1.z) / dist;
}


Point projectionPointLine(Point l1, Point l2, Point current)
{
  normalize(l1, l2);
  double a1 = crossProduct(l1, current, l1, l2);
  return Point(a1*(l2.x - l1.x), a1*(l2.y - l1.y), a1*(l2.z - l1.z));
}

double calculateAngle(Point ax1, Point ax2, Point current, Point target)
{
  double distR = distPointAxis(ax1, ax2, current);
  double distF = distPointAxis(ax1, ax2, target);

  double a = distR * distR + distF * distF;
  Point pointOnAxis = projectionPointLine(ax1, ax2, current);

  Point copyCurrent = current;
  normalize(pointOnAxis, copyCurrent);

  Point prependicular(pointOnAxis.x + (ax2.y - ax1.y)*(current.z - pointOnAxis.z) - (ax2.z - ax1.z)*(current.y - pointOnAxis.y),
    pointOnAxis.y + (ax2.x - ax1.x)*(current.z - pointOnAxis.z) - (ax2.z - ax1.z)*(current.x - pointOnAxis.x),
    pointOnAxis.z + (ax2.x - ax1.x)*(current.y - pointOnAxis.y) - (ax2.y - ax1.y)*(current.x - pointOnAxis.x));
  normalize(pointOnAxis, prependicular);
  double b = 2 * distR * crossProduct(pointOnAxis, target, pointOnAxis, copyCurrent);
  double c = 2 * distR * crossProduct(pointOnAxis, target, pointOnAxis, prependicular);
  std::default_random_engine re;
  std::uniform_real_distribution<double> dist(0, 2 * 3.14);
  if (b*b + c*c < 0.0001)
    return dist(re);
  else
    return asin(c / sqrt(b*b + c*c));
}

void rotatePointAroundAxis(Point& prev, double angle, Point ax1, Point ax2)
{
  normalize(ax1, ax2);
  double newX = (ax1.x*(ax2.y * ax2.y + ax2.z * ax2.z) - ax2.x*(ax1.y*ax2.y + ax1.z*ax2.z - prev.x*ax2.x - prev.y*ax2.y - prev.z*ax2.z)) * (1 - cos(angle)) +
    prev.x * cos(angle) + (-ax1.z*ax2.y + ax1.y*ax2.z - prev.y*ax2.z + ax2.y*prev.z) * sin(angle);
  double newY = (ax1.y*(ax2.x * ax2.x + ax2.z * ax2.z) - ax2.y*(ax1.x*ax2.x + ax1.z*ax2.z - prev.x*ax2.x - prev.y*ax2.y - prev.z*ax2.z)) * (1 - cos(angle)) +
    prev.y * cos(angle) + (ax1.z*ax2.x - ax1.x*ax2.z + prev.x*ax2.z - ax2.x*prev.z) * sin(angle);
  double newZ = (ax1.z*(ax2.y * ax2.y + ax2.x * ax2.x) - ax2.z*(ax1.y*ax2.y + ax1.x*ax2.x - prev.x*ax2.x - prev.y*ax2.y - prev.z*ax2.z)) * (1 - cos(angle)) +
    prev.z * cos(angle) + (-ax1.y*ax2.x + ax1.x*ax2.y - prev.x*ax2.y + ax2.x*prev.y) * sin(angle);
  prev.x = newX;
  prev.y = newY;
  prev.z = newZ;
}


void CCD(std::vector<Atom>& loopAtoms, Point& target)
{
  //Fix first 3 atoms, even 4
  for (int i = 3; i < loopAtoms.size() - 1; ++i)
  {
    if (loopAtoms[i].type == "C" && loopAtoms[i + 1].type == "N")
    {

      continue;
    }
    if (distance(target, loopAtoms[loopAtoms.size() - 1].coord) < 0.001)
    {
      continue;
    }
    double angle = calculateAngle(loopAtoms[i].coord, loopAtoms[i + 1].coord, loopAtoms[loopAtoms.size() - 1].coord, target);
    for (int j = i + 1; i < loopAtoms.size(); ++i)
    {
      rotatePointAroundAxis(loopAtoms[j].coord, angle, loopAtoms[i].coord, loopAtoms[i + 1].coord);
    }
  }
}


int main(int argc, char* argv[])
{
  std::string filename = argv[1];
  Point target(atof(argv[2]), atof(argv[3]), atof(argv[4]));
  std::ifstream in(filename);
  std::vector<Atom> loopAtoms;
  std::vector<std::string> filePDB;

  std::string temp = "";
  while (in >> temp)
  {
    if (temp.substr(0, 4) == "ATOM")
    {
      loopAtoms.push_back(Atom(temp));
    }
    filePDB.push_back(temp);
  }
  
  CCD(loopAtoms, target);

  std::ofstream out("output.pdb");
  int index = 0;
  for (auto d : filePDB)
  {
    if (d.substr(0, 4) == "ATOM")
    {
      out << loopAtoms[index].prefix << std::setw(7) << loopAtoms[index].coord.x << std::setw(7) << loopAtoms[index].coord.y << std::setw(7) << loopAtoms[index].coord.z << loopAtoms[index].suffix << std::endl;
      index++;
    }
    else
    {
      out << d << std::endl;
    }
  }
  return 0;
}