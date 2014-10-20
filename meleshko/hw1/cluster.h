#include <unordered_set>


struct cluster
{
public:
  cluster(int n);
  int getSize();
  std::unordered_set<int> getElems();
  void merge(cluster b);
private:
  int size;
  std::unordered_set<int> elems;
};