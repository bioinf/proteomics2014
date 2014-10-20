#include "cluster.h"

cluster::cluster(int n)
{
  size = 1;
  elems.insert(n);
}

int cluster::getSize()
{
  return size;
}

std::unordered_set<int> cluster::getElems()
{
  return elems;
}


void cluster::merge(cluster b)
{
  size += b.getSize();
  for (auto d : b.getElems())
  {
    elems.insert(d);
  }
}
