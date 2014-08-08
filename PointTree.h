#ifndef PointTree_h
#define PointTree_h

#include "Export.h"
#include <vector>

struct Node;

/// a k-d Tree for maintaining unique lists of points
class EXPORT PointTree
{
public:
  PointTree() : Tree(0) {}
  ~PointTree();

  /**
  given a set of interleaved coordinates build the tree
  list of coordnates must not be deleted until after
  this class
  */
  void Initialize(const std::vector<double> &coords);

  /// return true if the given point is already in the tree
  bool Exists(const double *point) const;

private:
  Node *Tree;
};

#endif
