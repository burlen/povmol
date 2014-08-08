#include "PointTree.h"
#include <cstddef>
#include <vector>
#include <algorithm>
#include <cmath>

using std::vector;

/**
Datastructure to work with a point stored indirectly
*/
struct Point
{
  Point() : Points(0), Id(0) {}
  Point(const double *points, size_t id) : Points(points), Id(id) {}
  ~Point() {}
  double operator[](size_t i) const { return Coordinates()[i]; }
  const double *Coordinates() const { return Points + 3*Id; }
  const double *Points;
  size_t Id;
};

/**
Node in the tree, stores partiton direction and partiton point
*/
struct Node
{
  Node() : Direction(0), LeftSubtree(0), RightSubtree(0) {}
  Node(size_t d) : Direction(d), LeftSubtree(0), RightSubtree(0) {}
  ~Node() { delete LeftSubtree; delete RightSubtree; }
  const Node *Subtree(const double *point) const
  {
    if (point[Direction] < Partition[Direction]) return LeftSubtree;
    return RightSubtree;
  }
  size_t Direction;
  Point Partition;
  Node *LeftSubtree;
  Node *RightSubtree;
};

/**
Compares points in one coordinate dircetion
*/
struct DirectionLess
{
  DirectionLess() : D(0) {}
  DirectionLess(size_t d) : D(d) {}
  bool operator()(const Point &lhs, const Point &rhs) const
  { return lhs[D] < rhs[D]; }
  size_t D;
};

// **************************************************************************
inline
bool equal(double a, double b, double eps) { return fabs(a - b) < eps; }

// **************************************************************************
void Build(Node *node, vector<Point> &points, size_t depth)
{
  size_t dir = depth%3;
  // choose partiton
  std::sort(points.begin(), points.end(), DirectionLess(dir));
  size_t partition = points.size()/2;
  node->Partition = points[partition];
  node->Direction = dir;
  // build left subtree
  if (partition)
    {
    node->LeftSubtree = new Node;
    vector<Point> subtreePoints(points.begin(), points.begin()+partition);
    Build(node->LeftSubtree, subtreePoints, depth+1);
    }
  // build right subtree
  if (partition+1 < points.size())
    {
    node->RightSubtree = new Node;
    vector<Point> subtreePoints(points.begin()+partition+1, points.end());
    Build(node->RightSubtree, subtreePoints, depth+1);
    }
}

// --------------------------------------------------------------------------
PointTree::~PointTree()
{
  delete Tree;
}

// --------------------------------------------------------------------------
void PointTree::Initialize(const vector<double> &coords)
{
  delete Tree;
  Tree = 0;
  size_t nPoints = coords.size()/3;
  if (nPoints)
    {
    const double *pCoords = &coords[0];
    vector<Point> points(nPoints);
    for (size_t i = 0; i<nPoints; ++i)
      {
      points[i] = Point(pCoords, i);
      }
    Tree = new Node;
    Build(Tree, points, 0);
    }
}

// --------------------------------------------------------------------------
bool PointTree::Exists(const double *point) const
{
  const Node *node = Tree;
  while (node)
    {
    const double *nodePt = node->Partition.Coordinates();
    if ( equal(nodePt[0], point[0], 1e-5)
      && equal(nodePt[1], point[1], 1e-5)
      && equal(nodePt[2], point[2], 1e-5) )
      {
      return true;
      }
    node = node->Subtree(point);
    }
  return false;
}
