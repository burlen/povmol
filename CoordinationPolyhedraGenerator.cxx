#include "CoordinationPolyhedraGenerator.h"

#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkIntArray.h>

#include <set>
#include <vector>
#include <utility>
#include <iostream>

using std::map;
using std::vector;
using std::pair;
using std::cerr;
using std::endl;

extern "C" {
#include <hull/hull.h>
};

// point iterator state
static site gPointListAt;
static site gPointListBegin;
static site gPointListEnd;
static vector<vtkIdType> *gRealPointIds;

// **************************************************************************
// set the points to comute comvex hull on. Initializes the iterator.
static void set_points(vector<double> &points)
{
  // intialize the iterator
  gPointListAt = &points[0];
  gPointListBegin = &points[0];
  gPointListEnd = &points[0] + points.size();
}

// **************************************************************************
// set real the point ids
static void set_real_point_ids(vector<vtkIdType> &ids)
{
  gRealPointIds = &ids;
}

// **************************************************************************
// point iterator
extern "C" { static site get_next_point(void)
{
  if (gPointListAt < gPointListEnd)
    {
    site current = gPointListAt;
    gPointListAt += 3;
    return current;
    }
  return NULL;
}
}

// **************************************************************************
// given a point return the point id
extern "C" { static long get_point_id(site p)
{
  return (p - gPointListBegin)/3;
}
}

// **************************************************************************
// given a point return the point id
extern "C" { static long get_real_point_id(site p)
{
  return gRealPointIds->at(get_point_id(p));
}
}

// output state for copying convex hull
// into vtk data structures
static vtkDoubleArray *gHullPoints;
static map<vtkIdType,vtkIdType> *gHullPointIds;
static vtkCellArray *gHullTriangles;
static vtkIntArray *gHullLabels;
static int gHullLabel;

// **************************************************************************
static void set_hull_label(int id)
{
  gHullLabel = id;
}

// **************************************************************************
static void set_hull_labels(vtkIntArray *ids)
{
  gHullLabels = ids;
}

// **************************************************************************
// set the vtk array hull coordinates are coppied into
static void set_hull_points(vtkDoubleArray *da)
{
  gHullPoints = da;
}

// **************************************************************************
// set the set used to maintain unique point list as
// points are coppied from the convex hull into vtk data
// structures
static void set_hull_point_ids(map<vtkIdType,vtkIdType> &ids)
{
  gHullPointIds = &ids;
}

// **************************************************************************
// set the vtk array hull triangles are coppied into
static void set_hull_triangles(vtkCellArray *ia)
{
  gHullTriangles = ia;
}

// **************************************************************************
// copy the trianle into VTK arrays
extern "C" { static void *copy_triangle(simplex *s, void*)
{
  vtkIdType verts[3];
  for (int j=0; j<3; ++j)
    {
    double *pt = s->neigh[j].vert;
    vtkIdType id = get_real_point_id(pt);

    // attempt an insert, if it fails use the mapped value
    // if it succeeds insert the new point
    vtkIdType nextId = gHullPoints->GetNumberOfTuples();
    auto ret = gHullPointIds->insert(pair<vtkIdType,vtkIdType>(id, gHullPoints->GetNumberOfTuples()));
    if (ret.second)
      {
      // new point
      gHullPoints->InsertNextTupleValue(pt);
      // label it
      gHullLabels->InsertNextValue(gHullLabel);
      }
    else
      {
      // using existing mapped value
      nextId = ret.first->second;
      }
    verts[j] = nextId;

    #ifdef CoordinationPolyhedraGeneratorDEBUG
    cerr << id << " " << nextId << " "
      << pt[0] << " " << pt[1] << " " << pt[2] << endl;
    #endif
    }
  gHullTriangles->InsertNextCell(3, verts);
  return NULL;
}
}

// --------------------------------------------------------------------------
CoordinationPolyhedraGenerator::CoordinationPolyhedraGenerator() :
    Points(0),
    Positions(0),
    Triangles(0),
    Labels(0),
    MinNumberOfPoints(4)
{
  this->Initialize();
  init_hull();
}

// --------------------------------------------------------------------------
CoordinationPolyhedraGenerator::~CoordinationPolyhedraGenerator()
{
  this->Clear();
}

// --------------------------------------------------------------------------
void CoordinationPolyhedraGenerator::Initialize()
{
  this->Triangles = vtkCellArray::New();
  this->Positions = vtkDoubleArray::New();
  this->Positions->SetNumberOfComponents(3);
  this->Points = vtkPoints::New();
  this->Points->SetData(this->Positions);
  this->Positions->Delete();
  this->Labels = vtkIntArray::New();
  this->Labels->SetName("label ids");
}

// --------------------------------------------------------------------------
void CoordinationPolyhedraGenerator::Clear()
{
  if (this->Points) this->Points->Delete();
  if (this->Triangles) this->Triangles->Delete();
  if (this->Labels) this->Labels->Delete();
  this->PointIds.clear();
}

// --------------------------------------------------------------------------
int CoordinationPolyhedraGenerator::Insert(
        int label,
        vector<double> &points,
        vector<vtkIdType> &ids)
{
  #ifdef CoordinationPolyhedraGeneratorDEBUG
  cerr << "CoordinationPolyhedraGenerator::Insert" << endl;
  #endif

  // avoid degenerate cases
  if (ids.size() < static_cast<size_t>(this->MinNumberOfPoints))
    {
    return 0;
    }

  set_points(points);
  set_real_point_ids(ids);

  simplex *hull = build_convex_hull(get_next_point, get_point_id, 3, 0);

  set_hull_points(this->Positions);
  set_hull_point_ids(this->PointIds);
  set_hull_triangles(this->Triangles);
  set_hull_labels(this->Labels);
  set_hull_label(label);

  visit_hull(hull, copy_triangle);

  return 1;
}
