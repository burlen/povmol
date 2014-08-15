#ifndef CoordinationPolyhedra_h
#define CoordinationPolyhedra_h

class vtkPoints;
class vtkDoubleArray;
class vtkCellArray;
class vtkIntArray;

#include <vtkType.h>
#include <map>
#include <vector>

/**
Helper class to generate coordination polyhedra for
crystal structures.
*/
class CoordinationPolyhedraGenerator
{
public:
  CoordinationPolyhedraGenerator();
  ~CoordinationPolyhedraGenerator();

  /**
  Clear out internal reosurces, release memory etc.
  */
  void Clear();

  /**
  Allocate interanl data structures
  */
  void Initialize();

  /**
  Insert a polyhedra for the given points. Returns
  0 if the polyhedra was not generated,
  */
  int Insert(
        int label,
        std::vector<double> &points,
        std::vector<vtkIdType> &ids);

  /**
  Returns a VTK data array containing all of the points
  for all of the polyhedra.
  */
  vtkPoints *GetPoints(){ return this->Points; }
  vtkCellArray *GetTriangles(){ return this->Triangles; }
  vtkIntArray *GetLabels(){ return this->Labels; }

  /**
  Skip generation when given fewer than this number
  of input points.
  */
  void SetMinNumberOfPoints(int n){ this->MinNumberOfPoints = n; }

private:
  vtkPoints *Points;
  vtkDoubleArray *Positions;
  vtkCellArray *Triangles;
  vtkIntArray *Labels;
  std::map<vtkIdType,vtkIdType> PointIds;
  int MinNumberOfPoints;
};

#endif
