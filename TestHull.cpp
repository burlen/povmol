#include "CoordinationPolyhedraGenerator.h"

#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>

#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char **argv)
{
  vector<double> pts {
       0.0, 0.0, 0.0,
      -1.0, 0.0, 0.0,
       1.0, 0.0, 0.0,
       0.0, -1.0, 0.0,
       0.0,  1.0, 0.0,
       0.0, 0.0, -1.0,
       0.0, 0.0,  1.0
      };

  vector<vtkIdType> ids {
      0, 1, 2, 3, 4, 5, 6
      };

  CoordinationPolyhedraGenerator *gen = new CoordinationPolyhedraGenerator;
  gen->Insert(pts, ids);

  vtkPolyData *pd = vtkPolyData::New();
  pd->SetPolys(gen->GetTriangles());
  pd->SetPoints(gen->GetPoints());

  vtkPolyDataWriter *pdw = vtkPolyDataWriter::New();
  pdw->SetInputData(pd);
  pdw->SetFileName("gen.vtk");
  pdw->Write();

  pd->Delete();
  pdw->Delete();

  delete gen;

  return 0;
}
