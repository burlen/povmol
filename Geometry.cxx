#include "Geometry.h"
#include "Utility.h"
#include "AtomicProperties.h"
#include "Math3D.h"

#include <vtkMolecule.h>
#include <vtkPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>

#include <vector>
using std::vector;

// ---------------------------------------------------------------------------
void BuildAxes(
    vtkPolyData *data,
    const vector<double> &axes)
{
  if (axes.size() != 9)
    {
    pError(cerr) << "invalid axes points" << endl;
    return;
    }

  vtkDoubleArray *pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents(3);
  pts->SetNumberOfTuples(4);
  double *pPts = pts->GetPointer(0);
  pPts[0] = 0.0;
  pPts[1] = 0.0;
  pPts[2] = 0.0;
  pPts += 3;
  for (size_t i = 0; i<9; ++i)
    {
    pPts[i] = axes[i];
    }

  vtkPoints *points = vtkPoints::New();
  points->SetData(pts);
  pts->Delete();

  data->SetPoints(points);
  points->Delete();

  vtkIdTypeArray *ptIds = vtkIdTypeArray::New();
  ptIds->SetNumberOfTuples(9);
  vtkIdType *pPtIds = ptIds->GetPointer(0);
  pPtIds[0] = 2;
  pPtIds[1] = 0;
  pPtIds[2] = 1;
  pPtIds[3] = 2;
  pPtIds[4] = 0;
  pPtIds[5] = 2;
  pPtIds[6] = 2;
  pPtIds[7] = 0;
  pPtIds[8] = 3;

  vtkCellArray *cells = vtkCellArray::New();
  cells->SetCells(3, ptIds);
  ptIds->Delete();

  data->SetLines(cells);

  vtkUnsignedCharArray *colors = vtkUnsignedCharArray::New();
  colors->SetNumberOfComponents(3);
  colors->SetNumberOfTuples(3);
  colors->SetName("colors");
  unsigned char *pColors = colors->GetPointer(0);
  memset(pColors, 0, 9);
  pColors[0] = 255;
  pColors[4] = 255;
  pColors[8] = 255;
  data->GetCellData()->AddArray(colors);
  data->GetCellData()->SetActiveScalars("colors");
  colors->Delete();
}

// ---------------------------------------------------------------------------
void BuildMolecule(
    vtkMolecule *molecule,
    vector<double> &positions,
    vector<unsigned short> &numbers,
    int detectionMode,
    double proximityFactor)
{
  molecule->Initialize();
  // atoms
  unsigned long natoms = numbers.size();
  for (unsigned long i=0; i<natoms; ++i)
    {
    if (!numbers[i]) continue;
    unsigned long ii = 3*i;
    molecule->AppendAtom(numbers[i], positions[ii], positions[ii+1], positions[ii+2]);
    }
  // bonds
  for (unsigned int i=0; i<natoms; ++i)
    {
    if (!numbers[i]) continue;
    unsigned int ii = 3*i;
    vector<unsigned long> bonds;
    for (unsigned int j=i+1; j<natoms; ++j)
      {
      if (!numbers[j]) continue;
      unsigned int jj = 3*j;
      double r[3];
      Math3D::sub(r, &positions[jj], &positions[ii]);
      if ( Math3D::length(r)
        < proximityFactor*( AtomicProperties::Radius(numbers[i], detectionMode)
                + AtomicProperties::Radius(numbers[j], detectionMode) ) )
        {
        bonds.push_back(j);
        }
      }
    unsigned long nTargets = bonds.size();
    for (unsigned long j=0; j<nTargets; ++j)
      {
      molecule->AppendBond(i, bonds[j]);
      }
    }
}
