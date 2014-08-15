#include "Utility.h"
#include "Math3D.h"
#include "Geometry.h"
#include "AtomicProperties.h"
#include "AtomicPropertiesBondDetector.h"

#include <vtkMolecule.h>
#include <vtkStructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>

#include <iostream>
#include <vector>
#include <string>

using std::vector;
using std::string;
using std::cerr;
using std::endl;

namespace {

/**
Read and transform atoms from XSF file, store positions in vector
of doubles ax and forces in af. points and forces are stores as
3 tuples. store atomic numbers for each atom in an.
*/
int ReadXSFAtoms(
    ifstream &file,
    vector<double> &ax,
    vector<double> &af,
    vector<unsigned short> &an);

/**
Read and transform grid from XSF file.
*/
int ReadXSFGrid(ifstream &file, vtkStructuredGrid *sgrid);

// ---------------------------------------------------------------------------
int ReadXSFAtoms(
    ifstream &file,
    vector<double> &ax,
    vector<double> &af,
    vector<unsigned short> &an)
{
  string header;
  file >> header;
  if (!file.good() || (header != "PRIMVEC"))
    {
    cerr << "ERROR: PRIMVEC not found" << endl;
    return -3;
    }

  double A[3] = {1.0, 0.0, 0.0};
  double B[3] = {0.0, 1.0, 0.0};
  double C[3] = {0.0, 0.0, 1.0};

  file
    >> A[0] >> A[1] >> A[2]
    >> B[0] >> B[1] >> B[2]
    >> C[0] >> C[1] >> C[2];
  if (!file.good())
    {
    cerr << "ERROR: PRIMVEC values not found" << endl;
    return -4;
    }
  double a[3];
  Math3D::scale(a,A,0.5);

  double b[3];
  Math3D::scale(b,B,0.5);

  double r0[3];
  Math3D::add(r0,a,b);

  double alpha[4];
  Math3D::defplane(alpha,A,r0);

  double beta[4];
  Math3D::defplane(beta,B,r0);


  file >> header;
  if (!file.good() || (header != "PRIMCOORD" ))
    {
    cerr << "ERROR: PRIMCOORD not found" << endl;
    return -5;
    }

  unsigned long nMols = 0;
  unsigned long nMolsPerLine = 1;
  file >> nMols >> nMolsPerLine;
  if (!file.good())
    {
    cerr << "ERROR: number of molecules not found" << endl;
    return -6;
    }

  // atoms
  unsigned long npts = 3*nMols;
  ax.resize(npts);
  af.resize(npts);
  an.resize(nMols);
  for (unsigned long j=0; j<nMols; ++j)
    {
    unsigned long jj = 3*j;
    file
      >> an[j]
      >> ax[jj] >> ax[jj+1] >> ax[jj+2]
      >> af[jj] >> af[jj+1] >> af[jj+2];
    if (!file.good())
      {
      cerr << "ERROR: bad atom data " << jj << endl;
      return -7;
      }

    double *pax = &ax[jj];
    if ((Math3D::evalplane(alpha,pax) < 0.0) && (Math3D::evalplane(beta,pax) < 0.0))
      {
      Math3D::add(pax,pax,a);
      Math3D::add(pax,pax,b);
      }
    else
    if ((Math3D::evalplane(alpha,pax) >= 0.0) && (Math3D::evalplane(beta,pax) >= 0.0))
      {
      Math3D::sub(pax,pax,a);
      Math3D::sub(pax,pax,b);
      }
    else
    if ((Math3D::evalplane(alpha,pax) < 0.0) && (Math3D::evalplane(beta,pax) >= 0.0))
      {
      Math3D::add(pax,pax,a);
      Math3D::sub(pax,pax,b);
      }
    else
    if ((Math3D::evalplane(alpha,pax) >= 0.0) && (Math3D::evalplane(beta,pax) < 0.0))
      {
      Math3D::sub(pax,pax,a);
      Math3D::add(pax,pax,b);
      }
    else
      {
      cerr << "ERROR: not in any case" << endl;
      }
    }
  return 0;
}

// ---------------------------------------------------------------------------
int ReadXSFGrid(ifstream &file, vtkStructuredGrid *sgrid)
{
  string header;
  file >> header;
  if (!file.good() || (header != "BEGIN_BLOCK_DATAGRID_3D"))
    {
    cerr << "ERROR: BEGIN_BLOCK_DATAGRID_3D not found" << endl;
    return -8;
    }

  string comment;
  file >> comment;
  if (!file.good())
    {
    cerr << "ERROR: comment not found" << endl;
    return -9;
    }

  string grid;
  file >> grid;
  if (!file.good() || (grid.find("BEGIN_DATAGRID_3D")))
    {
    cerr << "ERROR: BEGIN_DATAGRID_3D not found" << endl;
    return -10;
    }
  string arrayName = grid.substr(18);

  unsigned long nx[4] = {0, 0, 0, 0};
  double O[3] = {0.0, 0.0, 0.0};
  double A[3] = {1.0, 0.0, 0.0};
  double B[3] = {0.0, 1.0, 0.0};
  double C[3] = {0.0, 0.0, 1.0};

  file
    >> nx[0] >> nx[1] >> nx[2]
    >> O[0] >> O[1] >> O[2]
    >> A[0] >> A[1] >> A[2]
    >> B[0] >> B[1] >> B[2]
    >> C[0] >> C[1] >> C[2];

  nx[3] = nx[0]*nx[1];
  unsigned long gnx[3] = {nx[0]/2, nx[1]/2, nx[2]/2};

  double da[3];
  Math3D::sub(da,A,O);
  Math3D::scale(da,da,1.0/nx[0]);

  double db[3];
  Math3D::sub(db,B,O);
  Math3D::scale(db,db,1.0/nx[1]);

  double dc[3];
  Math3D::sub(dc,C,O);
  Math3D::scale(dc,dc,1.0/nx[2]);

  unsigned long n = nx[0]*nx[1]*nx[2];

  vtkPoints *points = vtkPoints::New();
  vtkDoubleArray *pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents(3);
  pts->SetNumberOfTuples(n);
  double *pPoints = pts->GetPointer(0);
  points->SetData(pts);
  pts->Delete();

  vtkDoubleArray *array = vtkDoubleArray::New();
  array->SetName(arrayName.c_str());
  array->SetNumberOfTuples(n);
  double *arrayData = array->GetPointer(0);

  unsigned long q = 0;
  for (unsigned long k=0; k<nx[2]; ++k)
    {
    unsigned long t = (k+gnx[2])%(nx[2]-1);
    for (unsigned long j=0; j<nx[1]; ++j)
      {
      unsigned long s =(j+gnx[1])%(nx[1]-1);
      for (unsigned long i=0; i<nx[0]; ++i,++q)
        {
        unsigned long r =(i+gnx[0])%(nx[0]-1);
        file >> arrayData[r + s*nx[0] + t*nx[3]];
        if (!file.good())
          {
          cerr
            << "ERROR: failed to read data at "
            << i << ", " << j << ", " << k << endl;
          return -11;
          }
        double pt[3] = {
            O[0]+i*da[0]+j*db[0]+k*dc[0],
            O[1]+i*da[1]+j*db[1]+k*dc[1],
            O[2]+i*da[2]+j*db[2]+k*dc[2]
            };

        unsigned long qq = 3*q;
        Math3D::copy(pPoints+qq, pt);
        }
      }
    }

  sgrid->SetDimensions(nx[0],nx[1],nx[2]);
  sgrid->SetPoints(points);
  sgrid->GetPointData()->AddArray(array);

  points->Delete();
  array->Delete();

  file >> header;
  if (!file.good() || (header != "END_DATAGRID_3D"))
    {
    cerr << "ERROR: END_DATAGRID_3D not found" << endl;
    return -11;
    }

  file >> header;
  if (!file.good() || (header != "END_BLOCK_DATAGRID_3D"))
    {
    cerr << "ERROR: END_BLOCK_DATAGRID_3D not found" << endl;
    return -12;
    }

  return 0;
}

};

namespace XSF
{
// ---------------------------------------------------------------------------
int Read(const string &fileName, vtkMolecule *molecule, vtkStructuredGrid *sgrid)
{
  ifstream file(fileName.c_str());
  if (!file.good())
    {
    cerr << "ERROR: failed to open " << fileName << endl;
    return -1;
    }

  string header;
  file >> header;
  if (!file.good() || (header != "CRYSTAL"))
    {
    cerr << "ERROR: CRYSTAL not found" << endl;
    return -2;
    }

  // read atoms
  vector<double> ax;
  vector<double> af;
  vector<unsigned short> an;
  vector<string> al;
  vector<int> ali;
  vector<bool> as;
  vector<bool> acs;
  vector<bool> gh;
  if (ReadXSFAtoms(file, ax, af, an))
    {
    cerr << "ERROR: reading atoms" << endl;
    return -3;
    }

  // construct molecule
  AtomicPropertiesBondDetector *detector = new AtomicPropertiesBondDetector;
  detector->SetDetectionMode(AtomicProperties::ATOMIC);
  detector->SetTolerance(1.05);

  al.resize(an.size());
  ali.resize(an.size());
  as.resize(an.size(), true);
  acs.resize(an.size(), false);
  gh.resize(an.size(), false);

  BuildMolecule(molecule, ax, an, al, ali, as, acs, gh, detector);
  delete detector;

  if (ReadXSFGrid(file, sgrid))
    {
    cerr << "ERROR: reading grid" << endl;
    return -5;
    }

  return 0;
}
};
