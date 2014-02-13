/*=========================================================================

  Program:   Visualization Toolkit

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkLightCollection.h"
#include "vtkLight.h"
#include "vtkCamera.h"
#include "vtkCMLMoleculeReader.h"
#include "vtkMolecule.h"
#include "vtkMoleculeMapper2.h"
#include "vtkMoleculeToBondStickFilter.h"
#include "vtkNew.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkDataSetWriter.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkStructuredGrid.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <cfloat>

using namespace std;

namespace {

void initbounds(double *b)
{
  b[0] =  DBL_MAX;
  b[1] = -DBL_MAX;
  b[2] =  DBL_MAX;
  b[3] = -DBL_MAX;
  b[4] =  DBL_MAX;
  b[5] = -DBL_MAX;
}

// ---------------------------------------------------------------------------
void bounds(double *b, double *a)
{
  b[0] = min(b[0],a[0]);
  b[1] = max(b[1],a[0]);
  b[2] = min(b[2],a[1]);
  b[3] = max(b[3],a[1]);
  b[4] = min(b[4],a[2]);
  b[5] = max(b[5],a[2]);
}

// ---------------------------------------------------------------------------
double frac(double a, double b, double f)
{
  return (b-a)*f+a;
}

// ---------------------------------------------------------------------------
double *abs(double *b, double *a)
{
  b[0] = fabs(a[0]);
  b[1] = fabs(a[1]);
  b[2] = fabs(a[2]);
  return b;
}

// ---------------------------------------------------------------------------
double *scale(double *c, double *a, double b)
{
  c[0] = a[0] * b;
  c[1] = a[1] * b;
  c[2] = a[2] * b;
  return c;
}

// ---------------------------------------------------------------------------
double *sub(double *c, double *b, double *a)
{
  c[0] = b[0] - a[0];
  c[1] = b[1] - a[1];
  c[2] = b[2] - a[2];
  return c;
}

// ---------------------------------------------------------------------------
double *add(double *c, double *a, double *b)
{
  c[0] = b[0] + a[0];
  c[1] = b[1] + a[1];
  c[2] = b[2] + a[2];
  return c;
}

// ---------------------------------------------------------------------------
double dot(double *a, double *b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// ---------------------------------------------------------------------------
double length(double *a, double *b)
{
    double r[3];
    sub(r,b,a);
    return sqrt(dot(r,r));
}

// ---------------------------------------------------------------------------
double length(double *a)
{
    return sqrt(dot(a,a));
}

// ---------------------------------------------------------------------------
double *copy(double *a, double *b)
{
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
  return a;
}

// ---------------------------------------------------------------------------
double *normalize(double *b, double *a)
{
  double r = length(a);
  b[0] = a[0]/r;
  b[1] = a[1]/r;
  b[2] = a[2]/r;
  return b;
}

// ---------------------------------------------------------------------------
void defplane(double *abcd, double *n, double *r0)
{
  double nn[3];
  copy(abcd,normalize(nn,n));
  abcd[3] = -dot(nn,r0);
}

// ---------------------------------------------------------------------------
double evalplane(double *abcd, double *xyz)
{
  return dot(abcd,xyz)+abcd[3];
}

/*
  iabr[1][1] = 0.0;
  iabr[8][8] = 0.0;
  iabr[49][49] = 0.0;
  iabr[1][8] = 1.025770337248545;
  iabr[8][1] = 1.025770337248545;
  iabr[1][49] = 1.739025683089;
  iabr[49][1] = 1.739025683089;
  iabr[8][49] = 2.165492231328947;
  iabr[49][8] = 2.165492231328947;
*/

enum {
    ATOMIC,
    IONIC,
    COVALENT,
    VANDERWAALS,
    CRYSTAL
    };

double atomradius(unsigned short an, int mode=ATOMIC)
{
// Nu  Elem  Atomic R  Ionic R  Covalent R  Van-der-Waals R  "Crystal" R
// 1   H     0.53      0.25     0.37        1.20             0.10
// 7   N     0.56      0.65     0.75        1.55             0.30
// 49  In    1.56      1.55     1.44        1.93             0.94
  static map<int,double> radius;
  static bool init = false;
  if (!init)
    {
    init = true;
    switch(mode)
      {
      case ATOMIC:
        radius[1]  = 0.53;
        radius[7]  = 0.56;
        radius[49] = 1.56;
        break;
      case IONIC:
        radius[1]  = 0.25;
        radius[7]  = 0.65;
        radius[49] = 1.55;
        break;
      case COVALENT:
        radius[1]  = 0.37;
        radius[7]  = 0.75;
        radius[49] = 1.44;
        break;
      case VANDERWAALS:
        radius[1]  = 1.20;
        radius[7]  = 1.55;
        radius[49] = 1.93;
        break;
      case CRYSTAL:
        radius[1]  = 0.10;
        radius[7]  = 0.30;
        radius[49] = 0.94;
      }
    }
  return radius[an];
}


// ---------------------------------------------------------------------------
int BuildVTKMolecule(
    vtkMolecule *molecule,
    vector<double> &ax,
    vector<unsigned short> &an)
{

  molecule->Initialize();
  // atoms
  double bondr = 0.0;
  unsigned long natoms = an.size();
  for (unsigned long i=0; i<natoms; ++i)
    {
    if (!an[i]) continue;
    unsigned long ii = 3*i;
    molecule->AppendAtom(an[i], ax[ii], ax[ii+1], ax[ii+2]);
    /*else
      {
      bondr = length(&ax[ii]);
      }*/
    }
  // bonds
  for (unsigned int i=0; i<natoms; ++i)
    {
    if (!an[i]) continue;
    unsigned int ii = 3*i;
    vector<unsigned long> bonds;
    for (unsigned int j=i+1; j<natoms; ++j)
      {
      if (!an[j]) continue;
      unsigned int jj = 3*j;
      double r[3];
      sub(r,&ax[jj],&ax[ii]);
      if (length(r) < 1.05*(atomradius(an[i])+atomradius(an[j])))
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
  return 0;
}

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
  scale(a,A,0.5);

  double b[3];
  scale(b,B,0.5);

  double r0[3];
  add(r0,a,b);

  double alpha[4];
  defplane(alpha,A,r0);

  double beta[4];
  defplane(beta,B,r0);


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
  double R = 0.0;
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
    if ((evalplane(alpha,pax) < 0.0) && (evalplane(beta,pax) < 0.0))
      {
      add(pax,pax,a);
      add(pax,pax,b);
      }
    else
    if ((evalplane(alpha,pax) >= 0.0) && (evalplane(beta,pax) >= 0.0))
      {
      sub(pax,pax,a);
      sub(pax,pax,b);
      }
    else
    if ((evalplane(alpha,pax) < 0.0) && (evalplane(beta,pax) >= 0.0))
      {
      add(pax,pax,a);
      sub(pax,pax,b);
      }
    else
    if ((evalplane(alpha,pax) >= 0.0) && (evalplane(beta,pax) < 0.0))
      {
      sub(pax,pax,a);
      add(pax,pax,b);
      }
    else
      {
      cerr << "ERROR: not in any case" << endl;
      }
    }
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

  unsigned long nx[4] = {0, 0, 0, 0.0};
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
  sub(da,A,O);
  scale(da,da,1.0/nx[0]);

  double db[3];
  sub(db,B,O);
  scale(db,db,1.0/nx[1]);

  double dc[3];
  sub(dc,C,O);
  scale(dc,dc,1.0/nx[2]);

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
        copy(pPoints+qq, pt);
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

// ---------------------------------------------------------------------------
int ReadXSFFile(const string &fileName, vtkMolecule *molecule, vtkStructuredGrid *sgrid)
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
  if (ReadXSFAtoms(file, ax, af, an))
    {
    cerr << "ERROR: reading atoms" << endl;
    return -3;
    }

  // construct molecule
  if (BuildVTKMolecule(molecule,ax,an))
    {
    cerr << "ERROR: constructing vtk molecule" << endl;
    return -4;
    }

  if (ReadXSFGrid(file, sgrid))
    {
    cerr << "ERROR: reading grid" << endl;
    return -5;
    }

  return 0;
}

// Define interaction style
class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
  static KeyPressInteractorStyle* New();
  vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);

  KeyPressInteractorStyle() :
    Molecule(NULL),
    Grid(NULL),
    Mapper(NULL),
    Renderer(NULL),
    RenderWindow(NULL),
    POVFile(NULL)
  {}

  ~KeyPressInteractorStyle()
    {
    this->SetMolecule(NULL);
    this->SetGrid(NULL);
    this->SetMapper(NULL);
    this->SetRenderer(NULL);
    this->SetRenderWindow(NULL);
    }
  // Description:
  // Set the objects
  vtkSetObjectMacro(Molecule, vtkMolecule);
  vtkSetObjectMacro(Grid, vtkStructuredGrid);
  vtkSetObjectMacro(Mapper, vtkMoleculeMapper2);
  vtkSetObjectMacro(Renderer, vtkRenderer);
  vtkSetObjectMacro(RenderWindow, vtkRenderWindow);

  // Description:
  // Set the output file name
  vtkSetStringMacro(POVFile);

  // Description:
  // Hanlde key press events
  virtual void OnKeyPress()
    {
    // Get the keypress
    vtkRenderWindowInteractor *rwi = this->Interactor;
    std::string key = rwi->GetKeySym();

    if (key == "W")
      {
      ofstream povf(this->POVFile);
      povf << this->Mapper->GetPOVRayStream() << endl;
      cerr << "wrote " << this->POVFile << endl;

      vtkMoleculeToBondStickFilter *bf = vtkMoleculeToBondStickFilter::New();
      bf->SetInputData(this->Molecule);


      //vtkPolyData *bonds = this->Mapper->GetBonds();
      vtkPolyDataWriter *w = vtkPolyDataWriter::New();
      //w->SetInputData(bonds);
      w->SetInputConnection(bf->GetOutputPort());
      w->SetFileName("bonds.vtk");
      w->Write();
      w->Delete();
      //bonds->Delete();
      bf->Delete();
      cerr << "wrote bonds to bonds.vtk" << endl;

      vtkPolyData *atoms = this->Mapper->GetAtoms();
      w = vtkPolyDataWriter::New();
      w->SetInputData(atoms);
      w->SetFileTypeToBinary();
      w->SetFileName("atoms.vtk");
      w->Write();
      w->Delete();
      atoms->Delete();
      cerr << "wrote atoms to atoms.vtk" << endl;

      vtkDataSetWriter *sw = vtkDataSetWriter::New();
      sw->SetInputData(this->Grid);
      sw->SetFileName("grid.vtk");
      sw->Write();
      sw->Delete();
      cerr << "wrote grid.vtk" << endl;

      return;
      }
    else
    if (key == "A")
      {
      double r = this->Mapper->GetAtomicRadiusScaleFactor();
      cerr << "Enter atom radius scale factor (" << r << ") :";
      cin >> r;
      this->Mapper->SetAtomicRadiusScaleFactor(r);
      this->Mapper->Update();
      return;
      }
    else
    if (key == "B")
      {
      double t = this->Mapper->GetBondRadius();
      cerr << "Enter bond radius (" << t << ") :";
      cin >> t;
      this->Mapper->SetBondRadius(t);
      this->Mapper->Update();
      return;
      }
    else
    if (key == "S")
      {
      this->Mapper->SetBondColorModeToSingleColor();
      cerr << "set bonds to single color" << endl;
      return;
      }
    else
    if (key == "D")
      {
      this->Mapper->SetBondColorModeToDiscreteByAtom();
      cerr << "set bonds to discrete color" << endl;
      return;
      }
    else
    if (key == "I")
      {
      string file;
      cerr << "Enter file name :";
      cin >> file;
      vtkWindowToImageFilter *i=vtkWindowToImageFilter::New();
      i->SetInput(this->RenderWindow);
      vtkPNGWriter *w=vtkPNGWriter::New();
      w->SetFileName(file.c_str());
      w->SetInputConnection(i->GetOutputPort());
      w->Write();
      w->Delete();
      i->Delete();
      }
    else
    if (key == "L")
      {
      cerr << "Enter light intensity :" << endl;
      double i = 1.0;
      cin >> i;
      vtkLightCollection *lc = this->Renderer->GetLights();
      lc->InitTraversal();
      vtkLight *l;
      while (l = lc->GetNextItem())
        {
        l->SetIntensity(i);
        }
      }
    else
    if (key == "h")
      {
      cerr << "help:" << endl
        << "W - write POV file" << endl
        << "B - enter bond radius" << endl
        << "A - enter atom radius scale factor" << endl
        << "S - single color for bonds" << endl
        << "D - discrete color for bonds" << endl
        << endl;
      return;
      }

    // Forward other events
    vtkInteractorStyleTrackballCamera::OnKeyPress();
    }

private:
  vtkMolecule *Molecule;
  vtkStructuredGrid *Grid;
  vtkMoleculeMapper2 *Mapper;
  vtkRenderer *Renderer;
  vtkRenderWindow *RenderWindow;
  char *POVFile;

};

// --------------------------------------------------------------------------
vtkStandardNewMacro(KeyPressInteractorStyle);

// --------------------------------------------------------------------------
void fixPath(string &path)
{
#ifdef _WIN32
  size_t n = path.size();
  for (size_t i=0; i<n; ++i)
    {
    if (path[i] == '/')
      {
      path[i] = '\\';
      }
    }
#endif
}

};

// **************************************************************************
int main(int argc, char *argv[])
{
  if (argc != 3)
    {
    cerr << "Error usage: povmol input.cml output.pov" << endl;
    return -1;
    }
  string input(argv[1]);
  fixPath(input);

  string output(argv[2]);
  fixPath(output);

  /*
  vtkNew<vtkCMLMoleculeReader> cmlSource;
  cmlSource->SetFileName(input.c_str());
  cmlSource->Update();
  */
  vtkMolecule *molecule = vtkMolecule::New();
  vtkStructuredGrid *grid = vtkStructuredGrid::New();
  if (ReadXSFFile(input, molecule, grid))
    {
    cerr << "ERROR: I/O error" << endl;
    return -1;
    }
  cerr << "read " << input << endl;

  /* while looking for interesting molecules ...
  if (cmlSource->GetOutput()->GetNumberOfAtoms() < 30)
    {
    return -2;
    }*/

  vtkNew<vtkMoleculeMapper2> molmapper;
  molmapper->SetInputData(molecule);
  molecule->Delete();
  //molmapper->SetInputConnection(cmlSource->GetOutputPort());
  molmapper->UseBallAndStickSettings();
  molmapper->SetAtomicRadiusScaleFactor(0.25);
  molmapper->SetBondRadius(0.15);
  molmapper->SetPOVRayStreaming(true);

  vtkNew<vtkActor> actor;
  actor->GetProperty()->SetAmbient(0.1);
  actor->GetProperty()->SetSpecular(0.5);
  actor->GetProperty()->SetDiffuse(0.2);
  actor->SetMapper(molmapper.GetPointer());

  vtkNew<vtkRenderer> ren;
  ren->SetBackground2(0.08,0.08,0.08);
  ren->SetBackground2(0.87,0.87,0.87);
  ren->GradientBackgroundOn();
  vtkNew<vtkRenderWindow> win;
  win->AddRenderer(ren.GetPointer());
  vtkNew<vtkRenderWindowInteractor> iren;
  iren->SetRenderWindow(win.GetPointer());

  vtkNew<KeyPressInteractorStyle> style;
  style->SetMolecule(molecule);
  style->SetGrid(grid);
  style->SetMapper(molmapper.GetPointer());
  style->SetPOVFile(output.c_str());
  style->SetRenderer(ren.GetPointer());
  style->SetRenderWindow(win.GetPointer());
  grid->Delete();

  iren->SetInteractorStyle(style.GetPointer());
  style->SetCurrentRenderer(ren.GetPointer());

  ren->AddActor(actor.GetPointer());

  ren->SetBackground(0.0,0.0,0.0);
  win->SetSize(450,450);
  win->Render();

  win->GetInteractor()->Initialize();
  win->GetInteractor()->Start();

  return 0;
}
