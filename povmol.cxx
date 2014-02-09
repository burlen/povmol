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
#include "vtkCamera.h"
#include "vtkCMLMoleculeReader.h"
#include "vtkMolecule.h"
#include "vtkMoleculeMapper2.h"
#include "vtkNew.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.cxx"

#include <string>
#include <fstream>
#include <iostream>

using namespace std;


namespace {

// Define interaction style
class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
  static KeyPressInteractorStyle* New();
  vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);

  KeyPressInteractorStyle() : Mapper(NULL), POVFile(NULL) {}

  // Description:
  // Set the mapper
  vtkSetObjectMacro(Mapper, vtkMoleculeMapper2);

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

    if (key == "w")
      {
      cerr << "wrote " << this->POVFile << endl;
      ofstream povf(this->POVFile);
      povf << this->Mapper->GetPOVRayStream() << endl;
      }
    else
    if (key == "a")
      {
      double r = this->Mapper->GetAtomicRadiusScaleFactor();
      cerr << "Enter atom radius scale factor (" << r << ") :";
      cin >> r;
      this->Mapper->SetAtomicRadiusScaleFactor(r);
      this->Mapper->Update();
      }
    else
    if (key == "b")
      {
      double t = this->Mapper->GetBondRadius();
      cerr << "Enter bond radius (" << t << ") :";
      cin >> t;
      this->Mapper->SetBondRadius(t);
      this->Mapper->Update();
      }
    else
    if (key == "s")
      {
      this->Mapper->SetBondColorModeToSingleColor();
      }
    else
    if (key == "d")
      {
      this->Mapper->SetBondColorModeToDiscreteByAtom();
      }
    else
    if (key == "p")
      {
      vtkPolyData *bonds = this->Mapper->GetBonds();
      vtkPolyDataWriter *w = vtkPolyDataWriter::New();
      w->SetInputData(bonds);
      w->SetFileName("bonds.vtk");
      w->Write();
      w->Delete();
      bonds->Delete();

      vtkPolyData *atoms = this->Mapper->GetAtoms();
      w = vtkPolyDataWriter::New();
      w->SetInputData(atoms);
      w->SetFileName("atoms.vtk");
      w->Write();
      w->Delete();
      atoms->Delete();
      }
    else
    if (key == "h")
      {
      cerr << "help:" << endl
        << "w - write POV file" << endl
        << "b - enter bond radius" << endl
        << "a - enter atom radius scale factor" << endl
        << "s - single color for bonds" << endl
        << "d - discrete color for bonds" << endl
        << "p - write poly data" << endl
        << endl;
      }

    // Forward events
    vtkInteractorStyleTrackballCamera::OnKeyPress();
    }

private:
  vtkMoleculeMapper2 *Mapper;
  char *POVFile;

};

// --------------------------------------------------------------------------
vtkStandardNewMacro(KeyPressInteractorStyle);

// **************************************************************************
void fixPath(string &path)
{
  size_t n = path.size();
  for (size_t i=0; i<n; ++i)
    {
    if (path[i] == '/')
      {
      path[i] = '\\';
      }
    }
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

  vtkNew<vtkCMLMoleculeReader> cmlSource;
  cmlSource->SetFileName(input.c_str());

  cmlSource->Update();
  if (cmlSource->GetOutput()->GetNumberOfAtoms() < 30)
    {
    return -2;
    }

  vtkNew<vtkMoleculeMapper2> molmapper;
  molmapper->SetInputConnection(cmlSource->GetOutputPort());
  molmapper->UseBallAndStickSettings();
  molmapper->SetAtomicRadiusScaleFactor(0.0075);
  molmapper->SetBondRadius(0.0035);
  molmapper->SetPOVRayStreaming(true);

  vtkNew<vtkActor> actor;
  actor->SetMapper(molmapper.GetPointer());

  vtkNew<vtkRenderer> ren;
  vtkNew<vtkRenderWindow> win;
  win->AddRenderer(ren.GetPointer());
  vtkNew<vtkRenderWindowInteractor> iren;
  iren->SetRenderWindow(win.GetPointer());

  vtkNew<KeyPressInteractorStyle> style;
  style->SetMapper(molmapper.GetPointer());
  style->SetPOVFile(output.c_str());

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
