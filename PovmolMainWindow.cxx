#include <PovmolMainWindow.h>

#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkLightCollection.h"
#include "vtkLight.h"
#include "vtkCamera.h"
#include "vtkCMLMoleculeReader.h"
#include "vtkCIFMoleculeReader.h"
#include "vtkPolyDataMapper.h"
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
#include "vtkAxesActor.h"
#include "vtkOrientationMarkerWidget.h"
#include "vtkCamera.h"

#include "Utility.h"
#include "Math3D.h"

#include <string>
#include <iostream>
#include <vector>
#include <cfloat>
#include <sstream>

using std::vector;
using std::cerr;
using std::endl;
using std::ostringstream;

#include <QFileDialog>
#include <QMessageBox>
#include <QAction>
#include <QDoubleValidator>
#include <QListWidgetItem>

#define PovmolMainWindowDEBUG

// --------------------------------------------------------------------------
PovmolMainWindow::PovmolMainWindow()
        :
    Renderer(0),
    Reader(0),
    MoleculeMapper(0)
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::PovmolMainWindow" << endl;
  #endif
  this->Ui = new Ui_PovmolMainWindowUi;
  this->Ui->setupUi(this);
  this->Ui->BondRadius->setValidator(new QDoubleValidator(this->Ui->BondRadius));
  this->Ui->AtomRadiusFactor->setValidator(new QDoubleValidator(this->Ui->AtomRadiusFactor));

  //his->Ui->ViewWidget->GetRenderWindow()->SetSize(

  this->Renderer = vtkRenderer::New();
  this->Ui->ViewWidget->GetRenderWindow()->AddRenderer(this->Renderer);
  this->Renderer->Delete();
  this->Renderer->SetBackground2(0.08,0.08,0.08);
  this->Renderer->SetBackground2(0.87,0.87,0.87);
  this->Renderer->GradientBackgroundOn();
  this->UpdateLightIntensity();

  vtkAxesActor *axes = vtkAxesActor::New();
  vtkOrientationMarkerWidget *axesWidget = vtkOrientationMarkerWidget::New();
  axesWidget->SetOutlineColor(0.9300, 0.5700, 0.1300);
  axesWidget->SetOrientationMarker(axes);
  axesWidget->SetInteractor(this->Ui->ViewWidget->GetInteractor());
  axesWidget->SetViewport(0.0, 0.0, 0.3, 0.3);
  axesWidget->SetEnabled(1);
  axesWidget->InteractiveOn();

  this->Ui->ViewWidget->GetRenderWindow()->Render();

  connect(this->Ui->OpenAction, SIGNAL(triggered()), this, SLOT(OpenFile()));
  connect(this->Ui->WritePOVAction, SIGNAL(triggered()), this, SLOT(WritePOV()));
  connect(this->Ui->WriteImageAction, SIGNAL(triggered()), this, SLOT(WriteImage()));
  connect(this->Ui->WriteGeometryAction, SIGNAL(triggered()), this, SLOT(WriteVTK()));
  connect(this->Ui->AtomRadiusFactor, SIGNAL(editingFinished()), this, SLOT(UpdateAtomRadiusFactor()));
  connect(this->Ui->BondRadius, SIGNAL(editingFinished()), this, SLOT(UpdateBondRadius()));
  connect(this->Ui->BondDetectionMode, SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateBondDetectionMode(int)));
  connect(this->Ui->BondProximityFactor, SIGNAL(editingFinished()), this, SLOT(UpdateBondProximityFactor()));
  connect(this->Ui->BondColorAtoms, SIGNAL(toggled(bool)), this, SLOT(UpdateBondColorMode()));
  connect(this->Ui->BondColorSingle, SIGNAL(toggled(bool)), this, SLOT(UpdateBondColorMode()));
  connect(this->Ui->LightIntensity, SIGNAL(valueChanged(int)), this, SLOT(UpdateLightIntensity()));
  connect(this->Ui->ActiveTransforms, SIGNAL(itemChanged(QListWidgetItem*)), this, SLOT(UpdateActiveTransforms()));
  connect(this->Ui->ViewDownX, SIGNAL(released()), this, SLOT(ViewDownX()));
  connect(this->Ui->ViewDownY, SIGNAL(released()), this, SLOT(ViewDownY()));
  connect(this->Ui->ViewDownZ, SIGNAL(released()), this, SLOT(ViewDownZ()));
  connect(this->Ui->ViewUpX, SIGNAL(released()), this, SLOT(ViewUpX()));
  connect(this->Ui->ViewUpY, SIGNAL(released()), this, SLOT(ViewUpY()));
  connect(this->Ui->ViewUpZ, SIGNAL(released()), this, SLOT(ViewUpZ()));
  
  //connect(this->Ui->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));
}

// --------------------------------------------------------------------------
PovmolMainWindow::~PovmolMainWindow()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::~PovmolMainWindow" << endl;
  #endif
  this->Initialize();
  //delete this->Ui;
}

// --------------------------------------------------------------------------
void PovmolMainWindow::OpenFile()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::OpenFile" << endl;
  #endif
  string fileName = QFileDialog::getOpenFileName(
        this, tr("Open"), "", tr("Molecule Files (*.cifpp *.XSF)")).toStdString();

  if (fileName != "")
    {
    this->Initialize();

    if (vtkCIFMoleculeReader::CanReadFile(fileName.c_str()))
      {
      vtkCIFMoleculeReader *cifReader = vtkCIFMoleculeReader::New();
      cifReader->SetFileName(fileName.c_str());
      if (!cifReader->GetInitialized())
        {
        ostringstream oss;
        oss << "CIF Reader failed to open " << fileName;
        QMessageBox::warning(this, tr("Error"), tr(oss.str().c_str()));
        return;
        }
      this->Reader = cifReader;
      }
    else
      {
      ostringstream oss;
      oss << "No this->Reader could open " << fileName;
      QMessageBox::warning(this, tr("Error"), tr(oss.str().c_str()));
      return;
      }

    this->BuildPipeline();
    this->ShowTransforms();
    }
}

// --------------------------------------------------------------------------
void PovmolMainWindow::ShowTransforms()
{
  size_t n = this->Reader->GetNumberOfTransforms();
  for (size_t i=0; i<n; ++i)
    {
    QListWidgetItem *item = new QListWidgetItem;
    item->setData(Qt::DisplayRole, this->Reader->GetTransformLabel(i));
    item->setData(Qt::CheckStateRole, Qt::Checked);
    this->Ui->ActiveTransforms->addItem(item);
    }
}

// --------------------------------------------------------------------------
void PovmolMainWindow::UpdateActiveTransforms()
{
  int n = this->Ui->ActiveTransforms->count();
  for (int i=0; i<n; ++i)
    {
    QListWidgetItem *item = this->Ui->ActiveTransforms->item(i);
    if (item->checkState() == Qt::Checked)
      {
      this->Reader->ActivateTransform(i);
      }
    else
      {
      this->Reader->DeactivateTransform(i);
      }
    }
  this->Render();
}

// --------------------------------------------------------------------------
void PovmolMainWindow::Initialize()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::Initialize" << endl;
  #endif
  if (this->Reader)
    {
    this->Reader->Delete();
    this->Reader = NULL;
    }

  if (this->MoleculeMapper)
    {
    this->MoleculeMapper->Delete();
    this->MoleculeMapper = NULL;
    }

  this->Ui->ActiveTransforms->clear();
}

// --------------------------------------------------------------------------
void PovmolMainWindow::BuildPipeline()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::BuildPipeline" << endl;
  #endif
  this->Reader->SetBondDetectionMode(this->Ui->BondDetectionMode->currentIndex());


  this->MoleculeMapper = vtkMoleculeMapper2::New();
  this->MoleculeMapper->SetInputConnection(this->Reader->GetOutputPort(0));
  this->MoleculeMapper->UseBallAndStickSettings();
  this->MoleculeMapper->SetAtomicRadiusScaleFactor(this->Ui->AtomRadiusFactor->text().toDouble());
  this->MoleculeMapper->SetBondRadius(this->Ui->BondRadius->text().toDouble());
  this->MoleculeMapper->SetPOVRayStreaming(true);

  vtkActor *molActor = vtkActor::New();
  molActor->GetProperty()->SetAmbient(0.1);
  molActor->GetProperty()->SetSpecular(0.5);
  molActor->GetProperty()->SetDiffuse(0.2);
  molActor->SetMapper(this->MoleculeMapper);

  this->Renderer->AddActor(molActor);

  vtkPolyDataMapper *axesMapper = vtkPolyDataMapper::New();
  axesMapper->SetInputConnection(this->Reader->GetOutputPort(1));
  axesMapper->SetScalarModeToUseCellData();
  axesMapper->ScalarVisibilityOn();

  vtkActor *axesActor = vtkActor::New();
  axesActor->SetMapper(axesMapper);
  axesMapper->Delete();

  this->Renderer->AddActor(axesActor);
  axesActor->Delete();

  this->Renderer->Render();
  this->Renderer->ResetCamera();
  this->Renderer->Render();
}

// --------------------------------------------------------------------------
void PovmolMainWindow::WriteVTK()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::WriteVTK" << endl;
  #endif
  string fileName = QFileDialog::getSaveFileName(
    this, tr("Write Geometry"), "", tr("VTK Legacy (*.vtk)")).toStdString();

  if (fileName.c_str())
    {
    vtkMoleculeToBondStickFilter *bf = vtkMoleculeToBondStickFilter::New();
    bf->SetInputData(this->MoleculeMapper->GetInput());

    //vtkPolyData *bonds = this->MoleculeMapper->GetBonds();
    vtkPolyDataWriter *w = vtkPolyDataWriter::New();
    //w->SetInputData(bonds);
    w->SetInputConnection(bf->GetOutputPort());
    w->SetFileName("bonds.vtk");
    w->Write();
    w->Delete();
    //bonds->Delete();
    bf->Delete();
    cerr << "wrote bonds to bonds.vtk" << endl;

    vtkPolyData *atoms = this->MoleculeMapper->GetAtoms();
    w = vtkPolyDataWriter::New();
    w->SetInputData(atoms);
    w->SetFileTypeToBinary();
    w->SetFileName("atoms.vtk");
    w->Write();
    w->Delete();
    atoms->Delete();
    cerr << "wrote atoms to atoms.vtk" << endl;
    }
}

// --------------------------------------------------------------------------
void PovmolMainWindow::WritePOV()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::WritePOV" << endl;
  #endif
  string fileName = QFileDialog::getSaveFileName(
    this, tr("Export"), "", tr("POVRay (*.pov)")).toStdString();

  if (fileName != "")
    {
    ofstream povf(fileName.c_str());
    povf << this->MoleculeMapper->GetPOVRayStream() << endl;
    cerr << "wrote " << fileName << endl;
    }
}

// --------------------------------------------------------------------------
void PovmolMainWindow::WriteImage()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::WriteImage" << endl;
  #endif
  string fileName = QFileDialog::getSaveFileName(
    this, tr("Save Screenshot"), "", tr("Images (*.png)")).toStdString();

  if (fileName != "")
    {
    vtkWindowToImageFilter *i=vtkWindowToImageFilter::New();
    i->SetInput(this->Ui->ViewWidget->GetRenderWindow());
    vtkPNGWriter *w = vtkPNGWriter::New();
    w->SetFileName(fileName.c_str());
    w->SetInputConnection(i->GetOutputPort());
    w->Write();
    w->Delete();
    i->Delete();
    }
}

// --------------------------------------------------------------------------
void PovmolMainWindow::UpdateBondDetectionMode(int id)
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::UpdateBondDetectionMode" << endl;
  #endif
  if (this->Reader)
    {
    this->Reader->SetBondDetectionMode(id);
    this->Render();
    }
}

// --------------------------------------------------------------------------
void PovmolMainWindow::UpdateAtomRadiusFactor()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::UpdateAtomRadiusFactor" << endl;
  #endif
  if (this->MoleculeMapper)
    {
    this->MoleculeMapper->SetAtomicRadiusScaleFactor(this->Ui->AtomRadiusFactor->text().toDouble());
    this->Render();
    }
}

// --------------------------------------------------------------------------
void PovmolMainWindow::UpdateBondRadius()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::UpdateBondRadius" << endl;
  #endif
  if (this->MoleculeMapper)
    {
    this->MoleculeMapper->SetBondRadius(this->Ui->BondRadius->text().toDouble());
    this->Render();
    }

}

// --------------------------------------------------------------------------
void PovmolMainWindow::UpdateBondProximityFactor()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::UpdateBondProximityFactor" << endl;
  #endif
  if (this->MoleculeMapper)
    {
    this->Reader->SetBondProximityFactor(this->Ui->BondProximityFactor->text().toDouble());
    this->Render();
    }

}

// --------------------------------------------------------------------------
void PovmolMainWindow::UpdateBondColorMode()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::UpdateBondColorMode" << endl;
  #endif
  if (this->Ui->BondColorSingle->isChecked())
    {
    this->MoleculeMapper->SetBondColorModeToSingleColor();
    }
  else
    {
    this->MoleculeMapper->SetBondColorModeToDiscreteByAtom();
    }
  this->Render();
}

// --------------------------------------------------------------------------
void PovmolMainWindow::UpdateLightIntensity()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::UpdateLightIntensity" << endl;
  #endif
  double intensity
    = static_cast<double>(this->Ui->LightIntensity->value()) /
         static_cast<double>(this->Ui->LightIntensity->maximum());

  vtkLightCollection *lc = this->Renderer->GetLights();
  lc->InitTraversal();
  vtkLight *l;
  while (l = lc->GetNextItem())
    {
    l->SetIntensity(intensity);
    }

  this->Render();
}

// --------------------------------------------------------------------------
void PovmolMainWindow::Render()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::Render" << endl;
  #endif
  this->Ui->ViewWidget->GetRenderWindow()->Render();
}

// --------------------------------------------------------------------------
void PovmolMainWindow::ViewDownX()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::ViewDownX" << endl;
  #endif
  vtkCamera *cam = this->Renderer->GetActiveCamera();
  double d = Math3D::length(cam->GetPosition());
  double pos[3] = {d, 0.0, 0.0};
  Math3D::add(pos, pos, cam->GetFocalPoint());
  cam->SetPosition(pos);
  Renderer->ResetCamera();
  this->Render();
}

// --------------------------------------------------------------------------
void PovmolMainWindow::ViewUpX()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::ViewUpX" << endl;
  #endif
  vtkCamera *cam = this->Renderer->GetActiveCamera();
  double d = Math3D::length(cam->GetPosition());
  double pos[3] = {-d, 0.0, 0.0};
  Math3D::add(pos, pos, cam->GetFocalPoint());
  cam->SetPosition(pos);
  cam->SetViewUp(0.0, 1.0, 0.0);
  this->Render();
}

// --------------------------------------------------------------------------
void PovmolMainWindow::ViewDownY()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::ViewDownY" << endl;
  #endif
  vtkCamera *cam = this->Renderer->GetActiveCamera();
  double d = Math3D::length(cam->GetPosition());
  double pos[3] = {0.0, d, 0.0};
  Math3D::add(pos, pos, cam->GetFocalPoint());
  cam->SetPosition(pos);
  cam->SetViewUp(0.0, 0.0, 1.0);
  this->Render();
}

// --------------------------------------------------------------------------
void PovmolMainWindow::ViewUpY()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::ViewUpY" << endl;
  #endif
  vtkCamera *cam = this->Renderer->GetActiveCamera();
  double d = Math3D::length(cam->GetPosition());
  double pos[3] = {0.0, -d, 0.0};
  Math3D::add(pos, pos, cam->GetFocalPoint());
  cam->SetPosition(pos);
  cam->SetViewUp(0.0, 0.0, 1.0);
  this->Render();
}

// --------------------------------------------------------------------------
void PovmolMainWindow::ViewDownZ()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::ViewDownZ" << endl;
  #endif
  vtkCamera *cam = this->Renderer->GetActiveCamera();
  double d = Math3D::length(cam->GetPosition());
  double pos[3] = {0.0, 0.0, d};
  Math3D::add(pos, pos, cam->GetFocalPoint());
  cam->SetPosition(pos);
  cam->SetViewUp(0.0, 1.0, 0.0);
  this->Render();
}

// --------------------------------------------------------------------------
void PovmolMainWindow::ViewUpZ()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::ViewUpZ" << endl;
  #endif
  vtkCamera *cam = this->Renderer->GetActiveCamera();
  double d = Math3D::length(cam->GetPosition());
  double pos[3] = {0.0, 0.0, -d};
  Math3D::add(pos, pos, cam->GetFocalPoint());
  cam->SetPosition(pos);
  cam->SetViewUp(0.0, 1.0, 0.0);
  this->Render();
}

  /*
// --------------------------------------------------------------------------
void PovmolMainWindow::closeEvent(QCloseEvent *event)
{
  event->ignore();
  if (maybeSave())
    {
    writeSettings();
    event->accept();
    }
  else
    {
    event->ignore();
    }
}
  */
