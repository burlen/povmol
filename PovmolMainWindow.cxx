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
#include "vtkLookupTable.h"
#include "vtkPOVExporter.h"

#include "Utility.h"
#include "Math3D.h"
#include "AtomicPropertiesBondDetector.h"
#include "TableBasedBondDetector.h"

#include <string>
#include <iostream>
#include <vector>
#include <cfloat>
#include <sstream>

using std::vector;
using std::cerr;
using std::endl;
using std::ostringstream;
using std::shared_ptr;
using std::make_shared;
using std::static_pointer_cast;

#include <QFileDialog>
#include <QColorDialog>
#include <QMessageBox>
#include <QAction>
#include <QDoubleValidator>
#include <QListWidgetItem>
#include <QColor>
#include <QPixmap>
#include <QIcon>

#include "PovmolBondTableDialog.h"
#include "PovmolLookupTableDialog.h"

#define PovmolMainWindowDEBUG

namespace {
// **************************************************************************
QIcon makeColorSwatch(QColor color)
{
  QPixmap patch(24, 24);
  patch.fill(color);
  QIcon swatch(patch);
  return swatch;
}
};

// --------------------------------------------------------------------------
PovmolMainWindow::PovmolMainWindow()
        :
    Renderer(0),
    Reader(0),
    MoleculeMapper(0),
    CoordinationSiteActor(0),
    TableDetector(make_shared<TableBasedBondDetector>()),
    PropertiesDetector(make_shared<AtomicPropertiesBondDetector>()),
    DepthPeelingEnabled(false)
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::PovmolMainWindow" << endl;
  #endif
  this->Ui = new Ui_PovmolMainWindowUi;
  this->Ui->setupUi(this);

  this->Ui->BondRadius->setValidator(new QDoubleValidator(this->Ui->BondRadius));
  this->Ui->AtomRadiusFactor->setValidator(new QDoubleValidator(this->Ui->AtomRadiusFactor));

  this->Ui->DuplicateAMinus->setValidator(new QIntValidator(this->Ui->DuplicateAMinus));
  this->Ui->DuplicateAPlus->setValidator(new QIntValidator(this->Ui->DuplicateAPlus));
  this->Ui->DuplicateBMinus->setValidator(new QIntValidator(this->Ui->DuplicateBMinus));
  this->Ui->DuplicateBPlus->setValidator(new QIntValidator(this->Ui->DuplicateBPlus));
  this->Ui->DuplicateCMinus->setValidator(new QIntValidator(this->Ui->DuplicateCMinus));
  this->Ui->DuplicateCPlus->setValidator(new QIntValidator(this->Ui->DuplicateCPlus));

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
  connect(this->Ui->EditBondTableAction, SIGNAL(triggered()), this, SLOT(EditBondTable()));
  connect(this->Ui->EditLookupTableAction, SIGNAL(triggered()), this, SLOT(EditLookupTable()));
  connect(this->Ui->AtomRadiusFactor, SIGNAL(editingFinished()), this, SLOT(UpdateAtomRadiusFactor()));
  connect(this->Ui->BondRadius, SIGNAL(editingFinished()), this, SLOT(UpdateBondRadius()));
  connect(this->Ui->BondDetectionMode, SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateBondDetectionMode(int)));
  connect(this->Ui->BondProximityFactor, SIGNAL(editingFinished()), this, SLOT(UpdateBondProximityFactor()));
  connect(this->Ui->BondColorAtoms, SIGNAL(toggled(bool)), this, SLOT(UpdateBondColorMode()));
  connect(this->Ui->BondColorSingle, SIGNAL(toggled(bool)), this, SLOT(UpdateBondColorMode()));
  connect(this->Ui->BondColor, SIGNAL(released()), this, SLOT(UpdateBondColor()));
  connect(this->Ui->LightIntensity, SIGNAL(valueChanged(int)), this, SLOT(UpdateLightIntensity()));
  connect(this->Ui->ActiveTransforms, SIGNAL(itemChanged(QListWidgetItem*)), this, SLOT(UpdateActiveTransforms()));
  connect(this->Ui->ActiveSites, SIGNAL(itemChanged(QListWidgetItem*)), this, SLOT(UpdateSites()));
  connect(this->Ui->CoordinationSites, SIGNAL(toggled(bool)), this, SLOT(EnableCoordinationSites(bool)));
  connect(this->Ui->ActiveCoordinationSites, SIGNAL(itemChanged(QListWidgetItem*)), this, SLOT(UpdateCoordinationSites()));
  connect(this->Ui->ActiveCoordinationSites, SIGNAL(itemDoubleClicked(QListWidgetItem*)), this, SLOT(UpdateCoordinationSiteColor(QListWidgetItem*)));
  connect(this->Ui->DuplicateAMinus, SIGNAL(editingFinished()), this, SLOT(UpdateDuplicates()));
  connect(this->Ui->DuplicateBMinus, SIGNAL(editingFinished()), this, SLOT(UpdateDuplicates()));
  connect(this->Ui->DuplicateCMinus, SIGNAL(editingFinished()), this, SLOT(UpdateDuplicates()));
  connect(this->Ui->DuplicateAPlus, SIGNAL(editingFinished()), this, SLOT(UpdateDuplicates()));
  connect(this->Ui->DuplicateBPlus, SIGNAL(editingFinished()), this, SLOT(UpdateDuplicates()));
  connect(this->Ui->DuplicateCPlus, SIGNAL(editingFinished()), this, SLOT(UpdateDuplicates()));
  connect(this->Ui->GhostBonds, SIGNAL(toggled(bool)), this, SLOT(UpdateGhostBonds(bool)));
  connect(this->Ui->ViewDownX, SIGNAL(released()), this, SLOT(ViewDownX()));
  connect(this->Ui->ViewDownY, SIGNAL(released()), this, SLOT(ViewDownY()));
  connect(this->Ui->ViewDownZ, SIGNAL(released()), this, SLOT(ViewDownZ()));
  connect(this->Ui->ViewUpX, SIGNAL(released()), this, SLOT(ViewUpX()));
  connect(this->Ui->ViewUpY, SIGNAL(released()), this, SLOT(ViewUpY()));
  connect(this->Ui->ViewUpZ, SIGNAL(released()), this, SLOT(ViewUpZ()));
  connect(this->Ui->CoordinationSiteProperties, SIGNAL(Translucent(bool)), this, SLOT(EnableDepthPeeling(bool)));
  connect(this->Ui->CoordinationSiteProperties, SIGNAL(Modified()), this, SLOT(Render()));
  connect(this->Ui->SiteProperties, SIGNAL(Modified()), this, SLOT(Render()));

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

    this->UpdateBondDetectionMode(this->Ui->BondDetectionMode->currentIndex());
    this->BuildPipeline();
    this->ShowTransforms();
    this->ShowSites();
    this->ShowCoordinationSites();
    }
}

// ----------------------------------------------------------------------------
void PovmolMainWindow::EnableDepthPeeling(bool enabled)
{
  if (this->DepthPeelingEnabled == enabled) return;
  this->DepthPeelingEnabled = enabled;
  // 1. Use a render window with alpha bits (as initial value is 0 (false)):
  // 2. Force to not pick a framebuffer with a multisample buffer
  // 3. Choose to use depth peeling (if supported) (initial value is 0 (false)):
  // 4. Set depth peeling parameters
  vtkRenderWindow *renderWindow = this->Ui->ViewWidget->GetRenderWindow();
  renderWindow->SetAlphaBitPlanes(enabled);
  renderWindow->SetMultiSamples(enabled?0:8);

  vtkRenderer *renderer = this->Renderer;
  renderer->SetUseDepthPeeling(enabled);
  renderer->SetMaximumNumberOfPeels(200);
  renderer->SetOcclusionRatio(0.1);
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
void PovmolMainWindow::ShowSites()
{
  size_t n = this->Reader->GetNumberOfSites();
  for (size_t i=0; i<n; ++i)
    {
    QListWidgetItem *item = new QListWidgetItem;
    item->setData(Qt::DisplayRole, this->Reader->GetSiteLabel(i));
    item->setData(Qt::CheckStateRole, Qt::Checked);
    this->Ui->ActiveSites->addItem(item);
    }
}

// --------------------------------------------------------------------------
void PovmolMainWindow::UpdateSites()
{
  if (this->Reader)
    {
    int n = this->Ui->ActiveSites->count();
    for (int i=0; i<n; ++i)
      {
      QListWidgetItem *item = this->Ui->ActiveSites->item(i);
      if (item->checkState() == Qt::Checked)
        {
        this->Reader->ActivateSite(i);
        }
      else
        {
        this->Reader->DeactivateSite(i);
        }
      }
    this->Render();
    }
}

// --------------------------------------------------------------------------
void PovmolMainWindow::ShowCoordinationSites()
{
  size_t n = this->Reader->GetNumberOfSites();
  for (size_t i=0; i<n; ++i)
    {
    QListWidgetItem *item = new QListWidgetItem;
    item->setData(Qt::DisplayRole, this->Reader->GetSiteLabel(i));
    item->setData(Qt::CheckStateRole, Qt::Unchecked);
    item->setData(Qt::UserRole  , 255);
    item->setData(Qt::UserRole+1, 255);
    item->setData(Qt::UserRole+2, 255);
    item->setIcon(makeColorSwatch(QColor(255,255,255)));
    this->Ui->ActiveCoordinationSites->addItem(item);
    }


  vtkLookupTable *lut = vtkLookupTable::New();
  lut->SetNumberOfTableValues(n);

  this->UpdateCoordinationSiteColors(lut);

  vtkPolyDataMapper *mapper
     = dynamic_cast<vtkPolyDataMapper*>(this->CoordinationSiteActor->GetMapper());
  mapper->SetScalarModeToUsePointData();
  mapper->SelectColorArray("label ids");
  mapper->SetColorModeToMapScalars();
  mapper->ScalarVisibilityOn();
  mapper->SetScalarRange(0,n-1);
  mapper->SetLookupTable(lut);
  lut->Delete();
}

// --------------------------------------------------------------------------
void PovmolMainWindow::UpdateCoordinationSites()
{
  if (this->Reader)
    {
    int n = this->Ui->ActiveCoordinationSites->count();
    for (int i=0; i<n; ++i)
      {
      QListWidgetItem *item = this->Ui->ActiveCoordinationSites->item(i);
      if (item->checkState() == Qt::Checked)
        {
        this->Reader->ActivateCoordinationSite(i);
        }
      else
        {
        this->Reader->DeactivateCoordinationSite(i);
        }
      }
    this->Render();
    }
}

// --------------------------------------------------------------------------
void PovmolMainWindow::UpdateCoordinationSiteColors(vtkLookupTable *lut)
{
  int nColors = this->Ui->ActiveCoordinationSites->count();
  for (int i=0; i<nColors; ++i)
    {
    QListWidgetItem *site = this->Ui->ActiveCoordinationSites->item(i);

    double r = site->data(Qt::UserRole  ).toDouble()/255.0;
    double g = site->data(Qt::UserRole+1).toDouble()/255.0;
    double b = site->data(Qt::UserRole+2).toDouble()/255.0;

    lut->SetTableValue(i, r, g, b, 1.0);
    }
  lut->Modified();
}

// --------------------------------------------------------------------------
void PovmolMainWindow::UpdateCoordinationSiteColor(QListWidgetItem *site)
{
  QColor currentColor(
    site->data(Qt::UserRole  ).toInt(),
    site->data(Qt::UserRole+1).toInt(),
    site->data(Qt::UserRole+2).toInt());

  QColor newColor = QColorDialog::getColor(currentColor, this);
  if (newColor.isValid())
    {
    site->setIcon(makeColorSwatch(newColor));
    site->setData(Qt::UserRole  , newColor.red());
    site->setData(Qt::UserRole+1, newColor.green());
    site->setData(Qt::UserRole+2, newColor.blue());
    }

  vtkPolyDataMapper *mapper
    = dynamic_cast<vtkPolyDataMapper*>(this->CoordinationSiteActor->GetMapper());

  vtkLookupTable *lut
    = dynamic_cast<vtkLookupTable*>(mapper->GetLookupTable());

  this->UpdateCoordinationSiteColors(lut);

  this->Render();
}

// --------------------------------------------------------------------------
void PovmolMainWindow::EnableCoordinationSites(bool v)
{
  if (this->Reader)
    {
    if (v)
      {
      this->Reader->SetGenerateCoordinationSites(1);
      }
    else
      {
      this->Reader->SetGenerateCoordinationSites(0);
      }
    this->Render();
    }
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

  if (this->CoordinationSiteActor)
    {
    this->CoordinationSiteActor->Delete();
    this->CoordinationSiteActor = NULL;
    }

  this->Ui->ActiveTransforms->clear();
}

// --------------------------------------------------------------------------
void PovmolMainWindow::BuildPipeline()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::BuildPipeline" << endl;
  #endif

  // molecule
  this->MoleculeMapper = vtkMoleculeMapper2::New();
  this->MoleculeMapper->SetInputConnection(this->Reader->GetOutputPort(0));
  this->MoleculeMapper->UseBallAndStickSettings();
  this->MoleculeMapper->SetAtomicRadiusScaleFactor(this->Ui->AtomRadiusFactor->text().toDouble());
  this->MoleculeMapper->SetBondRadius(this->Ui->BondRadius->text().toDouble());
  this->MoleculeMapper->SetPOVRayStreaming(true);

  unsigned char color[3];
  this->MoleculeMapper->GetBondColor(color);
  this->Ui->BondColor->setIcon(makeColorSwatch(QColor(color[0], color[1], color[2])));

  vtkActor *molActor = vtkActor::New();
  molActor->SetMapper(this->MoleculeMapper);

  this->Renderer->AddActor(molActor);
  molActor->Delete();

  this->Ui->SiteProperties->SetProperty(molActor->GetProperty());

  // coordination polyhedra
  vtkPolyDataMapper *polyhedraMapper = vtkPolyDataMapper::New();
  polyhedraMapper->SetInputConnection(this->Reader->GetOutputPort(1));

  this->CoordinationSiteActor = vtkActor::New();
  this->Ui->CoordinationSiteProperties->SetProperty(this->CoordinationSiteActor->GetProperty());

  //this->CoordinationSiteActor->GetProperty()->EdgeVisibilityOn();
  this->CoordinationSiteActor->SetMapper(polyhedraMapper);
  polyhedraMapper->Delete();

  this->Renderer->AddActor(this->CoordinationSiteActor);

  // axes
  vtkPolyDataMapper *axesMapper = vtkPolyDataMapper::New();
  axesMapper->SetInputConnection(this->Reader->GetOutputPort(2));
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
    this, tr("Bond Geometry"), "", tr("VTK Legacy (*.vtk)")).toStdString();

  if (fileName != "")
    {
    vtkMoleculeToBondStickFilter *bf = vtkMoleculeToBondStickFilter::New();
    bf->SetInputData(this->MoleculeMapper->GetInput());

    vtkPolyDataWriter *w = vtkPolyDataWriter::New();
    w->SetInputConnection(bf->GetOutputPort());
    w->SetFileName(fileName.c_str());
    w->SetFileTypeToBinary();
    w->Write();
    w->Delete();
    bf->Delete();
    }

  fileName = QFileDialog::getSaveFileName(
    this, tr("Atom Geometry"), "", tr("VTK Legacy (*.vtk)")).toStdString();

  if (fileName != "")
    {
    vtkPolyData *atoms = this->MoleculeMapper->GetAtoms();
    vtkPolyDataWriter *w = vtkPolyDataWriter::New();
    w->SetInputData(atoms);
    w->SetFileTypeToBinary();
    w->SetFileName(fileName.c_str());
    w->Write();
    w->Delete();
    atoms->Delete();
    }

  fileName = QFileDialog::getSaveFileName(
    this, tr("Coordination Polyhedra Geometry"), "", tr("VTK Legacy (*.vtk)")).toStdString();

  if (fileName != "")
    {
    vtkPolyDataWriter *w = vtkPolyDataWriter::New();
    w->SetInputConnection(this->Reader->GetOutputPort(1));
    w->SetFileTypeToBinary();
    w->SetFileName(fileName.c_str());
    w->Write();
    w->Delete();
    }
}

// --------------------------------------------------------------------------
void PovmolMainWindow::WritePOV()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::WritePOV" << endl;
  #endif

  // molecule
  string fileName = QFileDialog::getSaveFileName(
    this, tr("Export molecule"), "", tr("POVRay (*.pov)")).toStdString();

  if (fileName != "")
    {
    ofstream povf(fileName.c_str());
    povf << this->MoleculeMapper->GetPOVRayStream() << endl;
    cerr << "wrote " << fileName << endl;
    }

  // the rest of the geometry
  fileName = QFileDialog::getSaveFileName(
    this, tr("Export geometry"), "", tr("POVRay (*.pov)")).toStdString();

  if (fileName != "")
    {
    vtkPOVExporter *e = vtkPOVExporter::New();
    e->SetFileName(fileName.c_str());
    e->SetRenderWindow(this->Ui->ViewWidget->GetRenderWindow());
    e->Write();
    e->Delete();
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
    if (id == TABLE)
      {
      this->Reader->SetBondDetector(static_pointer_cast<BondDetector>(this->TableDetector));
      }
    else
      {
      this->PropertiesDetector->SetDetectionMode(id);
      this->Reader->SetBondDetector(static_pointer_cast<BondDetector>(this->PropertiesDetector));
      }
    this->Reader->BondDetectorModified();
    this->Render();
    }
}

// --------------------------------------------------------------------------
void PovmolMainWindow::UpdateBondProximityFactor()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::UpdateBondProximityFactor" << endl;
  #endif
  if (this->Reader)
    {
    this->PropertiesDetector->SetTolerance(this->Ui->BondProximityFactor->text().toDouble());
    this->Reader->BondDetectorModified();
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
void PovmolMainWindow::UpdateBondColor()
{
  unsigned char color[3];
  this->MoleculeMapper->GetBondColor(color);
  QColor currentColor(color[0], color[1], color[2]);
  QColor newColor = QColorDialog::getColor(currentColor, this);
  if (newColor.isValid())
    {
    this->Ui->BondColor->setIcon(makeColorSwatch(newColor));
    this->MoleculeMapper->SetBondColor(newColor.red(), newColor.green(), newColor.blue());
    this->Render();
    }
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
  while ((l = lc->GetNextItem()))
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

// --------------------------------------------------------------------------
void PovmolMainWindow::EditBondTable()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::EditBondTable" << endl;
  #endif

  PovmolBondTableDialog *dialog = new PovmolBondTableDialog(this);

  TableBasedBondDetector *detector
    = dynamic_cast<TableBasedBondDetector*>(this->TableDetector.get());

  int nRows = detector->GetNumberOfRows();
  for (int i=0; i<nRows; ++i)
    {
    dialog->AddRow(
        detector->GetSourceType(i),
        detector->GetDestinationType(i),
        detector->GetMinLength(i),
        detector->GetMaxLength(i));
    }

  if (dialog->exec() == QDialog::Accepted)
    {
    detector->ClearTable();
    int nRows = dialog->GetNumberOfRows();
    for (int i=0; i<nRows; ++i)
      {
      detector->InsertTableEntry(
            dialog->GetSourceType(i),
            dialog->GetDestinationType(i),
            dialog->GetMinLength(i),
            dialog->GetMaxLength(i));
      }
    int id = this->Ui->BondDetectionMode->currentIndex();
    if (id == TABLE)
      {
      this->UpdateBondDetectionMode(id);
      }
    }
}

// --------------------------------------------------------------------------
void PovmolMainWindow::EditLookupTable()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::EditLookupTable" << endl;
  #endif
  if (this->MoleculeMapper)
    {
    vtkLookupTable *atomLut = this->MoleculeMapper->GetAtomLookupTable();
    PovmolLookupTableDialog *dialog = new PovmolLookupTableDialog(this);
    dialog->Initialize(atomLut);
    if (dialog->exec() == QDialog::Accepted)
      {
      dialog->CopyTo(atomLut);
      atomLut->Modified();

      //vtkLookupTable *bondLut = this->MoleculeMapper->GetBondLookupTable();
      //dialog->CopyTo(bondLut);
      //bondLut->Modified();

      this->MoleculeMapper->LookupTableModified();
      this->Render();
      }
    }
}

// --------------------------------------------------------------------------
void PovmolMainWindow::UpdateDuplicates()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::UpdateDuplicates" << endl;
  #endif

  if (this->Reader)
    {
    this->Reader->SetDuplicateAMinus(this->Ui->DuplicateAMinus->text().toInt());
    this->Reader->SetDuplicateBMinus(this->Ui->DuplicateBMinus->text().toInt());
    this->Reader->SetDuplicateCMinus(this->Ui->DuplicateCMinus->text().toInt());
    this->Reader->SetDuplicateAPlus(this->Ui->DuplicateAPlus->text().toInt());
    this->Reader->SetDuplicateBPlus(this->Ui->DuplicateBPlus->text().toInt());
    this->Reader->SetDuplicateCPlus(this->Ui->DuplicateCPlus->text().toInt());
    this->Render();
    }
}

// --------------------------------------------------------------------------
void PovmolMainWindow::UpdateGhostBonds(bool v)
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::UpdateGhostBonds" << endl;
  #endif
  if (this->Reader)
    {
    if (v)
      {
      this->Reader->SetGenerateGhostBonds(1);
      }
    else
      {
      this->Reader->SetGenerateGhostBonds(0);
      }
    this->Render();
    }
}

/*
// --------------------------------------------------------------------------
void PovmolMainWindow::UpdateCoordinationSiteAlpha()
{
  #ifdef PovmolMainWindowDEBUG
  cerr << ":::::PovmolMainWindow::UpdateCoordinationSiteAlpha" << endl;
  #endif
  double alpha
    = static_cast<double>(this->Ui->CoordinationSiteAlpha->value()) /
         static_cast<double>(this->Ui->CoordinationSiteAlpha->maximum());

  if (alpha >= 0.9999999)
    {
    this->EnableDepthPeeling(false);
    }
  else
    {
    this->EnableDepthPeeling(true);
    }

  this->CoordinationSiteActor->GetProperty()->SetOpacity(alpha);

  this->Render();
}
*/
