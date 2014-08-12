#ifndef povmolMainWindow_h
#define povmolMainWindow_h

#include <QMainWindow>
#include "ui_PovmolMainWindowUi.h"
#include <memory>

class vtkRenderer;
class vtkCIFMoleculeReader;
class vtkMoleculeMapper2;
class BondDetector;
class TableBasedBondDetector;
class AtomicPropertiesBondDetector;

class PovmolMainWindow : public QMainWindow
{
private:
  Q_OBJECT
public:
  PovmolMainWindow();
  ~PovmolMainWindow();

protected:
  // detection modes
  enum {ATOMIC, IONIC, COVALENT, VANDERWAALS, CRYSTAL, TABLE};

private slots:
  void UpdateBondDetectionMode(int);
  void UpdateAtomRadiusFactor();
  void UpdateBondRadius();
  void UpdateBondProximityFactor();
  void UpdateBondColorMode();
  void UpdateActiveTransforms();
  void UpdateDuplicates();
  void OpenFile();
  void WritePOV();
  void WriteVTK();
  void WriteImage();
  void Initialize();
  void Render();
  void BuildPipeline();
  void UpdateLightIntensity();
  void ShowTransforms();
  void ViewDownX();
  void ViewDownY();
  void ViewDownZ();
  void ViewUpX();
  void ViewUpY();
  void ViewUpZ();
  void EditBondTable();

private:
  Ui_PovmolMainWindowUi *Ui;
  vtkRenderer *Renderer;
  vtkCIFMoleculeReader *Reader;
  vtkMoleculeMapper2 *MoleculeMapper;
  std::shared_ptr<TableBasedBondDetector> TableDetector;
  std::shared_ptr<AtomicPropertiesBondDetector> PropertiesDetector;
};

#endif
