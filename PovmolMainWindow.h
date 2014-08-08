#ifndef povmolMainWindow_h
#define povmolMainWindow_h

#include <QMainWindow>
#include "ui_PovmolMainWindowUi.h"

class vtkRenderer;
class vtkCIFMoleculeReader;
class vtkMoleculeMapper2;

class PovmolMainWindow : public QMainWindow
{
private:
  Q_OBJECT
public:
  PovmolMainWindow();
  ~PovmolMainWindow();

protected:
  //virtual void closeEvent(QCloseEvent *event);

private slots:
  void UpdateBondDetectionMode(int);
  void UpdateAtomRadiusFactor();
  void UpdateBondRadius();
  void UpdateBondProximityFactor();
  void UpdateBondColorMode();
  void OpenFile();
  void WritePOV();
  void WriteVTK();
  void WriteImage();
  void Initialize();
  void Render();
  void BuildPipeline();
  void UpdateLightIntensity();
  void ShowTransforms();
  void UpdateActiveTransforms();
  void ViewDownX();
  void ViewDownY();
  void ViewDownZ();
  void ViewUpX();
  void ViewUpY();
  void ViewUpZ();


private:
  Ui_PovmolMainWindowUi *Ui;
  vtkRenderer *Renderer;
  vtkCIFMoleculeReader *Reader;
  vtkMoleculeMapper2 *MoleculeMapper;
};

#endif
