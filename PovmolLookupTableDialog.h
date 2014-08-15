#ifndef PovmolLookupTableDialog_h
#define PovmolLookupTableDialog_h

#include <QDialog>
#include "ui_PovmolLookupTableDialog.h"
#include <vector>
#include <QColor>

class vtkLookupTable;

class PovmolLookupTableDialog : public QDialog
{
private:
  Q_OBJECT
public:
  explicit PovmolLookupTableDialog(QWidget *parent = 0);
  ~PovmolLookupTableDialog();

  void Initialize(vtkLookupTable *lut);
  void CopyTo(vtkLookupTable *lut);

public slots:
  void EditColorEntry(QListWidgetItem *colorEntry);

private:
  Ui_PovmolLookupTableDialog *Ui;
  std::vector<QColor> Colors;
};

#endif

