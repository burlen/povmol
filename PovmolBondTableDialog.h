#ifndef PovmolBondTableDialog_h
#define PovmolBondTableDialog_h

#include <QDialog>
#include "ui_PovmolBondTableDialogUi.h"
class QStandardItemModel;

//#include "AtomTypeDelegate.h"

class PovmolBondTableDialog : public QDialog
{
private:
  Q_OBJECT
public:
  explicit PovmolBondTableDialog(QWidget *parent = 0);
  ~PovmolBondTableDialog();

  int GetNumberOfRows() const;
  unsigned short GetSourceType(int row) const;
  unsigned short GetDestinationType(int row) const;
  double GetMinLength(int row) const;
  double GetMaxLength(int row) const;

  void AddRow(
      unsigned short sourceType,
      unsigned short destType,
      double minLength,
      double maxLength);

public slots:
  void AddRow();
  void RemoveRow();
  void ResizeTable();

private:
  Ui_PovmolBondTableDialogUi *Ui;
  QStandardItemModel *BondTableModel;
  //AtomTypeDelegate *AtomTypes;
};

#endif
