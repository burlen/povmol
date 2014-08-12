#include "PovmolBondTableDialog.h"

#include <QStandardItemModel>
#include <QStandardItem>
#include <QStringListModel>
#include <QStringList>
#include <QList>
#include <QHeaderView>

#include "AtomicProperties.h"
#include "AtomTypeDelegate.h"
#include "BondLengthDelegate.h"

#include <string>
#include <iostream>

using std::cerr;
using std::endl;
using std::string;

// --------------------------------------------------------------------------
PovmolBondTableDialog::PovmolBondTableDialog(QWidget *parent)
      : QDialog(parent)
{
  this->Ui = new Ui_PovmolBondTableDialogUi;
  this->Ui->setupUi(this);

  QStringList headerLabels;
  headerLabels << tr("Source") << tr("Destination") << tr("Min") << tr("Max");

  this->BondTableModel = new QStandardItemModel(0, 4, this);
  this->BondTableModel->setHorizontalHeaderLabels(headerLabels);

  this->Ui->BondTableView->setModel(this->BondTableModel);
  this->Ui->BondTableView->setItemDelegateForColumn(0, new AtomTypeDelegate(this));
  this->Ui->BondTableView->setItemDelegateForColumn(1, new AtomTypeDelegate(this));
  this->Ui->BondTableView->setItemDelegateForColumn(2, new BondLengthDelegate(this));
  this->Ui->BondTableView->setItemDelegateForColumn(3, new BondLengthDelegate(this));
  this->Ui->BondTableView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  this->Ui->BondTableView->verticalHeader()->setVisible(false);
  this->ResizeTable();

  connect(this->Ui->AddRow, SIGNAL(released()), this, SLOT(AddRow()));
  connect(this->Ui->RemoveRow, SIGNAL(released()), this, SLOT(RemoveRow()));
  connect(this->BondTableModel, SIGNAL(itemChanged(QStandardItem*)), this, SLOT(ResizeTable()));

}

// --------------------------------------------------------------------------
PovmolBondTableDialog::~PovmolBondTableDialog()
{
  cerr << ":::::PovmolBondTableDialog::~PovmolBondTableDialog" << endl;
}

// --------------------------------------------------------------------------
void PovmolBondTableDialog::ResizeTable()
{
  this->Ui->BondTableView->resizeColumnsToContents();
  int width = 5;
  for (int i=0; i<4; ++i)
    {
    width += this->Ui->BondTableView->columnWidth(i);
    }
  this->Ui->BondTableView->setFixedWidth(width);
}

// --------------------------------------------------------------------------
void PovmolBondTableDialog::AddRow()
{
  QList<QStandardItem*> row;
  row << new QStandardItem("");
  row << new QStandardItem("");
  row << new QStandardItem("0.0");
  row << new QStandardItem("0.0");
  this->BondTableModel->appendRow(row);
  this->ResizeTable();
}

// --------------------------------------------------------------------------
void PovmolBondTableDialog::AddRow(
      unsigned short sourceType,
      unsigned short destType,
      double minLength,
      double maxLength)
{
  QList<QStandardItem*> row;
  row << new QStandardItem(AtomicProperties::Symbol(sourceType).c_str());
  row << new QStandardItem(AtomicProperties::Symbol(destType).c_str());
  row << new QStandardItem(QString("%1").arg(minLength));
  row << new QStandardItem(QString("%2").arg(maxLength));
  this->BondTableModel->appendRow(row);
  this->ResizeTable();
}

// --------------------------------------------------------------------------
void PovmolBondTableDialog::RemoveRow()
{
  QModelIndex idx = this->Ui->BondTableView->currentIndex();
  this->BondTableModel->removeRow(idx.row());
}

// --------------------------------------------------------------------------
int PovmolBondTableDialog::GetNumberOfRows() const
{
  return this->BondTableModel->rowCount();
}

// --------------------------------------------------------------------------
unsigned short PovmolBondTableDialog::GetSourceType(int row) const
{
  string symbol = this->BondTableModel->item(row, 0)->text().toStdString();
  return AtomicProperties::Number(symbol);
}

// --------------------------------------------------------------------------
unsigned short PovmolBondTableDialog::GetDestinationType(int row) const
{
  string symbol = this->BondTableModel->item(row, 1)->text().toStdString();
  return AtomicProperties::Number(symbol);
}

// --------------------------------------------------------------------------
double PovmolBondTableDialog::GetMinLength(int row) const
{
  return this->BondTableModel->item(row, 2)->text().toDouble();
}

// --------------------------------------------------------------------------
double PovmolBondTableDialog::GetMaxLength(int row) const
{
  return this->BondTableModel->item(row, 3)->text().toDouble();
}
