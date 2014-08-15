#include "PovmolLookupTableDialog.h"

#include "AtomicProperties.h"

#include <vtkLookupTable.h>
#include <vtkFloatArray.h>

#include <QListWidget>
#include <QListWidgetItem>
#include <QPixmap>
#include <QIcon>
#include <QColor>
#include <QColorDialog>

#include <iostream>
using std::cerr;
using std::endl;


namespace {
// **************************************************************************
QIcon makeColorSwatch(QColor color)
{
  QPixmap patch(24, 24);
  patch.fill(color);
  QIcon swatch(patch);
  return swatch;
}

// **************************************************************************
void Print(
    ostream &os,
    unsigned char r,
    unsigned char g,
    unsigned char b,
    unsigned char a)
{
  os
    << setw(4) << int(r)
    << setw(4) << int(g)
    << setw(4) << int(b)
    << setw(4) << int(a);
}

// **************************************************************************
void Print(ostream &os, QColor color)
{
  os
    << setw(4) << color.red()
    << setw(4) << color.green()
    << setw(4) << color.blue()
    << setw(4) << color.alpha();
}

// **************************************************************************
void Print(ostream &os, unsigned char *rgba)
{
  os
    << setw(4) << int(rgba[0])
    << setw(4) << int(rgba[1])
    << setw(4) << int(rgba[2])
    << setw(4) << int(rgba[3]);
}
};

// --------------------------------------------------------------------------
PovmolLookupTableDialog::PovmolLookupTableDialog(QWidget *parent) :
    QDialog(parent)
{
  this->Ui = new Ui_PovmolLookupTableDialog;
  this->Ui->setupUi(this);

  connect(this->Ui->ColorList, SIGNAL(itemDoubleClicked(QListWidgetItem*)),
       this, SLOT(EditColorEntry(QListWidgetItem*)));
}

// --------------------------------------------------------------------------
PovmolLookupTableDialog::~PovmolLookupTableDialog()
{
}

// --------------------------------------------------------------------------
void PovmolLookupTableDialog::Initialize(vtkLookupTable *lut)
{
  //size_t nColors = lut->GetNumberOfTableValues();
  size_t nColors = AtomicProperties::MaxNumber();
  this->Colors.reserve(nColors);
  unsigned char *colors = lut->GetTable()->GetPointer(0);
  for (size_t i=1; i<nColors; ++i)
    {
    size_t ii = 4*i;
    unsigned char r = colors[ii];
    unsigned char g = colors[ii+1];
    unsigned char b = colors[ii+2];
    unsigned char a = colors[ii+3];

    QColor color(r,g,b,a);

    QListWidgetItem *colorEntry
       = new QListWidgetItem(
            makeColorSwatch(color),
            tr(AtomicProperties::Symbol(i).c_str()));

    colorEntry->setData(Qt::UserRole  , color.red());
    colorEntry->setData(Qt::UserRole+1, color.green());
    colorEntry->setData(Qt::UserRole+2, color.blue());
    colorEntry->setData(Qt::UserRole+3, color.alpha());

    this->Ui->ColorList->addItem(colorEntry);

    cerr << setw(4) << i << setw(4) << AtomicProperties::Symbol(i);
       Print(cerr, color); cerr << endl;
    }
  this->Ui->ColorList->sortItems();
}

// --------------------------------------------------------------------------
void PovmolLookupTableDialog::CopyTo(vtkLookupTable *lut)
{
  size_t nColorEntries = this->Ui->ColorList->count();
  unsigned char *colors = lut->GetTable()->GetPointer(0);
  for (size_t i=0; i<nColorEntries; ++i)
    {
    QListWidgetItem *colorEntry = this->Ui->ColorList->item(i);

    unsigned short number
      = AtomicProperties::Number(colorEntry->text().toStdString());

    size_t ii = 4*number;

    colors[ii  ] = colorEntry->data(Qt::UserRole  ).toInt();
    colors[ii+1] = colorEntry->data(Qt::UserRole+1).toInt();
    colors[ii+2] = colorEntry->data(Qt::UserRole+2).toInt();
    colors[ii+3] = colorEntry->data(Qt::UserRole+3).toInt();

    cerr << setw(4) << number << setw(4) << AtomicProperties::Symbol(number);
       Print(cerr, colors+ii); cerr << endl;
    }
}

// --------------------------------------------------------------------------
void PovmolLookupTableDialog::EditColorEntry(QListWidgetItem *colorEntry)
{
  QColor currentColor(
    colorEntry->data(Qt::UserRole  ).toInt(),
    colorEntry->data(Qt::UserRole+1).toInt(),
    colorEntry->data(Qt::UserRole+2).toInt(),
    colorEntry->data(Qt::UserRole+3).toInt());

  unsigned short number
    = AtomicProperties::Number(colorEntry->text().toStdString());

  cerr << setw(4) << AtomicProperties::Symbol(number);
      Print(cerr, currentColor); cerr << endl;

  QColor newColor = QColorDialog::getColor(currentColor, this);
  if (newColor.isValid())
    {
    colorEntry->setIcon(makeColorSwatch(newColor));
    colorEntry->setData(Qt::UserRole  , newColor.red());
    colorEntry->setData(Qt::UserRole+1, newColor.green());
    colorEntry->setData(Qt::UserRole+2, newColor.blue());
    colorEntry->setData(Qt::UserRole+3, newColor.alpha());

    cerr << setw(4) << number << setw(4) << AtomicProperties::Symbol(number);
        Print(cerr, newColor); cerr << endl;
    }
}
