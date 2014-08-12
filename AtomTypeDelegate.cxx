#include "AtomTypeDelegate.h"

#include <QWidget>
#include <QComboBox>
#include <QStringList>
#include <QAbstractItemModel>
#include <QModelIndex>

#include <AtomicProperties.h>

// --------------------------------------------------------------------------
AtomTypeDelegate::AtomTypeDelegate(QObject *parent)
      : QStyledItemDelegate(parent)
{
}

// --------------------------------------------------------------------------
AtomTypeDelegate::~AtomTypeDelegate()
{
}

// --------------------------------------------------------------------------
QWidget *AtomTypeDelegate::createEditor(
      QWidget *parent,
      const QStyleOptionViewItem &option,
      const QModelIndex &index)
      const
{
  QStringList types;
  types << "";
  unsigned short nTypes = AtomicProperties::MaxNumber();
  for (unsigned short i=1; i<=nTypes; ++i)
    {
    types << AtomicProperties::Symbol(i).c_str();
    }
  types.sort();

  QComboBox *typeCombo = new QComboBox(parent);
  typeCombo->addItems(types);

  return typeCombo;
}

// --------------------------------------------------------------------------
void AtomTypeDelegate::setEditorData(
      QWidget *editor,
      const QModelIndex &index)
      const
{
  QComboBox *typeCombo = dynamic_cast<QComboBox*>(editor);
  QString symbol = index.data().toString();
  unsigned short number = 0;
  if (symbol != "")
    {
    number = AtomicProperties::Number(symbol.toStdString());
    }
  typeCombo->setCurrentIndex(number);
}

// --------------------------------------------------------------------------
void AtomTypeDelegate::setModelData(
      QWidget *editor,
      QAbstractItemModel *model,
      const QModelIndex & index)
      const
{
  QComboBox *typeCombo = dynamic_cast<QComboBox*>(editor);
  model->setData(index, typeCombo->currentText());
}

// --------------------------------------------------------------------------
void AtomTypeDelegate::updateEditorGeometry(
      QWidget *editor,
      const QStyleOptionViewItem &option,
      const QModelIndex &index)
      const
{
  editor->setGeometry(option.rect);
}
