#include "BondLengthDelegate.h"

#include <QWidget>
#include <QLineEdit>
#include <QDoubleValidator>
#include <QStringList>
#include <QAbstractItemModel>
#include <QModelIndex>

#include <AtomicProperties.h>

// --------------------------------------------------------------------------
BondLengthDelegate::BondLengthDelegate(QObject *parent)
      : QStyledItemDelegate(parent)
{
}

// --------------------------------------------------------------------------
BondLengthDelegate::~BondLengthDelegate()
{
}

// --------------------------------------------------------------------------
QWidget *BondLengthDelegate::createEditor(
      QWidget *parent,
      const QStyleOptionViewItem &option,
      const QModelIndex &index)
      const
{
  QDoubleValidator *dv = new QDoubleValidator(parent);
  dv->setBottom(0.0);

  QLineEdit *edit = new QLineEdit(parent);
  edit->setValidator(dv);

  return edit;
}

// --------------------------------------------------------------------------
void BondLengthDelegate::setEditorData(
      QWidget *editor,
      const QModelIndex &index)
      const
{
  QLineEdit *edit = dynamic_cast<QLineEdit*>(editor);
  edit->setText(index.data().toString());
}

// --------------------------------------------------------------------------
void BondLengthDelegate::setModelData(
      QWidget *editor,
      QAbstractItemModel *model,
      const QModelIndex & index)
      const
{
  QLineEdit *edit = dynamic_cast<QLineEdit*>(editor);
  model->setData(index, edit->text());
}

// --------------------------------------------------------------------------
void BondLengthDelegate::updateEditorGeometry(
      QWidget *editor,
      const QStyleOptionViewItem &option,
      const QModelIndex &index)
      const
{
  editor->setGeometry(option.rect);
}
