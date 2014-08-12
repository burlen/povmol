#ifndef BondLengthDelegate_h
#define BondLengthDelegate_h

#include <QStyledItemDelegate>

class QWidget;
class QModelIndex;
class QSytleOptionViewItem;
class QAbstractItemModel;

class BondLengthDelegate : public QStyledItemDelegate
{
private:
  Q_OBJECT
public:
  explicit BondLengthDelegate(QObject *parent = 0);
  ~BondLengthDelegate();

  QWidget *createEditor(
        QWidget *parent,
        const QStyleOptionViewItem &option,
        const QModelIndex &index) const;

  void setEditorData(
        QWidget *editor,
        const QModelIndex &index) const;

  void setModelData(
        QWidget *editor,
        QAbstractItemModel *model,
        const QModelIndex & index) const;

  void updateEditorGeometry(
        QWidget *editor,
        const QStyleOptionViewItem &option,
        const QModelIndex &index) const;

};

#endif
