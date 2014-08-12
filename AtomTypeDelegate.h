#ifndef AtomTypeDelegate_h
#define AtomTypeDelegate_h

#include <QStyledItemDelegate>

class QWidget;
class QModelIndex;
class QSytleOptionViewItem;
class QAbstractItemModel;

class AtomTypeDelegate : public QStyledItemDelegate
{
private:
  Q_OBJECT
public:
  explicit AtomTypeDelegate(QObject *parent = 0);
  ~AtomTypeDelegate();

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
