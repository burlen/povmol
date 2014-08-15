#ifndef QVTKPropertyWidget_h
#define QVTKPropertyWidget_h

#include <QWidget>
#include "ui_QVTKPropertyWidget.h"

class vtkProperty;

class QVTKPropertyWidget : public QWidget
{
private:
  Q_OBJECT;
public:
  QVTKPropertyWidget(QWidget *parent = 0);
  ~QVTKPropertyWidget();

  void SetProperty(vtkProperty *prop);

signals:
  void Translucent(bool v);
  void Modified();

public slots:
  void UpdateAmbient();
  void UpdateDiffuse();
  void UpdateSpecular();
  void UpdateSpecularPower();
  void UpdateOpacity();

private:
  Ui_QVTKPropertyWidget *Ui;
  vtkProperty *Prop;
};

#endif
