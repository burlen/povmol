#include "QVTKPropertyWidget.h"

#include <vtkProperty.h>

// --------------------------------------------------------------------------
QVTKPropertyWidget::QVTKPropertyWidget(QWidget *parent) :
    QWidget(parent),
    Prop(0)
{
  this->Ui = new Ui_QVTKPropertyWidget;
  this->Ui->setupUi(this);

  connect(this->Ui->Ambient, SIGNAL(valueChanged(int)), this, SLOT(UpdateAmbient()));
  connect(this->Ui->Diffuse, SIGNAL(valueChanged(int)), this, SLOT(UpdateDiffuse()));
  connect(this->Ui->Specular, SIGNAL(valueChanged(int)), this, SLOT(UpdateSpecular()));
  connect(this->Ui->SpecularPower, SIGNAL(valueChanged(int)), this, SLOT(UpdateSpecularPower()));
  connect(this->Ui->Opacity, SIGNAL(valueChanged(int)), this, SLOT(UpdateOpacity()));
}

// --------------------------------------------------------------------------
QVTKPropertyWidget::~QVTKPropertyWidget()
{
  this->SetProperty(0);
}

// --------------------------------------------------------------------------
void QVTKPropertyWidget::SetProperty(vtkProperty *prop)
{
  if (prop == this->Prop)
   {
   return;
   }

  if (this->Prop)
    {
    this->Prop->Delete();
    }

  this->Prop = prop;

  if (this->Prop)
    {
    this->Prop->Register(0);
    this->Ui->Ambient->setValue(this->Prop->GetAmbient()*this->Ui->Ambient->maximum());
    this->Ui->Diffuse->setValue(this->Prop->GetDiffuse()*this->Ui->Diffuse->maximum());
    this->Ui->Specular->setValue(this->Prop->GetSpecular()*this->Ui->Specular->maximum());
    this->Ui->SpecularPower->setValue(this->Prop->GetSpecularPower()*this->Ui->SpecularPower->maximum());
    this->Ui->Opacity->setValue(this->Prop->GetOpacity()*this->Ui->Opacity->maximum());
    }
}

// --------------------------------------------------------------------------
void QVTKPropertyWidget::UpdateAmbient()
{
  if (!this->Prop) return;

  double val
    = static_cast<double>(this->Ui->Ambient->value()) /
         static_cast<double>(this->Ui->Ambient->maximum());

  this->Prop->SetAmbient(val);

  emit Modified();
}

// --------------------------------------------------------------------------
void QVTKPropertyWidget::UpdateDiffuse()
{
  if (!this->Prop) return;

  double val
    = static_cast<double>(this->Ui->Diffuse->value()) /
         static_cast<double>(this->Ui->Diffuse->maximum());

  this->Prop->SetDiffuse(val);

  emit Modified();
}

// --------------------------------------------------------------------------
void QVTKPropertyWidget::UpdateSpecular()
{
  if (!this->Prop) return;

  double val
    = static_cast<double>(this->Ui->Specular->value()) /
         static_cast<double>(this->Ui->Specular->maximum());

  this->Prop->SetSpecular(val);

  emit Modified();
}

// --------------------------------------------------------------------------
void QVTKPropertyWidget::UpdateSpecularPower()
{
  if (!this->Prop) return;

  double val
    = static_cast<double>(this->Ui->SpecularPower->value()) /
         static_cast<double>(this->Ui->SpecularPower->maximum());

  this->Prop->SetSpecularPower(val);

  emit Modified();
}

// --------------------------------------------------------------------------
void QVTKPropertyWidget::UpdateOpacity()
{
  if (!this->Prop) return;

  double val
    = static_cast<double>(this->Ui->Opacity->value()) /
         static_cast<double>(this->Ui->Opacity->maximum());

  this->Prop->SetOpacity(val);

  emit Translucent(val < 0.999999);
  emit Modified();
}

