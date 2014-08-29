#include "vtkVector.h"
#include "vtkVectorOperators.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include <iostream>
using std::ostream;
#include <string>
using std::string;

// Description:
// Point primative for POV Ray.
class POVRayPoint : public vtkVector3f
{
public:
  POVRayPoint(float v)
    {
    this->Data[0] = this->Data[1] = this->Data[2] = v;
    }

  POVRayPoint(vtkVector3f &v)
    {
    this->Data[0] = v[0];
    this->Data[1] = v[1];
    this->Data[2] = v[2];
    }

  POVRayPoint(float v0, float v1, float v2)
    {
    this->Data[0] = v0;
    this->Data[1] = v1;
    this->Data[2] = v2;
    }
};

// Description:
// Color primative for POV Ray.
class POVRayColor : public vtkVector<float,5>
{
public:
  POVRayColor() : Name("DefaultColor")
    {
    this->Data[0] = this->Data[1] = this->Data[2]
      = this->Data[3] = this->Data[4] = 0.0f;
    }

  POVRayColor(string name, float v, float f=0.0f, float t=0.0f) : Name(name)
    {
    this->Data[0] = this->Data[1] = this->Data[2] = v;
    this->Data[3] = f;
    this->Data[4] = t;
    }

  POVRayColor(string name, const vtkVector3f &v) : Name(name)
    {
    this->Data[0] = v[0];
    this->Data[1] = v[1];
    this->Data[2] = v[2];
    this->Data[3] = 0.0f;
    this->Data[4] = 0.0f;
    }

  POVRayColor(
        string name,
        unsigned char v0,
        unsigned char v1,
        unsigned char v2,
        float f=0.0f,
        float t=0.0f)
        :
    Name(name)
    {
    this->Data[0] = v0/255.0f;
    this->Data[1] = v1/255.0f;
    this->Data[2] = v2/255.0f;
    this->Data[3] = f;
    this->Data[4] = t;
    }

  POVRayColor(
        string name,
        double v0,
        double v1,
        double v2,
        float f=0.0f,
        float t=0.0f)
        :
    Name(name)
    {
    this->Data[0] = v0;
    this->Data[1] = v1;
    this->Data[2] = v2;
    this->Data[3] = f;
    this->Data[4] = t;
    }

  POVRayColor(
      string name,
      float v0,
      float v1,
      float v2,
      float f=0.0f,
      float t=0.0f)
        :
    Name(name)
    {
    this->Data[0] = v0;
    this->Data[1] = v1;
    this->Data[2] = v2;
    this->Data[3] = f;
    this->Data[4] = t;
    }
  string Name;
};


// Description:
// Finish primative for POV Ray.
class POVRayFinish
{
public:
  POVRayFinish()
          :
      Name("DefaultFinish"),
      Ambient(1.0f), Diffuse(1.0f), Specular(1.0f),
      Roughness(0.0f),
      Phong(1.0f), PhongSize(1.0f),
      Reflection(0.0f),
      Metal(false)
      {}

  POVRayFinish(
      string name,
      float ambient,
      float diffuse,
      float specular,
      float roughness,
      float phong,
      float phongSize,
      float reflection,
      bool metal)
          :
      Name(name),
      Ambient(ambient),
      Diffuse(diffuse),
      Specular(specular),
      Roughness(roughness),
      Phong(phong),
      PhongSize(phongSize),
      Reflection(reflection),
      Metal(metal)
      {}

  string Name;
  float Ambient;
  float Diffuse;
  float Specular;
  float Roughness;
  float Phong;
  float PhongSize;
  float Reflection;
  bool Metal;
};

// Description:
// Cylinder primative for POV Ray.
class POVRayCylinder
{
public:
  POVRayCylinder(
      vtkVector3f &point1,
      vtkVector3f &point2,
      float radius,
      POVRayColor &color,
      POVRayFinish &finish)
          :
      Point1(point1),
      Point2(point2),
      Radius(radius),
      Color(color),
      Finish(finish)
      {}

  POVRayPoint Point1;
  POVRayPoint Point2;
  float Radius;
  POVRayColor Color;
  POVRayFinish Finish;
};


// Description:
// Sphere primative for POV Ray.
class POVRaySphere
{
public:
  POVRaySphere(
      const string &name,
      vtkVector3f &center,
      float radius,
      POVRayColor &color,
      POVRayFinish &finish)
          :
      Name(name),
      Center(center),
      Radius(radius),
      Color(color),
      Finish(finish)
      {}

  string Name;
  POVRayPoint Center;
  float Radius;
  POVRayColor Color;
  POVRayFinish Finish;
};


// Description:
// Camera primative for POV Ray.
class POVRayCamera
{
public:
  POVRayCamera(vtkCamera *cam) : Position(0.0f), FocalPoint(0.0f), ViewAngle(30.0f)
    {
    double pt[3];
    cam->GetPosition(pt);
    this->Position[0] = pt[0];
    this->Position[1] = pt[1];
    this->Position[2] = pt[2];

    cam->GetFocalPoint(pt);
    this->FocalPoint[0] = pt[0];
    this->FocalPoint[1] = pt[1];
    this->FocalPoint[2] = pt[2];

    this->ViewAngle = cam->GetViewAngle();
    }

  POVRayPoint Position;
  POVRayPoint FocalPoint;
  float ViewAngle;

};

// Description:
// Light primative for POV Ray.
class POVRayLight
{
public:
  POVRayLight(vtkLight *light) : Position(0.0f), Color("lightColor", 1.0f)
    {
    double v[3];
    light->GetPosition(v);
    this->Position[0] = v[0];
    this->Position[1] = v[1];
    this->Position[2] = v[2];
    /*
    light->GetAmbientColor(v);
    this->Color[0] = v[0];
    this->Color[1] = v[1];
    this->Color[2] = v[2];
    */
    }

  POVRayPoint Position;
  POVRayColor Color;
};

// Description:
// ostream overloads for the POV Ray primatives.
ostream &operator<<(ostream &os, POVRayPoint &v);
ostream &operator<<(ostream &os, POVRayColor &c);
ostream &operator<<(ostream &os, POVRayFinish &f);
ostream &operator<<(ostream &os, POVRayCylinder &c);
ostream &operator<<(ostream &os, POVRaySphere &s);
ostream &operator<<(ostream &os, POVRayLight &light);
ostream &operator<<(ostream &os, POVRayCamera &cam);

vtkVector3f &operator+=(vtkVector3f &left, vtkVector3f &right);
vtkVector3f &operator*=(vtkVector3f &left, float right);
