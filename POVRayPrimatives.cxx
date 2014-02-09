#include "POVRayPrimatives.h"

//----------------------------------------------------------------------------
ostream &operator<<(ostream &os, POVRayPoint &v)
{
  os << "<" << v[0] << ", " << v[1] << ", " << v[2] << ">";
  return os;
}

//----------------------------------------------------------------------------
ostream &operator<<(ostream &os, POVRayColor &c)
{
  os
    << "#declare " << c.Name << "=rgbft <"
    << c[0] << ", " << c[1] << ", " << c[2] << ", " << c[3] << ", " << c[4]
    << ">;";
  return os;
}

//----------------------------------------------------------------------------
ostream &operator<<(ostream &os, POVRayFinish &f)
{
  os
    << "#declare " << f.Name << "=finish {" << endl
    << " ambient " << f.Ambient << endl
    << " diffuse " << f.Diffuse << endl
    << " specular " << f.Specular << endl
    << " roughness " << f.Roughness << endl
    << " phong " << f.Phong << endl
    << " phong_size " << f.PhongSize << endl
    << " reflection ";
  if (f.Metal)
    {
    os << "{ " << f.Reflection << " metallic }";
    }
  else
    {
    os << f.Reflection;
    }
  os << endl << "}";
  return os;
}

//----------------------------------------------------------------------------
ostream &operator<<(ostream &os, POVRayCylinder &c)
{
  os
    << "cylinder { "
    << c.Point1 << ", " << c.Point2 << ", " << c.Radius << endl
    << "pigment { " << c.Color.Name << "}" << endl
    << "finish { " << c.Finish.Name << "}" << endl
    << "}";
  return os;
}

//----------------------------------------------------------------------------
ostream &operator<<(ostream &os, POVRaySphere &s)
{
  os
    << "sphere { "
    << s.Center << ", " << s.Radius << endl
    << "pigment { " << s.Color.Name << "}" << endl
    << "finish { " << s.Finish.Name << "}" << endl
    << "}";
  return os;
}

//----------------------------------------------------------------------------
ostream &operator<<(ostream &os, POVRayLight &light)
{

  os << "light_source { " << light.Position << " color <1.0, 1.0, 1.0> }";
  return os;
}

//----------------------------------------------------------------------------
ostream &operator<<(ostream &os, POVRayCamera &cam)
{
  os
    << "camera { location " << cam.Position
    << " look_at " << cam.FocalPoint
    << " angle " << cam.ViewAngle << " }";
  return os;
}

//----------------------------------------------------------------------------
vtkVector3f &operator+=(vtkVector3f &left, vtkVector3f &right)
{
  left[0] += right[0];
  left[1] += right[1];
  left[2] += right[2];
  return left;
}

//----------------------------------------------------------------------------
vtkVector3f &operator*=(vtkVector3f &left, float right)
{
  left[0] *= right;
  left[1] *= right;
  left[2] *= right;
  return left;
}
