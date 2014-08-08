#include "Math3D.h"

#include <cmath>
#include <limits>
#include <algorithm>

namespace Math3D
{

// ---------------------------------------------------------------------------
void initbounds(double *b)
{
  b[0] =  std::numeric_limits<double>::max();
  b[1] = -std::numeric_limits<double>::max();
  b[2] =  std::numeric_limits<double>::max();
  b[3] = -std::numeric_limits<double>::max();
  b[4] =  std::numeric_limits<double>::max();
  b[5] = -std::numeric_limits<double>::max();
}

// ---------------------------------------------------------------------------
void bounds(double *b, double *a)
{
  b[0] = std::min(b[0],a[0]);
  b[1] = std::max(b[1],a[0]);
  b[2] = std::min(b[2],a[1]);
  b[3] = std::max(b[3],a[1]);
  b[4] = std::min(b[4],a[2]);
  b[5] = std::max(b[5],a[2]);
}

// ---------------------------------------------------------------------------
double frac(double a, double b, double f)
{
  return (b-a)*f+a;
}

// ---------------------------------------------------------------------------
double *abs(double *b, double *a)
{
  b[0] = fabs(a[0]);
  b[1] = fabs(a[1]);
  b[2] = fabs(a[2]);
  return b;
}

// ---------------------------------------------------------------------------
double *scale(double *a, double b)
{
  a[0] *= b;
  a[1] *= b;
  a[2] *= b;
  return a;
}

// ---------------------------------------------------------------------------
double *scale(double *c, double *a, double b)
{
  c[0] = a[0] * b;
  c[1] = a[1] * b;
  c[2] = a[2] * b;
  return c;
}

// ---------------------------------------------------------------------------
double *sub(double *c, double *b, double *a)
{
  c[0] = b[0] - a[0];
  c[1] = b[1] - a[1];
  c[2] = b[2] - a[2];
  return c;
}

// ---------------------------------------------------------------------------
double *add(double *c, double *a, double *b)
{
  c[0] = b[0] + a[0];
  c[1] = b[1] + a[1];
  c[2] = b[2] + a[2];
  return c;
}

// ---------------------------------------------------------------------------
double dot(double *a, double *b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// ---------------------------------------------------------------------------
double length(double *a, double *b)
{
    double r[3];
    sub(r,b,a);
    return sqrt(dot(r,r));
}

// ---------------------------------------------------------------------------
double length(double *a)
{
    return sqrt(dot(a,a));
}

// ---------------------------------------------------------------------------
double *copy(double *a, const double *b)
{
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
  return a;
}

// ---------------------------------------------------------------------------
double *normalize(double *b, double *a)
{
  double r = length(a);
  b[0] = a[0]/r;
  b[1] = a[1]/r;
  b[2] = a[2]/r;
  return b;
}

// ---------------------------------------------------------------------------
void defplane(double *abcd, double *n, double *r0)
{
  double nn[3];
  copy(abcd,normalize(nn,n));
  abcd[3] = -dot(nn,r0);
}

// ---------------------------------------------------------------------------
double evalplane(double *abcd, double *xyz)
{
  return dot(abcd,xyz) + abcd[3];
}

// ---------------------------------------------------------------------------
double *reflect(double *p, double *n)
{
  double a = n[0];
  double b = n[1];
  double c = n[2];
  double M[12] = {
      1.0-2.0*a*a,    -2.0*a*b,    -2.0*a*c, 0.0,
         -2.0*a*b, 1.0-2.0*b*b,    -2.0*b*c, 0.0,
         -2.0*a*c,    -2.0*b*c, 1.0-2.0*c*c, 0.0
      };
  return transform(M, p);
}

// ---------------------------------------------------------------------------
double *transform(const double *M, double *a)
{
  double b[4] = {a[0], a[1], a[2]};
  a[0] = M[0]*b[0] + M[1]*b[1] + M[2]*b[2] + M[3];
  a[1] = M[4]*b[0] + M[5]*b[1] + M[6]*b[2] + M[7];
  a[2] = M[8]*b[0] + M[9]*b[1] + M[10]*b[2] + M[11];
  return a;
}


// ---------------------------------------------------------------------------
double deg2rad(double deg)
{
  return M_PI/180.0*deg;
}

// ---------------------------------------------------------------------------
double *rotatei(double *vec, double rad)
{
  double s = sin(rad);
  double c = cos(rad);
  double y = vec[1];
  double z = vec[2];
  vec[1] = c*y - s*z;
  vec[2] = s*y + c*z;
  return vec;
}

// ---------------------------------------------------------------------------
double *rotatej(double *vec, double rad)
{
  double s = sin(rad);
  double c = cos(rad);
  double x = vec[0];
  double z = vec[2];
  vec[0] = c*x + s*z;
  vec[2] = -s*x + c*z;
  return vec;
}

// ---------------------------------------------------------------------------
double *rotatek(double *vec, double rad)
{
  double s = sin(rad);
  double c = cos(rad);
  double x = vec[0];
  double y = vec[1];
  vec[0] = c*x - s*y;
  vec[1] = s*x + c*y;
  return vec;
}

// ---------------------------------------------------------------------------
double *rotate(double *vec, const double *axis, double rad)
{
  double r = sqrt(axis[0]*axis[0]+axis[1]*axis[1]);
  double t = atan(axis[1]/axis[0]);
  double p = atan(r/axis[2]);
  rotatek(vec, -t);
  rotatej(vec, -p);
  rotatek(vec, rad);
  rotatej(vec, p);
  rotatek(vec, t);
  return vec;
}

};
