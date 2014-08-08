#ifndef pointMath_h
#define pointMath_h

#include "Export.h"

namespace Math3D
{

/**
initialize bounds vector(6x1)
*/
EXPORT void initbounds(double *b);

/**
given initialized bounds vector(6x1) b include point a(3x1)
*/
EXPORT void bounds(double *b, double *a);

/**
(b-a)*f + a
*/
EXPORT double frac(double a, double b, double f);

/**
compute component wise absolute value of vector(3x1) a
store in vector b. return pointer to b
*/
EXPORT double *abs(double *b, double *a);

/**
multiply vector(3x1) a by value b store the result in c
return a pointer to c
*/
EXPORT double *scale(double *c, double *a, double b);

/**
multiply vector(3x1) a inplace by value b
return a pointer to a
*/
EXPORT double *scale(double *a, double b);

/**
subtract two vectors(3x1) a and b store the result in c
return a pointer to c
*/
EXPORT double *sub(double *c, double *b, double *a);
EXPORT inline double *sub(double *b, double *a){ return sub(b, b, a); }

/**
add two vectors(3x1) a and b store the result in c
return a pointer to c
*/
EXPORT double *add(double *c, double *a, double *b);
EXPORT inline double *add(double *a, double *b){ return add(a, a, b); }

/**
compute dot product of two vectors(3x1)
*/
EXPORT double dot(double *a, double *b);

/**
compute the length of vectors(3x1) b-a
*/
EXPORT double length(double *a, double *b);

/**
compute length of vector(3x1) a
*/
EXPORT double length(double *a);

/**
copy 3x1 b into a return a pointer to a
*/
EXPORT double *copy(double *a, const double *b);

/**
nnormalize the vector a and store in b. return a
pointer to b
*/
EXPORT double *normalize(double *b, double *a);

/**
given a point(3x1) r0 an a normal vector(3x1) compute
coeficients of the plane equation and store in abcd
*/
EXPORT void defplane(double *abcd, double *n, double *r0);

/**
evaluate a point(3x1) in the equation of a plane
given by abcd
*/
EXPORT double evalplane(double *abcd, double *xyz);

/**
reflect a point p through a plane passing through the origin
whose normal is n. n must be normalized.
*/
EXPORT double *reflect(double *p, double *n);

/**
transform a point(3x1) in place a using 3x4 matrix
return pointer to a
*/
EXPORT double *transform(const double *M, double *a);

/**
convert from degrees to radians
*/
EXPORT double deg2rad(double deg);

/**
rotate point(3x1) about coordinate axis
*/
EXPORT double *rotatei(double *vec, double rad);
EXPORT double *rotatej(double *vec, double rad);
EXPORT double *rotatek(double *vec, double rad);

/**
rotate a point(3x1) about arbitrary axis
*/
EXPORT double *rotate(double *vec, const double *axis, double rad);

};

#endif
