#include "AtomicPropertiesBondDetector.h"
#include "AtomicProperties.h"
#include "Math3D.h"

// --------------------------------------------------------------------------
bool AtomicPropertiesBondDetector::Connected(
    unsigned short a1,
    double *pos1,
    unsigned short a2,
    double *pos2)
    const
{
  double tmp[3];
  return ( Math3D::length(Math3D::sub(tmp, pos2, pos1))
    < this->Tolerance*( AtomicProperties::Radius(a1, this->DetectionMode)
            + AtomicProperties::Radius(a2, this->DetectionMode) ) );
}
