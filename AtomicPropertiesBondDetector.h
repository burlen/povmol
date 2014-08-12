#ifndef AtomicPropertiesBondDetector_h
#define AtomicPropertiesBondDetector_h

#include "BondDetector.h"

#include "AtomicProperties.h"

class AtomicPropertiesBondDetector : public BondDetector
{
public:
  AtomicPropertiesBondDetector() :
    Tolerance(0.05),
    DetectionMode(AtomicProperties::ATOMIC)
    {}

  virtual ~AtomicPropertiesBondDetector() {}

  /**
  Sets minimun difference to consider connected
  */
  void SetTolerance(double tol)
  { this->Tolerance = tol; }

  /**
  Set the mode used for detections. Can be one of
  ATOMIC, IONIC, COVALENT, VANDERWAALS, CRYSTAL
  defined in AtomicProeprties.h
  */
  void SetDetectionMode(int mode)
  { this->DetectionMode = mode; }

  /**
  Return true if the two atoms are close enought to
  be connected in the corrent mode.
  */
  virtual bool Connected(
      unsigned short a1,
      double *pos1,
      unsigned short a2,
      double *pos2) const;

private:
  double Tolerance;
  int DetectionMode;
};

#endif
