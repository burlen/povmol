#ifndef BondDetector_h
#define BondDetector_h

/// Interface to bond detector classes
class BondDetector
{
public:
  virtual ~BondDetector() {}
  /**
  Return true if the two atoms are close enough to
  be connected
  */
  virtual bool Connected(
      unsigned short a1,
      double *pos1,
      unsigned short a2,
      double *pos2) const = 0;
};

#endif
