#ifndef TableBasedBondDetector_h
#define TableBasedBondDetector_h

#include "BondDetector.h"
#include <map>
#include <utility>

/// Interface to bond detector classes
class TableBasedBondDetector : public BondDetector
{
public:
  virtual ~TableBasedBondDetector() {}

  /**
  Return true if the two atoms are close enough to
  be connected
  */
  virtual bool Connected(
      unsigned short a1,
      double *pos1,
      unsigned short a2,
      double *pos2) const;

  /**
  Insert entries for source->destination and destination->source
  bonds. Overwrites existing entry or adds a new one.
  */
  void InsertTableEntry(
      unsigned short source,
      unsigned short destination,
      double minLenghth,
      double maxLength);

  /**
  Remove the entries in the table.
  */
  void ClearTable(){ this->Table.clear(); }

  /**
  Access cells in the table by row,
  */
  int GetNumberOfRows() const;
  unsigned short GetSourceType(int row) const;
  unsigned short GetDestinationType(int row) const;
  double GetMinLength(int row) const;
  double GetMaxLength(int row) const;

private:
  typedef std::map<unsigned int,std::pair<double,double > > TableType;
  TableType Table;
};

#endif
