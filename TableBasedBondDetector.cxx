#include "TableBasedBondDetector.h"

#include "Math3D.h"
#include "Utility.h"
#include "AtomicProperties.h"

#include <iostream>

using std::cerr;
using std::endl;
using std::map;
using std::pair;

namespace {

// ***************************************************************************
unsigned int log2(unsigned int v)
{
  unsigned int r = 0; // r will be lg(v)
  while (v >>= 1) { ++r; }
  return r;
}

// ***************************************************************************
unsigned long makeKey(unsigned short source, unsigned short destination)
{
  unsigned int shortBits = 1 << log2(sizeof(unsigned short)*8);
  unsigned long key = source;
  key <<= shortBits;
  key |= destination;
  return key;
}

// **************************************************************************
unsigned short sourceFromKey(unsigned long key)
{
  unsigned int shortBits = 1 << log2(sizeof(unsigned short)*8);
  key >>= shortBits;
  return key;
}

// **************************************************************************
unsigned short destinationFromKey(unsigned long key)
{
  unsigned int shortBits = 1 << log2(sizeof(unsigned short)*8);
  unsigned long mask = (1 << shortBits) - 1;
  key &= mask;
  return key;
}

};

// ---------------------------------------------------------------------------
bool TableBasedBondDetector::Connected(
    unsigned short a1,
    double *pos1,
    unsigned short a2,
    double *pos2) const
{
  double tmp[3];
  double d = Math3D::length(Math3D::sub(tmp, pos2, pos1));

  // skip expensive table lookup if these atoms are
  // further away than the largest bond distance in the
  // table.
  if (d > this->MaxBondLength)
    {
    return false;
    }

  unsigned long key = makeKey(a1, a2);

  // check that this combination is in the table
  // if not assume that there's no desire to generate
  // bonds between these two types of atoms
  if (!this->Table.count(key))
    {
    return false;
    }

  const pair<double, double> &range = this->Table.at(key);

  return ( (d >= range.first) && (d <= range.second));
}

// ---------------------------------------------------------------------------
void TableBasedBondDetector::InsertTableEntry(
    unsigned short source,
    unsigned short destination,
    double minLength,
    double maxLength)
{
  if (this->MaxBondLength < maxLength)
    {
    this->MaxBondLength = maxLength;
    }

  pair<double,double> range(minLength, maxLength);
  this->Table[makeKey(source, destination)] = range;
  this->Table[makeKey(destination, source)] = range;
}

// ---------------------------------------------------------------------------
int TableBasedBondDetector::GetNumberOfRows() const
{
  // 2 because of duplicate to handle src->dest and dest->src
  return static_cast<int>(this->Table.size()/2);
}

// ---------------------------------------------------------------------------
unsigned short TableBasedBondDetector::GetSourceType(int row) const
{
  row *= 2; // 2 because of duplicate to handle src->dest and dest->src
  auto it = this->Table.begin();
  for (int i=0; i<row; ++i) ++it;
  return sourceFromKey(it->first);
}

// ---------------------------------------------------------------------------
unsigned short TableBasedBondDetector::GetDestinationType(int row) const
{
  row *= 2; // 2 because of duplicate to handle src->dest and dest->src
  auto it = this->Table.begin();
  for (int i=0; i<row; ++i) ++it;
  return destinationFromKey(it->first);
}

// ---------------------------------------------------------------------------
double TableBasedBondDetector::GetMinLength(int row) const
{
  row *= 2; // 2 because of duplicate to handle src->dest and dest->src
  auto it = this->Table.begin();
  for (int i=0; i<row; ++i) ++it;
  return it->second.first;
}

// ---------------------------------------------------------------------------
double TableBasedBondDetector::GetMaxLength(int row) const
{
  row *= 2; // 2 because of duplicate to handle src->dest and dest->src
  auto it = this->Table.begin();
  for (int i=0; i<row; ++i) ++it;
  return it->second.second;
}
// ---------------------------------------------------------------------------
void TableBasedBondDetector::ClearTable()
{
  this->MaxBondLength = 0.0;
  this->Table.clear();
}
