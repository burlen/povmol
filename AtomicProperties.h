#ifndef atomicProperties_h
#define atomicProperties_h

#include "Export.h"
#include <string>

namespace AtomicProperties
{

enum {
    ATOMIC,
    IONIC,
    COVALENT,
    VANDERWAALS,
    CRYSTAL
    };

/**
Compute a radius of interaction that could be
used to identify molecular bonds. Mode aregument
determiines which type of radius is computed.
Supported modes are: ATOMIC, IONIC, COVALENT,
VANDERWAALS, CRYSTAL
*/
EXPORT double Radius(unsigned short number, int mode=ATOMIC);
EXPORT double Radius(std::string &sym, int mode=ATOMIC);

/**
Given a symbol from the periodic table convert to
its atomic number.
*/
EXPORT unsigned short Number(const std::string &sym);

/**
Given an atomic number get its symbol in the
periodic table
*/
EXPORT std::string Symbol(unsigned short number);

};

#endif
