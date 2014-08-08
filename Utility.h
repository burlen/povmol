#ifndef utility_h
#define utility_h

#include "Export.h"

#include <iostream>
#include <string>

#define pError(s) \
  s << "Error: " << __FILE__ << ":" << __LINE__ << std::endl

/**
convert unix path for windows
*/
EXPORT
void fixPath(std::string &path);

/**
given a file name return extension
*/
EXPORT
std::string extension(std::string filename);

#endif
