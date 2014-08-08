#include "Utility.h"

using std::string;

// --------------------------------------------------------------------------
void fixPath(string &path)
{
#ifndef _WIN32
  (void)path;
#else
  size_t n = path.size();
  for (size_t i=0; i<n; ++i)
    {
    if (path[i] == '/')
      {
      path[i] = '\\';
      }
    }
#endif
}

// --------------------------------------------------------------------------
string extension(string filename)
{
  size_t at = filename.rfind('.');
  string ext = filename.substr(at);
  return ext;
}
