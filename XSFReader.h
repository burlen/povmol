#ifndef XSF_h
#define XSF_h

#include "Export.h"

namespace XSF
{
/**
Read atoms and grid data from XSF file. Periodic transformations
are applied.
*/
EXPORT
int Read(
    const string &fileName,
    vtkMolecule *molecule,
    vtkStructuredGrid *sgrid);
};

#endif
