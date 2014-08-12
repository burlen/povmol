#ifndef geometry_h
#define geometry_h

#include "Export.h"
#include <vector>

class vtkPolyData;
class vtkMolecule;
class BondDetector;

EXPORT void BuildMolecule(
    vtkMolecule *molecule,
    std::vector<double> &positions,
    std::vector<unsigned short> &numbers,
    const BondDetector *detector);

EXPORT void BuildAxes(
    vtkPolyData *data,
    const std::vector<double> &axes);

#endif
