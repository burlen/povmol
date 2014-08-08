#ifndef geometry_h
#define geometry_h

#include "Export.h"
#include <vector>

class vtkPolyData;
class vtkMolecule;

EXPORT void BuildMolecule(
    vtkMolecule *molecule,
    std::vector<double> &positions,
    std::vector<unsigned short> &numbers,
    int detectionMode,
    double proximityFactor);

EXPORT void BuildAxes(
    vtkPolyData *data,
    const std::vector<double> &axes);

#endif
