#ifndef geometry_h
#define geometry_h

#include "Export.h"
#include <vector>
#include <string>

class vtkPolyData;
class vtkMolecule;
class BondDetector;

EXPORT void BuildMolecule(
    vtkMolecule *molecule,
    std::vector<double> &positions,
    std::vector<unsigned short> &numbers,
    std::vector<std::string> &labels,
    std::vector<int> &labelIds,
    std::vector<bool> &sites,
    std::vector<bool> &coordinationSites,
    std::vector<bool> &ghosts,
    const BondDetector *detector);

EXPORT void CopyActiveSites(
    vtkMolecule *source,
    bool ghostBonds,
    vtkMolecule *dest);

EXPORT void BuildAxes(
    vtkPolyData *data,
    const std::vector<double> &axes);

EXPORT void BuildPolyhedra(
    vtkMolecule *molecule,
    vtkPolyData *polyhedra);

#endif
