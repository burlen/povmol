/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCIFMoleculeReader.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkCIFMoleculeReader - Read a CIF file and output a
// vtkMolecule object
// .SECTION Description
#ifndef vtkCIFMoleculeReader_h
#define vtkCIFMoleculeReader_h

#include "vtkDomainsChemistryModule.h" // For export macro
#include "vtkMoleculeAlgorithm.h"

#include <string>
#include <vector>

class vtkMolecule;
class vtkPolyData;

class VTKDOMAINSCHEMISTRY_EXPORT vtkCIFMoleculeReader : public vtkMoleculeAlgorithm
{
public:
  static vtkCIFMoleculeReader *New();
  vtkTypeMacro(vtkCIFMoleculeReader,vtkMoleculeAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Return true if the reader can read the file.
  static
  bool CanReadFile(const char *fileName);

  // Description:
  // Get/Set the name of the CIF file
  void SetFileName(const char *fileName);
  const char *GetFileName(){ return this->FileName.c_str(); }

  // Description:
  // Set/Get mode by which bonds are detected
  enum {
    ATOMIC,
    IONIC,
    COVALENT,
    VANDERWAALS,
    CRYSTAL
    };
  vtkSetMacro(BondDetectionMode, int);
  vtkGetMacro(BondDetectionMode, int);

  vtkSetMacro(BondProximityFactor, double);
  vtkGetMacro(BondProximityFactor, double);

  vtkGetMacro(Initialized, int);

  size_t GetNumberOfTransforms(){ return this->Transforms.size(); }
  const char *GetTransformLabel(size_t i){ return this->TransformLabels[i].c_str(); }
  void ActivateTransform(size_t i);
  void DeactivateTransform(size_t i);
  void DeacivateTransforms();

protected:
  vtkCIFMoleculeReader();
  ~vtkCIFMoleculeReader();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int FillOutputPortInformation(int, vtkInformation*);

private:
  std::string FileName;
  bool Initialized;

  // primative cell
  std::vector<double> CellLengths;
  std::vector<double> CellAngles;
  std::vector<double> CellAxes;

  // crystal structure etc
  std::vector<double> BasisPositions;
  std::vector<unsigned short> BasisTypes;
  std::vector<std::vector<double> > Transforms;
  std::vector<std::string> TransformLabels;
  std::vector<bool> ActiveTransforms;
  std::vector<double> Positions;
  std::vector<unsigned short> Types;

  int BondDetectionMode;
  double BondProximityFactor;

private:
  vtkCIFMoleculeReader(const vtkCIFMoleculeReader&);  // Not implemented.
  void operator=(const vtkCIFMoleculeReader&);  // Not implemented.
};

#endif
