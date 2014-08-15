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
#include <memory>

class vtkMolecule;
class vtkPolyData;
class BondDetector;

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
  // Set/Get the bond detector
  void SetBondDetector(const std::shared_ptr<BondDetector> &detector);
  void BondDetectorModified();

  // Description:
  // return true when a file has been successfully opened
  // and get'ers return valid results.
  vtkGetMacro(Initialized, int);

  // Description:
  // Return the number of symmetry transforms in the file
  size_t GetNumberOfTransforms(){ return this->Transforms.size(); }

  // Description:
  // Get the i'th symmetry transform's label
  const char *GetTransformLabel(size_t i){ return this->TransformLabels[i].c_str(); }

  // Description:
  // Activate/Deactive the i'th symmetry transform
  void ActivateTransform(size_t i);
  void DeactivateTransform(size_t i);

  // Description:
  // Deactivate all symmetry transforms.
  void DeacivateTransforms();


  // Description:
  // Return the number of symmetry transforms in the file
  size_t GetNumberOfSites(){ return this->BasisLabels.size(); }

  // Description:
  // Get the i'th symmetry transform's label
  const char *GetSiteLabel(size_t i){ return this->BasisLabels[i].c_str(); }

  // Description:
  // Activate/Deactive the i'th symmetry transform
  void ActivateSite(size_t i);
  void DeactivateSite(size_t i);

  // Description:
  // Deactivate all symmetry transforms.
  void DeacivateSites();

  // Description:
  // Activate/Deactive the i'th symmetry transform
  void ActivateCoordinationSite(size_t i);
  void DeactivateCoordinationSite(size_t i);

  // Description:
  // Deactivate all symmetry transforms.
  void DeacivateCoordinationSites();

  // Description:
  // Enable/Disable coordination site generation and rendering.
  vtkSetMacro(GenerateCoordinationSites, int);
  vtkGetMacro(GenerateCoordinationSites, int);

  // Description:
  // Enable/Disable ghost bonds
  vtkSetMacro(GenerateGhostBonds, int);
  vtkGetMacro(GenerateGhostBonds, int);

  // Description:
  // Set the number of duplicate in primitive vector
  // directions
  vtkSetMacro(DuplicateAPlus, unsigned int);
  vtkSetMacro(DuplicateBPlus, unsigned int);
  vtkSetMacro(DuplicateCPlus, unsigned int);
  vtkSetMacro(DuplicateAMinus, unsigned int);
  vtkSetMacro(DuplicateBMinus, unsigned int);
  vtkSetMacro(DuplicateCMinus, unsigned int);
  vtkGetMacro(DuplicateAPlus, unsigned int);
  vtkGetMacro(DuplicateBPlus, unsigned int);
  vtkGetMacro(DuplicateCPlus, unsigned int);
  vtkGetMacro(DuplicateAMinus, unsigned int);
  vtkGetMacro(DuplicateBMinus, unsigned int);
  vtkGetMacro(DuplicateCMinus, unsigned int);

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
  std::vector<std::string> BasisLabels;
  std::vector<int> BasisLabelIds;
  std::vector<bool> BasisSites;
  std::vector<bool> BasisCoordinationSites;
  std::vector<std::vector<double> > Transforms;
  std::vector<std::string> TransformLabels;
  std::vector<bool> ActiveTransforms;
  std::vector<double> Positions;
  std::vector<unsigned short> Types;
  std::vector<std::string> Labels;
  std::vector<int> LabelIds;
  std::vector<bool> Sites;
  std::vector<bool> CoordinationSites;

  std::shared_ptr<BondDetector> Detector;

  unsigned int DuplicateAPlus;
  unsigned int DuplicateBPlus;
  unsigned int DuplicateCPlus;
  unsigned int DuplicateAMinus;
  unsigned int DuplicateBMinus;
  unsigned int DuplicateCMinus;

  int GenerateCoordinationSites;
  int GenerateGhostBonds;

private:
  vtkCIFMoleculeReader(const vtkCIFMoleculeReader&);  // Not implemented.
  void operator=(const vtkCIFMoleculeReader&);  // Not implemented.
};

#endif
