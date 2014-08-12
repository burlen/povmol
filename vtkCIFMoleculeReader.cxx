/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCIFMoleculeReader.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/
#include "vtkCIFMoleculeReader.h"

#include "vtkDataObject.h"
#include "vtkExecutive.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkMolecule.h"
#include "vtkPolyData.h"

#include "AtomicProperties.h"
#include "Math3D.h"
#include "Utility.h"
#include "Geometry.h"
#include "PointTree.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <limits>
#include <iterator>
#include <algorithm>

using std::vector;
using std::cerr;
using std::endl;
using std::istream;
using std::string;
using std::setw;
using std::ifstream;
using std::numeric_limits;
using std::back_inserter;
using std::shared_ptr;

#define vtkCIFMoleculeReaderDEBUG

namespace {

// **************************************************************************
void Print(
    ostream &os,
    const vector<double> &positions,
    const vector<unsigned short> &types)
{
  size_t n = types.size();
  for (size_t i=0; i<n; ++i)
    {
    size_t ii = 3*i;
    os
      << setw(4) << i << setw(6) << types[i] << setw(6) << AtomicProperties::Symbol(types[i])
      << setw(15) << positions[ii] << setw(15) << positions[ii+1] << setw(15) << positions[ii+2]
      << endl;
    }
}

// **************************************************************************
int ParseTransforms(
      istream &file,
      vector<vector<double> > &transforms,
      vector<string> &transformLabels)
{
  size_t n = 0;
  if (!(file >> n))
    {
    pError(cerr) << "length not found." << endl;
    return -1;
    }
  transforms.resize(n);
  transformLabels.resize(n);
  for (size_t i = 0; i<n; ++i)
    {
    string &label = transformLabels[i];
    if (!std::getline(file >> std::ws, label))
      {
      pError(cerr) << i << " label not found" << endl;
      return -1;
      }
    cerr << "  " << label << endl;
    vector<double> &T = transforms[i];
    T.resize(12);
    for (size_t q = 0; q<12; ++q)
      {
      if (!(file >> T[q]))
        {
        pError(cerr) << i << ", " << q << " invalid enry" << endl;
        return -1;
        }
      }
    }
  return 0;
}

// **************************************************************************
int ParseTypes(istream &file, vector<unsigned short> &numbers)
{
  size_t n = 0;
  if (!(file >> n))
    {
    pError(cerr) << "length not found." << endl;
    return -1;
    }
  numbers.resize(n);
  for (size_t i = 0; i<n; ++i)
    {
    string type;
    if (!(file >> type))
      {
      pError(cerr) << i << " invalid type" << endl;
      return -1;
      }
    numbers[i] = AtomicProperties::Number(type);
    }
  return 0;
}

// **************************************************************************
int ParsePositions(istream &file, vector<double> &positions)
{
  size_t n = 0;
  if (!(file >> n))
    {
    pError(cerr) << "length not found." << endl;
    return -1;
    }
  positions.resize(3*n);
  for (size_t i = 0; i<n; ++i)
    {
    if ( !(file >> positions[3*i])
      || !(file >> positions[3*i+1])
      || !(file >> positions[3*i+2]) )
      {
      pError(cerr) << i << " invalid coordinates" << endl;
      return -1;
      }
    }
  return 0;
}

// **************************************************************************
bool Exists(double *pt, const vector<double> &pts, double tol)
{
  //return false;
  size_t n = pts.size()/3;
  for (size_t i=0; i<n; ++i)
    {
    size_t ii=3*i;
    if ( (fabs(pts[ii  ]-pt[0]) < tol)
      && (fabs(pts[ii+1]-pt[1]) < tol)
      && (fabs(pts[ii+2]-pt[2]) < tol) )
      {
      return true;
      }
    }
  return false;
}

// **************************************************************************
void NormalizeCoordinate(vector<double> &points)
{
  double coordMax[3] = {
    numeric_limits<double>::min(),
    numeric_limits<double>::min(),
    numeric_limits<double>::min()
    };
  double coordMin[3] = {
    numeric_limits<double>::max(),
    numeric_limits<double>::max(),
    numeric_limits<double>::max()
    };
  size_t n = points.size()/3;
  for (size_t i=0; i<n; ++i)
    {
    double *pt = &points[3*i];
    if (pt[0] < coordMin[0]) coordMin[0] = pt[0];
    if (pt[1] < coordMin[1]) coordMin[1] = pt[1];
    if (pt[2] < coordMin[2]) coordMin[2] = pt[2];
    if (pt[0] > coordMax[0]) coordMax[0] = pt[0];
    if (pt[1] > coordMax[1]) coordMax[1] = pt[1];
    if (pt[2] > coordMax[2]) coordMax[2] = pt[2];
    }
  double coordMaxMinDiff[3] = {
    coordMax[0] - coordMin[0],
    coordMax[1] - coordMin[1],
    coordMax[2] - coordMin[2]
    };
  for (size_t i=0; i<n; ++i)
    {
    double *pt = &points[3*i];
    pt[0] = (pt[0] - coordMin[0]) / coordMaxMinDiff[0];
    pt[1] = (pt[1] - coordMin[1]) / coordMaxMinDiff[1];
    pt[2] = (pt[2] - coordMin[2]) / coordMaxMinDiff[2];
    }
}

// **************************************************************************
void NormalizeCoordinate(vector<double> &points, int dir)
{
  double coordMax = numeric_limits<double>::min();
  double coordMin = numeric_limits<double>::max();
  size_t n = points.size()/3;
  for (size_t i=0; i<n; ++i)
    {
    double *pt = &points[3*i];
    if (pt[dir] < coordMin) coordMin = pt[dir];
    if (pt[dir] > coordMax) coordMax = pt[dir];
    }
  double coordMaxMinDiff = coordMax - coordMin;
  for (size_t i=0; i<n; ++i)
    {
    double *pt = &points[3*i];
    pt[dir] = (pt[dir] - coordMin) / coordMaxMinDiff;
    }
}

// **************************************************************************
void DuplicatePeriodicPositions(
      vector<double> &points,
      vector<unsigned short> &numbers,
      int dir,
      double tol)
{
  size_t nPoints = numbers.size();
  for (size_t i=0; i<nPoints; ++i)
    {
    size_t ii = 3*i;
    if (fabs(points[ii+dir]) < tol)
      {
      double tmp[3] = {points[ii], points[ii+1], points[ii+2]};
      tmp[dir] = 1.0;
      // position
      points.push_back(tmp[0]);
      points.push_back(tmp[1]);
      points.push_back(tmp[2]);
      // type
      numbers.push_back(numbers[i]);
      }
    }
}

// **************************************************************************
void Mirror(
      vector<double> &points,
      vector<unsigned short> &numbers,
      double *n,
      double *p)
{
  vector<double> mirroredPoints;
  mirroredPoints.reserve(points.size());
  size_t nPoints = points.size()/3;
  for (size_t i=0; i<nPoints; ++i)
    {
    double pt[3];
    Math3D::copy(pt, &points[3*i]);
    Math3D::sub(pt, p);
    Math3D::reflect(pt, n);
    Math3D::add(pt, p);
    if (!Exists(pt, points, 1e-5))
      {
      // position
      mirroredPoints.push_back(pt[0]);
      mirroredPoints.push_back(pt[1]);
      mirroredPoints.push_back(pt[2]);
      // type
      numbers.push_back(numbers[i]);
      }
    }
  std::copy(
      mirroredPoints.begin(),
      mirroredPoints.end(),
      back_inserter(points));
}

// **************************************************************************
void ApplyPeriodicBC(double *pt)
{
  if (pt[0] < 0.0) pt[0] += 1.0;
  else if (pt[0] > 1.0) pt[0] -= 1.0;

  if (pt[1] < 0.0) pt[1] += 1.0;
  else if (pt[1] > 1.0) pt[1] -= 1.0;

  if (pt[2] < 0.0) pt[2] += 1.0;
  else if (pt[2] > 1.0) pt[2] -= 1.0;
}

// **************************************************************************
void ApplyPeriodicBC(vector<double> &points, vector<unsigned short> &numbers)
{
  vector<double> newPoints;
  vector<unsigned short> newNumbers;

  size_t nPoints = numbers.size();
  for (size_t i=0; i<nPoints; ++i)
    {
    size_t ii = 3*i;
    double *pt = &points[ii];

    ApplyPeriodicBC(pt);

    if (!Exists(pt, newPoints, 1e-5))
      {
      newPoints.push_back(pt[0]);
      newPoints.push_back(pt[1]);
      newPoints.push_back(pt[2]);

      newNumbers.push_back(numbers[i]);
      }
    }

  points.swap(newPoints);
  numbers.swap(newNumbers);
}

// **************************************************************************
void ApplyTransform(
      vector<double> &points,
      vector<unsigned short> &numbers,
      const vector<double> &basisPoints,
      const vector<unsigned short> &basisNumbers,
      const vector<double> &transform)
{
//#define USE_KDTREE
#ifdef USE_KDTREE
  PointTree kdTree;
  kdTree.Initialize(points);
#endif

  size_t nPoints = numbers.size();
  size_t nBasisPoints = basisNumbers.size();
  size_t nOut = nPoints + nBasisPoints;
  points.reserve(3*nOut);
  numbers.reserve(nOut);

  for (size_t i=0; i<nBasisPoints; ++i)
    {
    size_t ii = 3*i;
    double pt[3];
    Math3D::transform(&transform[0], Math3D::copy(pt, &basisPoints[ii]));
    //ApplyPeriodicBC(pt);
#ifdef USE_KDTREE
    if (!kdTree.Exists(pt))
#else
    if (!Exists(pt, points, 1e-5))
#endif
      {
      // save position
      points.push_back(pt[0]);
      points.push_back(pt[1]);
      points.push_back(pt[2]);
      // save type
      numbers.push_back(basisNumbers[i]);
      }
    }
}

// **************************************************************************
void ApplyTransforms(
      vector<double> &points,
      vector<unsigned short> &numbers,
      const vector<double> &basisPoints,
      const vector<unsigned short> &basisNumbers,
      const vector<vector<double> > &transforms,
      const vector<string> &transformLabels,
      const vector<bool> &activeTransforms)
{
  size_t nTransforms = transforms.size();
  for (size_t i = 0; i<nTransforms; ++i)
    {
    if (activeTransforms[i])
      {
      cerr << "transform " << transformLabels[i] << endl;
      ApplyTransform(points, numbers, basisPoints, basisNumbers, transforms[i]);
      }
    }
}


// **************************************************************************
void ComputePrimitiveCell(
       vector<double> &primvec,
       const vector<double> &lengths,
       const vector<double> &angles)
{
  const double &A = lengths[0];
  const double &B = lengths[1];
  const double &C = lengths[2];

  const double &alpha = Math3D::deg2rad(angles[0]);
  const double &beta = Math3D::deg2rad(angles[1]);
  const double &gamma = Math3D::deg2rad(angles[2]);

  primvec.resize(9, 0.0);
  double *a = &primvec[0];
  double *b = &primvec[3];
  double *c = &primvec[6];

  // a
  a[0] = A;

  // b
  b[0] = B;
  Math3D::rotatek(b, gamma);

  // c
  c[0] = C;
  Math3D::rotatek(c, gamma);
  Math3D::rotatei(c, alpha);
  Math3D::rotate(c, b, M_PI/2.0 - beta);

  cerr << "a = [" << a[0] << ", " << a[1] << ", " << a[2] << "]" << endl;
  cerr << "b = [" << b[0] << ", " << b[1] << ", " << b[2] << "]" << endl;
  cerr << "c = [" << c[0] << ", " << c[1] << ", " << c[2] << "]" << endl;
}

// **************************************************************************
void ComputePrimitiveCellPositions(
      vector<double> &points,
      const vector<double> &primvec)
{
  const double *a = &primvec[0];
  const double *b = &primvec[3];
  const double *c = &primvec[6];

  size_t nPts = points.size()/3;
  for (size_t i = 0; i<nPts; ++i)
    {
    size_t ii = 3*i;
    double *pt = &points[ii];

    double p[3] = {a[0], a[1], a[2]};
    double q[3] = {b[0], b[1], b[2]};
    double r[3] = {c[0], c[1], c[2]};

    Math3D::scale(p, pt[0]);
    Math3D::scale(q, pt[1]);
    Math3D::scale(r, pt[2]);
    Math3D::add(pt, p, q);
    Math3D::add(pt, pt, r);
    }
}

// **************************************************************************
void CopyTranslate(
      const vector<double> &basisPoints,
      const vector<unsigned short> &basisNumbers,
      double *offset,
      vector<double> &points,
      vector<unsigned short> &numbers)
{
  size_t nPts = basisPoints.size()/3;
  for (size_t i = 0; i<nPts; ++i)
    {
    double pt[3];
    Math3D::copy(pt, &basisPoints[3*i]);
    Math3D::add(pt, offset);

    if (!Exists(pt, points, 1e-5))
      {
      points.push_back(pt[0]);
      points.push_back(pt[1]);
      points.push_back(pt[2]);

      numbers.push_back(basisNumbers[i]);
      }
    }
}

// **************************************************************************
void ComputeDuplicates(
      vector<double> &points,
      vector<unsigned short> &numbers,
      const vector<double> &primvec,
      unsigned int nAPlus,
      unsigned int nAMinus,
      unsigned int nBPlus,
      unsigned int nBMinus,
      unsigned int nCPlus,
      unsigned int nCMinus)
{
  vector<double> basisPoints(points);
  vector<unsigned short> basisNumbers(numbers);

  const double *a = &primvec[0];
  const double *b = &primvec[3];
  const double *c = &primvec[6];

  // a hat
  for (unsigned int i=1; i<=nAPlus; ++i)
    {
    double offset[3];
    Math3D::copy(offset, a);
    Math3D::scale(offset, static_cast<double>(i));
    CopyTranslate(basisPoints, basisNumbers, offset, points, numbers);
    }

  // -a hat
  for (unsigned int i=1; i<=nAMinus; ++i)
    {
    double offset[3];
    Math3D::copy(offset, a);
    Math3D::scale(offset, -static_cast<double>(i));
    CopyTranslate(basisPoints, basisNumbers, offset, points, numbers);
    }

  if (nAMinus || nAPlus)
    {
    basisPoints.assign(points.begin(), points.end());
    basisNumbers.assign(numbers.begin(), numbers.end());
    }

  // b hat
  for (unsigned int i=1; i<=nBPlus; ++i)
    {
    double offset[3];
    Math3D::copy(offset, b);
    Math3D::scale(offset, static_cast<double>(i));
    CopyTranslate(basisPoints, basisNumbers, offset, points, numbers);
    }

  // -b hat
  for (unsigned int i=1; i<=nBMinus; ++i)
    {
    double offset[3];
    Math3D::copy(offset, b);
    Math3D::scale(offset, -static_cast<double>(i));
    CopyTranslate(basisPoints, basisNumbers, offset, points, numbers);
    }

  if (nBMinus || nBPlus)
    {
    basisPoints.assign(points.begin(), points.end());
    basisNumbers.assign(numbers.begin(), numbers.end());
    }

  // c hat
  for (unsigned int i=1; i<=nCPlus; ++i)
    {
    double offset[3];
    Math3D::copy(offset, c);
    Math3D::scale(offset, static_cast<double>(i));
    CopyTranslate(basisPoints, basisNumbers, offset, points, numbers);
    }

  // -c hat
  for (unsigned int i=1; i<=nCMinus; ++i)
    {
    double offset[3];
    Math3D::copy(offset, c);
    Math3D::scale(offset, -static_cast<double>(i));
    CopyTranslate(basisPoints, basisNumbers, offset, points, numbers);
    }
}

// **************************************************************************
int Parse(
    istream &file,
    vector<double> &lengths,
    vector<double> &angles,
    vector<double> &positions,
    vector<unsigned short> &numbers,
    vector<vector<double> > &transforms,
    vector<string> &transformLabels)
{
  string type, ver;
  file >> type >> ver;
  if (type != "CIFPP")
    {
    pError(cerr) << "read failed, this is not a cifpp file" << endl;
    return -1;
    }

  lengths.resize(3, 0.0);
  angles.resize(3, 0.0);

  string label;
  while (file && file >> label)
    {
    cerr << "parsing " << label << endl;
    if ((label == "_cell_length_a") && !(file >> lengths[0]))
      {
      pError(cerr) << "Failed to read " << label << endl;
      return -1;
      }
    else
    if ((label == "_cell_length_b") && !(file >> lengths[1]))
      {
      pError(cerr) << "Failed to read " << label << endl;
      return -1;
      }
    else
    if ((label == "_cell_length_c") && !(file >> lengths[2]))
      {
      pError(cerr) << "Failed to read " << label << endl;
      return -1;
      }
    else
    if ((label == "_cell_angle_alpha") && !(file >> angles[0]))
      {
      pError(cerr) << "Failed to read " << label << endl;
      return -1;
      }
    else
    if ((label == "_cell_angle_beta") && !(file >> angles[1]))
      {
      pError(cerr) << "Failed to read " << label << endl;
      return -1;
      }
    else
    if ((label == "_cell_angle_gamma") && !(file >> angles[2]))
      {
      pError(cerr) << "Failed to read " << label << endl;
      return -1;
      }
    else
    if ( (label == "_symmetry_transforms")
      && ParseTransforms(file, transforms, transformLabels) )
      {
      pError(cerr) << "Failed to read transforms" << endl;
      return -1;
      }
    else
    if ( (label == "_atom_site_type_symbol")
      && ParseTypes(file, numbers) )
      {
      pError(cerr) << "Failed to read types" << endl;
      return -1;
      }
    else
    if ( (label == "_atom_site_fract_xyz")
      && ParsePositions(file, positions) )
      {
      pError(cerr) << "Failed to read positions" << endl;
      return -1;
      }
    }
  return 0;
}
};

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkCIFMoleculeReader);

//----------------------------------------------------------------------------
vtkCIFMoleculeReader::vtkCIFMoleculeReader() :
      Initialized(false),
      DuplicateAPlus(0),
      DuplicateBPlus(0),
      DuplicateCPlus(0),
      DuplicateAMinus(0),
      DuplicateBMinus(0),
      DuplicateCMinus(0)
{
  #ifdef vtkCIFMoleculeReaderDEBUG
  cerr << "=====vtkCIFMoleculeReader::vtkCIFMoleculeReader" << endl;
  #endif
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(2);
}

//----------------------------------------------------------------------------
vtkCIFMoleculeReader::~vtkCIFMoleculeReader()
{
  #ifdef vtkCIFMoleculeReaderDEBUG
  cerr << "=====vtkCIFMoleculeReader::~vtkCIFMoleculeReader" << endl;
  #endif
  this->Detector = nullptr;
}


//----------------------------------------------------------------------------
void vtkCIFMoleculeReader::SetBondDetector(const shared_ptr<BondDetector> &detector)
{
  if (detector.get() == this->Detector.get())
    {
    return;
    }
  this->Detector = detector;
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkCIFMoleculeReader::BondDetectorModified()
{
  this->Modified();
}

//----------------------------------------------------------------------------
int vtkCIFMoleculeReader::FillOutputPortInformation(
      int port,
      vtkInformation *info)
{
  #ifdef vtkCIFMoleculeReaderDEBUG
  cerr << "=====vtkCIFMoleculeReader::FillOutputPortInformation" << endl;
  #endif
  switch (port)
    {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMolecule");
      break;
    case 1:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
      break;
    default:
      vtkErrorMacro("Invalid output port " << port);
      break;
    }
  return 1;
}

//----------------------------------------------------------------------------
bool vtkCIFMoleculeReader::CanReadFile(const char *fileName)
{
  #ifdef vtkCIFMoleculeReaderDEBUG
  cerr << "=====vtkCIFMoleculeReader::CanReadFile" << endl;
  #endif
  string ext = extension(fileName);
  if (ext == ".cifpp") return true;
  return false;
}

//----------------------------------------------------------------------------
void vtkCIFMoleculeReader::DeacivateTransforms()
{
  #ifdef vtkCIFMoleculeReaderDEBUG
  cerr << "=====vtkCIFMoleculeReader::DeacivateTransforms" << endl;
  #endif
  size_t n = this->Transforms.size();
  for (size_t i=0; i<n; ++i)
    {
    this->DeactivateTransform(i);
    }
}

//----------------------------------------------------------------------------
void vtkCIFMoleculeReader::ActivateTransform(size_t i)
{
  #ifdef vtkCIFMoleculeReaderDEBUG
  cerr << "=====vtkCIFMoleculeReader::ActivateTransform " << i << endl;
  #endif
  if (!this->ActiveTransforms[i])
    {
    this->ActiveTransforms[i] = true;
    this->Modified();
    }
}
//----------------------------------------------------------------------------
void vtkCIFMoleculeReader::DeactivateTransform(size_t i)
{
  #ifdef vtkCIFMoleculeReaderDEBUG
  cerr << "=====vtkCIFMoleculeReader::DeactivateTransform " << i << endl;
  #endif
  if (this->ActiveTransforms[i])
    {
    this->ActiveTransforms[i] = false;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkCIFMoleculeReader::SetFileName(const char *cFileName)
{
  #ifdef vtkCIFMoleculeReaderDEBUG
  cerr << "=====vtkCIFMoleculeReader::SetFileName" << endl;
  #endif
  string fileName(cFileName);
  if (cFileName == this->FileName)
    {
    return;
    }

  this->FileName = cFileName;
  this->Initialized = false;
  this->CellLengths.clear();
  this->CellLengths.resize(3, 0.0);
  this->CellAngles.clear();
  this->CellAngles.resize(3, 0.0);
  this->BasisPositions.clear();
  this->BasisTypes.clear();
  this->Transforms.clear();
  this->TransformLabels.clear();
  this->DuplicateAPlus = 0;
  this->DuplicateBPlus = 0;
  this->DuplicateCPlus = 0;
  this->DuplicateAMinus = 0;
  this->DuplicateBMinus = 0;
  this->DuplicateCMinus = 0;

  ifstream file(this->FileName.c_str());
  if (!file)
    {
    vtkErrorMacro("failed to open " << this->FileName);
    return;
    }

  if (Parse(
        file,
        this->CellLengths,
        this->CellAngles,
        this->BasisPositions,
        this->BasisTypes,
        this->Transforms,
        this->TransformLabels))
    {
    vtkErrorMacro("failed to parse " <<this->FileName);
    }

  ::Print(cerr, this->BasisPositions, this->BasisTypes);

  this->ActiveTransforms.clear();
  this->ActiveTransforms.resize(this->Transforms.size(), true);

  this->CellAxes.clear();
  this->CellAxes.resize(9, 0.0);
  ComputePrimitiveCell(this->CellAxes, this->CellLengths, this->CellAngles);

  this->Initialized = true;
  this->Modified();
}

//----------------------------------------------------------------------------
int vtkCIFMoleculeReader::RequestData(
  vtkInformation *request,
  vtkInformationVector **inInfos,
  vtkInformationVector *outInfo)
{
  #ifdef vtkCIFMoleculeReaderDEBUG
  cerr << "=====vtkCIFMoleculeReader::RequestData" << endl;
  #endif
  int port = request->Get(vtkExecutive::FROM_OUTPUT_PORT());
//  if (port == 0)
    {
    vtkMolecule *molecule
      = vtkMolecule::SafeDownCast(vtkDataObject::GetData(outInfo, 0));

    if (!molecule)
      {
      vtkErrorMacro("Empty molecule output");
      return 1;
      }

    this->Positions.clear();
    this->Types.clear();
    ApplyTransforms(
          this->Positions,
          this->Types,
          this->BasisPositions,
          this->BasisTypes,
          this->Transforms,
          this->TransformLabels,
          this->ActiveTransforms);

    //ApplyPeriodicBC(this->Positions, this->Types);

    // TODO -- works for this particualr case but not sure
    // how is this encoded in the file
    //double n[3] = {0.0, 1.0, 0.0};
    //double p[3] = {0.0, 0.5, 0.0};
    //Mirror(this->Positions, this->Types, n, p);
    // Duplicating across the periodic boundary solves it!

    DuplicatePeriodicPositions(this->Positions, this->Types, 0, 1.0e-3);
    DuplicatePeriodicPositions(this->Positions, this->Types, 1, 1.0e-3);
    DuplicatePeriodicPositions(this->Positions, this->Types, 2, 1.0e-3);

    cerr << "transformed positions" << endl;
    ::Print(cerr, this->Positions, this->Types);
    cerr << endl;

    ComputePrimitiveCellPositions(
          this->Positions,
          this->CellAxes);

    cerr << "primative cell positions" << endl;
    ::Print(cerr, this->Positions, this->Types);
    cerr << endl;

    ComputeDuplicates(
        this->Positions,
        this->Types,
        this->CellAxes,
        this->DuplicateAPlus,
        this->DuplicateAMinus,
        this->DuplicateBPlus,
        this->DuplicateBMinus,
        this->DuplicateCPlus,
        this->DuplicateCMinus);

    BuildMolecule(
        molecule,
        this->Positions,
        this->Types,
        this->Detector.get());

    #if 0 && defined(vtkCIFMoleculeReaderDEBUG)
    cerr << "port 0" << endl;
    molecule->Print(cerr);
    cerr << endl;
    #endif
    }
//  else
//  if (port == 1)
    {
    vtkPolyData *axes
      = vtkPolyData::SafeDownCast(vtkDataObject::GetData(outInfo, 1));

    if (!axes)
      {
      vtkErrorMacro("Empty axes output");
      return 1;
      }

    BuildAxes(axes, this->CellAxes);

    #if 0 && defined(vtkCIFMoleculeReaderDEBUG)
    cerr << "port 1" << endl;
    axes->Print(cerr);
    cerr << endl;
    #endif
    }
//  else
//    {
//    vtkErrorMacro("invalid port " << port);
//    }

  return 1;
}

//----------------------------------------------------------------------------
void vtkCIFMoleculeReader::PrintSelf(ostream& os, vtkIndent indent)
{
  #ifdef vtkCIFMoleculeReaderDEBUG
  cerr << "=====vtkCIFMoleculeReader::PrintSelf" << endl;
  #endif
  this->Superclass::PrintSelf(os,indent);
}
