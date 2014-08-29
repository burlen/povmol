/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMoleculeMapper2.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkMoleculeMapper2.h"

#include "vtkActor.h"
#include "vtkColor.h"
#include "vtkCylinderSource.h"
#include "vtkDataArray.h"
#include "vtkEventForwarderCommand.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkGlyph3DMapper.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkLookupTable.h"
#include "vtkMath.h"
#include "vtkMolecule.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkPeriodicTable.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSelection.h"
#include "vtkSelectionNode.h"
#include "vtkSphereSource.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTrivialProducer.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkVector.h"
#include "vtkVectorOperators.h"

#include "vtkLightCollection.h"
#include "vtkLight.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"

#include <string>
using std::string;
#include <sstream>
using std::ostringstream;
#include <utility>
using std::pair;
using std::map;

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkMoleculeMapper2);

//----------------------------------------------------------------------------
vtkMoleculeMapper2::vtkMoleculeMapper2()
  : RenderAtoms(true),
    AtomicRadiusType(VDWRadius),
    AtomicRadiusScaleFactor(0.3),
    RenderBonds(true),
    BondColorMode(DiscreteByAtom),
    UseMultiCylindersForBonds(true),
    BondRadius(0.075),
    UniformBondColor("BondColor", 1.0f),
    AtomFinish(),
    BondFinish(),
    POVRayStreaming(false)
{
  // Initialize ivars:
  this->BondColor[0] = this->BondColor[1] = this->BondColor[2] = 50;

  // Setup glyph sources
  vtkNew<vtkSphereSource> sphere;
  sphere->SetRadius(1.0);
  sphere->SetPhiResolution(50);
  sphere->SetThetaResolution(50);
  sphere->Update();
  this->AtomGlyphMapper->SetSourceConnection(sphere->GetOutputPort());

  vtkNew<vtkCylinderSource> cylinder;
  cylinder->SetRadius(1.0);
  cylinder->SetResolution(20);
  cylinder->SetHeight(1.0);
  cylinder->Update();
  // Rotate the glyph so that the cylinder is aligned with the x-axis,
  // rather than the y-axis. This makes glyph orientation much easier.
  vtkNew<vtkTransform> cylXform;
  cylXform->RotateWXYZ(90, 0.0, 0.0, 1.0);
  vtkNew<vtkTransformPolyDataFilter> cylXformFilter;
  cylXformFilter->SetInputConnection(cylinder->GetOutputPort());
  cylXformFilter->SetTransform(cylXform.GetPointer());
  cylXformFilter->Update();
  this->BondGlyphMapper->SetSourceConnection(cylXformFilter->GetOutputPort());

  // Setup glyph mappers
  vtkNew<vtkLookupTable> lut;
  this->PeriodicTable->GetDefaultLUT(lut.GetPointer());
  this->AtomGlyphMapper->SetLookupTable(lut.GetPointer());
  this->AtomGlyphMapper->SetScalarRange
    (0, this->PeriodicTable->GetNumberOfElements());
  this->AtomGlyphMapper->SetColorModeToMapScalars();
  this->AtomGlyphMapper->SetScalarModeToUsePointFieldData();
  this->AtomGlyphMapper->SetScaleModeToScaleByMagnitude();
  this->BondGlyphMapper->SetScaleModeToScaleByVectorComponents();
  // Bond color mode is setup during updates

  // Forward commands to instance mappers
  vtkNew<vtkEventForwarderCommand> cb;
  cb->SetTarget(this);

  this->AtomGlyphMapper->AddObserver(vtkCommand::StartEvent, cb.GetPointer());
  this->AtomGlyphMapper->AddObserver(vtkCommand::EndEvent, cb.GetPointer());
  this->AtomGlyphMapper->AddObserver(vtkCommand::ProgressEvent,
                                     cb.GetPointer());

  this->BondGlyphMapper->AddObserver(vtkCommand::StartEvent, cb.GetPointer());
  this->BondGlyphMapper->AddObserver(vtkCommand::EndEvent, cb.GetPointer());
  this->BondGlyphMapper->AddObserver(vtkCommand::ProgressEvent,
                                     cb.GetPointer());

  // Connect the trivial producers to forward the glyph polydata
  this->AtomGlyphPointOutput->SetOutput(this->AtomGlyphPolyData.GetPointer());
  this->AtomGlyphMapper->SetInputConnection
    (this->AtomGlyphPointOutput->GetOutputPort());

  this->BondGlyphPointOutput->SetOutput(this->BondGlyphPolyData.GetPointer());
  this->BondGlyphMapper->SetInputConnection
    (this->BondGlyphPointOutput->GetOutputPort());

  // Force the glyph data to be generated on the next render:
  this->GlyphDataInitialized = false;
}

//----------------------------------------------------------------------------
vtkMoleculeMapper2::~vtkMoleculeMapper2()
{
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::SetInputData(vtkMolecule *input)
{
  this->SetInputDataInternal(0, input);
}

//----------------------------------------------------------------------------
vtkMolecule *vtkMoleculeMapper2::GetInput()
{
  return vtkMolecule::SafeDownCast(this->GetExecutive()->GetInputData(0, 0));
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::UseBallAndStickSettings()
{
  this->SetRenderAtoms(true);
  this->SetRenderBonds(true);
  this->SetAtomicRadiusType( VDWRadius );
  this->SetAtomicRadiusScaleFactor( 0.3 );
  this->SetBondColorMode( DiscreteByAtom );
  this->SetUseMultiCylindersForBonds( true );
  this->SetBondRadius( 0.075 );
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::UseVDWSpheresSettings()
{
  this->SetRenderAtoms(true);
  this->SetRenderBonds(true);
  this->SetAtomicRadiusType( VDWRadius );
  this->SetAtomicRadiusScaleFactor( 1.0 );
  this->SetBondColorMode( DiscreteByAtom );
  this->SetUseMultiCylindersForBonds( true );
  this->SetBondRadius( 0.075 );
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::UseLiquoriceStickSettings()
{
  this->SetRenderAtoms(true);
  this->SetRenderBonds(true);
  this->SetAtomicRadiusType( UnitRadius );
  this->SetAtomicRadiusScaleFactor( 0.15 );
  this->SetBondColorMode( DiscreteByAtom );
  this->SetUseMultiCylindersForBonds( false );
  this->SetBondRadius( 0.15 );
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::UseFastSettings()
{
  this->SetRenderAtoms(true);
  this->SetRenderBonds(true);
  this->SetAtomicRadiusType( UnitRadius );
  this->SetAtomicRadiusScaleFactor( 0.60 );
  this->SetBondColorMode( SingleColor );
  this->SetBondColor( 50, 50, 50 );
  this->SetUseMultiCylindersForBonds( false );
  this->SetBondRadius( 0.075 );
}

//----------------------------------------------------------------------------
vtkLookupTable *vtkMoleculeMapper2::GetAtomLookupTable()
{
  return
    dynamic_cast<vtkLookupTable*>(this->AtomGlyphMapper->GetLookupTable());
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::SetAtomLookupTable(vtkLookupTable *lut)
{
  if (dynamic_cast<vtkLookupTable*>(this->AtomGlyphMapper->GetLookupTable()) == lut)
    {
    return;
    }
  this->AtomGlyphMapper->SetLookupTable(lut);
  this->Modified();
}

//----------------------------------------------------------------------------
vtkLookupTable *vtkMoleculeMapper2::GetBondLookupTable()
{
  return
    dynamic_cast<vtkLookupTable*>(this->BondGlyphMapper->GetLookupTable());
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::SetBondLookupTable(vtkLookupTable *lut)
{
  if (dynamic_cast<vtkLookupTable*>(this->BondGlyphMapper->GetLookupTable()) == lut)
    {
    return;
    }
  this->BondGlyphMapper->SetLookupTable(lut);
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::LookupTableModified()
{
  this->Modified();
}

//----------------------------------------------------------------------------
const char * vtkMoleculeMapper2::GetAtomicRadiusTypeAsString()
{
  switch (this->AtomicRadiusType)
    {
    case CovalentRadius:
      return "CovalentRadius";
    case VDWRadius:
      return "VDWRadius";
    case UnitRadius:
      return "UnitRadius";
    default:
      return "Invalid";
    }
}

//----------------------------------------------------------------------------
const char * vtkMoleculeMapper2::GetBondColorModeAsString()
{
  switch (this->BondColorMode)
    {
    case SingleColor:
      return "SingleColor";
    case DiscreteByAtom:
      return "DiscreteByAtom";
    default:
      return "Invalid";
    }
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::GetSelectedAtomsAndBonds(vtkSelection *selection,
                                                 vtkIdTypeArray *atomIds,
                                                 vtkIdTypeArray *bondIds)
{
  // Sanity check
  if (selection == NULL || (atomIds == NULL && bondIds == NULL) )
    {
    return;
    }

  // Clear the inputs
  if (atomIds != NULL)
    {
    atomIds->Reset();
    }
  if (bondIds != NULL)
    {
    bondIds->Reset();
    }

  const vtkIdType numAtoms = this->GetInput()->GetNumberOfAtoms();
  const vtkIdType numBonds = this->GetInput()->GetNumberOfBonds();
  const vtkIdType numAtomsAndBonds = numAtoms + numBonds;

  // Find selection node that we're interested in:
  const vtkIdType numNodes = selection->GetNumberOfNodes();
  for (vtkIdType nodeId = 0; nodeId < numNodes; ++nodeId)
    {
    vtkSelectionNode *node = selection->GetNode(nodeId);

    // Check if the mapper is this instance of MoleculeMapper2
    vtkActor *selActor = vtkActor::SafeDownCast(
               node->GetProperties()->Get(vtkSelectionNode::PROP()));
    if (selActor && (selActor->GetMapper() == this))
      {
      // Separate the selection ids into atoms and bonds
      vtkIdTypeArray *selIds = vtkIdTypeArray::SafeDownCast(
            node->GetSelectionList());
      if (selIds)
        {
        vtkIdType numIds = selIds->GetNumberOfTuples();
        for (vtkIdType i = 0; i < numIds; ++i)
          {
          vtkIdType curId = selIds->GetValue(i);
          if (atomIds != NULL && curId < numAtoms) // atoms
            {
            atomIds->InsertNextValue(curId);
            }
          else if (bondIds != NULL && curId < numAtomsAndBonds)// bonds
            {
            // Remove offset
            curId -= numAtoms;
            bondIds->InsertNextValue(curId);
            }
          }
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::Render(vtkRenderer *ren, vtkActor *act )
{
  // If we add more rendering backend (e.g. point sprites), add a switch here.
  this->GlyphRender(ren, act);
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::GlyphRender(vtkRenderer *ren, vtkActor *act)
{
  // Update cached polydata if needed
  this->UpdateGlyphPolyData();

  // update POVRay output if needed
  if (this->POVRayStreaming)
    {
    this->POVRaySStream.str("");

    POVRayCamera cam(ren->GetActiveCamera());
    this->POVRaySStream << cam << endl;

    vtkLight *l;
    vtkLightCollection *lights = ren->GetLights();
    lights->InitTraversal();
    while ((l = lights->GetNextItem()))
      {
      POVRayLight light(l);
      this->POVRaySStream << light << endl;
      }

    this->WriteColorDefs(this->POVRaySStream);

     this->POVRaySStream << "#declare Molecule = union {" << endl;
     if (this->RenderAtoms)
       {
       this->WriteAtoms(this->POVRaySStream);
       }
     if (this->RenderBonds)
       {
       this->WriteBonds(this->POVRaySStream);
       }
    this->POVRaySStream << "}" << endl;
    }

  // Pass rendering call on
  if (this->RenderAtoms)
    {
    this->AtomGlyphMapper->Render(ren, act);
    }

  if (this->RenderBonds)
    {
    this->BondGlyphMapper->Render(ren, act);
    }
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::UpdateGlyphPolyData()
{
  vtkMolecule *molecule = this->GetInput();

  if (!this->GlyphDataInitialized || (
        (molecule->GetMTime() > this->AtomGlyphPolyData->GetMTime() ||
         this->GetMTime() > this->AtomGlyphPolyData->GetMTime()) &&
        this->RenderAtoms))
    {
    this->UpdateAtomGlyphPolyData();
    }

  if (!this->GlyphDataInitialized || (
        (molecule->GetMTime() > this->BondGlyphPolyData->GetMTime() ||
         this->GetMTime() > this->BondGlyphPolyData->GetMTime()) &&
        this->RenderBonds))
    {
    this->UpdateBondGlyphPolyData();
    }

  this->GlyphDataInitialized = true;
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::WriteColorDefs(ostream &os)
{
  // finishes
  this->AtomFinish
    = POVRayFinish("AtomFinish",0.1, 0.55, 1, 0.001, 0.75, 1, 0.85, true);

  os << this->AtomFinish << endl;

  this->BondFinish
    = POVRayFinish("BondFinish",0.1, 0.1, 1, 0.001, 0.75, 1, 0.95, true);

  os << this->BondFinish << endl;

  // for fixed bond color
  this->UniformBondColor
    = POVRayColor("BondColor", this->BondColor[0], this->BondColor[1], this->BondColor[2]);

  os << this->UniformBondColor << endl;

  // atom based colors
  vtkMolecule *molecule = this->GetInput();
  const vtkIdType numAtoms = molecule->GetNumberOfAtoms();
  vtkUnsignedShortArray *atomicNums = molecule->GetAtomicNumberArray();

  vtkLookupTable *lut
    = dynamic_cast<vtkLookupTable*>(this->AtomGlyphMapper->GetLookupTable());

  this->AtomColors.clear();
  for (vtkIdType i = 0; i < numAtoms; ++i)
    {
    int atomicNum = atomicNums->GetValue(i);
    double mappedColor[3];
    lut->GetColor(static_cast<double>(atomicNum), mappedColor);
    ostringstream oss;
    oss << "Atom" << atomicNum << "Color";
    POVRayColor atomColor(oss.str(), mappedColor[0], mappedColor[1], mappedColor[2]);
    pair<int,POVRayColor> val(atomicNum,atomColor);
    pair<map<int,POVRayColor>::iterator,bool> insval = this->AtomColors.insert(val);
    if (insval.second)
      {
      os << atomColor << endl;
      }
    }

}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::WriteBonds(ostream &os)
{

  vtkMolecule *molecule = this->GetInput();
  const vtkIdType numBonds = molecule->GetNumberOfBonds();

  // Generate the scale, orientation, and position of each cylinder
  for (vtkIdType bondInd = 0; bondInd<numBonds; ++bondInd)
    {
    vtkBond bond = molecule->GetBond(bondInd);
    unsigned short bondOrder = bond.GetOrder();

    vtkVector3f pos1 = bond.GetBeginAtom().GetPosition();
    vtkVector3f pos2 = bond.GetEndAtom().GetPosition();

    vtkVector3f bondVec = pos2 - pos1;
    bondVec.Normalize();

    vtkVector3f pos3(pos1);
    pos3 += pos2;
    pos3 *= 0.5f;

    unsigned short atomicNumbers[2];
    atomicNumbers[0] = bond.GetBeginAtom().GetAtomicNumber();
    atomicNumbers[1] = bond.GetEndAtom().GetAtomicNumber();

    vtkVector3f delta;
    delta.Set(0.0, 0.0, 0.0);

    vtkVector3f initialDisp;
    initialDisp.Set(0.0, 0.0, 0.0);

    // Set up delta step vector and bond radius from bond order:
    const static vtkVector3f unitZ(0.0f, 0.0f, 1.0f);
    if (this->UseMultiCylindersForBonds)
      {
      const float deltaLength = this->BondRadius * 2.6;
      switch (bondOrder)
        {
        case 1:
        default:
          delta.Set(0.0, 0.0, 0.0);
          initialDisp.Set(0.0, 0.0, 0.0);
          break;
        case 2:
          delta = bondVec.Cross(unitZ).Normalized();
          delta[0] *= deltaLength;
          delta[1] *= deltaLength;
          delta[2] *= deltaLength;
          initialDisp.Set(delta[0]*(-0.5), delta[1]*(-0.5), delta[2]*(-0.5));
          break;
        case 3:
          delta = bondVec.Cross(unitZ).Normalized();
          delta[0] *= deltaLength;
          delta[1] *= deltaLength;
          delta[2] *= deltaLength;
          initialDisp.Set(-delta[0], -delta[1], -delta[2]);
          break;
        }
      pos1[0] += initialDisp[0];
      pos1[1] += initialDisp[1];
      pos1[2] += initialDisp[2];

      pos2[0] += initialDisp[0];
      pos2[1] += initialDisp[1];
      pos2[2] += initialDisp[2];

      pos3[0] += initialDisp[0];
      pos3[1] += initialDisp[1];
      pos3[2] += initialDisp[2];
      }

    // For each bond order, add a point to the glyph points, translate
    // by delta, and repeat.
    for (unsigned short iter = 0; iter<bondOrder; ++iter)
      {
      // Single color mode adds a single cylinder, while
      // DiscreteByAtom adds two differently colored and positioned
      // cylinders.
      switch (this->BondColorMode)
        {
        case SingleColor:
          {
          POVRayCylinder bond(pos1, pos2, this->BondRadius, this->UniformBondColor, this->BondFinish);
          os << bond << endl;
          }
          break;
        default:
        case DiscreteByAtom:
          {
          POVRayCylinder bond0(pos1, pos3, this->BondRadius, this->AtomColors[atomicNumbers[0]], this->BondFinish);
          POVRayCylinder bond1(pos3, pos2, this->BondRadius, this->AtomColors[atomicNumbers[1]], this->BondFinish);
          os << bond0 << endl << bond1 << endl;
          }
        }
      // Prepare for next multicylinder
      if (this->UseMultiCylindersForBonds && bondOrder != 1)
        {
        pos1 += delta;
        pos2 += delta;
        pos3 += delta;
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::WriteAtoms(ostream &os)
{
  vtkMolecule *molecule = this->GetInput();
  const vtkIdType numAtoms = molecule->GetNumberOfAtoms();
  vtkPoints *atomCenters = molecule->GetAtomicPositionArray();
  vtkUnsignedShortArray *atomicNums = molecule->GetAtomicNumberArray();

  // generate the atoms
  for (vtkIdType i = 0; i < numAtoms; ++i)
    {
    int atomicNum = atomicNums->GetValue(i);

    float radius = 1.0f;
    switch (this->AtomicRadiusType)
      {
      default:
        vtkWarningMacro(<<"Unknown radius type: " << this->AtomicRadiusType
                        <<". Falling back to 'VDWRadius' ("<<VDWRadius<<").");
      case VDWRadius:
        radius
          = this->AtomicRadiusScaleFactor
          * this->PeriodicTable->GetVDWRadius(atomicNums->GetValue(i));
        break;

      case CovalentRadius:
        radius
          = this->AtomicRadiusScaleFactor
          * this->PeriodicTable->GetCovalentRadius(atomicNums->GetValue(i));
        break;
      case UnitRadius:
        radius = this->AtomicRadiusScaleFactor;
        break;
      }

    ostringstream oss;
    oss << "Atom" << atomicNum;

    // TODO -- template on pt type
    double pt[3];
    atomCenters->GetPoint(i, pt);
    vtkVector3f pos(pt[0], pt[1], pt[2]);

    POVRaySphere atom(oss.str(), pos, radius, this->AtomColors[atomicNum], this->AtomFinish);
    os << atom << endl;
    }
}

//----------------------------------------------------------------------------
// Generate scale and position information for each atom sphere
void vtkMoleculeMapper2::UpdateAtomGlyphPolyData()
{
  this->AtomGlyphPolyData->Initialize();

  vtkMolecule *molecule = this->GetInput();
  const vtkIdType numAtoms = molecule->GetNumberOfAtoms();

  vtkUnsignedShortArray *atomicNums = molecule->GetAtomicNumberArray();
  this->AtomGlyphPolyData->GetPointData()->AddArray(atomicNums);
  this->AtomGlyphPolyData->SetPoints(molecule->GetAtomicPositionArray());
  this->AtomGlyphMapper->SelectColorArray("Atomic Numbers");

  vtkNew<vtkFloatArray> scaleFactors;
  scaleFactors->SetNumberOfComponents(1);
  scaleFactors->SetName("Scale Factors");
  scaleFactors->Allocate(numAtoms);

  switch (this->AtomicRadiusType)
    {
    default:
      vtkWarningMacro(<<"Unknown radius type: " << this->AtomicRadiusType
                      <<". Falling back to 'VDWRadius' ("<<VDWRadius<<").");
    case VDWRadius:
      for (vtkIdType i = 0; i < numAtoms; ++i)
        {
        scaleFactors->InsertNextValue(this->AtomicRadiusScaleFactor *
          this->PeriodicTable->GetVDWRadius(atomicNums->GetValue(i)));
        }
      break;
    case CovalentRadius:
      for (vtkIdType i = 0; i < numAtoms; ++i)
        {
        scaleFactors->InsertNextValue(this->AtomicRadiusScaleFactor *
          this->PeriodicTable->GetCovalentRadius(atomicNums->GetValue(i)));
        }
      break;
    case UnitRadius:
      for (vtkIdType i = 0; i < numAtoms; ++i)
        {
        scaleFactors->InsertNextValue(this->AtomicRadiusScaleFactor);
        }
      break;
    }

  this->AtomGlyphPolyData->GetPointData()->AddArray(scaleFactors.GetPointer());
  this->AtomGlyphMapper->SetScaleArray("Scale Factors");
}

//----------------------------------------------------------------------------
// Generate position, scale, and orientation vectors for each bond cylinder
void vtkMoleculeMapper2::UpdateBondGlyphPolyData()
{
  this->BondGlyphPolyData->Initialize();

  vtkMolecule *molecule = this->GetInput();
  const vtkIdType numBonds = molecule->GetNumberOfBonds();

  // For selection ID offset:
  const vtkIdType numAtoms = molecule->GetNumberOfAtoms();

  // Create arrays
  vtkNew<vtkPoints> cylCenters;
  vtkNew<vtkFloatArray> cylScales;
  vtkNew<vtkFloatArray> orientationVectors;
  // Since vtkHardwareSelector won't distinguish between the internal instance
  // mappers of vtkMoleculeMapper2, use a custom selection ID range. This also
  // fixes the issue of bonds that are colored-by-atom, as these are rendered
  // as two glyphs.
  vtkNew<vtkIdTypeArray> selectionIds;

  // Setup arrays
  // vtkPoints uses three components by default
  cylScales->SetNumberOfComponents(3);
  orientationVectors->SetNumberOfComponents(3);
  selectionIds->SetNumberOfComponents(1);

  // Name arrays
  // Can't set name for points
  cylScales->SetName("Scale Factors");
  orientationVectors->SetName("Orientation Vectors");
  selectionIds->SetName("Selection Ids");

  // Allocate memory -- find out how many cylinders are needed
  vtkIdType numCylinders = numBonds;
  // Up to three cylinders per bond if multicylinders are enabled:
  if (this->UseMultiCylindersForBonds)
    {
    numCylinders *= 3;
    }
  // If DiscreteByAtom coloring is used, each cylinder is represented
  // by two individual cylinders
  if (this->BondColorMode == DiscreteByAtom)
    {
    numCylinders *= 2;
    }

  // Allocate memory. Multiply numCylinders by number of components in array.
  cylCenters->Allocate(3*numCylinders);
  cylScales->Allocate(3*numCylinders);
  orientationVectors->Allocate(3*numCylinders);
  selectionIds->Allocate(numCylinders);

  // Add arrays to BondGlyphPolyData
  this->BondGlyphPolyData->SetPoints(cylCenters.GetPointer());
  this->BondGlyphPolyData->GetPointData()->AddArray(cylScales.GetPointer());
  this->BondGlyphPolyData->GetPointData()->AddArray(
    orientationVectors.GetPointer());
  this->BondGlyphPolyData->GetPointData()->AddArray(
        selectionIds.GetPointer());

  // Set up coloring mode
  vtkDataArray *cylColors = 0;
  switch(this->BondColorMode)
    {
    case SingleColor:
      cylColors = vtkUnsignedCharArray::New();
      cylColors->SetNumberOfComponents(3);
      cylColors->Allocate(3*numCylinders);
      cylColors->SetName("Colors");
      this->BondGlyphPolyData->GetPointData()->SetScalars(cylColors);
      this->BondGlyphMapper->SetColorModeToDefault();
      this->BondGlyphMapper->SetScalarModeToUsePointData();
      break;
    default:
    case DiscreteByAtom:
      cylColors = vtkUnsignedShortArray::New();
      cylColors->SetNumberOfComponents(1);
      cylColors->Allocate(numCylinders);
      cylColors->SetName("Colors");
      this->BondGlyphPolyData->GetPointData()->SetScalars(cylColors);
      vtkNew<vtkLookupTable> lut;
      lut->DeepCopy(this->AtomGlyphMapper->GetLookupTable());
      //this->PeriodicTable->GetDefaultLUT(lut.GetPointer());
      this->BondGlyphMapper->SetLookupTable(lut.GetPointer());
      //is->BondGlyphMapper->SetLookupTable(this->AtomGlyphMapper->GetLookupTable());
      this->BondGlyphMapper->SetScalarRange
        (0, this->PeriodicTable->GetNumberOfElements());
      this->BondGlyphMapper->SetScalarModeToUsePointData();
      this->AtomGlyphMapper->SetColorModeToMapScalars();
      break;
    }

  // Set up pointers to the specific color arrays
  vtkUnsignedCharArray *singleColorArray =
    vtkUnsignedCharArray::SafeDownCast(cylColors);
  vtkUnsignedShortArray *discreteColorArray =
    vtkUnsignedShortArray::SafeDownCast(cylColors);

  // Remove reference to cylColors
  cylColors->Delete();

  // Declare some variables for later
  unsigned short bondOrder;
  float bondLength;
  // Since the input cylinder is oriented along the the z axis, the
  // scale vector should be [radius, radius, bondLength]
  vtkVector3f scale;
  // Current cylinder's selection id
  vtkIdType selectionId;
  // Distance between multicylinder surfaces is approx. 1/3 of the
  // diameter:
  const float deltaLength = this->BondRadius * 2.6;
  // Vector between centers cylinders in a multibond:
  vtkVector3f delta;
  delta.Set(0.0, 0.0, 0.0);
  // The initial displacement when generating a multibond
  vtkVector3f initialDisp;
  initialDisp.Set(0.0, 0.0, 0.0);
  // The geometric center of the bond
  vtkVector3f bondCenter;
  // The center of the current cylinder
  vtkVector3f cylinderCenter;
  // Used in DiscreteByAtom color mode:
  vtkVector3f halfCylinderCenter;
  // The positions of the start and end atoms
  vtkVector3f pos1, pos2;
  // Array containing the atomic numbers of {begin, end} atoms
  unsigned short atomicNumbers[2];
  // Normalized vector pointing along bond from begin->end atom
  vtkVector3f bondVec;
  // Unit z vector -- used for multicylinder orientation
  const static vtkVector3f unitZ (0.0, 0.0, 1.0);
  // Can't use InsertNextTuple(unsigned char *) in a
  // vtkUnsignedCharArray. This float array is used instead.
  // Initialize with bondColor for SingleColor mode
  vtkColor3f bondColorf (static_cast<float>(this->BondColor[0]),
                         static_cast<float>(this->BondColor[1]),
                         static_cast<float>(this->BondColor[2]));

  // Generate the scale, orientation, and position of each cylinder
  for (vtkIdType bondInd = 0; bondInd < numBonds; ++bondInd)
    {
    selectionId = numAtoms + bondInd; // mixing 1 and 0 indexed ids on purpose
    // Extract bond info
    vtkBond bond = molecule->GetBond(bondInd);
    bondOrder = bond.GetOrder();
    pos1 = bond.GetBeginAtom().GetPosition();
    pos2 = bond.GetEndAtom().GetPosition();
    atomicNumbers[0] = bond.GetBeginAtom().GetAtomicNumber();
    atomicNumbers[1] = bond.GetEndAtom().GetAtomicNumber();

    // Compute additional bond info
    // - Normalized vector in direction of bond
    bondVec = pos2 - pos1;
    bondLength = bondVec.Normalize();
    // - Center of bond for translation
    // TODO vtkVector scalar multiplication
//    bondCenter = (pos1 + pos2) * 0.5;
    bondCenter[0] = (pos1[0] + pos2[0]) * 0.5;
    bondCenter[1] = (pos1[1] + pos2[1]) * 0.5;
    bondCenter[2] = (pos1[2] + pos2[2]) * 0.5;
    // end vtkVector TODO

    // Set up delta step vector and bond radius from bond order:
    if (this->UseMultiCylindersForBonds)
      {
      switch (bondOrder)
        {
        case 1:
        default:
          delta.Set(0.0, 0.0, 0.0);
          initialDisp.Set(0.0, 0.0, 0.0);
          break;
        case 2:
          // TODO vtkVector scalar multiplication
//          delta = bondVec.Cross(unitZ).Normalized() * deltaLength;
//          initialDisp = delta * (-0.5);
          delta = bondVec.Cross(unitZ).Normalized();
          delta[0] *= deltaLength;
          delta[1] *= deltaLength;
          delta[2] *= deltaLength;
          initialDisp.Set(delta[0]*(-0.5), delta[1]*(-0.5), delta[2]*(-0.5));
          // End vtkVector TODO
          break;
        case 3:
          // TODO vtkVector scalar multiplication, negation
//          delta = bondVec.Cross(unitZ).Normalized() * deltaLength;
//          initialDisp = -delta;
          delta = bondVec.Cross(unitZ).Normalized();
          delta[0] *= deltaLength;
          delta[1] *= deltaLength;
          delta[2] *= deltaLength;
          initialDisp.Set(-delta[0], -delta[1], -delta[2]);
          // End vtkVector TODO
          break;
        }
      }

    // Set up cylinder scale factors
    switch (this->BondColorMode)
      {
      case SingleColor:
        scale.Set( bondLength, this->BondRadius, this->BondRadius);
        break;
      default:
      case DiscreteByAtom:
        scale.Set( 0.5 * bondLength, this->BondRadius, this->BondRadius);
        break;
      }

    if (this->UseMultiCylindersForBonds)
      {
      cylinderCenter = bondCenter + initialDisp;
      }
    else
      {
      cylinderCenter = bondCenter;
      }

    // For each bond order, add a point to the glyph points, translate
    // by delta, and repeat.
    for (unsigned short iter = 0; iter < bondOrder; ++iter)
      {
      // Single color mode adds a single cylinder, while
      // DiscreteByAtom adds two differently colored and positioned
      // cylinders.
      switch (this->BondColorMode)
        {
        case SingleColor:
          cylCenters->InsertNextPoint(cylinderCenter.GetData());
          cylScales->InsertNextTuple(scale.GetData());
          singleColorArray->InsertNextTuple(bondColorf.GetData());
          orientationVectors->InsertNextTuple(bondVec.GetData());
          selectionIds->InsertNextValue(selectionId);
          break;
        default:
        case DiscreteByAtom:
          // Cache some scaling factors
          const float quarterLength = 0.25 * bondLength;

          // Add cylinder for begin atom:
          // TODO vtkVector subtraction, scalar multiplication
//          halfCylinderCenter = cylinderCenter - (bondVec * quarterLength);
          halfCylinderCenter[0] = cylinderCenter[0] -
              (bondVec[0] * quarterLength);
          halfCylinderCenter[1] = cylinderCenter[1] -
              (bondVec[1] * quarterLength);
          halfCylinderCenter[2] = cylinderCenter[2] -
              (bondVec[2] * quarterLength);
          // end vtkVector TODO

          cylCenters->InsertNextPoint(halfCylinderCenter.GetData());
          cylScales->InsertNextTuple(scale.GetData());
          discreteColorArray->InsertNextValue(atomicNumbers[0]);
          orientationVectors->InsertNextTuple(bondVec.GetData());
          selectionIds->InsertNextValue(selectionId);

          // Add cylinder for begin atom:
          // TODO vtkVector addition, scalar multiplication
//          halfCylinderCenter = cylinderCenter + (bondVec * quarterLength);
          halfCylinderCenter[0] = cylinderCenter[0] +
              (bondVec[0] * quarterLength);
          halfCylinderCenter[1] = cylinderCenter[1] +
              (bondVec[1] * quarterLength);
          halfCylinderCenter[2] = cylinderCenter[2] +
              (bondVec[2] * quarterLength);
          // end vtkVector TODO

          cylCenters->InsertNextPoint(halfCylinderCenter.GetData());
          cylScales->InsertNextTuple(scale.GetData());
          discreteColorArray->InsertNextValue(atomicNumbers[1]);
          orientationVectors->InsertNextTuple(bondVec.GetData());
          selectionIds->InsertNextValue(selectionId);
        }

      // Prepare for next multicylinder
      if (this->UseMultiCylindersForBonds && bondOrder != 1)
        {
        // TODO vtkVector in-place addition
//        cylinderCenter += delta;
        cylinderCenter[0] += delta[0];
        cylinderCenter[1] += delta[1];
        cylinderCenter[2] += delta[2];
        // end vtkVector TODO
        }
      }
    }

  // Free up some space
  this->BondGlyphPolyData->Squeeze();

  // Setup glypher
  this->BondGlyphMapper->SetScaleArray("Scale Factors");
  this->BondGlyphMapper->SetOrientationArray("Orientation Vectors");
  this->BondGlyphMapper->SetSelectionIdArray("Selection Ids");
  this->BondGlyphMapper->UseSelectionIdsOn();
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::ReleaseGraphicsResources(vtkWindow *w)
{
  this->AtomGlyphMapper->ReleaseGraphicsResources(w);
  this->BondGlyphMapper->ReleaseGraphicsResources(w);
}

double *vtkMoleculeMapper2::GetBounds()
{
  vtkMolecule *input = this->GetInput();
  if (!input)
    {
    vtkMath::UninitializeBounds(this->Bounds);
    }
  else
    {
    if (!this->Static)
      {
      this->Update();
      }
    input->GetBounds(this->Bounds);
    // Pad bounds by 3 Angstrom to contain spheres, etc
    this->Bounds[0] -= 3.0;
    this->Bounds[1] += 3.0;
    this->Bounds[2] -= 3.0;
    this->Bounds[3] += 3.0;
    this->Bounds[4] -= 3.0;
    this->Bounds[5] += 3.0;
    }
  return this->Bounds;
}

//----------------------------------------------------------------------------
int vtkMoleculeMapper2::FillInputPortInformation(int vtkNotUsed(port),
                                                vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMolecule");
  return 1;
}

//----------------------------------------------------------------------------
void vtkMoleculeMapper2::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "AtomGlyphMapper:\n";
  this->AtomGlyphMapper->PrintSelf(os, indent.GetNextIndent());

  os << indent << "BondGlyphMapper:\n";
  this->BondGlyphMapper->PrintSelf(os, indent.GetNextIndent());
}
