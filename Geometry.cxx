#include "Geometry.h"
#include "Utility.h"
#include "AtomicProperties.h"
#include "Math3D.h"
#include "BondDetector.h"
#include "CoordinationPolyhedraGenerator.h"

#include <vtkMolecule.h>
#include <vtkPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkStringArray.h>
#include <vtkDataSetAttributes.h>
#include <vtkInEdgeIterator.h>
#include <vtkVertexListIterator.h>
#include <vtkGraphEdge.h>
#include <vtkPointData.h>

#include <vtkPolyDataWriter.h>

#include <vector>
#include <string>
#include <map>
#include <utility>

using std::vector;
using std::string;
using std::map;
using std::pair;

// ---------------------------------------------------------------------------
void BuildAxes(
    vtkPolyData *data,
    const vector<double> &axes)
{
  if (axes.size() != 9)
    {
    pError(cerr) << "invalid axes points" << endl;
    return;
    }

  vtkDoubleArray *pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents(3);
  pts->SetNumberOfTuples(4);
  double *pPts = pts->GetPointer(0);
  pPts[0] = 0.0;
  pPts[1] = 0.0;
  pPts[2] = 0.0;
  pPts += 3;
  for (size_t i = 0; i<9; ++i)
    {
    pPts[i] = axes[i];
    }

  vtkPoints *points = vtkPoints::New();
  points->SetData(pts);
  pts->Delete();

  data->SetPoints(points);
  points->Delete();

  vtkIdTypeArray *ptIds = vtkIdTypeArray::New();
  ptIds->SetNumberOfTuples(9);
  vtkIdType *pPtIds = ptIds->GetPointer(0);
  pPtIds[0] = 2;
  pPtIds[1] = 0;
  pPtIds[2] = 1;
  pPtIds[3] = 2;
  pPtIds[4] = 0;
  pPtIds[5] = 2;
  pPtIds[6] = 2;
  pPtIds[7] = 0;
  pPtIds[8] = 3;

  vtkCellArray *cells = vtkCellArray::New();
  cells->SetCells(3, ptIds);
  ptIds->Delete();

  data->SetLines(cells);

  vtkUnsignedCharArray *colors = vtkUnsignedCharArray::New();
  colors->SetNumberOfComponents(3);
  colors->SetNumberOfTuples(3);
  colors->SetName("colors");
  unsigned char *pColors = colors->GetPointer(0);
  memset(pColors, 0, 9);
  pColors[0] = 255;
  pColors[4] = 255;
  pColors[8] = 255;
  data->GetCellData()->AddArray(colors);
  data->GetCellData()->SetActiveScalars("colors");
  colors->Delete();
}

// ---------------------------------------------------------------------------
void BuildMolecule(
    vtkMolecule *molecule,
    vector<double> &positions,
    vector<unsigned short> &numbers,
    vector<string> &labels,
    vector<int> &labelIds,
    vector<bool> &sites,
    vector<bool> &coordinationSites,
    vector<bool> &ghosts,
    const BondDetector *detector)
{
  molecule->Initialize();

  // atoms + labels
  unsigned long natoms = numbers.size();

  vtkStringArray *atomLabels = vtkStringArray::New();
  atomLabels->SetNumberOfTuples(natoms);
  atomLabels->SetName("labels");
  molecule->GetVertexData()->AddArray(atomLabels);
  atomLabels->Delete();

  vtkIntArray *atomLabelIds = vtkIntArray::New();
  atomLabelIds->SetNumberOfTuples(natoms);
  atomLabelIds->SetName("label ids");
  molecule->GetVertexData()->AddArray(atomLabelIds);
  atomLabelIds->Delete();

  vtkIntArray *atomSites = vtkIntArray::New();
  atomSites->SetNumberOfTuples(natoms);
  atomSites->SetName("sites");
  molecule->GetVertexData()->AddArray(atomSites);
  atomSites->Delete();

  vtkIntArray *atomCoordinationSites = vtkIntArray::New();
  atomCoordinationSites->SetNumberOfTuples(natoms);
  atomCoordinationSites->SetName("coordination sites");
  molecule->GetVertexData()->AddArray(atomCoordinationSites);
  atomCoordinationSites->Delete();

  vtkIntArray *atomGhosts = vtkIntArray::New();
  atomGhosts->SetNumberOfTuples(natoms);
  atomGhosts->SetName("ghosts");
  molecule->GetVertexData()->AddArray(atomGhosts);
  atomGhosts->Delete();

  for (unsigned long i=0; i<natoms; ++i)
    {
    unsigned long ii = 3*i;

    vtkAtom atom = molecule->AppendAtom(
        numbers[i], positions[ii], positions[ii+1], positions[ii+2]);

    vtkIdType id = atom.GetId();

    atomLabels->SetValue(id, labels[i].c_str());
    atomLabelIds->SetValue(id, labelIds[i]);
    atomSites->SetValue(id, static_cast<int>(sites[i]));
    atomCoordinationSites->SetValue(id, static_cast<int>(coordinationSites[i]));
    atomGhosts->SetValue(id, static_cast<int>(ghosts[i]));
    }

  // bonds
  for (unsigned int i=0; i<natoms; ++i)
    {
    unsigned int ii = 3*i;
    vector<unsigned long> bonds;
    for (unsigned int j=i+1; j<natoms; ++j)
      {
      unsigned int jj = 3*j;
      if (detector->Connected(numbers[i], &positions[ii], numbers[j], &positions[jj]))
        {
        bonds.push_back(j);
        }
      }
    unsigned long nTargets = bonds.size();
    for (unsigned long j=0; j<nTargets; ++j)
      {
      molecule->AppendBond(i, bonds[j]);
      }
    }
}

// ---------------------------------------------------------------------------
void CopyActiveSites(vtkMolecule *source, bool ghostBonds, vtkMolecule *dest)
{
  dest->Initialize();

  vtkIntArray *sites
    = dynamic_cast<vtkIntArray*>(source->GetVertexData()->GetArray("sites"));
  if (!sites)
    {
    pError(cerr) << "no sites array" << endl;
    return;
    }
  int *pSites = sites->GetPointer(0);

  vtkIntArray *ghosts
    = dynamic_cast<vtkIntArray*>(source->GetVertexData()->GetArray("ghosts"));
  if (!ghosts)
    {
    pError(cerr) << "no ghost array" << endl;
    return;
    }
  int *pGhosts = ghosts->GetPointer(0);

  // map used to track which atoms have been coppied
  // and prevent duplicates
  map<vtkIdType, vtkIdType> atomMap;

  // traverse all bonds
  vtkIdType nBonds = source->GetNumberOfBonds();
  for (vtkIdType i = 0; i< nBonds; ++i)
    {
    vtkBond bond = source->GetBond(i);
    vtkIdType sid = bond.GetBeginAtomId();
    vtkIdType did = bond.GetEndAtomId();
    if ((pSites[sid] && pSites[did]))
      {
      // both sites active
      if ( (!pGhosts[sid] && !pGhosts[did])
        || (ghostBonds && (!pGhosts[sid] &&  pGhosts[did])) )
        {
        // copy source atom
        vtkIdType oSid = -1;
        if (!atomMap.count(sid))
          {
          vtkAtom a = source->GetAtom(sid);
          vtkAtom b = dest->AppendAtom(a.GetAtomicNumber(), a.GetPosition());
          oSid = b.GetId();
          atomMap.insert(pair<vtkIdType,vtkIdType>(sid, oSid));
          }
        else
          {
          oSid = atomMap[sid];
          }
        // copy dest atom
        vtkIdType oDid = -1;
        if (!atomMap.count(did))
          {
          vtkAtom a = source->GetAtom(did);
          vtkAtom b = dest->AppendAtom(a.GetAtomicNumber(), a.GetPosition());
          oDid = b.GetId();
          atomMap.insert(pair<vtkIdType,vtkIdType>(did, oDid));
          }
        else
          {
          oDid = atomMap[did];
          }
        // copy bond
        dest->AppendBond(oSid, oDid, bond.GetOrder());
        }
      }
    }
}

// ---------------------------------------------------------------------------
void BuildPolyhedra(vtkMolecule *molecule, vtkPolyData *polyhedra)
{
  vtkIntArray *sites
    = dynamic_cast<vtkIntArray*>(molecule->GetVertexData()->GetArray("coordination sites"));
  if (!sites)
    {
    pError(cerr) << "no sites array" << endl;
    return;
    }
  int *pSites = sites->GetPointer(0);

  vtkIntArray *ghosts
    = dynamic_cast<vtkIntArray*>(molecule->GetVertexData()->GetArray("ghosts"));
  if (!ghosts)
    {
    pError(cerr) << "no ghost array" << endl;
    return;
    }
  int *pGhosts = ghosts->GetPointer(0);

  vtkIntArray *labels
    = dynamic_cast<vtkIntArray*>(molecule->GetVertexData()->GetArray("label ids"));
  if (!labels)
    {
    pError(cerr) << "no labels array" << endl;
    return;
    }
  int *pLabels = labels->GetPointer(0);

  CoordinationPolyhedraGenerator *generator = new CoordinationPolyhedraGenerator;
  generator->SetMinNumberOfPoints(7);

  vtkPoints *sitePoints = molecule->GetPoints();

  vtkInEdgeIterator *edgeIter = vtkInEdgeIterator::New();

  vtkVertexListIterator *vertIter = vtkVertexListIterator::New();
  vertIter->SetGraph(molecule);

  vector<double> points;
  vector<vtkIdType> ids;

  // visit each site
  while (vertIter->HasNext())
    {
    vtkIdType vertId = vertIter->Next();

    // skip inactive sites and ghost nodes
    if (!pSites[vertId] || pGhosts[vertId])
      {
      continue;
      }

    molecule->GetInEdges(vertId, edgeIter);

    points.clear();
    ids.clear();

    double pt[3];
    sitePoints->GetPoint(vertId, pt);

    points.push_back(pt[0]);
    points.push_back(pt[1]);
    points.push_back(pt[2]);

    ids.push_back(vertId);

    // accumulate sites bonded to
    while (edgeIter->HasNext())
      {
      vtkGraphEdge *edge = edgeIter->NextGraphEdge();
      vtkIdType sourceId = edge->GetSource();

      sitePoints->GetPoint(sourceId, pt);

      points.push_back(pt[0]);
      points.push_back(pt[1]);
      points.push_back(pt[2]);

      ids.push_back(sourceId);
      }

    // generate coordination polyhedra
    if (generator->Insert(pLabels[vertId], points, ids))
      {
      // flip the ghost flag on ids that are in
      // a coordination polyhedra
      size_t nIds = ids.size();
      for (size_t i =0; i<nIds; ++i)
        {
        pGhosts[ids[i]] = false;
        }
      }
    }

  // build the output
  polyhedra->Initialize();
  polyhedra->SetPolys(generator->GetTriangles());
  polyhedra->SetPoints(generator->GetPoints());
  polyhedra->GetPointData()->SetScalars(generator->GetLabels());

  polyhedra->Print(cerr);

  #ifdef GeometryDEBUG
  vtkPolyDataWriter *w = vtkPolyDataWriter::New();
  w->SetInputData(polyhedra);
  w->SetFileName("polyhedra.vtk");
  w->Write();
  w->Delete();
  #endif

  delete generator;
}
