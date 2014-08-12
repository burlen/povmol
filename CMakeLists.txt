cmake_minimum_required(VERSION 2.8)
PROJECT(2pov)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(Qt4 REQUIRED)
include(${QT_USE_FILE})

set(CXX_SOURCES
  AtomicProperties.cxx
  AtomicPropertiesBondDetector.cxx
  Geometry.cxx
  Math3D.cxx
  PointTree.cxx
  POVRayPrimatives.cxx
  Utility.cxx
  vtkMoleculeMapper2.cxx
  vtkCIFMoleculeReader.cxx
  XSFReader.cxx
  AtomTypeDelegate.cxx
  BondLengthDelegate.cxx
  PovmolBondTableDialog.cxx
  PovmolMainWindow.cxx
  TableBasedBondDetector.cxx
  povmol.cxx
  )

QT4_WRAP_UI(UI_SOURCES
  PovmolMainWindowUi.ui
  PovmolBondTableDialogUi.ui
  )

QT4_WRAP_CPP(MOC_SOURCES
  PovmolMainWindow.h
  PovmolBondTableDialog.h
  AtomTypeDelegate.h
  BondLengthDelegate.h
  )

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  )

add_executable(povmol
  ${CXX_SOURCES}
  ${MOC_SOURCES}
  ${UI_SOURCES}
  )

target_link_libraries(povmol
  ${VTK_LIBRARIES}
  ${QT_LIBRARIES}
  )
