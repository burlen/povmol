cmake_minimum_required(VERSION 2.8)
PROJECT(2pov)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(Qt4 REQUIRED)
include(${QT_USE_FILE})

add_subdirectory(hull)

set(CXX_SOURCES
  AtomicProperties.cxx
  AtomicPropertiesBondDetector.cxx
  CoordinationPolyhedraGenerator.cxx
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
  PovmolLookupTableDialog.cxx
  QVTKPropertyWidget.cxx
  TableBasedBondDetector.cxx
  povmol.cxx
  )

QT4_WRAP_UI(UI_SOURCES
  PovmolBondTableDialogUi.ui
  PovmolMainWindowUi.ui
  PovmolLookupTableDialog.ui
  QVTKPropertyWidget.ui
  )

QT4_WRAP_CPP(MOC_SOURCES
  AtomTypeDelegate.h
  BondLengthDelegate.h
  PovmolBondTableDialog.h
  PovmolMainWindow.h
  PovmolLookupTableDialog.h
  QVTKPropertyWidget.h
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
  hull
  )

if (BUILD_TESTING)
  add_executable(TestHull
    CoordinationPolyhedraGenerator.cxx
    TestHull.cpp
    )
  target_link_libraries(TestHull
    ${VTK_LIBRARIES}
    hull
    )
endif()
