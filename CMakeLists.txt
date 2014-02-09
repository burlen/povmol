cmake_minimum_required(VERSION 2.8)
PROJECT(povmol)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set(CXX_SOURCES
  vtkMoleculeMapper2.cxx
  POVRayPrimatives.cxx
  povmol.cxx
  )

add_executable(povmol ${CXX_SOURCES})
target_link_libraries(povmol ${VTK_LIBRARIES})
