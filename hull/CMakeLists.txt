cmake_minimum_required(VERSION 2.8)
project(hull)

set(LIB_SOURCES
  hull.c
  ch.c
  io.c
  rand.c
  pointops.c
  fg.c
  )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(hull STATIC ${LIB_SOURCES})

if (BUILD_TESTING)
  add_executable(TestHull TestHull.cpp)
  target_link_libraries(TestHull hull)
endif()
