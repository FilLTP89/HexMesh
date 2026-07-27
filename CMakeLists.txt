cmake_minimum_required(VERSION 3.16)
project(HexMesh LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# ---------------------------------------------------------------------------
# MPI
# ---------------------------------------------------------------------------
find_package(MPI REQUIRED)

# ---------------------------------------------------------------------------
# HDF5 (mesh output). To use a custom build (e.g. a serial-only install kept
# separate from a parallel system package) point CMake at it with
# -DHDF5_ROOT=/path/to/hdf5 instead of editing this file.
# ---------------------------------------------------------------------------
find_package(HDF5 REQUIRED COMPONENTS CXX HL)

# ---------------------------------------------------------------------------
# glib (required by GTS)
# ---------------------------------------------------------------------------
find_package(PkgConfig REQUIRED)
pkg_check_modules(GLIB REQUIRED IMPORTED_TARGET glib-2.0)

# ---------------------------------------------------------------------------
# GTS (GNU Triangulated Surface library). It ships a pkg-config file but no
# CMake package config, so it is looked up under a single GTS_DIR prefix,
# editable from ccmake, defaulting to the common from-source install path.
# ---------------------------------------------------------------------------
set(GTS_DIR "/usr/local" CACHE PATH "GTS (GNU Triangulated Surface library) install prefix")
find_path(GTS_INCLUDE_DIR NAMES gts.h HINTS "${GTS_DIR}/include")
find_library(GTS_LIBRARY NAMES gts HINTS "${GTS_DIR}/lib")
if(NOT GTS_INCLUDE_DIR OR NOT GTS_LIBRARY)
  message(FATAL_ERROR
    "GTS not found under GTS_DIR='${GTS_DIR}'. Install it (see README.md) "
    "and point GTS_DIR at its install prefix, e.g. with ccmake.")
endif()
add_library(HexMesh::GTS UNKNOWN IMPORTED)
set_target_properties(HexMesh::GTS PROPERTIES
  IMPORTED_LOCATION "${GTS_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${GTS_INCLUDE_DIR}")
target_link_libraries(HexMesh::GTS INTERFACE PkgConfig::GLIB)

# ---------------------------------------------------------------------------
# libsc (octree/hash-table support library from the p4est project). Built
# with autotools and installs no CMake or pkg-config metadata, so, like GTS,
# it is a single cache variable pointing at the --prefix used to install it.
# ---------------------------------------------------------------------------
set(SC_DIR "" CACHE PATH "libsc install prefix (the --prefix used when building libsc)")
find_path(SC_INCLUDE_DIR NAMES sc.h HINTS "${SC_DIR}/include")
find_library(SC_LIBRARY NAMES sc HINTS "${SC_DIR}/lib")
if(NOT SC_INCLUDE_DIR OR NOT SC_LIBRARY)
  message(FATAL_ERROR
    "libsc not found. Set SC_DIR (currently '${SC_DIR}') to its install "
    "prefix, e.g. with ccmake.")
endif()
add_library(HexMesh::SC UNKNOWN IMPORTED)
set_target_properties(HexMesh::SC PROPERTIES
  IMPORTED_LOCATION "${SC_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${SC_INCLUDE_DIR}")

# ---------------------------------------------------------------------------
# CGAL + Boost, used by the topography/bathymetry surface-intersection code
# (src/intercept_surface.cpp). CGAL's own package config pulls in GMP/MPFR.
# ---------------------------------------------------------------------------
find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)

# ---------------------------------------------------------------------------
# Mesquite (optional). No active code path calls into it today -- mesh
# optimization is commented out in HexMesh.cpp, and its glue file
# (src/mesquite_interface_cpp) is intentionally kept out of the build. Turn
# this on and point MESQUITE_DIR at a checkout of
# git@github.com:FilLTP89/mesquite.git (built in-tree; its src/libmesquite.so
# is used directly, no `make install` needed) to revive that feature.
# ---------------------------------------------------------------------------
option(HEXMESH_ENABLE_MESQUITE "Link against Mesquite (unused by any active code path today)" OFF)
if(HEXMESH_ENABLE_MESQUITE)
  set(MESQUITE_DIR "" CACHE PATH "Mesquite source/build tree")
  find_library(MESQUITE_LIBRARY NAMES mesquite HINTS "${MESQUITE_DIR}/src")
  if(NOT MESQUITE_LIBRARY)
    message(FATAL_ERROR
      "HEXMESH_ENABLE_MESQUITE=ON but libmesquite was not found under "
      "MESQUITE_DIR='${MESQUITE_DIR}'.")
  endif()
  file(GLOB_RECURSE _mesquite_headers "${MESQUITE_DIR}/src/*.hpp" "${MESQUITE_DIR}/src/*.h")
  set(_mesquite_include_dirs "")
  foreach(_hdr ${_mesquite_headers})
    get_filename_component(_dir "${_hdr}" DIRECTORY)
    list(APPEND _mesquite_include_dirs "${_dir}")
  endforeach()
  list(REMOVE_DUPLICATES _mesquite_include_dirs)
  add_library(HexMesh::Mesquite UNKNOWN IMPORTED)
  set_target_properties(HexMesh::Mesquite PROPERTIES
    IMPORTED_LOCATION "${MESQUITE_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${_mesquite_include_dirs}")
endif()

# ---------------------------------------------------------------------------
# hexmesh executable
# ---------------------------------------------------------------------------
file(GLOB HEXMESH_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_executable(hexmesh ${HEXMESH_SOURCES} HexMesh.cpp)

target_include_directories(hexmesh PRIVATE include ${HDF5_INCLUDE_DIRS})
target_link_libraries(hexmesh PRIVATE
  MPI::MPI_CXX
  HexMesh::GTS
  HexMesh::SC
  CGAL::CGAL
  Boost::boost
  ${HDF5_LIBRARIES}
)
if(HEXMESH_ENABLE_MESQUITE)
  target_link_libraries(hexmesh PRIVATE HexMesh::Mesquite)
endif()
