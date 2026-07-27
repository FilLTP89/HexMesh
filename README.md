# HexMesh

## PRESENTATION

This project is focused on the development of a scalable hexahedral mesh generator for large domains based on octree and 27-tree structures. It allows for consideration of topography, bathymetry and coastlines, as well as water bodies and basins for geophysical applications.

The library was initiated at NACAD (Universidad Federal do Rio de Janeiro/COPPE, Brazil) and developed by LMA (Laboratoire de Mécanique et d'Acoustique, UMR 7031 AMU - CNRS - Centrale Marseille). This is fork of the original repository, mainly maintained and developed by LMPS (Laboratoire de Mécanique Paris-Saclay, UMR 9026 - Université Paris-Saclay, CentraleSupélec, CNRS, ENS Paris-Saclay).


* main developer: [Lúcio de Abreu Corrêa](mailto:labcorrea@gmail.com)
* contact (forked repository): [Filippo Gatti](mailto:filippo.gatti@centralesupelec.fr) (LMPS)
* initial developer: José Camata (NACAD), no longer maintaining the project
* contributors (by order of first commit): L. A. Corrêa (LMA), R. Cottereau (LMA), F. Gatti (LMPS, on forked repository)

HexMesh itself is written in C++/MPI. The Matlab routines that used to prepare GTS topography/bathymetry files have been retired from this repository; that preparation is now handled by the Python API in the companion `easyrisk-post` project.

## REFERENCES

If you use the library, please cite the following paper:
1. J. Camata, A. Coutinho. Parallel implementation and performance analysis of a linear octree finite element mesh generation scheme, _ Concurrency and Computation: Practice and Experience _ (2013), pp. 826-842. (http://dx.doi.org/10.1002/cpe.2869)

## DEPENDENCIES

HexMesh links against the following libraries. Where an Ubuntu/Debian package exists it is listed; the others (GTS, libsc, and optionally Mesquite) ship no distro package and are expected to be built from source, each into its own install prefix that you then point CMake at (see [COMPILATION](#compilation)).

| Library | Purpose | How to get it |
| --- | --- | --- |
| MPI | mesh generation is distributed across ranks | `sudo apt install libopenmpi-dev openmpi-bin` |
| HDF5 (C++ API) | mesh/output serialization (`H5Cpp.h`) | `sudo apt install libhdf5-openmpi-dev`, or a custom build (e.g. a serial-only install) pointed to via `HDF5_ROOT` |
| glib-2.0 | required by GTS | `sudo apt install libglib2.0-dev` |
| GTS (GNU Triangulated Surface library) | topography/bathymetry surface representation | build from source: <https://gts.sourceforge.net> (no Debian dev package ships a `.pc`/headers pair usable here) |
| libsc | octree/hash-table containers, from the p4est project | build from source: <https://github.com/cburstedde/libsc> |
| CGAL | exact-predicate surface/segment intersection (`src/intercept_surface.cpp`) | `sudo apt install libcgal-dev`, or build from source and expose it as an environment module (see below) |
| Boost | `boost::variant`/`static_visitor` used alongside CGAL | `sudo apt install libboost-dev` |
| Mesquite (optional, off by default) | mesh smoothing/optimization backend | not wired into any active code path today; build the fork at <git@github.com:FilLTP89/mesquite.git> only if you intend to revive `MeshOptimization` |

GTS and libsc both need GMP/MPFR-independent, plain C toolchains; libgmp-dev/libmpfr-dev are pulled in transitively by CGAL and don't need to be installed by hand.

### Building GTS

Download the source from <https://gts.sourceforge.net> (tarball release, or the repository linked from there), then:

```
>> cd gts
>> ./autogen.sh && ./configure --prefix=/usr/local
>> make && sudo make install
```

### Building libsc

```
>> git clone https://github.com/cburstedde/libsc.git
>> cd libsc
>> ./bootstrap
>> mkdir build && cd build
>> ../configure --prefix=$HOME/libsc/local --enable-mpi
>> make -j && make install
```

### Building CGAL as an environment module

CGAL is header-only (`make install` just copies headers and CMake package files, no compilation), so building it from source is cheap. This installs it into its own versioned prefix and exposes it via [environment-modules](https://modules.sourceforge.net/) instead of a system package, matching how `hdf5`/`fftw` are already handled under `~/modules` on this machine:

```
>> git clone --depth 1 --branch v6.2 https://github.com/CGAL/cgal.git
>> cmake -S cgal -B cgal/build -DCMAKE_INSTALL_PREFIX=$HOME/.local/cgal-6.2
>> cmake --install cgal/build
```

Then add a modulefile at `~/modules/cgal/6.2`:

```tcl
#%Module1.0#####################################################################

proc ModulesHelp { } {
    puts stderr "CGAL 6.2 (header-only) user install"
}

module-whatis "CGAL 6.2 header-only user installation"

set root $env(HOME)/.local/cgal-6.2

prepend-path CPATH             $root/include
prepend-path INCLUDE           $root/include
prepend-path CMAKE_PREFIX_PATH $root

# Picked up automatically by find_package(CGAL) via CMake's CMP0074 policy
setenv CGAL_ROOT $root
```

Then, before configuring HexMesh:

```
>> module load cgal/6.2
```

## COMPILATION

HexMesh is built with CMake (>= 3.16); there is no more hand-edited, per-machine Makefile. A normal out-of-source build looks like:

```
>> mkdir build && cd build
>> cmake .. -DSC_DIR=<libsc prefix> -DHDF5_ROOT=<HDF5 prefix>
>> make -j
```

`SC_DIR` has no default and must always be set (see below); `HDF5_ROOT` should be set explicitly too rather than relying on auto-detection: CMake's `FindHDF5` looks for an `h5cc`/`h5c++` compiler wrapper on `PATH` first, and if a `conda`/other user environment shadows the system one with a broken wrapper (e.g. one that shells out to a compiler binary that no longer exists), configuration fails with `Could NOT find HDF5 (missing: HDF5_INCLUDE_DIRS)` even though a perfectly good HDF5 is installed. Passing `HDF5_ROOT` makes it search that prefix's headers/libraries directly instead.

`cmake ..` will otherwise fail with a clear message if a dependency can't be found automatically (glib, CGAL, Boost, MPI are located with `find_package`/`pkg-config`; GTS and libsc are not, see below). Use `ccmake ..` (or `cmake-gui ..`) to interactively set the relevant cache variables instead of passing them all on the command line:

```
>> ccmake ..
```

Cache variables of interest:

* `GTS_DIR` (default `/usr/local`) — install prefix used when building GTS.
* `SC_DIR` (default empty, must be set) — `--prefix` used when building libsc.
* `HDF5_ROOT` — set this (as a normal CMake variable, e.g. `-DHDF5_ROOT=...`) to point at a non-system HDF5 build instead of the one `find_package(HDF5)` would pick up automatically.
* `HEXMESH_ENABLE_MESQUITE` (default `OFF`) and `MESQUITE_DIR` — turn Mesquite linkage on and point it at a Mesquite source/build tree; leave off otherwise.
* `CMAKE_BUILD_TYPE` (default `RelWithDebInfo`).

Example non-interactive configure with everything pinned explicitly:

```
>> cmake .. \
     -DGTS_DIR=/usr/local \
     -DSC_DIR=$HOME/libsc/local \
     -DHDF5_ROOT=/opt/hdf5-seq \
     -DCMAKE_BUILD_TYPE=Release
>> make -j
```

The resulting binary is `build/hexmesh`.

## USE

To create the mesh, run (from the build directory, or adjust the path to `hexmesh`):

```
>> mpirun -np <nb_proc> ./hexmesh <refinement_level> <gts_topofile> <outmesh_tag> <n_pml_layers> <pml_length> [gts_bathyfile]
```

* `<nb_proc>` — number of MPI processes used to create the mesh (each process writes its own output file).
* `<refinement_level>` — number of level refinements of the 27-tree structure (an integer corresponding to a power of 3: 4, 5, 6, ...).
* `<gts_topofile>` — path to the GTS file for the topography surface.
* `<outmesh_tag>` — output mesh file prefix.
* `<n_pml_layers>` — number of PML element layers extruded on each side of the domain.
* `<pml_length>` — physical thickness (in meters) of the PML region on each side of the domain.
* `[gts_bathyfile]` — optional path to the GTS file for the bathymetry surface. When omitted, the bathymetry-cutting steps (element interception, node projection, pillowing) are skipped entirely and the soil mesh (with its PML) is generated directly from the topography surface.

`gts_topofile`/`gts_bathyfile` are produced from STL files with the `stl2gts` tool (installed alongside GTS); the STL files themselves are generated by the Python API in `easyrisk-post`, which replaces the old `mainSRTM.m` Matlab pipeline.

By default (this will be made more general later), the depth of the mesh is the larger dimension of the two horizontal dimensions of the topography file; and the depths at which the refinements occur are set in function `hexa_tree_cube`, in `hexa.cpp`.
