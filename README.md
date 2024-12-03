# HexMesh

## PRESENTATION

This project is focused on the development of a scalable hexahedral mesh generator for large domains based on octree and 27-tree structures. It allows for consideration of topography, bathymetry and coastlines, as well as water bodies and basins for geophysical applications.

The library was initiated at NACAD (Universidad Federal do Rio de Janeiro/COPPE, Brazil) and developed by LMA (Laboratoire de Mécanique et d'Acoustique, UMR 7031 AMU - CNRS - Centrale Marseille). This is fork of the original repository, mainly maintained and developed by LMPS (Laboratoire de Mécanique Paris-Saclay, UMR 9026 - Université Paris-Saclay, CentraleSupélec, CNRS, ENS Paris-Saclay).


* contact : [José Camata](mailto:camata@nacad.ufrj.br) [Lúcio de Abreu Corrêa](mailto:labcorrea@gmail.com)
* contributors (by order of first commit): L. A. Corrêa (LMA), R. Cottereau (LMA), F. Gatti (LMPS, on forked repository)

HexMesh is written in C++/MPI. Additional routines are written in Matlab for the preparation of GTS topography files.
 
## REFERENCES

If you use the library, please cite the following paper:
1. J. Camata, A. Coutinho. Parallel implementation and performance analysis of a linear octree finite element mesh generation scheme, _ Concurrency and Computation: Practice and Experience _ (2013), pp. 826-842. (http://dx.doi.org/10.1002/cpe.2869)

## INSTALLATION

Before using the software, the following libraries need to be installed and available:
1. gts (https://gts.sourceforge.net)
2. libsc (https://github.com/cburstedde/libsc)
3. hdf5 (https://www.hdfgroup.org/download-hdf5)
4. mesquite (https://www.mesquiteproject.org/Installation.html)

## COMPILATION

Depending on the OS you are using, modify the paths for GTS_LIB, SC_LIB, HDF5_DIR, MESQUITE_DIR and GLIB_INCLUDE in Makefile

Compile with (replace OS by Linux or mac):

```
>> make -f Makefile
```

## USE

To prepare the geometry files, modify the headers in `mainSRTM.m` (in particular choose the bounding box in latitude/longitude) and run in Matlab:

```
>> mainSRTM
```

You need an internet connection to download the surface topography, bathymetric coastlines files (no connection needed if they are already available on your computer). The output files are a topography STL file <stl_topofile> and a bathymetry STL file  <stl_bathyfile>. These files should be transformed to GTS using `stl2gts` command.

To create the mesh, you should run (in a Terminal from the directory $(HEXMESH)):

```
>> mpirun -np <nb_proc> ./hexmesh <refinement_level> <gts_topofile> <gts_bathyfile> <tag_mesh>
```

where <nb_proc> is an integer specifying the number of processes used to create the mesh (each process creates its own VTK file), and <refinement_level> is an integer specifying the number of level refinements of the 27-tree structure! For now, <refinement_level> is an integer (corresponding to a power of 3): 4, 5, 6, ...<gts_topofile> and <gts_bathyfile> are the paths to the GTS files corresponding to topography and bathymetry. <tag_mesh> corresponds to the output mesh file prefix.

To create a mesh without bathymetry cut, you should run (in a Terminal from the directory $(HEXMESH))

```
>> mpirun -np <nb_proc> ./hexmesh <refinement_level> <gts_topofile> <gts_topofile> <tag_mesh>
```

By default (this will be made more general later), the depth of the mesh is the larger dimension of the two horizontal dimensions of the topography file; and the depths at which the refinements occur at set in function `hexa_tree_cube`, in `hexa.cpp` (line 190).
