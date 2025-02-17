/*
 * File:   HexMesh.cpp
 * Author: Lucio de Abreu Correa
 *
 * Created on April 28, 2021, 12.00 PM
 */

#include <cstdlib>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>

#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>
#include "hexa.h"
#include "hilbert.h"
// #include "refinement.h"
#include <ctime>
#include <chrono>

// void GetMeshFromSurface(hexa_tree_t *tree, const char *surface_topo, std::vector<double> &coords);
// // void GetInterceptedElements(hexa_tree_t *mesh, std::vector<double> &coords, std::vector<int> &elements_ids, const char *surface_bathy);
// void CheckOctreeTemplate(hexa_tree_t *mesh, const std::vector<double> &coords, std::vector<int> &elements_ids, bool flag);
// void ApplyOctreeTemplate(hexa_tree_t *mesh, std::vector<double> &coords, std::vector<int> &elements_ids);

// int CheckTemplate(hexa_tree_t *mesh, const std::vector<double> &coords, std::vector<int> &elements_ids, bool flag);
// void CutTemplate(hexa_tree_t *mesh, std::vector<double> &coords, std::vector<int> &elements_ids);

// void Apply_material(hexa_tree_t *mesh, std::vector<double> &coords, std::vector<int> &element_ids, const char *surface_bathy);
// void Adjust_material(hexa_tree_t *mesh);
// // void AddPMLElements(hexa_tree_t* mesh);
// void ExtrudePMLElements(hexa_tree_t *mesh, std::vector<double> &coords);

// void IdentifyTemplate(hexa_tree_t *mesh, const std::vector<double> &coords, std::vector<int> &elements_ids);
// void MovingNodes(hexa_tree_t *mesh, std::vector<double> &coords, std::vector<int> &nodes_b_mat);
// void MeshOpt(hexa_tree_t *mesh, std::vector<double> &coords, std::vector<int> material_fixed_nodes);
// void UntagleMesh(hexa_tree_t *mesh, std::vector<double> &coords, std::vector<int> material_fixed_nodes);

int main(int argc, char *argv[])
{

  hexa_tree_t mesh;

  std::vector<double> coords;
  std::vector<int> element_ids;
  std::vector<int> nodes_b_mat;
  auto start = std::chrono::steady_clock::now();
  int l = atoi(argv[1]);

  // mpi init
  hexa_init(argc, argv, &mesh);
  // set the initial number of elements in x,y,z
  hexa_tree_init(&mesh, l);
  // build the referene mesh
  hexa_tree_cube(&mesh);

  // deal with the mpi com
  hexa_mesh(&mesh);

  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the initialization %lld millisecond(s).\n", elapsed.count());
  std::cout << "Time in the initialization " << elapsed.count() << " millisecond(s)." << std::endl;

  // const char *bathy = argv[2];
  const char *topo = argv[2];
  const char *outmesh = argv[3];
  printf("Loading files: \n \t %s \n", topo);

  start = std::chrono::steady_clock::now();
  // Note that here we use a gts file.
  // There is a tool called stl2gts that convert STL files to GTS.
  // It is installed together with the gts library.
  // create the geometrical mesh
  GetMeshFromSurface(&mesh, topo, coords);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the GetMeshFromSurface %lld millisecond(s).\n", elapsed.count());
  std::cout << "Time in the GetMeshFromSurface " << elapsed.count() << " millisecond(s)." << std::endl;

  // find the elements intercepted by the bathy
  // start = std::chrono::steady_clock::now();
  //GetInterceptedElements(&mesh, coords, element_ids, bathy);
  //printf(" Elements intercepted: %lld\n\n", element_ids.size());
  // elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  // fprintf(mesh.profile, "Time in the GetInterceptedElements %lld millisecond(s).\n", elapsed.count());
  // std::cout << "Time in GetInterceptedElements " << elapsed.count() << " millisecondsecond(s)." << std::endl;

  // apply a deformation in the mesh to fit the bathy
  start = std::chrono::steady_clock::now();
  printf(" Project nodes to the surface\n\n");
  MovingNodes(&mesh, coords, nodes_b_mat);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the MovingNodes %lld millisecond(s).\n", elapsed.count());
  std::cout << "Time in MovingNodes " << elapsed.count() << " millisecond(s)." << std::endl;

  // apply material
  start = std::chrono::steady_clock::now();
  printf(" Applying material \n\n");
  element_ids.clear();
  // Apply_material(&mesh, coords, bathy);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the Apply_material %lld millisecond(s).\n", elapsed.count());
  std::cout << "Time in Apply_material " << elapsed.count() << " millisecond(s)." << std::endl;
  
  // do the pillow
  start = std::chrono::steady_clock::now();
  printf(" Applying pillowing process\n\n");
  PillowingInterface(&mesh, coords, nodes_b_mat);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the PillowingInterface %lld millisecond(s).\n", elapsed.count());
  std::cout << "Time in PillowingInterface " << elapsed.count() << " millisecond(s)." << std::endl;

  element_ids.clear();

  start = std::chrono::steady_clock::now();
  printf(" Extrude elements\n\n");
  ExtrudePMLElements(&mesh, coords);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the ExtrudePMLElements %lld millisecond(s).\n", elapsed.count());
  std::cout << "Time in ExtrudePMLElements " << elapsed.count() << " millisecond(s)." << std::endl;

  // clean vectors
  std::vector<int>().swap(element_ids);
  std::vector<int>().swap(nodes_b_mat);

  printf(" Writing output files \n\n");
  // hexa_mesh_write_vtk(&mesh, "mesh", &coords);
  // hexa_mesh_write_msh(&mesh, "mesh", &coords);
  hexa_mesh_write_h5(&mesh, outmesh, coords);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in Writing output files %lld millisecond(s).\n", elapsed.count());
  std::cout << "Time in Writing output files " << elapsed.count() << " millisecond(s)." << std::endl;

  start = std::chrono::steady_clock::now();
  printf(" Cleaning variables \n\n");

  hexa_tree_destroy_nocut(&mesh);
  hexa_finalize(&mesh);
  std::vector<double>().swap(coords);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the Cleaning variables %lld millisecond(s).\n", elapsed.count());
  std::cout << "Time in Cleaning variables " << elapsed.count() << " millisecond(s)." << std::endl;

  return 0;
}
