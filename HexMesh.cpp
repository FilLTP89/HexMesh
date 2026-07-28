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
#include <ctime>
#include <chrono>

namespace
{

// Settings read from the HexMesh.input-style config file. Sentinel negative
// values mark a numeric setting as "not yet seen" so missing required keys
// can be reported together after the whole file has been read.
struct HexMeshConfig
{
  int ref = -1;
  std::string topo;
  std::string outmesh;
  int n_pml_layers = -1;
  double pml_length = -1.0;
  std::string bathy; // empty => no bathymetry cut
};

std::string Trim(const std::string &s)
{
  size_t b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos)
    return "";
  size_t e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

HexMeshConfig ReadHexMeshConfig(const std::string &path)
{
  std::ifstream file(path);
  if (!file.is_open())
  {
    std::cerr << "Error: could not open config file '" << path << "'" << std::endl;
    exit(EXIT_FAILURE);
  }

  HexMeshConfig cfg;
  std::string line;
  int line_no = 0;
  while (std::getline(file, line))
  {
    ++line_no;
    size_t comment = line.find('#');
    if (comment != std::string::npos)
      line = line.substr(0, comment);
    line = Trim(line);
    if (line.empty())
      continue;

    size_t eq = line.find('=');
    if (eq == std::string::npos)
    {
      std::cerr << "Warning: " << path << ":" << line_no
                 << ": ignoring malformed line '" << line << "' (expected 'key = value')" << std::endl;
      continue;
    }
    std::string key = Trim(line.substr(0, eq));
    std::string value = Trim(line.substr(eq + 1));

    try
    {
      if (key == "ref")
        cfg.ref = std::stoi(value);
      else if (key == "topo")
        cfg.topo = value;
      else if (key == "outmesh")
        cfg.outmesh = value;
      else if (key == "n_pml_layers")
        cfg.n_pml_layers = std::stoi(value);
      else if (key == "pml_length")
        cfg.pml_length = std::stod(value);
      else if (key == "bathy")
        cfg.bathy = value;
      else
        std::cerr << "Warning: " << path << ":" << line_no << ": unknown key '" << key << "'" << std::endl;
    }
    catch (const std::exception &)
    {
      std::cerr << "Error: " << path << ":" << line_no
                 << ": invalid value '" << value << "' for key '" << key << "'" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  std::vector<std::string> missing;
  if (cfg.ref < 0)
    missing.push_back("ref");
  if (cfg.topo.empty())
    missing.push_back("topo");
  if (cfg.outmesh.empty())
    missing.push_back("outmesh");
  if (cfg.n_pml_layers < 0)
    missing.push_back("n_pml_layers");
  if (cfg.pml_length < 0)
    missing.push_back("pml_length");
  if (!missing.empty())
  {
    std::cerr << "Error: " << path << ": missing required setting(s):";
    for (const auto &m : missing)
      std::cerr << " " << m;
    std::cerr << std::endl;
    exit(EXIT_FAILURE);
  }

  return cfg;
}

} // namespace

int main(int argc, char *argv[])
{

  hexa_tree_t mesh;

  std::vector<double> coords;
  std::vector<int> element_ids;
  std::vector<int> nodes_b_mat;
  auto start = std::chrono::steady_clock::now();

  std::string config_path = (argc > 1) ? argv[1] : "HexMesh.input";
  HexMeshConfig cfg = ReadHexMeshConfig(config_path);

  // mpi init
  hexa_init(argc, argv, &mesh);
  // set the initial number of elements in x,y,z
  hexa_tree_init(&mesh, cfg.ref);
  // build the referene mesh
  hexa_tree_cube(&mesh);

  // deal with the mpi com
  hexa_mesh(&mesh);

  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the initialization %lld millisecond(s).\n", elapsed.count());
  std::cout << "Time in the initialization " << elapsed.count() << " millisecond(s)." << std::endl;

  const char *topo = cfg.topo.c_str();
  const char *outmesh = cfg.outmesh.c_str();
  int n_pml_layers = cfg.n_pml_layers;
  double pml_length = cfg.pml_length;
  const char *bathy = cfg.bathy.empty() ? nullptr : cfg.bathy.c_str();
  if (bathy)
  {
    printf("Loading files:\n \t %s \n \t %s \n", bathy, topo);
  }
  else
  {
    printf("Loading files:\n \t %s \n", topo);
  }

  start = std::chrono::steady_clock::now();
  // Note that here we use a gts file.
  // There is a tool called stl2gts that convert STL files to GTS.
  // It is installed together with the gts library.
  // create the geometrical mesh
  GetMeshFromSurface(&mesh, topo, coords);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the GetMeshFromSurface %lld millisecond(s).\n", elapsed.count());
  std::cout << "Time in the GetMeshFromSurface " << elapsed.count() << " millisecond(s)." << std::endl;

  if (bathy)
  {
    // find the elements intercepted by the bathy
    start = std::chrono::steady_clock::now();
    GetInterceptedElements(&mesh, coords, element_ids, bathy);
    printf(" Elements intercepted: %lld\n\n", element_ids.size());
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
    fprintf(mesh.profile, "Time in the GetInterceptedElements %lld millisecond(s).\n", elapsed.count());
    std::cout << "Time in GetInterceptedElements " << elapsed.count() << " millisecondsecond(s)." << std::endl;

    // apply a deformation in the mesh to fit the bathy
    start = std::chrono::steady_clock::now();
    printf(" Project nodes to the surface\n\n");
    MovingNodes(&mesh, coords, nodes_b_mat);
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
    fprintf(mesh.profile, "Time in the MovingNodes %lld millisecond(s).\n", elapsed.count());
    std::cout << "Time in MovingNodes " << elapsed.count() << " millisecond(s)." << std::endl;
  }

  // apply material
  start = std::chrono::steady_clock::now();
  printf(" Applying material \n\n");
  element_ids.clear();
  Apply_material(&mesh, coords);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the Apply_material %lld millisecond(s).\n", elapsed.count());
  std::cout << "Time in Apply_material " << elapsed.count() << " millisecond(s)." << std::endl;

  if (bathy)
  {
    // pillowing
    start = std::chrono::steady_clock::now();
    printf(" Applying pillowing process\n\n");
    PillowingInterface(&mesh, coords, nodes_b_mat);
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
    fprintf(mesh.profile, "Time in the PillowingInterface %lld millisecond(s).\n", elapsed.count());
    std::cout << "Time in PillowingInterface " << elapsed.count() << " millisecond(s)." << std::endl;
  }

  // // opt mesh
  // start = std::chrono::steady_clock::now();
  // printf(" Mesh Optimization\n\n");
  // // MeshOptimization(&mesh, coords, nodes_b_mat);
  // elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  // fprintf(mesh.profile, "Time in the MeshOptimization %lld millisecond(s).\n", elapsed.count());
  // std::cout << "Time in MeshOptimization " << elapsed.count() << " millisecond(s)." << std::endl;

  start = std::chrono::steady_clock::now();
  printf(" Extrude elements\n\n");
  ExtrudePMLElements(&mesh, coords, n_pml_layers, pml_length);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the ExtrudePMLElements %lld millisecond(s).\n", elapsed.count());
  std::cout << "Time in ExtrudePMLElements " << elapsed.count() << " millisecond(s)." << std::endl;

  printf(" Writing output files \n\n");
  hexa_mesh_write_h5(&mesh, outmesh, coords);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in Writing output files %lld millisecond(s).\n", elapsed.count());
  std::cout << "Time in Writing output files " << elapsed.count() << " millisecond(s)." << std::endl;

  start = std::chrono::steady_clock::now();
  printf(" Cleaning variables \n\n");

  hexa_tree_destroy(&mesh);
  std::vector<double>().swap(coords);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the Cleaning variables %lld millisecond(s).\n", elapsed.count());
  std::cout << "Time in Cleaning variables " << elapsed.count() << " millisecond(s)." << std::endl;

  // hexa_finalize() closes mesh.profile and calls MPI_Finalize(), so it must
  // run last, after every other use of mesh.profile/MPI above.
  hexa_finalize(&mesh);

  return 0;
}
