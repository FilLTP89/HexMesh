/* 
 * File:   main.cpp
 * Author: camata
 *
 * Created on March 19, 2015, 1:42 PM
 */

#include <cstdlib>
#include <vector>

#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>
#include "hexa.h"
#include "hilbert.h"

//void GetMeshFromSurface(hexa_tree_t* tree, const char* surface_bathy, const char* surface_topo, const char* coastline, std::vector<double>& coords);
void GetMeshFromSurface(hexa_tree_t* tree, const char* surface_topo, std::vector<double>& coords);
//void GetInterceptedElements(hexa_tree_t* tree, std::vector<double>& coords, sc_array_t* intercepted_elements);
//void ApllyTemplate(hexa_tree_t* mesh, std::vector<double>& coords, sc_array_t* intercepted_elements);

void GetInterceptedElements(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, const char* surface_bathy);
void ApllyTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids);
void CheckTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids);


void AddPMLElements(hexa_tree_t* mesh);

/*
 * 
 */

int main(int argc, char** argv) {

    hexa_tree_t mesh;

    std::vector<double> coords;
    std::vector<int> element_ids;
    int l = atoi(argv[1]);

    hexa_init(argc, argv, &mesh);

    hexa_tree_init(&mesh, l);
    hexa_tree_cube(&mesh);


    //hexa_debug_face_hanging(&mesh);
    
    hexa_mesh(&mesh);
    GetMeshFromSurface(&mesh, "./input/Topo.gts", coords);
    //GetMeshFromSurface(&mesh, "./input/Bathymetry.gts", "./input/Topo.gts", "./input/Argostoli_coastline.dat", coords);
    //GetMeshFromSurface(&mesh, "./input/01.gts", "./input/Topo.gts", "./input/Argostoli_coastline.dat", coords);


    GetInterceptedElements(&mesh, coords, element_ids, "./input/Bathymetry.gts");
    //GetInterceptedElements(&mesh, coords, element_ids, "./input/template3_1.gts");
    //GetInterceptedElements(&mesh, coords, element_ids, "./input/template1_z.gts");


    //printf(" Elements intercepted: %d\n", element_ids.size());

    CheckTemplate(&mesh, coords, element_ids);
    //ApllyTemplate(&mesh, coords, element_ids);
    //TO DO
    //AddPMLElements(&mesh);
    
    hexa_mesh_write_vtk(&mesh, "mesh", &coords);

    //hexa_mesh_write_unv(&mesh, "mesh", &coords);
    //hexa_mesh_write_vtk(&mesh,"template", NULL);
    //hexa_mesh_write_unv(&mesh,"teste", NULL);


    hexa_tree_destroy(&mesh);
    hexa_finalize(&mesh);
    return 0;
}

void add_pml_elements(hexa_tree_t *mesh) {

}
