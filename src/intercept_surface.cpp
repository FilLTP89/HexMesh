
#include <gts.h>
#include <glib.h>
#include <vector>
#include <iostream>
using namespace std;
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>

#include "hexa.h"
#include "hilbert.h"
#include "coastline.h"

void Coastline_Reader(const char* file, SetOfCoastLine* scl);

static unsigned nDims = 2;
static unsigned nBits = 32;

typedef struct {
    bitmask_t coord[2];
    int node_id;
} node_in_edge_t;

/*
unsigned edge_hash_fn(const void *v, const void *u) {
    const node_in_edge_t *q = (const node_in_edge_t*) v;
    bitmask_t index;
    index = hilbert_c2i(2, 32, q->coord);
    return (unsigned) index;
}
 * 
 */

unsigned edge_hash_fn(const void *v, const void *u) {
    const node_in_edge_t *q = (const node_in_edge_t*) v;
    uint32_t a, b, c;

    a = (uint32_t) q->coord[0];
    b = (uint32_t) q->coord[1];
    c = (uint32_t) 0;
    sc_hash_mix(a, b, c);
    sc_hash_final(a, b, c);
    return (unsigned) c;
}

int edge_equal_fn(const void *v, const void *u, const void *w) {
    const node_in_edge_t *e1 = (const node_in_edge_t*) v;
    const node_in_edge_t *e2 = (const node_in_edge_t*) u;

    return (unsigned) ((e1->coord[0] == e2->coord[0]) && (e2->coord[1] == e2->coord[1]));

}

void shift_n_position_from_index(sc_array_t *elements, uint32_t n, uint32_t index);

int EdgeVerticesMap[12][2] = {
    {0, 1}, // Edge 0
    {1, 2}, // Edge 1
    {2, 3}, //      2
    {3, 0},
    {0, 4},
    {1, 5},
    {2, 6},
    {3, 7},
    {4, 5},
    {5, 6},
    {6, 7},
    {7, 4}
};

int EdgeVerticesMap_surf_diagonal[12][2] = {
    {0, 5},
    {1, 4},
    {3, 6},
    {2, 7},
    {0, 7},
    {3, 4},
    {1, 6},
    {2, 5},
    {4, 6},
    {5, 7},
    {0, 2},
    {1, 3}
};

int EdgeVerticesMap_vol_diagonal[4][2] = {
    {0, 6},
    {1, 7},
    {2, 4},
    {3, 5}
};

int FaceEdgesMap[6][4] = {
    {4, 11, 7, 3},
    {5, 1, 6, 9},
    {0, 5, 8, 4},
    {2, 6, 10, 7},
    {8, 9, 10, 11},
    {0, 1, 2, 3}
};

int EdgeEdgeMap[12][4] = {
    {3, 4, 1, 5}, // Edge 0
    {0, 5, 2, 6}, // Edge 1
    {1, 6, 3, 7}, //      2
    {2, 7, 0, 4,}, //3  
    {0, 3, 8, 11}, //4 
    {1, 0, 9, 8}, //5
    {2, 1, 10, 9}, //6
    {3, 2, 11, 10}, //7
    {11, 4, 9, 5}, //8
    {8, 5, 10, 6}, //9 
    {9, 6, 11, 7}, //10 
    {10, 7, 8, 9}//11
};

// Read the gts file format and create a gts surface.

GtsSurface* SurfaceRead(const char* fname) {
    FILE *gts_file;
    GtsSurface *s;
    GtsPoint *p;
    GtsFile *fp;

    gts_file = fopen(fname, "r");
    if (!gts_file) return NULL;
    fp = gts_file_new(gts_file);
    s = gts_surface_new(gts_surface_class(),
            gts_face_class(),
            gts_edge_class(),
            gts_vertex_class());
    if (gts_surface_read(s, fp)) {
        fputs("file on standard input is not a valid GTS file\n", stderr);
        fprintf(stderr, "stdin:%d:%d: %s\n", fp->line, fp->pos, fp->error);
        return NULL; /* failure */
    }

    gts_file_destroy(fp);
    return s;

}

// Compute the distance between a point and a triangle.

gdouble distance(GtsPoint *p, gpointer bounded) {
    GtsTriangle *t = (GtsTriangle*) bounded;
    return gts_point_triangle_distance(p, t);
}

int pnpoly(int nvert, std::vector<Point2D>& vert, double testx, double testy) {
    int i, j, c = 0;
    for (i = 0, j = nvert - 1; i < nvert; j = i++) {
        if (((vert[i].y > testy) != (vert[j].y > testy)) &&
                (testx < (vert[j].x - vert[i].x) * (testy - vert[i].y) / (vert[j].y - vert[i].y) + vert[i].x))
            c = !c;
    }
    return c;
}

int AddPointOnEdge(int* nodes, sc_hash_array_t* hash, int &npoints, GtsPoint *p, std::vector<double> &coords);

//void GetMeshFromSurface(hexa_tree_t* tree, const char* surface_bathy, const char* surface_topo, const char* coastline, vector<double>& coords) {

void GetMeshFromSurface(hexa_tree_t* tree, const char* surface_topo, vector<double>& coords) {

    GtsSurface *s;
    GtsBBox *box;
    GNode *t;
    GtsPoint *p;
    int nx, ny, nz;
    double dx, dy, dz;
    sc_array_t *nodes = &tree->nodes;
    sc_array_t *elements = &tree->elements;

    //hexa_tree_t tree_t;


    // Note that here we use a gts file.
    // There is a tool called stl2gts that convert STL files to GTS.
    // It is installed together with the gts library.
    tree->gdata.s = SurfaceRead(surface_topo);
    //tree->gdata.s = SurfaceRead(surface_bathy);
    //tree_t.gdata.s = SurfaceRead(surface_topo);

    FILE *fout = fopen("surface.dat", "w");
    gts_surface_print_stats(tree->gdata.s, fout);
    fclose(fout);

    // Get the surface bounding box
    tree->gdata.bbox = gts_bbox_surface(gts_bbox_class(), tree->gdata.s);
    //tree_t.gdata.bbox = gts_bbox_surface(gts_bbox_class(), tree_t.gdata.s);
    if (tree->mpi_rank == 0) {
        printf("Bounding box: \n");
        printf(" x ranges from %f to %f\n", tree->gdata.bbox->x1, tree->gdata.bbox->x2);
        printf(" y ranges from %f to %f\n", tree->gdata.bbox->y1, tree->gdata.bbox->y2);
        printf(" z ranges from %f to %f\n", tree->gdata.bbox->z1, tree->gdata.bbox->z2);
    }

    //double factor = 0.01;
    double x_factor = (tree->gdata.bbox->x2 - tree->gdata.bbox->x1)*0.01;
    double y_factor = (tree->gdata.bbox->y2 - tree->gdata.bbox->y1)*0.01;

    tree->gdata.bbox->x1 += x_factor;
    tree->gdata.bbox->y1 += y_factor;

    tree->gdata.bbox->x2 -= x_factor;
    tree->gdata.bbox->y2 -= y_factor;

    double Lx = (tree->gdata.bbox->x2 - tree->gdata.bbox->x1);
    double Ly = (tree->gdata.bbox->y2 - tree->gdata.bbox->y1);
    double zmin = ((Lx < Ly) ? -Lx : -Ly);
    //double zmin = (tree->gdata.bbox->x1<tree->gdata.bbox->y1)?tree->gdata.bbox->x1: tree->gdata.bbox->y1;

    // Get grid-spacing at x and y direction
    dx = (tree->gdata.bbox->x2 - tree->gdata.bbox->x1) / (double) tree->ncellx;
    dy = (tree->gdata.bbox->y2 - tree->gdata.bbox->y1) / (double) tree->ncelly;

    coords.resize(nodes->elem_count * 3);

    // Build the bounding box tree
    tree->gdata.bbt = gts_bb_tree_surface(tree->gdata.s);
    //tree_t.gdata.bbt = gts_bb_tree_surface(tree_t.gdata.s);

    p = gts_point_new(gts_point_class(), 0.0, 0.0, tree->gdata.bbox->z2);

    //GSList *list = gts_bb_tree_overlap(t, box);
    // Call the coastline reader
    //SetOfCoastLine scl;
    //Coastline_Reader(coastline, &scl);

    for (int i = 0; i < nodes->elem_count; ++i) {
        octant_node_t* n = (octant_node_t*) sc_array_index(nodes, i);
        p->x = tree->gdata.bbox->x1 + n->x*dx;
        p->y = tree->gdata.bbox->y1 + n->y*dy;
        //double xp = tree->gdata.bbox->x1 + n->x*dx;
        //double yp = tree->gdata.bbox->y1 + n->y*dy;
        double d;
        double zmax;
        //for (int ii = 0; ii < scl.npart - 1; ++ii) {
        //    int tes = pnpoly(scl.coastlines[ii].npoints, scl.coastlines[ii].points, xp, yp);
        //    if (tes == 0) {
        d = gts_bb_tree_point_distance(tree->gdata.bbt, p, distance, NULL);
        zmax = tree->gdata.bbox->z2 - d;
        //        zmax = 0;
        //    } else {
        //        d = gts_bb_tree_point_distance(tree_t.gdata.bbt, p, distance, NULL);
        //        zmax = tree_t.gdata.bbox->z2 - d;
        //        break;
        //    }
        //}

        dz = (zmax - zmin) / (double) tree->ncellz;
        double z = zmax - (n->z) * dz;

        coords[i * 3 + 0] = p->x;
        coords[i * 3 + 1] = p->y;
        coords[i * 3 + 2] = z;
    }
}

void GetInterceptedElements(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, const char* surface_bathy) {
    sc_array_t *elements = &mesh->elements;
    GtsBBox *box;

    //hexa_tree_t bathy;

    mesh->gdata.s = SurfaceRead(surface_bathy);
    mesh->gdata.bbox = gts_bbox_surface(gts_bbox_class(), mesh->gdata.s);
    mesh->gdata.bbt = gts_bb_tree_surface(mesh->gdata.s);

    box = gts_bbox_new(gts_bbox_class(), 0, 0, 0, 0, 1, 1, 1);

    for (int iel = 0; iel < elements->elem_count; ++iel) {

        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

        box->x1 = box->y1 = box->z1 = 1.0E10;
        box->x2 = box->y2 = box->z2 = -1.0E10;

        for (int i = 0; i < 8; ++i) {
            octant_node_t* node = &elem->nodes[i];
            int id = node->id;
            double x = coords[id * 3];
            double y = coords[id * 3 + 1];
            double z = coords[id * 3 + 2];
            box->x1 = (x < box->x1) ? x : box->x1;
            box->y1 = (y < box->y1) ? y : box->y1;
            box->z1 = (z < box->z1) ? z : box->z1;
            box->x2 = (x > box->x2) ? x : box->x2;
            box->y2 = (y > box->y2) ? y : box->y2;
            box->z2 = (z > box->z2) ? z : box->z2;
        }
        elem->pad = 0;

        if (gts_bb_tree_is_overlapping(mesh->gdata.bbt, box)) {
            elements_ids.push_back(iel);
            elem->pad = 1;
        }
    }
}

/*

 1. Para cada octante, 

 1.a vefificar se duas faces paralelas nao sao interceptadas pela
 superficie.
   => Aplica-se o template 1: divide elemento ao meio gerando dois novos elementos

 1.b: quatro faces inteceptadas e as restantes nao sao paralelas
   => Aplica-se o template 2:
           - Adicionar ponto na intersecao seguindo a projecao diagonal do ponto da esquina do elemento, 
           -  Tres novos elementos são criados. 
            Obs.: Possibilidade de formacao de triangulos.  Solucao: mover pontos.

 1.c: tres faces interceptadas: Veja template3 hexMesh/doc.


 */


GtsPoint* SegmentTriangleIntersection(GtsSegment * s, GtsTriangle * t) {
    GtsPoint * A, * B, * C, * D, * E;
    gint ABCE, ABCD, ADCE, ABDE, BCDE;
    GtsEdge * AB, * BC, * CA;
    gdouble a, b, c;

    g_return_val_if_fail(s != NULL, NULL);
    g_return_val_if_fail(t != NULL, NULL);
    //if( s == NULL) return NULL;
    //if (t == NULL) return NULL;

    gts_triangle_vertices_edges(t, NULL,
            (GtsVertex **) & A,
            (GtsVertex **) & B,
            (GtsVertex **) & C,
            &AB, &BC, &CA);
    D = GTS_POINT(s->v1);
    E = GTS_POINT(s->v2);

    ABCE = gts_point_orientation_3d_sos(A, B, C, E);
    ABCD = gts_point_orientation_3d_sos(A, B, C, D);
    if (ABCE < 0 || ABCD > 0) {
        GtsPoint * tmpp;
        gint tmp;

        tmpp = E;
        E = D;
        D = tmpp;
        tmp = ABCE;
        ABCE = ABCD;
        ABCD = tmp;
    }
    if (ABCE < 0 || ABCD > 0)
        return NULL;
    ADCE = gts_point_orientation_3d_sos(A, D, C, E);
    if (ADCE < 0)
        return NULL;
    ABDE = gts_point_orientation_3d_sos(A, B, D, E);
    if (ABDE < 0)
        return NULL;
    BCDE = gts_point_orientation_3d_sos(B, C, D, E);
    if (BCDE < 0)
        return NULL;
    a = gts_point_orientation_3d(A, B, C, E);
    b = gts_point_orientation_3d(A, B, C, D);
    if (a != b) {
        c = a / (a - b);
        return gts_point_new(gts_point_class(),
                E->x + c * (D->x - E->x),
                E->y + c * (D->y - E->y),
                E->z + c * (D->z - E->z));
    }
    /* D and E are contained within ABC */
#ifdef DEBUG
    fprintf(stderr,
            "segment: %p:%s triangle: %p:%s intersection\n"
            "D and E contained in ABC\n",
            s, GTS_NEDGE(s)->name, t, GTS_NFACE(t)->name);
#endif /* DEBUG */  
    g_assert(a == 0.0);
    return gts_point_new(gts_point_class(),
            (E->x + D->x) / 2.,
            (E->y + D->y) / 2.,
            (E->z + D->z) / 2.);
}

//void ApllyTemplate(hexa_tree_t* mesh, std::vector<double>& coords, sc_array_t* intercepted_elements) {

void ApllyTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids) {
    bool clamped = true;
    GtsSegment * segments[12];
    bool face_intecepted[6];
    GtsPoint * point[12];
    int Edge2GNode[12][2];
    int conn_p[4];
    int original_conn[8];
    FILE * fdbg;

    int id_offset = 0;

    fdbg = fopen("intercepted_faces.dbg", "w");

    sc_hash_array_t* hash_nodes = sc_hash_array_new(sizeof (node_in_edge_t), edge_hash_fn, edge_equal_fn, &clamped);

    for (int iel = 0; iel < elements_ids.size(); ++iel) {
        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
        elem->pad = -1;

        for (int i = 0; i < 8; i++) original_conn[i] = elem->nodes[i].id;

        //double z = (coords[original_conn[0]+2] + coords[original_conn[2]+2])*0.5;
        //if(z >= 0.0) continue;

        bool edge_list[12];
        int ed_cont = 0;

         for (int edge = 0; edge < 12; ++edge) {
            point[edge] = NULL;
            int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
            int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;

            //fprintf(fdbg, " Nodes : %d, %d\n", node1, node2);

            Edge2GNode[edge][0] = node1 <= node2 ? node1 : node2;
            Edge2GNode[edge][1] = node1 >= node2 ? node1 : node2;

            GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
            GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

            segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
            GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
            GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
            edge_list[edge] = false;
            if (list == NULL) continue;
            while (list) {
                GtsBBox *b = GTS_BBOX(list->data);
                point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
                if (point[edge]) {
                    //fprintf(fdbg, "edge, %d\n ", edge);
                    edge_list[edge] = true;
                    ed_cont++;
                    break;
                }
                list = list->next;
            }
            // 
            octant_node_t *gnode1 = (octant_node_t*) sc_array_index(&mesh->nodes, node1);
            octant_node_t *gnode2 = (octant_node_t*) sc_array_index(&mesh->nodes, node2);

            if (edge_list[edge] && gnode1->pad != 1 && gnode2->pad != 1) {

                //GtsPoint *p0 = point[edge];
                GtsPoint *p0 = gts_point_new(gts_point_class(),
                        coords[node1 * 3],
                        coords[node1 * 3 + 1],
                        coords[node1 * 3 + 2]);
                GtsPoint *p1 = gts_point_new(gts_point_class(),
                        coords[node2 * 3],
                        coords[node2 * 3 + 1],
                        coords[node2 * 3 + 2]);


                double d_c1 = gts_point_distance(point[edge], p0);
                double d_c2 = gts_point_distance(point[edge], p1);

                //printf("El: %d\n", elements_ids[iel]);
                //printf("nodes  : %d, %d\n", node1, node2);
                //printf("point : %f, %f, %f\n", point[edge]->x, point[edge]->y, point[edge]->z);
                //printf("no 0 : %f, %f, %f\n", coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
                //printf("no 1 : %f, %f, %f\n", coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);
                //printf("d_c0: %f, d_c1: %f\n", d_c0, d_c1);

                if (edge == 4 || edge == 5 || edge == 6 || edge == 7) {
                    if (coords[node2 * 3 + 2] == 0) {
                        gnode1->pad = 1;
                        coords[node1 * 3 + 2] = point[edge]->z;
                    } else if (coords[node1 * 3 + 2] == 0) {
                        gnode2->pad = 1;
                        coords[node2 * 3 + 2] = point[edge]->z;
                    } else {
                        if (d_c1 >= d_c2) {
                            coords[node2 * 3 + 2] = point[edge]->z;
                            gnode2->pad = 1;

                        } else {
                            coords[node1 * 3 + 2] = point[edge]->z;
                            gnode1->pad = 1;
                        }
                    }

                } else {
                    if (d_c1 >= d_c2) {
                        coords[node2 * 3] = point[edge]->x;
                        coords[node2 * 3 + 1] = point[edge]->y;
                        gnode2->pad = 1;

                    } else {
                        coords[node1 * 3] = point[edge]->x;
                        coords[node1 * 3 + 1] = point[edge]->y;
                        gnode1->pad = 1;
                    }
                }



            }

        }

        //check parallel faces
        for (int face = 0; face < 6; ++face) {
            face_intecepted[face] = false;
            for (int fe = 0; fe < 4; ++fe) {
                int edge = FaceEdgesMap[face][fe];
                if (point[edge] != NULL) {
                    face_intecepted[face] = true;
                    break;
                }
            }
            if (face_intecepted[face]) continue;
        }

        fprintf(fdbg, "Faces:  Elem. %d: %d %d %d %d %d %d\n", elements_ids[iel], face_intecepted[0],
                face_intecepted[1], face_intecepted[2],
                face_intecepted[3], face_intecepted[4], face_intecepted[5]);
        fprintf(fdbg, "Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
                edge_list[1], edge_list[2], edge_list[3],
                edge_list[4], edge_list[5], edge_list[6],
                edge_list[7], edge_list[8], edge_list[9],
                edge_list[10], edge_list[11], ed_cont);

        int n_parallel_faces = (face_intecepted[0] && face_intecepted[1]) +
                (face_intecepted[2] && face_intecepted[3]) +
                (face_intecepted[4] && face_intecepted[5]);

        if (n_parallel_faces == 2) {
            // Apply template 1.
            elem->pad = 1;

            int cut_edge[4];
            int connec_order_in[8];
            int connec_order_out1[8];
            int connec_order_out2[8];
#if 0

            if (((face_intecepted[0]) && (face_intecepted[1])) && (!face_intecepted[2] && !face_intecepted[3] &&
                    face_intecepted[4] && face_intecepted[5])) {
                //int cut_edge[4] = {3, 2, 11, 9};
                //int connec_order_in[8] = {0, 1, 4, 5, 2, 3, 6, 7};
                //int connec_order_out1[8] = {0, 1, 4, 5, 3, 2, 6, 7};
                //int connec_order_out2[8] = {0, 1, 4, 5, 2, 3, 6, 7};

                cut_edge[0] = 3;
                cut_edge[1] = 1;
                cut_edge[2] = 11;
                cut_edge[3] = 9;

                connec_order_in[0] = 0;
                connec_order_in[1] = 1;
                connec_order_in[2] = 4;
                connec_order_in[3] = 5;
                connec_order_in[4] = 2;
                connec_order_in[5] = 3;
                connec_order_in[6] = 6;
                connec_order_in[7] = 7;

                connec_order_out1[0] = 0;
                connec_order_out1[1] = 1;
                connec_order_out1[2] = 4;
                connec_order_out1[3] = 5;
                connec_order_out1[4] = 3;
                connec_order_out1[5] = 2;
                connec_order_out1[6] = 7;
                connec_order_out1[7] = 6;

                connec_order_out2[0] = 0;
                connec_order_out2[1] = 1;
                connec_order_out2[2] = 4;
                connec_order_out2[3] = 5;
                connec_order_out2[4] = 2;
                connec_order_out2[5] = 3;
                connec_order_out2[6] = 6;
                connec_order_out2[7] = 7;

            } else if ((face_intecepted[2]) && (face_intecepted[3]) && (!face_intecepted[0] && !face_intecepted[1] &&
                    face_intecepted[4] && face_intecepted[5])) {
                //int cut_edge[4] = {0, 2, 8, 10};
                //int connec_order_in[8] = {0, 3, 4, 7, 1, 2, 5, 6};
                //int connec_order_out1[8] = {0, 3, 4, 7, 1, 2, 5, 6};
                //int connec_order_out2[8] = {0, 3, 4, 7, 1, 2, 5, 6};

                cut_edge[0] = 0;
                cut_edge[1] = 2;
                cut_edge[2] = 8;
                cut_edge[3] = 10;

                connec_order_in[0] = 0;
                connec_order_in[1] = 3;
                connec_order_in[2] = 4;
                connec_order_in[3] = 7;
                connec_order_in[4] = 1;
                connec_order_in[5] = 2;
                connec_order_in[6] = 5;
                connec_order_in[7] = 6;

                connec_order_out1[0] = 0;
                connec_order_out1[1] = 3;
                connec_order_out1[2] = 4;
                connec_order_out1[3] = 7;
                connec_order_out1[4] = 1;
                connec_order_out1[5] = 2;
                connec_order_out1[6] = 5;
                connec_order_out1[7] = 6;

                connec_order_out2[0] = 0;
                connec_order_out2[1] = 3;
                connec_order_out2[2] = 4;
                connec_order_out2[3] = 7;
                connec_order_out2[4] = 1;
                connec_order_out2[5] = 2;
                connec_order_out2[6] = 5;
                connec_order_out2[7] = 6;

            } else if ((!face_intecepted[4]) && (!face_intecepted[5]) && (face_intecepted[0] && face_intecepted[1] &&
                    face_intecepted[2] && face_intecepted[3])) {
                //int cut_edge[4] = {4, 5, 6, 7};
                //int connec_order_in[8] = {0, 1, 2, 3, 4, 5, 6, 7};
                //int connec_order_out1[8] = {0, 1, 2, 3, 4, 5, 6, 7};
                //int connec_order_out2[8] = {0, 1, 2, 3, 4, 5, 6, 7};

                cut_edge[0] = 4;
                cut_edge[1] = 5;
                cut_edge[2] = 6;
                cut_edge[3] = 7;

                connec_order_in[0] = 0;
                connec_order_in[1] = 1;
                connec_order_in[2] = 2;
                connec_order_in[3] = 3;
                connec_order_in[4] = 4;
                connec_order_in[5] = 5;
                connec_order_in[6] = 6;
                connec_order_in[7] = 7;

                connec_order_out1[0] = 0;
                connec_order_out1[1] = 1;
                connec_order_out1[2] = 2;
                connec_order_out1[3] = 3;
                connec_order_out1[4] = 4;
                connec_order_out1[5] = 5;
                connec_order_out1[6] = 6;
                connec_order_out1[7] = 7;

                connec_order_out2[0] = 0;
                connec_order_out2[1] = 1;
                connec_order_out2[2] = 2;
                connec_order_out2[3] = 3;
                connec_order_out2[4] = 4;
                connec_order_out2[5] = 5;
                connec_order_out2[6] = 6;
                connec_order_out2[7] = 7;

            }

            GtsPoint *p0 = point[cut_edge[0]];
            GtsPoint *p1 = point[cut_edge[1]];
            GtsPoint *p2 = point[cut_edge[2]];
            GtsPoint *p3 = point[cut_edge[3]];

            g_assert(p0 != NULL);
            g_assert(p1 != NULL);
            g_assert(p2 != NULL);
            g_assert(p3 != NULL);

            conn_p[0] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p0, coords);
            conn_p[1] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p1, coords);
            conn_p[2] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p2, coords);
            conn_p[3] = AddPointOnEdge(Edge2GNode[cut_edge[3]], hash_nodes, mesh->local_n_nodes, p3, coords);

            octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

            elem1->nodes[connec_order_out1[0]].id = original_conn[connec_order_in[0]];
            elem1->nodes[connec_order_out1[1]].id = original_conn[connec_order_in[1]];
            elem1->nodes[connec_order_out1[2]].id = original_conn[connec_order_in[2]];
            elem1->nodes[connec_order_out1[3]].id = original_conn[connec_order_in[3]];

            elem1->nodes[connec_order_out1[4]].id = conn_p[0];
            elem1->nodes[connec_order_out1[5]].id = conn_p[1];
            elem1->nodes[connec_order_out1[6]].id = conn_p[2];
            elem1->nodes[connec_order_out1[7]].id = conn_p[3];

            octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);
            id_offset++;

            elem2->nodes[connec_order_out2[0]].id = conn_p[0];
            elem2->nodes[connec_order_out2[1]].id = conn_p[1];
            elem2->nodes[connec_order_out2[2]].id = conn_p[2];
            elem2->nodes[connec_order_out2[3]].id = conn_p[3];

            elem2->nodes[connec_order_out2[4]].id = original_conn[connec_order_in[4]];
            elem2->nodes[connec_order_out2[5]].id = original_conn[connec_order_in[5]];
            elem2->nodes[connec_order_out2[6]].id = original_conn[connec_order_in[6]];
            elem2->nodes[connec_order_out2[7]].id = original_conn[connec_order_in[7]];

            elem2->pad = elem->pad;
            elem2->level = elem->level;

            /*

                        if (((face_intecepted[0]) && (face_intecepted[1])) && (!face_intecepted[2] && !face_intecepted[3] &&
                                face_intecepted[4] && face_intecepted[5])) {

                            // New points are created splitting edges 0,2,8,10
                            // Two new elements are created
                            //  NP1 = Point splitting Edge1 

                            printf("primeiro\n");

                            GtsPoint *p0 = point[3];
                            GtsPoint *p1 = point[1];
                            GtsPoint *p2 = point[11];
                            GtsPoint *p3 = point[9];

                            g_assert(p0 != NULL);
                            g_assert(p1 != NULL);
                            g_assert(p2 != NULL);
                            g_assert(p3 != NULL);

                            conn_p[0] = AddPointOnEdge(Edge2GNode[3], hash_nodes, mesh->local_n_nodes, p0, coords);
                            conn_p[1] = AddPointOnEdge(Edge2GNode[1], hash_nodes, mesh->local_n_nodes, p1, coords);
                            conn_p[2] = AddPointOnEdge(Edge2GNode[11], hash_nodes, mesh->local_n_nodes, p2, coords);
                            conn_p[3] = AddPointOnEdge(Edge2GNode[9], hash_nodes, mesh->local_n_nodes, p3, coords);

                            // Remove the parent element form the list and
                            // introduce two new children elements
                            // Elem 1:
                            // N0, NP1, NP2, N3, N4, NP3, NP4, N7
                            // Elem 2:
                            // NP1, N1, N2, NP2, NP3, N5, N6, NP4
                            // Reajustando o array de elementos.
                            // Inserindo mais um elemento entre os indices elements_ids[i] e 
                            // elements_ids[i]+1

                            octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elem->id + id_offset);

                            elem1->nodes[0].id = original_conn[0];
                            elem1->nodes[1].id = original_conn[1];
                            elem1->nodes[4].id = original_conn[4];
                            elem1->nodes[5].id = original_conn[5];

                            elem1->nodes[3].id = conn_p[0];
                            elem1->nodes[2].id = conn_p[1];
                            elem1->nodes[6].id = conn_p[3];
                            elem1->nodes[7].id = conn_p[2];

                            //shift_n_position_from_index(&mesh->elements, 1, elem->id);

                            //octant_t *elem2 = (octant_t*) sc_array_index(&mesh->elements, elem->id + id_offset + 1);
                            octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);
                            id_offset++;

                            elem2->nodes[0].id = conn_p[0];
                            elem2->nodes[1].id = conn_p[1];
                            elem2->nodes[4].id = conn_p[2];
                            elem2->nodes[5].id = conn_p[3];
                            elem2->nodes[2].id = original_conn[2];
                            elem2->nodes[3].id = original_conn[3];
                            elem2->nodes[6].id = original_conn[6];
                            elem2->nodes[7].id = original_conn[7];

                            elem2->pad = elem->pad;
                            elem2->level = elem->level;


                        } else if ((face_intecepted[2]) && (face_intecepted[3]) && (!face_intecepted[0] && !face_intecepted[1] &&
                                face_intecepted[4] && face_intecepted[5])) {
                            // New points are created splitting edges 0,2,8,10
                            // Two new elements are created
                            //  NP1 = Point splitting Edge1  
                            GtsPoint *p0 = point[0];
                            GtsPoint *p1 = point[2];
                            GtsPoint *p2 = point[8];
                            GtsPoint *p3 = point[10];

                            printf("segundo %d\n",point[0]);

                            g_assert(p0 != NULL);
                            g_assert(p1 != NULL);
                            g_assert(p2 != NULL);
                            g_assert(p3 != NULL);

                            conn_p[0] = AddPointOnEdge(Edge2GNode[0 ], hash_nodes, mesh->local_n_nodes, p0, coords);
                            conn_p[1] = AddPointOnEdge(Edge2GNode[2 ], hash_nodes, mesh->local_n_nodes, p1, coords);
                            conn_p[2] = AddPointOnEdge(Edge2GNode[8 ], hash_nodes, mesh->local_n_nodes, p2, coords);
                            conn_p[3] = AddPointOnEdge(Edge2GNode[10], hash_nodes, mesh->local_n_nodes, p3, coords);

                            // Remove the parent element form the list and
                            // introduce two new children elements
                            // Elem 1:
                            // N0, NP1, NP2, N3, N4, NP3, NP4, N7
                            // Elem 2:
                            // NP1, N1, N2, NP2, NP3, N5, N6, NP4
                            // Reajustando o array de elementos.
                            // Inserindo mais um elemento entre os indices elements_ids[i] e 
                            // elements_ids[i]+1

                            octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elem->id + id_offset);

                            elem1->nodes[0].id = original_conn[0];
                            elem1->nodes[1].id = conn_p[0];
                            elem1->nodes[2].id = conn_p[1];
                            elem1->nodes[3].id = original_conn[3];
                            elem1->nodes[4].id = original_conn[4];
                            elem1->nodes[5].id = conn_p[2];
                            elem1->nodes[6].id = conn_p[3];
                            elem1->nodes[7].id = original_conn[7];

                            //shift_n_position_from_index(&mesh->elements, 1, elem->id);

                            //octant_t *elem2 = (octant_t*) sc_array_index(&mesh->elements, elem->id + id_offset + 1);
                            octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);
                            id_offset++;

                            elem2->nodes[0].id = conn_p[0];
                            elem2->nodes[1].id = original_conn[1];
                            elem2->nodes[2].id = original_conn[2];
                            elem2->nodes[3].id = conn_p[1];
                            elem2->nodes[4].id = conn_p[2];
                            elem2->nodes[5].id = original_conn[5];
                            elem2->nodes[6].id = original_conn[6];
                            elem2->nodes[7].id = conn_p[3];

                            elem2->pad = elem->pad;
                            elem2->level = elem->level;

                        } else if ((!face_intecepted[4]) && (!face_intecepted[5]) && (face_intecepted[0] && face_intecepted[1] &&
                                face_intecepted[2] && face_intecepted[3])) {
                            // New points are created splitting edges 0,2,8,10
                            // Two new elements are created
                            //  NP1 = Point splitting Edge1  
                            GtsPoint *p0 = point[4];
                            GtsPoint *p1 = point[5];
                            GtsPoint *p2 = point[6];
                            GtsPoint *p3 = point[7];

                            g_assert(p0 != NULL);
                            g_assert(p1 != NULL);
                            g_assert(p2 != NULL);
                            g_assert(p3 != NULL);

                            // Remove the parent element form the list and
                            // introduce two new children elements
                            // Elem 1:
                            // N0, NP1, NP2, N3, N4, NP3, NP4, N7
                            // Elem 2:
                            // NP1, N1, N2, NP2, NP3, N5, N6, NP4
                            // Reajustando o array de elementos.
                            // Inserindo mais um elemento entre os indices elements_ids[i] e 
                            // elements_ids[i]+1
                            conn_p[0] = AddPointOnEdge(Edge2GNode[4 ], hash_nodes, mesh->local_n_nodes, p0, coords);
                            conn_p[1] = AddPointOnEdge(Edge2GNode[5 ], hash_nodes, mesh->local_n_nodes, p1, coords);
                            conn_p[2] = AddPointOnEdge(Edge2GNode[6 ], hash_nodes, mesh->local_n_nodes, p2, coords);
                            conn_p[3] = AddPointOnEdge(Edge2GNode[7 ], hash_nodes, mesh->local_n_nodes, p3, coords);

                            fprintf(fdbg, " Conn_p %d, %d, %d, %d\n", conn_p[0], conn_p[1], conn_p[2], conn_p[3]);
                            octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elem->id + id_offset);

                            elem1->nodes[0].id = original_conn[0];
                            elem1->nodes[1].id = original_conn[1];
                            elem1->nodes[2].id = original_conn[2];
                            elem1->nodes[3].id = original_conn[3];
                            elem1->nodes[4].id = conn_p[0];
                            elem1->nodes[5].id = conn_p[1];
                            elem1->nodes[6].id = conn_p[2];
                            elem1->nodes[7].id = conn_p[3];

                            //shift_n_position_from_index(&mesh->elements, 1, elem->id);

                            //octant_t *elem2 = (octant_t*) sc_array_index(&mesh->elements, elem->id + id_offset + 1);
                            octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);
                            id_offset++;

                            elem2->nodes[0].id = conn_p[0];
                            elem2->nodes[1].id = conn_p[1];
                            elem2->nodes[2].id = conn_p[2];
                            elem2->nodes[3].id = conn_p[3];
                            elem2->nodes[4].id = original_conn[4];
                            elem2->nodes[5].id = original_conn[5];
                            elem2->nodes[6].id = original_conn[6];
                            elem2->nodes[7].id = original_conn[7];

                            elem2->pad = elem->pad;
                            elem2->level = elem->level;
                        }
             */
#endif

            for (int edge = 0; edge < 12; edge++) {

                if (point[edge]) gts_object_destroy(GTS_OBJECT(point[edge]));
                point[edge] = NULL;
            }

        }

        //check the diagonals in the surface and find the intersections
        GtsSegment * segments_s[12];
        GtsPoint * point_s[12];
        int Edge2GNode_s[12][2];

        for (int edge = 0; edge < 12; ++edge) {
            point_s[edge] = NULL;
            int node1 = elem->nodes[EdgeVerticesMap_surf_diagonal[edge][0]].id;
            int node2 = elem->nodes[EdgeVerticesMap_surf_diagonal[edge][1]].id;

            //fprintf(fdbg, " Nodes diagonal : %d, %d\n", node1, node2);

            Edge2GNode_s[edge][0] = node1 <= node2 ? node1 : node2;
            Edge2GNode_s[edge][1] = node1 >= node2 ? node1 : node2;

            GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
            GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

            segments_s[edge] = gts_segment_new(gts_segment_class(), v1, v2);
            GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments_s[edge]);
            GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);

            if (list == NULL) continue;
            while (list) {
                GtsBBox *b = GTS_BBOX(list->data);
                point_s[edge] = SegmentTriangleIntersection(segments_s[edge], GTS_TRIANGLE(b->bounded));
                if (point_s[edge]) {
                    //fprintf(fdbg, " edges diagonal : %d\n", edge);
                    break;
                }
                list = list->next;
            }
        }



        if (n_parallel_faces == 1) {
            // Apply template 1.
            elem->pad = 2;



#if 0
            int conn_t2[6];
            int cut_edge[4];
            int cut_edge_s[2];
            int connec_order_in1[8];
            int connec_order_out1[8];
            int connec_order_in2[8];
            int connec_order_out2[8];
            int connec_order_in3[8];
            int connec_order_out3[8];

            //if (((face_intecepted[0]) && (face_intecepted[1])) && (!face_intecepted[2] && !face_intecepted[3] &&
            //      face_intecepted[4] && face_intecepted[5])) {
            //int cut_edge[6] = {1, 3, 6, 7};
            //int connec_order_in[8] = {0, 1, 4, 5, 2, 3, 6, 7};
            //int connec_order_out1[8] = {0, 1, 4, 5, 3, 2, 6, 7};
            //int connec_order_out2[8] = {0, 1, 4, 5, 2, 3, 6, 7};
            //}

            //Edge 0
            if (((face_intecepted[0]) && (face_intecepted[1])) && (!face_intecepted[3] && !face_intecepted[4] &&
                    face_intecepted[2] && face_intecepted[5])) {
                printf("Entrou na 0!");

                cut_edge[0] = 3;
                cut_edge[1] = 1;
                cut_edge[2] = 4;
                cut_edge[3] = 5;

                cut_edge_s[0] = 4;
                cut_edge_s[1] = 6;

                connec_order_in1[0] = 2;
                connec_order_in1[1] = 3;
                connec_order_in1[2] = 6;
                connec_order_in1[3] = 7;
                connec_order_in1[4] = 0;
                connec_order_in1[5] = 1;
                connec_order_in1[6] = 4;
                connec_order_in1[7] = 5;

                connec_order_out1[0] = 2;
                connec_order_out1[1] = 3;
                connec_order_out1[2] = 6;
                connec_order_out1[3] = 7;
                connec_order_out1[4] = 0;
                connec_order_out1[5] = 1;
                connec_order_out1[6] = 4;
                connec_order_out1[7] = 5;

                connec_order_in2[0] = 0;
                connec_order_in2[1] = 1;
                connec_order_in2[2] = 2;
                connec_order_in2[3] = 3;
                connec_order_in2[4] = 4;
                connec_order_in2[5] = 5;
                connec_order_in2[6] = 6;
                connec_order_in2[7] = 7;

                connec_order_out2[0] = 2;
                connec_order_out2[1] = 3;
                connec_order_out2[2] = 5;
                connec_order_out2[3] = 4;
                connec_order_out2[4] = 4;
                connec_order_out2[5] = 5;
                connec_order_out2[6] = 6;
                connec_order_out2[7] = 7;

                connec_order_in3[0] = 0;
                connec_order_in3[1] = 1;
                connec_order_in3[2] = 2;
                connec_order_in3[3] = 3;
                connec_order_in3[4] = 4;
                connec_order_in3[5] = 5;
                connec_order_in3[6] = 6;
                connec_order_in3[7] = 7;

                connec_order_out3[0] = 0;
                connec_order_out3[1] = 1;
                connec_order_out3[2] = 1;
                connec_order_out3[3] = 0;
                connec_order_out3[4] = 2;
                connec_order_out3[5] = 3;
                connec_order_out3[6] = 5;
                connec_order_out3[7] = 4;
            }//Edge 1
            else if (((face_intecepted[1]) && (face_intecepted[2])) && (!face_intecepted[0] && !face_intecepted[4] &&
                    face_intecepted[3] && face_intecepted[5])) {
                printf("Entrou na 1!");

                cut_edge[0] = 0;
                cut_edge[1] = 2;
                cut_edge[2] = 5;
                cut_edge[3] = 6;

                cut_edge_s[0] = 1;
                cut_edge_s[1] = 3;

                connec_order_in1[0] = 0;
                connec_order_in1[1] = 3;
                connec_order_in1[2] = 4;
                connec_order_in1[3] = 7;
                connec_order_in1[4] = 1;
                connec_order_in1[5] = 2;
                connec_order_in1[6] = 5;
                connec_order_in1[7] = 6;

                connec_order_out1[0] = 0;
                connec_order_out1[1] = 3;
                connec_order_out1[2] = 4;
                connec_order_out1[3] = 7;
                connec_order_out1[4] = 0;
                connec_order_out1[5] = 1;
                connec_order_out1[6] = 4;
                connec_order_out1[7] = 5;

                connec_order_in2[0] = 0;
                connec_order_in2[1] = 1;
                connec_order_in2[2] = 2;
                connec_order_in2[3] = 3;
                connec_order_in2[4] = 4;
                connec_order_in2[5] = 5;
                connec_order_in2[6] = 6;
                connec_order_in2[7] = 7;

                connec_order_out2[0] = 4;
                connec_order_out2[1] = 2;
                connec_order_out2[2] = 3;
                connec_order_out2[3] = 5;
                connec_order_out2[4] = 4;
                connec_order_out2[5] = 5;
                connec_order_out2[6] = 6;
                connec_order_out2[7] = 7;

                connec_order_in3[0] = 1;
                connec_order_in3[1] = 2;
                connec_order_in3[2] = 0;
                connec_order_in3[3] = 3;
                connec_order_in3[4] = 4;
                connec_order_in3[5] = 5;
                connec_order_in3[6] = 6;
                connec_order_in3[7] = 7;

                connec_order_out3[0] = 1;
                connec_order_out3[1] = 2;
                connec_order_out3[2] = 0;
                connec_order_out3[3] = 1;
                connec_order_out3[4] = 4;
                connec_order_out3[5] = 2;
                connec_order_out3[6] = 3;
                connec_order_out3[7] = 5;
            }//Edge 2
            else if (((face_intecepted[0]) && (face_intecepted[1])) && (!face_intecepted[2] && !face_intecepted[4] &&
                    face_intecepted[3] && face_intecepted[5])) {
                printf("Entrou na 2!");

                cut_edge[0] = 3;
                cut_edge[1] = 1;
                cut_edge[2] = 7;
                cut_edge[3] = 6;

                cut_edge_s[0] = 5;
                cut_edge_s[1] = 7;

                connec_order_in1[0] = 0;
                connec_order_in1[1] = 1;
                connec_order_in1[2] = 4;
                connec_order_in1[3] = 5;
                connec_order_in1[4] = 2;
                connec_order_in1[5] = 3;
                connec_order_in1[6] = 6;
                connec_order_in1[7] = 7;

                connec_order_out1[0] = 0;
                connec_order_out1[1] = 1;
                connec_order_out1[2] = 4;
                connec_order_out1[3] = 5;
                connec_order_out1[4] = 1;
                connec_order_out1[5] = 0;
                connec_order_out1[6] = 5;
                connec_order_out1[7] = 4;

                connec_order_in2[0] = 0;
                connec_order_in2[1] = 1;
                connec_order_in2[2] = 2;
                connec_order_in2[3] = 3;
                connec_order_in2[4] = 4;
                connec_order_in2[5] = 5;
                connec_order_in2[6] = 6;
                connec_order_in2[7] = 7;

                connec_order_out2[0] = 4;
                connec_order_out2[1] = 5;
                connec_order_out2[2] = 3;
                connec_order_out2[3] = 2;
                connec_order_out2[4] = 4;
                connec_order_out2[5] = 5;
                connec_order_out2[6] = 6;
                connec_order_out2[7] = 7;

                connec_order_in3[0] = 2;
                connec_order_in3[1] = 3;
                connec_order_in3[2] = 0;
                connec_order_in3[3] = 1;
                connec_order_in3[4] = 4;
                connec_order_in3[5] = 5;
                connec_order_in3[6] = 6;
                connec_order_in3[7] = 7;

                connec_order_out3[0] = 2;
                connec_order_out3[1] = 3;
                connec_order_out3[2] = 0;
                connec_order_out3[3] = 1;
                connec_order_out3[4] = 4;
                connec_order_out3[5] = 5;
                connec_order_out3[6] = 3;
                connec_order_out3[7] = 2;
            }//Edge 3
            else if (((face_intecepted[0]) && (face_intecepted[2])) && (!face_intecepted[1] && !face_intecepted[4] &&
                    face_intecepted[3] && face_intecepted[5])) {
                printf("Entrou na 3!");

                cut_edge[0] = 0;
                cut_edge[1] = 2;
                cut_edge[2] = 4;
                cut_edge[3] = 7;

                cut_edge_s[0] = 0;
                cut_edge_s[1] = 2;

                connec_order_in1[0] = 1;
                connec_order_in1[1] = 2;
                connec_order_in1[2] = 5;
                connec_order_in1[3] = 6;
                connec_order_in1[4] = 0;
                connec_order_in1[5] = 3;
                connec_order_in1[6] = 4;
                connec_order_in1[7] = 7;

                connec_order_out1[0] = 1;
                connec_order_out1[1] = 2;
                connec_order_out1[2] = 5;
                connec_order_out1[3] = 6;
                connec_order_out1[4] = 0;
                connec_order_out1[5] = 1;
                connec_order_out1[6] = 4;
                connec_order_out1[7] = 5;

                connec_order_in2[0] = 0;
                connec_order_in2[1] = 1;
                connec_order_in2[2] = 2;
                connec_order_in2[3] = 3;
                connec_order_in2[4] = 4;
                connec_order_in2[5] = 5;
                connec_order_in2[6] = 6;
                connec_order_in2[7] = 7;

                connec_order_out2[0] = 2;
                connec_order_out2[1] = 4;
                connec_order_out2[2] = 5;
                connec_order_out2[3] = 3;
                connec_order_out2[4] = 4;
                connec_order_out2[5] = 5;
                connec_order_out2[6] = 6;
                connec_order_out2[7] = 7;

                connec_order_in3[0] = 0;
                connec_order_in3[1] = 3;
                connec_order_in3[2] = 1;
                connec_order_in3[3] = 2;
                connec_order_in3[4] = 4;
                connec_order_in3[5] = 5;
                connec_order_in3[6] = 6;
                connec_order_in3[7] = 7;

                connec_order_out3[0] = 0;
                connec_order_out3[1] = 3;
                connec_order_out3[2] = 0;
                connec_order_out3[3] = 1;
                connec_order_out3[4] = 2;
                connec_order_out3[5] = 4;
                connec_order_out3[6] = 5;
                connec_order_out3[7] = 3;

            }//Edge 4
            else if (((face_intecepted[0]) && (face_intecepted[2])) && (!face_intecepted[1] && !face_intecepted[3] &&
                    face_intecepted[4] && face_intecepted[5])) {
                printf("Entrou na 4!");

                cut_edge[0] = 0;
                cut_edge[1] = 3;
                cut_edge[2] = 11;
                cut_edge[3] = 8;

                cut_edge_s[0] = 10;
                cut_edge_s[1] = 8;

                connec_order_in1[0] = 1;
                connec_order_in1[1] = 2;
                connec_order_in1[2] = 5;
                connec_order_in1[3] = 6;
                connec_order_in1[4] = 0;
                connec_order_in1[5] = 3;
                connec_order_in1[6] = 4;
                connec_order_in1[7] = 7;

                connec_order_out1[0] = 1;
                connec_order_out1[1] = 2;
                connec_order_out1[2] = 5;
                connec_order_out1[3] = 6;
                connec_order_out1[4] = 0;
                connec_order_out1[5] = 4;
                connec_order_out1[6] = 3;
                connec_order_out1[7] = 5;

                connec_order_in2[0] = 0;
                connec_order_in2[1] = 1;
                connec_order_in2[2] = 4;
                connec_order_in2[3] = 5;
                connec_order_in2[4] = 2;
                connec_order_in2[5] = 3;
                connec_order_in2[6] = 6;
                connec_order_in2[7] = 7;

                connec_order_out2[0] = 1;
                connec_order_out2[1] = 4;
                connec_order_out2[2] = 2;
                connec_order_out2[3] = 5;
                connec_order_out2[4] = 2;
                connec_order_out2[5] = 3;
                connec_order_out2[6] = 6;
                connec_order_out2[7] = 7;

                connec_order_in3[0] = 0;
                connec_order_in3[1] = 4;
                connec_order_in3[2] = 1;
                connec_order_in3[3] = 2;
                connec_order_in3[4] = 3;
                connec_order_in3[5] = 5;
                connec_order_in3[6] = 6;
                connec_order_in3[7] = 7;

                connec_order_out3[0] = 0;
                connec_order_out3[1] = 4;
                connec_order_out3[2] = 0;
                connec_order_out3[3] = 4;
                connec_order_out3[4] = 1;
                connec_order_out3[5] = 3;
                connec_order_out3[6] = 5;
                connec_order_out3[7] = 2;

            }//Edge 5
            else if (((face_intecepted[1]) && (face_intecepted[2])) && (!face_intecepted[0] && !face_intecepted[3] &&
                    face_intecepted[4] && face_intecepted[5])) {
                printf("Entrou na 5!");

                cut_edge[0] = 1;
                cut_edge[1] = 0;
                cut_edge[2] = 8;
                cut_edge[3] = 9;

                cut_edge_s[0] = 11;
                cut_edge_s[1] = 9;

                connec_order_in1[0] = 0;
                connec_order_in1[1] = 3;
                connec_order_in1[2] = 4;
                connec_order_in1[3] = 7;
                connec_order_in1[4] = 1;
                connec_order_in1[5] = 2;
                connec_order_in1[6] = 5;
                connec_order_in1[7] = 6;

                connec_order_out1[0] = 0;
                connec_order_out1[1] = 3;
                connec_order_out1[2] = 4;
                connec_order_out1[3] = 7;
                connec_order_out1[4] = 1;
                connec_order_out1[5] = 4;
                connec_order_out1[6] = 2;
                connec_order_out1[7] = 5;

                connec_order_in2[0] = 0;
                connec_order_in2[1] = 1;
                connec_order_in2[2] = 4;
                connec_order_in2[3] = 5;
                connec_order_in2[4] = 2;
                connec_order_in2[5] = 3;
                connec_order_in2[6] = 6;
                connec_order_in2[7] = 7;

                connec_order_out2[0] = 4;
                connec_order_out2[1] = 0;
                connec_order_out2[2] = 5;
                connec_order_out2[3] = 3;
                connec_order_out2[4] = 2;
                connec_order_out2[5] = 3;
                connec_order_out2[6] = 6;
                connec_order_out2[7] = 7;

                connec_order_in3[0] = 1;
                connec_order_in3[1] = 5;
                connec_order_in3[2] = 0;
                connec_order_in3[3] = 2;
                connec_order_in3[4] = 3;
                connec_order_in3[5] = 4;
                connec_order_in3[6] = 6;
                connec_order_in3[7] = 7;

                connec_order_out3[0] = 1;
                connec_order_out3[1] = 5;
                connec_order_out3[2] = 1;
                connec_order_out3[3] = 0;
                connec_order_out3[4] = 4;
                connec_order_out3[5] = 2;
                connec_order_out3[6] = 3;
                connec_order_out3[7] = 5;

            }//Edge 6
            else if (((face_intecepted[1]) && (face_intecepted[3])) && (!face_intecepted[0] && !face_intecepted[2] &&
                    face_intecepted[4] && face_intecepted[5])) {
                printf("Entrou na 6!");

                cut_edge[0] = 2;
                cut_edge[1] = 10;
                cut_edge[2] = 1;
                cut_edge[3] = 9;

                cut_edge_s[0] = 10;
                cut_edge_s[1] = 8;

                connec_order_in1[0] = 0;
                connec_order_in1[1] = 1;
                connec_order_in1[2] = 4;
                connec_order_in1[3] = 5;
                connec_order_in1[4] = 2;
                connec_order_in1[5] = 3;
                connec_order_in1[6] = 6;
                connec_order_in1[7] = 7;

                connec_order_out1[0] = 0;
                connec_order_out1[1] = 1;
                connec_order_out1[2] = 4;
                connec_order_out1[3] = 5;
                connec_order_out1[4] = 2;
                connec_order_out1[5] = 4;
                connec_order_out1[6] = 3;
                connec_order_out1[7] = 5;

                connec_order_in2[0] = 1;
                connec_order_in2[1] = 2;
                connec_order_in2[2] = 5;
                connec_order_in2[3] = 6;
                connec_order_in2[4] = 0;
                connec_order_in2[5] = 3;
                connec_order_in2[6] = 4;
                connec_order_in2[7] = 7;

                connec_order_out2[0] = 4;
                connec_order_out2[1] = 0;
                connec_order_out2[2] = 5;
                connec_order_out2[3] = 1;
                connec_order_out2[4] = 0;
                connec_order_out2[5] = 3;
                connec_order_out2[6] = 4;
                connec_order_out2[7] = 7;

                connec_order_in3[0] = 2;
                connec_order_in3[1] = 6;
                connec_order_in3[2] = 0;
                connec_order_in3[3] = 1;
                connec_order_in3[4] = 3;
                connec_order_in3[5] = 4;
                connec_order_in3[6] = 5;
                connec_order_in3[7] = 7;

                connec_order_out3[0] = 2;
                connec_order_out3[1] = 6;
                connec_order_out3[2] = 4;
                connec_order_out3[3] = 2;
                connec_order_out3[4] = 0;
                connec_order_out3[5] = 5;
                connec_order_out3[6] = 3;
                connec_order_out3[7] = 1;

            }//Edge 7
            else if (((face_intecepted[0]) && (face_intecepted[3])) && (!face_intecepted[1] && !face_intecepted[2] &&
                    face_intecepted[4] && face_intecepted[5])) {
                printf("Entrou na 7!");

                cut_edge[0] = 3;
                cut_edge[1] = 11;
                cut_edge[2] = 2;
                cut_edge[3] = 10;

                cut_edge_s[0] = 11;
                cut_edge_s[1] = 9;

                connec_order_in1[0] = 0;
                connec_order_in1[1] = 1;
                connec_order_in1[2] = 4;
                connec_order_in1[3] = 5;
                connec_order_in1[4] = 2;
                connec_order_in1[5] = 3;
                connec_order_in1[6] = 6;
                connec_order_in1[7] = 7;

                connec_order_out1[0] = 0;
                connec_order_out1[1] = 1;
                connec_order_out1[2] = 4;
                connec_order_out1[3] = 5;
                connec_order_out1[4] = 4;
                connec_order_out1[5] = 0;
                connec_order_out1[6] = 5;
                connec_order_out1[7] = 1;

                connec_order_in2[0] = 0;
                connec_order_in2[1] = 3;
                connec_order_in2[2] = 4;
                connec_order_in2[3] = 7;
                connec_order_in2[4] = 1;
                connec_order_in2[5] = 2;
                connec_order_in2[6] = 5;
                connec_order_in2[7] = 6;

                connec_order_out2[0] = 4;
                connec_order_out2[1] = 2;
                connec_order_out2[2] = 5;
                connec_order_out2[3] = 3;
                connec_order_out2[4] = 1;
                connec_order_out2[5] = 2;
                connec_order_out2[6] = 5;
                connec_order_out2[7] = 6;

                connec_order_in3[0] = 3;
                connec_order_in3[1] = 7;
                connec_order_in3[2] = 0;
                connec_order_in3[3] = 1;
                connec_order_in3[4] = 2;
                connec_order_in3[5] = 4;
                connec_order_in3[6] = 5;
                connec_order_in3[7] = 6;

                connec_order_out3[0] = 3;
                connec_order_out3[1] = 7;
                connec_order_out3[2] = 0;
                connec_order_out3[3] = 4;
                connec_order_out3[4] = 2;
                connec_order_out3[5] = 1;
                connec_order_out3[6] = 5;
                connec_order_out3[7] = 3;

            }//Edge 8
            else if (((face_intecepted[0]) && (face_intecepted[1])) && (!face_intecepted[3] && !face_intecepted[5] &&
                    face_intecepted[2] && face_intecepted[4])) {
                printf("Entrou na 8!");

                cut_edge[0] = 11;
                cut_edge[1] = 9;
                cut_edge[2] = 4;
                cut_edge[3] = 5;

                cut_edge_s[0] = 5;
                cut_edge_s[1] = 7;

                connec_order_in1[0] = 2;
                connec_order_in1[1] = 3;
                connec_order_in1[2] = 6;
                connec_order_in1[3] = 7;
                connec_order_in1[4] = 0;
                connec_order_in1[5] = 1;
                connec_order_in1[6] = 4;
                connec_order_in1[7] = 5;

                connec_order_out1[0] = 2;
                connec_order_out1[1] = 3;
                connec_order_out1[2] = 6;
                connec_order_out1[3] = 7;
                connec_order_out1[4] = 4;
                connec_order_out1[5] = 5;
                connec_order_out1[6] = 0;
                connec_order_out1[7] = 1;

                connec_order_in2[0] = 4;
                connec_order_in2[1] = 5;
                connec_order_in2[2] = 6;
                connec_order_in2[3] = 7;
                connec_order_in2[4] = 0;
                connec_order_in2[5] = 1;
                connec_order_in2[6] = 2;
                connec_order_in2[7] = 3;

                connec_order_out2[0] = 2;
                connec_order_out2[1] = 3;
                connec_order_out2[2] = 5;
                connec_order_out2[3] = 4;
                connec_order_out2[4] = 0;
                connec_order_out2[5] = 1;
                connec_order_out2[6] = 2;
                connec_order_out2[7] = 3;

                connec_order_in3[0] = 4;
                connec_order_in3[1] = 5;
                connec_order_in3[2] = 0;
                connec_order_in3[3] = 1;
                connec_order_in3[4] = 2;
                connec_order_in3[5] = 3;
                connec_order_in3[6] = 6;
                connec_order_in3[7] = 7;

                connec_order_out3[0] = 4;
                connec_order_out3[1] = 5;
                connec_order_out3[2] = 2;
                connec_order_out3[3] = 3;
                connec_order_out3[4] = 5;
                connec_order_out3[5] = 4;
                connec_order_out3[6] = 1;
                connec_order_out3[7] = 0;


            }//Edge 9
            else if (((face_intecepted[1]) && (face_intecepted[2])) && (!face_intecepted[0] && !face_intecepted[5] &&
                    face_intecepted[3] && face_intecepted[4])) {
                printf("Entrou na 9!");

                cut_edge[0] = 8;
                cut_edge[1] = 10;
                cut_edge[2] = 5;
                cut_edge[3] = 6;

                cut_edge_s[0] = 0;
                cut_edge_s[1] = 2;

                connec_order_in1[0] = 0;
                connec_order_in1[1] = 3;
                connec_order_in1[2] = 4;
                connec_order_in1[3] = 7;
                connec_order_in1[4] = 1;
                connec_order_in1[5] = 2;
                connec_order_in1[6] = 5;
                connec_order_in1[7] = 6;

                connec_order_out1[0] = 0;
                connec_order_out1[1] = 3;
                connec_order_out1[2] = 4;
                connec_order_out1[3] = 7;
                connec_order_out1[4] = 4;
                connec_order_out1[5] = 5;
                connec_order_out1[6] = 0;
                connec_order_out1[7] = 1;

                connec_order_in2[0] = 4;
                connec_order_in2[1] = 5;
                connec_order_in2[2] = 6;
                connec_order_in2[3] = 7;
                connec_order_in2[4] = 0;
                connec_order_in2[5] = 1;
                connec_order_in2[6] = 2;
                connec_order_in2[7] = 3;

                connec_order_out2[0] = 4;
                connec_order_out2[1] = 2;
                connec_order_out2[2] = 3;
                connec_order_out2[3] = 5;
                connec_order_out2[4] = 0;
                connec_order_out2[5] = 1;
                connec_order_out2[6] = 2;
                connec_order_out2[7] = 3;

                connec_order_in3[0] = 5;
                connec_order_in3[1] = 6;
                connec_order_in3[2] = 0;
                connec_order_in3[3] = 1;
                connec_order_in3[4] = 2;
                connec_order_in3[5] = 3;
                connec_order_in3[6] = 4;
                connec_order_in3[7] = 7;

                connec_order_out3[0] = 5;
                connec_order_out3[1] = 6;
                connec_order_out3[2] = 4;
                connec_order_out3[3] = 2;
                connec_order_out3[4] = 3;
                connec_order_out3[5] = 5;
                connec_order_out3[6] = 0;
                connec_order_out3[7] = 1;

            }//Edge 10
            else if (((face_intecepted[0]) && (face_intecepted[1])) && (!face_intecepted[2] && !face_intecepted[5] &&
                    face_intecepted[3] && face_intecepted[4])) {
                printf("Entrou na 10!");

                cut_edge[0] = 9;
                cut_edge[1] = 11;
                cut_edge[2] = 6;
                cut_edge[3] = 7;

                cut_edge_s[0] = 6;
                cut_edge_s[1] = 4;

                connec_order_in1[0] = 0;
                connec_order_in1[1] = 1;
                connec_order_in1[2] = 4;
                connec_order_in1[3] = 5;
                connec_order_in1[4] = 2;
                connec_order_in1[5] = 3;
                connec_order_in1[6] = 6;
                connec_order_in1[7] = 7;

                connec_order_out1[0] = 0;
                connec_order_out1[1] = 1;
                connec_order_out1[2] = 4;
                connec_order_out1[3] = 5;
                connec_order_out1[4] = 4;
                connec_order_out1[5] = 5;
                connec_order_out1[6] = 0;
                connec_order_out1[7] = 1;

                connec_order_in2[0] = 4;
                connec_order_in2[1] = 5;
                connec_order_in2[2] = 6;
                connec_order_in2[3] = 7;
                connec_order_in2[4] = 0;
                connec_order_in2[5] = 1;
                connec_order_in2[6] = 2;
                connec_order_in2[7] = 3;

                connec_order_out2[0] = 5;
                connec_order_out2[1] = 4;
                connec_order_out2[2] = 2;
                connec_order_out2[3] = 3;
                connec_order_out2[4] = 0;
                connec_order_out2[5] = 1;
                connec_order_out2[6] = 2;
                connec_order_out2[7] = 3;

                connec_order_in3[0] = 6;
                connec_order_in3[1] = 7;
                connec_order_in3[2] = 0;
                connec_order_in3[3] = 1;
                connec_order_in3[4] = 2;
                connec_order_in3[5] = 3;
                connec_order_in3[6] = 4;
                connec_order_in3[7] = 5;

                connec_order_out3[0] = 6;
                connec_order_out3[1] = 7;
                connec_order_out3[2] = 5;
                connec_order_out3[3] = 4;
                connec_order_out3[4] = 2;
                connec_order_out3[5] = 3;
                connec_order_out3[6] = 1;
                connec_order_out3[7] = 0;

            }//Edge 11
            else if (((face_intecepted[0]) && (face_intecepted[2])) && (!face_intecepted[1] && !face_intecepted[5] &&
                    face_intecepted[3] && face_intecepted[4])) {
                printf("Entrou na 11!");

                cut_edge[0] = 8;
                cut_edge[1] = 10;
                cut_edge[2] = 4;
                cut_edge[3] = 7;

                cut_edge_s[0] = 1;
                cut_edge_s[1] = 3;

                connec_order_in1[0] = 1;
                connec_order_in1[1] = 2;
                connec_order_in1[2] = 5;
                connec_order_in1[3] = 6;
                connec_order_in1[4] = 0;
                connec_order_in1[5] = 3;
                connec_order_in1[6] = 4;
                connec_order_in1[7] = 7;

                connec_order_out1[0] = 1;
                connec_order_out1[1] = 2;
                connec_order_out1[2] = 5;
                connec_order_out1[3] = 6;
                connec_order_out1[4] = 4;
                connec_order_out1[5] = 5;
                connec_order_out1[6] = 0;
                connec_order_out1[7] = 1;

                connec_order_in2[0] = 4;
                connec_order_in2[1] = 5;
                connec_order_in2[2] = 6;
                connec_order_in2[3] = 7;
                connec_order_in2[4] = 0;
                connec_order_in2[5] = 1;
                connec_order_in2[6] = 2;
                connec_order_in2[7] = 3;

                connec_order_out2[0] = 2;
                connec_order_out2[1] = 4;
                connec_order_out2[2] = 5;
                connec_order_out2[3] = 3;
                connec_order_out2[4] = 0;
                connec_order_out2[5] = 1;
                connec_order_out2[6] = 2;
                connec_order_out2[7] = 3;

                connec_order_in3[0] = 4;
                connec_order_in3[1] = 7;
                connec_order_in3[2] = 0;
                connec_order_in3[3] = 1;
                connec_order_in3[4] = 2;
                connec_order_in3[5] = 3;
                connec_order_in3[6] = 5;
                connec_order_in3[7] = 6;

                connec_order_out3[0] = 4;
                connec_order_out3[1] = 7;
                connec_order_out3[2] = 2;
                connec_order_out3[3] = 4;
                connec_order_out3[4] = 5;
                connec_order_out3[5] = 3;
                connec_order_out3[6] = 0;
                connec_order_out3[7] = 1;

            }

            GtsPoint *p0 = point[cut_edge[0]];
            GtsPoint *p1 = point[cut_edge[1]];
            GtsPoint *p2 = point[cut_edge[2]];
            GtsPoint *p3 = point[cut_edge[3]];

            GtsPoint *p4 = point_s[cut_edge_s[0]];
            GtsPoint *p5 = point_s[cut_edge_s[1]];

            g_assert(p0 != NULL);
            g_assert(p1 != NULL);
            g_assert(p2 != NULL);
            g_assert(p3 != NULL);
            g_assert(p4 != NULL);
            g_assert(p5 != NULL);

            conn_t2[0] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p0, coords);
            conn_t2[1] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p1, coords);
            conn_t2[2] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p2, coords);
            conn_t2[3] = AddPointOnEdge(Edge2GNode[cut_edge[3]], hash_nodes, mesh->local_n_nodes, p3, coords);

            // add 2 extra points in the surface
            conn_t2[4] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[0]], hash_nodes, mesh->local_n_nodes, p4, coords);
            conn_t2[5] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[1]], hash_nodes, mesh->local_n_nodes, p5, coords);

            octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

            elem1->nodes[connec_order_in1[0]].id = original_conn[connec_order_out1[0]];
            elem1->nodes[connec_order_in1[1]].id = original_conn[connec_order_out1[1]];
            elem1->nodes[connec_order_in1[2]].id = original_conn[connec_order_out1[2]];
            elem1->nodes[connec_order_in1[3]].id = original_conn[connec_order_out1[3]];

            elem1->nodes[connec_order_in1[4]].id = conn_t2[connec_order_out1[4]];
            elem1->nodes[connec_order_in1[5]].id = conn_t2[connec_order_out1[5]];
            elem1->nodes[connec_order_in1[6]].id = conn_t2[connec_order_out1[6]];
            elem1->nodes[connec_order_in1[7]].id = conn_t2[connec_order_out1[7]];

            octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

            elem2->nodes[connec_order_in2[0]].id = conn_t2[connec_order_out2[0]];
            elem2->nodes[connec_order_in2[1]].id = conn_t2[connec_order_out2[1]];
            elem2->nodes[connec_order_in2[2]].id = conn_t2[connec_order_out2[2]];
            elem2->nodes[connec_order_in2[3]].id = conn_t2[connec_order_out2[3]];

            elem2->nodes[connec_order_in2[4]].id = original_conn[connec_order_out2[4]];
            elem2->nodes[connec_order_in2[5]].id = original_conn[connec_order_out2[5]];
            elem2->nodes[connec_order_in2[6]].id = original_conn[connec_order_out2[6]];
            elem2->nodes[connec_order_in2[7]].id = original_conn[connec_order_out2[7]];

            elem2->pad = 3;
            elem2->level = elem->level;

            octant_t* elem3 = (octant_t*) sc_array_push(&mesh->elements);

            elem3->nodes[connec_order_in3[0]].id = original_conn[connec_order_out3[0]];
            elem3->nodes[connec_order_in3[1]].id = original_conn[connec_order_out3[1]];

            elem3->nodes[connec_order_in3[2]].id = conn_t2[connec_order_out3[2]];
            elem3->nodes[connec_order_in3[3]].id = conn_t2[connec_order_out3[3]];
            elem3->nodes[connec_order_in3[4]].id = conn_t2[connec_order_out3[4]];
            elem3->nodes[connec_order_in3[5]].id = conn_t2[connec_order_out3[5]];
            elem3->nodes[connec_order_in3[6]].id = conn_t2[connec_order_out3[6]];
            elem3->nodes[connec_order_in3[7]].id = conn_t2[connec_order_out3[7]];

            elem3->pad = 4;
            elem3->level = elem->level;

#endif
            for (int edge = 0; edge < 12; edge++) {
                if (point_s[edge]) gts_object_destroy(GTS_OBJECT(point_s[edge]));
                point_s[edge] = NULL;
            }

            for (int edge = 0; edge < 12; edge++) {
                if (point[edge]) gts_object_destroy(GTS_OBJECT(point[edge]));
                point[edge] = NULL;
            }

        }



        //check the diagonals in the volume and find the intersections
        GtsSegment * segments_v[4];
        GtsPoint * point_v[4];
        int Edge2GNode_v[4][2];

        for (int edge = 0; edge < 4; ++edge) {
            point_v[edge] = NULL;
            int node1 = elem->nodes[EdgeVerticesMap_vol_diagonal[edge][0]].id;
            int node2 = elem->nodes[EdgeVerticesMap_vol_diagonal[edge][1]].id;

            //fprintf(fdbg, " Nodes diagonal vol : %d, %d\n", node1, node2);

            Edge2GNode_v[edge][0] = node1 <= node2 ? node1 : node2;
            Edge2GNode_v[edge][1] = node1 >= node2 ? node1 : node2;

            GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
            GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

            segments_v[edge] = gts_segment_new(gts_segment_class(), v1, v2);
            GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments_v[edge]);
            GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);

            if (list == NULL) continue;
            while (list) {
                GtsBBox *b = GTS_BBOX(list->data);
                point_v[edge] = SegmentTriangleIntersection(segments_v[edge], GTS_TRIANGLE(b->bounded));
                if (point_v[edge]) {
                    //fprintf(fdbg, " edges diagonal vol: %d\n", edge);
                    //           printf(" edges diagonal vol: %d\n", edge);

                    break;
                }
                list = list->next;
            }
        }


        //        printf("n_parallel_faces %d, ed_cont, %d\n", n_parallel_faces, ed_cont);
        if (n_parallel_faces == 0 && ed_cont == 3) {
            // Apply template 1.
            elem->pad = 5;

#if 0

            int conn_t2[7];
            int cut_edge[3];
            int cut_edge_s[3];
            int cut_edge_v[1];
            int connec_order_in1[8];
            int connec_order_out1[8];
            int connec_order_in2[8];
            int connec_order_out2[8];
            int connec_order_in3[8];
            int connec_order_out3[8];
            int connec_order_in4[8];
            int connec_order_out4[8];


            //Corner 0
            if ((face_intecepted[0]) && (!face_intecepted[1]) && (!face_intecepted[3] && !face_intecepted[4] &&
                    face_intecepted[2] && face_intecepted[5])) {
                printf("Entrou no 0!");

                cut_edge[0] = 4;
                cut_edge[1] = 0;
                cut_edge[2] = 3;

                cut_edge_s[0] = 0;
                cut_edge_s[1] = 4;
                cut_edge_s[2] = 10;

                cut_edge_v[0] = 0;

                connec_order_in1[0] = 4;
                connec_order_in1[1] = 5;
                connec_order_in1[2] = 6;
                connec_order_in1[3] = 7;
                connec_order_in1[4] = 0;
                connec_order_in1[5] = 1;
                connec_order_in1[6] = 2;
                connec_order_in1[7] = 3;

                connec_order_out1[0] = 4;
                connec_order_out1[1] = 5;
                connec_order_out1[2] = 6;
                connec_order_out1[3] = 7;
                connec_order_out1[4] = 0;
                connec_order_out1[5] = 3;
                connec_order_out1[6] = 6;
                connec_order_out1[7] = 4;

                connec_order_in2[0] = 0;
                connec_order_in2[1] = 3;
                connec_order_in2[2] = 4;
                connec_order_in2[3] = 7;
                connec_order_in2[4] = 1;
                connec_order_in2[5] = 2;
                connec_order_in2[6] = 5;
                connec_order_in2[7] = 6;

                connec_order_out2[0] = 1;
                connec_order_out2[1] = 5;
                connec_order_out2[2] = 3;
                connec_order_out2[3] = 6;
                connec_order_out2[4] = 1;
                connec_order_out2[5] = 2;
                connec_order_out2[6] = 5;
                connec_order_out2[7] = 6;

                connec_order_in3[0] = 0;
                connec_order_in3[1] = 1;
                connec_order_in3[2] = 4;
                connec_order_in3[3] = 5;
                connec_order_in3[4] = 2;
                connec_order_in3[5] = 3;
                connec_order_in3[6] = 6;
                connec_order_in3[7] = 7;

                connec_order_out3[0] = 2;
                connec_order_out3[1] = 5;
                connec_order_out3[2] = 4;
                connec_order_out3[3] = 6;
                connec_order_out3[4] = 2;
                connec_order_out3[5] = 3;
                connec_order_out3[6] = 6;
                connec_order_out3[7] = 7;

                connec_order_in4[0] = 0;
                connec_order_in4[1] = 1;
                connec_order_in4[2] = 2;
                connec_order_in4[3] = 3;
                connec_order_in4[4] = 4;
                connec_order_in4[5] = 5;
                connec_order_in4[6] = 6;
                connec_order_in4[7] = 7;

                connec_order_out4[0] = 0;
                connec_order_out4[1] = 1;
                connec_order_out4[2] = 5;
                connec_order_out4[3] = 2;
                connec_order_out4[4] = 0;
                connec_order_out4[5] = 3;
                connec_order_out4[6] = 6;
                connec_order_out4[7] = 4;
            }

            //Corner 1
            if ((face_intecepted[1]) && (!face_intecepted[0]) && (!face_intecepted[3] && !face_intecepted[4] &&
                    face_intecepted[2] && face_intecepted[5])) {
                printf("Entrou no 1!");

                cut_edge[0] = 5;
                cut_edge[1] = 1;
                cut_edge[2] = 0;

                cut_edge_s[0] = 6;
                cut_edge_s[1] = 1;
                cut_edge_s[2] = 11;

                cut_edge_v[0] = 1;

                connec_order_in1[0] = 4;
                connec_order_in1[1] = 5;
                connec_order_in1[2] = 6;
                connec_order_in1[3] = 7;
                connec_order_in1[4] = 0;
                connec_order_in1[5] = 1;
                connec_order_in1[6] = 2;
                connec_order_in1[7] = 3;

                connec_order_out1[0] = 4;
                connec_order_out1[1] = 5;
                connec_order_out1[2] = 6;
                connec_order_out1[3] = 7;
                connec_order_out1[4] = 4;
                connec_order_out1[5] = 0;
                connec_order_out1[6] = 3;
                connec_order_out1[7] = 6;

                connec_order_in2[0] = 0;
                connec_order_in2[1] = 1;
                connec_order_in2[2] = 4;
                connec_order_in2[3] = 5;
                connec_order_in2[4] = 2;
                connec_order_in2[5] = 3;
                connec_order_in2[6] = 6;
                connec_order_in2[7] = 7;

                connec_order_out2[0] = 5;
                connec_order_out2[1] = 1;
                connec_order_out2[2] = 6;
                connec_order_out2[3] = 3;
                connec_order_out2[4] = 2;
                connec_order_out2[5] = 3;
                connec_order_out2[6] = 6;
                connec_order_out2[7] = 7;

                connec_order_in3[0] = 1;
                connec_order_in3[1] = 2;
                connec_order_in3[2] = 5;
                connec_order_in3[3] = 6;
                connec_order_in3[4] = 0;
                connec_order_in3[5] = 3;
                connec_order_in3[6] = 4;
                connec_order_in3[7] = 7;

                connec_order_out3[0] = 2;
                connec_order_out3[1] = 5;
                connec_order_out3[2] = 4;
                connec_order_out3[3] = 6;
                connec_order_out3[4] = 0;
                connec_order_out3[5] = 3;
                connec_order_out3[6] = 4;
                connec_order_out3[7] = 7;

                connec_order_in4[0] = 1;
                connec_order_in4[1] = 0;
                connec_order_in4[2] = 2;
                connec_order_in4[3] = 3;
                connec_order_in4[4] = 4;
                connec_order_in4[5] = 5;
                connec_order_in4[6] = 6;
                connec_order_in4[7] = 7;

                connec_order_out4[0] = 1;
                connec_order_out4[1] = 2;
                connec_order_out4[2] = 1;
                connec_order_out4[3] = 5;
                connec_order_out4[4] = 4;
                connec_order_out4[5] = 0;
                connec_order_out4[6] = 3;
                connec_order_out4[7] = 6;

            }

            //Corner 2
            if ((face_intecepted[1]) && (!face_intecepted[0]) && (!face_intecepted[2] && !face_intecepted[4] &&
                    face_intecepted[3] && face_intecepted[5])) {
                printf("Entrou no 2!");

                cut_edge[0] = 6;
                cut_edge[1] = 2;
                cut_edge[2] = 1;

                cut_edge_s[0] = 3;
                cut_edge_s[1] = 7;
                cut_edge_s[2] = 10;

                cut_edge_v[0] = 2;

                connec_order_in1[0] = 4;
                connec_order_in1[1] = 5;
                connec_order_in1[2] = 6;
                connec_order_in1[3] = 7;
                connec_order_in1[4] = 1;
                connec_order_in1[5] = 2;
                connec_order_in1[6] = 3;
                connec_order_in1[7] = 4;

                connec_order_out1[0] = 4;
                connec_order_out1[1] = 5;
                connec_order_out1[2] = 6;
                connec_order_out1[3] = 7;
                connec_order_out1[4] = 4;
                connec_order_out1[5] = 0;
                connec_order_out1[6] = 3;
                connec_order_out1[7] = 6;

                connec_order_in2[0] = 1;
                connec_order_in2[1] = 2;
                connec_order_in2[2] = 5;
                connec_order_in2[3] = 6;
                connec_order_in2[4] = 0;
                connec_order_in2[5] = 3;
                connec_order_in2[6] = 4;
                connec_order_in2[7] = 7;

                connec_order_out2[0] = 5;
                connec_order_out2[1] = 1;
                connec_order_out2[2] = 6;
                connec_order_out2[3] = 3;
                connec_order_out2[4] = 0;
                connec_order_out2[5] = 3;
                connec_order_out2[6] = 4;
                connec_order_out2[7] = 7;

                connec_order_in3[0] = 2;
                connec_order_in3[1] = 3;
                connec_order_in3[2] = 6;
                connec_order_in3[3] = 7;
                connec_order_in3[4] = 0;
                connec_order_in3[5] = 1;
                connec_order_in3[6] = 4;
                connec_order_in3[7] = 5;

                connec_order_out3[0] = 2;
                connec_order_out3[1] = 5;
                connec_order_out3[2] = 4;
                connec_order_out3[3] = 6;
                connec_order_out3[4] = 0;
                connec_order_out3[5] = 1;
                connec_order_out3[6] = 4;
                connec_order_out3[7] = 5;

                connec_order_in4[0] = 2;
                connec_order_in4[1] = 0;
                connec_order_in4[2] = 1;
                connec_order_in4[3] = 3;
                connec_order_in4[4] = 4;
                connec_order_in4[5] = 5;
                connec_order_in4[6] = 6;
                connec_order_in4[7] = 7;

                connec_order_out4[0] = 2;
                connec_order_out4[1] = 5;
                connec_order_out4[2] = 2;
                connec_order_out4[3] = 1;
                connec_order_out4[4] = 6;
                connec_order_out4[5] = 4;
                connec_order_out4[6] = 0;
                connec_order_out4[7] = 3;

            }

            //Corner 3
            if ((face_intecepted[0]) && (!face_intecepted[2]) && (!face_intecepted[1] && !face_intecepted[4] &&
                    face_intecepted[3] && face_intecepted[5])) {
                printf("Entrou no 3!");

                cut_edge[0] = 7;
                cut_edge[1] = 3;
                cut_edge[2] = 2;

                cut_edge_s[0] = 5;
                cut_edge_s[1] = 2;
                cut_edge_s[2] = 11;

                cut_edge_v[0] = 3;

                connec_order_in1[0] = 4;
                connec_order_in1[1] = 5;
                connec_order_in1[2] = 6;
                connec_order_in1[3] = 7;
                connec_order_in1[4] = 0;
                connec_order_in1[5] = 1;
                connec_order_in1[6] = 2;
                connec_order_in1[7] = 3;

                connec_order_out1[0] = 4;
                connec_order_out1[1] = 5;
                connec_order_out1[2] = 6;
                connec_order_out1[3] = 7;
                connec_order_out1[4] = 3;
                connec_order_out1[5] = 6;
                connec_order_out1[6] = 4;
                connec_order_out1[7] = 0;

                connec_order_in2[0] = 2;
                connec_order_in2[1] = 3;
                connec_order_in2[2] = 6;
                connec_order_in2[3] = 7;
                connec_order_in2[4] = 0;
                connec_order_in2[5] = 1;
                connec_order_in2[6] = 4;
                connec_order_in2[7] = 5;

                connec_order_out2[0] = 5;
                connec_order_out2[1] = 1;
                connec_order_out2[2] = 6;
                connec_order_out2[3] = 3;
                connec_order_out2[4] = 0;
                connec_order_out2[5] = 1;
                connec_order_out2[6] = 4;
                connec_order_out2[7] = 5;

                connec_order_in3[0] = 0;
                connec_order_in3[1] = 3;
                connec_order_in3[2] = 4;
                connec_order_in3[3] = 7;
                connec_order_in3[4] = 1;
                connec_order_in3[5] = 2;
                connec_order_in3[6] = 5;
                connec_order_in3[7] = 6;

                connec_order_out3[0] = 5;
                connec_order_out3[1] = 2;
                connec_order_out3[2] = 6;
                connec_order_out3[3] = 4;
                connec_order_out3[4] = 1;
                connec_order_out3[5] = 2;
                connec_order_out3[6] = 5;
                connec_order_out3[7] = 6;

                connec_order_in4[0] = 3;
                connec_order_in4[1] = 0;
                connec_order_in4[2] = 1;
                connec_order_in4[3] = 2;
                connec_order_in4[4] = 4;
                connec_order_in4[5] = 5;
                connec_order_in4[6] = 6;
                connec_order_in4[7] = 7;

                connec_order_out4[0] = 3;
                connec_order_out4[1] = 1;
                connec_order_out4[2] = 5;
                connec_order_out4[3] = 2;
                connec_order_out4[4] = 3;
                connec_order_out4[5] = 6;
                connec_order_out4[6] = 4;
                connec_order_out4[7] = 0;

            }


            //Corner 4
            if ((face_intecepted[0]) && (!face_intecepted[1]) && (!face_intecepted[3] && !face_intecepted[5] &&
                    face_intecepted[2] && face_intecepted[4])) {
                printf("Entrou no 4!");

                cut_edge[0] = 4;
                cut_edge[1] = 8;
                cut_edge[2] = 11;

                cut_edge_s[0] = 1;
                cut_edge_s[1] = 5;
                cut_edge_s[2] = 8;

                cut_edge_v[0] = 2;

                connec_order_in1[0] = 0;
                connec_order_in1[1] = 1;
                connec_order_in1[2] = 2;
                connec_order_in1[3] = 3;
                connec_order_in1[4] = 4;
                connec_order_in1[5] = 5;
                connec_order_in1[6] = 6;
                connec_order_in1[7] = 7;

                connec_order_out1[0] = 0;
                connec_order_out1[1] = 1;
                connec_order_out1[2] = 2;
                connec_order_out1[3] = 3;
                connec_order_out1[4] = 0;
                connec_order_out1[5] = 3;
                connec_order_out1[6] = 6;
                connec_order_out1[7] = 4;

                connec_order_in2[0] = 0;
                connec_order_in2[1] = 3;
                connec_order_in2[2] = 4;
                connec_order_in2[3] = 7;
                connec_order_in2[4] = 1;
                connec_order_in2[5] = 2;
                connec_order_in2[6] = 5;
                connec_order_in2[7] = 6;

                connec_order_out2[0] = 3;
                connec_order_out2[1] = 6;
                connec_order_out2[2] = 1;
                connec_order_out2[3] = 5;
                connec_order_out2[4] = 1;
                connec_order_out2[5] = 2;
                connec_order_out2[6] = 5;
                connec_order_out2[7] = 6;

                connec_order_in3[0] = 0;
                connec_order_in3[1] = 1;
                connec_order_in3[2] = 4;
                connec_order_in3[3] = 5;
                connec_order_in3[4] = 2;
                connec_order_in3[5] = 3;
                connec_order_in3[6] = 6;
                connec_order_in3[7] = 7;

                connec_order_out3[0] = 4;
                connec_order_out3[1] = 6;
                connec_order_out3[2] = 2;
                connec_order_out3[3] = 5;
                connec_order_out3[4] = 2;
                connec_order_out3[5] = 3;
                connec_order_out3[6] = 6;
                connec_order_out3[7] = 7;

                connec_order_in4[0] = 4;
                connec_order_in4[1] = 0;
                connec_order_in4[2] = 1;
                connec_order_in4[3] = 2;
                connec_order_in4[4] = 3;
                connec_order_in4[5] = 5;
                connec_order_in4[6] = 6;
                connec_order_in4[7] = 7;

                connec_order_out4[0] = 4;
                connec_order_out4[1] = 0;
                connec_order_out4[2] = 3;
                connec_order_out4[3] = 6;
                connec_order_out4[4] = 4;
                connec_order_out4[5] = 1;
                connec_order_out4[6] = 5;
                connec_order_out4[7] = 2;

            }

            //Corner 5
            if ((face_intecepted[1]) && (!face_intecepted[0]) && (!face_intecepted[3] && !face_intecepted[5] &&
                    face_intecepted[2] && face_intecepted[4])) {
                printf("Entrou no 5!");

                cut_edge[0] = 5;
                cut_edge[1] = 9;
                cut_edge[2] = 8;

                cut_edge_s[0] = 7;
                cut_edge_s[1] = 0;
                cut_edge_s[2] = 9;

                cut_edge_v[0] = 3;

                connec_order_in1[0] = 0;
                connec_order_in1[1] = 1;
                connec_order_in1[2] = 2;
                connec_order_in1[3] = 3;
                connec_order_in1[4] = 4;
                connec_order_in1[5] = 5;
                connec_order_in1[6] = 6;
                connec_order_in1[7] = 7;

                connec_order_out1[0] = 0;
                connec_order_out1[1] = 1;
                connec_order_out1[2] = 2;
                connec_order_out1[3] = 3;
                connec_order_out1[4] = 4;
                connec_order_out1[5] = 0;
                connec_order_out1[6] = 3;
                connec_order_out1[7] = 6;

                connec_order_in2[0] = 0;
                connec_order_in2[1] = 1;
                connec_order_in2[2] = 4;
                connec_order_in2[3] = 5;
                connec_order_in2[4] = 2;
                connec_order_in2[5] = 3;
                connec_order_in2[6] = 6;
                connec_order_in2[7] = 7;

                connec_order_out2[0] = 6;
                connec_order_out2[1] = 3;
                connec_order_out2[2] = 5;
                connec_order_out2[3] = 1;
                connec_order_out2[4] = 2;
                connec_order_out2[5] = 3;
                connec_order_out2[6] = 6;
                connec_order_out2[7] = 7;

                connec_order_in3[0] = 1;
                connec_order_in3[1] = 2;
                connec_order_in3[2] = 5;
                connec_order_in3[3] = 6;
                connec_order_in3[4] = 0;
                connec_order_in3[5] = 3;
                connec_order_in3[6] = 4;
                connec_order_in3[7] = 7;

                connec_order_out3[0] = 4;
                connec_order_out3[1] = 6;
                connec_order_out3[2] = 2;
                connec_order_out3[3] = 5;
                connec_order_out3[4] = 0;
                connec_order_out3[5] = 3;
                connec_order_out3[6] = 4;
                connec_order_out3[7] = 7;

                connec_order_in4[0] = 5;
                connec_order_in4[1] = 0;
                connec_order_in4[2] = 1;
                connec_order_in4[3] = 2;
                connec_order_in4[4] = 3;
                connec_order_in4[5] = 4;
                connec_order_in4[6] = 6;
                connec_order_in4[7] = 7;

                connec_order_out4[0] = 5;
                connec_order_out4[1] = 4;
                connec_order_out4[2] = 0;
                connec_order_out4[3] = 3;
                connec_order_out4[4] = 6;
                connec_order_out4[5] = 2;
                connec_order_out4[6] = 1;
                connec_order_out4[7] = 5;

            }

            //Corner 6
            if ((face_intecepted[1]) && (!face_intecepted[0]) && (!face_intecepted[2] && !face_intecepted[5] &&
                    face_intecepted[3] && face_intecepted[4])) {
                printf("Entrou no 6!");


                cut_edge[0] = 6;
                cut_edge[1] = 10;
                cut_edge[2] = 9;

                cut_edge_s[0] = 2;
                cut_edge_s[1] = 6;
                cut_edge_s[2] = 8;

                cut_edge_v[0] = 0;

                connec_order_in1[0] = 0;
                connec_order_in1[1] = 1;
                connec_order_in1[2] = 2;
                connec_order_in1[3] = 3;
                connec_order_in1[4] = 4;
                connec_order_in1[5] = 5;
                connec_order_in1[6] = 6;
                connec_order_in1[7] = 7;

                connec_order_out1[0] = 0;
                connec_order_out1[1] = 1;
                connec_order_out1[2] = 2;
                connec_order_out1[3] = 3;
                connec_order_out1[4] = 6;
                connec_order_out1[5] = 4;
                connec_order_out1[6] = 0;
                connec_order_out1[7] = 3;

                connec_order_in2[0] = 1;
                connec_order_in2[1] = 2;
                connec_order_in2[2] = 5;
                connec_order_in2[3] = 6;
                connec_order_in2[4] = 0;
                connec_order_in2[5] = 3;
                connec_order_in2[6] = 4;
                connec_order_in2[7] = 7;

                connec_order_out2[0] = 6;
                connec_order_out2[1] = 3;
                connec_order_out2[2] = 5;
                connec_order_out2[3] = 1;
                connec_order_out2[4] = 0;
                connec_order_out2[5] = 3;
                connec_order_out2[6] = 4;
                connec_order_out2[7] = 7;

                connec_order_in3[0] = 2;
                connec_order_in3[1] = 3;
                connec_order_in3[2] = 6;
                connec_order_in3[3] = 7;
                connec_order_in3[4] = 0;
                connec_order_in3[5] = 1;
                connec_order_in3[6] = 4;
                connec_order_in3[7] = 5;

                connec_order_out3[0] = 4;
                connec_order_out3[1] = 6;
                connec_order_out3[2] = 2;
                connec_order_out3[3] = 5;
                connec_order_out3[4] = 0;
                connec_order_out3[5] = 1;
                connec_order_out3[6] = 4;
                connec_order_out3[7] = 5;

                connec_order_in4[0] = 6;
                connec_order_in4[1] = 0;
                connec_order_in4[2] = 1;
                connec_order_in4[3] = 2;
                connec_order_in4[4] = 3;
                connec_order_in4[5] = 4;
                connec_order_in4[6] = 5;
                connec_order_in4[7] = 7;

                connec_order_out4[0] = 6;
                connec_order_out4[1] = 6;
                connec_order_out4[2] = 4;
                connec_order_out4[3] = 0;
                connec_order_out4[4] = 3;
                connec_order_out4[5] = 5;
                connec_order_out4[6] = 2;
                connec_order_out4[7] = 1;


            }

            //Corner 7
            if ((face_intecepted[0]) && (!face_intecepted[2]) && (!face_intecepted[1] && !face_intecepted[5] &&
                    face_intecepted[3] && face_intecepted[4])) {
                printf("Entrou no 7!");


                cut_edge[0] = 7;
                cut_edge[1] = 11;
                cut_edge[2] = 10;

                cut_edge_s[0] = 4;
                cut_edge_s[1] = 3;
                cut_edge_s[2] = 9;

                cut_edge_v[0] = 1;

                connec_order_in1[0] = 0;
                connec_order_in1[1] = 1;
                connec_order_in1[2] = 2;
                connec_order_in1[3] = 3;
                connec_order_in1[4] = 4;
                connec_order_in1[5] = 5;
                connec_order_in1[6] = 6;
                connec_order_in1[7] = 7;

                connec_order_out1[0] = 0;
                connec_order_out1[1] = 1;
                connec_order_out1[2] = 2;
                connec_order_out1[3] = 3;
                connec_order_out1[4] = 3;
                connec_order_out1[5] = 6;
                connec_order_out1[6] = 4;
                connec_order_out1[7] = 0;

                connec_order_in2[0] = 2;
                connec_order_in2[1] = 3;
                connec_order_in2[2] = 6;
                connec_order_in2[3] = 7;
                connec_order_in2[4] = 0;
                connec_order_in2[5] = 1;
                connec_order_in2[6] = 4;
                connec_order_in2[7] = 5;

                connec_order_out2[0] = 6;
                connec_order_out2[1] = 3;
                connec_order_out2[2] = 5;
                connec_order_out2[3] = 1;
                connec_order_out2[4] = 0;
                connec_order_out2[5] = 1;
                connec_order_out2[6] = 4;
                connec_order_out2[7] = 5;

                connec_order_in3[0] = 0;
                connec_order_in3[1] = 3;
                connec_order_in3[2] = 4;
                connec_order_in3[3] = 7;
                connec_order_in3[4] = 1;
                connec_order_in3[5] = 2;
                connec_order_in3[6] = 5;
                connec_order_in3[7] = 6;

                connec_order_out3[0] = 6;
                connec_order_out3[1] = 4;
                connec_order_out3[2] = 5;
                connec_order_out3[3] = 2;
                connec_order_out3[4] = 1;
                connec_order_out3[5] = 2;
                connec_order_out3[6] = 5;
                connec_order_out3[7] = 6;

                connec_order_in4[0] = 7;
                connec_order_in4[1] = 0;
                connec_order_in4[2] = 1;
                connec_order_in4[3] = 2;
                connec_order_in4[4] = 3;
                connec_order_in4[5] = 4;
                connec_order_in4[6] = 5;
                connec_order_in4[7] = 6;

                connec_order_out4[0] = 7;
                connec_order_out4[1] = 3;
                connec_order_out4[2] = 6;
                connec_order_out4[3] = 4;
                connec_order_out4[4] = 0;
                connec_order_out4[5] = 1;
                connec_order_out4[6] = 5;
                connec_order_out4[7] = 2;

            }

            GtsPoint *p0 = point[cut_edge[0]];
            GtsPoint *p1 = point[cut_edge[1]];
            GtsPoint *p2 = point[cut_edge[2]];

            GtsPoint *p3 = point_s[cut_edge_s[0]];
            GtsPoint *p4 = point_s[cut_edge_s[1]];
            GtsPoint *p5 = point_s[cut_edge_s[2]];

            GtsPoint *p6 = point_v[cut_edge_v[0]];

            g_assert(p0 != NULL);
            g_assert(p1 != NULL);
            g_assert(p2 != NULL);
            g_assert(p3 != NULL);
            g_assert(p4 != NULL);
            g_assert(p5 != NULL);
            g_assert(p6 != NULL);

            conn_t2[0] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p0, coords);
            conn_t2[1] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p1, coords);
            conn_t2[2] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p2, coords);

            // add 3 extra points in the surface
            conn_t2[3] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[0]], hash_nodes, mesh->local_n_nodes, p3, coords);
            conn_t2[4] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[1]], hash_nodes, mesh->local_n_nodes, p4, coords);
            conn_t2[5] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[2]], hash_nodes, mesh->local_n_nodes, p5, coords);

            // add 1 extra points in the volume
            conn_t2[6] = AddPointOnEdge(Edge2GNode_s[cut_edge_v[0]], hash_nodes, mesh->local_n_nodes, p6, coords);

            octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

            elem1->nodes[connec_order_in1[0]].id = original_conn[connec_order_out1[0]];
            elem1->nodes[connec_order_in1[1]].id = original_conn[connec_order_out1[1]];
            elem1->nodes[connec_order_in1[2]].id = original_conn[connec_order_out1[2]];
            elem1->nodes[connec_order_in1[3]].id = original_conn[connec_order_out1[3]];

            elem1->nodes[connec_order_in1[4]].id = conn_t2[connec_order_out1[4]];
            elem1->nodes[connec_order_in1[5]].id = conn_t2[connec_order_out1[5]];
            elem1->nodes[connec_order_in1[6]].id = conn_t2[connec_order_out1[6]];
            elem1->nodes[connec_order_in1[7]].id = conn_t2[connec_order_out1[7]];

            octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

            elem2->nodes[connec_order_in2[0]].id = conn_t2[connec_order_out2[0]];
            elem2->nodes[connec_order_in2[1]].id = conn_t2[connec_order_out2[1]];
            elem2->nodes[connec_order_in2[2]].id = conn_t2[connec_order_out2[2]];
            elem2->nodes[connec_order_in2[3]].id = conn_t2[connec_order_out2[3]];

            elem2->nodes[connec_order_in2[4]].id = original_conn[connec_order_out2[4]];
            elem2->nodes[connec_order_in2[5]].id = original_conn[connec_order_out2[5]];
            elem2->nodes[connec_order_in2[6]].id = original_conn[connec_order_out2[6]];
            elem2->nodes[connec_order_in2[7]].id = original_conn[connec_order_out2[7]];

            elem2->pad = 6;
            elem2->level = elem->level;

            octant_t* elem3 = (octant_t*) sc_array_push(&mesh->elements);

            elem3->nodes[connec_order_in3[0]].id = conn_t2[connec_order_out3[0]];
            elem3->nodes[connec_order_in3[1]].id = conn_t2[connec_order_out3[1]];
            elem3->nodes[connec_order_in3[2]].id = conn_t2[connec_order_out3[2]];
            elem3->nodes[connec_order_in3[3]].id = conn_t2[connec_order_out3[3]];

            elem3->nodes[connec_order_in3[4]].id = original_conn[connec_order_out3[4]];
            elem3->nodes[connec_order_in3[5]].id = original_conn[connec_order_out3[5]];
            elem3->nodes[connec_order_in3[6]].id = original_conn[connec_order_out3[6]];
            elem3->nodes[connec_order_in3[7]].id = original_conn[connec_order_out3[7]];

            elem3->pad = 7;
            elem3->level = elem->level;

            octant_t* elem4 = (octant_t*) sc_array_push(&mesh->elements);

            elem4->nodes[connec_order_in4[0]].id = original_conn[connec_order_out4[0]];

            elem4->nodes[connec_order_in4[1]].id = conn_t2[connec_order_out4[1]];
            elem4->nodes[connec_order_in4[2]].id = conn_t2[connec_order_out4[2]];
            elem4->nodes[connec_order_in4[3]].id = conn_t2[connec_order_out4[3]];
            elem4->nodes[connec_order_in4[4]].id = conn_t2[connec_order_out4[4]];
            elem4->nodes[connec_order_in4[5]].id = conn_t2[connec_order_out4[5]];
            elem4->nodes[connec_order_in4[6]].id = conn_t2[connec_order_out4[6]];
            elem4->nodes[connec_order_in4[7]].id = conn_t2[connec_order_out4[7]];

            elem4->pad = 8;
            elem4->level = elem->level;


#endif




            for (int edge = 0; edge < 4; edge++) {
                if (point_v[edge]) gts_object_destroy(GTS_OBJECT(point_v[edge]));
                point_v[edge] = NULL;
            }

            for (int edge = 0; edge < 12; edge++) {
                if (point_s[edge]) gts_object_destroy(GTS_OBJECT(point_s[edge]));
                point_s[edge] = NULL;
            }

            for (int edge = 0; edge < 12; edge++) {
                if (point[edge]) gts_object_destroy(GTS_OBJECT(point[edge]));
                point[edge] = NULL;
            }

        }

    //if (elem->pad == -1) {
    //    elem->pad = ed_cont;
    //}
        
    }

    mesh->local_n_elements = mesh->elements.elem_count;
    //mesh->local_n_nodes = coords.size()/3;

    //MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    if (mesh->mpi_rank == 0) {
        printf("Total number of elements: %lld\n", mesh->local_n_elements);
        printf("Total number of nodes: %lld\n", mesh->local_n_nodes);
    }


    fclose(fdbg);
    sc_hash_array_destroy(hash_nodes);


}

/*
 *
 */
int AddPointOnEdge(int* nodes, sc_hash_array_t* hash, int &npoints, GtsPoint *p, std::vector<double> &coords) {
    size_t position;
    node_in_edge_t *r;
    node_in_edge_t key;
    key.coord[0] = nodes[0];
    key.coord[1] = nodes[1];

    r = (node_in_edge_t*) sc_hash_array_insert_unique(hash, &key, &position);
    if (r != NULL) {
        r->coord[0] = key.coord[0];
        r->coord[1] = key.coord[1];
        r->node_id = npoints;
        npoints++;
        coords.push_back(p->x);
        coords.push_back(p->y);
        coords.push_back(p->z);
        return r->node_id;
    } else {
        r = (node_in_edge_t*) sc_array_index(&hash->a, position);
        return r->node_id;
    }

}

void shift_n_position_from_index(sc_array_t *elements, uint32_t n, uint32_t index) {
    uint32_t incount = elements->elem_count;

    sc_array_resize(elements, incount + n);

    for (int i = incount - 1; i < index; i--) {
        octant_t *new_pos = (octant_t*) sc_array_index(elements, i + n);
        octant_t *old_pos = (octant_t*) sc_array_index(elements, i);
        new_pos->level = old_pos->level;
        new_pos->pad = old_pos->pad;
        new_pos->x = old_pos->x;
        new_pos->y = old_pos->y;
        new_pos->z = old_pos->z;

        //new_pos->nodes->id = old_pos->nodes->id

        for (int i = 0; i < 8; i++) {
            new_pos->nodes[i].id = old_pos->nodes[i].id;
            new_pos->nodes[i].x = old_pos->nodes[i].x;
            new_pos->nodes[i].y = old_pos->nodes[i].y;
            new_pos->nodes[i].z = old_pos->nodes[i].z;
        }

        //memcpy(new_pos->nodes, old_pos->nodes, 8 * sizeof (octant_node_t));
    }
}
