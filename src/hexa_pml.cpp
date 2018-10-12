
#include <gts.h>
#include <gts.h>
#include <glib.h>
#include <vector>
#include <iostream>
using namespace std;
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>

#include "hexa.h"
#include "pml.h"
#include "hilbert.h"

/*
typedef struct {
	bitmask_t coord[2];
	int node_id;
} node_in_edge_t;

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
*/
int8_t SetNodePML(hexa_tree_t* tree, octant_node_t* node) {
	int8_t pml_id = 0;
	if (node->x == 0) pml_id |= PML_X0;
	if (node->x == tree->ncellx) pml_id |= PML_X1;
	if (node->y == 0) pml_id |= PML_Y0;
	if (node->y == tree->ncelly) pml_id |= PML_Y1;
	//if(node->z == 0)               pml_id |= PML_Z0;
	if (node->z == tree->max_z) pml_id |= PML_Z1;
	return pml_id;
}

inline int isX0(int8_t pml_id) {
	return ((pml_id & PML_X0) == PML_X0);
}

inline int isX1(int8_t pml_id) {
	return ((pml_id & PML_X1) == PML_X1);
}

inline int isY0(int8_t pml_id) {
	return ((pml_id & PML_Y0) == PML_Y0);
}

inline int isY1(int8_t pml_id) {
	return ((pml_id & PML_Y1) == PML_Y1);
}

inline int isZ0(int8_t pml_id) {
	return ((pml_id & PML_Z0) == PML_Z0);
}

inline int isZ1(int8_t pml_id) {
	return ((pml_id & PML_Z1) == PML_Z1);
}

void SetElemPML(hexa_tree_t* tree, octant_t *elem) {

	int8_t pml_id[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	for (int i = 0; i < 8; ++i) {
		pml_id[i] = SetNodePML(tree, &elem->nodes[i]);
	}

	for (int face = 0; face < 6; ++face) {
		int no1 = face_map[face][0];
		int no2 = face_map[face][1];
		int no3 = face_map[face][2];
		int no4 = face_map[face][3];

		if (isX0(pml_id[no1]) && isX0(pml_id[no2]) && isX0(pml_id[no3]) && isX0(pml_id[no4]))
			elem->pml_id |= PML_X0;

		if (isX1(pml_id[no1]) && isX1(pml_id[no2]) && isX1(pml_id[no3]) && isX1(pml_id[no4]))
			elem->pml_id |= PML_X1;

		if (isY0(pml_id[no1]) && isY0(pml_id[no2]) && isY0(pml_id[no3]) && isY0(pml_id[no4]))
			elem->pml_id |= PML_Y0;

		if (isY1(pml_id[no1]) && isY1(pml_id[no2]) && isY1(pml_id[no3]) && isY1(pml_id[no4]))
			elem->pml_id |= PML_Y1;

		if (isZ0(pml_id[no1]) && isZ0(pml_id[no2]) && isZ0(pml_id[no3]) && isZ0(pml_id[no4]))
			elem->pml_id |= PML_Z0;

		if (isZ1(pml_id[no1]) && isZ1(pml_id[no2]) && isZ1(pml_id[no3]) && isZ1(pml_id[no4]))
			elem->pml_id |= PML_Z1;

	}

}

/*
inline void SetPMLMask(int8_t* mask, int8_t pml_id)
{

    mask[PML_CORNER_X0Y0Z0] = (((pml_id&PML_X0)==PML_X0) && ((pml_id&PML_Y0)==PML_Y0) && ((pml_id&PML_Z0)==PML_Z0));
    mask[PML_CORNER_X1Y0Z0] = (((pml_id&PML_X1)==PML_X1) && ((pml_id&PML_Y0)==PML_Y0) && ((pml_id&PML_Z0)==PML_Z0));
    mask[PML_CORNER_X0Y1Z0] = (((pml_id&PML_X0)==PML_X0) && ((pml_id&PML_Y1)==PML_Y1) && ((pml_id&PML_Z0)==PML_Z0));
    mask[PML_CORNER_X1Y1Z0] = (((pml_id&PML_X1)==PML_X1) && ((pml_id&PML_Y1)==PML_Y1) && ((pml_id&PML_Z0)==PML_Z0));
    mask[PML_CORNER_X0Y0Z1] = (((pml_id&PML_X0)==PML_X0) && ((pml_id&PML_Y0)==PML_Y0) && ((pml_id&PML_Z1)==PML_Z1));
    mask[PML_CORNER_X1Y0Z1] = (((pml_id&PML_X1)==PML_X1) && ((pml_id&PML_Y0)==PML_Y0) && ((pml_id&PML_Z1)==PML_Z1));
    mask[PML_CORNER_X0Y1Z1] = (((pml_id&PML_X0)==PML_X0) && ((pml_id&PML_Y1)==PML_Y1) && ((pml_id&PML_Z1)==PML_Z1));
    mask[PML_CORNER_X1Y1Z1] = (((pml_id&PML_X1)==PML_X1) && ((pml_id&PML_Y1)==PML_Y1) && ((pml_id&PML_Z1)==PML_Z1));

    mask[PML_EDGE_Z0_X0] = ((pml_id&PML_X0==PML_X0) && (pml_id&PML_Z0==PML_Z0));
    mask[PML_EDGE_Z0_X1] = ((pml_id&PML_X1==PML_X1) && (pml_id&PML_Z0==PML_Z0));
    mask[PML_EDGE_Z0_Y0] = ((pml_id&PML_Y0==PML_Y0) && (pml_id&PML_Z0==PML_Z0));
    mask[PML_EDGE_Z0_Y1] = ((pml_id&PML_Y1==PML_Y1) && (pml_id&PML_Z0==PML_Z0));

    mask[PML_EDGE_X0_Y0] = ((pml_id&PML_X0==PML_X0) && (pml_id&PML_Y0==PML_Y0));
    mask[PML_EDGE_X0_Y1] = ((pml_id&PML_X0==PML_X1) && (pml_id&PML_Y1==PML_Y1));
    mask[PML_EDGE_X1_Y0] = ((pml_id&PML_X1==PML_X1) && (pml_id&PML_Y0==PML_Y0));
    mask[PML_EDGE_X1_Y1] = ((pml_id&PML_X1==PML_X1) && (pml_id&PML_Y1==PML_Y1));

    mask[PML_EDGE_Z1_X0] = ((pml_id&PML_X0==PML_X0) && (pml_id&PML_Z1==PML_Z1));
    mask[PML_EDGE_Z1_X1] = ((pml_id&PML_X1==PML_X1) && (pml_id&PML_Z1==PML_Z1));
    mask[PML_EDGE_Z1_Y0] = ((pml_id&PML_Y0==PML_Y0) && (pml_id&PML_Z1==PML_Z1));
    mask[PML_EDGE_Z1_Y1] = ((pml_id&PML_Y1==PML_Y1) && (pml_id&PML_Z1==PML_Z1));

    mask[PML_FACE_X0] = (pml_id&PML_X0==PML_X0);
    mask[PML_FACE_X1] = (pml_id&PML_X1==PML_X1);
    mask[PML_FACE_Y0] = (pml_id&PML_Y0==PML_Y0);
    mask[PML_FACE_Y1] = (pml_id&PML_Y1==PML_Y1);
    mask[PML_FACE_Z0] = (pml_id&PML_Z0==PML_Z0);
    mask[PML_FACE_Z1] = (pml_id&PML_Z1==PML_Z1);

}
 */

inline void SetPMLMask(int8_t* mask, int8_t pml_id) {

	mask[PML_CORNER_X0Y0Z0] = (isX0(pml_id) && isY0(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X1Y0Z0] = (isX1(pml_id) && isY0(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X0Y1Z0] = (isX0(pml_id) && isY1(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X1Y1Z0] = (isX1(pml_id) && isY1(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X0Y0Z1] = (isX0(pml_id) && isY0(pml_id) && isZ1(pml_id));
	mask[PML_CORNER_X1Y0Z1] = (isX1(pml_id) && isY0(pml_id) && isZ1(pml_id));
	mask[PML_CORNER_X0Y1Z1] = (isX0(pml_id) && isY1(pml_id) && isZ1(pml_id));
	mask[PML_CORNER_X1Y1Z1] = (isX1(pml_id) && isY1(pml_id) && isZ1(pml_id));

	mask[PML_EDGE_Z0_X0] = (isX0(pml_id) && isZ0(pml_id));
	mask[PML_EDGE_Z0_X1] = (isX1(pml_id) && isZ0(pml_id));
	mask[PML_EDGE_Z0_Y0] = (isY0(pml_id) && isZ0(pml_id));
	mask[PML_EDGE_Z0_Y1] = (isY1(pml_id) && isZ0(pml_id));

	mask[PML_EDGE_X0_Y0] = (isX0(pml_id) && isY0(pml_id));
	mask[PML_EDGE_X0_Y1] = (isX0(pml_id) && isY1(pml_id));
	mask[PML_EDGE_X1_Y0] = (isX1(pml_id) && isY0(pml_id));
	mask[PML_EDGE_X1_Y1] = (isX1(pml_id) && isY1(pml_id));

	mask[PML_EDGE_Z1_X0] = (isX0(pml_id) && isZ1(pml_id));
	mask[PML_EDGE_Z1_X1] = (isX1(pml_id) && isZ1(pml_id));
	mask[PML_EDGE_Z1_Y0] = (isY0(pml_id) && isZ1(pml_id));
	mask[PML_EDGE_Z1_Y1] = (isY1(pml_id) && isZ1(pml_id));

	mask[PML_FACE_X0] = isX0(pml_id);
	mask[PML_FACE_X1] = isX1(pml_id);
	mask[PML_FACE_Y0] = isY0(pml_id);
	mask[PML_FACE_Y1] = isY1(pml_id);
	mask[PML_FACE_Z0] = isZ0(pml_id);
	mask[PML_FACE_Z1] = isZ1(pml_id);
}

inline void SetPMLMask_corner(int8_t* mask, int8_t pml_id) {

	mask[PML_CORNER_X0Y0Z0] = (isX0(pml_id) && isY0(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X1Y0Z0] = (isX1(pml_id) && isY0(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X0Y1Z0] = (isX0(pml_id) && isY1(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X1Y1Z0] = (isX1(pml_id) && isY1(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X0Y0Z1] = (isX0(pml_id) && isY0(pml_id) && isZ1(pml_id));
	mask[PML_CORNER_X1Y0Z1] = (isX1(pml_id) && isY0(pml_id) && isZ1(pml_id));
	mask[PML_CORNER_X0Y1Z1] = (isX0(pml_id) && isY1(pml_id) && isZ1(pml_id));
	mask[PML_CORNER_X1Y1Z1] = (isX1(pml_id) && isY1(pml_id) && isZ1(pml_id));

}

inline void SetPMLMask_edge(int8_t* mask, int8_t pml_id) {

	mask[PML_EDGE_Z0_X0] = (isX0(pml_id) && isZ0(pml_id));
	mask[PML_EDGE_Z0_X1] = (isX1(pml_id) && isZ0(pml_id));
	mask[PML_EDGE_Z0_Y0] = (isY0(pml_id) && isZ0(pml_id));
	mask[PML_EDGE_Z0_Y1] = (isY1(pml_id) && isZ0(pml_id));

	mask[PML_EDGE_X0_Y0] = (isX0(pml_id) && isY0(pml_id));
	mask[PML_EDGE_X0_Y1] = (isX0(pml_id) && isY1(pml_id));
	mask[PML_EDGE_X1_Y0] = (isX1(pml_id) && isY0(pml_id));
	mask[PML_EDGE_X1_Y1] = (isX1(pml_id) && isY1(pml_id));

	mask[PML_EDGE_Z1_X0] = (isX0(pml_id) && isZ1(pml_id));
	mask[PML_EDGE_Z1_X1] = (isX1(pml_id) && isZ1(pml_id));
	mask[PML_EDGE_Z1_Y0] = (isY0(pml_id) && isZ1(pml_id));
	mask[PML_EDGE_Z1_Y1] = (isY1(pml_id) && isZ1(pml_id));

}

inline void SetPMLMask_face(int8_t* mask, int8_t pml_id) {

	mask[PML_FACE_X0] = isX0(pml_id);
	mask[PML_FACE_X1] = isX1(pml_id);
	mask[PML_FACE_Y0] = isY0(pml_id);
	mask[PML_FACE_Y1] = isY1(pml_id);
	mask[PML_FACE_Z0] = isZ0(pml_id);
	mask[PML_FACE_Z1] = isZ1(pml_id);
}

void AddPMLElements(hexa_tree_t* mesh) {

	int8_t mask[NPML];
	int32_t npmls[NPML] = {0};

	int32_t X_pml = 1;
	int32_t Y_pml = 1;
	int32_t Z_pml = 1;
	mesh->total_n_mat = 18;
	sc_array_t *elements = &mesh->elements;
	int nelem = elements->elem_count;
	for (int i = 0; i < nelem; ++i) {
		octant_t* elem = (octant_t*) sc_array_index(elements, i);
		elem->pml_id = 0;

		SetElemPML(mesh, elem);
		SetPMLMask(mask, elem->pml_id);
		elem->n_mat = 18;

		if (mask[PML_CORNER_X0Y0Z0]) {

			npmls[PML_CORNER_X0Y0Z0]++;

			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_CORNER_X0Y0Z0;
			pml_e->n_mat = 50;

			pml_e->nodes[6] = elem->nodes[0];

			pml_e->nodes[0].x = elem->nodes[0].x - X_pml;
			pml_e->nodes[0].y = elem->nodes[0].y - Y_pml;
			pml_e->nodes[0].z = elem->nodes[0].z - Z_pml;

			pml_e->nodes[1].x = elem->nodes[0].x;
			pml_e->nodes[1].y = elem->nodes[0].y - Y_pml;
			pml_e->nodes[1].z = elem->nodes[0].z - Z_pml;

			pml_e->nodes[2].x = elem->nodes[0].x;
			pml_e->nodes[2].y = elem->nodes[0].y;
			pml_e->nodes[2].z = elem->nodes[0].z - Z_pml;


			pml_e->nodes[3].x = elem->nodes[0].x - X_pml;
			pml_e->nodes[3].y = elem->nodes[0].y;
			pml_e->nodes[3].z = elem->nodes[0].z - Z_pml;

			pml_e->nodes[4].x = elem->nodes[0].x - X_pml;
			pml_e->nodes[4].y = elem->nodes[0].y - Y_pml;
			pml_e->nodes[4].z = elem->nodes[0].z;

			pml_e->nodes[5].x = elem->nodes[0].x;
			pml_e->nodes[5].y = elem->nodes[0].y - Y_pml;
			pml_e->nodes[5].z = elem->nodes[0].z;

			pml_e->nodes[7].x = elem->nodes[0].x - X_pml;
			pml_e->nodes[7].y = elem->nodes[0].y;
			pml_e->nodes[7].z = elem->nodes[0].z;

			pml_e->x = elem->x-1;
			pml_e->y = elem->y-1;
			pml_e->z = elem->z-1;

			//continue;
		}

		if (mask[PML_CORNER_X1Y0Z0]) {

			npmls[PML_CORNER_X1Y0Z0]++;

			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_CORNER_X1Y0Z0;
			elem->n_mat = 50;

			pml_e->nodes[7] = elem->nodes[1];

			pml_e->nodes[0].x = elem->nodes[1].x;
			pml_e->nodes[0].y = elem->nodes[1].y - Y_pml;
			pml_e->nodes[0].z = elem->nodes[1].z - Z_pml;

			pml_e->nodes[1].x = elem->nodes[1].x + X_pml;
			pml_e->nodes[1].y = elem->nodes[1].y - Y_pml;
			pml_e->nodes[1].z = elem->nodes[1].z - Z_pml;

			pml_e->nodes[2].x = elem->nodes[1].x + X_pml;
			pml_e->nodes[2].y = elem->nodes[1].y;
			pml_e->nodes[2].z = elem->nodes[1].z - Z_pml;


			pml_e->nodes[3].x = elem->nodes[1].x;
			pml_e->nodes[3].y = elem->nodes[1].y;
			pml_e->nodes[3].z = elem->nodes[1].z - Z_pml;

			pml_e->nodes[4].x = elem->nodes[1].x;
			pml_e->nodes[4].y = elem->nodes[1].y - Y_pml;
			pml_e->nodes[4].z = elem->nodes[1].z;

			pml_e->nodes[5].x = elem->nodes[1].x + X_pml;
			pml_e->nodes[5].y = elem->nodes[1].y - Y_pml;
			pml_e->nodes[5].z = elem->nodes[1].z;

			pml_e->nodes[6].x = elem->nodes[1].x + X_pml;
			pml_e->nodes[6].y = elem->nodes[1].y;
			pml_e->nodes[6].z = elem->nodes[1].z;

			pml_e->x = elem->x+1;
			pml_e->y = elem->y-1;
			pml_e->z = elem->z-1;

			//continue;
		}

		if (mask[PML_CORNER_X0Y1Z0]) {
			npmls[PML_CORNER_X0Y1Z0]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_CORNER_X0Y1Z0;
			pml_e->n_mat = 50;

			pml_e->nodes[5] = elem->nodes[3];

			pml_e->nodes[0].x = elem->nodes[3].x - X_pml;
			pml_e->nodes[0].y = elem->nodes[3].y;
			pml_e->nodes[0].z = elem->nodes[3].z - Z_pml;

			pml_e->nodes[1].x = elem->nodes[3].x;
			pml_e->nodes[1].y = elem->nodes[3].y;
			pml_e->nodes[1].z = elem->nodes[3].z - Z_pml;

			pml_e->nodes[2].x = elem->nodes[3].x;
			pml_e->nodes[2].y = elem->nodes[3].y + Y_pml;
			pml_e->nodes[2].z = elem->nodes[3].z - Z_pml;


			pml_e->nodes[3].x = elem->nodes[3].x - X_pml;
			pml_e->nodes[3].y = elem->nodes[3].y + Y_pml;
			pml_e->nodes[3].z = elem->nodes[3].z - Z_pml;

			pml_e->nodes[4].x = elem->nodes[3].x - X_pml;
			pml_e->nodes[4].y = elem->nodes[3].y;
			pml_e->nodes[4].z = elem->nodes[3].z;

			pml_e->nodes[6].x = elem->nodes[3].x;
			pml_e->nodes[6].y = elem->nodes[3].y + Y_pml;
			pml_e->nodes[6].z = elem->nodes[3].z;

			pml_e->nodes[7].x = elem->nodes[3].x - X_pml;
			pml_e->nodes[7].y = elem->nodes[3].y + Y_pml;
			pml_e->nodes[7].z = elem->nodes[3].z;

			pml_e->x = elem->x-1;
			pml_e->y = elem->y+1;
			pml_e->z = elem->z-1;

			//continue;
		}

		if (mask[PML_CORNER_X1Y1Z0]) {
			npmls[PML_CORNER_X1Y1Z0]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_CORNER_X1Y1Z0;
			pml_e->n_mat = 50;

			pml_e->nodes[4] = elem->nodes[2];

			pml_e->nodes[0].x = elem->nodes[2].x;
			pml_e->nodes[0].y = elem->nodes[2].y;
			pml_e->nodes[0].z = elem->nodes[2].z - Z_pml;

			pml_e->nodes[1].x = elem->nodes[2].x + X_pml;
			pml_e->nodes[1].y = elem->nodes[2].y;
			pml_e->nodes[1].z = elem->nodes[2].z - Z_pml;

			pml_e->nodes[2].x = elem->nodes[2].x + X_pml;
			pml_e->nodes[2].y = elem->nodes[2].y + Y_pml;
			pml_e->nodes[2].z = elem->nodes[2].z - Z_pml;


			pml_e->nodes[3].x = elem->nodes[2].x;
			pml_e->nodes[3].y = elem->nodes[2].y + Y_pml;
			pml_e->nodes[3].z = elem->nodes[2].z - Z_pml;

			pml_e->nodes[5].x = elem->nodes[2].x + X_pml;
			pml_e->nodes[5].y = elem->nodes[2].y;
			pml_e->nodes[5].z = elem->nodes[2].z;

			pml_e->nodes[6].x = elem->nodes[2].x + X_pml;
			pml_e->nodes[6].y = elem->nodes[2].y + Y_pml;
			pml_e->nodes[6].z = elem->nodes[2].z;

			pml_e->nodes[7].x = elem->nodes[2].x;
			pml_e->nodes[7].y = elem->nodes[2].y + Y_pml;
			pml_e->nodes[7].z = elem->nodes[2].z;

			pml_e->x = elem->x+1;
			pml_e->y = elem->y+1;
			pml_e->z = elem->z-1;

			//continue;
		}

		if (mask[PML_CORNER_X0Y0Z1]) {
			npmls[PML_CORNER_X0Y0Z1]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_CORNER_X0Y0Z1;
			pml_e->n_mat = 1;

			pml_e->nodes[2] = elem->nodes[4];

			pml_e->nodes[0].x = elem->nodes[4].x - X_pml;
			pml_e->nodes[0].y = elem->nodes[4].y - Y_pml;
			pml_e->nodes[0].z = elem->nodes[4].z;

			pml_e->nodes[1].x = elem->nodes[4].x;
			pml_e->nodes[1].y = elem->nodes[4].y - Y_pml;
			pml_e->nodes[1].z = elem->nodes[4].z;

			pml_e->nodes[3].x = elem->nodes[4].x - X_pml;
			pml_e->nodes[3].y = elem->nodes[4].y;
			pml_e->nodes[3].z = elem->nodes[4].z;


			pml_e->nodes[4].x = elem->nodes[4].x - Y_pml;
			pml_e->nodes[4].y = elem->nodes[4].y - Y_pml;
			pml_e->nodes[4].z = elem->nodes[4].z + Z_pml;

			pml_e->nodes[5].x = elem->nodes[4].x;
			pml_e->nodes[5].y = elem->nodes[4].y - Y_pml;
			pml_e->nodes[5].z = elem->nodes[4].z + Z_pml;

			pml_e->nodes[6].x = elem->nodes[4].x;
			pml_e->nodes[6].y = elem->nodes[4].y;
			pml_e->nodes[6].z = elem->nodes[4].z + Z_pml;

			pml_e->nodes[7].x = elem->nodes[4].x - X_pml;
			pml_e->nodes[7].y = elem->nodes[4].y;
			pml_e->nodes[7].z = elem->nodes[4].z + Z_pml;

			pml_e->x = elem->x-1;
			pml_e->y = elem->y-1;
			pml_e->z = elem->z+1;

			//continue;
		}

		if (mask[PML_CORNER_X1Y0Z1]) {
			npmls[PML_CORNER_X1Y0Z1]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_CORNER_X1Y0Z1;
			pml_e->n_mat = 2;

			pml_e->nodes[3] = elem->nodes[5];

			pml_e->nodes[0].x = elem->nodes[5].x;
			pml_e->nodes[0].y = elem->nodes[5].y - Y_pml;
			pml_e->nodes[0].z = elem->nodes[5].z;

			pml_e->nodes[1].x = elem->nodes[5].x + X_pml;
			pml_e->nodes[1].y = elem->nodes[5].y - Y_pml;
			pml_e->nodes[1].z = elem->nodes[5].z;

			pml_e->nodes[2].x = elem->nodes[5].x + X_pml;
			pml_e->nodes[2].y = elem->nodes[5].y;
			pml_e->nodes[2].z = elem->nodes[5].z;


			pml_e->nodes[4].x = elem->nodes[5].x;
			pml_e->nodes[4].y = elem->nodes[5].y - Y_pml;
			pml_e->nodes[4].z = elem->nodes[5].z + Z_pml;

			pml_e->nodes[5].x = elem->nodes[5].x + X_pml;
			pml_e->nodes[5].y = elem->nodes[5].y - Y_pml;
			pml_e->nodes[5].z = elem->nodes[5].z + Z_pml;

			pml_e->nodes[6].x = elem->nodes[5].x + X_pml;
			pml_e->nodes[6].y = elem->nodes[5].y;
			pml_e->nodes[6].z = elem->nodes[5].z + Z_pml;

			pml_e->nodes[7].x = elem->nodes[5].x;
			pml_e->nodes[7].y = elem->nodes[5].y;
			pml_e->nodes[7].z = elem->nodes[5].z + Z_pml;

			pml_e->x = elem->x+1;
			pml_e->y = elem->y-1;
			pml_e->z = elem->z+1;

			//continue;
		}

		if (mask[PML_CORNER_X0Y1Z1]) {
			npmls[PML_CORNER_X0Y1Z1]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_CORNER_X0Y1Z1;
			pml_e->n_mat = 3;

			pml_e->nodes[1] = elem->nodes[7];

			pml_e->nodes[0].x = elem->nodes[7].x - X_pml;
			pml_e->nodes[0].y = elem->nodes[7].y;
			pml_e->nodes[0].z = elem->nodes[7].z;

			pml_e->nodes[2].x = elem->nodes[7].x;
			pml_e->nodes[2].y = elem->nodes[7].y + Y_pml;
			pml_e->nodes[2].z = elem->nodes[7].z;

			pml_e->nodes[3].x = elem->nodes[7].x - X_pml;
			pml_e->nodes[3].y = elem->nodes[7].y + Y_pml;
			pml_e->nodes[3].z = elem->nodes[7].z;


			pml_e->nodes[4].x = elem->nodes[7].x - X_pml;
			pml_e->nodes[4].y = elem->nodes[7].y;
			pml_e->nodes[4].z = elem->nodes[7].z + Z_pml;

			pml_e->nodes[5].x = elem->nodes[7].x;
			pml_e->nodes[5].y = elem->nodes[7].y;
			pml_e->nodes[5].z = elem->nodes[7].z + Z_pml;

			pml_e->nodes[6].x = elem->nodes[7].x;
			pml_e->nodes[6].y = elem->nodes[7].y + Y_pml;
			pml_e->nodes[6].z = elem->nodes[7].z + Z_pml;

			pml_e->nodes[7].x = elem->nodes[7].x - X_pml;
			pml_e->nodes[7].y = elem->nodes[7].y + Y_pml;
			pml_e->nodes[7].z = elem->nodes[7].z + Z_pml;

			pml_e->x = elem->x-1;
			pml_e->y = elem->y+1;
			pml_e->z = elem->z+1;

			//continue;
		}

		if (mask[PML_CORNER_X1Y1Z1]) {
			npmls[PML_CORNER_X1Y1Z1]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_CORNER_X1Y1Z1;
			pml_e->n_mat = 4;

			pml_e->nodes[0] = elem->nodes[6];

			pml_e->nodes[1].x = elem->nodes[6].x + X_pml;
			pml_e->nodes[1].y = elem->nodes[6].y;
			pml_e->nodes[1].z = elem->nodes[6].z;

			pml_e->nodes[2].x = elem->nodes[6].x + X_pml;
			pml_e->nodes[2].y = elem->nodes[6].y + Y_pml;
			pml_e->nodes[2].z = elem->nodes[6].z;

			pml_e->nodes[3].x = elem->nodes[6].x;
			pml_e->nodes[3].y = elem->nodes[6].y + Y_pml;
			pml_e->nodes[3].z = elem->nodes[6].z;


			pml_e->nodes[4].x = elem->nodes[6].x;
			pml_e->nodes[4].y = elem->nodes[6].y;
			pml_e->nodes[4].z = elem->nodes[6].z + Z_pml;

			pml_e->nodes[5].x = elem->nodes[6].x + X_pml;
			pml_e->nodes[5].y = elem->nodes[6].y;
			pml_e->nodes[5].z = elem->nodes[6].z + Z_pml;

			pml_e->nodes[6].x = elem->nodes[6].x + X_pml;
			pml_e->nodes[6].y = elem->nodes[6].y + Y_pml;
			pml_e->nodes[6].z = elem->nodes[6].z + Z_pml;

			pml_e->nodes[7].x = elem->nodes[6].x;
			pml_e->nodes[7].y = elem->nodes[6].y + Y_pml;
			pml_e->nodes[7].z = elem->nodes[6].z + Z_pml;

			pml_e->x = elem->x+1;
			pml_e->y = elem->y+1;
			pml_e->z = elem->z+1;

			//continue;
		}

		if (mask[PML_EDGE_Z0_X0]) {
			npmls[PML_EDGE_Z0_X0]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_EDGE_Z0_X0;
			pml_e->n_mat = 50;

			pml_e->nodes[5] = elem->nodes[0];
			pml_e->nodes[6] = elem->nodes[3];

			pml_e->nodes[0].x = elem->nodes[0].x - X_pml;
			pml_e->nodes[0].y = elem->nodes[0].y;
			pml_e->nodes[0].z = elem->nodes[0].z - Z_pml;

			pml_e->nodes[1].x = elem->nodes[0].x;
			pml_e->nodes[1].y = elem->nodes[0].y;
			pml_e->nodes[1].z = elem->nodes[0].z - Z_pml;

			pml_e->nodes[4].x = elem->nodes[0].x - X_pml;
			pml_e->nodes[4].y = elem->nodes[0].y;
			pml_e->nodes[4].z = elem->nodes[0].z;


			pml_e->nodes[2].x = elem->nodes[3].x;
			pml_e->nodes[2].y = elem->nodes[3].y;
			pml_e->nodes[2].z = elem->nodes[3].z - Z_pml;

			pml_e->nodes[3].x = elem->nodes[3].x - X_pml;
			pml_e->nodes[3].y = elem->nodes[3].y;
			pml_e->nodes[3].z = elem->nodes[3].z - Z_pml;

			pml_e->nodes[7].x = elem->nodes[3].x - X_pml;
			pml_e->nodes[7].y = elem->nodes[3].y;
			pml_e->nodes[7].z = elem->nodes[3].z;

			pml_e->x = elem->x-1;
			pml_e->y = elem->y;
			pml_e->z = elem->z-1;

			//continue;
		}

		if (mask[PML_EDGE_Z0_X1]) {
			npmls[PML_EDGE_Z0_X1]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_EDGE_Z0_X1;
			pml_e->n_mat = 50;

			pml_e->nodes[4] = elem->nodes[1];
			pml_e->nodes[7] = elem->nodes[2];

			pml_e->nodes[0].x = elem->nodes[1].x;
			pml_e->nodes[0].y = elem->nodes[1].y;
			pml_e->nodes[0].z = elem->nodes[1].z - Z_pml;

			pml_e->nodes[1].x = elem->nodes[1].x + X_pml;
			pml_e->nodes[1].y = elem->nodes[1].y;
			pml_e->nodes[1].z = elem->nodes[1].z - Z_pml;

			pml_e->nodes[5].x = elem->nodes[1].x + X_pml;
			pml_e->nodes[5].y = elem->nodes[1].y;
			pml_e->nodes[5].z = elem->nodes[1].z;


			pml_e->nodes[2].x = elem->nodes[2].x + X_pml;
			pml_e->nodes[2].y = elem->nodes[2].y;
			pml_e->nodes[2].z = elem->nodes[2].z - Z_pml;

			pml_e->nodes[3].x = elem->nodes[2].x;
			pml_e->nodes[3].y = elem->nodes[2].y;
			pml_e->nodes[3].z = elem->nodes[2].z - Z_pml;

			pml_e->nodes[6].x = elem->nodes[2].x + X_pml;
			pml_e->nodes[6].y = elem->nodes[2].y;
			pml_e->nodes[6].z = elem->nodes[2].z;

			pml_e->x = elem->x+1;
			pml_e->y = elem->y;
			pml_e->z = elem->z-1;

			//continue;
		}

		if (mask[PML_EDGE_Z0_Y0]) {
			npmls[PML_EDGE_Z0_Y0]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_EDGE_Z0_Y0;
			pml_e->n_mat = 50;


			pml_e->nodes[6] = elem->nodes[1];
			pml_e->nodes[7] = elem->nodes[0];

			pml_e->nodes[0].x = elem->nodes[0].x;
			pml_e->nodes[0].y = elem->nodes[0].y - Y_pml;
			pml_e->nodes[0].z = elem->nodes[0].z - Z_pml;

			pml_e->nodes[3].x = elem->nodes[0].x;
			pml_e->nodes[3].y = elem->nodes[0].y;
			pml_e->nodes[3].z = elem->nodes[0].z - Z_pml;

			pml_e->nodes[4].x = elem->nodes[0].x;
			pml_e->nodes[4].y = elem->nodes[0].y - Y_pml;
			pml_e->nodes[4].z = elem->nodes[0].z;


			pml_e->nodes[1].x = elem->nodes[1].x;
			pml_e->nodes[1].y = elem->nodes[1].y - Y_pml;
			pml_e->nodes[1].z = elem->nodes[1].z - Z_pml;

			pml_e->nodes[2].x = elem->nodes[1].x;
			pml_e->nodes[2].y = elem->nodes[1].y;
			pml_e->nodes[2].z = elem->nodes[1].z - Z_pml;

			pml_e->nodes[5].x = elem->nodes[1].x;
			pml_e->nodes[5].y = elem->nodes[1].y - Y_pml;
			pml_e->nodes[5].z = elem->nodes[1].z;

			pml_e->x = elem->x;
			pml_e->y = elem->y-1;
			pml_e->z = elem->z-1;

			//continue;
		}

		if (mask[PML_EDGE_Z0_Y1]) {
			npmls[PML_EDGE_Z0_Y1]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_EDGE_Z0_Y1;
			pml_e->n_mat = 50;

			pml_e->nodes[5] = elem->nodes[2];
			pml_e->nodes[4] = elem->nodes[3];

			pml_e->nodes[1].x = elem->nodes[2].x;
			pml_e->nodes[1].y = elem->nodes[2].y;
			pml_e->nodes[1].z = elem->nodes[2].z - Z_pml;

			pml_e->nodes[2].x = elem->nodes[2].x;
			pml_e->nodes[2].y = elem->nodes[2].y + Y_pml;
			pml_e->nodes[2].z = elem->nodes[2].z - Z_pml;

			pml_e->nodes[6].x = elem->nodes[2].x;
			pml_e->nodes[6].y = elem->nodes[2].y + Y_pml;
			pml_e->nodes[6].z = elem->nodes[2].z;


			pml_e->nodes[0].x = elem->nodes[3].x;
			pml_e->nodes[0].y = elem->nodes[3].y;
			pml_e->nodes[0].z = elem->nodes[3].z - Z_pml;

			pml_e->nodes[3].x = elem->nodes[3].x;
			pml_e->nodes[3].y = elem->nodes[3].y + Y_pml;
			pml_e->nodes[3].z = elem->nodes[3].z - Z_pml;

			pml_e->nodes[7].x = elem->nodes[3].x;
			pml_e->nodes[7].y = elem->nodes[3].y + Y_pml;
			pml_e->nodes[7].z = elem->nodes[3].z;

			pml_e->x = elem->x;
			pml_e->y = elem->y+1;
			pml_e->z = elem->z-1;

			//continue;
		}

		if (mask[PML_EDGE_X0_Y0]) {
			npmls[PML_EDGE_X0_Y0]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_EDGE_X0_Y0;
			pml_e->n_mat = 5;

			pml_e->nodes[2] = elem->nodes[0];
			pml_e->nodes[6] = elem->nodes[4];

			pml_e->nodes[0].x = elem->nodes[0].x - X_pml;
			pml_e->nodes[0].y = elem->nodes[0].y - Y_pml;
			pml_e->nodes[0].z = elem->nodes[0].z;

			pml_e->nodes[1].x = elem->nodes[0].x;
			pml_e->nodes[1].y = elem->nodes[0].y - Y_pml;
			pml_e->nodes[1].z = elem->nodes[0].z;

			pml_e->nodes[3].x = elem->nodes[0].x - X_pml;
			pml_e->nodes[3].y = elem->nodes[0].y;
			pml_e->nodes[3].z = elem->nodes[0].z;


			pml_e->nodes[4].x = elem->nodes[4].x - X_pml;
			pml_e->nodes[4].y = elem->nodes[4].y - Y_pml;
			pml_e->nodes[4].z = elem->nodes[4].z;

			pml_e->nodes[5].x = elem->nodes[4].x;
			pml_e->nodes[5].y = elem->nodes[4].y - Y_pml;
			pml_e->nodes[5].z = elem->nodes[4].z;

			pml_e->nodes[7].x = elem->nodes[4].x - X_pml;
			pml_e->nodes[7].y = elem->nodes[4].y;
			pml_e->nodes[7].z = elem->nodes[4].z;

			pml_e->x = elem->x-1;
			pml_e->y = elem->y-1;
			pml_e->z = elem->z;

			//continue;
		}

		if (mask[PML_EDGE_X0_Y1]) {
			npmls[PML_EDGE_X0_Y1]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_EDGE_X0_Y1;
			pml_e->n_mat = 6;

			pml_e->nodes[1] = elem->nodes[3];
			pml_e->nodes[5] = elem->nodes[7];

			pml_e->nodes[0].x = elem->nodes[3].x - X_pml;
			pml_e->nodes[0].y = elem->nodes[3].y;
			pml_e->nodes[0].z = elem->nodes[3].z;

			pml_e->nodes[2].x = elem->nodes[3].x;
			pml_e->nodes[2].y = elem->nodes[3].y + Y_pml;
			pml_e->nodes[2].z = elem->nodes[3].z;

			pml_e->nodes[3].x = elem->nodes[3].x - X_pml;
			pml_e->nodes[3].y = elem->nodes[3].y + Y_pml;
			pml_e->nodes[3].z = elem->nodes[3].z;


			pml_e->nodes[4].x = elem->nodes[7].x - X_pml;
			pml_e->nodes[4].y = elem->nodes[7].y;
			pml_e->nodes[4].z = elem->nodes[7].z;

			pml_e->nodes[6].x = elem->nodes[7].x;
			pml_e->nodes[6].y = elem->nodes[7].y + Y_pml;
			pml_e->nodes[6].z = elem->nodes[7].z;

			pml_e->nodes[7].x = elem->nodes[7].x - X_pml;
			pml_e->nodes[7].y = elem->nodes[7].y + Y_pml;
			pml_e->nodes[7].z = elem->nodes[7].z;

			pml_e->x = elem->x-1;
			pml_e->y = elem->y+1;
			pml_e->z = elem->z;

			//continue;
		}

		if (mask[PML_EDGE_X1_Y0]) {
			npmls[PML_EDGE_X1_Y0]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_EDGE_X1_Y0;
			pml_e->n_mat = 7;


			pml_e->nodes[3] = elem->nodes[1];
			pml_e->nodes[7] = elem->nodes[5];

			pml_e->nodes[0].x = elem->nodes[1].x;
			pml_e->nodes[0].y = elem->nodes[1].y - Y_pml;
			pml_e->nodes[0].z = elem->nodes[1].z;

			pml_e->nodes[1].x = elem->nodes[1].x + X_pml;
			pml_e->nodes[1].y = elem->nodes[1].y - Y_pml;
			pml_e->nodes[1].z = elem->nodes[1].z;

			pml_e->nodes[2].x = elem->nodes[1].x + X_pml;
			pml_e->nodes[2].y = elem->nodes[1].y;
			pml_e->nodes[2].z = elem->nodes[1].z;


			pml_e->nodes[4].x = elem->nodes[5].x;
			pml_e->nodes[4].y = elem->nodes[5].y - Y_pml;
			pml_e->nodes[4].z = elem->nodes[5].z;

			pml_e->nodes[5].x = elem->nodes[5].x + X_pml;
			pml_e->nodes[5].y = elem->nodes[5].y - Y_pml;
			pml_e->nodes[5].z = elem->nodes[5].z;

			pml_e->nodes[6].x = elem->nodes[5].x + X_pml;
			pml_e->nodes[6].y = elem->nodes[5].y;
			pml_e->nodes[6].z = elem->nodes[5].z;

			pml_e->x = elem->x+1;
			pml_e->y = elem->y-1;
			pml_e->z = elem->z;

			//continue;
		}

		if (mask[PML_EDGE_X1_Y1]) {
			npmls[PML_EDGE_X1_Y1]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->n_mat = 8;
			pml_e->pml_id = PML_EDGE_X1_Y1;

			pml_e->nodes[0] = elem->nodes[2];
			pml_e->nodes[4] = elem->nodes[6];

			pml_e->nodes[1].x = elem->nodes[2].x + X_pml;
			pml_e->nodes[1].y = elem->nodes[2].y;
			pml_e->nodes[1].z = elem->nodes[2].z;

			pml_e->nodes[2].x = elem->nodes[2].x + X_pml;
			pml_e->nodes[2].y = elem->nodes[2].y + Y_pml;
			pml_e->nodes[2].z = elem->nodes[2].z;

			pml_e->nodes[3].x = elem->nodes[2].x;
			pml_e->nodes[3].y = elem->nodes[2].y + Y_pml;
			pml_e->nodes[3].z = elem->nodes[2].z;


			pml_e->nodes[5].x = elem->nodes[6].x + X_pml;
			pml_e->nodes[5].y = elem->nodes[6].y;
			pml_e->nodes[5].z = elem->nodes[6].z;

			pml_e->nodes[6].x = elem->nodes[6].x + X_pml;
			pml_e->nodes[6].y = elem->nodes[6].y + Y_pml;
			pml_e->nodes[6].z = elem->nodes[6].z;

			pml_e->nodes[7].x = elem->nodes[6].x;
			pml_e->nodes[7].y = elem->nodes[6].y + Y_pml;
			pml_e->nodes[7].z = elem->nodes[6].z;

			pml_e->x = elem->x+1;
			pml_e->y = elem->y+1;
			pml_e->z = elem->z;

			//continue;
		}

		if (mask[PML_EDGE_Z1_X0]) {
			npmls[PML_EDGE_Z1_X0]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->n_mat = 9;
			pml_e->pml_id = PML_EDGE_Z1_X0;

			pml_e->nodes[1] = elem->nodes[4];
			pml_e->nodes[2] = elem->nodes[7];

			pml_e->nodes[0].x = elem->nodes[4].x - X_pml;
			pml_e->nodes[0].y = elem->nodes[4].y;
			pml_e->nodes[0].z = elem->nodes[4].z;

			pml_e->nodes[5].x = elem->nodes[4].x;
			pml_e->nodes[5].y = elem->nodes[4].y;
			pml_e->nodes[5].z = elem->nodes[4].z + Z_pml;

			pml_e->nodes[4].x = elem->nodes[4].x - X_pml;
			pml_e->nodes[4].y = elem->nodes[4].y;
			pml_e->nodes[4].z = elem->nodes[4].z + Z_pml;


			pml_e->nodes[3].x = elem->nodes[7].x - X_pml;
			pml_e->nodes[3].y = elem->nodes[7].y;
			pml_e->nodes[3].z = elem->nodes[7].z;

			pml_e->nodes[6].x = elem->nodes[7].x;
			pml_e->nodes[6].y = elem->nodes[7].y;
			pml_e->nodes[6].z = elem->nodes[7].z + Z_pml;

			pml_e->nodes[7].x = elem->nodes[7].x - X_pml;
			pml_e->nodes[7].y = elem->nodes[7].y;
			pml_e->nodes[7].z = elem->nodes[7].z + Z_pml;

			pml_e->x = elem->x-1;
			pml_e->y = elem->y;
			pml_e->z = elem->z+1;

			//continue;
		}

		if (mask[PML_EDGE_Z1_X1]) {
			npmls[PML_EDGE_Z1_X1]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->n_mat = 10;
			pml_e->pml_id = PML_EDGE_Z1_X1;

			pml_e->nodes[0] = elem->nodes[5];
			pml_e->nodes[3] = elem->nodes[6];

			pml_e->nodes[1].x = elem->nodes[5].x + X_pml;
			pml_e->nodes[1].y = elem->nodes[5].y;
			pml_e->nodes[1].z = elem->nodes[5].z;

			pml_e->nodes[5].x = elem->nodes[5].x + X_pml;
			pml_e->nodes[5].y = elem->nodes[5].y;
			pml_e->nodes[5].z = elem->nodes[5].z + Z_pml;

			pml_e->nodes[4].x = elem->nodes[5].x;
			pml_e->nodes[4].y = elem->nodes[5].y;
			pml_e->nodes[4].z = elem->nodes[5].z + Z_pml;


			pml_e->nodes[2].x = elem->nodes[6].x + X_pml;
			pml_e->nodes[2].y = elem->nodes[6].y;
			pml_e->nodes[2].z = elem->nodes[6].z;

			pml_e->nodes[6].x = elem->nodes[6].x + X_pml;
			pml_e->nodes[6].y = elem->nodes[6].y;
			pml_e->nodes[6].z = elem->nodes[6].z + Z_pml;

			pml_e->nodes[7].x = elem->nodes[6].x;
			pml_e->nodes[7].y = elem->nodes[6].y;
			pml_e->nodes[7].z = elem->nodes[6].z + Z_pml;

			pml_e->x = elem->x+1;
			pml_e->y = elem->y;
			pml_e->z = elem->z+1;

			//continue;
		}

		if (mask[PML_EDGE_Z1_Y0]) {
			npmls[PML_EDGE_Z1_Y0]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->n_mat = 11;
			pml_e->pml_id = PML_EDGE_Z1_Y0;

			pml_e->nodes[2] = elem->nodes[5];
			pml_e->nodes[3] = elem->nodes[4];

			pml_e->nodes[0].x = elem->nodes[4].x;
			pml_e->nodes[0].y = elem->nodes[4].y - Y_pml;
			pml_e->nodes[0].z = elem->nodes[4].z;

			pml_e->nodes[4].x = elem->nodes[4].x;
			pml_e->nodes[4].y = elem->nodes[4].y - Y_pml;
			pml_e->nodes[4].z = elem->nodes[4].z + Z_pml;

			pml_e->nodes[7].x = elem->nodes[4].x;
			pml_e->nodes[7].y = elem->nodes[4].y;
			pml_e->nodes[7].z = elem->nodes[4].z + Z_pml;

			pml_e->nodes[1].x = elem->nodes[5].x;
			pml_e->nodes[1].y = elem->nodes[5].y - Y_pml;
			pml_e->nodes[1].z = elem->nodes[5].z;

			pml_e->nodes[5].x = elem->nodes[5].x;
			pml_e->nodes[5].y = elem->nodes[5].y - Y_pml;
			pml_e->nodes[5].z = elem->nodes[5].z + Z_pml;

			pml_e->nodes[6].x = elem->nodes[5].x;
			pml_e->nodes[6].y = elem->nodes[5].y;
			pml_e->nodes[6].z = elem->nodes[5].z + Z_pml;

			pml_e->x = elem->x;
			pml_e->y = elem->y-1;
			pml_e->z = elem->z+1;

			//continue;
		}

		if (mask[PML_EDGE_Z1_Y1]) {
			npmls[PML_EDGE_Z1_Y1]++;

			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->n_mat = 12;
			pml_e->pml_id = PML_EDGE_Z1_Y1;

			pml_e->nodes[0] = elem->nodes[7];
			pml_e->nodes[1] = elem->nodes[6];

			pml_e->nodes[3].x = elem->nodes[7].x;
			pml_e->nodes[3].y = elem->nodes[7].y + Y_pml;
			pml_e->nodes[3].z = elem->nodes[7].z;

			pml_e->nodes[4].x = elem->nodes[7].x;
			pml_e->nodes[4].y = elem->nodes[7].y;
			pml_e->nodes[4].z = elem->nodes[7].z + Z_pml;

			pml_e->nodes[7].x = elem->nodes[7].x;
			pml_e->nodes[7].y = elem->nodes[7].y + Y_pml;
			pml_e->nodes[7].z = elem->nodes[7].z + Z_pml;

			pml_e->nodes[2].x = elem->nodes[6].x;
			pml_e->nodes[2].y = elem->nodes[6].y + Y_pml;
			pml_e->nodes[2].z = elem->nodes[6].z;

			pml_e->nodes[5].x = elem->nodes[6].x;
			pml_e->nodes[5].y = elem->nodes[6].y;
			pml_e->nodes[5].z = elem->nodes[6].z + Z_pml;

			pml_e->nodes[6].x = elem->nodes[6].x;
			pml_e->nodes[6].y = elem->nodes[6].y + Y_pml;
			pml_e->nodes[6].z = elem->nodes[6].z + Z_pml;

			pml_e->x = elem->x;
			pml_e->y = elem->y+1;
			pml_e->z = elem->z+1;

			//continue;
		}

		if (mask[PML_FACE_X0]) {
			npmls[PML_FACE_X0]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->n_mat = 13;
			pml_e->pml_id = PML_FACE_X0;

			pml_e->nodes[0].x = elem->nodes[0].x - X_pml;
			pml_e->nodes[0].y = elem->nodes[0].y;
			pml_e->nodes[0].z = elem->nodes[0].z;

			pml_e->nodes[3].x = elem->nodes[3].x - X_pml;
			pml_e->nodes[3].y = elem->nodes[3].y;
			pml_e->nodes[3].z = elem->nodes[3].z;

			pml_e->nodes[4].x = elem->nodes[4].x - X_pml;
			pml_e->nodes[4].y = elem->nodes[4].y;
			pml_e->nodes[4].z = elem->nodes[4].z;

			pml_e->nodes[7].x = elem->nodes[7].x - X_pml;
			pml_e->nodes[7].y = elem->nodes[7].y;
			pml_e->nodes[7].z = elem->nodes[7].z;

			pml_e->nodes[1] = elem->nodes[0];
			pml_e->nodes[2] = elem->nodes[3];
			pml_e->nodes[5] = elem->nodes[4];
			pml_e->nodes[6] = elem->nodes[7];

			pml_e->x = elem->x-1;
			pml_e->y = elem->y;
			pml_e->z = elem->z;

			//continue;
		}

		if (mask[PML_FACE_X1]) {
			npmls[PML_FACE_X1]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->n_mat = 14;
			pml_e->pml_id = PML_FACE_X1;

			pml_e->nodes[1].x = elem->nodes[1].x + X_pml;
			pml_e->nodes[1].y = elem->nodes[1].y;
			pml_e->nodes[1].z = elem->nodes[1].z;

			pml_e->nodes[2].x = elem->nodes[2].x + X_pml;
			pml_e->nodes[2].y = elem->nodes[2].y;
			pml_e->nodes[2].z = elem->nodes[2].z;

			pml_e->nodes[5].x = elem->nodes[5].x + X_pml;
			pml_e->nodes[5].y = elem->nodes[5].y;
			pml_e->nodes[5].z = elem->nodes[5].z;

			pml_e->nodes[6].x = elem->nodes[6].x + X_pml;
			pml_e->nodes[6].y = elem->nodes[6].y;
			pml_e->nodes[6].z = elem->nodes[6].z;

			pml_e->nodes[0] = elem->nodes[1];
			pml_e->nodes[3] = elem->nodes[2];
			pml_e->nodes[4] = elem->nodes[5];
			pml_e->nodes[7] = elem->nodes[6];

			pml_e->x = elem->x+1;
			pml_e->y = elem->y;
			pml_e->z = elem->z;

			//continue;
		}

		if (mask[PML_FACE_Y0]) {
			npmls[PML_FACE_Y0]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->n_mat = 15;
			pml_e->pml_id = PML_FACE_Y0;

			pml_e->nodes[0].x = elem->nodes[0].x;
			pml_e->nodes[0].y = elem->nodes[0].y - Y_pml;
			pml_e->nodes[0].z = elem->nodes[0].z;

			pml_e->nodes[1].x = elem->nodes[1].x;
			pml_e->nodes[1].y = elem->nodes[1].y - Y_pml;
			pml_e->nodes[1].z = elem->nodes[1].z;

			pml_e->nodes[5].x = elem->nodes[5].x;
			pml_e->nodes[5].y = elem->nodes[5].y - Y_pml;
			pml_e->nodes[5].z = elem->nodes[5].z;

			pml_e->nodes[4].x = elem->nodes[4].x;
			pml_e->nodes[4].y = elem->nodes[4].y - Y_pml;
			pml_e->nodes[4].z = elem->nodes[4].z;

			pml_e->nodes[2] = elem->nodes[1];
			pml_e->nodes[3] = elem->nodes[0];
			pml_e->nodes[6] = elem->nodes[5];
			pml_e->nodes[7] = elem->nodes[4];

			pml_e->x = elem->x;
			pml_e->y = elem->y-1;
			pml_e->z = elem->z;

			//continue;
		}

		if (mask[PML_FACE_Y1]) {
			npmls[PML_FACE_Y1]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->n_mat = 16;
			pml_e->pml_id = PML_FACE_Y1;


			pml_e->nodes[3].x = elem->nodes[3].x;
			pml_e->nodes[3].y = elem->nodes[3].y + Y_pml;
			pml_e->nodes[3].z = elem->nodes[3].z;

			pml_e->nodes[2].x = elem->nodes[2].x;
			pml_e->nodes[2].y = elem->nodes[2].y + Y_pml;
			pml_e->nodes[2].z = elem->nodes[2].z;

			pml_e->nodes[6].x = elem->nodes[6].x;
			pml_e->nodes[6].y = elem->nodes[6].y + Y_pml;
			pml_e->nodes[6].z = elem->nodes[6].z;

			pml_e->nodes[7].x = elem->nodes[7].x;
			pml_e->nodes[7].y = elem->nodes[7].y + Y_pml;
			pml_e->nodes[7].z = elem->nodes[7].z;

			pml_e->nodes[1] = elem->nodes[2];
			pml_e->nodes[0] = elem->nodes[3];
			pml_e->nodes[5] = elem->nodes[6];
			pml_e->nodes[4] = elem->nodes[7];

			pml_e->x = elem->x;
			pml_e->y = elem->y+1;
			pml_e->z = elem->z;

			//continue;
		}

		if (mask[PML_FACE_Z0]) {
			npmls[PML_FACE_Z0]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->pml_id = PML_FACE_Z0;
			pml_e->n_mat = 20;

			pml_e->nodes[0].x = elem->nodes[0].x;
			pml_e->nodes[0].y = elem->nodes[0].y;
			pml_e->nodes[0].z = elem->nodes[0].z - Z_pml;

			pml_e->nodes[1].x = elem->nodes[1].x;
			pml_e->nodes[1].y = elem->nodes[1].y;
			pml_e->nodes[1].z = elem->nodes[1].z - Z_pml;

			pml_e->nodes[2].x = elem->nodes[2].x;
			pml_e->nodes[2].y = elem->nodes[2].y;
			pml_e->nodes[2].z = elem->nodes[2].z - Z_pml;

			pml_e->nodes[3].x = elem->nodes[3].x;
			pml_e->nodes[3].y = elem->nodes[3].y;
			pml_e->nodes[3].z = elem->nodes[3].z - Z_pml;

			pml_e->nodes[4] = elem->nodes[0];
			pml_e->nodes[5] = elem->nodes[1];
			pml_e->nodes[6] = elem->nodes[2];
			pml_e->nodes[7] = elem->nodes[3];

			pml_e->x = elem->x;
			pml_e->y = elem->y;
			pml_e->z = elem->z-1;

			//continue;
		}

		if (mask[PML_FACE_Z1]) {
			npmls[PML_FACE_Z1]++;
			octant_t* pml_e = (octant_t*) sc_array_push(elements);

			elem->pml_id = 0;
			pml_e->n_mat = 17;
			pml_e->pml_id = PML_FACE_Z1;

			pml_e->nodes[4].x = elem->nodes[4].x;
			pml_e->nodes[4].y = elem->nodes[4].y;
			pml_e->nodes[4].z = elem->nodes[4].z + Z_pml;

			pml_e->nodes[5].x = elem->nodes[5].x;
			pml_e->nodes[5].y = elem->nodes[5].y;
			pml_e->nodes[5].z = elem->nodes[5].z + Z_pml;

			pml_e->nodes[6].x = elem->nodes[6].x;
			pml_e->nodes[6].y = elem->nodes[6].y;
			pml_e->nodes[6].z = elem->nodes[6].z + Z_pml;

			pml_e->nodes[7].x = elem->nodes[7].x;
			pml_e->nodes[7].y = elem->nodes[7].y;
			pml_e->nodes[7].z = elem->nodes[7].z + Z_pml;

			pml_e->nodes[0] = elem->nodes[4];
			pml_e->nodes[1] = elem->nodes[5];
			pml_e->nodes[2] = elem->nodes[6];
			pml_e->nodes[3] = elem->nodes[7];

			pml_e->x = elem->x;
			pml_e->y = elem->y;
			pml_e->z = elem->z+1;

			//continue;
		}
	}
}


void ExtrudePMLElements(hexa_tree_t* mesh) {

	double X_pml = 50000;
	double Y_pml = 50000;
	double Z_pml = 50000;
	int layers_x = 1;
	int layers_y = 1;
	int layers_z = 1;

	sc_array_t *elements = &mesh->elements;
	int nelem = elements->elem_count;

	for (int i = 0; i < nelem; ++i) {
		octant_t* elem = (octant_t*) sc_array_index(elements, i);

		if(elem->pml_id==PML_CORNER_X0Y0Z0){

		}else if(elem->pml_id==PML_CORNER_X1Y0Z0){

		}else if(elem->pml_id==PML_CORNER_X0Y1Z0){

		}else if(elem->pml_id==PML_CORNER_X1Y1Z0){

		}else if(elem->pml_id==PML_CORNER_X0Y0Z1){

		}else if(elem->pml_id==PML_CORNER_X1Y0Z1){

		}else if(elem->pml_id==PML_CORNER_X0Y1Z1){

		}else if(elem->pml_id==PML_CORNER_X1Y1Z1){

		}else if(elem->pml_id==PML_EDGE_Z0_X0){

		}else if(elem->pml_id==PML_EDGE_Z0_X1){

		}else if(elem->pml_id==PML_EDGE_Z0_Y0){

		}else if(elem->pml_id==PML_EDGE_Z0_Y1){

		}else if(elem->pml_id==PML_EDGE_X0_Y0){

		}else if(elem->pml_id==PML_EDGE_X0_Y1){

		}else if(elem->pml_id==PML_EDGE_X1_Y0){

		}else if(elem->pml_id==PML_EDGE_X1_Y1){

		}else if(elem->pml_id==PML_EDGE_Z1_X0){

		}else if(elem->pml_id==PML_EDGE_Z1_X1){

		}else if(elem->pml_id==PML_EDGE_Z1_Y0){

		}else if(elem->pml_id==PML_EDGE_Z1_Y1){

		}else if(elem->pml_id==PML_FACE_X0){

		}else if(elem->pml_id==PML_FACE_X1){

		}else if(elem->pml_id==PML_FACE_Y0){

		}else if(elem->pml_id==PML_FACE_Y1){

		}else if(elem->pml_id==PML_FACE_Z0){

		}else if(elem->pml_id==PML_FACE_Z1){

			for(int n_l = 0; n_l < layers_z; ++n_l ){

				octant_t* pml_e = (octant_t*) sc_array_push(elements);

				pml_e->pml_id = PML_FACE_Z1;

				pml_e->nodes[4].x = elem->nodes[4].x;
				pml_e->nodes[4].y = elem->nodes[4].y;
				pml_e->nodes[4].z = elem->nodes[4].z + 1 + n_l;

				pml_e->nodes[5].x = elem->nodes[5].x;
				pml_e->nodes[5].y = elem->nodes[5].y;
				pml_e->nodes[5].z = elem->nodes[5].z + 1 + n_l;

				pml_e->nodes[6].x = elem->nodes[6].x;
				pml_e->nodes[6].y = elem->nodes[6].y;
				pml_e->nodes[6].z = elem->nodes[6].z + 1 + n_l;

				pml_e->nodes[7].x = elem->nodes[7].x;
				pml_e->nodes[7].y = elem->nodes[7].y;
				pml_e->nodes[7].z = elem->nodes[7].z + 1 + n_l;

				pml_e->nodes[0] = elem->nodes[4];
				pml_e->nodes[1] = elem->nodes[5];
				pml_e->nodes[2] = elem->nodes[6];
				pml_e->nodes[3] = elem->nodes[7];

				pml_e->x = elem->x;
				pml_e->y = elem->y;
				pml_e->z = elem->z+1+n_l;


				/*
				pml_e->nodes[0].id = conn_p[0];
				pml_e->nodes[1].id = conn_p[0];
				pml_e->nodes[2].id = conn_p[0];
				pml_e->nodes[3].id = conn_p[0];

				pml_e->nodes[4].id = conn_p[0];
				pml_e->nodes[5].id = conn_p[1];
				pml_e->nodes[6].id = conn_p[2];
				pml_e->nodes[7].id = conn_p[3];
*/


			}


		}

	}
}

/*
int AddPointOnEdge_b(int* nodes, sc_hash_array_t* hash, int &npoints, GtsPoint *p, std::vector<double> &coords) {
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
*/
