
#ifndef MAT_TREE_H
#define MAT_TREE_H

#include "xsvit.h"

enum {
    SET_MAT_SPHERE = 0,
    SET_MAT_VOXEL = 1000,
    SET_MAT_CYLINDER = 2000,
    SET_MAT_CONE = 3000,
    SET_MAT_RCONE = 4000,
    SET_MAT_GWYDD = 5000,
    SET_MAT_MESH = 6000,
    SET_MAT_TETRAHEDRON = 100000
};

void mt_tree_assemble_row_text(XGControls *xgc, gchar *buff, gint buff_size, gint node_id, gint leaf_pos);
void mt_remeber_expanded_rows_and_selection(XGControls *xgc);
void mt_restore_expanded_rows_and_selection(XGControls *xgc);
void mt_create_tree(XGControls *xgc);

#endif  /* MAT_TREE_H */