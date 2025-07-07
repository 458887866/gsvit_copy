
#ifndef PAR_TREE_H
#define PAR_TREE_H

#include "global.h"

enum {
    SET_POOL = 1,
    SET_BASIC = 2,
    SET_SF = 3,
    SET_TSF = 4,
    SET_TSFF = 5,
    SET_LTSF = 6,
    SET_LTSFF = 7,
    SET_BND = 8,
    SET_MEDIUM = 9,
    SET_OUT = 10,
    SET_NFFF = 11,
    SET_PNFFF = 12,
    SET_PSOURCE = 100,
    SET_POUT = 200,
    SET_IOUT = 300,
    SET_IIOUT = 400,
    SET_COUT = 500,
    SET_SOUT = 600,
    SET_FOUT = 700,    
    SET_GROW = 800,
    SET_ROUGHNESS = 900,
    SET_SPECTRAL = 1000,
    SET_EXPRESSION = 1100,
    SET_NFAREA = 99900,
    SET_NFFFP = 100000,
    SET_PNFAREA = 199900,
    SET_PNFFFP = 200000
};

/* Number Of Leaves */
#define NOL_MEDIUM_MAT_PROPS        4   // number of tree node MEDIUM non-volatile properties: voxel, vector, material mode and material smoothing */
#define NOL_OUTPUT_PROPS            1   // general output file
#define NOL_NFFF_PROPS              7   // NFFF box, NFFF skips
#define NOL_PNFFF_PROPS             3   // PNFFF box, PNFFF skips

void pt_assemble_row_text(XGControls *xgc, gchar *buff, gint buff_size, gint node_id, gint leaf_pos);
void pt_remeber_expanded_rows_and_selection(XGControls *xgc);
void pt_restore_expanded_rows_and_selection(XGControls *xgc);
void pt_create_tree(XGControls *xgc);

#endif  /* PAR_TREE_H */