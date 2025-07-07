
#ifndef PAR_CONTROLS_H
#define PAR_CONTROLS_H

#include "global.h"

enum {
    ROW_POOL_INDEX = 0,
    ROW_BASIC_INDEX = 1,
    ROW_SOURCES_INDEX = 2,
    ROW_BOUNDARY_INDEX = 3,
    ROW_MEDIA_INDEX = 4,
    ROW_OUTPUTS_INDEX = 5,
    ROW_NFFF_INDEX = 6,
    ROW_PNFFF_INDEX = 7,

    ROW_POOL_SIZE = 100,
    ROW_POOL_SPACING = 101,

    ROW_BASIC_STEPS = 200,
    ROW_BASIC_VERBOSE = 201,
    ROW_BASIC_THREADS = 202,
    ROW_BASIC_UGPU = 203,
    ROW_BASIC_GPU = 204,

    ROW_SOURCES_POINT = 300,
    ROW_SOURCES_SF = 301,
    ROW_SOURCES_TSF = 302,
    ROW_SOURCES_TSFF = 303,
    ROW_SOURCES_LTSF = 304,
    ROW_SOURCES_LTSFF = 305,

    //ROW_BOUNDARY_X0 = 301,
};

void par_controls_changed(XGControls *xgc);
void par_controls_row_select_process(XGControls *xgc, gint id, GtkTreePath *path);
void par_controls_sensitive(XGControls *xgc, gint id);

#endif  /* PAR_CONTROLS_H */