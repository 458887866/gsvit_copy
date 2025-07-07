
#ifndef MAT_CONTROLS_H
#define MAT_CONTROLS_H

#include "global.h"

enum {
    ROW_MAT_OBJECT_INDEX = 0,    
};

void mat_controls_sensitive(XGControls *xgc, gint id);
void mat_controls_row_select_process(XGControls *xgc, gint id, GtkTreePath *path);
void mat_controls_changed(XGControls *xgc);

#endif  /* MAT_CONTROLS_H */