
#ifndef WRITE_H
#define WRITE_H

#include "global.h"

gchar* get_temporary_par_filename(XGControls *xgc);
gchar* get_temporary_mat_filename(XGControls *xgc);
gchar* planestring(gint i, gint j, gint k);
gchar* scpstring(SvSumType component);
gchar* vcpstring(SvOutputVolumeType component);
gchar* cpstring(SvCompType component);

gboolean write_parfile(XGControls *xgc);
gboolean write_matfile(XGControls *xgc);

gboolean write_textbuffer_to_parfile(XGControls *xgc);
gboolean write_textbuffer_to_matfile(XGControls *xgc);

#endif  /* WRITE_H */