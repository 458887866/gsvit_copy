
#ifndef XSVIT_H
#define XSVIT_H

#include "global.h"


/* menu File */
void file_new_cb(GtkWidget *widget, XGControls *xgc);
void file_open_cb(GtkWidget *widget, XGControls *xgc);
void file_save_as_cb(GtkWidget *widget, XGControls *xgc);
void file_save_cb(GtkWidget *widget, XGControls *xgc);
void file_save_mat_as_cb(GtkWidget *widget, XGControls *xgc);
void recent_chooser_item_activated_cb(GtkWidget *widget, XGControls *xgc);
void preferences_cb(GtkWidget *widget, XGControls *xgc);
void quit_cb(GtkWidget *widget, XGControls *xgc);

/* menu Edit parameters */
void par_add_src_point_cb(GtkWidget *widget, XGControls *xgc);
void par_add_src_sf_cb(GtkWidget *widget, XGControls *xgc);
void par_add_src_tsf_cb(GtkWidget *widget, XGControls *xgc);
void par_add_src_tsff_cb(GtkWidget *widget, XGControls *xgc);
void par_add_src_ltsf_cb(GtkWidget *widget, XGControls *xgc);
void par_add_src_ltsff_cb(GtkWidget *widget, XGControls *xgc);
void par_add_out_point_cb(GtkWidget *widget, XGControls *xgc);
void par_add_out_image_cb(GtkWidget *widget, XGControls *xgc);
void par_add_out_plane_cb(GtkWidget *widget, XGControls *xgc);
void par_add_out_volume_cb(GtkWidget *widget, XGControls *xgc);
void par_add_out_sum_cb(GtkWidget *widget, XGControls *xgc);
void par_add_out_force_cb(GtkWidget *widget, XGControls *xgc);
void par_add_grow_cb(GtkWidget *widget, XGControls *xgc);
void par_add_roughness_cb(GtkWidget *widget, XGControls *xgc);
void par_add_spectral_cb(GtkWidget *widget, XGControls *xgc);
void par_add_expression_cb(GtkWidget *widget, XGControls *xgc);
void par_add_nfffp_cb(GtkWidget *widget, XGControls *xgc);
void par_add_nfffa_cb(GtkWidget *widget, XGControls *xgc);
void par_add_pnfffp_cb(GtkWidget *widget, XGControls *xgc);
void par_add_pnfffa_cb(GtkWidget *widget, XGControls *xgc);
void par_row_remove_cb(GtkButton *button, XGControls *xgc);

/* menu Edit material objects */
void mat_add_sphere_cb(GtkWidget *widget, XGControls *xgc);
void mat_add_box_cb(GtkWidget *widget, XGControls *xgc);
void mat_add_cylinder_cb(GtkWidget *widget, XGControls *xgc);
void mat_add_cone_cb(GtkWidget *widget, XGControls *xgc);
void mat_add_rcone_cb(GtkWidget *widget, XGControls *xgc);
void mat_add_tetrahedron_cb(GtkWidget *widget, XGControls *xgc);
void mat_row_remove_cb(GtkButton *button, XGControls *xgc);

/* menu Execute */
void gsvit_run_cb(GtkWidget *widget, XGControls *xgc);
void gsvit_stop_cb(GtkWidget *widget, XGControls *xgc);
void gwydd_cb(GtkWidget *widget, XGControls *xgc);

/* menu Help */
void help_x_cb(GtkWidget *widget, XGControls *xgc);
void help_g_cb(GtkWidget *widget, XGControls *xgc);
void about_cb(GtkWidget *widget, XGControls *xgc);

/* toolbar */
void check_files_consistency_cb(GtkWidget *widget, XGControls *xgc);


#endif  /* XSVIT_H */
