
/*
 *  Copyright (C) 2013 Petr Klapetek
 *  E-mail: klapetek@gwyddion.net.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 USA
 */


 /*  svit.c :
  *  main file calling all the algorithms, including main() function
  */

#include <gtk/gtk.h>
#include <gtk/gtkgl.h>
#include <locale.h>
#include <gdk/gdkkeysyms.h>

#ifdef G_OS_WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <wait.h>
#endif

#include <stdio.h>
#include <glib.h>
#include <glib/gstdio.h>
#include <stdlib.h>
#include <string.h>
#if defined (_WIN32) && defined (_MSC_VER)
#include "../config_win.h"
#else
#include "../config.h"
#endif
#include "xsvit.h"
#include "xs_dialogs.h"
#include "xs_layer_treeview.h"
#include "xs_par_tree.h"
#include "xs_par_controls.h"
#include "xs_mat_tree.h"
#include "xs_mat_controls.h"
#include "xs_main_menu.h"
#include "xs_scene.h"
#include <math.h>
#include "constants.h"
#include "messages.h"

#include <app/gwyapp.h>
#include <libgwymodule/gwymodule.h>

/* temporary includes */
#include "write.h"


#define ENTRY_FILENAME_CHARS        16  // number of char in GTK entry box


/* controls - par file */
static void par_controls_changed_cb(XGControls *xgc);
static gboolean on_key_press_par_tree_cb (GtkWidget *widget, GdkEventKey *event, XGControls *xgc);

/* controls - mat file */
static void mat_controls_changed_cb(XGControls *xgc);
static gboolean on_key_press_mat_tree_cb (GtkWidget *widget, GdkEventKey *event, XGControls *xgc);

gboolean update_graph(gpointer data);
gboolean update_image(gpointer data);

static void set_all_controls_invisible(XGControls *xgc);
static void set_all_mat_controls_invisible(XGControls *xgc);
static void set_controls_visible(XGControls *xgc, gint id);
static void mat_set_controls_visible(XGControls *xgc, gint id);
static void set_computation_controls_sensitive(XGControls *xgc);
static void ss_choose_source_file_cb(GtkWidget *widget, XGControls *xgc);
static void so_choose_output_file_cb(GtkWidget *widget, XGControls *xgc);
static void nfffp_choose_output_file_cb(GtkWidget *widget, XGControls *xgc);
static void material_choose_voxel_file_cb(GtkWidget *widget, XGControls *xgc);
static void material_choose_vector_file_cb(GtkWidget *widget, XGControls *xgc);

static void file_load(XGControls *xgc, gchar *filename, gchar* uri);
static void file_save(XGControls *xgc);


void puthelp()
{
    fprintf(stderr, "XSvit - graphical fronted for running GSvit FDTD solver\n");
    fprintf(stderr, "Usage: ./xsvit\n\n");
    fprintf(stderr, "Report bugs to %s\n", PACKAGE_BUGREPORT);
}

/*
void edit_undo_cb(GtkWidget *widget, XGControls *xgc)
{
   xgc->undo = 0;
   gtk_widget_set_sensitive(GTK_WIDGET(xgc->undo_tool_item), FALSE);
   gtk_widget_set_sensitive(GTK_WIDGET(xgc->undo_menu_item), FALSE);
   gtk_widget_set_sensitive(GTK_WIDGET(xgc->redo_tool_item), TRUE);
   gtk_widget_set_sensitive(GTK_WIDGET(xgc->redo_menu_item), TRUE);

}

void edit_redo_cb(GtkWidget *widget, XGControls *xgc)
{
   xgc->undo = 1;
   gtk_widget_set_sensitive(GTK_WIDGET(xgc->undo_tool_item), TRUE);
   gtk_widget_set_sensitive(GTK_WIDGET(xgc->undo_menu_item), TRUE);
   gtk_widget_set_sensitive(GTK_WIDGET(xgc->redo_tool_item), FALSE);
   gtk_widget_set_sensitive(GTK_WIDGET(xgc->redo_menu_item), FALSE);

}*/


/*static inline gboolean
ggwy_data_field_inside(GwyDataField *data_field, gint i, gint j)
{
    if (i >= 0 && j >= 0 && i < data_field->xres && j < data_field->yres)
        return TRUE;
    else
        return FALSE;
}*/


void mat_buffer_changed(XGControls *xgc)
{
    gchar *filename = NULL;
    gint type, testout, xres, yres, zres, j, k;
    gdouble xreal, yreal, zreal, val;
    GtkWidget *dialog;
    GtkTextIter   start;
    GtkTextIter   end;
    GError* glib_error = 0;
    FILE *fr = NULL;
    gint ret;
    SvSphere sx;
    SvVoxel vx;
    SvCylinder cx;
    SvCone cnx;
    SvRCone rcnx;
    SvTetrahedron tx;
    SvGwydd gx;
    gint ngtot = 0;
    SvMatProp mat;
    GLfloat xshift, yshift, zshift, xmult, ymult, zmult;
    gchar buff[256], matfile[256];
    gchar filebase[256];
    gint attribute_pos, attribute_val;
    gint i, m, gwydd_channel, isobject;
    gdouble mask;
    GtkTreeIter   iter, child;
    SvDCube *buf;
    gdouble min, max;    

    if (xgc->curdir)
        filename = g_build_filename(xgc->curdir, "xsvit_temporary.mat", NULL);
    else {
#ifdef G_OS_WIN32
        filename = g_build_filename(g_get_home_dir(), ".xsvit", "xsvit_temporary.mat", NULL);
#else
        gint ft;
        if (!(ft = g_file_open_tmp("xsvit.matfile.XXXXXX", &filename, &glib_error))) {
            dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                                            GTK_DIALOG_DESTROY_WITH_PARENT,
                                            GTK_MESSAGE_ERROR,
                                            GTK_BUTTONS_CLOSE,
                                            "Error opening temporary file: %s",
                                            g_strdup(glib_error->message));
            gtk_dialog_run(GTK_DIALOG(dialog));
            gtk_widget_destroy(dialog);

        }
        close(ft);
#endif
    }

    if (filename) {
        gtk_text_buffer_get_start_iter(xgc->tb_mat, &start);
        gtk_text_buffer_get_end_iter(xgc->tb_mat, &end);
        if (!g_file_set_contents(filename, gtk_text_buffer_get_text(xgc->tb_mat, &start, &end, TRUE), -1, &glib_error)) {
            dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                                            GTK_DIALOG_DESTROY_WITH_PARENT,
                                            GTK_MESSAGE_ERROR,
                                            GTK_BUTTONS_CLOSE,
                                            "Error writing to temporary file: %s",
                                            g_strdup(glib_error->message));
            gtk_dialog_run(GTK_DIALOG(dialog));
            gtk_widget_destroy(dialog);
        }        

        { //do this always?? 

            xgc->data.set_mat.nmats = 0;
            xgc->data.set_mat.materials = (gchar **)g_malloc(100 * sizeof(gchar *)); //FIXME: limited size array
            xgc->data.set_mat.overpos = (gint *)g_malloc(100 * sizeof(gint));
            for (i = 0; i < 100; i++)
                xgc->data.set_mat.overpos[i] = -1;

            if ((testout = test_script_file(filename)) == SV_MATMODE_DATABASE) {
                if (override_matfile_from_database(filename, 600e-9, 600e-9, 600e-9, 12345, &(xgc->data.set_mat.nmats), xgc->data.set_mat.materials, xgc->data.set_mat.overpos, TRUE)) { //FIXME, here should be set->sm.localfile
                    g_snprintf(buff, sizeof(buff), "Error parsing material from database.");
                    gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
                } else {
                    g_snprintf(buff, sizeof(buff), "Material loaded from database (%d entries listed).", xgc->data.set_mat.nmats);
                    gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);

                    /*for (i=0; i<xgc->nmats; i++)
                        printf("Material %d: %s to object %d\n", i, xgc->materials[i], xgc->overpos[i]); */
                }

                //FIXME! use real wavelengths
                g_snprintf(matfile, sizeof(matfile), "tmp_matfile_%06d", 12345);
                filename = g_strdup(matfile);
            }
            if (testout == SV_MATMODE_ERROR) {
                printf("Error in material file\n");

                g_snprintf(buff, sizeof(buff), "Error parsing material file.");
                gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);

                return;
            }


            fr = fopen(filename, "r");
            if (fr == NULL) {
                fprintf(stderr, "Error: cannot open material vector file for loading.\n");
                g_snprintf(buff, sizeof(buff), "Error parsing material file, see console for output");
                //gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
            } else {
                /* alloc data*/
                alloc_set_par_mat(&xgc->data.set, &xgc->data.set_mat, FALSE);

                while ((ret = fscanf(fr, "%d", &type)) != EOF) { //TODO use proper locale handling functions here
                    if (ret != 1) {
                        fprintf(stderr, "Error parsing material file\n");
                        goto exit;
                    }

                    switch (type) {
                        case SV_ENTITY_SPHERE:
                            init_settings_mat_prop(&sx.mat);

                            if (scan_point(fr, sx.pnt1, sx.pnt1 + 1, sx.pnt1 + 2)) {
                                fprintf(stderr, "Error: cannot scan sphere center from material file\n");
                                goto exit;
                            }
                            if (scan_radius(fr, &sx.radius)) {
                                fprintf(stderr, "Error: cannot scan sphere radius from material file\n");
                                goto exit;
                            }
                            if (scan_material(fr, &sx.mat, 1)) {
                                fprintf(stderr, "Error: cannot scan sphere material from material file\n");
                                goto exit;
                            }
                            //if (sx.mat.type!=0) sv_pool_add_material(mp, &sx.mat);
                            //sx.n = ngtot++;
                            g_array_append_val(xgc->data.set_mat.spheres, sx);
                            break;

                        case SV_ENTITY_CYLINDER:
                            if (scan_point(fr, cx.pnt1, cx.pnt1 + 1, cx.pnt1 + 2)) {
                                fprintf(stderr, "Error: cannot scan cylinder center1 from material file\n");
                                goto exit;
                            }
                            if (scan_point(fr, cx.pnt2, cx.pnt2 + 1, cx.pnt2 + 2)) {
                                fprintf(stderr, "Error: cannot scan cylinder center2 from material file\n");
                                goto exit;
                            }
                            if (scan_radius(fr, &cx.radius)) {
                                fprintf(stderr, "Error: cannot scan cylinder radius from material file\n");
                                goto exit;
                            }
                            if (scan_material(fr, &cx.mat, 1)) {
                                fprintf(stderr, "Error: cannot scan cylinder material from material file\n");
                                goto exit;
                            }
                            //if (cx.mat.type!=0) sv_pool_add_material(mp, &cx.mat);
                            cx.n = ngtot++;
                            g_array_append_val(xgc->data.set_mat.cylinders, cx);
                            break;

                        case SV_ENTITY_CONE:
                            if (scan_point(fr, cnx.pnt1, cnx.pnt1 + 1, cnx.pnt1 + 2)) {
                                fprintf(stderr, "Error: cannot scan cone center1 from material file\n");
                                goto exit;
                            }
                            if (scan_point(fr, cnx.pnt2, cnx.pnt2 + 1, cnx.pnt2 + 2)) {
                                fprintf(stderr, "Error: cannot scan cone center2 from material file\n");
                                goto exit;
                            }
                            if (scan_radius(fr, &cnx.radius)) {
                                fprintf(stderr, "Error: cannot scan cone radius from material file\n");
                                goto exit;
                            }
                            if (scan_material(fr, &cnx.mat, 1)) {
                                fprintf(stderr, "Error: cannot scan cone material from material file\n");
                                goto exit;
                            }
                            //if (cnx.mat.type!=0) sv_pool_add_material(mp, &cnx.mat);
                            cnx.n = ngtot++;
                            g_array_append_val(xgc->data.set_mat.cones, cnx);
                            break;

                        case SV_ENTITY_RCONE:
                            if (scan_point(fr, rcnx.pnt1, rcnx.pnt1 + 1, rcnx.pnt1 + 2)) {
                                fprintf(stderr, "Error: cannot scan cut cone center1 from material file\n");
                                goto exit;
                            }
                            if (scan_point(fr, rcnx.pnt2, rcnx.pnt2 + 1, rcnx.pnt2 + 2)) {
                                fprintf(stderr, "Error: cannot scan cut cone center2 from material file\n");
                                goto exit;
                            }
                            if (scan_radius(fr, &rcnx.radius1)) {
                                fprintf(stderr, "Error: cannot scan cut cone radius1 from material file\n");
                                goto exit;
                            }
                            if (scan_radius(fr, &rcnx.radius2)) {
                                fprintf(stderr, "Error: cannot scan cut cone radius2 from material file\n");
                                goto exit;
                            }
                            if (scan_material(fr, &rcnx.mat, 1)) {
                                fprintf(stderr, "Error: cannot scan cut cone material from material file\n");
                                goto exit;
                            }
                            //if (rcnx.mat.type!=0) sv_pool_add_material(mp, &rcnx.mat);
                            rcnx.n = ngtot++;
                            g_array_append_val(xgc->data.set_mat.rcones, rcnx);
                            break;

                        case SV_ENTITY_VOXEL:
                            if (scan_point(fr, vx.pnt1, vx.pnt1 + 1, vx.pnt1 + 2)) { 
                                fprintf(stderr, "Error: cannot scan voxel start vertex from material file\n"); 
                                goto exit; 
                            }
                            if (scan_point(fr, vx.pnt2, vx.pnt2 + 1, vx.pnt2 + 2)) { 
                                fprintf(stderr, "Error: cannot scan voxel end vertex from material file\n"); 
                                goto exit; 
                            }
                            if (scan_material(fr, &vx.mat, 1)) {
                                fprintf(stderr, "Error: cannot scan voxel material from material file\n"); 
                                goto exit;
                            }
                            //if (vx.mat.type!=0) sv_pool_add_material(mp, &vx.mat);
                            vx.n = ngtot++;
                            g_array_append_val(xgc->data.set_mat.voxels, vx);
                            break;

                        case SV_ENTITY_TETRAHEDRON:
                            if (scan_point(fr, tx.pnt1, tx.pnt1 + 1, tx.pnt1 + 2)) {
                                fprintf(stderr, "Error: cannot scan tetrahedron point1 from material file\n");
                                goto exit;
                            }
                            if (scan_point(fr, tx.pnt2, tx.pnt2 + 1, tx.pnt2 + 2)) {
                                fprintf(stderr, "Error: cannot scan tetrahedron point2 from material file\n");
                                goto exit;
                            }
                            if (scan_point(fr, tx.pnt3, tx.pnt3 + 1, tx.pnt3 + 2)) {
                                fprintf(stderr, "Error: cannot scan tetrahedron point3 from material file\n");
                                goto exit;
                            }
                            if (scan_point(fr, tx.pnt4, tx.pnt4 + 1, tx.pnt4 + 2)) {
                                fprintf(stderr, "Error: cannot scan tetrahedron point4 from material file\n");
                                goto exit;
                            }
                            if (scan_material(fr, &tx.mat, 1)) {
                                fprintf(stderr, "Error: cannot scan tetrahedron material from material file\n");
                                goto exit;
                            }
                            //if (tx.mat.type!=0) sv_pool_add_material(mp, &tx.mat);
                            tx.setpart = 0;
                            tx.n = ngtot++;
                            g_array_append_val(xgc->data.set_mat.tetrahedrons, tx);
                            break;

                        case SV_ENTITY_TETGEN:
                            fscanf(fr, "%s", filebase);             //TODO this should be locale dependent
                            if (scan_int(fr, &attribute_pos)) {
                                fprintf(stderr, "Error: cannot scan tetrahedral mesh attribute_number from material file\n");
                                goto exit;
                            }
                            if (scan_int(fr, &attribute_val)) {
                                fprintf(stderr, "Error: cannot scan tetrahedral mesh material_index from material file\n");
                                goto exit;
                            }
                            if (scan_point(fr, &xshift, &yshift, &zshift)) {
                                fprintf(stderr, "Error: cannot scan tetrahedral mesh shift from material file\n");
                                goto exit;
                            }
                            if (scan_point(fr, &xmult, &ymult, &zmult)) {
                                fprintf(stderr, "Error: cannot scan tetrahedral mesh mult from material file\n");
                                goto exit;
                            }
                            if (scan_material(fr, &mat, 1)) {
                                fprintf(stderr, "Error: cannot scan tetrahedral mesh material from material file\n");
                                goto exit;
                            }
                            //if (mat.type!=0) sv_pool_add_material(mp, &mat);

                            xgc->data.set_mat.tetgen_attribute_pos[xgc->data.set_mat.ntetgens] = attribute_pos;
                            xgc->data.set_mat.tetgen_attribute_val[xgc->data.set_mat.ntetgens] = attribute_val;
                            xgc->data.set_mat.tetgen_start[xgc->data.set_mat.ntetgens] = ngtot;
                            xgc->data.set_mat.tetgen_xshift[xgc->data.set_mat.ntetgens] = xshift;
                            xgc->data.set_mat.tetgen_yshift[xgc->data.set_mat.ntetgens] = yshift;
                            xgc->data.set_mat.tetgen_zshift[xgc->data.set_mat.ntetgens] = zshift;
                            xgc->data.set_mat.tetgen_xmult[xgc->data.set_mat.ntetgens] = xmult;
                            xgc->data.set_mat.tetgen_ymult[xgc->data.set_mat.ntetgens] = ymult;
                            xgc->data.set_mat.tetgen_zmult[xgc->data.set_mat.ntetgens] = zmult;
                            g_snprintf(xgc->data.set_mat.tetgen_filebase[xgc->data.set_mat.ntetgens], 100, "%s", filebase);
                            if (load_tets(filebase, xgc->data.set_mat.tetrahedrons, attribute_pos, attribute_val, xshift, yshift, zshift, xmult, ymult, zmult, 3, &mat, &ngtot) != 0) {
                                g_snprintf(buff, sizeof(buff), "Error parsing tetrahedral mesh from %s", filebase);
                                gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
                                goto exit;
                            } else {
                                g_snprintf(buff, sizeof(buff), "Loaded tetrahedral mesh from %s", filebase);
                                gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
                            }
                            xgc->data.set_mat.tetgen_end[xgc->data.set_mat.ntetgens] = ngtot;
                            xgc->data.set_mat.ntetgens++;
                            break;

                        case SV_ENTITY_GWYDDION:
                            fscanf(fr, "%s", filebase);             //TODO this should be locale dependent
                            g_snprintf(gx.filebase, 256, "%s", filebase);
                            if (scan_int(fr, &gwydd_channel)) {
                                fprintf(stderr, "Error: cannot scan Gwyddion field channel from material file\n");
                                goto exit;
                            }
                            gx.channel = gwydd_channel;
                            if (scan_int(fr, &gx.mask)) {
                                fprintf(stderr, "Error: cannot scan Gwyddion field mask from material file\n");
                                goto exit;
                            }
                            if (scan_int(fr, &gx.i)) {
                                fprintf(stderr, "Error: cannot scan Gwyddion field i from material file\n");
                                goto exit;
                            }
                            if (scan_int(fr, &gx.j)) {
                                fprintf(stderr, "Error: cannot scan Gwyddion field j from material file\n");
                                goto exit;
                            }
                            if (scan_int(fr, &gx.k)) {
                                fprintf(stderr, "Error: cannot scan Gwyddion field k from material file\n");
                                goto exit;
                            }
                            if (scan_point(fr, &gx.xoffset, &gx.yoffset, &gx.zoffset)) {
                                fprintf(stderr, "Error: cannot scan Gwyddion field offset from material file\n");
                                goto exit;
                            }
                            if (scan_int(fr, &gx.depth)) {
                                fprintf(stderr, "Error: cannot scan Gwyddion field depth from material file\n");
                                goto exit;
                            }
                            if (scan_material(fr, &gx.mat, 1)) {
                                fprintf(stderr, "Error: cannot scan Gwyddion field material from material file\n");
                                goto exit;
                            }
                            //if (mat.type!=0) sv_pool_add_material(mp, &gx.mat);

                            load_gwydd(filebase, &gx, 0, gwydd_channel);

                            if (gx.dfield == NULL || (!GWY_IS_DATA_FIELD(gx.dfield))) {
                                fprintf(stderr, "Error: no valid Gwyddion datafield\n");
                            } else {
                                //                              printf("loaded %d %d\n", gwy_data_field_get_xres(gx.dfield),  gwy_data_field_get_yres(gx.dfield));
                                gx.n = ngtot++;
                                g_array_append_val(xgc->data.set_mat.gwydds, gx);

                                if (xgc->data.set_mat.gwyddata) 
                                    sv_dcube_free(xgc->data.set_mat.gwyddata);
                                xgc->data.set_mat.gwyddata = sv_dcube_new(xgc->data.set.sp.xres, xgc->data.set.sp.yres, xgc->data.set.sp.zres,
                                                                          xgc->data.set.sp.xres, xgc->data.set.sp.yres, xgc->data.set.sp.zres, 1);
                                buf = sv_dcube_new(xgc->data.set.sp.xres, xgc->data.set.sp.yres, xgc->data.set.sp.zres,
                                                   xgc->data.set.sp.xres, xgc->data.set.sp.yres, xgc->data.set.sp.zres, 1);

                                xres = xgc->data.set.sp.xres;
                                xreal = xgc->data.set.sp.dx*xres;

                                yres = xgc->data.set.sp.yres;
                                yreal = xgc->data.set.sp.dy*yres;

                                zres = xgc->data.set.sp.zres;
                                zreal = xgc->data.set.sp.dz*zres;

                                min = gwy_data_field_get_min(gx.dfield);
                                max = gwy_data_field_get_max(gx.dfield);

                                if (gx.i != -1) {
                                    xgc->data.set_mat.gwydd_ymin[xgc->data.set_mat.gwydds->len - 1] = (-gx.yoffset) / xgc->data.set.sp.dy;
                                    xgc->data.set_mat.gwydd_zmin[xgc->data.set_mat.gwydds->len - 1] = (-gx.zoffset) / xgc->data.set.sp.dz;
                                    xgc->data.set_mat.gwydd_ymax[xgc->data.set_mat.gwydds->len - 1] = ((gdouble)yres*xgc->data.set.sp.dy - gx.yoffset) / xgc->data.set.sp.dy;
                                    xgc->data.set_mat.gwydd_zmax[xgc->data.set_mat.gwydds->len - 1] = ((gdouble)zres*xgc->data.set.sp.dz - gx.zoffset) / xgc->data.set.sp.dz;
                                    xgc->data.set_mat.gwydd_xmin[xgc->data.set_mat.gwydds->len - 1] = gx.i + (gdouble)xres / xreal*min;
                                    xgc->data.set_mat.gwydd_xmax[xgc->data.set_mat.gwydds->len - 1] = gx.i + (gdouble)xres / xreal*max;
                                }
                                if (gx.j != -1) {
                                    xgc->data.set_mat.gwydd_xmin[xgc->data.set_mat.gwydds->len - 1] = (-gx.xoffset) / xgc->data.set.sp.dx;
                                    xgc->data.set_mat.gwydd_zmin[xgc->data.set_mat.gwydds->len - 1] = (-gx.zoffset) / xgc->data.set.sp.dz;
                                    xgc->data.set_mat.gwydd_xmax[xgc->data.set_mat.gwydds->len - 1] = ((gdouble)xres*xgc->data.set.sp.dx - gx.xoffset) / xgc->data.set.sp.dx;
                                    xgc->data.set_mat.gwydd_zmax[xgc->data.set_mat.gwydds->len - 1] = ((gdouble)zres*xgc->data.set.sp.dz - gx.zoffset) / xgc->data.set.sp.dz;
                                    xgc->data.set_mat.gwydd_ymin[xgc->data.set_mat.gwydds->len - 1] = gx.j + (gdouble)yres / yreal*min;
                                    xgc->data.set_mat.gwydd_ymax[xgc->data.set_mat.gwydds->len - 1] = gx.j + (gdouble)yres / yreal*max;
                                }
                                if (gx.k != -1) {
                                    xgc->data.set_mat.gwydd_xmin[xgc->data.set_mat.gwydds->len - 1] = (-gx.xoffset) / xgc->data.set.sp.dx;
                                    xgc->data.set_mat.gwydd_ymin[xgc->data.set_mat.gwydds->len - 1] = (-gx.yoffset) / xgc->data.set.sp.dy;
                                    xgc->data.set_mat.gwydd_xmax[xgc->data.set_mat.gwydds->len - 1] = ((gdouble)xres*xgc->data.set.sp.dx - gx.xoffset) / xgc->data.set.sp.dx;
                                    xgc->data.set_mat.gwydd_ymax[xgc->data.set_mat.gwydds->len - 1] = ((gdouble)yres*xgc->data.set.sp.dy - gx.yoffset) / xgc->data.set.sp.dy;
                                    xgc->data.set_mat.gwydd_zmin[xgc->data.set_mat.gwydds->len - 1] = gx.k + (gdouble)zres / zreal*min;
                                    xgc->data.set_mat.gwydd_zmax[xgc->data.set_mat.gwydds->len - 1] = gx.k + (gdouble)zres / zreal*max;
                                }
                                xgc->data.set_mat.gwydd_nvx[xgc->data.set_mat.gwydds->len - 1] = 0;

                                for (i = 0; i < xres; i++) { //x
                                    for (j = 0; j < yres; j++) { //y
                                        for (k = 0; k < zres; k++) { //z
                                            mask = 0;

                                            if (gx.i != -1 && ggwy_data_field_inside(gx.dfield, (gint)(j - gx.yoffset / xgc->data.set.sp.dy), (gint)(k - gx.zoffset / xgc->data.set.sp.dz))) {
                                                val = gx.i + (gdouble)xres / xreal*(gwy_data_field_get_dval_real(gx.dfield, (gdouble)j*xgc->data.set.sp.dy - gx.yoffset, (gdouble)k*xgc->data.set.sp.dz - gx.zoffset, GWY_INTERPOLATION_BILINEAR) - gx.xoffset);
                                                if (gx.mask == 1 && gx.mfield)
                                                    mask = gwy_data_field_get_dval_real(gx.mfield, (gdouble)j*xgc->data.set.sp.dy - gx.xoffset, (gdouble)k*xgc->data.set.sp.dz - gx.yoffset, GWY_INTERPOLATION_ROUND);
                                                else if (gx.mask == -1 && gx.mfield)
                                                    mask = fabs(gwy_data_field_get_dval_real(gx.mfield, (gdouble)j*xgc->data.set.sp.dy - gx.xoffset, (gdouble)k*xgc->data.set.sp.dz - gx.yoffset, GWY_INTERPOLATION_ROUND) - 1);

                                                if (gx.depth > 0 && mask == 0 && i > val && i < (val + gx.depth)) {
                                                    buf->data[i][j][k] = (gdouble)i / (gdouble)xres;
                                                    if (xgc->data.set_mat.gwydd_xmin[xgc->data.set_mat.gwydds->len - 1] > i)
                                                        xgc->data.set_mat.gwydd_xmin[xgc->data.set_mat.gwydds->len - 1] = i;
                                                    if (xgc->data.set_mat.gwydd_xmax[xgc->data.set_mat.gwydds->len - 1] < i)
                                                        xgc->data.set_mat.gwydd_xmax[xgc->data.set_mat.gwydds->len - 1] = i;
                                                    xgc->data.set_mat.gwydd_nvx[xgc->data.set_mat.gwydds->len - 1]++;
                                                }
                                            }

                                            if (gx.j != -1 && ggwy_data_field_inside(gx.dfield, (gint)(i - gx.xoffset / xgc->data.set.sp.dx), (gint)(k - gx.zoffset / xgc->data.set.sp.dz))) {
                                                val = gx.j + (gdouble)yres / yreal*(gwy_data_field_get_dval_real(gx.dfield, (gdouble)i*xgc->data.set.sp.dx - gx.xoffset, (gdouble)k*xgc->data.set.sp.dz - gx.zoffset, GWY_INTERPOLATION_BILINEAR) - gx.yoffset);
                                                if (gx.mask == 1 && gx.mfield)
                                                    mask = gwy_data_field_get_dval_real(gx.mfield, (gdouble)i*xgc->data.set.sp.dx - gx.xoffset, (gdouble)k*xgc->data.set.sp.dz - gx.yoffset, GWY_INTERPOLATION_ROUND);
                                                else if (gx.mask == -1 && gx.mfield)
                                                    mask = fabs(gwy_data_field_get_dval_real(gx.mfield, (gdouble)i*xgc->data.set.sp.dx - gx.xoffset, (gdouble)k*xgc->data.set.sp.dz - gx.yoffset, GWY_INTERPOLATION_ROUND) - 1);


                                                if (gx.depth > 0 && mask == 0 && j > val && j < (val + gx.depth)) {
                                                    buf->data[i][j][k] = (gdouble)j / (gdouble)yres;
                                                    if (xgc->data.set_mat.gwydd_ymin[xgc->data.set_mat.gwydds->len - 1] > k)
                                                        xgc->data.set_mat.gwydd_ymin[xgc->data.set_mat.gwydds->len - 1] = j;
                                                    if (xgc->data.set_mat.gwydd_ymax[xgc->data.set_mat.gwydds->len - 1] < k)
                                                        xgc->data.set_mat.gwydd_ymax[xgc->data.set_mat.gwydds->len - 1] = j;
                                                    xgc->data.set_mat.gwydd_nvx[xgc->data.set_mat.gwydds->len - 1]++;
                                                }
                                            }

                                            if (gx.k != -1 && ggwy_data_field_inside(gx.dfield, (gint)(i - gx.xoffset / xgc->data.set.sp.dx), (gint)(j - gx.yoffset / xgc->data.set.sp.dy))) {
                                                val = gx.k + (gdouble)zres / zreal*(gwy_data_field_get_dval_real(gx.dfield, (gdouble)i*xgc->data.set.sp.dx - gx.xoffset, (gdouble)j*xgc->data.set.sp.dy - gx.yoffset, GWY_INTERPOLATION_BILINEAR) - gx.zoffset);
                                                if (gx.mask == 1 && gx.mfield)
                                                    mask = gwy_data_field_get_dval_real(gx.mfield, (gdouble)i*xgc->data.set.sp.dx - gx.xoffset, (gdouble)j*xgc->data.set.sp.dy - gx.yoffset, GWY_INTERPOLATION_ROUND);
                                                else if (gx.mask == -1 && gx.mfield)
                                                    mask = fabs(gwy_data_field_get_dval_real(gx.mfield, (gdouble)i*xgc->data.set.sp.dx - gx.xoffset, (gdouble)j*xgc->data.set.sp.dy - gx.yoffset, GWY_INTERPOLATION_ROUND) - 1);

                                                //                                                printf("%d %g %g\n", k, val, (val+gx.depth));
                                                if (gx.depth > 0 && k > val && k < (val + gx.depth)) {
                                                    buf->data[i][j][k] = (gdouble)k / (gdouble)zres;
                                                    if (xgc->data.set_mat.gwydd_zmin[xgc->data.set_mat.gwydds->len - 1] > k) xgc->data.set_mat.gwydd_zmin[xgc->data.set_mat.gwydds->len - 1] = k;
                                                    if (xgc->data.set_mat.gwydd_zmax[xgc->data.set_mat.gwydds->len - 1] < k) xgc->data.set_mat.gwydd_zmax[xgc->data.set_mat.gwydds->len - 1] = k;
                                                    xgc->data.set_mat.gwydd_nvx[xgc->data.set_mat.gwydds->len - 1]++;
                                                }
                                            }
                                        }
                                    }
                                }
                                /*simplify*/

                                for (i = 1; i < xres - 1; i++) { //x
                                    for (j = 1; j < yres - 1; j++) { //y
                                        for (k = 1; k < zres - 1; k++) { //z
                                            if (((gx.i != -1) && (((j / 2) % 10) == 0 || ((k / 2) % 10) == 0))) {
                                                if (buf->data[i][j][k] > 0 &&
                                                    !(buf->data[i - 1][j][k] > 0 && buf->data[i + 1][j][k] > 0
                                                      && buf->data[i][j - 1][k] > 0 && buf->data[i][j + 1][k] > 0
                                                      && buf->data[i][j][k - 1] > 0 && buf->data[i][j][k + 1] > 0))
                                                    xgc->data.set_mat.gwyddata->data[i][j][k] = buf->data[i][j][k];
                                            } else if (((gx.j != -1) && (((i / 2) % 10) == 0 || ((k / 2) % 10) == 0))) {
                                                if (buf->data[i][j][k] > 0 &&
                                                    !(buf->data[i - 1][j][k] > 0 && buf->data[i + 1][j][k] > 0
                                                      && buf->data[i][j - 1][k] > 0 && buf->data[i][j + 1][k] > 0
                                                      && buf->data[i][j][k - 1] > 0 && buf->data[i][j][k + 1] > 0))
                                                    xgc->data.set_mat.gwyddata->data[i][j][k] = buf->data[i][j][k];
                                            } else if (((gx.k != -1) && (((i / 2) % 10) == 0 || ((j / 2) % 10) == 0))) {
                                                if (buf->data[i][j][k] > 0 &&
                                                    !(buf->data[i - 1][j][k] > 0 && buf->data[i + 1][j][k] > 0
                                                      && buf->data[i][j - 1][k] > 0 && buf->data[i][j + 1][k] > 0
                                                      && buf->data[i][j][k - 1] > 0 && buf->data[i][j][k + 1] > 0))
                                                    xgc->data.set_mat.gwyddata->data[i][j][k] = buf->data[i][j][k];
                                            }
                                        }
                                    }
                                }
                                sv_dcube_free(buf);
                            }
                            break;


                        default:
                            fprintf(stderr, "Error: unknown object type in material file (%d)\n", type);
                            fclose(fr);
                            fr = NULL;
                            return;
                            break;
                    }
                }

                fclose(fr);
                fr = NULL;

                gtk_tree_store_clear(xgc->ts_mat);
                gtk_tree_store_append(xgc->ts_mat, &iter, NULL);
                if (xgc->data.set.sm.in_vector_filename)
                    g_snprintf(buff, sizeof(buff), "File: %s", xgc->data.set.sm.in_vector_filename);
                else
                    g_snprintf(buff, sizeof(buff), "File: undefined");


                gtk_tree_store_set(xgc->ts_mat, &iter, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, -1, COLUMN_SHOW_TOGGLE, FALSE, -1);


                for (m = 0; m < (gint)MIN(xgc ->data.set_mat.spheres->len, TREE_MAXENT); m++) {
                    sx = g_array_index(xgc->data.set_mat.spheres, SvSphere, m);

                    gtk_tree_store_append(xgc->ts_mat, &child, &iter);
                    g_snprintf(buff, sizeof(buff), "Sphere radius %g at %g, %g, %g", sx.radius, sx.pnt1[0], sx.pnt1[1], sx.pnt1[2]);
                    gtk_tree_store_set(xgc->ts_mat, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MAT_SPHERE + m, COLUMN_SHOW_TOGGLE, TRUE, -1);
                }
                for (m = 0; m < (gint)MIN(xgc->data.set_mat.voxels->len, TREE_MAXENT); m++) {
                    vx = g_array_index(xgc->data.set_mat.voxels, SvVoxel, m);

                    gtk_tree_store_append(xgc->ts_mat, &child, &iter);
                    g_snprintf(buff, sizeof(buff), "Box at (%g %g %g), (%g %g %g)", vx.pnt1[0], vx.pnt1[1], vx.pnt1[2], vx.pnt2[0], vx.pnt2[1], vx.pnt2[2]);
                    gtk_tree_store_set(xgc->ts_mat, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MAT_VOXEL + m, COLUMN_SHOW_TOGGLE, TRUE, -1);
                }
                for (m = 0; m < (gint)MIN(xgc->data.set_mat.cylinders->len, TREE_MAXENT); m++) {
                    cx = g_array_index(xgc->data.set_mat.cylinders, SvCylinder, m);

                    gtk_tree_store_append(xgc->ts_mat, &child, &iter);
                    g_snprintf(buff, sizeof(buff), "Cylinder radius %g at  (%g %g %g), (%g %g %g)", cx.radius, cx.pnt1[0], cx.pnt1[1], cx.pnt1[2], cx.pnt2[0], cx.pnt2[1], cx.pnt2[2]);
                    gtk_tree_store_set(xgc->ts_mat, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MAT_CYLINDER + m, COLUMN_SHOW_TOGGLE, TRUE, -1);
                }
                for (m = 0; m < (gint)MIN(xgc->data.set_mat.cones->len, TREE_MAXENT); m++) {
                    cnx = g_array_index(xgc->data.set_mat.cones, SvCone, m);

                    gtk_tree_store_append(xgc->ts_mat, &child, &iter);
                    g_snprintf(buff, sizeof(buff), "Cone radius %g at (%g %g %g), (%g %g %g)", cnx.radius, cnx.pnt1[0], cnx.pnt1[1], cnx.pnt1[2], cnx.pnt2[0], cnx.pnt2[1], cnx.pnt2[2]);
                    gtk_tree_store_set(xgc->ts_mat, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MAT_CONE + m, COLUMN_SHOW_TOGGLE, TRUE, -1);
                }
                for (m = 0; m < (gint)MIN(xgc->data.set_mat.rcones->len, TREE_MAXENT); m++) {
                    rcnx = g_array_index(xgc->data.set_mat.rcones, SvRCone, m);

                    gtk_tree_store_append(xgc->ts_mat, &child, &iter);
                    g_snprintf(buff, sizeof(buff), "Cut cone %g...%g at (%g %g %g), (%g %g %g)", rcnx.radius1, rcnx.radius2, rcnx.pnt1[0], rcnx.pnt1[1], rcnx.pnt1[2], rcnx.pnt2[0], rcnx.pnt2[1], rcnx.pnt2[2]);
                    gtk_tree_store_set(xgc->ts_mat, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MAT_RCONE + m, COLUMN_SHOW_TOGGLE, TRUE, -1);
                }
                for (m = 0; m < (gint)MIN(xgc->data.set_mat.tetrahedrons->len, 100 * TREE_MAXENT); m++) { //FIXME, this should count only listed tetrahedrons, not tetgens!
                    tx = g_array_index(xgc->data.set_mat.tetrahedrons, SvTetrahedron, m);
                    isobject = 0;
                    for (i = 0; i < xgc->data.set_mat.ntetgens; i++) {
                        if (tx.n >= xgc->data.set_mat.tetgen_start[i] && tx.n < xgc->data.set_mat.tetgen_end[i])
                            isobject = 1;
                    }
                    if (!isobject) {
                        gtk_tree_store_append(xgc->ts_mat, &child, &iter);
                        g_snprintf(buff, sizeof(buff), "Tetrahedron starting at %g, %g, %g", tx.pnt1[0], tx.pnt1[1], tx.pnt1[2]);
                        gtk_tree_store_set(xgc->ts_mat, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MAT_TETRAHEDRON + m, COLUMN_SHOW_TOGGLE, TRUE, -1);
                    }
                }
                for (m = 0; m < MIN(xgc->data.set_mat.ntetgens, TREE_MAXENT); m++) {
                    gtk_tree_store_append(xgc->ts_mat, &child, &iter);
                    g_snprintf(buff, sizeof(buff), "mesh (%s) at %g, %g, %g", xgc->data.set_mat.tetgen_filebase[m], xgc->data.set_mat.tetgen_xshift[m], xgc->data.set_mat.tetgen_yshift[m], xgc->data.set_mat.tetgen_zshift[m]);
                    gtk_tree_store_set(xgc->ts_mat, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MAT_MESH + m, COLUMN_SHOW_TOGGLE, TRUE, -1);
                }
                for (m = 0; m < (gint)MIN(xgc->data.set_mat.gwydds->len, TREE_MAXENT); m++) {
                    gx = g_array_index(xgc->data.set_mat.gwydds, SvGwydd, m);

                    gtk_tree_store_append(xgc->ts_mat, &child, &iter);

                    if (xgc->data.set_mat.gwydd_nvx[xgc->data.set_mat.gwydds->len - 1] == 0)
                        g_snprintf(buff, sizeof(buff), "Gwyddion field missed: pos (%d, %d, %d), span (%d %d %d) to (%d %d %d) vx", gx.i, gx.j, gx.k,
                        (gint)xgc->data.set_mat.gwydd_xmin[m], (gint)xgc->data.set_mat.gwydd_ymin[m], (gint)xgc->data.set_mat.gwydd_zmin[m],
                                   (gint)xgc->data.set_mat.gwydd_xmax[m], (gint)xgc->data.set_mat.gwydd_ymax[m], (gint)xgc->data.set_mat.gwydd_zmax[m]);
                    else
                        g_snprintf(buff, sizeof(buff), "Gwyddion field at %d, %d, %d, spanning (%d %d %d) to (%d %d %d) vx", gx.i, gx.j, gx.k,
                        (gint)xgc->data.set_mat.gwydd_xmin[m], (gint)xgc->data.set_mat.gwydd_ymin[m], (gint)xgc->data.set_mat.gwydd_zmin[m],
                                   (gint)xgc->data.set_mat.gwydd_xmax[m], (gint)xgc->data.set_mat.gwydd_ymax[m], (gint)xgc->data.set_mat.gwydd_zmax[m]);
                    gtk_tree_store_set(xgc->ts_mat, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MAT_GWYDD + m, COLUMN_SHOW_TOGGLE, TRUE, -1);
                }

                g_snprintf(buff, sizeof(buff), "Ready");
                gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);

                gtk_widget_queue_draw(xgc->view_scene);
            }
        }

        xgc->data.set_mat.ngtot = ngtot;
    }

exit:
    if (NULL != filename)
        g_free(filename);
    if (NULL != fr)
        fclose(fr);    

    return;
} /* mat_buffer_changed */


static void
cb_child_watch(GPid pid, gint status, XGControls *xgc)
{
    /* Close pid */
    gchar buff[256];    

    if (xgc->timeout_func)
        g_source_remove(xgc->timeout_func);
    if (xgc->watch)
        g_source_remove(xgc->watch);
    if (xgc->data.pid)
        g_spawn_close_pid(xgc->data.pid);

    xgc->watch = 0;

#ifdef G_OS_WIN32
    if (xgc->usekiller && g_file_test("terminate.txt", G_FILE_TEST_EXISTS))
        g_remove("terminate.txt");
#endif

    while (gtk_events_pending()) {
        gtk_main_iteration();
    }

#ifdef G_OS_WIN32
    xgc->data.pid = NULL;
#endif

    xgc->timeout_func = 0;

#ifdef G_OS_WIN32
    CloseHandle(xgc->data.pid);
#endif

    set_computation_controls_sensitive(xgc);

    update_graph(xgc);    

    g_snprintf(buff, sizeof(buff), "GSvit successfully finished.");
    gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
    while (gtk_events_pending()) {
        gtk_main_iteration();
    }
}

static gboolean
cb_out_watch(GIOChannel *channel, GIOCondition cond, XGControls *xgc)
{
    gchar *string, buff[50];
    gsize  size;
    GtkTextIter end;

    //printf("watch! %d\n", xgc->watch);

    if (xgc->watch == 0)
        return TRUE;


    if (cond == G_IO_HUP) {
        g_io_channel_unref(channel);
        return FALSE;
    }

    g_io_channel_read_line(channel, &string, &size, NULL, NULL);
    if (string == NULL || g_utf8_strlen(string, 100) <= 0) {
        //g_snprintf(buff, sizeof(buff), "No string to show in output.");
        //gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff); 
    } else {
        //gtk_text_buffer_insert_at_cursor(xgc->tb_out, string, -1);

        gtk_text_buffer_get_end_iter(xgc->tb_out, &end);
        gtk_text_buffer_insert(xgc->tb_out, &end, string, -1);

        g_free(string);
        gtk_text_buffer_get_end_iter(xgc->tb_out, &end);
        gtk_text_buffer_place_cursor(xgc->tb_out, &end);
        gtk_text_view_scroll_to_iter(GTK_TEXT_VIEW(xgc->tv_out), &end, 0, 0, 0, 0);

        g_snprintf(buff, sizeof(buff), "Running GSvit. Computing..");
        gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
    }
    while (gtk_events_pending() && xgc->watch != 0) {
        gtk_main_iteration();
    }
    return TRUE;
}


gboolean update_image(gpointer data)
{
    XGControls *xgc = (XGControls *)data;

    gtk_widget_queue_draw(xgc->view_scene);

    //printf("update image from %s\n", xgc->data.set.so.outfile);
    return TRUE;
}



gboolean update_graph(gpointer data)
{
    XGControls *xgc = (XGControls *)data;
    gint pos;
    gint step;
    gfloat ex, ey, ez, hx, hy, hz;
    gdouble x[10000];
    gdouble y[10000];
    gint n = 0;
    gchar buff[256];
    FILE *fr;

#ifndef G_OS_WIN32
    gint status;
    pid_t return_pid;
#else
    gint res = 0;
#endif

    // printf("update graph\n");
    if (xgc->data.pid) {

#ifndef G_OS_WIN32
        return_pid = waitpid(xgc->data.pid, &status, WNOHANG);
        if (return_pid == -1 || return_pid == xgc->data.pid) {
#else
        GetExitCodeProcess(xgc->data.pid, &res);
        if (res != STILL_ACTIVE) {
#endif	
            //printf("Error checking status or GSvit finished\n");
            if (xgc->timeout_func)
                g_source_remove(xgc->timeout_func);
            if (xgc->watch)
                g_source_remove(xgc->watch);
            if (xgc->data.pid)
                g_spawn_close_pid(xgc->data.pid);

            xgc->watch = 0;
#ifdef G_OS_WIN32
            xgc->data.pid = NULL;
#endif
            xgc->timeout_func = 0;

            gtk_widget_set_sensitive(GTK_WIDGET(xgc->run_tool_item), TRUE);
            gtk_widget_set_sensitive(GTK_WIDGET(xgc->run_menu_item), TRUE);
        }
    }

    pos = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->graph_combo));

    if (pos < 0) {
        //printf("Unsuitable selection of graph in combo\n");
        return TRUE;
    }


    //if (pos>=0) {
    //   printf("updating %d [%s]\n", pos, xgc->filestoshow[pos]);
    //} else printf("nothing to update\n");

    fr = fopen(xgc->filestoshow[pos], "r");
    if (fr) {
        if (xgc->formattoshow[pos] == SV_COMP_ALL) {
            n = 0;
            while (fscanf(fr, "%d", &step) != EOF) {
                fscanf(fr, "%f", &ex);
                fscanf(fr, "%f", &ey);
                fscanf(fr, "%f", &ez);
                fscanf(fr, "%f", &hx);
                fscanf(fr, "%f", &hy);
                fscanf(fr, "%f", &hz);
                x[n] = step;
                y[n] = (gdouble)ex*(gdouble)ex + (gdouble)ey*(gdouble)ey + (gdouble)ez*(gdouble)ez;
                n++;
            }
            fclose(fr);
        } else if (xgc->formattoshow[pos] == SV_COMP_EX || xgc->formattoshow[pos] == SV_COMP_EY || xgc->formattoshow[pos] == SV_COMP_EZ || xgc->formattoshow[pos] == SV_COMP_HX || xgc->formattoshow[pos] == SV_COMP_HY || xgc->formattoshow[pos] == SV_COMP_HZ) {
            n = 0;
            while (fscanf(fr, "%d", &step) != EOF) {
                fscanf(fr, "%f", &ex);
                x[n] = step;
                y[n] = ex;
                n++;
            }
            fclose(fr);
        }
        if (xgc->formattoshow[pos] == 99) {
            n = 0;
            while (fscanf(fr, "%d", &step) != EOF) {
                fscanf(fr, "%f", &ex);
                fscanf(fr, "%f", &ey);
                fscanf(fr, "%f", &ez);
                x[n] = step;
                y[n] = (gdouble)ex*(gdouble)ex + (gdouble)ey*(gdouble)ey + (gdouble)ez*(gdouble)ez;
                n++;
            }
            fclose(fr);
        }
        if (xgc->formattoshow[pos] == 100) {
            n = 0;
            while (fscanf(fr, "%d", &step) != EOF) {
                fscanf(fr, "%f", &ex);
                fscanf(fr, "%f", &ey);
                fscanf(fr, "%f", &ez);
                fscanf(fr, "%f", &ex);
                fscanf(fr, "%f", &ey);
                fscanf(fr, "%f", &ez);
                x[n] = step;
                y[n] = sqrt((gdouble)ex*(gdouble)ex + (gdouble)ey*(gdouble)ey + (gdouble)ez*(gdouble)ez);
                n++;
            }
            fclose(fr);
        }
        g_object_set(xgc->gcmodel1, "description", xgc->filestoshow[pos], "mode", GWY_GRAPH_CURVE_LINE, NULL);

        if (xgc->gcmodel1 && n > 0)
            gwy_graph_curve_model_set_data(xgc->gcmodel1, x, y, n);

        g_snprintf(buff, sizeof(buff), "Computing, updating graph %s", xgc->filestoshow[pos]);
        gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);

        while (gtk_events_pending()) {
            gtk_main_iteration();
        }
    } else {
        g_snprintf(buff, sizeof(buff), "Error: cannot update graph %s", xgc->filestoshow[pos]);
        gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
    }

    return TRUE;
}
#define MAGIC "GWYO"
#define MAGIC2 "GWYP"
#define MAGIC_SIZE (sizeof(MAGIC)-1)

void outfile_monitor_cb(GFileMonitor *monitor, GFile *file,
                        GFile *other_file,
                        GFileMonitorEvent event_type,
                        XGControls *xgc)
{
    gsize size = 0;
    gsize pos = 0;
    guchar *buffer = NULL;
    gchar *sbuf = NULL;
    GError *err = NULL;
    GObject *object = NULL;
    gchar data_key[30];
    gchar matchkey[100];
    gint i, k, *ids, last;


    if (event_type == G_FILE_MONITOR_EVENT_CHANGES_DONE_HINT) {
        //printf("Loading Gwyddion file %s\n", xgc->data.set.so.outfile);

        if (!g_file_get_contents(xgc->data.set.so.outfile, &sbuf, &size, &err)) {
            fprintf(stderr, "Error: cannot load file\n");
            for (i = 0; i < 100; i++)
                xgc->outputfield[i] = NULL;
            return;
        }
        buffer = (guchar *)sbuf;
        if (size < MAGIC_SIZE || (memcmp(buffer, MAGIC, MAGIC_SIZE) && memcmp(buffer, MAGIC2, MAGIC_SIZE))) {
            fprintf(stderr, "Error: wrong file type.\n");
            g_free(buffer);
            for (i = 0; i < 100; i++)
                xgc->outputfield[i] = NULL;
            return;
        }
        if (!memcmp(buffer, MAGIC, MAGIC_SIZE)) {
            fprintf(stderr, "Error: File in old format, now unsupported.\n");
            g_free(buffer);
            for (i = 0; i < 100; i++)
                xgc->outputfield[i] = NULL;
            return;
        } else
            object = gwy_serializable_deserialize(buffer + MAGIC_SIZE, size - MAGIC_SIZE, &pos);

        g_free(buffer);


        if (!object) {
            fprintf(stderr, "Error: Deserialization failed.\n");
            for (i = 0; i < 100; i++)
                xgc->outputfield[i] = NULL;
            return;
        }
        if (!GWY_IS_CONTAINER(object)) {
            fprintf(stderr, "Error: Deserialization resulted in something unexpected.\n");
            gwy_object_unref(object);
            for (i = 0; i < 100; i++)
                xgc->outputfield[i] = NULL;
            return;
        }


        gwy_object_unref(xgc->container);
        xgc->container = GWY_CONTAINER(object);

        gwy_app_data_browser_add(xgc->container);

        xgc->image_selected = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->image_combo));

        for (k = 0; k < xgc->data.set.so.nimgs; k++) {
            /*if (xgc->image_selected>=0) {*/
            g_snprintf(matchkey, sizeof(matchkey), "%s*", xgc->data.set.so.imgs[k].filebase);
            ids = gwy_app_data_browser_find_data_by_title(xgc->container,
                                                          matchkey);



            //if (ids[0]==-1) printf("nothing matches %s\n", matchkey);
            i = 0;
            last = -1;
            while (ids[i] != -1) {
                //printf("Matching id: %d for %s\n", ids[i], matchkey);
                last = ids[i];
                i++;
            };
            //                    printf("key %d will be used for id key %s\n", last, matchkey);

            if (last >= 0) {
                g_snprintf(data_key, sizeof(data_key), "/%i/data", last);

                xgc->outputfield[k] = gwy_container_get_object(xgc->container, g_quark_from_string(data_key));

                /*redraw 3d view*/
                gtk_widget_queue_draw(xgc->view_scene);
            }

        }
        gwy_app_data_browser_remove(xgc->container);
    }
}

void check_files_consistency_cb(GtkWidget *widget, XGControls *xgc)
{

    GtkWidget *dialog, *label, *content_area;
    GString *message;
    gboolean ok;

    message = g_string_new("<b>Consistency check:</b>\n");

    message = sv_check(&(xgc->data.set), message, &ok);

    dialog = gtk_dialog_new_with_buttons("Consistency check",
                                         GTK_WINDOW(xgc->toplevel),
                                         GTK_DIALOG_DESTROY_WITH_PARENT,
                                         GTK_STOCK_OK,
                                         GTK_RESPONSE_NONE,
                                         NULL);

    gtk_window_set_default_size(GTK_WINDOW(dialog), 200, 100);

    content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    label = gtk_label_new(NULL);
    gtk_label_set_markup(GTK_LABEL(label), message->str);
    /* Ensure that the dialog box is destroyed when the user responds. */
    g_signal_connect_swapped(dialog, "response", G_CALLBACK(gtk_widget_destroy), dialog);
    /* Add the label, and show everything we've added to the dialog. */
    gtk_container_add(GTK_CONTAINER(content_area), label);
    gtk_widget_show_all(dialog);
}


void gsvit_run_cb(GtkWidget *widget, XGControls *xgc)
{
    GPid pid;

    pt_remeber_expanded_rows_and_selection(xgc);
    mt_remeber_expanded_rows_and_selection(xgc);
    write_parfile(xgc);
    pt_restore_expanded_rows_and_selection(xgc);
    mt_restore_expanded_rows_and_selection(xgc);

#ifndef G_OS_WIN32
    gchar *argv[] = {g_strdup(xgc->data.gsvit_location), g_strdup(xgc->data.parfilename), NULL};
#endif

#ifdef G_OS_WIN32
    FILE *fw;
    gchar *argv[] = {g_strdup(xgc->data.gsvit_location), g_strdup(xgc->data.parfilename), "--killer", "terminate.txt", NULL};
#endif

    gint out;
    /* gint  err; */
    GIOChannel *out_ch;
    /* GIOChannel *err_ch; */
    gboolean   ret;
    GtkWidget *dialog;
    GError *error = NULL;

    /* Spawn child process */

    if (!xgc->data.gsvit_location) {
        dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                                        GTK_DIALOG_DESTROY_WITH_PARENT,
                                        GTK_MESSAGE_ERROR,
                                        GTK_BUTTONS_CLOSE,
                                        "GSvit location not specified.\nPlease check GSvit executable location in File/Preferences dialog.");
        gtk_dialog_run(GTK_DIALOG(dialog));
        gtk_widget_destroy(dialog);

        return;
    }
    if (!xgc->data.parfilename) {
        dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                                        GTK_DIALOG_DESTROY_WITH_PARENT,
                                        GTK_MESSAGE_ERROR,
                                        GTK_BUTTONS_CLOSE,
                                        "Parameter file not specified.\nPlease save parameter file.");
        gtk_dialog_run(GTK_DIALOG(dialog));
        gtk_widget_destroy(dialog);

        return;
    }


    //printf("%s     %s %s\n", g_get_current_dir(), xgc->data.gsvit_location, xgc->data.parfilename);

    file_save(xgc);

#ifdef G_OS_WIN32
    //printf("writing to dir %s %s\n", g_get_current_dir(),  xgc->data.parfilename);
    if (fw = fopen("terminate.txt", "w")) {
        //if (fw = fopen(argv[3], "w")) {
        fprintf(fw, "0\n");
        fclose(fw);
        xgc->usekiller = TRUE;
    } else
        xgc->usekiller = FALSE;
#endif

    /* do not redirect stderr until we know why it breaks cow on Linux */
    ret = g_spawn_async_with_pipes(NULL, argv, NULL,
                                   G_SPAWN_DO_NOT_REAP_CHILD, NULL,
                                   NULL, &pid, NULL, &out, NULL, &error);

    if (ret == FALSE && error != NULL) {
        dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                                        GTK_DIALOG_DESTROY_WITH_PARENT,
                                        GTK_MESSAGE_ERROR,
                                        GTK_BUTTONS_CLOSE,
                                        "Running GSvit failed: %s\nPlease check GSvit executable location in File/Preferences dialog.",
                                        g_strdup(error->message));
        gtk_dialog_run(GTK_DIALOG(dialog));
        gtk_widget_destroy(dialog);

        return;
    }

    xgc->data.pid = pid;
    xgc->timeout_func = g_timeout_add_full(G_PRIORITY_HIGH, 1000, update_graph, xgc, NULL);

    /* Add watch function to catch termination of the process. This function
     * will clean any remnants of process. */
    g_child_watch_add(pid, (GChildWatchFunc)cb_child_watch, xgc);

    /* Create channels that will be used to read data from pipes. */
#ifdef G_OS_WIN32
    out_ch = g_io_channel_win32_new_fd(out);
    /* err_ch = g_io_channel_win32_new_fd(err); */
#else
    out_ch = g_io_channel_unix_new(out);
    /* err_ch = g_io_channel_unix_new(err); */
#endif

    /* Add watches to channels */
    xgc->watch = g_io_add_watch(out_ch, G_IO_IN | G_IO_HUP, (GIOFunc)cb_out_watch, xgc);

    /* /\* also watch the error channel, mit it into the output channel display *\/ */
    /* g_io_add_watch(err_ch, G_IO_IN | G_IO_HUP, (GIOFunc)cb_out_watch, xgc); */

    xgc->outputfile_monitor = g_file_monitor_file(g_file_new_for_path(xgc->data.set.so.outfile),
                                                  G_FILE_MONITOR_NONE,
                                                  NULL,
                                                  &error);

    if (error)
        printf("Error while connecting monitor to output file: %s\n", error->message);
    else
        g_signal_connect(xgc->outputfile_monitor, "changed", G_CALLBACK(outfile_monitor_cb), xgc);

    set_computation_controls_sensitive(xgc);
}


void gsvit_stop_cb(GtkWidget *widget, XGControls *xgc)
{    
    if ( (xgc == NULL) || (xgc->data.pid == 0) )
        return;

#ifdef G_OS_WIN32
    FILE *fw;
    if (xgc->usekiller) { //stupid workaround to prevent termination troubles.
        if (fw = fopen("terminate.txt", "w")) {
            fprintf(fw, "1\n");
            fclose(fw);
        } else
            TerminateProcess(xgc->data.pid, 0);
    } else
        TerminateProcess(xgc->data.pid, 0);    
#else
    kill(xgc->data.pid, SIGINT);
#endif
}

void gwydd_cb(GtkWidget *widget, XGControls *xgc)
{
    GPid pid;
    gchar *argv[] = {g_strdup(xgc->data.gwyddion_location), g_strdup(xgc->data.set.so.outfile), NULL};
    gchar buff[256];
    GError *error = NULL;
    GtkWidget *dialog;

    if (!xgc->data.gwyddion_location) {
        dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                                        GTK_DIALOG_DESTROY_WITH_PARENT,
                                        GTK_MESSAGE_ERROR,
                                        GTK_BUTTONS_CLOSE,
                                        "Gwyddion location not specified (see Preferences)");
        gtk_dialog_run(GTK_DIALOG(dialog));
        gtk_widget_destroy(dialog);

        return;
    }


    if (xgc->data.set.so.outfile) {
        if (xgc->data.gwyddion_location) {
            /* Spawn child process */
            if (!g_spawn_async(NULL, argv, NULL, G_SPAWN_DO_NOT_REAP_CHILD, NULL, NULL, &pid, &error)) {
                g_snprintf(buff, sizeof(buff), "Error starting Gwyddion: %s", error->message);
                gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
            } else {
                g_snprintf(buff, sizeof(buff), "Calling gwyddion on file %s.", xgc->data.set.so.outfile);
                gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
            }

        }
    } else {
        g_snprintf(buff, sizeof(buff), "Error starting Gwyddion: cannot find an output file to show.");
        gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
    }
}

static void
pars_page_switched_cb(XGControls *xgc, G_GNUC_UNUSED GtkNotebookPage *page, gint pagenum)
{
    GtkWidget *dialog;

    if (pagenum == 0) {     /* Parameters tab */
        if (xgc->par_file_success == FALSE) {
            if (xgc->par_tree_success == FALSE) {
                dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                                                GTK_DIALOG_DESTROY_WITH_PARENT,
                                                GTK_MESSAGE_WARNING,
                                                GTK_BUTTONS_OK,
                                                MSG_MB_ERROR_PARSING_PARFILE);
                gtk_dialog_run(GTK_DIALOG(dialog));
                gtk_widget_destroy(dialog);
            } else {
                dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                                                GTK_DIALOG_DESTROY_WITH_PARENT,
                                                GTK_MESSAGE_QUESTION,
                                                GTK_BUTTONS_YES_NO,
                                                "Error parsing parameter file. Reload settings?");
                if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_YES) {
                    /* reload settings from tree (data) */

                    write_parfile(xgc);
                    //pt_create_tree(xgc);
                    pt_restore_expanded_rows_and_selection(xgc);

                    //gtk_notebook_set_current_page(GTK_NOTEBOOK(xgc->param_notebook), 1);
                } else {                
                    /* create tree from error parfile */

                    clear_settings(&(xgc->data.set), FALSE);
                    pt_create_tree(xgc);
                    pt_restore_expanded_rows_and_selection(xgc);
                }
                gtk_widget_destroy(dialog);                
            }            
        } else {            
            pt_create_tree(xgc);
            pt_restore_expanded_rows_and_selection(xgc);
        }
    } else {                /* Parameter file tab */
        pt_remeber_expanded_rows_and_selection(xgc);
        if(xgc->par_tree_success)
            write_parfile(xgc);
    }
}

static void
material_page_switched_cb(XGControls *xgc, G_GNUC_UNUSED GtkNotebookPage *page, gint pagenum)
{
    GtkWidget *dialog;

    if (pagenum == 0) {     /* Parameters tab */
        if (xgc->mat_file_success == FALSE) {
            if (xgc->mat_tree_success == FALSE) {
                dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                    GTK_DIALOG_DESTROY_WITH_PARENT,
                    GTK_MESSAGE_WARNING,
                    GTK_BUTTONS_OK,
                    MSG_MB_ERROR_PARSING_PARFILE);
                gtk_dialog_run(GTK_DIALOG(dialog));
                gtk_widget_destroy(dialog);
            } else {
                dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                    GTK_DIALOG_DESTROY_WITH_PARENT,
                    GTK_MESSAGE_QUESTION,
                    GTK_BUTTONS_YES_NO,
                    "Error parsing material file. Reload settings?");
                if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_YES) {
                    /* reload settings from tree (data) */

                    write_matfile(xgc);
                    //pt_create_tree(xgc);
                    mt_restore_expanded_rows_and_selection(xgc);

                    //gtk_notebook_set_current_page(GTK_NOTEBOOK(xgc->param_notebook), 1);
                } else {
                    /* create tree from error parfile */

                    //clear_settings(&(xgc->data.set), FALSE);
                    //pt_create_tree(xgc);
                    //pt_restore_expanded_rows_and_selection(xgc);

                    alloc_set_mat(&(xgc->data.set_mat));
                    mt_create_tree(xgc);
                    mt_restore_expanded_rows_and_selection(xgc);
                }
                gtk_widget_destroy(dialog);
            }
        } else {
            mt_create_tree(xgc);
            mt_restore_expanded_rows_and_selection(xgc);
        }
    } else {                /* Material file tab */
        mt_remeber_expanded_rows_and_selection(xgc);
        if (xgc->mat_tree_success)
            write_matfile(xgc);
    }
}

static void
visible_toggle(GtkCellRendererToggle *cell, gchar *path_str, XGControls *xgc)
{
    GtkTreeModel *model = GTK_TREE_MODEL(xgc->ts_par);
    GtkTreeIter  iter;
    GtkTreePath *path = gtk_tree_path_new_from_string(path_str);
    gboolean checked;
    gchar *text;
    gint id;


    gtk_tree_model_get_iter(model, &iter, path);
    gtk_tree_model_get(model, &iter, COLUMN_PARAMETER, &text, COLUMN_CHECK, &checked, COLUMN_ID, &id, -1);

    //checked ^= 1;
    checked = !checked;
    //gtk_tree_store_set(GTK_TREE_STORE (model), &iter, COLUMN_CHECK, checked, -1);
    gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_CHECK, checked, -1);                  /* BUG: fixed */

    if (id == SET_POOL) {
        xgc->data.is_pool = checked;
    } else if (id == SET_SF)
        xgc->data.is_sf = checked;
    else if (id == SET_TSF)
        xgc->data.is_tsf = checked;
    else if (id == SET_TSFF)
        xgc->data.is_tsff = checked;
    else if (id == SET_LTSF)
        xgc->data.is_ltsf = checked;
    else if (id == SET_LTSFF)
        xgc->data.is_ltsff = checked;
    else if (id >= SET_PSOURCE && id < SET_POUT)
        xgc->data.is_psrc[id - SET_PSOURCE] = checked;
    else if (id >= SET_POUT && id < SET_IOUT)
        xgc->data.is_outpnt[id - SET_POUT] = checked;
    else if (id >= SET_IOUT && id < SET_IIOUT)
        xgc->data.is_outpln[id - SET_IOUT] = checked;
    else if (id >= SET_IIOUT && id < SET_COUT)
        xgc->data.is_outimg[id - SET_IIOUT] = checked;
    else /*if (id >= SET_COUT && id < SET_SOUT)
        xgc->data.is_outcub[id-SET_COUT] = checked; 
    else*/ if (id >= SET_SOUT && id < SET_FOUT)
        xgc->data.is_outsum[id - SET_SOUT] = checked;
    else if (id >= SET_FOUT && id < SET_GROW)
        xgc->data.is_outforce[id - SET_FOUT] = checked;
    else if (id == SET_BND) {
        if (g_strstr_len(text, -1, "x0:"))
            xgc->data.is_bx0 = checked;
        if (g_strstr_len(text, -1, "xn:"))
            xgc->data.is_bxn = checked;
        if (g_strstr_len(text, -1, "y0:"))
            xgc->data.is_by0 = checked;
        if (g_strstr_len(text, -1, "yn:"))
            xgc->data.is_byn = checked;
        if (g_strstr_len(text, -1, "z0:"))
            xgc->data.is_bz0 = checked;
        if (g_strstr_len(text, -1, "zn:"))
            xgc->data.is_bzn = checked;
        if (g_strstr_len(text, -1, "x0 periodic"))
            xgc->data.is_mbx0 = checked;
        if (g_strstr_len(text, -1, "xn periodic"))
            xgc->data.is_mbxn = checked;
        if (g_strstr_len(text, -1, "y0 periodic"))
            xgc->data.is_mby0 = checked;
        if (g_strstr_len(text, -1, "yn periodic"))
            xgc->data.is_mbyn = checked;
        if (g_strstr_len(text, -1, "z0 periodic"))
            xgc->data.is_mbz0 = checked;
        if (g_strstr_len(text, -1, "zn periodic"))
            xgc->data.is_mbzn = checked;
    } else if (id == SET_NFFF) {
        if (g_strstr_len(text, -1, "box"))
            xgc->data.is_nfff = checked;
        if (g_strstr_len(text, -1, "skip i0"))
            xgc->data.is_nfff_skipi0 = checked;
        if (g_strstr_len(text, -1, "skip in"))
            xgc->data.is_nfff_skipin = checked;
        if (g_strstr_len(text, -1, "skip j0"))
            xgc->data.is_nfff_skipj0 = checked;
        if (g_strstr_len(text, -1, "skip jn"))
            xgc->data.is_nfff_skipjn = checked;
        if (g_strstr_len(text, -1, "skip k0"))
            xgc->data.is_nfff_skipk0 = checked;
        if (g_strstr_len(text, -1, "skip kn"))
            xgc->data.is_nfff_skipkn = checked;
    } else if (id >= SET_NFFFP && id < SET_PNFAREA)
        xgc->data.is_nfff_point[id -= SET_NFFFP] = checked;
    else if (id == SET_PNFFF) {
        if (g_strstr_len(text, -1, "box"))
            xgc->data.is_pnfff = checked;
        if (g_strstr_len(text, -1, "skip k0"))
            xgc->data.is_pnfff_skipk0 = checked;
        if (g_strstr_len(text, -1, "skip kn"))
            xgc->data.is_pnfff_skipkn = checked;
    } else if (id >= SET_PNFFFP)
        xgc->data.is_pnfff_point[id -= SET_PNFFFP] = checked;


    gtk_widget_queue_draw(xgc->view_scene);

    gtk_tree_path_free(path);
}

static void
mat_toggled_cb(GtkCellRendererToggle *cell, gchar *path_str, XGControls *xgc)
{
    GtkTreeModel *model = GTK_TREE_MODEL(xgc->ts_mat);
    GtkTreeIter  iter;
    GtkTreePath *path = gtk_tree_path_new_from_string(path_str);
    gboolean checked;
    gchar *text;
    gint id, i, nmesh;

    gtk_tree_model_get_iter(model, &iter, path);
    //gtk_tree_model_get(model, &iter,  0, &text, 1, &checked, 2, &id, -1);
    gtk_tree_model_get(model, &iter, COLUMN_PARAMETER, &text, COLUMN_CHECK, &checked, COLUMN_ID, &id, -1);

    checked ^= 1;
    //gtk_tree_store_set(GTK_TREE_STORE (model), &iter, 1, checked, -1);
    gtk_tree_store_set(xgc->ts_mat, &iter, 1, checked, -1);

    if (id >= SET_MAT_SPHERE && id < SET_MAT_VOXEL)
        xgc->data.is_sphere[id - SET_MAT_SPHERE] = checked;
    else if (id >= SET_MAT_VOXEL && id < SET_MAT_CYLINDER)
        xgc->data.is_voxel[id - SET_MAT_VOXEL] = checked;
    else if (id >= SET_MAT_CYLINDER && id < SET_MAT_CONE)
        xgc->data.is_cylinder[id - SET_MAT_CYLINDER] = checked;
    else if (id >= SET_MAT_CONE && id < SET_MAT_RCONE)
        xgc->data.is_cone[id - SET_MAT_CONE] = checked;
    else if (id >= SET_MAT_RCONE && id < SET_MAT_GWYDD)
        xgc->data.is_rcone[id - SET_MAT_RCONE] = checked;
    else if (id >= SET_MAT_GWYDD && id < SET_MAT_MESH)
        xgc->data.is_gwydd[id - SET_MAT_GWYDD] = checked;
    else if (id >= SET_MAT_MESH && id < SET_MAT_TETRAHEDRON) {
        nmesh = id - SET_MAT_MESH;
        for (i = xgc->data.set_mat.tetgen_start[nmesh]; i < xgc->data.set_mat.tetgen_end[nmesh]; i++)
            xgc->data.is_tetrahedron[i] = checked;
    } else if (id >= SET_MAT_TETRAHEDRON)
        xgc->data.is_tetrahedron[id - SET_MAT_TETRAHEDRON] = checked;

    gtk_widget_queue_draw(xgc->view_scene);

    gtk_tree_path_free(path);
}


void tree_expand_row_and_select_item(XGControls *xgc, gint id_to_select)
{
    GtkTreePath *path = NULL;
    GtkTreePath *path_to_select = NULL;
    gchar *text;
    gint id, i, n_children;
    GtkTreeIter  iter, iter_parent;
    GtkTreeSelection *selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(xgc->tv_par));

    if ((id_to_select == SET_SF) || (id_to_select == SET_TSF) || (id_to_select == SET_TSFF) || (id_to_select == SET_LTSF) || (id_to_select == SET_LTSFF)
        || (id_to_select >= SET_PSOURCE && id_to_select < SET_POUT))
        path = gtk_tree_path_new_from_indices(ROW_SOURCES_INDEX, -1);           /* expand "Sources" row */
    else if (   (id_to_select >= SET_GROW && id_to_select < SET_ROUGHNESS)
             || (id_to_select >= SET_ROUGHNESS && id_to_select == SET_SPECTRAL) 
             || (id_to_select >= SET_SPECTRAL && id_to_select < SET_EXPRESSION)
             || (id_to_select >= SET_EXPRESSION && id_to_select < SET_NFAREA)
             )
        path = gtk_tree_path_new_from_indices(ROW_MEDIA_INDEX, -1);             /* expand "Media" row */
    else if (   (id_to_select >= SET_POUT && id_to_select < SET_IOUT)
             || (id_to_select >= SET_IIOUT && id_to_select < SET_COUT)
             || (id_to_select >= SET_IOUT && id_to_select < SET_IIOUT)
             || (id_to_select >= SET_COUT && id_to_select < SET_SOUT)
             || (id_to_select >= SET_SOUT && id_to_select < SET_FOUT)
             || (id_to_select >= SET_FOUT && id_to_select < SET_GROW)
             )
        path = gtk_tree_path_new_from_indices(ROW_OUTPUTS_INDEX, -1);           /* expand "Outputs" row */
    else if (   (id_to_select == SET_NFFF)
             || (id_to_select >= SET_NFFFP && id_to_select < SET_PNFAREA)      /* VII. 2. Near field to far field point */
             || (id_to_select >= SET_NFAREA && id_to_select < SET_NFFFP))     /* VII. 3. Near field to far field area */
        path = gtk_tree_path_new_from_indices(ROW_NFFF_INDEX, -1);              /* expand "NFFF" row */
    else if (   (id_to_select == SET_PNFFF)
             || (id_to_select >= SET_PNFAREA && id_to_select < SET_PNFFFP)     /* VIII. 3. Periodic near field to far field area */
             || (id_to_select >= SET_PNFFFP)) {                                /* VIII. 2. Periodic near field to far field point */
        path = gtk_tree_path_new_from_indices(ROW_PNFFF_INDEX, -1);             /* expand "PNFFF" row */
        //if (xgc->data.set.sf.nrs || xgc->data.set.sf.nareas)
        //    path = gtk_tree_path_new_from_indices(ROW_OUTPUTS_INDEX+2, -1);   /* expand "PNFFF" row */
        //else
        //    path = gtk_tree_path_new_from_indices(ROW_OUTPUTS_INDEX+1, -1);   /* expand "PNFFF" row */
    }

    if (NULL != path) {
        gtk_tree_view_expand_row(GTK_TREE_VIEW(xgc->tv_par), path, TRUE);

        /* select item */
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);
        n_children = gtk_tree_model_iter_n_children(GTK_TREE_MODEL(xgc->ts_par), &iter_parent);
        gtk_tree_model_iter_children(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent);
        for (i = 0; i < n_children; i++) {
            gtk_tree_model_get(GTK_TREE_MODEL(xgc->ts_par), &iter, COLUMN_PARAMETER, &text, COLUMN_ID, &id, -1);
            if (id == id_to_select) {
                if (id == SET_PSOURCE)
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, xgc->data.set.ss.npnts - 1);
                else if (id == SET_POUT)    /* point output */
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, xgc->data.set.so.npnts);
                else if (id == SET_IIOUT)   /* image output */
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, xgc->data.set.so.npnts + xgc->data.set.so.nimgs);
                else if (id == SET_IOUT)    /* plane output */
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, xgc->data.set.so.npnts + xgc->data.set.so.nimgs + xgc->data.set.so.nplns);
                else if (id == SET_COUT)     /* volume output */
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, xgc->data.set.so.npnts + xgc->data.set.so.nimgs + xgc->data.set.so.nplns + xgc->data.set.so.ncubs);
                else if (id == SET_SOUT)    /* sum output */
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, xgc->data.set.so.npnts + xgc->data.set.so.nimgs + xgc->data.set.so.nplns + xgc->data.set.so.ncubs + xgc->data.set.so.nsums);
                else if (id == SET_FOUT)    /* force output */
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, xgc->data.set.so.npnts + xgc->data.set.so.nimgs + xgc->data.set.so.nplns + xgc->data.set.so.ncubs + xgc->data.set.so.nsums + xgc->data.set.so.nforces);
                else if (id == SET_GROW)    /* material grow */
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, NOL_MEDIUM_MAT_PROPS + xgc->data.set.sm.ngrowths - 1);
                else if (id == SET_ROUGHNESS) /* material roughen */
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, NOL_MEDIUM_MAT_PROPS + xgc->data.set.sm.ngrowths + xgc->data.set.sm.nroughens - 1);
                else if (id == SET_SPECTRAL) /* material spectral */
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, NOL_MEDIUM_MAT_PROPS + xgc->data.set.sm.ngrowths + xgc->data.set.sm.nroughens + xgc->data.set.sm.nspectrals - 1);
                else if (id == SET_EXPRESSION) /* material expression */
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, NOL_MEDIUM_MAT_PROPS + xgc->data.set.sm.ngrowths + xgc->data.set.sm.nroughens + xgc->data.set.sm.nspectrals + xgc->data.set.sm.nexprs - 1);
                else if (id_to_select >= SET_NFFFP && id_to_select < SET_PNFAREA)   /* VII. 2. Near field to far field point */
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, NOL_NFFF_PROPS + xgc->data.set.sf.nrs - 1);
                else if (id_to_select >= SET_NFAREA && id_to_select < SET_NFFFP)    /* VII. 3. Near field to far field area */
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, NOL_NFFF_PROPS + xgc->data.set.sf.nrs + xgc->data.set.sf.nareas - 1);
                else if (id_to_select >= SET_PNFFFP)                                /* VIII. 2. Periodic near field to far field point */
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, NOL_PNFFF_PROPS + xgc->data.set.spf.nrs - 1);
                else if (id_to_select >= SET_PNFAREA && id_to_select < SET_PNFFFP)  /* VIII. 3. Periodic near field to far field area */
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, NOL_PNFFF_PROPS + xgc->data.set.spf.nrs + xgc->data.set.spf.nareas - 1);

                gtk_tree_selection_select_iter(selection, &iter);

                path_to_select = gtk_tree_model_get_path(GTK_TREE_MODEL(xgc->ts_par), &iter);
                gtk_tree_view_scroll_to_cell(GTK_TREE_VIEW(xgc->tv_par), path_to_select, NULL, TRUE, 0, 0);
                break;
            }
            gtk_tree_model_iter_next(GTK_TREE_MODEL(xgc->ts_par), &iter);
        }
    }
} /* tree_expand_row_and_select_item */

void tree_mat_expand_row_and_select_item(XGControls *xgc, gint id_to_select)
{
    GtkTreePath *path = NULL;
    gchar *text;
    gint id, i, n_children;
    GtkTreeIter  iter, iter_parent;
    GtkTreeSelection *selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(xgc->tv_mat));

    if ((id_to_select == SET_MAT_SPHERE) || (id_to_select == SET_MAT_VOXEL) || (id_to_select == SET_MAT_CYLINDER) || (id_to_select == SET_MAT_CONE) || (id_to_select == SET_MAT_RCONE)
        || (id_to_select == SET_MAT_GWYDD) || (id_to_select == SET_MAT_MESH) || (id_to_select == SET_MAT_TETRAHEDRON))
        path = gtk_tree_path_new_from_indices(ROW_MAT_OBJECT_INDEX, -1);           /* expand "File" row */
    

    if (NULL != path) {
        gtk_tree_view_expand_row(GTK_TREE_VIEW(xgc->tv_mat), path, TRUE);

        /* select item */
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_mat), &iter_parent, path);
        n_children = gtk_tree_model_iter_n_children(GTK_TREE_MODEL(xgc->ts_mat), &iter_parent);
        gtk_tree_model_iter_children(GTK_TREE_MODEL(xgc->ts_mat), &iter, &iter_parent);
        for (i = 0; i < n_children; i++) {
            gtk_tree_model_get(GTK_TREE_MODEL(xgc->ts_mat), &iter, COLUMN_PARAMETER, &text, COLUMN_ID, &id, -1);
            if (id == id_to_select) {
                if (id == SET_MAT_SPHERE)
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_mat), &iter, &iter_parent, xgc->data.set_mat.spheres->len - 1);
                else if (id == SET_MAT_VOXEL)
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_mat), &iter, &iter_parent, xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len - 1);
                else if (id == SET_MAT_CYLINDER)
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_mat), &iter, &iter_parent, xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + xgc->data.set_mat.cylinders->len - 1);
                else if (id == SET_MAT_CONE)    
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_mat), &iter, &iter_parent, xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + xgc->data.set_mat.cylinders->len + xgc->data.set_mat.cones->len - 1);
                else if (id == SET_MAT_RCONE)   
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_mat), &iter, &iter_parent, xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + xgc->data.set_mat.cylinders->len + xgc->data.set_mat.cones->len + xgc->data.set_mat.rcones->len - 1);
                else if (id == SET_MAT_TETRAHEDRON)
                    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_mat), &iter, &iter_parent, xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + xgc->data.set_mat.cylinders->len + xgc->data.set_mat.cones->len + xgc->data.set_mat.rcones->len + xgc->data.set_mat.tetrahedrons->len - 1);

                gtk_tree_selection_select_iter(selection, &iter);
                break;
            }
            gtk_tree_model_iter_next(GTK_TREE_MODEL(xgc->ts_mat), &iter);
        }
    }
} /* tree_mat_expand_row_and_select_item */


void file_load(XGControls *xgc, gchar *fname, gchar* uri)
{
    gchar *contents;
    gchar *filename;
    gchar buff[256];

    filename = g_canonicalize_filename(fname, NULL);

    /*delete par and mat filenames*/
    g_free(xgc->data.parfilename);
    xgc->data.parfilename = NULL;
    g_free(xgc->data.matfilename);
    xgc->data.matfilename = NULL;

    /*delete settings*/
    delete_point_sources(&(xgc->data.set));

    /*clear settings*/
    clear_settings(&(xgc->data.set), TRUE);
    clear_settings_mat(&(xgc->data.set_mat), TRUE);

    /*init data*/
    alloc_set_par_mat(&xgc->data.set, &xgc->data.set_mat, TRUE);

    if (g_file_get_contents(filename, &contents, NULL, NULL)) {
        xgc->curdir = g_path_get_dirname(filename);

        if (g_strcmp0(xgc->curdir, ".") == 0) {
            /* file name has no directory components */
            g_free(xgc->curdir);
            xgc->curdir = g_get_current_dir();
        }
        else {
            g_chdir(xgc->curdir);
        }

        /*skip path!*/
        xgc->data.parfilename = g_path_get_basename(filename);

        set_main_window_title(xgc);
        //g_snprintf(title, sizeof(title), "%s - XSvit", xgc->data.parfilename);
        //gtk_window_set_title(GTK_WINDOW(xgc->toplevel), title);        

        //printf("current directory: %s\n", g_get_current_dir());
        //printf("xgc->curdir: %s\n", xgc->curdir);
        //printf("xgc->data.parfilename: %s\n", xgc->data.parfilename);

        g_snprintf(buff, sizeof(buff), "Data successfully loaded from %s", filename);
        gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);

        gtk_text_buffer_set_text(xgc->tb_par, contents, -1);
        pt_create_tree(xgc);
        mt_create_tree(xgc);

        //reset_visibility(xgc);

        scene_reset_cb(xgc->view_scene, xgc);

        // 2018-01-30 commented - gtk_text_buffer_set_text() implies par_file_changed() via par_buffer_changed_cb()
        //par_file_changed(xgc, TRUE);        
    } else {
        file_new_cb(xgc->toplevel, xgc);

        g_snprintf(buff, sizeof(buff), "Error loading parameter file %s.", filename);
        gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
    }

    g_free(filename);
} /* file_load */

void file_save(XGControls *xgc)
{
    GtkTextIter   start;
    GtkTextIter   end;
    gchar buff[256] = {0};

    write_parfile(xgc);

    gtk_text_buffer_get_start_iter(xgc->tb_par, &start);
    gtk_text_buffer_get_end_iter(xgc->tb_par, &end);

    if (g_file_set_contents(xgc->data.parfilename, gtk_text_buffer_get_text(xgc->tb_par, &start, &end, TRUE), -1, NULL)) {
        if (xgc->data.matfilename) {
            gtk_text_buffer_get_start_iter(xgc->tb_mat, &start);
            gtk_text_buffer_get_end_iter(xgc->tb_mat, &end);
            if (g_file_set_contents(xgc->data.matfilename, gtk_text_buffer_get_text(xgc->tb_mat, &start, &end, TRUE), -1, NULL)) {
                g_snprintf(buff, sizeof(buff), "Parameter file %s and material file %s successfully saved.", xgc->data.parfilename, xgc->data.matfilename/*g_build_filename(xgc->curdir, xgc->data.matfilename, NULL)*/);
                gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
            } else {
                g_snprintf(buff, sizeof(buff), "Parameter file %s successfully saved, error saving material file %s!", xgc->data.parfilename, g_build_filename(xgc->curdir, xgc->data.matfilename, NULL));
                gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
            }
        } else {
            g_snprintf(buff, sizeof(buff), "Parameter file %s successfully saved.", xgc->data.parfilename);
            gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
        }
    }
} /* file_save */


void file_new_cb(GtkWidget *widget, XGControls *xgc)
{
    gint i;
    GtkTreePath *path;

    /*delete par and mat filenames*/
    g_free(xgc->data.parfilename);
    xgc->data.parfilename = NULL;
    g_free(xgc->data.matfilename);
    xgc->data.matfilename = NULL;

    // delete current directory filename
    g_free(xgc->curdir);
    xgc->curdir = NULL;

    /*delete settings*/
    delete_point_sources(&(xgc->data.set));

    /*clear settings*/
    clear_settings(&(xgc->data.set), TRUE);
    clear_settings_mat(&(xgc->data.set_mat), TRUE);

    /*init data*/
    alloc_set_par_mat(&xgc->data.set, &xgc->data.set_mat, TRUE);

    for (i = 0; i < 100; i++)
        xgc->outputfield[i] = NULL;

    /*xgc->zoomfactor = 1;
    xgc->viewposx = 0;
    xgc->viewposy = 0;*/
    scene_reset_cb(xgc->view_scene, xgc);

    if (xgc->data.set_mat.gwyddata)
        sv_dcube_free(xgc->data.set_mat.gwyddata);
    xgc->data.set_mat.gwyddata = NULL;

    write_parfile(xgc);
    pt_create_tree(xgc);

    //gtk_tree_store_clear(xgc->ts_mat);
    if (write_matfile(xgc)) {
        path = gtk_tree_path_new_from_indices(0, -1);
        gtk_tree_view_expand_row(GTK_TREE_VIEW(xgc->tv_mat), path, TRUE);
    }
    mt_create_tree(xgc);

    reset_par_elemets_visibility(xgc);

    set_main_window_title(xgc);
}


void file_open_cb(GtkWidget *widget, XGControls *xgc)
{
    GtkWidget *dialog;
    GtkFileFilter *filter = gtk_file_filter_new();
    gchar* filename;
    gchar* dirname;

    dialog = gtk_file_chooser_dialog_new("Open parameter file (*.par)",
        GTK_WINDOW(xgc->toplevel),
        GTK_FILE_CHOOSER_ACTION_OPEN,
        GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
        GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
        NULL);
    gtk_file_filter_set_name(filter, "GSvit parameter files (*.par)");
    gtk_file_filter_add_pattern(filter, "*.par");
    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filter);
    if (NULL != xgc->curdir)
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), xgc->curdir);
    else {
        dirname = g_build_filename(g_get_home_dir(), ".xsvit", NULL);
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), dirname);
        g_free(dirname);
    }

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

        file_load(xgc, filename, NULL);            

        g_free(filename);
    }
    gtk_widget_destroy(dialog);

    reset_par_elemets_visibility(xgc);
}

void file_save_as_cb(GtkWidget *widget, XGControls *xgc)
{
    GtkWidget *dialog;
    gchar *filename;
    GtkFileFilter *filter = gtk_file_filter_new();
    GtkTextIter   start;
    GtkTextIter   end;
    gchar buff[256] = {0};
    gchar *matfilename;
    gchar *dirname;

    dialog = gtk_file_chooser_dialog_new("Save parameter and material files",
        GTK_WINDOW(xgc->toplevel),
        GTK_FILE_CHOOSER_ACTION_SAVE,
        GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
        GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
        NULL);
    gtk_file_filter_set_name(filter, "GSvit parameter files (*.par)");
    gtk_file_filter_add_pattern(filter, "*.par");

    gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(dialog), TRUE);
    if (!xgc->data.parfilename) {
        dirname = g_build_filename(g_get_home_dir(), ".xsvit", NULL);
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), dirname);
        gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(dialog), "Untitled.par");
        g_free(dirname);
    }
    else
        gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(dialog), xgc->data.parfilename);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

        xgc->curdir = g_path_get_dirname(filename);
        g_chdir(xgc->curdir);

        gtk_text_buffer_get_start_iter(xgc->tb_par, &start);
        gtk_text_buffer_get_end_iter(xgc->tb_par, &end);
        xgc->data.parfilename = g_strdup(filename);

        if (g_file_set_contents(xgc->data.parfilename, gtk_text_buffer_get_text(xgc->tb_par, &start, &end, TRUE), -1, NULL)) {
            if (xgc->data.matfilename) {
                gtk_text_buffer_get_start_iter(xgc->tb_mat, &start);
                gtk_text_buffer_get_end_iter(xgc->tb_mat, &end);

                if (g_path_is_absolute(xgc->data.matfilename))
                    matfilename = g_strdup(xgc->data.matfilename);
                else
                    matfilename = g_build_filename(xgc->curdir, xgc->data.matfilename, NULL);

                if (g_file_set_contents(matfilename, gtk_text_buffer_get_text(xgc->tb_mat, &start, &end, TRUE), -1, NULL)) {
                    g_snprintf(buff, sizeof(buff), "Data successfully saved to %s, material to %s", xgc->data.parfilename, matfilename);
                    gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
                }
                else {
                    if (g_path_is_absolute(xgc->data.matfilename))
                        g_snprintf(buff, sizeof(buff), "Data successfully saved to %s, error saving material to %s", xgc->data.parfilename, matfilename);
                    gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
                }
                g_free(matfilename);
            }
            else {
                g_snprintf(buff, sizeof(buff), "Data successfully saved to %s", xgc->data.parfilename);
                gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
            }

            set_main_window_title(xgc);
        }
        g_free(filename);
    }

    gtk_widget_destroy(dialog);
}

void file_save_cb(GtkWidget *widget, XGControls *xgc)
{
    if (!xgc->data.parfilename)
        file_save_as_cb(widget, xgc);
    else {
        //remeber_expanded_rows_and_selection(xgc);
        //write_parfile(xgc);
        //restore_expanded_rows_and_selection(xgc);

        file_save(xgc);
    }
}

void file_save_mat_as_cb(GtkWidget *widget, XGControls *xgc)
{

    GtkWidget *dialog;
    gchar *filename = NULL;
    gchar *dirname = NULL;
    GtkTextIter   start;
    GtkTextIter   end;
    gchar buff[256];
    gboolean saveok = FALSE;

    dialog = gtk_file_chooser_dialog_new("Save material file as..",
        GTK_WINDOW(xgc->toplevel),
        GTK_FILE_CHOOSER_ACTION_SAVE,
        GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
        GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
        NULL);

    gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(dialog), TRUE);
    if (!xgc->data.matfilename) {
        dirname = g_build_filename(g_get_home_dir(), ".xsvit", NULL);
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), dirname);
        gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(dialog), "Untitled");
        g_free(dirname);
    }
    else
        gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(dialog), xgc->data.matfilename);


    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

        gtk_text_buffer_get_start_iter(xgc->tb_mat, &start);
        gtk_text_buffer_get_end_iter(xgc->tb_mat, &end);
        if (g_file_set_contents(filename, gtk_text_buffer_get_text(xgc->tb_mat, &start, &end, TRUE), -1, NULL))
            g_snprintf(buff, sizeof(buff), "Material file %s successfully saved.", filename);
        gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
        saveok = TRUE;
    }

    gtk_widget_destroy(dialog);

    if (saveok) {
        dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
            GTK_DIALOG_DESTROY_WITH_PARENT,
            GTK_MESSAGE_QUESTION,
            GTK_BUTTONS_YES_NO,
            "Update also material file location in parameter file?");
        if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_YES) {
            if (g_strcmp0(g_get_current_dir(), g_path_get_dirname(filename)) == 0) {
                xgc->data.set.sm.in_vector_filename = g_path_get_basename(filename);
            }
            else
                xgc->data.set.sm.in_vector_filename = filename;
            write_parfile(xgc);
        }
        gtk_widget_destroy(dialog);
    }
    g_free(filename);
}

void recent_chooser_item_activated_cb(GtkWidget *widget, XGControls *xgc)
{
    gchar *uri = gtk_recent_chooser_get_current_uri(GTK_RECENT_CHOOSER(widget));
    gchar *filename = g_filename_from_uri(uri, NULL, NULL);

    file_load(xgc, filename, uri);

    reset_par_elemets_visibility(xgc);

    g_free(uri);
    g_free(filename);
}

void preferences_cb(GtkWidget *widget, XGControls *xgc)
{
    xsv_get_preferences(xgc);
}

void quit_cb(GtkWidget *widget, XGControls *xgc)
{
    gsvit_stop_cb(widget, xgc);

    /* cleanup */
    destroy_layer_treeview();

    clear_settings(&(xgc->data.set), TRUE);
    clear_settings_mat(&(xgc->data.set_mat), TRUE);

    if (NULL != xgc->tmpfilename)
        remove(xgc->tmpfilename);

    if (NULL != xgc->tmpmatfilename)
        remove(xgc->tmpmatfilename);


    gtk_main_quit();
}

void
par_add_src_point_cb(GtkWidget *widget, XGControls *xgc)
{    
    init_settings(&(xgc->data.set), SV_SRCTYPE_POINT);
    xgc->data.is_psrc[xgc->data.set.ss.npnts-1] = TRUE;

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_PSOURCE);
    }
}

void
par_add_src_sf_cb(GtkWidget *widget, XGControls *xgc)
{
    /* only one SF object allowed */
    if (xgc->data.set.ss.sf.valid == TRUE) {
        tree_expand_row_and_select_item(xgc, SET_SF);
        return;
    }

    xgc->data.is_sf = TRUE;
    init_settings(&(xgc->data.set), SV_SRCTYPE_SF);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_SF);
    }
}

void
par_add_src_tsf_cb(GtkWidget *widget, XGControls *xgc)
{    
    /* only one TSF object allowed */
    if (xgc->data.set.ss.tsf.valid == TRUE) {  
        tree_expand_row_and_select_item(xgc, SET_TSF);
        return;
    }

    init_settings(&(xgc->data.set), SV_SRCTYPE_TSF);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_TSF);
    }
}

void
par_add_src_tsff_cb(GtkWidget *widget, XGControls *xgc)
{
    /* only one TSFF object allowed */
    if (xgc->data.set.ss.tsff.valid == TRUE) {  
        tree_expand_row_and_select_item(xgc, SET_TSFF);
        return;
    }

    init_settings(&(xgc->data.set), SV_SRCTYPE_TSFF);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_TSFF);        
    }
}

void
par_add_src_ltsf_cb(GtkWidget *widget, XGControls *xgc)
{
    /* only one LTSF object allowed */
    if (xgc->data.set.ss.ltsf.valid == TRUE) {  
        tree_expand_row_and_select_item(xgc, SET_LTSF);
        return;
    }

    init_settings(&(xgc->data.set), SV_SRCTYPE_LTSF);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_LTSF);
    }
}

void
par_add_src_ltsff_cb(GtkWidget *widget, XGControls *xgc)
{
    /* only one LTSFF object allowed */
    if (xgc->data.set.ss.ltsff.valid == TRUE) {  
        tree_expand_row_and_select_item(xgc, SET_LTSFF);
        return;
    }

    init_settings(&(xgc->data.set), SV_SRCTYPE_LTSFF);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_LTSFF);
    }
}

void
par_add_out_point_cb(GtkWidget *widget, XGControls *xgc)
{
    init_settings(&(xgc->data.set), SV_OUTTYPE_POINT);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_POUT);
    }
}

void
par_add_out_image_cb(GtkWidget *widget, XGControls *xgc)
{
    init_settings(&(xgc->data.set), SV_OUTTYPE_IMAGE);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_IIOUT);
    }
}

void
par_add_out_plane_cb(GtkWidget *widget, XGControls *xgc)
{
    init_settings(&(xgc->data.set), SV_OUTTYPE_PLANE);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_IOUT);
    }
}

void
par_add_out_sum_cb(GtkWidget *widget, XGControls *xgc)
{
    init_settings(&(xgc->data.set), SV_OUTTYPE_SUM);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_SOUT);
    }
}

void
par_add_out_force_cb(GtkWidget *widget, XGControls *xgc)
{
    init_settings(&(xgc->data.set), SV_OUTTYPE_FORCE);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_FOUT);
    }
}

void
par_add_out_volume_cb(GtkWidget *widget, XGControls *xgc)
{
    init_settings(&(xgc->data.set), SV_OUTTYPE_VOLUME);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_COUT);
    }
}

void
par_add_grow_cb(GtkWidget *widget, XGControls *xgc)
{
    init_settings(&(xgc->data.set), SV_TYPE_MATERIAL_GROW);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_GROW);
    }
}

void
par_add_roughness_cb(GtkWidget *widget, XGControls *xgc)
{
    init_settings(&(xgc->data.set), SV_TYPE_MATERIAL_ROUGHNESS);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_ROUGHNESS);    
    }
}

void
par_add_spectral_cb(GtkWidget *widget, XGControls *xgc)
{
    init_settings(&(xgc->data.set), SV_TYPE_MATERIAL_SPECTRAL);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_SPECTRAL);
    }
}

void
par_add_expression_cb(GtkWidget *widget, XGControls *xgc)
{
    init_settings(&(xgc->data.set), SV_TYPE_MATERIAL_EXPRESSION);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_EXPRESSION);
    }
}

void
par_add_nfffp_cb(GtkWidget *widget, XGControls *xgc)
{
    init_settings(&(xgc->data.set), SV_TYPE_NFFF_POINT);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_NFFFP);
    }
}

void
par_add_pnfffp_cb(GtkWidget *widget, XGControls *xgc)
{
    init_settings(&(xgc->data.set), SV_TYPE_PNFFF_POINT);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_PNFFFP);
    }
}

void
par_add_nfffa_cb(GtkWidget *widget, XGControls *xgc)
{
    init_settings(&(xgc->data.set), SV_TYPE_NFFF_AREA);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_NFAREA);    
    }
}

void
par_add_pnfffa_cb(GtkWidget *widget, XGControls *xgc)
{
    init_settings(&(xgc->data.set), SV_TYPE_PNFFF_AREA);

    if (write_parfile(xgc)) {
        pt_create_tree(xgc);
        tree_expand_row_and_select_item(xgc, SET_PNFAREA);
    }
}

void
mat_add_sphere_cb(GtkWidget *widget, XGControls *xgc)
{
    init_add_settings_mat(&(xgc->data.set_mat), SV_TYPE_MAT_SPHERE);

    if (write_matfile(xgc)) {
        mt_create_tree(xgc);
        tree_mat_expand_row_and_select_item(xgc, SET_MAT_SPHERE);
    }

    /*GtkTreePath *path;
    gboolean writeme = FALSE;

    writeme = xsv_get_object(xgc, xgc->data.set_mat.spheres->len, 0);

    if (writeme) {
        if (write_matfile(xgc)) {
            path = gtk_tree_path_new_from_indices(0, -1);
            gtk_tree_view_expand_row(GTK_TREE_VIEW(xgc->tv_mat), path, TRUE);
        }
    }*/
}

void
mat_add_box_cb(GtkWidget *widget, XGControls *xgc)
{
    init_add_settings_mat(&(xgc->data.set_mat), SV_TYPE_MAT_VOXEL);

    if (write_matfile(xgc)) {
        mt_create_tree(xgc);
        tree_mat_expand_row_and_select_item(xgc, SET_MAT_VOXEL);
    }

    /*GtkTreePath *path;
    gboolean writeme = FALSE;

    writeme = xsv_get_object(xgc, xgc->data.set_mat.voxels->len, 1);

    if (writeme) {
        if (write_matfile(xgc)) {
            path = gtk_tree_path_new_from_indices(0, -1);
            gtk_tree_view_expand_row(GTK_TREE_VIEW(xgc->tv_mat), path, TRUE);
        }
    }*/
}

void
mat_add_cylinder_cb(GtkWidget *widget, XGControls *xgc)
{
    init_add_settings_mat(&(xgc->data.set_mat), SV_TYPE_MAT_CYLINDER);

    if (write_matfile(xgc)) {
        mt_create_tree(xgc);
        tree_mat_expand_row_and_select_item(xgc, SET_MAT_CYLINDER);
    }
}

void
mat_add_cone_cb(GtkWidget *widget, XGControls *xgc)
{
    init_add_settings_mat(&(xgc->data.set_mat), SV_TYPE_MAT_CONE);

    if (write_matfile(xgc)) {
        mt_create_tree(xgc);
        tree_mat_expand_row_and_select_item(xgc, SET_MAT_CONE);
    }
}
void
mat_add_rcone_cb(GtkWidget *widget, XGControls *xgc)
{
    init_add_settings_mat(&(xgc->data.set_mat), SV_TYPE_MAT_RCONE);

    if (write_matfile(xgc)) {
        mt_create_tree(xgc);
        tree_mat_expand_row_and_select_item(xgc, SET_MAT_RCONE);
    }
}

/*
void
add_mat_gwydd_cb(GtkWidget *widget, XGControls *xgc)
{
    GtkTreePath *path;
    gboolean writeme = FALSE;

    writeme = xsv_get_object(xgc, xgc->data.gwydds->len, 5);

    if (writeme) {
        if (write_matfile(xgc)) {
            path = gtk_tree_path_new_from_indices(0, -1);
            gtk_tree_view_expand_row(GTK_TREE_VIEW(xgc->tv_mat), path, TRUE);
        }
    }
}

void
add_mat_mesh_cb(GtkWidget *widget, XGControls *xgc)
{
    GtkTreePath *path;
    gboolean writeme = FALSE;

    writeme = xsv_get_object(xgc, xgc->data.ntetgens, 6);

    if (writeme) {
        if (write_matfile(xgc)) {
            path = gtk_tree_path_new_from_indices(0, -1);
            gtk_tree_view_expand_row(GTK_TREE_VIEW(xgc->tv_mat), path, TRUE);
        }
    }
}*/

void
mat_add_tetrahedron_cb(GtkWidget *widget, XGControls *xgc)
{
    init_add_settings_mat(&(xgc->data.set_mat), SV_TYPE_MAT_TETRAHEDRON);

    if (write_matfile(xgc)) {
        mt_create_tree(xgc);
        tree_mat_expand_row_and_select_item(xgc, SET_MAT_TETRAHEDRON);
    }
}

static void
par_row_remove_process(XGControls *xgc)
{
    gint i, index;
    gchar *text;
    GtkTreePath *path;
    gboolean writeme = FALSE;    
    gint id = xgc->par_selid;    

    if (xgc->par_selid == -1)
        return;

    path = gtk_tree_model_get_path(GTK_TREE_MODEL(xgc->ts_par), &(xgc->par_seliter));

    ///////////////////////////////////////////////////////////`///////////////
    GtkTreeIter iter;
    GtkTreePath *path_to_select;
    gchar* path_string_to_select = NULL;

    //model = gtk_tree_view_get_model(xgc->tv_par);
    //model = GTK_TREE_MODEL(xgc->ts_par);
    GtkTreeSelection *selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(xgc->tv_par));

    if (gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter, path)) {
        gtk_tree_model_get(GTK_TREE_MODEL(xgc->ts_par), &iter, COLUMN_PARAMETER, &text, COLUMN_ID, &id, -1);
        
        path_to_select = gtk_tree_model_get_path(GTK_TREE_MODEL(xgc->ts_par), &iter);        
        path_string_to_select = gtk_tree_path_to_string(path_to_select);
    }

    /*if (TRUE == gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(xgc->ts_par), &iter, path)) {
        gtk_tree_model_get(GTK_TREE_MODEL(xgc->ts_par), &iter, COLUMN_PARAMETER, &text, COLUMN_ID, &id, -1);

        gtk_tree_model_iter_next(GTK_TREE_MODEL(xgc->ts_par), &iter);
        gtk_tree_model_get(GTK_TREE_MODEL(xgc->ts_par), &iter, COLUMN_PARAMETER, &text, COLUMN_ID, &id, -1);
    }*/
    //////////////////////////////////////////////////////////////////////////
    

    /*use id to call appropriate callback*/
    if (id == SET_SF) {
        g_free(xgc->data.set.ss.sf.source_filename);
        xgc->data.set.ss.sf.valid = FALSE;
        writeme = TRUE;
    } else if (id == SET_TSF) {
        xgc->data.set.ss.tsf.valid = FALSE;
        writeme = TRUE;
    } else if (id == SET_TSFF) {
        xgc->data.set.ss.tsff.valid = FALSE;
        writeme = TRUE;
    } else if (id == SET_LTSF) {
        xgc->data.set.ss.ltsf.valid = FALSE;
        for (i = 0; i < xgc->data.set.ss.ltsf.layered_count; i++) {
            g_free(xgc->data.set.ss.ltsf.layered_material[i]);
            xgc->data.set.ss.ltsf.layered_material[i] = NULL;
        }
        writeme = TRUE;
    } else if (id == SET_LTSFF) {
        xgc->data.set.ss.ltsff.valid = FALSE;
        for (i = 0; i < xgc->data.set.ss.ltsff.layered_count; i++) {
            g_free(xgc->data.set.ss.ltsff.layered_material[i]);
            xgc->data.set.ss.ltsff.layered_material[i] = NULL;
        }
        writeme = TRUE;
    } else if (id >= SET_PSOURCE && id < SET_POUT) {    /* point source */
        if (xgc->data.set.ss.npnts == 1)
            xgc->data.set.ss.npnts = 0;
        else {
            for (i = id - SET_PSOURCE; i < (xgc->data.set.ss.npnts - 1); i++) {
                xgc->data.set.ss.pnts[i].point_origin_position_i = xgc->data.set.ss.pnts[i + 1].point_origin_position_i;
                xgc->data.set.ss.pnts[i].point_origin_position_j = xgc->data.set.ss.pnts[i + 1].point_origin_position_j;
                xgc->data.set.ss.pnts[i].point_origin_position_k = xgc->data.set.ss.pnts[i + 1].point_origin_position_k;
                xgc->data.set.ss.pnts[i].point_origin_theta = xgc->data.set.ss.pnts[i + 1].point_origin_theta;
                xgc->data.set.ss.pnts[i].point_origin_phi = xgc->data.set.ss.pnts[i + 1].point_origin_phi;
                xgc->data.set.ss.pnts[i].source_mode = xgc->data.set.ss.pnts[i + 1].source_mode;
                g_free(xgc->data.set.ss.pnts[i].source_filename);
                xgc->data.set.ss.pnts[i].source_filename = xgc->data.set.ss.pnts[i + 1].source_filename;
                xgc->data.set.ss.pnts[i].source_wl = xgc->data.set.ss.pnts[i + 1].source_wl;
                xgc->data.set.ss.pnts[i].source_pulsewidth = xgc->data.set.ss.pnts[i + 1].source_pulsewidth;
                xgc->data.set.ss.pnts[i].source_amplitude = xgc->data.set.ss.pnts[i + 1].source_amplitude;

            }
            xgc->data.set.ss.npnts -= 1;
        }
        writeme = TRUE;
    } else if (id >= SET_POUT && id < SET_IOUT) {   /* point output */
        if (xgc->data.set.so.npnts == 1)
            xgc->data.set.so.npnts = 0;
        else {
            for (i = id - SET_POUT; i < (xgc->data.set.so.npnts - 1); i++) {
                xgc->data.set.so.pnts[i].step = xgc->data.set.so.pnts[i + 1].step;
                xgc->data.set.so.pnts[i].i = xgc->data.set.so.pnts[i + 1].i;
                xgc->data.set.so.pnts[i].j = xgc->data.set.so.pnts[i + 1].j;
                xgc->data.set.so.pnts[i].k = xgc->data.set.so.pnts[i + 1].k;
                xgc->data.set.so.pnts[i].component = xgc->data.set.so.pnts[i + 1].component;
                g_snprintf(xgc->data.set.so.pnts[i].filebase, 100, "%s", xgc->data.set.so.pnts[i + 1].filebase);
            }
            xgc->data.set.so.npnts -= 1;
        }
        writeme = TRUE;
    } else if (id >= SET_SOUT && id < SET_FOUT) {   /* sum output */
        if (xgc->data.set.so.nsums == 1)
            xgc->data.set.so.nsums = 0;
        else {
            for (i = id - SET_SOUT; i < (xgc->data.set.so.nsums - 1); i++) {
                xgc->data.set.so.sums[i].step = xgc->data.set.so.sums[i + 1].step;
                xgc->data.set.so.sums[i].component = xgc->data.set.so.sums[i + 1].component;
                xgc->data.set.so.sums[i].box_i0 = xgc->data.set.so.sums[i + 1].box_i0;
                xgc->data.set.so.sums[i].box_j0 = xgc->data.set.so.sums[i + 1].box_j0;
                xgc->data.set.so.sums[i].box_k0 = xgc->data.set.so.sums[i + 1].box_k0;
                xgc->data.set.so.sums[i].box_in = xgc->data.set.so.sums[i + 1].box_in;
                xgc->data.set.so.sums[i].box_jn = xgc->data.set.so.sums[i + 1].box_jn;
                xgc->data.set.so.sums[i].box_kn = xgc->data.set.so.sums[i + 1].box_kn;
                xgc->data.set.so.sums[i].layered_epsilon = xgc->data.set.so.sums[i + 1].layered_epsilon;
                xgc->data.set.so.sums[i].layered_mu = xgc->data.set.so.sums[i + 1].layered_mu;
                xgc->data.set.so.sums[i].layered_sigma = xgc->data.set.so.sums[i + 1].layered_sigma;
                xgc->data.set.so.sums[i].layered_sigast = xgc->data.set.so.sums[i + 1].layered_sigast;
                xgc->data.set.so.sums[i].stringbased = xgc->data.set.so.sums[i + 1].stringbased;
                g_snprintf(xgc->data.set.so.sums[i].string, 100, "%s", xgc->data.set.so.sums[i + 1].string);
                g_snprintf(xgc->data.set.so.sums[i].filename, 100, "%s", xgc->data.set.so.sums[i + 1].filename);
            }
            xgc->data.set.so.nsums -= 1;
        }
        writeme = TRUE;
    } else if (id >= SET_FOUT && id < SET_GROW) {    /* force output*/
        if (xgc->data.set.so.nforces == 1)
            xgc->data.set.so.nforces = 0;
        else {
            for (i = id - SET_FOUT; i < (xgc->data.set.so.nforces - 1); i++) {
                xgc->data.set.so.forces[i].step = xgc->data.set.so.forces[i + 1].step;
                xgc->data.set.so.forces[i].box_i0 = xgc->data.set.so.forces[i + 1].box_i0;
                xgc->data.set.so.forces[i].box_j0 = xgc->data.set.so.forces[i + 1].box_j0;
                xgc->data.set.so.forces[i].box_k0 = xgc->data.set.so.forces[i + 1].box_k0;
                xgc->data.set.so.forces[i].box_in = xgc->data.set.so.forces[i + 1].box_in;
                xgc->data.set.so.forces[i].box_jn = xgc->data.set.so.forces[i + 1].box_jn;
                xgc->data.set.so.forces[i].box_kn = xgc->data.set.so.forces[i + 1].box_kn;
                g_snprintf(xgc->data.set.so.forces[i].filename, 100, "%s", xgc->data.set.so.forces[i + 1].filename);
            }
            xgc->data.set.so.nforces -= 1;
        }
        writeme = TRUE;
    } else if (id >= SET_IOUT && id < SET_IIOUT) {  /* plane output */
        if (xgc->data.set.so.nplns == 1)
            xgc->data.set.so.nplns = 0;
        else {
            for (i = id - SET_IOUT; i < (xgc->data.set.so.nplns - 1); i++) {
                xgc->data.set.so.plns[i].step = xgc->data.set.so.plns[i + 1].step;
                xgc->data.set.so.plns[i].i = xgc->data.set.so.plns[i + 1].i;
                xgc->data.set.so.plns[i].j = xgc->data.set.so.plns[i + 1].j;
                xgc->data.set.so.plns[i].k = xgc->data.set.so.plns[i + 1].k;
                xgc->data.set.so.plns[i].component = xgc->data.set.so.plns[i + 1].component;
                g_snprintf(xgc->data.set.so.plns[i].filebase, 100, "%s", xgc->data.set.so.plns[i + 1].filebase);
            }
            xgc->data.set.so.nplns -= 1;
        }
        writeme = TRUE;
    } else if (id >= SET_IIOUT && id < SET_COUT) {  /* image output */
        if (xgc->data.set.so.nimgs == 1)
            xgc->data.set.so.nimgs = 0;
        else {
            for (i = id - SET_IIOUT; i < (xgc->data.set.so.nimgs - 1); i++) {
                xgc->data.set.so.imgs[i].step = xgc->data.set.so.imgs[i + 1].step;
                xgc->data.set.so.imgs[i].i = xgc->data.set.so.imgs[i + 1].i;
                xgc->data.set.so.imgs[i].j = xgc->data.set.so.imgs[i + 1].j;
                xgc->data.set.so.imgs[i].k = xgc->data.set.so.imgs[i + 1].k;
                xgc->data.set.so.imgs[i].component = xgc->data.set.so.imgs[i + 1].component;
                g_snprintf(xgc->data.set.so.imgs[i].filebase, 100, "%s", xgc->data.set.so.imgs[i + 1].filebase);
            }
            xgc->data.set.so.nimgs -= 1;
        }
        writeme = TRUE;
    } else if (id >= SET_COUT && id < SET_SOUT) {   /* volume output */
        if (xgc->data.set.so.ncubs == 1)
            xgc->data.set.so.ncubs = 0;
        else {
            for (i = id - SET_COUT; i < (xgc->data.set.so.ncubs - 1); i++) {
                xgc->data.set.so.cubs[i].step = xgc->data.set.so.cubs[i + 1].step;
                xgc->data.set.so.cubs[i].start = xgc->data.set.so.cubs[i + 1].start;
                xgc->data.set.so.cubs[i].stop = xgc->data.set.so.cubs[i + 1].stop;
                xgc->data.set.so.cubs[i].format = xgc->data.set.so.cubs[i + 1].format;
                xgc->data.set.so.cubs[i].component = xgc->data.set.so.cubs[i + 1].component;
                g_snprintf(xgc->data.set.so.cubs[i].filebase, 100, "%s", xgc->data.set.so.cubs[i + 1].filebase);
            }
            xgc->data.set.so.ncubs -= 1;
        }
        writeme = TRUE;
    } else if (id == SET_NFFF) {                        /* VII. 1. Near field to far field transform box */
        //xgc->data.set.sf.nrs = xgc->data.set.sf.nareas = 0;
        //writeme = TRUE;
    } else if (id == SET_PNFFF) {                       /* VIII. 1. Near field to far field transform box */   
        //xgc->data.set.spf.nrs = xgc->data.set.spf.nareas = 0;
        //writeme = TRUE;
    } else if (id >= SET_NFFFP && id < SET_PNFAREA) {   /* VII. 2. Near field to far field point */
        if (xgc->data.set.sf.nrs == 1)
            xgc->data.set.sf.nrs = 0;
        else {
            for (i = id - SET_NFFFP; i < (xgc->data.set.sf.nrs - 1); i++) {
                xgc->data.set.sf.ri[i] = xgc->data.set.sf.ri[i + 1];
                xgc->data.set.sf.rj[i] = xgc->data.set.sf.rj[i + 1];
                xgc->data.set.sf.rk[i] = xgc->data.set.sf.rk[i + 1];
                xgc->data.set.sf.individual[i] = xgc->data.set.sf.individual[i + 1];
                xgc->data.set.sf.source_filename[i] = g_strdup_printf("%s", xgc->data.set.sf.source_filename[i + 1]);
            }
            xgc->data.set.sf.nrs -= 1;
        }
        writeme = TRUE;
    } else if (id >= SET_PNFFFP) {                      /* VIII. 2. Periodic near field to far field point */
        if (xgc->data.set.spf.nrs == 1)
            xgc->data.set.spf.nrs = 0;
        else {
            for (i = id - SET_PNFFFP; i < (xgc->data.set.spf.nrs - 1); i++) {
                xgc->data.set.spf.ri[i] = xgc->data.set.spf.ri[i + 1];
                xgc->data.set.spf.rj[i] = xgc->data.set.spf.rj[i + 1];
                xgc->data.set.spf.rk[i] = xgc->data.set.spf.rk[i + 1];
                xgc->data.set.spf.individual[i] = xgc->data.set.spf.individual[i + 1];
                xgc->data.set.spf.source_filename[i] = g_strdup_printf("%s", xgc->data.set.spf.source_filename[i + 1]);
            }
            xgc->data.set.spf.nrs -= 1;
        }
        writeme = TRUE;
    } else if (id >= SET_NFAREA && id < SET_NFFFP) {    /* VII. 3. Near field to far field area */
        if (xgc->data.set.sf.nareas == 1)
            xgc->data.set.sf.nareas = 0;
        else {
            for (i = id - SET_NFAREA; i < (xgc->data.set.sf.nareas - 1); i++) {
                xgc->data.set.sf.area_thetares[i] = xgc->data.set.sf.area_thetares[i + 1];
                xgc->data.set.sf.area_phires[i] = xgc->data.set.sf.area_phires[i + 1];
                xgc->data.set.sf.area_radius[i] = xgc->data.set.sf.area_radius[i + 1];
                xgc->data.set.sf.area_thetafrom[i] = xgc->data.set.sf.area_thetafrom[i + 1];
                xgc->data.set.sf.area_phifrom[i] = xgc->data.set.sf.area_phifrom[i + 1];
                xgc->data.set.sf.area_thetato[i] = xgc->data.set.sf.area_thetato[i + 1];
                xgc->data.set.sf.area_phito[i] = xgc->data.set.sf.area_phito[i + 1];
                xgc->data.set.sf.area_savefile[i] = xgc->data.set.sf.area_savefile[i + 1];
            }
            xgc->data.set.sf.nareas -= 1;
        }
        writeme = TRUE;
    } else if (id >= SET_PNFAREA && id < SET_PNFFFP) {  /* VIII. 3. Periodic near field to far field area */
        if (xgc->data.set.spf.nareas == 1)
            xgc->data.set.spf.nareas = 0;
        else {
            for (i = id - SET_PNFAREA; i < (xgc->data.set.spf.nareas - 1); i++) {
                xgc->data.set.spf.area_thetares[i] = xgc->data.set.spf.area_thetares[i + 1];
                xgc->data.set.spf.area_phires[i] = xgc->data.set.spf.area_phires[i + 1];
                xgc->data.set.spf.area_radius[i] = xgc->data.set.spf.area_radius[i + 1];
                xgc->data.set.spf.area_thetafrom[i] = xgc->data.set.spf.area_thetafrom[i + 1];
                xgc->data.set.spf.area_phifrom[i] = xgc->data.set.spf.area_phifrom[i + 1];
                xgc->data.set.spf.area_thetato[i] = xgc->data.set.spf.area_thetato[i + 1];
                xgc->data.set.spf.area_phito[i] = xgc->data.set.spf.area_phito[i + 1];
                xgc->data.set.spf.area_savefile[i] = xgc->data.set.spf.area_savefile[i + 1];
            }
            xgc->data.set.spf.nareas -= 1;
        }
        writeme = TRUE;
    } else if (id >= SET_ROUGHNESS && id < SET_GROW) {    /* material roughen */
        if (xgc->data.set.sm.nroughens == 1)
            xgc->data.set.sm.nroughens = 0;
        else {
            for (i = id - SET_ROUGHNESS; i < (xgc->data.set.sm.nroughens - 1); i++) {
                xgc->data.set.sm.rough_radius_peak[i] = xgc->data.set.sm.rough_radius_peak[i + 1];
                xgc->data.set.sm.rough_radius_span[i] = xgc->data.set.sm.rough_radius_span[i + 1];
                xgc->data.set.sm.rough_iterations[i] = xgc->data.set.sm.rough_iterations[i + 1];
                xgc->data.set.sm.rough_probability[i] = xgc->data.set.sm.rough_probability[i + 1];
                xgc->data.set.sm.rough_matindex[i] = xgc->data.set.sm.rough_matindex[i + 1];
                xgc->data.set.sm.rough_voidindex[i] = xgc->data.set.sm.rough_voidindex[i + 1];
                xgc->data.set.sm.rough_seed[i] = xgc->data.set.sm.rough_seed[i + 1];
            }
            xgc->data.set.sm.nroughens -= 1;
        }
        writeme = TRUE;
    } else if (id >= SET_GROW && id < SET_NFAREA) {     /* material grow */
        if (xgc->data.set.sm.ngrowths == 1)
            xgc->data.set.sm.ngrowths = 0;
        else {
            for (i = id - SET_GROW; i < (xgc->data.set.sm.ngrowths - 1); i++) {
                xgc->data.set.sm.grow_i0[i] = xgc->data.set.sm.grow_i0[i + 1];
                xgc->data.set.sm.grow_in[i] = xgc->data.set.sm.grow_in[i + 1];
                xgc->data.set.sm.grow_j0[i] = xgc->data.set.sm.grow_j0[i + 1];
                xgc->data.set.sm.grow_jn[i] = xgc->data.set.sm.grow_jn[i + 1];
                xgc->data.set.sm.grow_k0[i] = xgc->data.set.sm.grow_k0[i + 1];
                xgc->data.set.sm.grow_kn[i] = xgc->data.set.sm.grow_kn[i + 1];
                xgc->data.set.sm.grow_skipi0[i] = xgc->data.set.sm.grow_skipi0[i + 1];
                xgc->data.set.sm.grow_skipj0[i] = xgc->data.set.sm.grow_skipj0[i + 1];
                xgc->data.set.sm.grow_skipk0[i] = xgc->data.set.sm.grow_skipk0[i + 1];
                xgc->data.set.sm.grow_skipin[i] = xgc->data.set.sm.grow_skipin[i + 1];
                xgc->data.set.sm.grow_skipjn[i] = xgc->data.set.sm.grow_skipjn[i + 1];
                xgc->data.set.sm.grow_skipkn[i] = xgc->data.set.sm.grow_skipkn[i + 1];
                xgc->data.set.sm.grow_addindex[i] = xgc->data.set.sm.grow_addindex[i + 1];
                xgc->data.set.sm.grow_attachindex[i] = xgc->data.set.sm.grow_attachindex[i + 1];
                xgc->data.set.sm.grow_subsampling[i] = xgc->data.set.sm.grow_subsampling[i + 1];
                xgc->data.set.sm.grow_seed[i] = xgc->data.set.sm.grow_seed[i + 1];
                xgc->data.set.sm.grow_nsteps[i] = xgc->data.set.sm.grow_nsteps[i + 1];
                xgc->data.set.sm.grow_mobility[i] = xgc->data.set.sm.grow_mobility[i + 1];
                xgc->data.set.sm.grow_probability[i] = xgc->data.set.sm.grow_probability[i + 1];
            }
            xgc->data.set.sm.ngrowths -= 1;
        }
        writeme = TRUE;
    }

    /*rewrite the temporary parfile and then reload it and interpret completely in order to assure consistency between its text and gui interpretation*/
    if (writeme) {
        if (write_parfile(xgc)) {
            pt_create_tree(xgc);
            index = gtk_tree_path_get_indices(path)[0];
            path = gtk_tree_path_new_from_indices(index, -1);
            gtk_tree_view_expand_row(GTK_TREE_VIEW(xgc->tv_par), path, TRUE);
            

            //selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(xgc->tv_par));
            //gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(xgc->ts_par), &iter, path_string_to_select);
            //gtk_tree_selection_select_iter(selection, &iter);

            selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(xgc->tv_par));
            if (gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(xgc->ts_par), &iter, path_string_to_select))
                gtk_tree_selection_select_iter(selection, &iter);
            else {
                gint n;
                GtkTreeIter iter_parent;  

                if (gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path)) {
                    n = gtk_tree_model_iter_n_children(GTK_TREE_MODEL(xgc->ts_par), &iter);
                    if (gtk_tree_model_iter_nth_child (GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, n-1)) {
                        gtk_tree_selection_select_iter(selection, &iter);
                    } else
                        gtk_tree_selection_select_iter(selection, &iter_parent);
                }
            }

            //pc_row_select_process(xgc, id_to_select, path_to_select);
        }
    }

    if (path_string_to_select != NULL)
        g_free(path_string_to_select);
} /* par_row_remove_process */

void par_row_remove_cb(GtkButton *button, XGControls *xgc)
{
    //gint id = xgc->par_selid;
    //GtkTreePath *path;
    //
    //if (xgc->par_selid == -1)
    //    return;
    //
    //path = gtk_tree_model_get_path(GTK_TREE_MODEL(xgc->ts_par), &(xgc->par_seliter));
    par_row_remove_process(xgc);
}

gboolean par_row_selected_cb(GtkTreeSelection *selection, XGControls *xgc)
{
    //    if (block_selection_change)
    //        return FALSE;

    gint id;
    //gboolean delete;
    GtkTreeIter iter;
    gchar *text;
    GtkTreeModel *model = GTK_TREE_MODEL(xgc->ts_par);
    GtkTreePath *path;
    //gboolean delete = FALSE;
    xgc->par_selid = -1;
    xgc->par_sel_removable = FALSE;

    set_all_controls_invisible(xgc);

    if (gtk_tree_selection_get_selected(selection, &model, &iter)) {
        gtk_tree_model_get(GTK_TREE_MODEL(xgc->ts_par), &iter, COLUMN_PARAMETER, &text, COLUMN_ID, &id, -1);

        //
//        if (xgc->par_selid == id)
//            return FALSE;
//        xgc->par_selid = -1;
//        ss_set_all_controls_invisible(xgc);
        //

        xgc->par_selid = id;
        xgc->par_seliter = iter;

        //
        if (xgc->par_selid < 0)
            return FALSE;

        path = gtk_tree_model_get_path(GTK_TREE_MODEL(xgc->ts_par), &(xgc->par_seliter));
        par_controls_row_select_process(xgc, id, path);
        //

        if (id == SET_POOL) {
            xgc->par_sel_removable = FALSE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id == SET_BASIC) {
            xgc->par_sel_removable = FALSE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id == SET_SF) {
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id == SET_TSF) {
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id == SET_TSFF) {
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id == SET_LTSF) {
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id == SET_LTSFF) {
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id >= SET_PSOURCE && id < SET_POUT) {    /* point source */
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id == SET_BND) {                         /* boundary conditions */
            xgc->par_sel_removable = FALSE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id == SET_MEDIUM) {                      /* media */
            xgc->par_sel_removable = FALSE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id == SET_OUT) {
            xgc->par_sel_removable = FALSE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id >= SET_POUT && id < SET_IOUT) {       /* point output */
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id >= SET_IIOUT && id < SET_COUT) {      /* image output */
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id >= SET_IOUT && id < SET_IIOUT) {      /* plane output */
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id >= SET_COUT && id < SET_SOUT) {       /* volume output */
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id >= SET_SOUT && id < SET_FOUT) {       /* sum output */
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id >= SET_FOUT && id < SET_GROW) {       /* force output */
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id >= SET_GROW && id < SET_ROUGHNESS) {    /* V. 2. Add growth modifier */
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id >= SET_ROUGHNESS && id < SET_SPECTRAL) { /* V. 3. Add roughness modifier */
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id >= SET_SPECTRAL && id < SET_EXPRESSION) { /* V. 4. Add spectral modifier */
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id >= SET_EXPRESSION && id < SET_NFAREA) {   /* V. 5. Add expression modifier */
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id == SET_NFFF) {                        /* VII. 1. Near field to far field transform box */
            xgc->par_sel_removable = FALSE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id == SET_PNFFF) {                       /* VIII. 1. Near field to far field transform box */
            xgc->par_sel_removable = FALSE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id >= SET_NFFFP && id < SET_PNFAREA) {   /* VII. 2. Near field to far field point */
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id >= SET_NFAREA && id < SET_NFFFP) {    /* VII. 3. Near field to far field area */
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id >= SET_PNFAREA && id < SET_PNFFFP) {  /* VIII. 3. Periodic near field to far field area */
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        } else if (id >= SET_PNFFFP) {                      /* VIII. 2. Periodic near field to far field point */
            xgc->par_sel_removable = TRUE;
            set_controls_visible(xgc, id);
            par_controls_sensitive(xgc, id);
        }
    }

    //gtk_widget_set_sensitive(xgc->par_button_remove, xgc->par_sel_removable);

    return TRUE;
}

void
par_row_activated_cb(GtkTreeView *treeview, GtkTreePath *path, GtkTreeViewColumn *col, XGControls *xgc)
{
    GtkTreeModel *model;
    GtkTreeIter   iter;

    model = gtk_tree_view_get_model(treeview);
    if (gtk_tree_model_get_iter(model, &iter, path)) {
        if (gtk_tree_view_row_expanded(treeview, path))
            gtk_tree_view_collapse_row(treeview, path);
        else
            gtk_tree_view_expand_row(treeview, path, TRUE);
    }
}

void
par_view_popup_menu_on_remove(GtkWidget *menuitem, XGControls *xgc)
{
    par_row_remove_process(xgc);
}

void
par_view_popup_menu(GtkWidget *treeview, GdkEventButton *event, XGControls *xgc)
{
    GtkWidget *menu, *menuitem;

    menu = gtk_menu_new();

    menuitem = gtk_menu_item_new_with_label("Remove");

    g_signal_connect(menuitem, "activate", G_CALLBACK(par_view_popup_menu_on_remove), xgc);

    gtk_widget_set_sensitive(menuitem, xgc->par_sel_removable);

    gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);

    gtk_widget_show_all(menu);

    /* Note: event can be NULL here when called from view_onPopupMenu;
    *  gdk_event_get_time() accepts a NULL argument */
    gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL,
                   (event != NULL) ? event->button : 0,
                   gdk_event_get_time((GdkEvent*)event));
}

gboolean
par_row_button_release_cb(GtkWidget* widget, GdkEventButton * event, XGControls *xgc)
{
    if ( (event->type == GDK_BUTTON_RELEASE) && (event->button == 3) ) {
#if(1)
        //////////////////////////////////////////////////////////////////////////
        //GtkTreePath *path;
        //gint id = xgc->par_selid;        
        //
        //if (xgc->par_selid == -1)
        //    return FALSE;
        //
        //path = gtk_tree_model_get_path(GTK_TREE_MODEL(xgc->ts_par), &(xgc->par_seliter));
        //row_delete_process(xgc, id, path);

        par_view_popup_menu(widget, event, xgc);
        //////////////////////////////////////////////////////////////////////////
#endif
#if(0)
        GtkTreeSelection *selection;
        selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(widget));

        if (gtk_tree_selection_count_selected_rows(selection) <= 1) {
            GtkTreePath *path;

            /* Get tree path for row that was clicked */
            if (gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(widget), (gint)event->x, (gint)event->y, &path, NULL, NULL, NULL))
            {
                gtk_tree_selection_unselect_all(selection);
                gtk_tree_selection_select_path(selection, path);

                

                gtk_tree_path_free(path);
            }
        }
#endif

        return TRUE;
    }

    return FALSE;
}


void
mat_row_activated_cb(GtkTreeView *treeview, GtkTreePath *path, GtkTreeViewColumn *col, XGControls *xgc)
{
    GtkTreeModel *model;
    GtkTreeIter   iter;

    model = gtk_tree_view_get_model(treeview);
    if (gtk_tree_model_get_iter(model, &iter, path)) {
        if (gtk_tree_view_row_expanded(treeview, path))
            gtk_tree_view_collapse_row(treeview, path);
        else
            gtk_tree_view_expand_row(treeview, path, TRUE);
    }

    /*gchar *text;
    gint id;
    GtkTreeIter  iter;

    gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_mat), &iter, arg1);
    gtk_tree_model_get(GTK_TREE_MODEL(xgc->ts_mat), &iter, COLUMN_PARAMETER, &text, COLUMN_ID, &id, -1);

    mat_row_edit_process(xgc, id, arg1);*/
}

static void
mat_row_remove_process(XGControls *xgc)
{
    gint index;
    gchar *text;
    GtkTreePath *path;
    gboolean writeme = FALSE;
    gint id = xgc->mat_selid;

    if (xgc->mat_selid == -1)
        return;

    path = gtk_tree_model_get_path(GTK_TREE_MODEL(xgc->ts_mat), &(xgc->mat_seliter));

    //////////////////////////////////////////////////////////////////////////
    GtkTreeIter iter;
    GtkTreePath *path_to_select;
    gchar* path_string_to_select = NULL;

    GtkTreeSelection *selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(xgc->tv_mat));

    if (gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_mat), &iter, path)) {
        gtk_tree_model_get(GTK_TREE_MODEL(xgc->ts_mat), &iter, COLUMN_PARAMETER, &text, COLUMN_ID, &id, -1);

        path_to_select = gtk_tree_model_get_path(GTK_TREE_MODEL(xgc->ts_mat), &iter);
        path_string_to_select = gtk_tree_path_to_string(path_to_select);
    }
    //////////////////////////////////////////////////////////////////////////    

    if (id >= SET_MAT_SPHERE && id < SET_MAT_VOXEL) {
        writeme = TRUE;
        xgc->data.set_mat.spheres = g_array_remove_index(xgc->data.set_mat.spheres, id - SET_MAT_SPHERE);
    } else if (id >= SET_MAT_VOXEL && id < SET_MAT_CYLINDER) {
        writeme = TRUE;
        xgc->data.set_mat.voxels = g_array_remove_index(xgc->data.set_mat.voxels, id - SET_MAT_VOXEL);
    } else if (id >= SET_MAT_CYLINDER && id < SET_MAT_CONE) {
        writeme = TRUE;
        xgc->data.set_mat.cylinders = g_array_remove_index(xgc->data.set_mat.cylinders, id - SET_MAT_CYLINDER);
    } else if (id >= SET_MAT_CONE && id < SET_MAT_RCONE) {
        writeme = TRUE;
        xgc->data.set_mat.cones = g_array_remove_index(xgc->data.set_mat.cones, id - SET_MAT_CONE);
    } else if (id >= SET_MAT_RCONE && id < SET_MAT_GWYDD) {
        writeme = TRUE;
        xgc->data.set_mat.rcones = g_array_remove_index(xgc->data.set_mat.rcones, id - SET_MAT_RCONE);
    }
    //    else if (id >= MAT_GWYDD && id<MAT_MESH) {
    //        writeme = TRUE; 
    //        xgc->data.gwydds = g_array_remove_index(xgc->data.gwydds, id-MAT_GWYDD);
    //    }
    //    else if (id >= MAT_MESH && id<MAT_TETRAHEDRON) 
    //        g_array_remove_index(xgc->data., id);
    else if (id >= SET_MAT_TETRAHEDRON) {
        writeme = TRUE;
        g_array_remove_index(xgc->data.set_mat.tetrahedrons, id - SET_MAT_TETRAHEDRON);
    }


    /*rewrite the temporary matfile and then reload it and interpret completely in order to assure consistency between its text and gui interpretation*/
    if (writeme) {
        if (write_matfile(xgc)) {
            path = gtk_tree_path_new_from_indices(0, -1);
            gtk_tree_view_expand_row(GTK_TREE_VIEW(xgc->tv_mat), path, TRUE);
        }
    }

    /*rewrite the temporary matfile and then reload it and interpret completely in order to assure consistency between its text and gui interpretation*/
    if (writeme) {
        if (write_matfile(xgc)) {
            mt_create_tree(xgc);
            index = gtk_tree_path_get_indices(path)[0];
            path = gtk_tree_path_new_from_indices(index, -1);
            gtk_tree_view_expand_row(GTK_TREE_VIEW(xgc->tv_mat), path, TRUE);

            selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(xgc->tv_mat));
            if (gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(xgc->ts_mat), &iter, path_string_to_select))
                gtk_tree_selection_select_iter(selection, &iter);
            else {
                gint n;
                GtkTreeIter iter_parent;

                if (gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_mat), &iter_parent, path)) {
                    n = gtk_tree_model_iter_n_children(GTK_TREE_MODEL(xgc->ts_mat), &iter);
                    if (gtk_tree_model_iter_nth_child (GTK_TREE_MODEL(xgc->ts_mat), &iter, &iter_parent, n-1)) {
                        gtk_tree_selection_select_iter(selection, &iter);
                    }
                    else
                        gtk_tree_selection_select_iter(selection, &iter_parent);
                }
            }

            //pc_row_select_process(xgc, id_to_select, path_to_select);
        }
    }
}

void
mat_view_popup_menu_on_remove(GtkWidget *menuitem, XGControls *xgc)
{
    mat_row_remove_process(xgc);
}

void
mat_view_popup_menu(GtkWidget *treeview, GdkEventButton *event, XGControls *xgc)
{
    GtkWidget *menu, *menuitem;

    menu = gtk_menu_new();

    menuitem = gtk_menu_item_new_with_label("Remove");

    g_signal_connect(menuitem, "activate", G_CALLBACK(mat_view_popup_menu_on_remove), xgc);

    gtk_widget_set_sensitive(menuitem, xgc->mat_sel_removable);

    gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);

    gtk_widget_show_all(menu);

    /* Note: event can be NULL here when called from view_onPopupMenu;
    *  gdk_event_get_time() accepts a NULL argument */
    gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL,
                   (event != NULL) ? event->button : 0,
                   gdk_event_get_time((GdkEvent*)event));
}

gboolean
mat_row_button_release_cb(GtkWidget* widget, GdkEventButton * event, XGControls *xgc)
{
    if ((event->type == GDK_BUTTON_RELEASE) && (event->button == 3)) {
        mat_view_popup_menu(widget, event, xgc);
        return TRUE;
    }

    return FALSE;
}

void mat_row_remove_cb(GtkButton *button, XGControls *xgc)
{    
    mat_row_remove_process(xgc);
}

gboolean mat_row_selected_cb(GtkTreeSelection *selection, XGControls *xgc)
{
    gint id;
    GtkTreeIter iter;
    gchar *text;
    GtkTreeModel *model = GTK_TREE_MODEL(xgc->ts_mat);
    GtkTreePath *path;
    xgc->mat_selid = -1;
    xgc->mat_sel_removable = FALSE;

    set_all_mat_controls_invisible(xgc);

    if (gtk_tree_selection_get_selected(selection, &model, &iter)) {
        gtk_tree_model_get(GTK_TREE_MODEL(xgc->ts_mat), &iter, COLUMN_PARAMETER, &text, COLUMN_ID, &id, -1);
        xgc->mat_selid = id;
        xgc->mat_seliter = iter;

        //
        if (xgc->mat_selid == -1)
            return FALSE;

        path = gtk_tree_model_get_path(GTK_TREE_MODEL(xgc->ts_mat), &(xgc->mat_seliter));
        mat_controls_row_select_process(xgc, id, path);
        //

        if (id >= SET_MAT_SPHERE && id < SET_MAT_VOXEL) {
            xgc->mat_sel_removable = TRUE;
            mat_set_controls_visible(xgc, id);
            mat_controls_sensitive(xgc, id);
        } else if (id >= SET_MAT_VOXEL && id < SET_MAT_CYLINDER) {
            xgc->mat_sel_removable = TRUE;
            mat_set_controls_visible(xgc, id);
            mat_controls_sensitive(xgc, id);
        } else if (id >= SET_MAT_CYLINDER && id < SET_MAT_CONE) {
            xgc->mat_sel_removable = TRUE;
            mat_set_controls_visible(xgc, id);
            mat_controls_sensitive(xgc, id);
        } else if (id >= SET_MAT_CONE && id < SET_MAT_RCONE) {
            xgc->mat_sel_removable = TRUE;
            mat_set_controls_visible(xgc, id);
            mat_controls_sensitive(xgc, id);
        } else if (id >= SET_MAT_RCONE && id < SET_MAT_GWYDD) {
            xgc->mat_sel_removable = TRUE;
            mat_set_controls_visible(xgc, id);
            mat_controls_sensitive(xgc, id);
        } else if (id >= SET_MAT_GWYDD && id < SET_MAT_MESH) {
            xgc->mat_sel_removable = TRUE;
            mat_set_controls_visible(xgc, id);
            mat_controls_sensitive(xgc, id);
        } else if (id >= SET_MAT_MESH && id < SET_MAT_TETRAHEDRON) {
            xgc->mat_sel_removable = TRUE;
            mat_set_controls_visible(xgc, id);
            mat_controls_sensitive(xgc, id);
        } else if (id >= SET_MAT_TETRAHEDRON) {
            xgc->mat_sel_removable = TRUE;
            mat_set_controls_visible(xgc, id);
            mat_controls_sensitive(xgc, id);
        }
    }

    return TRUE;
}

void
par_buffer_changed_cb(GtkTextBuffer *textbuffer, XGControls *xgc)
{
    //////////////////////////////////////////////////////////////////////////
    SvSet set;
    gchar msg_string[256];

    write_textbuffer_to_parfile(xgc);

    if (NULL == xgc->tmpfilename)
        xgc->tmpfilename = get_temporary_par_filename(xgc);

    if (xgc->tmpfilename) {
        clear_settings(&set, FALSE);

        if (parse_settings(xgc->tmpfilename, &set, FALSE)) {
            xgc->par_file_success = TRUE;

            gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, "");
        } else {
            xgc->par_file_success = FALSE;

            g_snprintf(msg_string, sizeof(msg_string), MSG_SB_ERROR_PARSING_PARFILE);
            gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, msg_string);
        }
    }
    //////////////////////////////////////////////////////////////////////////

    //pt_create_tree(xgc);
}

void
mat_buffer_changed_cb(GtkTextBuffer *textbuffer, XGControls *xgc)
{
    //mat_buffer_changed(xgc);

    SvSetMat set;
    gchar msg_string[256];

    write_textbuffer_to_matfile(xgc);

    if (NULL == xgc->tmpmatfilename)
        xgc->tmpmatfilename = get_temporary_mat_filename(xgc);

    if (xgc->tmpmatfilename) {
        clear_settings_mat(&set, FALSE);

        if (parse_settings_mat(xgc->tmpmatfilename, &xgc->data.set, &set, FALSE)) {
            xgc->mat_file_success = TRUE;

            gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, "");
        } else {
            xgc->mat_file_success = FALSE;

            g_snprintf(msg_string, sizeof(msg_string), MSG_SB_ERROR_PARSING_PARFILE);
            gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, msg_string);
        }
    }
}

void graph_changed_cb(GtkComboBox *combo, XGControls *xgc)
{
    update_graph(xgc);
}

void image_changed_cb(GtkComboBox *combo, XGControls *xgc)
{
    update_image(xgc);
}

void help_x_cb(GtkWidget *widget, XGControls *xgc)
{
    GtkWidget *dialog;
    gint ret = 0;


#ifdef G_OS_WIN32
    ret = (gint)ShellExecute(NULL, "open", "http://gsvit.net/wiki/doku.php/start:using_xsvit", NULL, NULL, SW_SHOWNORMAL);
    if (ret < 32)
        ret = 1;
    else
        ret = 0;
#endif

#ifndef G_OS_WIN32
    ret = system("xdg-open http://gsvit.net/xdocs17.php");

#endif


    if (ret != 0) {
        dialog = gtk_message_dialog_new_with_markup(GTK_WINDOW(xgc->toplevel),
                                                    GTK_DIALOG_DESTROY_WITH_PARENT,
                                                    GTK_MESSAGE_INFO,
                                                    GTK_BUTTONS_CLOSE,
                                                    "Error: cannot open link in a web browser.");

        gtk_dialog_run(GTK_DIALOG(dialog));
        gtk_widget_destroy(dialog);
    }
}

void help_g_cb(GtkWidget *widget, XGControls *xgc)
{
    GtkWidget *dialog;
    gint ret = 0;


#ifdef G_OS_WIN32
    ret = (gint)ShellExecute(NULL, "open", "http://gsvit.net/docs17.php", NULL, NULL, SW_SHOWNORMAL);
    if (ret < 32)
        ret = 1;
    else
        ret = 0;
#endif

#ifndef G_OS_WIN32
    ret = system("xdg-open http://gsvit.net/docs17.php");

#endif


    if (ret != 0) {
        dialog = gtk_message_dialog_new_with_markup(GTK_WINDOW(xgc->toplevel),
                                                    GTK_DIALOG_DESTROY_WITH_PARENT,
                                                    GTK_MESSAGE_INFO,
                                                    GTK_BUTTONS_CLOSE,
                                                    "Error: cannot open link in a web browser.");

        gtk_dialog_run(GTK_DIALOG(dialog));
        gtk_widget_destroy(dialog);
    }
}

void about_cb(GtkWidget *widget, XGControls *xgc)
{
    gchar buff[256];
    g_snprintf(buff, sizeof(buff), "Version: %s\n<b>XSvit</b> is a GUI for open source FDTD solver GSvit.", VERSION);

    GtkWidget *dialog = gtk_message_dialog_new_with_markup(GTK_WINDOW(xgc->toplevel),
        GTK_DIALOG_DESTROY_WITH_PARENT,
        GTK_MESSAGE_INFO,
        GTK_BUTTONS_CLOSE,
        buff);

    gtk_message_dialog_format_secondary_markup(GTK_MESSAGE_DIALOG(dialog), "Developed by P. Klapetek, various OS porting by M. Valtr.");
    gtk_window_set_title(GTK_WINDOW(dialog), "About XSvit");
    gtk_dialog_run(GTK_DIALOG(dialog));
    gtk_widget_destroy(dialog);
}



/* SIGNALS */

gboolean on_key_press_par_tree_cb(GtkWidget *widget, GdkEventKey *event, XGControls *xgc)
{
    switch (event->keyval)
    {
    case GDK_Delete:
        par_row_remove_process(xgc);
        break;
    default:
        return FALSE;
    }

    return TRUE;
}

gboolean on_key_press_mat_tree_cb(GtkWidget *widget, GdkEventKey *event, XGControls *xgc)
{
    switch (event->keyval)
    {
    case GDK_Delete:
        mat_row_remove_process(xgc);
        break;
    default:
        return FALSE;
    }

    return TRUE;
}

static void
set_all_controls_invisible(XGControls *xgc)
{
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), FALSE);

    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.cc.vbg_outer_comp_domain), FALSE);

    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.bac.vbg_outer_basic), FALSE);

    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_box), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_point_origin), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_source), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_ia), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_fs), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_ls), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_zmultiplier), FALSE);

    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.boc.vbg_outer_cvboundary), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.boc.vbg_outer_periodicboundary), FALSE);

    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.mc.vbg_outer_material_prop), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.mc.vbg_outer_material_grow), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.mc.vbg_outer_material_rough), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.mc.vbg_outer_material_spectral), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.mc.vbg_outer_material_expr), FALSE);

    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_general_output), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_point_output), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_image_output), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_plane_output), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_volume_output), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_sum_output), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_box_output), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_force_output), FALSE);

    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.fc.vbg_outer_nfff), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.fc.vbg_outer_nfffp), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.fc.vbg_outer_nfffa), FALSE);

    gtk_widget_set_visible(GTK_WIDGET(xgc->cs.fc.vbg_outer_pnfff), FALSE);
}


static void
set_all_mat_controls_invisible(XGControls *xgc)
{
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.label_name), FALSE);

    gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_sphere), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_box), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_cylinder), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_cone), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_rcone), FALSE);
    gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_tetrahedron), FALSE);

    gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_prop), FALSE);
}


static void
set_controls_visible(XGControls *xgc, gint id)
{
    if (id == SET_POOL) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Computational domain");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.cc.vbg_outer_comp_domain, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.cc.vbg_outer_comp_domain), TRUE);
    } else if (id == SET_BASIC) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Basic parameters");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.bac.vbg_outer_basic, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.bac.vbg_outer_basic), TRUE);
    } else if (id >= SET_PSOURCE && id < SET_POUT) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Point source");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_point_origin, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_source, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_point_origin), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_source), TRUE);
    } else if (id == SET_SF) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Scattered field source (SF)");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_source, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_ia, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_source), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_ia), TRUE);
    } else if (id == SET_TSF) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Total / Scattered field source (TSF)");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_box, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_source, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_ia, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_zmultiplier, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_box), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_source), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_ia), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_zmultiplier), TRUE);
    } else if (id == SET_TSFF) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Focused Total / Scattered field source (TSFF)");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_box, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_source, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_fs, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_box), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_source), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_fs), TRUE);
    } else if (id == SET_LTSF) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Layered Total / Scattered field source (LTSF)");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_box, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_source, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_ia, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_ls, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_zmultiplier, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_box), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_source), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_ia), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_ls), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_zmultiplier), TRUE);
    } else if (id == SET_LTSFF) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Layered Focused Total / Scattered field source (LTSFF)");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_box, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_source, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_fs, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.sc.vbg_outer_ls, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_box), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_source), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_fs), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.sc.vbg_outer_ls), TRUE);
    } else if (id == SET_OUT) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "General output");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.oc.vbg_outer_general_output, -1);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_general_output), TRUE);
    } else if (id >= SET_POUT && id < SET_IOUT) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Point output");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.oc.vbg_outer_point_output, -1);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_point_output), TRUE);
    } else if (id >= SET_IIOUT && id < SET_COUT) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Image output");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.oc.vbg_outer_image_output, -1);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_image_output), TRUE);
    } else if (id >= SET_IOUT && id < SET_IIOUT) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Plane output");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.oc.vbg_outer_plane_output, -1);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_plane_output), TRUE);
    } else if (id >= SET_COUT && id < SET_SOUT) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Volume output");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.oc.vbg_outer_volume_output, -1);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_volume_output), TRUE);
    } else if (id >= SET_SOUT && id < SET_FOUT) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Sum output");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.oc.vbg_outer_box_output, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.oc.vbg_outer_sum_output, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_box_output), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_sum_output), TRUE);
    } else if (id >= SET_FOUT && id < SET_GROW) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Force output");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.oc.vbg_outer_box_output, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.oc.vbg_outer_force_output, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_box_output), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.oc.vbg_outer_force_output), TRUE);
    } else if (id == SET_BND) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Boundary conditions");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.boc.vbg_outer_cvboundary, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.boc.vbg_outer_periodicboundary, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.boc.vbg_outer_cvboundary), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.boc.vbg_outer_periodicboundary), TRUE);
    } else if (id == SET_MEDIUM) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Material properties");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.mc.vbg_outer_material_prop, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.mc.vbg_outer_material_prop), TRUE);
    } else if (id >= SET_GROW && id < SET_ROUGHNESS) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Growth material modifier");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.mc.vbg_outer_material_grow, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.mc.vbg_outer_material_grow), TRUE);
    } else if (id >= SET_ROUGHNESS && id < SET_SPECTRAL) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Roughness material modifier");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.mc.vbg_outer_material_rough, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.mc.vbg_outer_material_rough), TRUE);
    } else if (id >= SET_SPECTRAL && id < SET_EXPRESSION) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Spectral material modifier");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.mc.vbg_outer_material_spectral, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.mc.vbg_outer_material_spectral), TRUE);
    } else if (id >= SET_EXPRESSION && id < SET_NFAREA) {
        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Expression material modifier");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.mc.vbg_outer_material_expr, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.mc.vbg_outer_material_expr), TRUE);
    } else if (id == SET_NFFF) {
        /* VII. 1. Near field to far field transform */

        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Near field to far field transform");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.fc.vbg_outer_nfff, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.fc.vbg_outer_nfff), TRUE);
    } else if (id >= SET_NFFFP && id < SET_PNFAREA) {
        /* VII. 2. Near field to far field point */

        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Far field point");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.fc.vbg_outer_nfffp, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.fc.vbg_outer_nfffp), TRUE);
    } else if (id >= SET_NFAREA && id < SET_NFFFP) {
        /* VII. 3. Near field to far field area */

        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Far field area output");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.fc.vbg_outer_nfffa, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.fc.vbg_outer_nfffa), TRUE);
    } else if (id == SET_PNFFF) {
        /* VIII. 1. Periodic near field to far field transform */

        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Periodic near field to far field transform");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.fc.vbg_outer_pnfff, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.fc.vbg_outer_pnfff), TRUE);
    } else if (id >= SET_PNFAREA && id < SET_PNFFFP) {
        /* VIII. 3. Periodic near field to far field area */

        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Periodic far field area output");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.fc.vbg_outer_nfffa, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.fc.vbg_outer_nfffa), TRUE);
    } else if (id >= SET_PNFFFP) {
        /* VIII. 2. Periodic near field to far field point */

        gtk_label_set_text(GTK_LABEL(xgc->cs.label_name), "Periodic far field point");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs.vbox_page), xgc->cs.fc.vbg_outer_nfffp, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs.fc.vbg_outer_nfffp), TRUE);
    }
}

static void
mat_set_controls_visible(XGControls *xgc, gint id)
{
    if (id >= SET_MAT_SPHERE && id < SET_MAT_VOXEL) {
        gtk_label_set_text(GTK_LABEL(xgc->cs_mat.label_name), "Sphere");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs_mat.vbox_page), xgc->cs_mat.oc.vbg_outer_sphere, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs_mat.vbox_page), xgc->cs_mat.oc.vbg_outer_prop, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_sphere), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_prop), TRUE);
    } else if (id >= SET_MAT_VOXEL && id < SET_MAT_CYLINDER) {
        gtk_label_set_text(GTK_LABEL(xgc->cs_mat.label_name), "Box");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs_mat.vbox_page), xgc->cs_mat.oc.vbg_outer_box, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs_mat.vbox_page), xgc->cs_mat.oc.vbg_outer_prop, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_box), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_prop), TRUE);
    } else if (id >= SET_MAT_CYLINDER && id < SET_MAT_CONE) {
        gtk_label_set_text(GTK_LABEL(xgc->cs_mat.label_name), "Cylinder");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs_mat.vbox_page), xgc->cs_mat.oc.vbg_outer_cylinder, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs_mat.vbox_page), xgc->cs_mat.oc.vbg_outer_prop, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_cylinder), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_prop), TRUE);
    } else if (id >= SET_MAT_CONE && id < SET_MAT_RCONE) {
        gtk_label_set_text(GTK_LABEL(xgc->cs_mat.label_name), "Cone");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs_mat.vbox_page), xgc->cs_mat.oc.vbg_outer_cone, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs_mat.vbox_page), xgc->cs_mat.oc.vbg_outer_prop, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_cone), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_prop), TRUE);
    } else if (id >= SET_MAT_RCONE && id < SET_MAT_GWYDD) {
        gtk_label_set_text(GTK_LABEL(xgc->cs_mat.label_name), "Cut cone");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs_mat.vbox_page), xgc->cs_mat.oc.vbg_outer_rcone, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs_mat.vbox_page), xgc->cs_mat.oc.vbg_outer_prop, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_rcone), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_prop), TRUE);
    } else if (id >= SET_MAT_GWYDD && id < SET_MAT_MESH) {
    } else if (id >= SET_MAT_MESH && id < SET_MAT_TETRAHEDRON) {
    } else if (id >= SET_MAT_TETRAHEDRON) {
        gtk_label_set_text(GTK_LABEL(xgc->cs_mat.label_name), "Tetrahedron");
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.label_name), TRUE);

        gtk_box_reorder_child(GTK_BOX(xgc->cs_mat.vbox_page), xgc->cs_mat.oc.vbg_outer_tetrahedron, -1);
        gtk_box_reorder_child(GTK_BOX(xgc->cs_mat.vbox_page), xgc->cs_mat.oc.vbg_outer_prop, -1);

        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_tetrahedron), TRUE);
        gtk_widget_set_visible(GTK_WIDGET(xgc->cs_mat.oc.vbg_outer_prop), TRUE);
    }
}

#if(0)
void cvbtype0_changed_cb(GtkComboBox *combo, XGControls *xgc)
{
    if (gtk_combo_box_get_active(combo) == 3)
        set_boundary_controls_sensitive(xgc, 0, TRUE);
    else
        set_boundary_controls_sensitive(xgc, 0, FALSE);

    par_controls_changed(xgc);
}
void cvbtype1_changed_cb(GtkComboBox *combo, XGControls *xgc)
{
    if (gtk_combo_box_get_active(combo) == 3)
        set_boundary_controls_sensitive(xgc, 1, TRUE);
    else
        set_boundary_controls_sensitive(xgc, 1, FALSE);

    par_controls_changed(xgc);
}
void cvbtype2_changed_cb(GtkComboBox *combo, XGControls *xgc)
{
    if (gtk_combo_box_get_active(combo) == 3)
        set_boundary_controls_sensitive(xgc, 2, TRUE);
    else
        set_boundary_controls_sensitive(xgc, 2, FALSE);

    par_controls_changed(xgc);
}
void cvbtype3_changed_cb(GtkComboBox *combo, XGControls *xgc)
{
    if (gtk_combo_box_get_active(combo) == 3)
        set_boundary_controls_sensitive(xgc, 3, TRUE);
    else
        set_boundary_controls_sensitive(xgc, 3, FALSE);

    par_controls_changed(xgc);
}
void cvbtype4_changed_cb(GtkComboBox *combo, XGControls *xgc)
{
    if (gtk_combo_box_get_active(combo) == 3)
        set_boundary_controls_sensitive(xgc, 4, TRUE);
    else
        set_boundary_controls_sensitive(xgc, 4, FALSE);

    par_controls_changed(xgc);
}
void cvbtype5_changed_cb(GtkComboBox *combo, XGControls *xgc)
{
    if (gtk_combo_box_get_active(combo) == 3)
        set_boundary_controls_sensitive(xgc, 5, TRUE);
    else
        set_boundary_controls_sensitive(xgc, 5, FALSE);

    par_controls_changed(xgc);
}
/* boundary sensitive functions */
#endif

static void 
set_computation_controls_sensitive(XGControls *xgc)
{
    gboolean running = (xgc->data.pid != 0);

    gtk_widget_set_sensitive(GTK_WIDGET(xgc->run_tool_item), !running);
    gtk_widget_set_sensitive(GTK_WIDGET(xgc->run_menu_item), !running);

    gtk_widget_set_sensitive(GTK_WIDGET(xgc->stop_tool_item), running);
    gtk_widget_set_sensitive(GTK_WIDGET(xgc->stop_menu_item), running);

    gtk_widget_set_sensitive(GTK_WIDGET(xgc->new_menu_item), !running);
    gtk_widget_set_sensitive(GTK_WIDGET(xgc->open_menu_item), !running);
    gtk_widget_set_sensitive(GTK_WIDGET(xgc->save_menu_item), !running);
    gtk_widget_set_sensitive(GTK_WIDGET(xgc->saveas_menu_item), !running);
    gtk_widget_set_sensitive(GTK_WIDGET(xgc->savematas_menu_item), !running);
    gtk_widget_set_sensitive(GTK_WIDGET(xgc->recentfiles_menu_item), !running);
}

/*
static void
fill_material_combo(GtkWidget **combo)
{
    GDir *dir;
    GError *error;
    const gchar *filename;
    gchar * path;
     
    path = get_spectra_path(NULL);

    dir = g_dir_open(path, 0, &error);
    //gtk_list_store_insert_with_values(*combo_list, NULL, -1, 0, SOURCE_LAYERED_MATERIAL, -1);
    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(*combo), MAT_PROP_MATERIAL_NONE);
    while ((filename = g_dir_read_name(dir))) {
        //printf("%s\n", filename);
        //gtk_list_store_insert_with_values(*combo_list, NULL, -1, 0, g_strdup(filename), -1);

        if (FALSE == (strcmp(filename, "Makefile.am") == 0 || strcmp(filename, "makefile.am") == 0))
            gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(*combo), g_strdup(filename));
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(*combo), MAT_PROP_MATERIAL_NONE_INDEX);

    g_dir_close(dir);
}
*/

/*  Properties Dialog Page Layout Conception

|-------------------------------------------------------|
|vbox_page                                              |
|                                                       |
||-----------------------------------------------------||
||vbg_outer_xxx (vbox group outer xxx)                 ||
||                                                     ||
|||---------------------------------------------------|||
|||hseparator_xxx (horizontal separator - optional)   |||
|||                                                   |||
|||hbg_xxx (hbox group)                               |||
|||                                                   |||
|||---------------------------------------------------|||
||||vbg_inner_xxx (vbox group inner xxx)             ||||
||||                                                 ||||
||||label_xxx (group label - optional                ||||
||||                                                 ||||
||||table_xxx                                        ||||
|||||---|---|---|---|                                ||||
|||||---|---|---|---|                                ||||
|||||---|---|---|---|                                ||||
||||-------------------------------------------------||||
|||---------------------------------------------------|||
||-----------------------------------------------------||
|-------------------------------------------------------|

*/

void create_vbox_group(GtkWidget **vbg_outer, 
                       GtkWidget **table,
                       const GtkWidget *vbox_page,
                       const gchar *label_text,
                       const gint rows, 
                       const gint columns)
{
    GtkWidget *hbg, *hseparator, *vbg_inner, *label;

    *vbg_outer = gtk_vbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox_page), *vbg_outer, TRUE, TRUE, 4);    

    if (NULL != label_text) {
        hseparator = gtk_hseparator_new();
        gtk_box_pack_start(GTK_BOX(*vbg_outer), hseparator, TRUE, TRUE, 4);
    }

    hbg = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(*vbg_outer), hbg, TRUE, TRUE, 4);
    
    vbg_inner = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbg_inner), 4);
    gtk_box_pack_start(GTK_BOX(hbg), vbg_inner, FALSE, FALSE, 0);

    if (NULL != label_text) {
        label = gtk_label_new(label_text);
        gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
        gtk_box_pack_start(GTK_BOX(vbg_inner), label, TRUE, TRUE, 4);
    }

    *table = gtk_table_new(rows, columns, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(*table), 2);
    gtk_table_set_col_spacings(GTK_TABLE(*table), 6);
    gtk_box_pack_start(GTK_BOX(vbg_inner), *table, TRUE, TRUE, 4);
}


int main(int argc, char *argv[])
{
    XGControls xgc;
//    GdkGLConfig *glconfig;

    GtkWidget *window;
    GtkWidget *vbox, *vbox2_left, *vbox2_right, *hbox, *hbox2, *page, *page_par_tree, *page_mat_tree, *vbox2_graph, *graph, *frame, *vbox2bottom_right;
    GtkWidget *view_par_tree, *view_par_file, *view_mat_file, *view_scene;
    GtkWidget *scroll_window_par_tree, *scroll_window_par_values, *scroll_window_par_file, *scroll_window_mat_tree, *scroll_window_mat_values, *scroll_window_mat_file, *scroll_window_output;
    GtkWidget *vbox3, *vbox_par, *vbox_mat;
    GtkCellRenderer *renderer;
    gint col_offset;
    GtkTreeViewColumn *column;
    GtkTreeSelection *selection;

    GtkWidget *menubar = NULL;
    GtkWidget *toolbar = NULL;
    GtkWidget *notebook;
    gchar *filename = NULL;;

    GKeyFileFlags flags;
    GError *error = NULL;

    /* set locale */
    setlocale(LC_ALL, "en_US.utf-8");
    gtk_disable_setlocale();

#if GLIB_MINOR_VERSION < 32
    g_thread_init(NULL);
#endif
#if GLIB_MINOR_VERSION < 36
    g_type_init();
#endif

    gtk_init(&argc, &argv);
    gtk_gl_init(&argc, &argv);
    gwy_process_type_init();


#ifdef G_OS_WIN32
    xgc.config_path = g_build_filename(g_get_home_dir(), ".xsvit", "xsvit.ini", NULL);
#endif

#ifndef G_OS_WIN32
    xgc.config_path = g_build_filename(g_get_home_dir(), ".xsvit", "xsvit.ini", NULL);
#endif


    // Create a new GKeyFile object and a bitwise list of flags. 
    xgc.keyfile = g_key_file_new();
    flags = G_KEY_FILE_KEEP_COMMENTS | G_KEY_FILE_KEEP_TRANSLATIONS;


    // Load the GKeyFile from keyfile.conf or return. 
    if (!g_key_file_load_from_file(xgc.keyfile, xgc.config_path, flags, &error)) {
        printf("No configuration file found (%s)\n", error->message);

        filename = g_build_filename(g_get_home_dir(), ".xsvit", NULL);
        if (!g_file_test(filename, G_FILE_TEST_EXISTS))
            g_mkdir(filename, 0777);
        g_free(filename);

        xgc.data.gsvit_location = NULL;
        xgc.data.gwyddion_location = NULL;
    } else {
        xgc.data.gsvit_location = g_key_file_get_string(xgc.keyfile, "Locations", "GSvit", NULL);
        xgc.data.gwyddion_location = g_key_file_get_string(xgc.keyfile, "Locations", "Gwyddion", NULL);
    }


    if (argc > 1) {
        if (!strcmp(argv[1], "--help") || gwy_strequal(argv[1], "-h")) {
            printf("XSvit is visualiser and editor of parameter files for GSvit open source FDTD solver\nUsage: ./xsvit or ./xsvit parameter_file\n");
            exit(0);
        }

        if (!strcmp(argv[1], "--version") || gwy_strequal(argv[1], "-v")) {
            printf("%s %s\n", PACKAGE_NAME, PACKAGE_VERSION);
            exit(0);
        }

        xgc.loadfile = g_strdup(argv[1]);
        //printf("Loading file %s\n", xgc.loadfile);
    } else
        xgc.loadfile = NULL;


    /*************   basic init  ************************/
    xgc.undo = 1;
    xgc.data.parfilename = NULL;
    xgc.data.matfilename = NULL;
    xgc.curdir = NULL;
    xgc.tmpfilename = NULL;
    xgc.tmpmatfilename = NULL;
    xgc.data.set_mat.gwyddata = NULL;
    xgc.container = NULL;
    xgc.xrot = xgc.yrot = 0;
    //xgc.data.set.sp.xres = xgc.data.set.sp.yres = xgc.data.set.sp.zres = 0;

    xgc.data.set.sp.xres = xgc.data.set.sp.yres = xgc.data.set.sp.zres = COMP_DOMAIN_SIZE;
    xgc.data.set.sp.dx = xgc.data.set.sp.dy = xgc.data.set.sp.dz = COMP_DOMAIN_SPACING;
    xgc.data.set.sc.nsteps = BASIC_STEPS;
    xgc.data.set.sc.verbose = BASIC_VERBOSE;
    xgc.data.set.sc.nthreads = BASIC_NTHREADS;
    xgc.data.set.sc.usegpu = BASIC_USEGPU;
    memset(xgc.data.set.sc.ugpu, 0, sizeof(gboolean)*MAX_GPUS);

    xgc.data.set_mat.spheres = g_array_new (FALSE, FALSE, sizeof (SvSphere));
    xgc.data.set_mat.cones = g_array_new (FALSE, FALSE, sizeof (SvCone));
    xgc.data.set_mat.rcones = g_array_new (FALSE, FALSE, sizeof (SvRCone));
    xgc.data.set_mat.cylinders = g_array_new (FALSE, FALSE, sizeof (SvCylinder));
    xgc.data.set_mat.voxels = g_array_new (FALSE, FALSE, sizeof (SvVoxel));
    xgc.data.set_mat.tetrahedrons = g_array_new (FALSE, FALSE, sizeof (SvTetrahedron));
    xgc.data.set_mat.gwydds = g_array_new (FALSE, FALSE, sizeof (SvGwydd));
    xgc.data.set_mat.ntetgens = 0;
    xgc.data.pid = 0;
    xgc.nfilestoshow = 0;
    xgc.timeout_func = 0;
    xgc.watch = 0;
    xgc.par_lock_controls_update = FALSE;    
    xgc.par_selid = -1;
    xgc.par_sel_removable = FALSE;
    xgc.par_file_success = FALSE;
    xgc.par_tree_success = FALSE;
    xgc.mat_lock_controls_update = FALSE;
    xgc.mat_selid = -1;
    xgc.mat_sel_removable = FALSE;
    xgc.mat_file_success = FALSE;
    xgc.mat_tree_success = FALSE;

    xgc.toplevel = window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_position(GTK_WINDOW(window), GTK_WIN_POS_CENTER);
    gtk_window_set_default_size(GTK_WINDOW(window), 1024, 600);

    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(window), vbox);

    xgc.data.is_nfff_point = (gboolean *)g_malloc(100000 * sizeof(gboolean));
    xgc.data.is_pnfff_point = (gboolean *)g_malloc(100000 * sizeof(gboolean));
    xgc.data.is_tetrahedron = (gboolean *)g_malloc(1000000 * sizeof(gboolean));

    /*xgc.data.sfb.ps_source_filename = NULL;
    xgc.data.sfb.sf_source_filename = NULL;
    xgc.data.sfb.tsf_source_filename = NULL;
    xgc.data.sfb.tsff_source_filename = NULL;
    xgc.data.sfb.ltsf_source_filename = NULL;
    xgc.data.sfb.ltsff_source_filename = NULL;*/

    reset_par_elemets_visibility(&xgc);
    

    /**************    main menubar and toolbar  *****************************/

    create_main_menubar_and_toolbar(&xgc, &menubar, &toolbar);
    gtk_box_pack_start(GTK_BOX(vbox), menubar, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(vbox), toolbar, FALSE, FALSE, 0);
       
    /********** data handling and visualisation ******************************************/
    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    /* left side */
    vbox2_left = gtk_vbox_new(FALSE, 0);
    //gtk_box_pack_start(GTK_BOX(hbox), vbox2_left, TRUE, TRUE, 5);


    /* parameter view notebook */
    vbox_par = gtk_vbox_new(FALSE, 0);
    //gtk_box_pack_start(GTK_BOX(vbox2_left), vbox_par, TRUE, TRUE, 5);


    xgc.param_notebook = notebook = gtk_notebook_new();
    //gtk_widget_set_size_request(notebook, 400, 350);
    gtk_box_pack_start(GTK_BOX(vbox_par), notebook, TRUE, TRUE, 0);

    /* parameter tree page */
    xgc.page_par_tree = page_par_tree = gtk_hbox_new(FALSE, 0);    
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page_par_tree, gtk_label_new(("Parameters")));

    //xgc.ts_par = gtk_tree_store_new(4, G_TYPE_STRING, G_TYPE_BOOLEAN, G_TYPE_INT, G_TYPE_BOOLEAN);
    xgc.ts_par = gtk_tree_store_new(NUM_COLUMNS, G_TYPE_STRING, G_TYPE_BOOLEAN, G_TYPE_INT, G_TYPE_BOOLEAN);

    xgc.tv_par = view_par_tree = gtk_tree_view_new_with_model(GTK_TREE_MODEL(xgc.ts_par));

    g_signal_connect(xgc.tv_par, "row-activated", G_CALLBACK(par_row_activated_cb), NULL);
    g_signal_connect(xgc.tv_par, "button-release-event", G_CALLBACK(par_row_button_release_cb), &xgc);
    g_signal_connect(xgc.tv_par, "key_press_event", G_CALLBACK (on_key_press_par_tree_cb), &xgc);

    memset(xgc.tv_par_expanded, FALSE, sizeof(xgc.tv_par_expanded));
    xgc.tv_par_path_select_string = FALSE;

    memset(xgc.tv_mat_expanded, FALSE, sizeof(xgc.tv_mat_expanded));
    xgc.tv_mat_path_select_string = FALSE;

    col_offset = gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(view_par_tree),
                                                             -1, "Parameter",
                                                             gtk_cell_renderer_text_new(),
                                                             "text", COLUMN_PARAMETER,
                                                             NULL);

    column = gtk_tree_view_get_column(GTK_TREE_VIEW(view_par_tree), col_offset - 1);
    gtk_tree_view_column_set_expand(GTK_TREE_VIEW_COLUMN(column), TRUE);

    renderer = gtk_cell_renderer_toggle_new();
    g_signal_connect(renderer, "toggled", G_CALLBACK(visible_toggle), &xgc);

    col_offset = gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(view_par_tree),
                                                             -1, "Visible",
                                                             renderer,
                                                             "active", COLUMN_CHECK,
                                                             "visible", COLUMN_SHOW_TOGGLE,
                                                             //"indicator-size", 3,
                                                             NULL);

    column = gtk_tree_view_get_column(GTK_TREE_VIEW(view_par_tree), col_offset - 1);
    gtk_tree_view_column_set_sizing(GTK_TREE_VIEW_COLUMN(column), GTK_TREE_VIEW_COLUMN_FIXED);
    gtk_tree_view_column_set_fixed_width(GTK_TREE_VIEW_COLUMN(column), 50);
    gtk_tree_view_column_set_alignment(GTK_TREE_VIEW_COLUMN(column), 1.0);
    gtk_tree_view_column_set_clickable(GTK_TREE_VIEW_COLUMN(column), TRUE);

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(view_par_tree));
    g_signal_connect(selection, "changed", G_CALLBACK(par_row_selected_cb), &xgc);

    scroll_window_par_tree = gtk_scrolled_window_new(NULL, NULL);
    //gtk_box_pack_start(GTK_BOX(page_par_tree), scroll_window_par_tree, TRUE, TRUE, 4);
    gtk_container_add(GTK_CONTAINER(scroll_window_par_tree), view_par_tree);    

    //////////////////////////////////////////////////////////////////////////        

    GtkWidget *table_comp_domain, *table_basic, *table_cvboundary, *table_periodicboundary;
    GtkWidget *table_box, *table_point_origin, *table_source, *table_ia, *table_fs, *table_zmultiplier;
    GtkWidget *table_general_output, *table_point_output, *table_image_output, *table_plane_output, *table_volume_output;
    GtkWidget *table_box_output, *table_sum_output, *table_force_output;
    GtkWidget *table_material_prop, *table_material_grow, *table_material_rough, *table_material_spectral, *table_material_expr;
    GtkWidget *table_nfff, *table_nfffp, *table_nfffa, *table_pnfff;

    GtkWidget *table_mat_sphere, *table_mat_box, *table_mat_cyl, *table_mat_cone, *table_mat_rcone, *table_mat_tthn, *table_mat_prop;

    GtkWidget *align, *label;
    GtkWidget *vbox_page, *treeview_packet_layered;
    GtkWidget *hpaned, *vpaned;

    /* vbox group and hbox group */
    GtkWidget *vbg_outer_ls, *vbg_inner_ls, *hbg_ls;   

    GtkWidget *hbox_point_origin_pos, *hbox_point_output_pos, *hbox_nfffp_pos;

    gint row_comp_domain, row_basic, row_cvboundary, row_periodicboundary, row_material_prop, row_material_grow, row_material_rough, row_material_spectral, row_material_expr;
    gint row_box, row_point_origin, row_source, row_ia, row_fs, row_zmultiplier;
    gint row_general_output, row_point_output, row_image_output, row_plane_output, row_volume_output, row_box_output, row_sum_output, row_force_output;
    gint row_nfff, row_nfffp, row_nfffa, row_pnfff;

    gint row_mat_sphere, row_mat_box, row_mat_cyl, row_mat_cone, row_mat_rcone, row_mat_tthn, row_mat_prop;

    gint cvboundary_typeval[6] = {CVBOUNDARY_TYPE, CVBOUNDARY_TYPE, CVBOUNDARY_TYPE, CVBOUNDARY_TYPE, CVBOUNDARY_TYPE, CVBOUNDARY_TYPE};
    gint cvboundary_depthval[6] = {CVBOUNDARY_DEPTH, CVBOUNDARY_DEPTH, CVBOUNDARY_DEPTH, CVBOUNDARY_DEPTH, CVBOUNDARY_DEPTH, CVBOUNDARY_DEPTH};
    gint cvboundary_mval[6] = {CVBOUNDARY_M, CVBOUNDARY_M, CVBOUNDARY_M, CVBOUNDARY_M, CVBOUNDARY_M, CVBOUNDARY_M};
    gdouble cvboundary_sigmaval[6] = {CVBOUNDARY_SIGMA, CVBOUNDARY_SIGMA, CVBOUNDARY_SIGMA, CVBOUNDARY_SIGMA, CVBOUNDARY_SIGMA, CVBOUNDARY_SIGMA};
    gdouble cvboundary_aval[6] = {CVBOUNDARY_A, CVBOUNDARY_A, CVBOUNDARY_A, CVBOUNDARY_A, CVBOUNDARY_A, CVBOUNDARY_A};
    gdouble cvboundary_kappaval[6] = {CVBOUNDARY_KAPPA, CVBOUNDARY_KAPPA, CVBOUNDARY_KAPPA, CVBOUNDARY_KAPPA, CVBOUNDARY_KAPPA, CVBOUNDARY_KAPPA};
    gboolean pboundary_mbval[6] = {CVBOUNDARY_MB, CVBOUNDARY_MB, CVBOUNDARY_MB, CVBOUNDARY_MB, CVBOUNDARY_MB, CVBOUNDARY_MB};
    gint pboundary_mbposval[6] = {CVBOUNDARY_MBPOSX0, CVBOUNDARY_MBPOSXN, CVBOUNDARY_MBPOSY0, CVBOUNDARY_MBPOSYN, CVBOUNDARY_MBPOSZ0, CVBOUNDARY_MBPOSZN};
    gint material_grow_inval = xgc.data.set.sp.xres - MATERIAL_GROW_INBORDER;
    gint material_grow_jnval = xgc.data.set.sp.yres - MATERIAL_GROW_JNBORDER;
    gint material_grow_knval = xgc.data.set.sp.zres - MATERIAL_GROW_KNBORDER;
    gint i;    

    /* properties page strip */
    xgc.cs.vbox_page = vbox_page = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox_page), 4);

    xgc.cs.sw_values = scroll_window_par_values = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scroll_window_par_values), GTK_SHADOW_NONE);
    //gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll), GTK_POLICY_NEVER, GTK_POLICY_ALWAYS);
    gtk_container_set_border_width(GTK_CONTAINER(scroll_window_par_values), 0);

    //gtk_box_pack_start(GTK_BOX(page_par_tree), scroll_window_par_values, TRUE, TRUE, 4);

    //gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroll_window), table);

    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_container_add(GTK_CONTAINER(align), vbox_page);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroll_window_par_values), align);

    gtk_widget_set_size_request(scroll_window_par_values, 50, -1);

    //xgc.cs.sc.label_name = gtk_label_new("Total/Scattered field source (TSF)");
    xgc.cs.label_name = gtk_label_new("");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.label_name), 0.0, 0.5);
    //gtk_table_attach(GTK_TABLE(table), xgc.cs.sc.label_name, 0, 1, row, row+1, GTK_EXPAND, 0, 0, 0);
    gtk_box_pack_start(GTK_BOX(vbox_page), xgc.cs.label_name, TRUE, TRUE, 4);    

    /* I. Computational domain */

    /* I. 1. Computational domain */

    create_vbox_group(&xgc.cs.cc.vbg_outer_comp_domain, &table_comp_domain, vbox_page, NULL, 2, 4);

    row_comp_domain = 0;

    //GtkWidget* ttt = tictactoe_new();
    //gtk_box_pack_start(GTK_BOX(vbox_page), ttt, FALSE, TRUE, 0);

    xgc.cs.cc.label_comp_domain_size = gtk_label_new("Total size (x, y, z) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.cc.label_comp_domain_size), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_comp_domain), xgc.cs.cc.label_comp_domain_size, 0, 2, row_comp_domain, row_comp_domain + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_comp_domain++;

    xgc.cs.cc.comp_domain_size_x = gtk_adjustment_new(COMP_DOMAIN_SIZE, 1, 10000, 1, 10, 0);
    xgc.cs.cc.comp_domain_size_x_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.cc.comp_domain_size_x), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.cc.comp_domain_size_x_spin), 5);
    gtk_table_attach(GTK_TABLE(table_comp_domain), xgc.cs.cc.comp_domain_size_x_spin, 0, 1, row_comp_domain, row_comp_domain + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.cc.comp_domain_size_x_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.cc.comp_domain_size_y = gtk_adjustment_new(COMP_DOMAIN_SIZE, 1, 10000, 1, 10, 0);
    xgc.cs.cc.comp_domain_size_y_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.cc.comp_domain_size_y), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.cc.comp_domain_size_y_spin), 5);
    gtk_table_attach(GTK_TABLE(table_comp_domain), xgc.cs.cc.comp_domain_size_y_spin, 1, 2, row_comp_domain, row_comp_domain + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.cc.comp_domain_size_y_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.cc.comp_domain_size_z = gtk_adjustment_new(COMP_DOMAIN_SIZE, 1, 10000, 1, 10, 0);
    xgc.cs.cc.comp_domain_size_z_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.cc.comp_domain_size_z), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.cc.comp_domain_size_z_spin), 5);
    gtk_table_attach(GTK_TABLE(table_comp_domain), xgc.cs.cc.comp_domain_size_z_spin, 2, 3, row_comp_domain, row_comp_domain + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.cc.comp_domain_size_z_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_comp_domain++;

    xgc.cs.cc.label_comp_domain_spacing = gtk_label_new("Voxel spacing (x, y, z) [μm]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.cc.label_comp_domain_spacing), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_comp_domain), xgc.cs.cc.label_comp_domain_spacing, 0, 2, row_comp_domain, row_comp_domain + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_comp_domain++;

    xgc.cs.cc.comp_domain_spacing_x = gtk_adjustment_new(COMP_DOMAIN_SPACING*1e6, 0, 1000, 1, 10, 0);    /* UI value in micrometers */
    xgc.cs.cc.comp_domain_spacing_x_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.cc.comp_domain_spacing_x), 1, 3);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.cc.comp_domain_spacing_x_spin), 5);
    gtk_table_attach(GTK_TABLE(table_comp_domain), xgc.cs.cc.comp_domain_spacing_x_spin, 0, 1, row_comp_domain, row_comp_domain + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.cc.comp_domain_spacing_x_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.cc.comp_domain_spacing_y = gtk_adjustment_new(COMP_DOMAIN_SPACING*1e6, 0, 1000, 1, 10, 0);    /* UI value in micrometers */
    xgc.cs.cc.comp_domain_spacing_y_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.cc.comp_domain_spacing_y), 1, 3);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.cc.comp_domain_spacing_y_spin), 5);
    gtk_table_attach(GTK_TABLE(table_comp_domain), xgc.cs.cc.comp_domain_spacing_y_spin, 1, 2, row_comp_domain, row_comp_domain + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.cc.comp_domain_spacing_y_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.cc.comp_domain_spacing_z = gtk_adjustment_new(COMP_DOMAIN_SPACING*1e6, 0, 1000, 1, 10, 0);    /* UI value in micrometers */
    xgc.cs.cc.comp_domain_spacing_z_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.cc.comp_domain_spacing_z), 1, 3);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.cc.comp_domain_spacing_z_spin), 5);
    gtk_table_attach(GTK_TABLE(table_comp_domain), xgc.cs.cc.comp_domain_spacing_z_spin, 2, 3, row_comp_domain, row_comp_domain + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.cc.comp_domain_spacing_z_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_comp_domain++;

    xgc.cs.cc.label_comp_domain_size_um = gtk_label_new("Total size (x, y, z) [μm]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.cc.label_comp_domain_size_um), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_comp_domain), xgc.cs.cc.label_comp_domain_size_um, 0, 2, row_comp_domain, row_comp_domain + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_comp_domain++;

    xgc.cs.cc.comp_domain_size_x_um = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.cc.comp_domain_size_x_um), 5);
    gtk_widget_set_sensitive(xgc.cs.cc.comp_domain_size_x_um, FALSE);
    gtk_table_attach(GTK_TABLE(table_comp_domain), xgc.cs.cc.comp_domain_size_x_um, 0, 1, row_comp_domain, row_comp_domain + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.cc.comp_domain_size_y_um = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.cc.comp_domain_size_y_um), 5);
    gtk_widget_set_sensitive(xgc.cs.cc.comp_domain_size_y_um, FALSE);
    gtk_table_attach(GTK_TABLE(table_comp_domain), xgc.cs.cc.comp_domain_size_y_um, 1, 2, row_comp_domain, row_comp_domain + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.cc.comp_domain_size_z_um = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.cc.comp_domain_size_z_um), 5);
    gtk_widget_set_sensitive(xgc.cs.cc.comp_domain_size_z_um, FALSE);
    gtk_table_attach(GTK_TABLE(table_comp_domain), xgc.cs.cc.comp_domain_size_z_um, 2, 3, row_comp_domain, row_comp_domain + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_comp_domain++;

    /* II. Basic parameters */

    /* II. 1. Basic parameters */

    create_vbox_group(&xgc.cs.bac.vbg_outer_basic, &table_basic, vbox_page, NULL, 2, 6);

    row_basic = 0;

    xgc.cs.bac.label_basic_nsteps = gtk_label_new("Number of steps:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.bac.label_basic_nsteps), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_basic), xgc.cs.bac.label_basic_nsteps, 0, 1, row_basic, row_basic + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs.bac.basic_nsteps = gtk_adjustment_new(BASIC_STEPS, 1, 1000000, 1, 10, 0);
    xgc.cs.bac.basic_nsteps_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.bac.basic_nsteps), 1, 0);
    gtk_table_attach(GTK_TABLE(table_basic), xgc.cs.bac.basic_nsteps_spin, 1, 2, row_basic, row_basic + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.bac.basic_nsteps_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_basic++;

    xgc.cs.bac.label_basic_verbose = gtk_label_new("Verbosity level:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.bac.label_basic_verbose), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_basic), xgc.cs.bac.label_basic_verbose, 0, 1, row_basic, row_basic + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs.bac.basic_verbose = gtk_adjustment_new(BASIC_VERBOSE, 0, 4, 1, 10, 0);
    xgc.cs.bac.basic_verbose_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.bac.basic_verbose), 1, 0);
    gtk_table_attach(GTK_TABLE(table_basic), xgc.cs.bac.basic_verbose_spin, 1, 2, row_basic, row_basic + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.bac.basic_verbose_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_basic++;

    xgc.cs.bac.label_basic_nthreads = gtk_label_new("Number of CPU cores (-1 for all):");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.bac.label_basic_nthreads), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_basic), xgc.cs.bac.label_basic_nthreads, 0, 2, row_basic, row_basic + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_basic++;

    xgc.cs.bac.basic_nthreads = gtk_adjustment_new(BASIC_NTHREADS, -1, 1000, 1, 10, 0);
    xgc.cs.bac.basic_nthreads_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.bac.basic_nthreads), 1, 0);
    gtk_table_attach(GTK_TABLE(table_basic), xgc.cs.bac.basic_nthreads_spin, 1, 2, row_basic, row_basic + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.bac.basic_nthreads_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_basic++;

    xgc.cs.bac.label_basic_ugpu = gtk_label_new("");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.bac.label_basic_ugpu), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_basic), xgc.cs.bac.label_basic_ugpu, 0, 1, row_basic, row_basic + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs.bac.basic_usegpu = gtk_check_button_new_with_label("Use GPU");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.bac.basic_usegpu), BASIC_USEGPU);
    gtk_table_attach(GTK_TABLE(table_basic), xgc.cs.bac.basic_usegpu, 1, 2, row_basic, row_basic + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.bac.basic_usegpu, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_basic++;

    xgc.cs.bac.label_basic_ugpu_index = gtk_label_new("GPU index:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.bac.label_basic_ugpu_index), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_basic), xgc.cs.bac.label_basic_ugpu_index, 0, 1, row_basic, row_basic + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs.bac.basic_ugpu_index = gtk_adjustment_new(BASIC_UGPUINDEX, 0, MAX_GPUS, 1, 10, 0);
    xgc.cs.bac.basic_ugpu_index_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.bac.basic_ugpu_index), 1, 0);
    gtk_table_attach(GTK_TABLE(table_basic), xgc.cs.bac.basic_ugpu_index_spin, 1, 2, row_basic, row_basic + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.bac.basic_ugpu_index_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    /* III. Sources */

    /* III. A. Source box group */

    create_vbox_group(&xgc.cs.sc.vbg_outer_box, &table_box, vbox_page, "Box:", 3, 9);
      
    row_box = 0;

    xgc.cs.sc.label_box_start = gtk_label_new("Start vertex (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_box_start), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.label_box_start, 0, 3, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_box++;

    xgc.cs.sc.box_i0 = gtk_adjustment_new(SOURCE_BOX_VERTEX0, 1, 10000, 1, 10, 0);
    xgc.cs.sc.box_i0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.box_i0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.box_i0_spin), 5);
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_i0_spin, 0, 1, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.box_i0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.sc.box_j0 = gtk_adjustment_new(SOURCE_BOX_VERTEX0, 1, 10000, 1, 10, 0);
    xgc.cs.sc.box_j0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.box_j0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.box_j0_spin), 5);
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_j0_spin, 1, 2, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.box_j0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.sc.box_k0 = gtk_adjustment_new(SOURCE_BOX_VERTEX0, 1, 10000, 1, 10, 0);
    xgc.cs.sc.box_k0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.box_k0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.box_k0_spin), 5);
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_k0_spin, 2, 3, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.box_k0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_box++;

    xgc.cs.sc.label_box_to = gtk_label_new("End vertex (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_box_to), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.label_box_to, 0, 3, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_box++;

    xgc.cs.sc.box_in = gtk_adjustment_new(SOURCE_BOX_VERTEXN, 1, 1000000, 1, 10, 0);
    xgc.cs.sc.box_in_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.box_in), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.box_in_spin), 5);
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_in_spin, 0, 1, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.box_in_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.sc.box_jn = gtk_adjustment_new(SOURCE_BOX_VERTEXN, 1, 1000000, 1, 10, 0);
    xgc.cs.sc.box_jn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.box_jn), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.box_jn_spin), 5);
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_jn_spin, 1, 2, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.box_jn_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.sc.box_kn = gtk_adjustment_new(SOURCE_BOX_VERTEXN, 1, 1000000, 1, 10, 0);
    xgc.cs.sc.box_kn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.box_kn), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.box_kn_spin), 5);
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_kn_spin, 2, 3, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.box_kn_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_box++;

    xgc.cs.sc.label_box_size_um = gtk_label_new("Total size (x, y, z) [μm]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_box_size_um), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.label_box_size_um, 0, 3, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_box++;

    xgc.cs.sc.box_size_x_um = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.box_size_x_um), 5);
    gtk_widget_set_sensitive(xgc.cs.sc.box_size_x_um, FALSE);
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_size_x_um, 0, 1, row_box, row_box + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.box_size_y_um = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.box_size_y_um), 5);
    gtk_widget_set_sensitive(xgc.cs.sc.box_size_y_um, FALSE);
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_size_y_um, 1, 2, row_box, row_box + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.box_size_z_um = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.box_size_z_um), 5);
    gtk_widget_set_sensitive(xgc.cs.sc.box_size_z_um, FALSE);
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_size_z_um, 2, 3, row_box, row_box + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_box++;

    xgc.cs.sc.label_box_boundary = gtk_label_new("Skip boundary:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_box_boundary), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.label_box_boundary, 0, 3, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_box++;

    xgc.cs.sc.box_boundary_skipi0 = gtk_check_button_new_with_label("i0");
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_boundary_skipi0, 0, 1, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.sc.box_boundary_skipi0), SOURCE_BOUNDARY_SKIP);
    g_signal_connect_swapped(xgc.cs.sc.box_boundary_skipi0, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.sc.box_boundary_skipj0 = gtk_check_button_new_with_label("j0");
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_boundary_skipj0, 1, 2, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.sc.box_boundary_skipj0), SOURCE_BOUNDARY_SKIP);
    g_signal_connect_swapped(xgc.cs.sc.box_boundary_skipj0, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.sc.box_boundary_skipk0 = gtk_check_button_new_with_label("k0");
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_boundary_skipk0, 2, 3, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.sc.box_boundary_skipk0), SOURCE_BOUNDARY_SKIP);
    g_signal_connect_swapped(xgc.cs.sc.box_boundary_skipk0, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_box++;

    xgc.cs.sc.box_boundary_skipin = gtk_check_button_new_with_label("in");
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_boundary_skipin, 0, 1, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.sc.box_boundary_skipin), SOURCE_BOUNDARY_SKIP);
    g_signal_connect_swapped(xgc.cs.sc.box_boundary_skipin, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.sc.box_boundary_skipjn = gtk_check_button_new_with_label("jn");
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_boundary_skipjn, 1, 2, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.sc.box_boundary_skipjn), SOURCE_BOUNDARY_SKIP);
    g_signal_connect_swapped(xgc.cs.sc.box_boundary_skipjn, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.sc.box_boundary_skipkn = gtk_check_button_new_with_label("kn");
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_boundary_skipkn, 2, 3, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.sc.box_boundary_skipkn), SOURCE_BOUNDARY_SKIP);
    g_signal_connect_swapped(xgc.cs.sc.box_boundary_skipkn, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_box++;    

    xgc.cs.sc.box_enable_skipdepth = gtk_check_button_new_with_label("Conformal");
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_enable_skipdepth, 0, 3, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.sc.box_enable_skipdepth), SOURCE_BOUNDARY_SKIP);
    g_signal_connect_swapped(xgc.cs.sc.box_enable_skipdepth, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_box++;

    xgc.cs.sc.label_box_skipdepth = gtk_label_new("Skip depth:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_box_skipdepth), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.label_box_skipdepth, 0, 2, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs.sc.box_skipdepth = gtk_adjustment_new(SOURCE_BOUNDARY_SKIP_DEPTH, -1, 1000000, 1, 10, 0);
    xgc.cs.sc.box_skipdepth_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.box_skipdepth), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.box_skipdepth_spin), 5);
    gtk_table_attach(GTK_TABLE(table_box), xgc.cs.sc.box_skipdepth_spin, 2, 3, row_box, row_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.box_skipdepth_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);


    /* III. 1. Point source properties */

    create_vbox_group(&xgc.cs.sc.vbg_outer_point_origin, &table_point_origin, vbox_page, "Point source properties:", 4, 2);

    row_point_origin = 0;

    xgc.cs.sc.label_point_origin_position = gtk_label_new("Source position (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_point_origin_position), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_point_origin), xgc.cs.sc.label_point_origin_position, 0, 3, row_point_origin, row_point_origin + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    hbox_point_origin_pos = gtk_hbox_new(FALSE, 0);

    row_point_origin++;

    xgc.cs.sc.point_origin_position_i = gtk_adjustment_new(SOURCE_POINT_ORIGIN_I, 1, 10000, 1, 10, 0);
    xgc.cs.sc.point_origin_position_i_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.point_origin_position_i), 1, 0);
    //gtk_table_attach(GTK_TABLE(table_point_origin), xgc.cs.sc.point_origin_position_i_spin, 0, 1, row_point_origin, row_point_origin + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_box_pack_start(GTK_BOX(hbox_point_origin_pos), xgc.cs.sc.point_origin_position_i_spin, TRUE, TRUE, 2);
    g_signal_connect_swapped(xgc.cs.sc.point_origin_position_i_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.sc.point_origin_position_j = gtk_adjustment_new(SOURCE_POINT_ORIGIN_J, 1, 10000, 1, 10, 0);
    xgc.cs.sc.point_origin_position_j_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.point_origin_position_j), 1, 0);
    //gtk_table_attach(GTK_TABLE(table_point_origin), xgc.cs.sc.point_origin_position_j_spin, 1, 2, row_point_origin, row_point_origin + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_box_pack_start(GTK_BOX(hbox_point_origin_pos), xgc.cs.sc.point_origin_position_j_spin, TRUE, TRUE, 2);
    g_signal_connect_swapped(xgc.cs.sc.point_origin_position_j_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.sc.point_origin_position_k = gtk_adjustment_new(SOURCE_POINT_ORIGIN_K, 1, 10000, 1, 10, 0);
    xgc.cs.sc.point_origin_position_k_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.point_origin_position_k), 1, 0);
    //gtk_table_attach(GTK_TABLE(table_point_origin), xgc.cs.sc.point_origin_position_k_spin, 2, 3, row_point_origin, row_point_origin + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_box_pack_start(GTK_BOX(hbox_point_origin_pos), xgc.cs.sc.point_origin_position_k_spin, TRUE, TRUE, 2);
    g_signal_connect_swapped(xgc.cs.sc.point_origin_position_k_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    gtk_table_attach(GTK_TABLE(table_point_origin), hbox_point_origin_pos, 0, 2, row_point_origin, row_point_origin + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_point_origin++;

    xgc.cs.sc.label_point_origin_theta = gtk_label_new("Theta [deg]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_point_origin_theta), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_point_origin), xgc.cs.sc.label_point_origin_theta, 0, 1, row_point_origin, row_point_origin + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs.sc.point_origin_theta = gtk_adjustment_new(SOURCE_POINT_ORIGIN_THETA, 0, 360, 1, 10, 0);
    xgc.cs.sc.point_origin_theta_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.point_origin_theta), 1, 3);
    gtk_table_attach(GTK_TABLE(table_point_origin), xgc.cs.sc.point_origin_theta_spin, 1, 2, row_point_origin, row_point_origin + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.point_origin_theta_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_point_origin++;

    xgc.cs.sc.label_point_origin_phi = gtk_label_new("Phi [deg]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_point_origin_phi), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_point_origin), xgc.cs.sc.label_point_origin_phi, 0, 1, row_point_origin, row_point_origin + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs.sc.point_origin_phi = gtk_adjustment_new(SOURCE_POINT_ORIGIN_PHI, 0, 360, 1, 10, 0);
    xgc.cs.sc.point_origin_phi_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.point_origin_phi), 1, 3);
    gtk_table_attach(GTK_TABLE(table_point_origin), xgc.cs.sc.point_origin_phi_spin, 1, 2, row_point_origin, row_point_origin + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.point_origin_phi_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_point_origin++;


    /* III. B. Source group */

    create_vbox_group(&xgc.cs.sc.vbg_outer_source, &table_source, vbox_page, "Source properties:", 4, 2);

    row_source = 0;

    xgc.cs.sc.label_source_mode = gtk_label_new("Source mode:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_source_mode), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_source), xgc.cs.sc.label_source_mode, 0, 1, row_source, row_source + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.source_mode = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.sc.source_mode), "Load from file");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.sc.source_mode), "Generate sine");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.sc.source_mode), "Generate pulse");
    gtk_combo_box_set_active(GTK_COMBO_BOX(xgc.cs.sc.source_mode), SOURCE_MODE);
    g_signal_connect_swapped(G_OBJECT(xgc.cs.sc.source_mode), "changed", G_CALLBACK(par_controls_changed_cb), &xgc);
    gtk_table_attach(GTK_TABLE(table_source), xgc.cs.sc.source_mode, 1, 2, row_source, row_source + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_source++;

    xgc.cs.sc.label_source_filename = gtk_label_new("Source file:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_source_filename), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_source), xgc.cs.sc.label_source_filename, 0, 1, row_source, row_source + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.source_filename = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.source_filename), ENTRY_FILENAME_CHARS);
    gtk_entry_set_text(GTK_ENTRY(xgc.cs.sc.source_filename), GENERAL_FILENAME);
    gtk_entry_set_editable(GTK_ENTRY(xgc.cs.sc.source_filename), FALSE);
    gtk_widget_set_can_focus(GTK_WIDGET(xgc.cs.sc.source_filename), TRUE);
    gtk_table_attach(GTK_TABLE(table_source), xgc.cs.sc.source_filename, 1, 2, row_source, row_source + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_source++;

    xgc.cs.sc.source_button_change = gtk_button_new_with_label("Change..");
    gtk_table_attach(GTK_TABLE(table_source), xgc.cs.sc.source_button_change, 1, 2, row_source, row_source + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect(xgc.cs.sc.source_button_change, "clicked", G_CALLBACK(ss_choose_source_file_cb), &xgc);

    row_source++;

    xgc.cs.sc.label_source_amplitude = gtk_label_new("Amplitude [V]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_source_amplitude), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_source), xgc.cs.sc.label_source_amplitude, 0, 1, row_source, row_source + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.source_amplitude = gtk_adjustment_new(SOURCE_AMPLITUDE, 0, 1000, 1, 10, 0);
    xgc.cs.sc.source_amplitude_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.source_amplitude), 1, 3);
    gtk_table_attach(GTK_TABLE(table_source), xgc.cs.sc.source_amplitude_spin, 1, 2, row_source, row_source + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.source_amplitude_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_source++;

    xgc.cs.sc.label_source_wavelength = gtk_label_new("Wavelength [μm]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_source_wavelength), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_source), xgc.cs.sc.label_source_wavelength, 0, 1, row_source, row_source + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.source_wavelength = gtk_adjustment_new(SOURCE_WAVELENGTH, -1, 1000, 1, 10, 0);
    xgc.cs.sc.source_wavelength_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.source_wavelength), 1, 3);
    gtk_table_attach(GTK_TABLE(table_source), xgc.cs.sc.source_wavelength_spin, 1, 2, row_source, row_source + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.source_wavelength_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_source++;

    xgc.cs.sc.label_source_pulsewidth = gtk_label_new("Pulse width [ts]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_source_pulsewidth), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_source), xgc.cs.sc.label_source_pulsewidth, 0, 1, row_source, row_source + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.source_pulsewidth = gtk_adjustment_new(SOURCE_PULSE_WIDTH, 1, 1000, 1, 10, 0);
    xgc.cs.sc.source_pulsewidth_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.source_pulsewidth), 1, 0);
    gtk_table_attach(GTK_TABLE(table_source), xgc.cs.sc.source_pulsewidth_spin, 1, 2, row_source, row_source + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.source_pulsewidth_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_source++;


    /* III. C. Incident angle group */

    create_vbox_group(&xgc.cs.sc.vbg_outer_ia, &table_ia, vbox_page, "Incident angle:", 3, 2);    

    row_ia = 0;

    xgc.cs.sc.label_ia_theta = gtk_label_new("Theta [deg]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_ia_theta), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_ia), xgc.cs.sc.label_ia_theta, 0, 1, row_ia, row_ia + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.ia_theta = gtk_adjustment_new(SOURCE_IA_THETADEG, 0, 360, 1, 10, 0);
    xgc.cs.sc.ia_theta_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.ia_theta), 1, 3);
    gtk_table_attach(GTK_TABLE(table_ia), xgc.cs.sc.ia_theta_spin, 1, 2, row_ia, row_ia + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.ia_theta_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_ia++;

    xgc.cs.sc.label_ia_phi = gtk_label_new("Phi [deg]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_ia_phi), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_ia), xgc.cs.sc.label_ia_phi, 0, 1, row_ia, row_ia + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.ia_phi = gtk_adjustment_new(SOURCE_IA_PHIDEG, 0, 360, 1, 10, 0);
    xgc.cs.sc.ia_phi_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.ia_phi), 1, 3);
    gtk_table_attach(GTK_TABLE(table_ia), xgc.cs.sc.ia_phi_spin, 1, 2, row_ia, row_ia + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.ia_phi_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_ia++;

    xgc.cs.sc.label_ia_psi = gtk_label_new("Psi [deg]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_ia_psi), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_ia), xgc.cs.sc.label_ia_psi, 0, 1, row_ia, row_ia + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.ia_psi = gtk_adjustment_new(SOURCE_IA_PSIDEG, 0, 360, 1, 10, 0);
    xgc.cs.sc.ia_psi_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.ia_psi), 1, 3);
    gtk_table_attach(GTK_TABLE(table_ia), xgc.cs.sc.ia_psi_spin, 1, 2, row_ia, row_ia + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.ia_psi_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_ia++;


    /* III. D. Focused source properties group */

    create_vbox_group(&xgc.cs.sc.vbg_outer_fs, &table_fs, vbox_page, "Focused source properties:", 5, 2);

    row_fs = 0;

    xgc.cs.sc.label_fs_thetamax = gtk_label_new("Theta max [deg]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_fs_thetamax), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_fs), xgc.cs.sc.label_fs_thetamax, 0, 1, row_fs, row_fs + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.fs_thetamax_deg = gtk_adjustment_new(SOURCE_FS_THETAMAXDEG, 0, 360, 1, 10, 0);
    xgc.cs.sc.fs_thetamax_deg_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.fs_thetamax_deg), 1, 3);
    gtk_table_attach(GTK_TABLE(table_fs), xgc.cs.sc.fs_thetamax_deg_spin, 1, 2, row_fs, row_fs + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.fs_thetamax_deg_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_fs++;

    xgc.cs.sc.label_fs_fdist = gtk_label_new("Focal distance [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_fs_fdist), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_fs), xgc.cs.sc.label_fs_fdist, 0, 1, row_fs, row_fs + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.fs_fdist = gtk_adjustment_new(SOURCE_FS_FDIST, 0, 360, 1, 10, 0);
    xgc.cs.sc.fs_fdist_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.fs_fdist), 1, 3);
    gtk_table_attach(GTK_TABLE(table_fs), xgc.cs.sc.fs_fdist_spin, 1, 2, row_fs, row_fs + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.fs_fdist_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_fs++;

    xgc.cs.sc.label_fs_polarisation = gtk_label_new("Polarisation [deg]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_fs_polarisation), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_fs), xgc.cs.sc.label_fs_polarisation, 0, 1, row_fs, row_fs + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.fs_pol_deg = gtk_adjustment_new(SOURCE_FS_POLARISATIONDEG, 0, 360, 1, 10, 0);
    xgc.cs.sc.fs_pol_deg_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.fs_pol_deg), 1, 3);
    gtk_table_attach(GTK_TABLE(table_fs), xgc.cs.sc.fs_pol_deg_spin, 1, 2, row_fs, row_fs + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.fs_pol_deg_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_fs++;

    xgc.cs.sc.label_fs_nip = gtk_label_new("N integration points:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_fs_nip), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_fs), xgc.cs.sc.label_fs_nip, 0, 1, row_fs, row_fs + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.fs_nip = gtk_adjustment_new(SOURCE_FS_NIP, 0, 360, 1, 10, 0);
    xgc.cs.sc.fs_nip_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.fs_nip), 1, 0);
    gtk_table_attach(GTK_TABLE(table_fs), xgc.cs.sc.fs_nip_spin, 1, 2, row_fs, row_fs + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.fs_nip_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_fs++;

    xgc.cs.sc.label_fs_mip = gtk_label_new("M integration points:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_fs_mip), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_fs), xgc.cs.sc.label_fs_mip, 0, 1, row_fs, row_fs + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.fs_mip = gtk_adjustment_new(SOURCE_FS_MIP, 0, 360, 1, 10, 0);
    xgc.cs.sc.fs_mip_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.fs_mip), 1, 0);
    gtk_table_attach(GTK_TABLE(table_fs), xgc.cs.sc.fs_mip_spin, 1, 2, row_fs, row_fs + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.fs_mip_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_fs++;

    /* III. E. Layered source group */

    xgc.cs.sc.vbg_outer_ls = vbg_outer_ls = gtk_vbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox_page), vbg_outer_ls, TRUE, TRUE, 4);

    xgc.cs.sc.hseparator_ls = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(vbg_outer_ls), xgc.cs.sc.hseparator_ls, TRUE, TRUE, 4);

    hbg_ls = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbg_outer_ls), hbg_ls, TRUE, TRUE, 4);

    vbg_inner_ls = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbg_inner_ls), 4);
    gtk_widget_set_size_request(vbg_inner_ls, 400, -1);
    gtk_box_pack_start(GTK_BOX(hbg_ls), vbg_inner_ls, FALSE, FALSE, 0);

    xgc.cs.sc.label_ls_prop = gtk_label_new("Layered source properties:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_ls_prop), 0.0, 0.5);
    gtk_box_pack_start(GTK_BOX(vbg_inner_ls), xgc.cs.sc.label_ls_prop, TRUE, TRUE, 4);

    treeview_packet_layered = NULL;
    //create_layer_treeview(&xgc.cs.sc.ls_prop, &treeview_packet_layered);
    create_layer_treeview((TreeViewWrapper**)&(xgc.cs.sc.ls_prop_wrapper), &treeview_packet_layered);
    gtk_widget_set_size_request(treeview_packet_layered, -1, 140);
    gtk_box_pack_start(GTK_BOX(vbg_inner_ls), treeview_packet_layered, TRUE, TRUE, 4);

    g_signal_new("cell-edited",
                 G_TYPE_OBJECT, G_SIGNAL_RUN_FIRST,
                 0, NULL, NULL,
                 g_cclosure_marshal_VOID__POINTER,
                 G_TYPE_NONE, 1, G_TYPE_POINTER);
    g_signal_connect_swapped(((TreeViewWrapper*)xgc.cs.sc.ls_prop_wrapper)->parent_treeview, "cell-edited", G_CALLBACK(par_controls_changed_cb), &xgc);

    /* III. F. Z multiplier group */

    create_vbox_group(&xgc.cs.sc.vbg_outer_zmultiplier, &table_zmultiplier, vbox_page, "Z multiplier:", 5, 3);

    row_zmultiplier = 0;    

    xgc.cs.sc.gaussian_mult_enable = gtk_check_button_new_with_label("Enable Gaussian multiplier");
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.gaussian_mult_enable, 0, 3, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.sc.gaussian_mult_enable), SOURCE_GAUSSIAN_MULT_ENABLE);
    g_signal_connect_swapped(xgc.cs.sc.gaussian_mult_enable, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_zmultiplier++;

    xgc.cs.sc.label_gaussian_mult_center = gtk_label_new("Center (i, j) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_gaussian_mult_center), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.label_gaussian_mult_center, 0, 1, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.gaussian_mult_center_pos_i = gtk_adjustment_new(COMP_DOMAIN_SIZE/2, 1, 10000, 1, 10, 0);
    xgc.cs.sc.gaussian_mult_center_pos_i_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.gaussian_mult_center_pos_i), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.gaussian_mult_center_pos_i_spin), 7);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.gaussian_mult_center_pos_i_spin, 1, 2, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.gaussian_mult_center_pos_i_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.sc.gaussian_mult_center_pos_j = gtk_adjustment_new(COMP_DOMAIN_SIZE / 2, 1, 10000, 1, 10, 0);
    xgc.cs.sc.gaussian_mult_center_pos_j_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.gaussian_mult_center_pos_j), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.gaussian_mult_center_pos_j_spin), 7);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.gaussian_mult_center_pos_j_spin, 2, 3, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.gaussian_mult_center_pos_j_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_zmultiplier++;

    xgc.cs.sc.label_gaussian_mult_radius = gtk_label_new("Radius (i, j) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_gaussian_mult_radius), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.label_gaussian_mult_radius, 0, 1, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.gaussian_mult_radius_i = gtk_adjustment_new(SOURCE_XXX_MULT_RADIUS, 1, 10000, 1, 10, 0);
    xgc.cs.sc.gaussian_mult_radius_i_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.gaussian_mult_radius_i), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.gaussian_mult_radius_i_spin), 7);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.gaussian_mult_radius_i_spin, 1, 2, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.gaussian_mult_radius_i_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.sc.gaussian_mult_radius_j = gtk_adjustment_new(SOURCE_XXX_MULT_RADIUS, 1, 10000, 1, 10, 0);
    xgc.cs.sc.gaussian_mult_radius_j_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.gaussian_mult_radius_j), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.gaussian_mult_radius_j_spin), 7);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.gaussian_mult_radius_j_spin, 2, 3, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.gaussian_mult_radius_j_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_zmultiplier++;

    xgc.cs.sc.radial_mult_enable = gtk_check_button_new_with_label("Enable radial multiplier");
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.radial_mult_enable, 0, 3, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.sc.radial_mult_enable), SOURCE_RADIAL_MULT_ENABLE);
    g_signal_connect_swapped(xgc.cs.sc.radial_mult_enable, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_zmultiplier++;

    xgc.cs.sc.label_radial_mult_center = gtk_label_new("Center (i, j) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_radial_mult_center), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.label_radial_mult_center, 0, 1, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.radial_mult_center_pos_i = gtk_adjustment_new(COMP_DOMAIN_SIZE / 2, 1, 10000, 1, 10, 0);
    xgc.cs.sc.radial_mult_center_pos_i_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.radial_mult_center_pos_i), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.radial_mult_center_pos_i_spin), 7);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.radial_mult_center_pos_i_spin, 1, 2, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.radial_mult_center_pos_i_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.sc.radial_mult_center_pos_j = gtk_adjustment_new(COMP_DOMAIN_SIZE / 2, 1, 10000, 1, 10, 0);
    xgc.cs.sc.radial_mult_center_pos_j_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.radial_mult_center_pos_j), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.radial_mult_center_pos_j_spin), 7);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.radial_mult_center_pos_j_spin, 2, 3, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.radial_mult_center_pos_j_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_zmultiplier++;

    xgc.cs.sc.label_radial_mult_radius = gtk_label_new("Radius (i, j) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_radial_mult_radius), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.label_radial_mult_radius, 0, 1, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.radial_mult_radius_i = gtk_adjustment_new(SOURCE_XXX_MULT_RADIUS, 1, 10000, 1, 10, 0);
    xgc.cs.sc.radial_mult_radius_i_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.radial_mult_radius_i), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.radial_mult_radius_i_spin), 7);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.radial_mult_radius_i_spin, 1, 2, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.radial_mult_radius_i_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.sc.radial_mult_radius_j = gtk_adjustment_new(SOURCE_XXX_MULT_RADIUS, 1, 10000, 1, 10, 0);
    xgc.cs.sc.radial_mult_radius_j_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.radial_mult_radius_j), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.radial_mult_radius_j_spin), 7);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.radial_mult_radius_j_spin, 2, 3, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.radial_mult_radius_j_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_zmultiplier++;

    xgc.cs.sc.fiber_mult_enable = gtk_check_button_new_with_label("Enable fiber multiplier");
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.fiber_mult_enable, 0, 3, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.sc.fiber_mult_enable), SOURCE_FIBER_MULT_ENABLE);
    g_signal_connect_swapped(xgc.cs.sc.fiber_mult_enable, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_zmultiplier++;

    xgc.cs.sc.label_fiber_mult_center = gtk_label_new("Center (i, j) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_fiber_mult_center), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.label_fiber_mult_center, 0, 1, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.fiber_mult_center_pos_i = gtk_adjustment_new(COMP_DOMAIN_SIZE / 2, 1, 10000, 1, 10, 0);
    xgc.cs.sc.fiber_mult_center_pos_i_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.fiber_mult_center_pos_i), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.fiber_mult_center_pos_i_spin), 7);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.fiber_mult_center_pos_i_spin, 1, 2, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.fiber_mult_center_pos_i_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.sc.fiber_mult_center_pos_j = gtk_adjustment_new(COMP_DOMAIN_SIZE / 2, 1, 10000, 1, 10, 0);
    xgc.cs.sc.fiber_mult_center_pos_j_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.fiber_mult_center_pos_j), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.fiber_mult_center_pos_j_spin), 7);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.fiber_mult_center_pos_j_spin, 2, 3, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.fiber_mult_center_pos_j_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_zmultiplier++;

    xgc.cs.sc.label_fiber_mult_radius = gtk_label_new("Radius [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_fiber_mult_radius), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.label_fiber_mult_radius, 0, 1, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.fiber_mult_radius = gtk_adjustment_new(SOURCE_XXX_MULT_RADIUS, 1, 10000, 1, 10, 0);
    xgc.cs.sc.fiber_mult_radius_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.fiber_mult_radius), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.fiber_mult_radius_spin), 7);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.fiber_mult_radius_spin, 1, 2, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.fiber_mult_radius_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_zmultiplier++;

    xgc.cs.sc.label_fiber_mult_cutoff = gtk_label_new("Cutoff [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_fiber_mult_cutoff), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.label_fiber_mult_cutoff, 0, 1, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.fiber_mult_cutoff = gtk_adjustment_new(SOURCE_FIBER_MULT_CUTOFF, 1, 10000, 1, 10, 0);
    xgc.cs.sc.fiber_mult_cutoff_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.fiber_mult_cutoff), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.fiber_mult_cutoff_spin), 7);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.fiber_mult_cutoff_spin, 1, 2, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.fiber_mult_cutoff_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_zmultiplier++;

    xgc.cs.sc.label_fiber_mult_epsilon_core = gtk_label_new("ε core:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_fiber_mult_epsilon_core), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.label_fiber_mult_epsilon_core, 0, 1, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.fiber_mult_epsilon_core = gtk_adjustment_new(SOURCE_FIBER_MULT_EPSILON_CORE, 1, 20, 1, 10, 0);
    xgc.cs.sc.fiber_mult_epsilon_core_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.fiber_mult_epsilon_core), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.fiber_mult_epsilon_core_spin), 4);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.fiber_mult_epsilon_core_spin, 1, 2, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.fiber_mult_epsilon_core_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_zmultiplier++;

    xgc.cs.sc.label_fiber_mult_epsilon_cladding = gtk_label_new("ε cladding:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.sc.label_fiber_mult_epsilon_cladding), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.label_fiber_mult_epsilon_cladding, 0, 1, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.sc.fiber_mult_epsilon_cladding = gtk_adjustment_new(SOURCE_FIBER_MULT_EPSILON_CLADDING, 1, 20, 1, 10, 0);
    xgc.cs.sc.fiber_mult_epsilon_cladding_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.sc.fiber_mult_epsilon_cladding), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.sc.fiber_mult_epsilon_cladding_spin), 4);
    gtk_table_attach(GTK_TABLE(table_zmultiplier), xgc.cs.sc.fiber_mult_epsilon_cladding_spin, 1, 2, row_zmultiplier, row_zmultiplier + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.sc.fiber_mult_epsilon_cladding_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    /* IV. Boundary conditions */

    /* IV. 1. Computational volume boundaries */

    create_vbox_group(&xgc.cs.boc.vbg_outer_cvboundary, &table_cvboundary, vbox_page, "Computational volume boundaries:", 19, 5);

    gchar *bnd[6] = { "x0 boundary:", "xn boundary:", 
                      "y0 boundary:", "yn boundary:", 
                      "z0 boundary:", "zn boundary:" };

    row_cvboundary = 0;

    for (i = 0; i < 6; i++) {
        xgc.cs.boc.label_btype[i] = gtk_label_new(bnd[i]);
        gtk_misc_set_alignment(GTK_MISC(xgc.cs.boc.label_btype[i]), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(table_cvboundary), xgc.cs.boc.label_btype[i], 0, 1, row_cvboundary, row_cvboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        xgc.cs.boc.btype[i] = gtk_combo_box_new_text();
        gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.boc.btype[i]), "None");
        gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.boc.btype[i]), "PEC");
        gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.boc.btype[i]), "Liao");
        gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.boc.btype[i]), "CPML");
        gtk_combo_box_set_active(GTK_COMBO_BOX(xgc.cs.boc.btype[i]), cvboundary_typeval[i]);

        g_signal_connect_swapped(G_OBJECT(xgc.cs.boc.btype[i]), "changed", G_CALLBACK(par_controls_changed_cb), &xgc);

        gtk_table_attach(GTK_TABLE(table_cvboundary), xgc.cs.boc.btype[i], 1, 2, row_cvboundary, row_cvboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        row_cvboundary++;

        xgc.cs.boc.label_cvCPMLparams = gtk_label_new("CPML parameters (depth, m, σ, a, κ):");
        gtk_misc_set_alignment(GTK_MISC(xgc.cs.boc.label_cvCPMLparams), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(table_cvboundary), xgc.cs.boc.label_cvCPMLparams, 0, 3, row_cvboundary, row_cvboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

        row_cvboundary++;

        xgc.cs.boc.bdepth[i] = gtk_adjustment_new(cvboundary_depthval[i], 1, 1000, 1, 10, 0);
        xgc.cs.boc.bdepthspin[i] = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.boc.bdepth[i]), 1, 0);
        gtk_table_attach(GTK_TABLE(table_cvboundary), xgc.cs.boc.bdepthspin[i], 0, 1, row_cvboundary, row_cvboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(xgc.cs.boc.bdepthspin[i], "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

        xgc.cs.boc.bm[i] = gtk_adjustment_new(cvboundary_mval[i], 1, 100, 1, 10, 0);
        xgc.cs.boc.bmspin[i] = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.boc.bm[i]), 1, 0);
        gtk_table_attach(GTK_TABLE(table_cvboundary), xgc.cs.boc.bmspin[i], 1, 2, row_cvboundary, row_cvboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(xgc.cs.boc.bmspin[i], "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

        xgc.cs.boc.bsigma[i] = gtk_adjustment_new(cvboundary_sigmaval[i], -1, 1000000, 1, 10, 0);
        xgc.cs.boc.bsigmaspin[i] = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.boc.bsigma[i]), 1, 3);
        gtk_table_attach(GTK_TABLE(table_cvboundary), xgc.cs.boc.bsigmaspin[i], 2, 3, row_cvboundary, row_cvboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(xgc.cs.boc.bsigmaspin[i], "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

        xgc.cs.boc.ba[i] = gtk_adjustment_new(cvboundary_aval[i], 0, 100, 1, 10, 0);
        xgc.cs.boc.baspin[i] = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.boc.ba[i]), 1, 3);
        gtk_table_attach(GTK_TABLE(table_cvboundary), xgc.cs.boc.baspin[i], 3, 4, row_cvboundary, row_cvboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(xgc.cs.boc.baspin[i], "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

        xgc.cs.boc.bkappa[i] = gtk_adjustment_new(cvboundary_kappaval[i], 1, 100, 1, 10, 0);
        xgc.cs.boc.bkappaspin[i] = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.boc.bkappa[i]), 1, 3);
        gtk_table_attach(GTK_TABLE(table_cvboundary), xgc.cs.boc.bkappaspin[i], 4, 5, row_cvboundary, row_cvboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
        g_signal_connect_swapped(xgc.cs.boc.bkappaspin[i], "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

        row_cvboundary++;
    }

    /* IV. 2. Periodic boundaries inside computational volume */

    create_vbox_group(&xgc.cs.boc.vbg_outer_periodicboundary, &table_periodicboundary, vbox_page, "Computational volume boundaries:", 3, 4);

    row_periodicboundary = 0;

    xgc.cs.boc.mb[0] = gtk_check_button_new_with_mnemonic("X0 periodic at ");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.boc.mb[0]), pboundary_mbval[0]);
    //g_signal_connect(G_OBJECT(xgc.cs.boc.mbx[0]), "clicked", G_CALLBACK(mbx_cb), &xgc);
    g_signal_connect_swapped(G_OBJECT(xgc.cs.boc.mb[0]), "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);
    gtk_table_attach(GTK_TABLE(table_periodicboundary), xgc.cs.boc.mb[0], 0, 1, row_periodicboundary, row_periodicboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.boc.mbpos[0] = gtk_adjustment_new(pboundary_mbposval[0], 1, 10000, 1, 10, 0);
    xgc.cs.boc.mbposspin[0] = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.boc.mbpos[0]), 1, 0);
    gtk_widget_set_sensitive(GTK_WIDGET(xgc.cs.boc.mbposspin[0]), pboundary_mbval[0]);
    gtk_table_attach(GTK_TABLE(table_periodicboundary), xgc.cs.boc.mbposspin[0], 1, 2, row_periodicboundary, row_periodicboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.boc.mbposspin[0], "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.boc.mb[1] = gtk_check_button_new_with_mnemonic("XN periodic at ");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.boc.mb[1]), pboundary_mbval[1]);
    //g_signal_connect(G_OBJECT(xgc.cs.boc.mb[1]), "clicked", G_CALLBACK(mb_cb), &xgc);
    g_signal_connect_swapped(G_OBJECT(xgc.cs.boc.mb[1]), "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);
    gtk_table_attach(GTK_TABLE(table_periodicboundary), xgc.cs.boc.mb[1], 2, 3, row_periodicboundary, row_periodicboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.boc.mbpos[1] = gtk_adjustment_new(pboundary_mbposval[1], 1, 10000, 1, 10, 0);
    xgc.cs.boc.mbposspin[1] = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.boc.mbpos[1]), 1, 0);
    gtk_widget_set_sensitive(GTK_WIDGET(xgc.cs.boc.mbposspin[1]), pboundary_mbval[1]);
    gtk_table_attach(GTK_TABLE(table_periodicboundary), xgc.cs.boc.mbposspin[1], 3, 4, row_periodicboundary, row_periodicboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.boc.mbposspin[1], "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_periodicboundary++;

    xgc.cs.boc.mb[2] = gtk_check_button_new_with_mnemonic("Y0 periodic at ");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.boc.mb[2]), pboundary_mbval[2]);
    //g_signal_connect(G_OBJECT(xgc.cs.boc.mb[2]), "clicked", G_CALLBACK(mb_cb), &xgc);
    g_signal_connect_swapped(G_OBJECT(xgc.cs.boc.mb[2]), "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);
    gtk_table_attach(GTK_TABLE(table_periodicboundary), xgc.cs.boc.mb[2], 0, 1, row_periodicboundary, row_periodicboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.boc.mbpos[2] = gtk_adjustment_new(pboundary_mbposval[2], 1, 10000, 1, 10, 0);
    xgc.cs.boc.mbposspin[2] = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.boc.mbpos[2]), 1, 0);
    gtk_widget_set_sensitive(GTK_WIDGET(xgc.cs.boc.mbposspin[2]), pboundary_mbval[2]);
    gtk_table_attach(GTK_TABLE(table_periodicboundary), xgc.cs.boc.mbposspin[2], 1, 2, row_periodicboundary, row_periodicboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.boc.mbposspin[2], "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.boc.mb[3] = gtk_check_button_new_with_mnemonic("YN periodic at ");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.boc.mb[3]), pboundary_mbval[3]);
    //g_signal_connect(G_OBJECT(xgc.cs.boc.mb[3]), "clicked", G_CALLBACK(mb_cb), &xgc);
    g_signal_connect_swapped(G_OBJECT(xgc.cs.boc.mb[3]), "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);
    gtk_table_attach(GTK_TABLE(table_periodicboundary), xgc.cs.boc.mb[3], 2, 3, row_periodicboundary, row_periodicboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.boc.mbpos[3] = gtk_adjustment_new(pboundary_mbposval[3], 1, 10000, 1, 10, 0);
    xgc.cs.boc.mbposspin[3] = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.boc.mbpos[3]), 1, 0);
    gtk_widget_set_sensitive(GTK_WIDGET(xgc.cs.boc.mbposspin[3]), pboundary_mbval[3]);
    gtk_table_attach(GTK_TABLE(table_periodicboundary), xgc.cs.boc.mbposspin[3], 3, 4, row_periodicboundary, row_periodicboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.boc.mbposspin[3], "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_periodicboundary++;

    xgc.cs.boc.mb[4] = gtk_check_button_new_with_mnemonic("Z0 periodic at ");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.boc.mb[4]), pboundary_mbval[4]);
    //g_signal_connect(G_OBJECT(xgc.cs.boc.mb[4]), "clicked", G_CALLBACK(mb_cb), &xgc);
    g_signal_connect_swapped(G_OBJECT(xgc.cs.boc.mb[4]), "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);
    gtk_table_attach(GTK_TABLE(table_periodicboundary), xgc.cs.boc.mb[4], 0, 1, row_periodicboundary, row_periodicboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.boc.mbpos[4] = gtk_adjustment_new(pboundary_mbposval[4], 1, 10000, 1, 10, 0);
    xgc.cs.boc.mbposspin[4] = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.boc.mbpos[4]), 1, 0);
    gtk_widget_set_sensitive(GTK_WIDGET(xgc.cs.boc.mbposspin[4]), pboundary_mbval[4]);
    gtk_table_attach(GTK_TABLE(table_periodicboundary), xgc.cs.boc.mbposspin[4], 1, 2, row_periodicboundary, row_periodicboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.boc.mbposspin[4], "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.boc.mb[5] = gtk_check_button_new_with_mnemonic("ZN periodic at ");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.boc.mb[5]), pboundary_mbval[5]);
    //g_signal_connect(G_OBJECT(xgc.cs.boc.mb[5]), "clicked", G_CALLBACK(mb_cb), &xgc);
    g_signal_connect_swapped(G_OBJECT(xgc.cs.boc.mb[5]), "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);
    gtk_table_attach(GTK_TABLE(table_periodicboundary), xgc.cs.boc.mb[5], 2, 3, row_periodicboundary, row_periodicboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.boc.mbpos[5] = gtk_adjustment_new(pboundary_mbposval[5], 1, 10000, 1, 10, 0);
    xgc.cs.boc.mbposspin[5] = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.boc.mbpos[5]), 1, 0);
    gtk_widget_set_sensitive(GTK_WIDGET(xgc.cs.boc.mbposspin[5]), pboundary_mbval[5]);
    gtk_table_attach(GTK_TABLE(table_periodicboundary), xgc.cs.boc.mbposspin[5], 3, 4, row_periodicboundary, row_periodicboundary + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.boc.mbposspin[5], "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    /* V. Media */

    /* V. 1. Material properties */

    create_vbox_group(&xgc.cs.mc.vbg_outer_material_prop, &table_material_prop, vbox_page, NULL, 8, 1);

    row_material_prop = 0;

    xgc.cs.mc.label_material_voxel_filename = gtk_label_new("Voxel by voxel material file:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_material_voxel_filename), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_prop), xgc.cs.mc.label_material_voxel_filename, 0, 2, row_material_prop, row_material_prop + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_material_prop++;

    xgc.cs.mc.material_voxel_filename = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.mc.material_voxel_filename), ENTRY_FILENAME_CHARS);
    gtk_entry_set_text(GTK_ENTRY(xgc.cs.mc.material_voxel_filename), GENERAL_FILENAME);
    gtk_entry_set_editable(GTK_ENTRY(xgc.cs.mc.material_voxel_filename), FALSE);
    gtk_widget_set_can_focus(GTK_WIDGET(xgc.cs.mc.material_voxel_filename), TRUE);
    gtk_table_attach(GTK_TABLE(table_material_prop), xgc.cs.mc.material_voxel_filename, 0, 2, row_material_prop, row_material_prop + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_material_prop++;

    xgc.cs.mc.material_voxel_button_change = gtk_button_new_with_label("Change..");
    gtk_table_attach(GTK_TABLE(table_material_prop), xgc.cs.mc.material_voxel_button_change, 0, 2, row_material_prop, row_material_prop + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect(xgc.cs.mc.material_voxel_button_change, "clicked", G_CALLBACK(material_choose_voxel_file_cb), &xgc);

    row_material_prop++;

    xgc.cs.mc.label_material_vector_filename = gtk_label_new("Vector material file:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_material_vector_filename), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_prop), xgc.cs.mc.label_material_vector_filename, 0, 2, row_material_prop, row_material_prop + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_material_prop++;

    xgc.cs.mc.material_vector_filename = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.mc.material_vector_filename), ENTRY_FILENAME_CHARS);
    //gtk_entry_set_text(GTK_ENTRY(xgc.cs.mc.material_vector_filename), GENERAL_FILENAME);
    gtk_entry_set_editable(GTK_ENTRY(xgc.cs.mc.material_vector_filename), FALSE);
    gtk_widget_set_can_focus(GTK_WIDGET(xgc.cs.mc.material_vector_filename), TRUE);
    gtk_table_attach(GTK_TABLE(table_material_prop), xgc.cs.mc.material_vector_filename, 0, 2, row_material_prop, row_material_prop + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_material_prop++;

    xgc.cs.mc.material_vector_button_change = gtk_button_new_with_label("Change..");
    gtk_table_attach(GTK_TABLE(table_material_prop), xgc.cs.mc.material_vector_button_change, 0, 2, row_material_prop, row_material_prop + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect(xgc.cs.mc.material_vector_button_change, "clicked", G_CALLBACK(material_choose_vector_file_cb), &xgc);

    row_material_prop++;

    xgc.cs.mc.label_material_modecheck = gtk_label_new("");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_material_modecheck), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_prop), xgc.cs.mc.label_material_modecheck, 0, 2, row_material_prop, row_material_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs.mc.material_modecheck = gtk_check_button_new_with_label("Material mode check");
    gtk_table_attach(GTK_TABLE(table_material_prop), xgc.cs.mc.material_modecheck, 0, 1, row_material_prop, row_material_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.mc.material_modecheck), MATERIAL_MODECHECK);
    g_signal_connect_swapped(xgc.cs.mc.material_modecheck, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_prop++;

    xgc.cs.mc.label_material_smoothsteps = gtk_label_new("Material smooth steps:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_material_smoothsteps), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_prop), xgc.cs.mc.label_material_smoothsteps, 0, 1, row_material_prop, row_material_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs.mc.material_smoothsteps = gtk_adjustment_new(MATERIAL_SMOOTHSTEPS, 0, 360, 1, 10, 0);
    xgc.cs.mc.material_smoothsteps_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.material_smoothsteps), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_prop), xgc.cs.mc.material_smoothsteps_spin, 1, 2, row_material_prop, row_material_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.material_smoothsteps_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_prop++;

    /* V. 2. Add growth modifier */

    create_vbox_group(&xgc.cs.mc.vbg_outer_material_grow, &table_material_grow, vbox_page, NULL, 13, 3);

    row_material_grow = 0;

    xgc.cs.mc.label_grow_span_from = gtk_label_new("Generator span from (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_grow_span_from), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.label_grow_span_from, 0, 3, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_material_grow++;

    xgc.cs.mc.grow_i0 = gtk_adjustment_new(MATERIAL_GROW_I0, 1, 10000, 1, 10, 0);
    xgc.cs.mc.grow_i0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.grow_i0), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_i0_spin, 0, 1, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_i0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.mc.grow_j0 = gtk_adjustment_new(MATERIAL_GROW_J0, 1, 10000, 1, 10, 0);
    xgc.cs.mc.grow_j0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.grow_j0), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_j0_spin, 1, 2, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_j0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.mc.grow_k0 = gtk_adjustment_new(MATERIAL_GROW_K0, 1, 10000, 1, 10, 0);
    xgc.cs.mc.grow_k0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.grow_k0), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_k0_spin, 2, 3, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_k0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_grow++;

    xgc.cs.mc.label_grow_span_to = gtk_label_new("Generator span to (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_grow_span_to), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.label_grow_span_to, 0, 3, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_material_grow++;

    xgc.cs.mc.grow_in = gtk_adjustment_new(material_grow_inval, 1, 10000, 1, 10, 0);
    xgc.cs.mc.grow_in_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.grow_in), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_in_spin, 0, 1, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_in_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.mc.grow_jn = gtk_adjustment_new(material_grow_jnval, 1, 10000, 1, 10, 0);
    xgc.cs.mc.grow_jn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.grow_jn), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_jn_spin, 1, 2, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_jn_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.mc.grow_kn = gtk_adjustment_new(material_grow_knval, 1, 10000, 1, 10, 0);
    xgc.cs.mc.grow_kn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.grow_kn), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_kn_spin, 2, 3, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_kn_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_grow++;

    xgc.cs.mc.grow_skipi0 = gtk_check_button_new_with_mnemonic("skip i0 face");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.mc.grow_skipi0), MATERIAL_GROW_SKIP);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_skipi0, 0, 1, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_skipi0, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.mc.grow_skipj0 = gtk_check_button_new_with_mnemonic("skip j0 face");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.mc.grow_skipj0), MATERIAL_GROW_SKIP);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_skipj0, 1, 2, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_skipj0, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.mc.grow_skipk0 = gtk_check_button_new_with_mnemonic("skip k0 face");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.mc.grow_skipk0), MATERIAL_GROW_SKIP);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_skipk0, 2, 3, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_skipk0, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_grow++;

    xgc.cs.mc.grow_skipin = gtk_check_button_new_with_mnemonic("skip in face");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.mc.grow_skipin), MATERIAL_GROW_SKIP);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_skipin, 0, 1, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_skipin, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.mc.grow_skipjn = gtk_check_button_new_with_mnemonic("skip jn face");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.mc.grow_skipjn), MATERIAL_GROW_SKIP);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_skipjn, 1, 2, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_skipjn, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.mc.grow_skipkn = gtk_check_button_new_with_mnemonic("skip kn face");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.mc.grow_skipkn), MATERIAL_GROW_SKIP);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_skipkn, 2, 3, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_skipkn, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_grow++;

    xgc.cs.mc.label_grow_attachindex = gtk_label_new("Material index where to attach:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_grow_attachindex), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.label_grow_attachindex, 0, 2, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.grow_attachindex = gtk_adjustment_new(MATERIAL_GROW_ATTACHINDEX, 0, 100, 1, 10, 0);
    xgc.cs.mc.grow_attachindex_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.grow_attachindex), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_attachindex_spin, 2, 3, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_attachindex_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_grow++;

    xgc.cs.mc.label_grow_addindex = gtk_label_new("Material index to grow:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_grow_addindex), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.label_grow_addindex, 0, 2, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.grow_addindex = gtk_adjustment_new(MATERIAL_GROW_ADDINDEX, 0, 100, 1, 10, 0);
    xgc.cs.mc.grow_addindex_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.grow_addindex), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_addindex_spin, 2, 3, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_addindex_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_grow++;

    xgc.cs.mc.label_grow_subsampling = gtk_label_new("Volume subsampling:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_grow_subsampling), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.label_grow_subsampling, 0, 2, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.grow_subsampling = gtk_adjustment_new(MATERIAL_GROW_SUBSAMPLING, 0, 100, 1, 10, 0);
    xgc.cs.mc.grow_subsampling_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.grow_subsampling), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_subsampling_spin, 2, 3, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_subsampling_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_grow++;

    xgc.cs.mc.label_grow_seed = gtk_label_new("Generator seed (-1 for random):");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_grow_seed), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.label_grow_seed, 0, 2, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.grow_seed = gtk_adjustment_new(MATERIAL_GROW_SEED, -1, 1000, 1, 10, 0);
    xgc.cs.mc.grow_seed_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.grow_seed), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_seed_spin, 2, 3, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_seed_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_grow++;

    xgc.cs.mc.label_grow_nsteps = gtk_label_new("Number of particles:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_grow_nsteps), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.label_grow_nsteps, 0, 2, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.grow_nsteps = gtk_adjustment_new(MATERIAL_GROW_NSTEPS, 0, 10000000, 1, 10, 0);
    xgc.cs.mc.grow_nsteps_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.grow_nsteps), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_nsteps_spin, 2, 3, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_nsteps_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_grow++;

    xgc.cs.mc.label_grow_mobility = gtk_label_new("Particle mobility:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_grow_mobility), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.label_grow_mobility, 0, 2, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.grow_mobility = gtk_adjustment_new(MATERIAL_GROW_MOBILITY, 0, 100, 1, 10, 0);
    xgc.cs.mc.grow_mobility_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.grow_mobility), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_mobility_spin, 2, 3, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_mobility_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_grow++;

    xgc.cs.mc.label_grow_probability = gtk_label_new("Relaxation probability:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_grow_probability), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.label_grow_probability, 0, 2, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.grow_probability = gtk_adjustment_new(MATERIAL_GROW_PROBABILITY, 0, 1.0, 0.01, 0.1, 0);
    xgc.cs.mc.grow_probability_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.grow_probability), 1, 3);
    gtk_table_attach(GTK_TABLE(table_material_grow), xgc.cs.mc.grow_probability_spin, 2, 3, row_material_grow, row_material_grow + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.grow_probability_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_grow++;

    /* V. 3. Add roughness modifier */

    create_vbox_group(&xgc.cs.mc.vbg_outer_material_rough, &table_material_rough, vbox_page, NULL, 7, 2);

    row_material_rough = 0;

    xgc.cs.mc.label_rough_matindex = gtk_label_new("Material index to add:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_rough_matindex), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_rough), xgc.cs.mc.label_rough_matindex, 0, 2, row_material_rough, row_material_rough + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.rough_matindex = gtk_adjustment_new(MATERIAL_ROUGHNESS_MATINDEX, 0, 100, 1, 10, 0);
    xgc.cs.mc.rough_matindex_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.rough_matindex), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_rough), xgc.cs.mc.rough_matindex_spin, 2, 3, row_material_rough, row_material_rough + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.rough_matindex_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_rough++;

    xgc.cs.mc.label_rough_voidindex = gtk_label_new("Void index:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_rough_voidindex), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_rough), xgc.cs.mc.label_rough_voidindex, 0, 2, row_material_rough, row_material_rough + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.rough_voidindex = gtk_adjustment_new(MATERIAL_ROUGHNESS_VOIDINDEX, 0, 100, 1, 10, 0);
    xgc.cs.mc.rough_voidindex_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.rough_voidindex), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_rough), xgc.cs.mc.rough_voidindex_spin, 2, 3, row_material_rough, row_material_rough + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.rough_voidindex_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_rough++;

    xgc.cs.mc.label_rough_iterations = gtk_label_new("Number of iterations:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_rough_iterations), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_rough), xgc.cs.mc.label_rough_iterations, 0, 2, row_material_rough, row_material_rough + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.rough_iterations = gtk_adjustment_new(MATERIAL_ROUGHNESS_ITERATIONS, 0, 100, 1, 10, 0);
    xgc.cs.mc.rough_iterations_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.rough_iterations), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_rough), xgc.cs.mc.rough_iterations_spin, 2, 3, row_material_rough, row_material_rough + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.rough_iterations_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_rough++;

    xgc.cs.mc.label_rough_seed = gtk_label_new("Generator seed (-1 for random):");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_rough_seed), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_rough), xgc.cs.mc.label_rough_seed, 0, 2, row_material_rough, row_material_rough + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.rough_seed = gtk_adjustment_new(MATERIAL_ROUGHNESS_SEED, -1, 1000, 1, 10, 0);
    xgc.cs.mc.rough_seed_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.rough_seed), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_rough), xgc.cs.mc.rough_seed_spin, 2, 3, row_material_rough, row_material_rough + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.rough_seed_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_rough++;

    xgc.cs.mc.label_rough_radiuspeak = gtk_label_new("Radius peak:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_rough_radiuspeak), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_rough), xgc.cs.mc.label_rough_radiuspeak, 0, 2, row_material_rough, row_material_rough + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.rough_radiuspeak = gtk_adjustment_new(MATERIAL_ROUGHNESS_RADIUSPEAK, 0, 1000, 1, 10, 0);
    xgc.cs.mc.rough_radiuspeak_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.rough_radiuspeak), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_rough), xgc.cs.mc.rough_radiuspeak_spin, 2, 3, row_material_rough, row_material_rough + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.rough_radiuspeak_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_rough++;

    xgc.cs.mc.label_rough_radiusspan = gtk_label_new("Radius span:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_rough_radiusspan), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_rough), xgc.cs.mc.label_rough_radiusspan, 0, 2, row_material_rough, row_material_rough + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.rough_radiusspan = gtk_adjustment_new(MATERIAL_ROUGHNESS_RADIUSSPAN, 0, 1000, 1, 10, 0);
    xgc.cs.mc.rough_radiusspan_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.rough_radiusspan), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_rough), xgc.cs.mc.rough_radiusspan_spin, 2, 3, row_material_rough, row_material_rough + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.rough_radiusspan_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_rough++;

    xgc.cs.mc.label_rough_probability = gtk_label_new("Addition probability:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_rough_probability), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_rough), xgc.cs.mc.label_rough_probability, 0, 2, row_material_rough, row_material_rough + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.rough_probability = gtk_adjustment_new(MATERIAL_ROUGHNESS_PROBABILITY, 0, 1.0, 0.01, 0.1, 0);
    xgc.cs.mc.rough_probability_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.rough_probability), 1, 3);
    gtk_table_attach(GTK_TABLE(table_material_rough), xgc.cs.mc.rough_probability_spin, 2, 3, row_material_rough, row_material_rough + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.rough_probability_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_rough++;

    /* V. 4. Add spectral modifier */

    create_vbox_group(&xgc.cs.mc.vbg_outer_material_spectral, &table_material_spectral, vbox_page, NULL, 4, 2);

    row_material_spectral = 0;

    xgc.cs.mc.label_spectral_sigma = gtk_label_new("Variance [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_spectral_sigma), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_spectral), xgc.cs.mc.label_spectral_sigma, 0, 1, row_material_spectral, row_material_spectral + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.spectral_sigma = gtk_adjustment_new(MATERIAL_SPECTRAL_SIGMA, 0, 1000.0f, 1, 10, 0);
    xgc.cs.mc.spectral_sigma_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.spectral_sigma), 1, 1);
    gtk_table_attach(GTK_TABLE(table_material_spectral), xgc.cs.mc.spectral_sigma_spin, 1, 2, row_material_spectral, row_material_spectral + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.spectral_sigma_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_spectral++;

    xgc.cs.mc.label_spectral_t = gtk_label_new("Correlation length [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_spectral_t), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_spectral), xgc.cs.mc.label_spectral_t, 0, 1, row_material_spectral, row_material_spectral + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.spectral_t = gtk_adjustment_new(MATERIAL_SPECTRAL_T, 0.1f, 1000.0f, 1, 10, 0);
    xgc.cs.mc.spectral_t_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.spectral_t), 1, 1);
    gtk_table_attach(GTK_TABLE(table_material_spectral), xgc.cs.mc.spectral_t_spin, 1, 2, row_material_spectral, row_material_spectral + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.spectral_t_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_spectral++;

    xgc.cs.mc.label_spectral_matindex = gtk_label_new("Material index:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_spectral_matindex), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_spectral), xgc.cs.mc.label_spectral_matindex, 0, 1, row_material_spectral, row_material_spectral + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.spectral_matindex = gtk_adjustment_new(MATERIAL_SPECTRAL_MATINDEX, 0, 100, 1, 10, 0);
    xgc.cs.mc.spectral_matindex_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.spectral_matindex), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_spectral), xgc.cs.mc.spectral_matindex_spin, 1, 2, row_material_spectral, row_material_spectral + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.spectral_matindex_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_spectral++;

    xgc.cs.mc.label_spectral_seed = gtk_label_new("Generator seed (-1 for random):");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_spectral_seed), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_spectral), xgc.cs.mc.label_spectral_seed, 0, 1, row_material_spectral, row_material_spectral + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.spectral_seed = gtk_adjustment_new(MATERIAL_SPECTRAL_SEED, -1, 1000, 1, 10, 0);
    xgc.cs.mc.spectral_seed_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.spectral_seed), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_spectral), xgc.cs.mc.spectral_seed_spin, 1, 2, row_material_spectral, row_material_spectral + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.spectral_seed_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_spectral++;

    /* V. 5. Add expression modifier */

    create_vbox_group(&xgc.cs.mc.vbg_outer_material_expr, &table_material_expr, vbox_page, NULL, 9, 3);

    row_material_expr = 0;

    xgc.cs.mc.label_expr_span_from = gtk_label_new("Generator span from (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_expr_span_from), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.label_expr_span_from, 0, 3, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_material_expr++;

    xgc.cs.mc.expr_i0 = gtk_adjustment_new(MATERIAL_EXPR_I0, 1, 10000, 1, 10, 0);
    xgc.cs.mc.expr_i0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.expr_i0), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.mc.expr_i0_spin), 0);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.expr_i0_spin, 0, 1, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.expr_i0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.mc.expr_j0 = gtk_adjustment_new(MATERIAL_EXPR_I0, 1, 10000, 1, 10, 0);
    xgc.cs.mc.expr_j0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.expr_j0), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.mc.expr_j0_spin), 0);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.expr_j0_spin, 1, 2, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.expr_j0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.mc.expr_k0 = gtk_adjustment_new(MATERIAL_EXPR_I0, 1, 10000, 1, 10, 0);
    xgc.cs.mc.expr_k0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.expr_k0), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.mc.expr_k0_spin), 0);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.expr_k0_spin, 2, 3, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.expr_k0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_expr++;

    xgc.cs.mc.label_expr_span_to = gtk_label_new("Generator span to (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_expr_span_to), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.label_expr_span_to, 0, 3, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_material_expr++;

    xgc.cs.mc.expr_in = gtk_adjustment_new(MATERIAL_EXPR_I0, 1, 10000, 1, 10, 0);
    xgc.cs.mc.expr_in_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.expr_in), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.mc.expr_in_spin), 0);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.expr_in_spin, 0, 1, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.expr_in_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.mc.expr_jn = gtk_adjustment_new(MATERIAL_EXPR_I0, 1, 10000, 1, 10, 0);
    xgc.cs.mc.expr_jn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.expr_jn), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.mc.expr_jn_spin), 0);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.expr_jn_spin, 1, 2, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.expr_jn_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.mc.expr_kn = gtk_adjustment_new(MATERIAL_EXPR_I0, 1, 10000, 1, 10, 0);
    xgc.cs.mc.expr_kn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.expr_kn), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.mc.expr_kn_spin), 0);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.expr_kn_spin, 2, 3, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.expr_kn_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_expr++;

    xgc.cs.mc.label_expr_matindex = gtk_label_new("Material index:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_expr_matindex), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.label_expr_matindex, 0, 2, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.expr_matindex = gtk_adjustment_new(MATERIAL_EXPR_MATINDEX, 0, 100, 1, 10, 0);
    xgc.cs.mc.expr_matindex_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.expr_matindex), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.expr_matindex_spin, 2, 3, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.expr_matindex_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_expr++;

    xgc.cs.mc.label_expr_voidindex = gtk_label_new("Void index:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_expr_voidindex), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.label_expr_voidindex, 0, 2, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.expr_voidindex = gtk_adjustment_new(MATERIAL_EXPR_MATINDEX, 0, 100, 1, 10, 0);
    xgc.cs.mc.expr_voidindex_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.expr_voidindex), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.expr_voidindex_spin, 2, 3, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.expr_voidindex_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_expr++;

    xgc.cs.mc.label_expr_maxdist = gtk_label_new("Maximum distance [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_expr_maxdist), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.label_expr_maxdist, 0, 2, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.expr_maxdist = gtk_adjustment_new(MATERIAL_EXPR_MATINDEX, 0.0, 1000, 1, 10, 0);
    xgc.cs.mc.expr_maxdist_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.expr_maxdist), 1, 1);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.expr_maxdist_spin, 2, 3, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.expr_maxdist_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_expr++;

    xgc.cs.mc.label_expr_distmode = gtk_label_new("Distance mode:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_expr_distmode), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.label_expr_distmode, 0, 2, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.mc.expr_distmode = gtk_adjustment_new(MATERIAL_EXPR_MATINDEX, 0, 10, 1, 10, 0);
    xgc.cs.mc.expr_distmode_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.mc.expr_distmode), 1, 0);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.expr_distmode_spin, 2, 3, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.expr_distmode_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_material_expr++;

    xgc.cs.mc.label_expr_expr = gtk_label_new("Expression:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.mc.label_expr_expr), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.label_expr_expr, 0, 2, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_material_expr++;

    xgc.cs.mc.expr_expr = gtk_entry_new();
    //gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.mc.expr_expr), MATERIAL_EXPR_EXPR_CHARS);
    gtk_entry_set_text(GTK_ENTRY(xgc.cs.mc.expr_expr), MATERIAL_EXPR_EXPR);
    gtk_entry_set_editable(GTK_ENTRY(xgc.cs.mc.expr_expr), TRUE);
    gtk_widget_set_can_focus(GTK_WIDGET(xgc.cs.mc.expr_expr), TRUE);
    gtk_table_attach(GTK_TABLE(table_material_expr), xgc.cs.mc.expr_expr, 0, 3, row_material_expr, row_material_expr + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.mc.expr_expr, "changed", G_CALLBACK(par_controls_changed_cb), &xgc);


    /* VI. Outputs */

    /* VI. 1. General output */

    create_vbox_group(&xgc.cs.oc.vbg_outer_general_output, &table_general_output, vbox_page, NULL, 2, 2);

    row_general_output = 0;

    xgc.cs.oc.label_general_output_filename = gtk_label_new("General output file:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_general_output_filename), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_general_output), xgc.cs.oc.label_general_output_filename, 0, 1, row_general_output, row_general_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.general_output_filename = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.oc.general_output_filename), ENTRY_FILENAME_CHARS);
    gtk_entry_set_text(GTK_ENTRY(xgc.cs.oc.general_output_filename), GENERAL_FILENAME);
    gtk_entry_set_editable(GTK_ENTRY(xgc.cs.oc.general_output_filename), FALSE);
    gtk_widget_set_can_focus(GTK_WIDGET(xgc.cs.oc.general_output_filename), TRUE);
    gtk_table_attach(GTK_TABLE(table_general_output), xgc.cs.oc.general_output_filename, 1, 2, row_general_output, row_general_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_general_output++;

    xgc.cs.oc.general_output_filename_button_change = gtk_button_new_with_label("Change..");
    gtk_table_attach(GTK_TABLE(table_general_output), xgc.cs.oc.general_output_filename_button_change, 1, 2, row_general_output, row_general_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect(xgc.cs.oc.general_output_filename_button_change, "clicked", G_CALLBACK(so_choose_output_file_cb), &xgc);

    row_general_output++;


    /* VI. 2. Point output properties */

    create_vbox_group(&xgc.cs.oc.vbg_outer_point_output, &table_point_output, vbox_page, NULL, 5, 2);

    row_point_output = 0;

    xgc.cs.oc.label_point_output_origin_position = gtk_label_new("Output position (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_point_output_origin_position), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_point_output), xgc.cs.oc.label_point_output_origin_position, 0, 2, row_point_output, row_point_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_point_output++;

    hbox_point_output_pos = gtk_hbox_new(FALSE, 0);    

    xgc.cs.oc.point_output_origin_position_i = gtk_adjustment_new(POINT_OUTPUT_ORIGIN_I, 1, 10000, 1, 10, 0);
    xgc.cs.oc.point_output_origin_position_i_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.point_output_origin_position_i), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.point_output_origin_position_i_spin), 0);
    //gtk_table_attach(GTK_TABLE(table_point_output), xgc.cs.oc.point_output_origin_position_i_spin, 0, 1, row_point_output, row_point_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_box_pack_start(GTK_BOX(hbox_point_output_pos), xgc.cs.oc.point_output_origin_position_i_spin, TRUE, TRUE, 2);
    g_signal_connect_swapped(xgc.cs.oc.point_output_origin_position_i_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.oc.point_output_origin_position_j = gtk_adjustment_new(POINT_OUTPUT_ORIGIN_J, 1, 10000, 1, 10, 0);
    xgc.cs.oc.point_output_origin_position_j_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.point_output_origin_position_j), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.point_output_origin_position_j_spin), 0);
    //gtk_table_attach(GTK_TABLE(table_point_output), xgc.cs.oc.point_output_origin_position_j_spin, 1, 2, row_point_output, row_point_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_box_pack_start(GTK_BOX(hbox_point_output_pos), xgc.cs.oc.point_output_origin_position_j_spin, TRUE, TRUE, 2);
    g_signal_connect_swapped(xgc.cs.oc.point_output_origin_position_j_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.oc.point_output_origin_position_k = gtk_adjustment_new(POINT_OUTPUT_ORIGIN_K, 1, 10000, 1, 10, 0);
    xgc.cs.oc.point_output_origin_position_k_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.point_output_origin_position_k), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.point_output_origin_position_k_spin), 0);
    //gtk_table_attach(GTK_TABLE(table_point_output), xgc.cs.oc.point_output_origin_position_k_spin, 2, 3, row_point_output, row_point_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    gtk_box_pack_start(GTK_BOX(hbox_point_output_pos), xgc.cs.oc.point_output_origin_position_k_spin, TRUE, TRUE, 2);
    g_signal_connect_swapped(xgc.cs.oc.point_output_origin_position_k_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    gtk_table_attach(GTK_TABLE(table_point_output), hbox_point_output_pos, 0, 2, row_point_output, row_point_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_point_output++;

    xgc.cs.oc.label_point_output_step = gtk_label_new("Output step:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_point_output_step), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_point_output), xgc.cs.oc.label_point_output_step, 0, 1, row_point_output, row_point_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.point_output_step = gtk_adjustment_new(OUTPUT_STEP, 0, 1000, 1, 10, 0);
    xgc.cs.oc.point_output_step_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.point_output_step), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.point_output_step_spin), 0);
    gtk_table_attach(GTK_TABLE(table_point_output), xgc.cs.oc.point_output_step_spin, 1, 2, row_point_output, row_point_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.point_output_step_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_point_output++;

    xgc.cs.oc.label_point_output_component = gtk_label_new("Output component:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_point_output_component), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_point_output), xgc.cs.oc.label_point_output_component, 0, 1, row_point_output, row_point_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.point_output_component = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.point_output_component), "Ex");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.point_output_component), "Ey");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.point_output_component), "Ez");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.point_output_component), "Hx");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.point_output_component), "Hy");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.point_output_component), "Hz");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.point_output_component), "All");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.point_output_component), "Currents");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.point_output_component), "Epsilon");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.point_output_component), "Sigma");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.point_output_component), "Mu");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.point_output_component), "Sigast");
    gtk_combo_box_set_active(GTK_COMBO_BOX(xgc.cs.oc.point_output_component), OUTPUT_COMPONENT);
    gtk_table_attach(GTK_TABLE(table_point_output), xgc.cs.oc.point_output_component, 1, 2, row_point_output, row_point_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.point_output_component, "changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_point_output++;

    xgc.cs.oc.label_point_output_filename = gtk_label_new("Output file:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_point_output_filename), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_point_output), xgc.cs.oc.label_point_output_filename, 0, 1, row_point_output, row_point_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    //row_point_output++;

    xgc.cs.oc.point_output_filename = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.oc.point_output_filename), ENTRY_FILENAME_CHARS);
    gtk_entry_set_text(GTK_ENTRY(xgc.cs.oc.point_output_filename), GENERAL_FILENAME);
    gtk_entry_set_editable(GTK_ENTRY(xgc.cs.oc.point_output_filename), FALSE);
    gtk_widget_set_can_focus(GTK_WIDGET(xgc.cs.oc.point_output_filename), TRUE);
    gtk_table_attach(GTK_TABLE(table_point_output), xgc.cs.oc.point_output_filename, 1, 2, row_point_output, row_point_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_point_output++;

    xgc.cs.oc.point_output_filename_button_change = gtk_button_new_with_label("Change..");
    gtk_table_attach(GTK_TABLE(table_point_output), xgc.cs.oc.point_output_filename_button_change, 1, 2, row_point_output, row_point_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect(xgc.cs.oc.point_output_filename_button_change, "clicked", G_CALLBACK(so_choose_output_file_cb), &xgc);

    row_point_output++;

    /* VI. 3. Image output properties */

    create_vbox_group(&xgc.cs.oc.vbg_outer_image_output, &table_image_output, vbox_page, NULL, 5, 2);

    row_image_output = 0;

    xgc.cs.oc.label_image_output_plane = gtk_label_new("Output plane:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_image_output_plane), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_image_output), xgc.cs.oc.label_image_output_plane, 0, 1, row_image_output, row_image_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.image_output_plane = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_plane), "YZ plane");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_plane), "XZ plane");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_plane), "XY plane");
    gtk_combo_box_set_active(GTK_COMBO_BOX(xgc.cs.oc.image_output_plane), OUTPUT_PLANE);
    gtk_table_attach(GTK_TABLE(table_image_output), xgc.cs.oc.image_output_plane, 1, 2, row_image_output, row_image_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.image_output_plane, "changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_image_output++;

    xgc.cs.oc.label_image_output_pos = gtk_label_new("Output position [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_image_output_pos), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_image_output), xgc.cs.oc.label_image_output_pos, 0, 1, row_image_output, row_image_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.image_output_pos = gtk_adjustment_new(OUTPUT_IMAGE_POSI, 0, 10000, 1, 10, 0);
    xgc.cs.oc.image_output_pos_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.image_output_pos), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.image_output_pos_spin), 0);
    gtk_table_attach(GTK_TABLE(table_image_output), xgc.cs.oc.image_output_pos_spin, 1, 2, row_image_output, row_image_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.image_output_pos_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_image_output++;

    xgc.cs.oc.label_image_output_step = gtk_label_new("Output step:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_image_output_step), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_image_output), xgc.cs.oc.label_image_output_step, 0, 1, row_image_output, row_image_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.image_output_step = gtk_adjustment_new(OUTPUT_STEP, 0, 1000, 1, 10, 0);
    xgc.cs.oc.image_output_step_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.image_output_step), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.image_output_step_spin), 0);
    gtk_table_attach(GTK_TABLE(table_image_output), xgc.cs.oc.image_output_step_spin, 1, 2, row_image_output, row_image_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.image_output_step_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_image_output++;

    xgc.cs.oc.label_image_output_component = gtk_label_new("Output component:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_image_output_component), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_image_output), xgc.cs.oc.label_image_output_component, 0, 1, row_image_output, row_image_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.image_output_component = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_component), "Ex");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_component), "Ey");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_component), "Ez");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_component), "Hx");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_component), "Hy");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_component), "Hz");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_component), "All");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_component), "Currents");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_component), "Epsilon");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_component), "Sigma");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_component), "Mu");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_component), "Sigast");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.image_output_component), "Material");
    gtk_combo_box_set_active(GTK_COMBO_BOX(xgc.cs.oc.image_output_component), OUTPUT_COMPONENT);
    gtk_table_attach(GTK_TABLE(table_image_output), xgc.cs.oc.image_output_component, 1, 2, row_image_output, row_image_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.image_output_component, "changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_image_output++;

    xgc.cs.oc.label_image_output_label = gtk_label_new("Output label:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_image_output_label), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_image_output), xgc.cs.oc.label_image_output_label, 0, 1, row_image_output, row_image_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.image_output_label = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.oc.image_output_label), ENTRY_FILENAME_CHARS);
    gtk_entry_set_text(GTK_ENTRY(xgc.cs.oc.image_output_label), GENERAL_FILENAME);
    gtk_entry_set_editable(GTK_ENTRY(xgc.cs.oc.image_output_label), TRUE);
    gtk_widget_set_can_focus(GTK_WIDGET(xgc.cs.oc.image_output_label), TRUE);
    gtk_table_attach(GTK_TABLE(table_image_output), xgc.cs.oc.image_output_label, 1, 2, row_image_output, row_image_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    //g_signal_connect_swapped(xgc.cs.oc.image_output_label, "activate", G_CALLBACK(controls_changed), &xgc);
    g_signal_connect_swapped(xgc.cs.oc.image_output_label, "changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_image_output++;


    /* VI. 4. Plane output properties */

    create_vbox_group(&xgc.cs.oc.vbg_outer_plane_output, &table_plane_output, vbox_page, NULL, 9, 2);

    row_plane_output = 0;

    xgc.cs.oc.label_plane_output_plane = gtk_label_new("Output plane:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_plane_output_plane), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.label_plane_output_plane, 0, 1, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.plane_output_plane = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_plane), "YZ plane");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_plane), "XZ plane");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_plane), "XY plane");
    gtk_combo_box_set_active(GTK_COMBO_BOX(xgc.cs.oc.plane_output_plane), OUTPUT_PLANE);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.plane_output_plane, 1, 2, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.plane_output_plane, "changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_plane_output++;

    xgc.cs.oc.label_plane_output_pos = gtk_label_new("Output position [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_plane_output_pos), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.label_plane_output_pos, 0, 1, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.plane_output_pos = gtk_adjustment_new(OUTPUT_PLANE_POSI, 0, 10000, 1, 10, 0);
    xgc.cs.oc.plane_output_pos_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.plane_output_pos), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.plane_output_pos_spin), 0);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.plane_output_pos_spin, 1, 2, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.plane_output_pos_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_plane_output++;

    xgc.cs.oc.label_plane_output_step = gtk_label_new("Output step:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_plane_output_step), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.label_plane_output_step, 0, 1, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.plane_output_step = gtk_adjustment_new(OUTPUT_STEP, 0, 1000, 1, 10, 0);
    xgc.cs.oc.plane_output_step_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.plane_output_step), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.plane_output_step_spin), 0);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.plane_output_step_spin, 1, 2, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.plane_output_step_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_plane_output++;

    xgc.cs.oc.label_plane_output_component = gtk_label_new("Output component:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_plane_output_component), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.label_plane_output_component, 0, 1, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.plane_output_component = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_component), "Ex");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_component), "Ey");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_component), "Ez");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_component), "Hx");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_component), "Hy");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_component), "Hz");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_component), "All");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_component), "Currents");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_component), "Epsilon");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_component), "Sigma");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_component), "Mu");
    gtk_combo_box_set_active(GTK_COMBO_BOX(xgc.cs.oc.plane_output_component), OUTPUT_COMPONENT);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.plane_output_component, 1, 2, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.plane_output_component, "changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_plane_output++;

    xgc.cs.oc.label_plane_output_filename = gtk_label_new("Output file:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_plane_output_filename), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.label_plane_output_filename, 0, 1, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.plane_output_filename = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.oc.plane_output_filename), ENTRY_FILENAME_CHARS);
    gtk_entry_set_text(GTK_ENTRY(xgc.cs.oc.plane_output_filename), GENERAL_FILENAME);
    gtk_entry_set_editable(GTK_ENTRY(xgc.cs.oc.plane_output_filename), FALSE);
    gtk_widget_set_can_focus(GTK_WIDGET(xgc.cs.oc.plane_output_filename), TRUE);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.plane_output_filename, 1, 2, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_plane_output++;

    xgc.cs.oc.plane_output_filename_button_change = gtk_button_new_with_label("Change..");
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.plane_output_filename_button_change, 1, 2, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect(xgc.cs.oc.plane_output_filename_button_change, "clicked", G_CALLBACK(so_choose_output_file_cb), &xgc);

    row_plane_output++;

    xgc.cs.oc.label_plane_output_format = gtk_label_new("Output format:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_plane_output_format), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.label_plane_output_format, 0, 1, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.plane_output_format = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_format), "Binary");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.plane_output_format), "Text ACSII");
    gtk_combo_box_set_active(GTK_COMBO_BOX(xgc.cs.oc.plane_output_format), OUTPUT_FORMAT);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.plane_output_format, 1, 2, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.plane_output_format, "changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_plane_output++;

    xgc.cs.oc.label_plane_output_start = gtk_label_new("Start at step:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_plane_output_start), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.label_plane_output_start, 0, 1, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.plane_output_start = gtk_adjustment_new(OUTPUT_START, 0, 100000, 1, 10, 0);
    xgc.cs.oc.plane_output_start_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.plane_output_start), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.plane_output_start_spin), 0);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.plane_output_start_spin, 1, 2, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.plane_output_start_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_plane_output++;

    xgc.cs.oc.label_plane_output_stop = gtk_label_new("Stop at step:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_plane_output_stop), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.label_plane_output_stop, 0, 1, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.plane_output_stop = gtk_adjustment_new(OUTPUT_STOP, 0, 100000, 1, 10, 0);
    xgc.cs.oc.plane_output_stop_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.plane_output_stop), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.plane_output_stop_spin), 0);
    gtk_table_attach(GTK_TABLE(table_plane_output), xgc.cs.oc.plane_output_stop_spin, 1, 2, row_plane_output, row_plane_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.plane_output_stop_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_plane_output++;

    /* B5. Volume output properties */

    create_vbox_group(&xgc.cs.oc.vbg_outer_volume_output, &table_volume_output, vbox_page, NULL, 7, 2);

    row_volume_output = 0;

    xgc.cs.oc.label_volume_output_step = gtk_label_new("Output step:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_volume_output_step), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_volume_output), xgc.cs.oc.label_volume_output_step, 0, 1, row_volume_output, row_volume_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.volume_output_step = gtk_adjustment_new(OUTPUT_STEP, 0, 1000, 1, 10, 0);
    xgc.cs.oc.volume_output_step_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.volume_output_step), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.volume_output_step_spin), 0);
    gtk_table_attach(GTK_TABLE(table_volume_output), xgc.cs.oc.volume_output_step_spin, 1, 2, row_volume_output, row_volume_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.volume_output_step_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_volume_output++;

    xgc.cs.oc.label_volume_output_component = gtk_label_new("Output component:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_volume_output_component), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_volume_output), xgc.cs.oc.label_volume_output_component, 0, 1, row_volume_output, row_volume_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.volume_output_component = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.volume_output_component), "Ex");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.volume_output_component), "Ey");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.volume_output_component), "Ez");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.volume_output_component), "Hx");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.volume_output_component), "Hy");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.volume_output_component), "Hz");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.volume_output_component), "All");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.volume_output_component), "Epsilon");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.volume_output_component), "Sigma");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.volume_output_component), "Mu");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.volume_output_component), "Material");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.volume_output_component), "Mat. type");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.volume_output_component), "Absorption");
    gtk_combo_box_set_active(GTK_COMBO_BOX(xgc.cs.oc.volume_output_component), OUTPUT_COMPONENT);
    gtk_table_attach(GTK_TABLE(table_volume_output), xgc.cs.oc.volume_output_component, 1, 2, row_volume_output, row_volume_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.volume_output_component, "changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_volume_output++;

    xgc.cs.oc.label_volume_output_filename = gtk_label_new("Output file:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_volume_output_filename), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_volume_output), xgc.cs.oc.label_volume_output_filename, 0, 1, row_volume_output, row_volume_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.volume_output_filename = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.oc.volume_output_filename), ENTRY_FILENAME_CHARS);
    gtk_entry_set_text(GTK_ENTRY(xgc.cs.oc.volume_output_filename), GENERAL_FILENAME);
    gtk_entry_set_editable(GTK_ENTRY(xgc.cs.oc.volume_output_filename), FALSE);
    gtk_widget_set_can_focus(GTK_WIDGET(xgc.cs.oc.volume_output_filename), TRUE);
    gtk_table_attach(GTK_TABLE(table_volume_output), xgc.cs.oc.volume_output_filename, 1, 2, row_volume_output, row_volume_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_volume_output++;

    xgc.cs.oc.volume_output_filename_button_change = gtk_button_new_with_label("Change..");
    gtk_table_attach(GTK_TABLE(table_volume_output), xgc.cs.oc.volume_output_filename_button_change, 1, 2, row_volume_output, row_volume_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect(xgc.cs.oc.volume_output_filename_button_change, "clicked", G_CALLBACK(so_choose_output_file_cb), &xgc);

    row_volume_output++;

    xgc.cs.oc.label_volume_output_format = gtk_label_new("Output format:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_volume_output_format), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_volume_output), xgc.cs.oc.label_volume_output_format, 0, 1, row_volume_output, row_volume_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.volume_output_format = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.volume_output_format), "Binary");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.volume_output_format), "Text ACSII");
    gtk_combo_box_set_active(GTK_COMBO_BOX(xgc.cs.oc.volume_output_format), OUTPUT_FORMAT);
    gtk_table_attach(GTK_TABLE(table_volume_output), xgc.cs.oc.volume_output_format, 1, 2, row_volume_output, row_volume_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.volume_output_format, "changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_volume_output++;

    xgc.cs.oc.label_volume_output_start = gtk_label_new("Start at step:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_volume_output_start), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_volume_output), xgc.cs.oc.label_volume_output_start, 0, 1, row_volume_output, row_volume_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.volume_output_start = gtk_adjustment_new(OUTPUT_START, 0, 100000, 1, 10, 0);
    xgc.cs.oc.volume_output_start_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.volume_output_start), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.volume_output_start_spin), 0);
    gtk_table_attach(GTK_TABLE(table_volume_output), xgc.cs.oc.volume_output_start_spin, 1, 2, row_volume_output, row_volume_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.volume_output_start_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_volume_output++;

    xgc.cs.oc.label_volume_output_stop = gtk_label_new("Stop at step:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_volume_output_stop), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_volume_output), xgc.cs.oc.label_volume_output_stop, 0, 1, row_volume_output, row_volume_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.volume_output_stop = gtk_adjustment_new(OUTPUT_STOP, 0, 100000, 1, 10, 0);
    xgc.cs.oc.volume_output_stop_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.volume_output_stop), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.volume_output_stop_spin), 0);
    gtk_table_attach(GTK_TABLE(table_volume_output), xgc.cs.oc.volume_output_stop_spin, 1, 2, row_volume_output, row_volume_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.volume_output_stop_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_volume_output++;

    /* VI. A. Box output group */

    create_vbox_group(&xgc.cs.oc.vbg_outer_box_output, &table_box_output, vbox_page, "Box:", 2, 4);

    row_box_output = 0;

    xgc.cs.oc.label_box_output_start = gtk_label_new("Start vertex (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_box_output_start), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_box_output), xgc.cs.oc.label_box_output_start, 0, 3, row_box_output, row_box_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_box_output++;

    xgc.cs.oc.box_output_i0 = gtk_adjustment_new(SOURCE_BOX_VERTEX0, 1, 10000, 1, 10, 0);
    xgc.cs.oc.box_output_i0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.box_output_i0), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.box_output_i0_spin), 0);
    gtk_table_attach(GTK_TABLE(table_box_output), xgc.cs.oc.box_output_i0_spin, 0, 1, row_box_output, row_box_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.box_output_i0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.oc.box_output_j0 = gtk_adjustment_new(SOURCE_BOX_VERTEX0, 1, 10000, 1, 10, 0);
    xgc.cs.oc.box_output_j0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.box_output_j0), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.box_output_j0_spin), 0);
    gtk_table_attach(GTK_TABLE(table_box_output), xgc.cs.oc.box_output_j0_spin, 1, 2, row_box_output, row_box_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.box_output_j0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.oc.box_output_k0 = gtk_adjustment_new(SOURCE_BOX_VERTEX0, 1, 10000, 1, 10, 0);
    xgc.cs.oc.box_output_k0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.box_output_k0), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.box_output_k0_spin), 0);
    gtk_table_attach(GTK_TABLE(table_box_output), xgc.cs.oc.box_output_k0_spin, 2, 3, row_box_output, row_box_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.box_output_k0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_box_output++;

    xgc.cs.oc.label_box_output_to = gtk_label_new("End vertex (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_box_output_to), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_box_output), xgc.cs.oc.label_box_output_to, 0, 3, row_box_output, row_box_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_box_output++;

    xgc.cs.oc.box_output_in = gtk_adjustment_new(SOURCE_BOX_VERTEXN, 1, 1000000, 1, 10, 0);
    xgc.cs.oc.box_output_in_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.box_output_in), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.box_output_in_spin), 0);
    gtk_table_attach(GTK_TABLE(table_box_output), xgc.cs.oc.box_output_in_spin, 0, 1, row_box_output, row_box_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.box_output_in_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.oc.box_output_jn = gtk_adjustment_new(SOURCE_BOX_VERTEXN, 1, 1000000, 1, 10, 0);
    xgc.cs.oc.box_output_jn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.box_output_jn), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.box_output_jn_spin), 0);
    gtk_table_attach(GTK_TABLE(table_box_output), xgc.cs.oc.box_output_jn_spin, 1, 2, row_box_output, row_box_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.box_output_jn_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.oc.box_output_kn = gtk_adjustment_new(SOURCE_BOX_VERTEXN, 1, 1000000, 1, 10, 0);
    xgc.cs.oc.box_output_kn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.box_output_kn), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.box_output_kn_spin), 0);
    gtk_table_attach(GTK_TABLE(table_box_output), xgc.cs.oc.box_output_kn_spin, 2, 3, row_box_output, row_box_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.box_output_kn_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_box_output++;

    xgc.cs.oc.label_box_output_size_um = gtk_label_new("Total size (x, y, z) [μm]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_box_output_size_um), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_box_output), xgc.cs.oc.label_box_output_size_um, 0, 3, row_box_output, row_box_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_box_output++;

    xgc.cs.oc.box_output_size_x_um = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.oc.box_output_size_x_um), 5);
    gtk_widget_set_sensitive(xgc.cs.oc.box_output_size_x_um, FALSE);
    gtk_table_attach(GTK_TABLE(table_box_output), xgc.cs.oc.box_output_size_x_um, 0, 1, row_box_output, row_box_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.box_output_size_y_um = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.oc.box_output_size_y_um), 5);
    gtk_widget_set_sensitive(xgc.cs.oc.box_output_size_y_um, FALSE);
    gtk_table_attach(GTK_TABLE(table_box_output), xgc.cs.oc.box_output_size_y_um, 1, 2, row_box_output, row_box_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.box_output_size_z_um = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.oc.box_output_size_z_um), 5);
    gtk_widget_set_sensitive(xgc.cs.oc.box_output_size_z_um, FALSE);
    gtk_table_attach(GTK_TABLE(table_box_output), xgc.cs.oc.box_output_size_z_um, 2, 3, row_box_output, row_box_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_box_output++;

    /* B7. Sum output properties */

    create_vbox_group(&xgc.cs.oc.vbg_outer_sum_output, &table_sum_output, vbox_page, "Sum output properties:", 2, 8);

    row_sum_output = 0;

    xgc.cs.oc.label_sum_output_epsilon = gtk_label_new("Selection ε:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_sum_output_epsilon), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_sum_output), xgc.cs.oc.label_sum_output_epsilon, 0, 1, row_sum_output, row_sum_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs.oc.sum_output_epsilon = gtk_adjustment_new(OUTPUT_SUM_EPS, 0, 100, 1, 10, 0);
    xgc.cs.oc.sum_output_epsilon_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.sum_output_epsilon), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.sum_output_epsilon_spin), 3);
    gtk_table_attach(GTK_TABLE(table_sum_output), xgc.cs.oc.sum_output_epsilon_spin, 1, 2, row_sum_output, row_sum_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.sum_output_epsilon_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_sum_output++;

    xgc.cs.oc.label_sum_output_mu = gtk_label_new("Selection μ:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_sum_output_mu), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_sum_output), xgc.cs.oc.label_sum_output_mu, 0, 1, row_sum_output, row_sum_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs.oc.sum_output_mu = gtk_adjustment_new(OUTPUT_SUM_MU, 0, 100, 1, 10, 0);
    xgc.cs.oc.sum_output_mu_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.sum_output_mu), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.sum_output_mu_spin), 3);
    gtk_table_attach(GTK_TABLE(table_sum_output), xgc.cs.oc.sum_output_mu_spin, 1, 2, row_sum_output, row_sum_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.sum_output_mu_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_sum_output++;

    xgc.cs.oc.label_sum_output_sigma = gtk_label_new("Selection σ:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_sum_output_sigma), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_sum_output), xgc.cs.oc.label_sum_output_sigma, 0, 1, row_sum_output, row_sum_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs.oc.sum_output_sigma = gtk_adjustment_new(OUTPUT_SUM_SIGMA, 0, 1e6, 1, 10, 0);
    xgc.cs.oc.sum_output_sigma_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.sum_output_sigma), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.sum_output_sigma_spin), 3);
    gtk_table_attach(GTK_TABLE(table_sum_output), xgc.cs.oc.sum_output_sigma_spin, 1, 2, row_sum_output, row_sum_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.sum_output_sigma_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_sum_output++;

    xgc.cs.oc.label_sum_output_sigast = gtk_label_new("Selection σ*:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_sum_output_sigast), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_sum_output), xgc.cs.oc.label_sum_output_sigast, 0, 1, row_sum_output, row_sum_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs.oc.sum_output_sigast = gtk_adjustment_new(OUTPUT_SUM_SIGAST, 0, 1e6, 1, 10, 0);
    xgc.cs.oc.sum_output_sigast_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.sum_output_sigast), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.sum_output_sigast_spin), 3);
    gtk_table_attach(GTK_TABLE(table_sum_output), xgc.cs.oc.sum_output_sigast_spin, 1, 2, row_sum_output, row_sum_output + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.sum_output_sigast_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_sum_output++;

    xgc.cs.oc.label_sum_output_step = gtk_label_new("Output step:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_sum_output_step), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_sum_output), xgc.cs.oc.label_sum_output_step, 0, 1, row_sum_output, row_sum_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.sum_output_step = gtk_adjustment_new(OUTPUT_STEP, 0, 1000, 1, 10, 0);
    xgc.cs.oc.sum_output_step_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.sum_output_step), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.sum_output_step_spin), 0);
    gtk_table_attach(GTK_TABLE(table_sum_output), xgc.cs.oc.sum_output_step_spin, 1, 2, row_sum_output, row_sum_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.sum_output_step_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_sum_output++;

    xgc.cs.oc.label_sum_output_component = gtk_label_new("Output component:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_sum_output_component), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_sum_output), xgc.cs.oc.label_sum_output_component, 0, 1, row_sum_output, row_sum_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.sum_output_component = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.sum_output_component), "Ex");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.sum_output_component), "Ey");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.sum_output_component), "Ez");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.sum_output_component), "All");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs.oc.sum_output_component), "Abs");
    gtk_combo_box_set_active(GTK_COMBO_BOX(xgc.cs.oc.sum_output_component), OUTPUT_COMPONENT);
    gtk_table_attach(GTK_TABLE(table_sum_output), xgc.cs.oc.sum_output_component, 1, 2, row_sum_output, row_sum_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.sum_output_component, "changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_sum_output++;

    xgc.cs.oc.label_sum_output_filename = gtk_label_new("Output file:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_sum_output_filename), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_sum_output), xgc.cs.oc.label_sum_output_filename, 0, 1, row_sum_output, row_sum_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.sum_output_filename = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.oc.sum_output_filename), ENTRY_FILENAME_CHARS);
    gtk_entry_set_text(GTK_ENTRY(xgc.cs.oc.sum_output_filename), GENERAL_FILENAME);
    gtk_entry_set_editable(GTK_ENTRY(xgc.cs.oc.sum_output_filename), FALSE);
    gtk_widget_set_can_focus(GTK_WIDGET(xgc.cs.oc.sum_output_filename), TRUE);
    gtk_table_attach(GTK_TABLE(table_sum_output), xgc.cs.oc.sum_output_filename, 1, 2, row_sum_output, row_sum_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_sum_output++;

    xgc.cs.oc.sum_output_filename_button_change = gtk_button_new_with_label("Change..");
    gtk_table_attach(GTK_TABLE(table_sum_output), xgc.cs.oc.sum_output_filename_button_change, 1, 2, row_sum_output, row_sum_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect(xgc.cs.oc.sum_output_filename_button_change, "clicked", G_CALLBACK(so_choose_output_file_cb), &xgc);

    /* B8. Force output properties */

    create_vbox_group(&xgc.cs.oc.vbg_outer_force_output, &table_force_output, vbox_page, "Force output properties:", 3, 2);

    row_force_output = 0;

    xgc.cs.oc.label_force_output_step = gtk_label_new("Output step:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_force_output_step), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_force_output), xgc.cs.oc.label_force_output_step, 0, 1, row_force_output, row_force_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.force_output_step = gtk_adjustment_new(OUTPUT_STEP, 0, 1000, 1, 10, 0);
    xgc.cs.oc.force_output_step_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.oc.force_output_step), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.oc.force_output_step_spin), 0);
    gtk_table_attach(GTK_TABLE(table_force_output), xgc.cs.oc.force_output_step_spin, 1, 2, row_force_output, row_force_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.oc.force_output_step_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_force_output++;

    xgc.cs.oc.label_force_output_filename = gtk_label_new("Output file:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.oc.label_force_output_filename), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_force_output), xgc.cs.oc.label_force_output_filename, 0, 1, row_force_output, row_force_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.oc.force_output_filename = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.oc.force_output_filename), ENTRY_FILENAME_CHARS);
    gtk_entry_set_text(GTK_ENTRY(xgc.cs.oc.force_output_filename), GENERAL_FILENAME);
    gtk_entry_set_editable(GTK_ENTRY(xgc.cs.oc.force_output_filename), FALSE);
    gtk_widget_set_can_focus(GTK_WIDGET(xgc.cs.oc.force_output_filename), TRUE);
    gtk_table_attach(GTK_TABLE(table_force_output), xgc.cs.oc.force_output_filename, 1, 4, row_force_output, row_force_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_force_output++;

    xgc.cs.oc.force_output_filename_button_change = gtk_button_new_with_label("Change..");
    gtk_table_attach(GTK_TABLE(table_force_output), xgc.cs.oc.force_output_filename_button_change, 1, 4, row_force_output, row_force_output + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect(xgc.cs.oc.force_output_filename_button_change, "clicked", G_CALLBACK(so_choose_output_file_cb), &xgc);

    /* VII. NFFF */

    /* VII. 1. Near field to far field transform box */

    create_vbox_group(&xgc.cs.fc.vbg_outer_nfff, &table_nfff, vbox_page, NULL, 18, 4);

    row_nfff = 0;

    xgc.cs.fc.label_nfff_box_from = gtk_label_new("Area span from: (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.fc.label_nfff_box_from), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.label_nfff_box_from, 0, 3, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_nfff++;

    xgc.cs.fc.nfff_box_i0 = gtk_adjustment_new(NFFF_BOX_VERTEX0, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfff_box_i0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_box_i0), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_box_i0_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_box_i0_spin, 0, 1, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_box_i0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_box_j0 = gtk_adjustment_new(NFFF_BOX_VERTEX0, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfff_box_j0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_box_j0), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_box_j0_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_box_j0_spin, 1, 2, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_box_j0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_box_k0 = gtk_adjustment_new(NFFF_BOX_VERTEX0, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfff_box_k0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_box_k0), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_box_k0_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_box_k0_spin, 2, 3, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_box_k0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfff++;

    xgc.cs.fc.label_nfff_box_to = gtk_label_new("Area span to: (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.fc.label_nfff_box_to), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.label_nfff_box_to, 0, 3, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_nfff++;

    xgc.cs.fc.nfff_box_in = gtk_adjustment_new(NFFF_BOX_VERTEXN, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfff_box_in_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_box_in), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_box_in_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_box_in_spin, 0, 1, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_box_in_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_box_jn = gtk_adjustment_new(NFFF_BOX_VERTEXN, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfff_box_jn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_box_jn), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_box_jn_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_box_jn_spin, 1, 2, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_box_jn_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_box_kn = gtk_adjustment_new(NFFF_BOX_VERTEXN, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfff_box_kn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_box_kn), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_box_kn_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_box_kn_spin, 2, 3, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_box_kn_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfff++;

    xgc.cs.fc.label_nfff_box_size_um = gtk_label_new("Total size (x, y, z) [μm]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.fc.label_nfff_box_size_um), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.label_nfff_box_size_um, 0, 3, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_nfff++;

    xgc.cs.fc.nfff_box_size_x_um = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.fc.nfff_box_size_x_um), 5);
    gtk_widget_set_sensitive(xgc.cs.fc.nfff_box_size_x_um, FALSE);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_box_size_x_um, 0, 1, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.fc.nfff_box_size_y_um = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.fc.nfff_box_size_y_um), 5);
    gtk_widget_set_sensitive(xgc.cs.fc.nfff_box_size_y_um, FALSE);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_box_size_y_um, 1, 2, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.fc.nfff_box_size_z_um = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.fc.nfff_box_size_z_um), 5);
    gtk_widget_set_sensitive(xgc.cs.fc.nfff_box_size_z_um, FALSE);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_box_size_z_um, 2, 3, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_nfff++;

    xgc.cs.fc.nfff_box_boundary_skipi0 = gtk_check_button_new_with_mnemonic("i0 skip in range (jmin, kmin, jmax, kmax):");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.fc.nfff_box_boundary_skipi0), NFFF_BOX_BOUNDARY_SKIP);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_box_boundary_skipi0, 0, 4, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_box_boundary_skipi0, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfff++;

    xgc.cs.fc.nfff_skipi0_jmin = gtk_adjustment_new(NFFF_SKIPI0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipi0_jmin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipi0_jmin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipi0_jmin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipi0_jmin_spin, 0, 1, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipi0_jmin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipi0_kmin = gtk_adjustment_new(NFFF_SKIPI0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipi0_kmin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipi0_kmin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipi0_kmin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipi0_kmin_spin, 1, 2, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipi0_kmin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipi0_jmax = gtk_adjustment_new(NFFF_SKIPI0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipi0_jmax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipi0_jmax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipi0_jmax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipi0_jmax_spin, 2, 3, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipi0_jmax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipi0_kmax = gtk_adjustment_new(NFFF_SKIPI0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipi0_kmax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipi0_kmax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipi0_kmax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipi0_kmax_spin, 3, 4, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipi0_kmax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfff++;

    xgc.cs.fc.nfff_box_boundary_skipin = gtk_check_button_new_with_mnemonic("in skip in range (jmin, kmin, jmax, kmax):");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.fc.nfff_box_boundary_skipin), NFFF_BOX_BOUNDARY_SKIP);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_box_boundary_skipin, 0, 4, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_box_boundary_skipin, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfff++;

    xgc.cs.fc.nfff_skipin_jmin = gtk_adjustment_new(NFFF_SKIPIN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipin_jmin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipin_jmin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipin_jmin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipin_jmin_spin, 0, 1, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipin_jmin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipin_kmin = gtk_adjustment_new(NFFF_SKIPIN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipin_kmin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipin_kmin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipin_kmin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipin_kmin_spin, 1, 2, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipin_kmin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipin_jmax = gtk_adjustment_new(NFFF_SKIPIN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipin_jmax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipin_jmax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipin_jmax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipin_jmax_spin, 2, 3, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipin_jmax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipin_kmax = gtk_adjustment_new(NFFF_SKIPIN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipin_kmax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipin_kmax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipin_kmax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipin_kmax_spin, 3, 4, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipin_kmax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfff++;

    xgc.cs.fc.nfff_box_boundary_skipj0 = gtk_check_button_new_with_mnemonic("j0 skip in range (imin, kmin, imax, kmax):");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.fc.nfff_box_boundary_skipj0), NFFF_BOX_BOUNDARY_SKIP);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_box_boundary_skipj0, 0, 4, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_box_boundary_skipj0, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfff++;

    xgc.cs.fc.nfff_skipj0_imin = gtk_adjustment_new(NFFF_SKIPJ0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipj0_imin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipj0_imin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipj0_imin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipj0_imin_spin, 0, 1, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipj0_imin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipj0_kmin = gtk_adjustment_new(NFFF_SKIPJ0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipj0_kmin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipj0_kmin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipj0_kmin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipj0_kmin_spin, 1, 2, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipj0_kmin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipj0_imax = gtk_adjustment_new(NFFF_SKIPJ0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipj0_imax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipj0_imax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipj0_imax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipj0_imax_spin, 2, 3, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipj0_imax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipj0_kmax = gtk_adjustment_new(NFFF_SKIPJ0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipj0_kmax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipj0_kmax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipj0_kmax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipj0_kmax_spin, 3, 4, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipj0_kmax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfff++;

    xgc.cs.fc.nfff_box_boundary_skipjn = gtk_check_button_new_with_mnemonic("jn skip in range (imin, kmin, imax, kmax):");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.fc.nfff_box_boundary_skipjn), NFFF_BOX_BOUNDARY_SKIP);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_box_boundary_skipjn, 0, 4, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_box_boundary_skipjn, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfff++;

    xgc.cs.fc.nfff_skipjn_imin = gtk_adjustment_new(NFFF_SKIPJN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipjn_imin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipjn_imin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipjn_imin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipjn_imin_spin, 0, 1, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipjn_imin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipjn_kmin = gtk_adjustment_new(NFFF_SKIPJN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipjn_kmin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipjn_kmin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipjn_kmin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipjn_kmin_spin, 1, 2, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipjn_kmin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipjn_imax = gtk_adjustment_new(NFFF_SKIPJN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipjn_imax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipjn_imax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipjn_imax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipjn_imax_spin, 2, 3, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipjn_imax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipjn_kmax = gtk_adjustment_new(NFFF_SKIPJN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipjn_kmax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipjn_kmax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipjn_kmax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipjn_kmax_spin, 3, 4, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipjn_kmax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfff++;

    xgc.cs.fc.nfff_box_boundary_skipk0 = gtk_check_button_new_with_mnemonic("k0 skip in range (imin, jmin, imax, jmax):");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.fc.nfff_box_boundary_skipk0), NFFF_BOX_BOUNDARY_SKIP);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_box_boundary_skipk0, 0, 4, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_box_boundary_skipk0, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfff++;

    xgc.cs.fc.nfff_skipk0_imin = gtk_adjustment_new(NFFF_SKIPK0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipk0_imin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipk0_imin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipk0_imin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipk0_imin_spin, 0, 1, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipk0_imin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipk0_jmin = gtk_adjustment_new(NFFF_SKIPK0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipk0_jmin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipk0_jmin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipk0_jmin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipk0_jmin_spin, 1, 2, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipk0_jmin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipk0_imax = gtk_adjustment_new(NFFF_SKIPK0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipk0_imax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipk0_imax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipk0_imax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipk0_imax_spin, 2, 3, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipk0_imax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipk0_jmax = gtk_adjustment_new(NFFF_SKIPK0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipk0_jmax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipk0_jmax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipk0_jmax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipk0_jmax_spin, 3, 4, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipk0_jmax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfff++;

    xgc.cs.fc.nfff_box_boundary_skipkn = gtk_check_button_new_with_mnemonic("kn skip in range (imin, jmin, imax, jmax):");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.fc.nfff_box_boundary_skipkn), NFFF_BOX_BOUNDARY_SKIP);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_box_boundary_skipkn, 0, 4, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_box_boundary_skipkn, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfff++;

    xgc.cs.fc.nfff_skipkn_imin = gtk_adjustment_new(NFFF_SKIPKN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipkn_imin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipkn_imin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipkn_imin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipkn_imin_spin, 0, 1, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipkn_imin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipkn_jmin = gtk_adjustment_new(NFFF_SKIPKN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipkn_jmin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipkn_jmin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipkn_jmin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipkn_jmin_spin, 1, 2, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipkn_jmin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipkn_imax = gtk_adjustment_new(NFFF_SKIPKN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipkn_imax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipkn_imax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipkn_imax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipkn_imax_spin, 2, 3, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipkn_imax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfff_skipkn_jmax = gtk_adjustment_new(NFFF_SKIPKN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.nfff_skipkn_jmax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfff_skipkn_jmax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfff_skipkn_jmax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfff), xgc.cs.fc.nfff_skipkn_jmax_spin, 3, 4, row_nfff, row_nfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfff_skipkn_jmax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfff++;

    /* VII. 2. Near field to far field point */

    create_vbox_group(&xgc.cs.fc.vbg_outer_nfffp, &table_nfffp, vbox_page, NULL, 4, 3);

    hbox_nfffp_pos = gtk_hbox_new(FALSE, 0);

    row_nfffp = 0;        

    xgc.cs.fc.label_nfffp_pos = gtk_label_new("Output position (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.fc.label_nfffp_pos), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_nfffp), xgc.cs.fc.label_nfffp_pos, 0, 3, row_nfffp, row_nfffp + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_nfffp++;

    xgc.cs.fc.nfffp_pos_i = gtk_adjustment_new(NFFFP_POS_I, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfffp_pos_i_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfffp_pos_i), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfffp_pos_i_spin), 0);
    //gtk_table_attach(GTK_TABLE(table_nfffp), xgc.cs.fc.nfffp_pos_i_spin, 0, 1, row_nfffp, row_nfffp + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_box_pack_start(GTK_BOX(hbox_nfffp_pos), xgc.cs.fc.nfffp_pos_i_spin, TRUE, TRUE, 2);
    g_signal_connect_swapped(xgc.cs.fc.nfffp_pos_i_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfffp_pos_j = gtk_adjustment_new(NFFFP_POS_J, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfffp_pos_j_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfffp_pos_j), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfffp_pos_j_spin), 0);
    //gtk_table_attach(GTK_TABLE(table_nfffp), xgc.cs.fc.nfffp_pos_j_spin, 1, 2, row_nfffp, row_nfffp + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_box_pack_start(GTK_BOX(hbox_nfffp_pos), xgc.cs.fc.nfffp_pos_j_spin, TRUE, TRUE, 2);
    g_signal_connect_swapped(xgc.cs.fc.nfffp_pos_j_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfffp_pos_k = gtk_adjustment_new(NFFFP_POS_K, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfffp_pos_k_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfffp_pos_k), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfffp_pos_k_spin), 0);
    //gtk_table_attach(GTK_TABLE(table_nfffp), xgc.cs.fc.nfffp_pos_k_spin, 2, 3, row_nfffp, row_nfffp + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_box_pack_start(GTK_BOX(hbox_nfffp_pos), xgc.cs.fc.nfffp_pos_k_spin, TRUE, TRUE, 2);
    g_signal_connect_swapped(xgc.cs.fc.nfffp_pos_k_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    gtk_table_attach(GTK_TABLE(table_nfffp), hbox_nfffp_pos, 0, 2, row_nfffp, row_nfffp + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_nfffp++;

    xgc.cs.fc.label_nfffp_filename = gtk_label_new("Output file:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.fc.label_nfffp_filename), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_nfffp), xgc.cs.fc.label_nfffp_filename, 0, 1, row_nfffp, row_nfffp + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.fc.nfffp_filename = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.fc.nfffp_filename), ENTRY_FILENAME_CHARS);
    gtk_entry_set_text(GTK_ENTRY(xgc.cs.fc.nfffp_filename), GENERAL_FILENAME);
    gtk_entry_set_editable(GTK_ENTRY(xgc.cs.fc.nfffp_filename), FALSE);
    gtk_widget_set_can_focus(GTK_WIDGET(xgc.cs.fc.nfffp_filename), TRUE);
    gtk_table_attach(GTK_TABLE(table_nfffp), xgc.cs.fc.nfffp_filename, 1, 2, row_nfffp, row_nfffp + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_nfffp++;

    xgc.cs.fc.nfffp_filename_button_change = gtk_button_new_with_label("Change..");
    gtk_table_attach(GTK_TABLE(table_nfffp), xgc.cs.fc.nfffp_filename_button_change, 1, 2, row_nfffp, row_nfffp + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect(xgc.cs.fc.nfffp_filename_button_change, "clicked", G_CALLBACK(nfffp_choose_output_file_cb), &xgc);

    row_nfffp++;

    /* VII. 3. Near field to far field area */

    create_vbox_group(&xgc.cs.fc.vbg_outer_nfffa, &table_nfffa, vbox_page, NULL, 6, 3);

    row_nfffa = 0;

    xgc.cs.fc.label_nfffa_thetares = gtk_label_new("Theta resolution [pixel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.fc.label_nfffa_thetares), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_nfffa), xgc.cs.fc.label_nfffa_thetares, 0, 1, row_nfffa, row_nfffa + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.fc.nfffa_thetares = gtk_adjustment_new(NFFFA_THETARES, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfffa_thetares_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfffa_thetares), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfffa_thetares_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfffa), xgc.cs.fc.nfffa_thetares_spin, 1, 2, row_nfffa, row_nfffa + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfffa_thetares_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfffa++;

    xgc.cs.fc.label_nfffa_phires = gtk_label_new("Phi resolution [pixel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.fc.label_nfffa_phires), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_nfffa), xgc.cs.fc.label_nfffa_phires, 0, 1, row_nfffa, row_nfffa + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.fc.nfffa_phires = gtk_adjustment_new(NFFFA_PHIRES, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfffa_phires_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfffa_phires), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfffa_phires_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfffa), xgc.cs.fc.nfffa_phires_spin, 1, 2, row_nfffa, row_nfffa + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfffa_phires_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfffa++;

    xgc.cs.fc.label_nfffa_radius = gtk_label_new("Radius [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.fc.label_nfffa_radius), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_nfffa), xgc.cs.fc.label_nfffa_radius, 0, 1, row_nfffa, row_nfffa + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.fc.nfffa_radius = gtk_adjustment_new(NFFFA_RADIUS, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfffa_radius_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfffa_radius), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfffa_radius_spin), 0);
    gtk_table_attach(GTK_TABLE(table_nfffa), xgc.cs.fc.nfffa_radius_spin, 1, 2, row_nfffa, row_nfffa + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfffa_radius_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfffa++;

    xgc.cs.fc.label_nfffa_thetafrom = gtk_label_new("Theta span [deg]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.fc.label_nfffa_thetafrom), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_nfffa), xgc.cs.fc.label_nfffa_thetafrom, 0, 1, row_nfffa, row_nfffa + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.fc.nfffa_thetafrom = gtk_adjustment_new(NFFFA_THETAFROMDEG, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfffa_thetafrom_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfffa_thetafrom), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfffa_thetafrom_spin), 3);
    gtk_table_attach(GTK_TABLE(table_nfffa), xgc.cs.fc.nfffa_thetafrom_spin, 1, 2, row_nfffa, row_nfffa + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfffa_thetafrom_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfffa_thetato = gtk_adjustment_new(NFFFA_THETATODEG, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfffa_thetato_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfffa_thetato), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfffa_thetato_spin), 3);
    gtk_table_attach(GTK_TABLE(table_nfffa), xgc.cs.fc.nfffa_thetato_spin, 2, 3, row_nfffa, row_nfffa + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfffa_thetato_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfffa++;

    xgc.cs.fc.label_nfffa_phifrom = gtk_label_new("Phi span [deg]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.fc.label_nfffa_phifrom), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_nfffa), xgc.cs.fc.label_nfffa_phifrom, 0, 1, row_nfffa, row_nfffa + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.fc.nfffa_phifrom = gtk_adjustment_new(NFFFA_PHIFROMDEG, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfffa_phifrom_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfffa_phifrom), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfffa_phifrom_spin), 3);
    gtk_table_attach(GTK_TABLE(table_nfffa), xgc.cs.fc.nfffa_phifrom_spin, 1, 2, row_nfffa, row_nfffa + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfffa_phifrom_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.nfffa_phito = gtk_adjustment_new(NFFFA_PHITODEG, 1, 10000, 1, 10, 0);
    xgc.cs.fc.nfffa_phito_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.nfffa_phito), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.nfffa_phito_spin), 3);
    gtk_table_attach(GTK_TABLE(table_nfffa), xgc.cs.fc.nfffa_phito_spin, 2, 3, row_nfffa, row_nfffa + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfffa_phito_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_nfffa++;

    xgc.cs.fc.nfffa_savefile = gtk_check_button_new_with_mnemonic("save individual outputs to file");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.fc.nfffa_savefile), NFFFA_SAVEFILE);
    gtk_table_attach(GTK_TABLE(table_nfffa), xgc.cs.fc.nfffa_savefile, 0, 2, row_nfffa, row_nfffa + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.nfffa_savefile, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    /* VIII. PNFFF */

    /* VIII. 1. Periodic near field to far field transform */

    create_vbox_group(&xgc.cs.fc.vbg_outer_pnfff, &table_pnfff, vbox_page, NULL, 14, 4);

    row_pnfff = 0;

    xgc.cs.fc.label_pnfff_box_from = gtk_label_new("Area span from: (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.fc.label_pnfff_box_from), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.label_pnfff_box_from, 0, 3, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_pnfff++;

    xgc.cs.fc.pnfff_box_i0 = gtk_adjustment_new(PNFFF_BOX_VERTEX0, 1, 10000, 1, 10, 0);
    xgc.cs.fc.pnfff_box_i0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_box_i0), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_box_i0_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_box_i0_spin, 0, 1, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_box_i0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.pnfff_box_j0 = gtk_adjustment_new(PNFFF_BOX_VERTEX0, 1, 10000, 1, 10, 0);
    xgc.cs.fc.pnfff_box_j0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_box_j0), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_box_j0_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_box_j0_spin, 1, 2, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_box_j0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.pnfff_box_k0 = gtk_adjustment_new(PNFFF_BOX_VERTEX0, 1, 10000, 1, 10, 0);
    xgc.cs.fc.pnfff_box_k0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_box_k0), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_box_k0_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_box_k0_spin, 2, 3, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_box_k0_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_pnfff++;

    xgc.cs.fc.label_pnfff_box_to = gtk_label_new("Area span to: (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.fc.label_pnfff_box_to), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.label_pnfff_box_to, 0, 3, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_pnfff++;

    xgc.cs.fc.pnfff_box_in = gtk_adjustment_new(PNFFF_BOX_VERTEXN, 1, 10000, 1, 10, 0);
    xgc.cs.fc.pnfff_box_in_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_box_in), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_box_in_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_box_in_spin, 0, 1, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_box_in_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.pnfff_box_jn = gtk_adjustment_new(PNFFF_BOX_VERTEXN, 1, 10000, 1, 10, 0);
    xgc.cs.fc.pnfff_box_jn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_box_jn), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_box_jn_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_box_jn_spin, 1, 2, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_box_jn_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.pnfff_box_kn = gtk_adjustment_new(PNFFF_BOX_VERTEXN, 1, 10000, 1, 10, 0);
    xgc.cs.fc.pnfff_box_kn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_box_kn), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_box_kn_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_box_kn_spin, 2, 3, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_box_kn_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_pnfff++;

    xgc.cs.fc.label_pnfff_box_size_um = gtk_label_new("Total size (x, y, z) [μm]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.fc.label_pnfff_box_size_um), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.label_pnfff_box_size_um, 0, 3, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_pnfff++;

    xgc.cs.fc.pnfff_box_size_x_um = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.fc.pnfff_box_size_x_um), 5);
    gtk_widget_set_sensitive(xgc.cs.fc.pnfff_box_size_x_um, FALSE);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_box_size_x_um, 0, 1, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.fc.pnfff_box_size_y_um = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.fc.pnfff_box_size_y_um), 5);
    gtk_widget_set_sensitive(xgc.cs.fc.pnfff_box_size_y_um, FALSE);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_box_size_y_um, 1, 2, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs.fc.pnfff_box_size_z_um = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs.fc.pnfff_box_size_z_um), 5);
    gtk_widget_set_sensitive(xgc.cs.fc.pnfff_box_size_z_um, FALSE);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_box_size_z_um, 2, 3, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_pnfff++;

    xgc.cs.fc.label_pnfff_integrationx = gtk_label_new("X periodic integration span:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.fc.label_pnfff_integrationx), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.label_pnfff_integrationx, 0, 3, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_pnfff++;

    xgc.cs.fc.pnfff_integration_xmin = gtk_adjustment_new(PNFFF_INTEGRATION_XMIN, -100, 100, 1, 10, 0);
    xgc.cs.fc.pnfff_integration_xmin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_integration_xmin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_integration_xmin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_integration_xmin_spin, 0, 1, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_integration_xmin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.pnfff_integration_xmax = gtk_adjustment_new(PNFFF_INTEGRATION_XMAX, -100, 100, 1, 10, 0);
    xgc.cs.fc.pnfff_integration_xmax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_integration_xmax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_integration_xmax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_integration_xmax_spin, 1, 2, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_integration_xmax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_pnfff++;

    xgc.cs.fc.label_pnfff_integrationy = gtk_label_new("Y periodic integration span:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs.fc.label_pnfff_integrationy), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.label_pnfff_integrationy, 0, 3, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_pnfff++;

    xgc.cs.fc.pnfff_integration_ymin = gtk_adjustment_new(PNFFF_INTEGRATION_YMIN, -100, 100, 1, 10, 0);
    xgc.cs.fc.pnfff_integration_ymin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_integration_ymin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_integration_ymin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_integration_ymin_spin, 0, 1, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_integration_ymin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.pnfff_integration_ymax = gtk_adjustment_new(PNFFF_INTEGRATION_YMAX, -100, 100, 1, 10, 0);
    xgc.cs.fc.pnfff_integration_ymax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_integration_ymax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_integration_ymax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_integration_ymax_spin, 1, 2, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_integration_ymax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_pnfff++;

    xgc.cs.fc.pnfff_box_boundary_skipk0 = gtk_check_button_new_with_mnemonic("k0 skip in range (imin, jmin, imax, jmax):");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.fc.pnfff_box_boundary_skipk0), PNFFF_BOX_BOUNDARY_SKIP);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_box_boundary_skipk0, 0, 4, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_box_boundary_skipk0, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_pnfff++;

    xgc.cs.fc.pnfff_skipk0_imin = gtk_adjustment_new(PNFFF_SKIPK0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.pnfff_skipk0_imin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_skipk0_imin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_skipk0_imin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_skipk0_imin_spin, 0, 1, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_skipk0_imin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.pnfff_skipk0_jmin = gtk_adjustment_new(PNFFF_SKIPK0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.pnfff_skipk0_jmin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_skipk0_jmin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_skipk0_jmin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_skipk0_jmin_spin, 1, 2, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_skipk0_jmin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.pnfff_skipk0_imax = gtk_adjustment_new(PNFFF_SKIPK0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.pnfff_skipk0_imax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_skipk0_imax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_skipk0_imax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_skipk0_imax_spin, 2, 3, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_skipk0_imax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.pnfff_skipk0_jmax = gtk_adjustment_new(PNFFF_SKIPK0, 0, 100000, 1, 10, 0);
    xgc.cs.fc.pnfff_skipk0_jmax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_skipk0_jmax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_skipk0_jmax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_skipk0_jmax_spin, 3, 4, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_skipk0_jmax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_pnfff++;

    xgc.cs.fc.pnfff_box_boundary_skipkn = gtk_check_button_new_with_mnemonic("kn skip in range (imin, jmin, imax, jmax):");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.cs.fc.pnfff_box_boundary_skipkn), PNFFF_BOX_BOUNDARY_SKIP);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_box_boundary_skipkn, 0, 4, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_box_boundary_skipkn, "clicked", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_pnfff++;

    xgc.cs.fc.pnfff_skipkn_imin = gtk_adjustment_new(PNFFF_SKIPKN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.pnfff_skipkn_imin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_skipkn_imin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_skipkn_imin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_skipkn_imin_spin, 0, 1, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_skipkn_imin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.pnfff_skipkn_jmin = gtk_adjustment_new(PNFFF_SKIPKN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.pnfff_skipkn_jmin_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_skipkn_jmin), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_skipkn_jmin_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_skipkn_jmin_spin, 1, 2, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_skipkn_jmin_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.pnfff_skipkn_imax = gtk_adjustment_new(PNFFF_SKIPKN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.pnfff_skipkn_imax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_skipkn_imax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_skipkn_imax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_skipkn_imax_spin, 2, 3, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_skipkn_imax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    xgc.cs.fc.pnfff_skipkn_jmax = gtk_adjustment_new(PNFFF_SKIPKN, 0, 100000, 1, 10, 0);
    xgc.cs.fc.pnfff_skipkn_jmax_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs.fc.pnfff_skipkn_jmax), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs.fc.pnfff_skipkn_jmax_spin), 0);
    gtk_table_attach(GTK_TABLE(table_pnfff), xgc.cs.fc.pnfff_skipkn_jmax_spin, 3, 4, row_pnfff, row_pnfff + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs.fc.pnfff_skipkn_jmax_spin, "value-changed", G_CALLBACK(par_controls_changed_cb), &xgc);

    row_pnfff++;



    //const char *build_str = "Version: " VERSION " " __DATE__ " " __TIME__;

    //////////////////////////////////////////////////////////////////////////

    /* material object page strip */
    xgc.cs_mat.vbox_page = vbox_page = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox_page), 4);

    xgc.cs_mat.sw_values = scroll_window_mat_values = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scroll_window_mat_values), GTK_SHADOW_NONE);
    //gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll), GTK_POLICY_NEVER, GTK_POLICY_ALWAYS);
    gtk_container_set_border_width(GTK_CONTAINER(scroll_window_mat_values), 0);

    //gtk_box_pack_start(GTK_BOX(page_par_tree), scroll_window_mat_values, TRUE, TRUE, 4);

    //gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroll_window), table);

    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_container_add(GTK_CONTAINER(align), vbox_page);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroll_window_mat_values), align);

    gtk_widget_set_size_request(scroll_window_mat_values, 50, -1);

    xgc.cs_mat.label_name = gtk_label_new("");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.label_name), 0.0, 0.5);
    gtk_box_pack_start(GTK_BOX(vbox_page), xgc.cs_mat.label_name, TRUE, TRUE, 4);

    /* I. Material objects */

    /* I. 1. Sphere properties */

    create_vbox_group(&xgc.cs_mat.oc.vbg_outer_sphere, &table_mat_sphere, vbox_page, "Object properties:", 2, 3);

    row_mat_sphere = 0;

    xgc.cs_mat.oc.label_sphere_center = gtk_label_new("Center position (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_sphere_center), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_sphere), xgc.cs_mat.oc.label_sphere_center, 0, 3, row_mat_sphere, row_mat_sphere + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_mat_sphere++;

    xgc.cs_mat.oc.sphere_center_i = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.sphere_center_i_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.sphere_center_i), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.sphere_center_i_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.sphere_center_i_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_sphere), xgc.cs_mat.oc.sphere_center_i_spin, 0, 1, row_mat_sphere, row_mat_sphere + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.sphere_center_i_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.sphere_center_j = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.sphere_center_j_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.sphere_center_j), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.sphere_center_j_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.sphere_center_j_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_sphere), xgc.cs_mat.oc.sphere_center_j_spin, 1, 2, row_mat_sphere, row_mat_sphere + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.sphere_center_j_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.sphere_center_k = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.sphere_center_k_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.sphere_center_k), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.sphere_center_k_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.sphere_center_k_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_sphere), xgc.cs_mat.oc.sphere_center_k_spin, 2, 3, row_mat_sphere, row_mat_sphere + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.sphere_center_k_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_sphere++;

    xgc.cs_mat.oc.label_sphere_radius = gtk_label_new("Radius [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_sphere_radius), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_sphere), xgc.cs_mat.oc.label_sphere_radius, 0, 2, row_mat_sphere, row_mat_sphere + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs_mat.oc.sphere_radius = gtk_adjustment_new(MAT_SPHERE_RADIUS, 1, 1000, 1, 10, 0);
    xgc.cs_mat.oc.sphere_radius_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.sphere_radius), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.sphere_radius_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.sphere_radius_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_sphere), xgc.cs_mat.oc.sphere_radius_spin, 2, 3, row_mat_sphere, row_mat_sphere + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.sphere_radius_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_sphere++;

    /* I. 2. Voxel properties */

    create_vbox_group(&xgc.cs_mat.oc.vbg_outer_box, &table_mat_box, vbox_page, "Object properties:", 4, 3);

    row_mat_box = 0;

    xgc.cs_mat.oc.label_box_from = gtk_label_new("Start vertex (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_box_from), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_box), xgc.cs_mat.oc.label_box_from, 0, 3, row_mat_box, row_mat_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_mat_box++;

    xgc.cs_mat.oc.box_i0 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.box_i0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.box_i0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.box_i0_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.box_i0_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_box), xgc.cs_mat.oc.box_i0_spin, 0, 1, row_mat_box, row_mat_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.box_i0_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.box_j0 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.box_j0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.box_j0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.box_j0_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.box_j0_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_box), xgc.cs_mat.oc.box_j0_spin, 1, 2, row_mat_box, row_mat_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.box_j0_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.box_k0 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.box_k0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.box_k0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.box_k0_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.box_k0_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_box), xgc.cs_mat.oc.box_k0_spin, 2, 3, row_mat_box, row_mat_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.box_k0_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_box++;

    xgc.cs_mat.oc.label_box_to = gtk_label_new("End vertex (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_box_to), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_box), xgc.cs_mat.oc.label_box_to, 0, 3, row_mat_box, row_mat_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_mat_box++;

    xgc.cs_mat.oc.box_in = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.box_in_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.box_in), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.box_in_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.box_in_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_box), xgc.cs_mat.oc.box_in_spin, 0, 1, row_mat_box, row_mat_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.box_in_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.box_jn = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.box_jn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.box_jn), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.box_jn_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.box_jn_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_box), xgc.cs_mat.oc.box_jn_spin, 1, 2, row_mat_box, row_mat_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.box_jn_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.box_kn = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.box_kn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.box_kn), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.box_kn_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.box_kn_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_box), xgc.cs_mat.oc.box_kn_spin, 2, 3, row_mat_box, row_mat_box + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.box_kn_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_box++;

    /* I. 3. Cylinder properties */

    create_vbox_group(&xgc.cs_mat.oc.vbg_outer_cylinder, &table_mat_cyl, vbox_page, "Object properties:", 2, 3);

    row_mat_cyl = 0;

    xgc.cs_mat.oc.label_cyl_from = gtk_label_new("Start vertex (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_cyl_from), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_cyl), xgc.cs_mat.oc.label_cyl_from, 0, 3, row_mat_cyl, row_mat_cyl + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_mat_cyl++;

    xgc.cs_mat.oc.cyl_i0 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.cyl_i0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.cyl_i0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.cyl_i0_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.cyl_i0_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_cyl), xgc.cs_mat.oc.cyl_i0_spin, 0, 1, row_mat_cyl, row_mat_cyl + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.cyl_i0_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.cyl_j0 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.cyl_j0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.cyl_j0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.cyl_j0_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.cyl_j0_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_cyl), xgc.cs_mat.oc.cyl_j0_spin, 1, 2, row_mat_cyl, row_mat_cyl + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.cyl_j0_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.cyl_k0 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.cyl_k0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.cyl_k0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.cyl_k0_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.cyl_k0_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_cyl), xgc.cs_mat.oc.cyl_k0_spin, 2, 3, row_mat_cyl, row_mat_cyl + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.cyl_k0_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_cyl++;

    xgc.cs_mat.oc.label_cyl_to = gtk_label_new("End vertex (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_cyl_to), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_cyl), xgc.cs_mat.oc.label_cyl_to, 0, 3, row_mat_cyl, row_mat_cyl + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_mat_cyl++;

    xgc.cs_mat.oc.cyl_in = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.cyl_in_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.cyl_in), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.cyl_in_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.cyl_in_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_cyl), xgc.cs_mat.oc.cyl_in_spin, 0, 1, row_mat_cyl, row_mat_cyl + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.cyl_in_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.cyl_jn = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.cyl_jn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.cyl_jn), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.cyl_jn_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.cyl_jn_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_cyl), xgc.cs_mat.oc.cyl_jn_spin, 1, 2, row_mat_cyl, row_mat_cyl + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.cyl_jn_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.cyl_kn = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.cyl_kn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.cyl_kn), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.cyl_kn_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.cyl_kn_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_cyl), xgc.cs_mat.oc.cyl_kn_spin, 2, 3, row_mat_cyl, row_mat_cyl + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.cyl_kn_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_cyl++;

    xgc.cs_mat.oc.label_cyl_radius = gtk_label_new("Radius [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_cyl_radius), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_cyl), xgc.cs_mat.oc.label_cyl_radius, 0, 2, row_mat_cyl, row_mat_cyl + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs_mat.oc.cyl_radius = gtk_adjustment_new(MAT_SPHERE_RADIUS, 1, 1000, 1, 10, 0);
    xgc.cs_mat.oc.cyl_radius_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.cyl_radius), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.cyl_radius_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.cyl_radius_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_cyl), xgc.cs_mat.oc.cyl_radius_spin, 2, 3, row_mat_cyl, row_mat_cyl + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.cyl_radius_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_cyl++;

    /* I. 4. Cone properties */

    create_vbox_group(&xgc.cs_mat.oc.vbg_outer_cone, &table_mat_cone, vbox_page, "Object properties:", 2, 3);

    row_mat_cone = 0;

    xgc.cs_mat.oc.label_cone_from = gtk_label_new("Start vertex (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_cone_from), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_cone), xgc.cs_mat.oc.label_cone_from, 0, 3, row_mat_cone, row_mat_cone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_mat_cone++;

    xgc.cs_mat.oc.cone_i0 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.cone_i0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.cone_i0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.cone_i0_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.cone_i0_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_cone), xgc.cs_mat.oc.cone_i0_spin, 0, 1, row_mat_cone, row_mat_cone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.cone_i0_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.cone_j0 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.cone_j0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.cone_j0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.cone_j0_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.cone_j0_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_cone), xgc.cs_mat.oc.cone_j0_spin, 1, 2, row_mat_cone, row_mat_cone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.cone_j0_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.cone_k0 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.cone_k0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.cone_k0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.cone_k0_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.cone_k0_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_cone), xgc.cs_mat.oc.cone_k0_spin, 2, 3, row_mat_cone, row_mat_cone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.cone_k0_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_cone++;

    xgc.cs_mat.oc.label_cone_to = gtk_label_new("End vertex (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_cone_to), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_cone), xgc.cs_mat.oc.label_cone_to, 0, 3, row_mat_cone, row_mat_cone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_mat_cone++;

    xgc.cs_mat.oc.cone_in = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.cone_in_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.cone_in), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.cone_in_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.cone_in_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_cone), xgc.cs_mat.oc.cone_in_spin, 0, 1, row_mat_cone, row_mat_cone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.cone_in_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.cone_jn = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.cone_jn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.cone_jn), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.cone_jn_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.cone_jn_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_cone), xgc.cs_mat.oc.cone_jn_spin, 1, 2, row_mat_cone, row_mat_cone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.cone_jn_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.cone_kn = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.cone_kn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.cone_kn), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.cone_kn_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.cone_kn_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_cone), xgc.cs_mat.oc.cone_kn_spin, 2, 3, row_mat_cone, row_mat_cone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.cone_kn_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_cone++;

    xgc.cs_mat.oc.label_cone_radius = gtk_label_new("Radius [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_cone_radius), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_cone), xgc.cs_mat.oc.label_cone_radius, 0, 2, row_mat_cone, row_mat_cone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs_mat.oc.cone_radius = gtk_adjustment_new(MAT_SPHERE_RADIUS, 1, 1000, 1, 10, 0);
    xgc.cs_mat.oc.cone_radius_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.cone_radius), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.cone_radius_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.cone_radius_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_cone), xgc.cs_mat.oc.cone_radius_spin, 2, 3, row_mat_cone, row_mat_cone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.cone_radius_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_cone++;

    /* I. 5. Cutcone properties */

    create_vbox_group(&xgc.cs_mat.oc.vbg_outer_rcone, &table_mat_rcone, vbox_page, "Object properties:", 2, 3);

    row_mat_rcone = 0;

    xgc.cs_mat.oc.label_rcone_from = gtk_label_new("Start vertex (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_rcone_from), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_rcone), xgc.cs_mat.oc.label_rcone_from, 0, 3, row_mat_rcone, row_mat_rcone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_mat_rcone++;

    xgc.cs_mat.oc.rcone_i0 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.rcone_i0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.rcone_i0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.rcone_i0_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.rcone_i0_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_rcone), xgc.cs_mat.oc.rcone_i0_spin, 0, 1, row_mat_rcone, row_mat_rcone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.rcone_i0_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.rcone_j0 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.rcone_j0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.rcone_j0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.rcone_j0_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.rcone_j0_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_rcone), xgc.cs_mat.oc.rcone_j0_spin, 1, 2, row_mat_rcone, row_mat_rcone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.rcone_j0_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.rcone_k0 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.rcone_k0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.rcone_k0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.rcone_k0_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.rcone_k0_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_rcone), xgc.cs_mat.oc.rcone_k0_spin, 2, 3, row_mat_rcone, row_mat_rcone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.rcone_k0_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_rcone++;

    xgc.cs_mat.oc.label_rcone_to = gtk_label_new("End vertex (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_rcone_to), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_rcone), xgc.cs_mat.oc.label_rcone_to, 0, 3, row_mat_rcone, row_mat_rcone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_mat_rcone++;

    xgc.cs_mat.oc.rcone_in = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.rcone_in_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.rcone_in), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.rcone_in_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.rcone_in_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_rcone), xgc.cs_mat.oc.rcone_in_spin, 0, 1, row_mat_rcone, row_mat_rcone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.rcone_in_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.rcone_jn = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.rcone_jn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.rcone_jn), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.rcone_jn_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.rcone_jn_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_rcone), xgc.cs_mat.oc.rcone_jn_spin, 1, 2, row_mat_rcone, row_mat_rcone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.rcone_jn_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.rcone_kn = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.rcone_kn_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.rcone_kn), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.rcone_kn_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.rcone_kn_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_rcone), xgc.cs_mat.oc.rcone_kn_spin, 2, 3, row_mat_rcone, row_mat_rcone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.rcone_kn_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_rcone++;

    xgc.cs_mat.oc.label_rcone_radius1 = gtk_label_new("Start radius [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_rcone_radius1), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_rcone), xgc.cs_mat.oc.label_rcone_radius1, 0, 2, row_mat_rcone, row_mat_rcone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs_mat.oc.rcone_radius1 = gtk_adjustment_new(MAT_SPHERE_RADIUS, 1, 1000, 1, 10, 0);
    xgc.cs_mat.oc.rcone_radius1_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.rcone_radius1), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.rcone_radius1_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.rcone_radius1_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_rcone), xgc.cs_mat.oc.rcone_radius1_spin, 2, 3, row_mat_rcone, row_mat_rcone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.rcone_radius1_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_rcone++;

    xgc.cs_mat.oc.label_rcone_radius2 = gtk_label_new("End radius [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_rcone_radius2), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_rcone), xgc.cs_mat.oc.label_rcone_radius2, 0, 2, row_mat_rcone, row_mat_rcone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    xgc.cs_mat.oc.rcone_radius2 = gtk_adjustment_new(MAT_SPHERE_RADIUS, 1, 1000, 1, 10, 0);
    xgc.cs_mat.oc.rcone_radius2_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.rcone_radius2), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.rcone_radius2_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.rcone_radius2_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_rcone), xgc.cs_mat.oc.rcone_radius2_spin, 2, 3, row_mat_rcone, row_mat_rcone + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.rcone_radius2_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_rcone++;


    /* I. 6. Tetrahedron properties */

    create_vbox_group(&xgc.cs_mat.oc.vbg_outer_tetrahedron, &table_mat_tthn, vbox_page, "Object properties:", 2, 3);

    row_mat_tthn = 0;

    xgc.cs_mat.oc.label_tthn_0 = gtk_label_new("Vertex (i, j, k) [voxel]:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_tthn_0), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_tthn), xgc.cs_mat.oc.label_tthn_0, 0, 3, row_mat_tthn, row_mat_tthn + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);

    row_mat_tthn++;

    xgc.cs_mat.oc.tthn_i0 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.tthn_i0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tthn_i0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tthn_i0_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tthn_i0_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_tthn), xgc.cs_mat.oc.tthn_i0_spin, 0, 1, row_mat_tthn, row_mat_tthn + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tthn_i0_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.tthn_j0 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.tthn_j0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tthn_j0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tthn_j0_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tthn_j0_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_tthn), xgc.cs_mat.oc.tthn_j0_spin, 1, 2, row_mat_tthn, row_mat_tthn + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tthn_j0_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.tthn_k0 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.tthn_k0_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tthn_k0), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tthn_k0_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tthn_k0_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_tthn), xgc.cs_mat.oc.tthn_k0_spin, 2, 3, row_mat_tthn, row_mat_tthn + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tthn_k0_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_tthn++;

    xgc.cs_mat.oc.tthn_i1 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.tthn_i1_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tthn_i1), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tthn_i1_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tthn_i1_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_tthn), xgc.cs_mat.oc.tthn_i1_spin, 0, 1, row_mat_tthn, row_mat_tthn + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tthn_i1_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.tthn_j1 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.tthn_j1_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tthn_j1), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tthn_j1_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tthn_j1_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_tthn), xgc.cs_mat.oc.tthn_j1_spin, 1, 2, row_mat_tthn, row_mat_tthn + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tthn_j1_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.tthn_k1 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.tthn_k1_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tthn_k1), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tthn_k1_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tthn_k1_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_tthn), xgc.cs_mat.oc.tthn_k1_spin, 2, 3, row_mat_tthn, row_mat_tthn + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tthn_k1_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_tthn++;

    xgc.cs_mat.oc.tthn_i2 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.tthn_i2_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tthn_i2), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tthn_i2_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tthn_i2_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_tthn), xgc.cs_mat.oc.tthn_i2_spin, 0, 1, row_mat_tthn, row_mat_tthn + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tthn_i2_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.tthn_j2 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.tthn_j2_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tthn_j2), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tthn_j2_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tthn_j2_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_tthn), xgc.cs_mat.oc.tthn_j2_spin, 1, 2, row_mat_tthn, row_mat_tthn + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tthn_j2_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.tthn_k2 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.tthn_k2_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tthn_k2), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tthn_k2_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tthn_k2_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_tthn), xgc.cs_mat.oc.tthn_k2_spin, 2, 3, row_mat_tthn, row_mat_tthn + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tthn_k2_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_tthn++;

    xgc.cs_mat.oc.tthn_i3 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.tthn_i3_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tthn_i3), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tthn_i3_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tthn_i3_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_tthn), xgc.cs_mat.oc.tthn_i3_spin, 0, 1, row_mat_tthn, row_mat_tthn + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tthn_i3_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.tthn_j3 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.tthn_j3_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tthn_j3), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tthn_j3_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tthn_j3_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_tthn), xgc.cs_mat.oc.tthn_j3_spin, 1, 2, row_mat_tthn, row_mat_tthn + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tthn_j3_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.tthn_k3 = gtk_adjustment_new(MAT_SPHERE_CENTER_I, -10000, 10000, 1, 10, 0);
    xgc.cs_mat.oc.tthn_k3_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tthn_k3), 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tthn_k3_spin), 7);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tthn_k3_spin), 1);
    gtk_table_attach(GTK_TABLE(table_mat_tthn), xgc.cs_mat.oc.tthn_k3_spin, 2, 3, row_mat_tthn, row_mat_tthn + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tthn_k3_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_tthn++;


    /* I. A. Material properties group */

    create_vbox_group(&xgc.cs_mat.oc.vbg_outer_prop, &table_mat_prop, vbox_page, "Material properties:", 14, 4);

    row_mat_prop = 0;

    xgc.cs_mat.oc.mm0 = gtk_radio_button_new_with_label(NULL, "voxel-by-voxel material properties (ε, μ, σ, σ*):");
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.mm0, 0, 4, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(GTK_TOGGLE_BUTTON(xgc.cs_mat.oc.mm0), "toggled", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_prop++;

    xgc.cs_mat.oc.epsilon = gtk_adjustment_new(MAT_PROP_EPS, 1, 100, 1, 10, 0);
    xgc.cs_mat.oc.epsilon_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.epsilon), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.epsilon_spin), 5);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.epsilon_spin), 2);
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.epsilon_spin, 0, 1, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.epsilon_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.mu = gtk_adjustment_new(MAT_PROP_MU, 1, 100, 1, 10, 0);
    xgc.cs_mat.oc.mu_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.mu), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.mu_spin), 5);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.mu_spin), 2);
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.mu_spin, 1, 2, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.mu_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.sigma = gtk_adjustment_new(MAT_PROP_SIGMA, 0, 9999999, 1, 10, 0);
    xgc.cs_mat.oc.sigma_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.sigma), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.sigma_spin), 9);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.sigma_spin), 2);
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.sigma_spin, 2, 3, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.sigma_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.sigast = gtk_adjustment_new(MAT_PROP_SIGAST, 0, 9999999, 1, 10, 0);
    xgc.cs_mat.oc.sigast_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.sigast), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.sigast_spin), 9);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.sigast_spin), 2);
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.sigast_spin, 3, 4, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.sigast_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_prop++;

    xgc.cs_mat.oc.mm1 = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(xgc.cs_mat.oc.mm0), "tabulated material properties (ε, μ, σ, σ*):");
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.mm1, 0, 4, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    //g_signal_connect_swapped(GTK_TOGGLE_BUTTON(xgc.cs_mat.oc.mm1), "toggled", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_prop++;

    xgc.cs_mat.oc.tepsilon = gtk_adjustment_new(MAT_PROP_EPS, 1, 100, 1, 10, 0);
    xgc.cs_mat.oc.tepsilon_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tepsilon), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tepsilon_spin), 5);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tepsilon_spin), 2);
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.tepsilon_spin, 0, 1, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tepsilon_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.tmu = gtk_adjustment_new(MAT_PROP_MU, 1, 100, 1, 10, 0);
    xgc.cs_mat.oc.tmu_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tmu), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tmu_spin), 5);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tmu_spin), 2);
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.tmu_spin, 1, 2, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tmu_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.tsigma = gtk_adjustment_new(MAT_PROP_SIGMA, 0, 9999999, 1, 10, 0);
    xgc.cs_mat.oc.tsigma_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tsigma), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tsigma_spin), 9);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tsigma_spin), 2);
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.tsigma_spin, 2, 3, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tsigma_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.tsigast = gtk_adjustment_new(MAT_PROP_SIGAST, 0, 9999999, 1, 10, 0);
    xgc.cs_mat.oc.tsigast_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.tsigast), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.tsigast_spin), 9);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.tsigast_spin), 2);
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.tsigast_spin, 3, 4, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.tsigast_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_prop++;

    xgc.cs_mat.oc.mm10 = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(xgc.cs_mat.oc.mm0), "perfect electric conductor");
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.mm10, 0, 4, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(GTK_TOGGLE_BUTTON(xgc.cs_mat.oc.mm10), "toggled", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_prop++;

    xgc.cs_mat.oc.mm4 = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(xgc.cs_mat.oc.mm0), "CP model based");
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.mm4, 0, 4, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    //g_signal_connect_swapped(GTK_TOGGLE_BUTTON(xgc.cs_mat.oc.mm4), "toggled", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_prop++;

    xgc.cs_mat.oc.label_eps_omega_nu = gtk_label_new(NULL);
    gtk_label_set_markup(GTK_LABEL(xgc.cs_mat.oc.label_eps_omega_nu), "ε, ω<sub>p</sub>, ν:");
    gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_eps_omega_nu), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.label_eps_omega_nu, 0, 1, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    
    xgc.cs_mat.oc.cepsilon = gtk_adjustment_new(MAT_PROP_CEPS, 0, 1e10, 1, 10, 0);
    xgc.cs_mat.oc.cepsilon_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.cepsilon), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.cepsilon_spin), 9);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.cepsilon_spin), 2);
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.cepsilon_spin, 1, 2, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.cepsilon_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.omegap = gtk_adjustment_new(MAT_PROP_OMEGA, 0, 1e10, 1, 10, 0);
    xgc.cs_mat.oc.omegap_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.omegap), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.omegap_spin), 9);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.omegap_spin), 2);
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.omegap_spin, 2, 3, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.omegap_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    xgc.cs_mat.oc.nu = gtk_adjustment_new(MAT_PROP_NU, 0, 1e10, 1, 10, 0);
    xgc.cs_mat.oc.nu_spin = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.nu), 1, 2);
    gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.nu_spin), 9);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.nu_spin), 2);
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.nu_spin, 3, 4, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
    g_signal_connect_swapped(xgc.cs_mat.oc.nu_spin, "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_prop++;

    for (i = 0; i < 2; i++) {
        xgc.cs_mat.oc.label_peak[i] = gtk_label_new("Peak parameters: A, φ, ω, γ:");
        gtk_misc_set_alignment(GTK_MISC(xgc.cs_mat.oc.label_peak[i]), 0.0, 0.5);
        gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.label_peak[i], 0, 4, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);        

        row_mat_prop++;

        xgc.cs_mat.oc.a[i] = gtk_adjustment_new(MAT_PROP_A, 2, 1e20, 1, 10, 0);
        xgc.cs_mat.oc.a_spin[i] = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.a[i]), 1, 2);
        gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.a_spin[i]), 9);
        gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.a_spin[i]), 2);
        gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.a_spin[i], 0, 1, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
        g_signal_connect_swapped(xgc.cs_mat.oc.a_spin[i], "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

        //xgc.cs_mat.oc.a[i] = gtk_adjustment_new(mat.cp3_a[i], 1, 1e20, 1, 10, 0);

        xgc.cs_mat.oc.ia_phi[i] = gtk_adjustment_new(MAT_PROP_PHI, 2, 1e20, 1, 10, 0);
        xgc.cs_mat.oc.ia_phi_spin[i] = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.ia_phi[i]), 1, 2);
        gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.ia_phi_spin[i]), 9);
        gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.ia_phi_spin[i]), 2);
        gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.ia_phi_spin[i], 1, 2, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
        g_signal_connect_swapped(xgc.cs_mat.oc.ia_phi_spin[i], "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

        //xgc.cs_mat.oc.ia_phi[i] = gtk_adjustment_new(mat.cp3_phi[i], 1, 1e20, 1, 10, 0);

        xgc.cs_mat.oc.omega[i] = gtk_adjustment_new(MAT_PROP_OMEGA, 2, 1e20, 1, 10, 0);
        xgc.cs_mat.oc.omega_spin[i] = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.omega[i]), 1, 2);
        gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.omega_spin[i]), 9);
        gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.omega_spin[i]), 2);
        gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.omega_spin[i], 2, 3, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
        g_signal_connect_swapped(xgc.cs_mat.oc.omega_spin[i], "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

        //xgc.cs_mat.oc.omega[i] = gtk_adjustment_new(mat.cp3_omega[i], 0, 1e20, 1, 10, 0);

        xgc.cs_mat.oc.gamma[i] = gtk_adjustment_new(MAT_PROP_GAMMA, 2, 1e20, 1, 10, 0);
        xgc.cs_mat.oc.gamma_spin[i] = gtk_spin_button_new(GTK_ADJUSTMENT(xgc.cs_mat.oc.gamma[i]), 1, 2);
        gtk_entry_set_width_chars(GTK_ENTRY(xgc.cs_mat.oc.gamma_spin[i]), 9);
        gtk_spin_button_set_digits(GTK_SPIN_BUTTON(xgc.cs_mat.oc.gamma_spin[i]), 2);
        gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.gamma_spin[i], 0, 1, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0, 0);
        g_signal_connect_swapped(xgc.cs_mat.oc.gamma_spin[i], "value-changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

        //xgc.cs_mat.oc.gamma[i] = gtk_adjustment_new(mat.cp3_gamma[i], 0, 1e20, 1, 10, 0);

        row_mat_prop++;
    }

    label = gtk_label_new("CP model algorithm:");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table_mat_prop), label, 0, 2, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    xgc.cs_mat.oc.metal_model = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs_mat.oc.metal_model), "RC");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs_mat.oc.metal_model), "ADE");
    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc.cs_mat.oc.metal_model), "PLRC");

    /*if (metalmodelval >= 4 && metalmodelval <= 6)
        gtk_combo_box_set_active(GTK_COMBO_BOX(xo.metalmodel), metalmodelval - 4);
    else
        gtk_combo_box_set_active(GTK_COMBO_BOX(xo.metalmodel), 1);*/
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.metal_model, 2, 4, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(G_OBJECT(xgc.cs_mat.oc.metal_model), "changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    row_mat_prop++;

    xgc.cs_mat.oc.mm99 = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(xgc.cs_mat.oc.mm0), "Loaded from database:");
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.mm99, 0, 2, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    //g_signal_connect_swapped(GTK_TOGGLE_BUTTON(xgc.cs_mat.oc.mm99), "toggled", G_CALLBACK(mat_controls_changed_cb), &xgc);

    //xgc.cs_mat.oc.overriden = gtk_entry_new();
    //gtk_entry_set_text(GTK_ENTRY(xgc.cs_mat.oc.overriden), MAT_PROP_MATERIAL_NONE);
    //gtk_entry_set_editable(GTK_ENTRY(xgc.cs_mat.oc.overriden), TRUE);
    //gtk_widget_set_can_focus(GTK_WIDGET(xgc.cs_mat.oc.overriden), TRUE);
    //gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.overriden, 2, 4, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    //g_signal_connect_swapped(G_OBJECT(xgc.cs_mat.oc.overriden), "changed", G_CALLBACK(mat_controls_changed_cb), &xgc);

    //xgc.cs_mat.oc.overriden = gtk_combo_box_new_text();
    xgc.cs_mat.oc.overriden = gtk_combo_box_text_new();
    xsv_fill_material_combo(&xgc.cs_mat.oc.overriden);    
    g_signal_connect_swapped(G_OBJECT(xgc.cs_mat.oc.overriden), "changed", G_CALLBACK(mat_controls_changed_cb), &xgc);
    gtk_table_attach(GTK_TABLE(table_mat_prop), xgc.cs_mat.oc.overriden, 2, 4, row_mat_prop, row_mat_prop + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row_mat_prop++;

    page = gtk_vbox_new(FALSE, 0);
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, gtk_label_new(("Parameter file")));    

    scroll_window_par_file = gtk_scrolled_window_new(NULL, NULL);
    gtk_box_pack_start(GTK_BOX(page), scroll_window_par_file, TRUE, TRUE, 4);

    xgc.tb_par = gtk_text_buffer_new(NULL);
    view_par_file = gtk_text_view_new_with_buffer(xgc.tb_par);
    gtk_container_add(GTK_CONTAINER(scroll_window_par_file), view_par_file);

    /* material view notebook */
    vbox_mat = gtk_vbox_new(FALSE, 0);
    //gtk_box_pack_start(GTK_BOX(vbox2_left), vbox_mat, TRUE, TRUE, 5);


    xgc.mat_notebook = notebook = gtk_notebook_new();
    gtk_box_pack_start(GTK_BOX(vbox_mat), notebook, TRUE, TRUE, 0);    

    /* parameter tree page */
    //page = gtk_vbox_new(FALSE, 0);
    xgc.page_mat_tree = page_mat_tree = gtk_hbox_new(FALSE, 0);
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page_mat_tree, gtk_label_new(("Vector objects")));
    
    xgc.ts_mat = gtk_tree_store_new(4, G_TYPE_STRING, G_TYPE_BOOLEAN, G_TYPE_INT, G_TYPE_BOOLEAN);

    xgc.tv_mat = gtk_tree_view_new_with_model(GTK_TREE_MODEL(xgc.ts_mat));

    g_signal_connect(xgc.tv_mat, "row-activated", G_CALLBACK(mat_row_activated_cb), NULL);
    g_signal_connect(xgc.tv_mat, "button-release-event", G_CALLBACK(mat_row_button_release_cb), &xgc);
    g_signal_connect(xgc.tv_mat, "key_press_event", G_CALLBACK (on_key_press_mat_tree_cb), &xgc);

    col_offset = gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(xgc.tv_mat), 
                                                             -1, "object", 
                                                             gtk_cell_renderer_text_new(), 
                                                             "text", COLUMN_PARAMETER, 
                                                             NULL);

    column = gtk_tree_view_get_column(GTK_TREE_VIEW(xgc.tv_mat), col_offset - 1);
    gtk_tree_view_column_set_expand(GTK_TREE_VIEW_COLUMN(column), TRUE);

    renderer = gtk_cell_renderer_toggle_new();
    g_signal_connect(renderer, "toggled", G_CALLBACK(mat_toggled_cb), &xgc);

    col_offset = gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(xgc.tv_mat),
                                                             -1, "Visible",
                                                             renderer,
                                                             "active", COLUMN_CHECK,
                                                             "visible", COLUMN_SHOW_TOGGLE,
                                                             //"indicator-size", 3,
                                                             NULL);

    column = gtk_tree_view_get_column(GTK_TREE_VIEW(xgc.tv_mat), col_offset - 1);
    gtk_tree_view_column_set_sizing(GTK_TREE_VIEW_COLUMN(column), GTK_TREE_VIEW_COLUMN_FIXED);
    gtk_tree_view_column_set_fixed_width(GTK_TREE_VIEW_COLUMN(column), 50);
    gtk_tree_view_column_set_alignment(GTK_TREE_VIEW_COLUMN(column), 1.0);
    gtk_tree_view_column_set_clickable(GTK_TREE_VIEW_COLUMN(column), TRUE);

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(xgc.tv_mat));
    g_signal_connect(selection, "changed", G_CALLBACK(mat_row_selected_cb), &xgc);


    scroll_window_mat_tree = gtk_scrolled_window_new(NULL, NULL);
    //gtk_box_pack_start(GTK_BOX(page), scroll_window_mat_view, TRUE, TRUE, 4);
    gtk_container_add(GTK_CONTAINER(scroll_window_mat_tree), xgc.tv_mat);



    page = gtk_vbox_new(FALSE, 0);
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, gtk_label_new(("Vector material file")));

    scroll_window_mat_file = gtk_scrolled_window_new(NULL, NULL);
    gtk_box_pack_start(GTK_BOX(page), scroll_window_mat_file, TRUE, TRUE, 4);

    xgc.tb_mat = gtk_text_buffer_new(NULL);
    view_mat_file = gtk_text_view_new_with_buffer(xgc.tb_mat);
    gtk_container_add(GTK_CONTAINER(scroll_window_mat_file), view_mat_file);

    /* right side */
    vbox2_right = gtk_vbox_new(FALSE, 0);
    //gtk_widget_set_size_request(vbox2_right, 300, 300);
    //gtk_box_pack_start(GTK_BOX(hbox), vbox2_right, TRUE, TRUE, 5);


    /* 3d view */
    view_scene = gtk_drawing_area_new();
    gtk_widget_set_size_request(view_scene, 50, 50);
    //gtk_box_pack_start(GTK_BOX(vbox2_right), view_scene, TRUE, TRUE, 5);

    /*setup gl*/
    if (scene_setup(view_scene) == FALSE)
        return 0;    

    g_signal_connect(view_scene, "configure-event", G_CALLBACK(scene_configure_cb), &xgc);
    g_signal_connect(view_scene, "expose-event", G_CALLBACK(scene_expose_cb), &xgc);

    gtk_widget_add_events(view_scene, GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK | GDK_POINTER_MOTION_MASK);

    g_signal_connect(G_OBJECT(view_scene), "motion-notify-event", G_CALLBACK(scene_motion_notify_cb), &xgc);
    g_signal_connect(G_OBJECT(view_scene), "button-press-event", G_CALLBACK(scene_button_press_cb), &xgc);
    g_signal_connect(G_OBJECT(view_scene), "button-release-event", G_CALLBACK(scene_button_release_cb), &xgc);
    g_signal_connect(G_OBJECT(view_scene), "scroll-event", G_CALLBACK(scene_scroll_cb), &xgc);    

    xgc.view_scene = view_scene;
    scene_reset_cb(xgc.toplevel, &xgc);

    /* 3d scene view output options + Graph view + Output view */
    vbox2bottom_right = gtk_vbox_new(FALSE, 0);
    //gtk_box_pack_start(GTK_BOX(vbox2_right), vbox2bottom_right, TRUE, TRUE, 0);

    /* 3d scene view output options */

    hbox2 = gtk_hbox_new(FALSE, 0);
    //gtk_box_pack_start(GTK_BOX(vbox2_right), hbox2, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox2bottom_right), hbox2, FALSE, FALSE, 5);


    xgc.image_check = gtk_check_button_new_with_mnemonic("Show cross-section (image):");

    gtk_box_pack_start(GTK_BOX(hbox2), xgc.image_check, FALSE, FALSE, 5);
    g_signal_connect(G_OBJECT(xgc.image_check), "toggled", G_CALLBACK(image_changed_cb), &xgc);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.image_check), TRUE);

    //xgc.image_combo = gtk_combo_box_new();
    xgc.image_combo = gtk_combo_box_text_new();
    g_signal_connect(G_OBJECT(xgc.image_combo), "changed", G_CALLBACK(image_changed_cb), &xgc);
    gtk_box_pack_start(GTK_BOX(hbox2), xgc.image_combo, TRUE, TRUE, 5);


    xgc.image_squarecheck = gtk_check_button_new_with_mnemonic("square values");
    g_signal_connect(G_OBJECT(xgc.image_squarecheck), "toggled", G_CALLBACK(image_changed_cb), &xgc);

    gtk_box_pack_start(GTK_BOX(hbox2), xgc.image_squarecheck, TRUE, TRUE, 5);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.image_squarecheck), TRUE);


    /* row of graphs with selection */

    hbox2 = gtk_hbox_new(FALSE, 0);
    //gtk_box_pack_start(GTK_BOX(vbox2_right), hbox2, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(vbox2bottom_right), hbox2, TRUE, TRUE, 5);

    vbox2_graph = gtk_vbox_new(FALSE, 0);
    //gtk_box_pack_start(GTK_BOX(hbox2), vbox2_graph, TRUE, TRUE, 0);

    vbox3 = gtk_vbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox2_graph), vbox3, FALSE, FALSE, 0);

    xgc.graph_check = gtk_check_button_new_with_mnemonic("Show graph:");

    gtk_box_pack_start(GTK_BOX(vbox3), xgc.graph_check, FALSE, FALSE, 5);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc.graph_check), TRUE);


    xgc.graph_combo = gtk_combo_box_new_text();
    g_signal_connect(G_OBJECT(xgc.graph_combo), "changed", G_CALLBACK(graph_changed_cb), &xgc);
    gtk_box_pack_start(GTK_BOX(vbox3), xgc.graph_combo, TRUE, TRUE, 5);


    xgc.gmodel1 = gwy_graph_model_new();
    xgc.gcmodel1 = gwy_graph_curve_model_new();
    gwy_graph_model_add_curve(xgc.gmodel1, xgc.gcmodel1);
    graph = gwy_graph_new(xgc.gmodel1);
    gtk_box_pack_start(GTK_BOX(vbox2_graph), graph, TRUE, TRUE, 5);

    frame = gtk_frame_new("GSvit output");
    gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_IN);
    //gtk_box_pack_start(GTK_BOX(hbox2), frame, TRUE, TRUE, 4);

    scroll_window_output = gtk_scrolled_window_new(NULL, NULL);    
    gtk_container_add(GTK_CONTAINER(frame), scroll_window_output);

    xgc.tb_out = gtk_text_buffer_new(NULL);
    //g_signal_connect(xgc.tb_mat, "changed", G_CALLBACK(mat_buffer_changed_cb), &xgc);
    xgc.tv_out = gtk_text_view_new_with_buffer(xgc.tb_out);
    gtk_container_add(GTK_CONTAINER(scroll_window_output), xgc.tv_out);


    //////////////////////////////////////////////////////////////////////////
    // Paned horizontal
    // left : Parameter view, Material object view
    // right: 3D scene view, 3D scene view options, Graph, Output
       
    hpaned = gtk_hpaned_new();
    gtk_box_pack_start(GTK_BOX(hbox), hpaned, TRUE, TRUE, 5);
    gtk_paned_set_position(GTK_PANED(hpaned), 800);

    gtk_paned_pack1(GTK_PANED(hpaned), vbox2_left, TRUE, FALSE);
    gtk_paned_pack2(GTK_PANED(hpaned), vbox2_right, TRUE, FALSE);
    //////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////
    // Paned vertical
    // top : 3D scene view
    // bottom: 3D scene view options, Graph, Output

    vpaned = gtk_vpaned_new();
    gtk_box_pack_start(GTK_BOX(vbox2_right), vpaned, TRUE, TRUE, 5);
    gtk_paned_set_position(GTK_PANED(vpaned), 300);

    gtk_paned_pack1(GTK_PANED(vpaned), view_scene, TRUE, FALSE);
    gtk_paned_pack2(GTK_PANED(vpaned), vbox2bottom_right, TRUE, FALSE);
    //////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////
    // Paned horizontal
    // left : Parameter tree view
    // right: Parameter values view

    hpaned = gtk_hpaned_new();
    gtk_box_pack_start(GTK_BOX(page_par_tree), hpaned, TRUE, TRUE, 5);
    gtk_paned_set_position(GTK_PANED(hpaned), 260);

    gtk_paned_pack1(GTK_PANED(hpaned), scroll_window_par_tree, TRUE, FALSE);
    gtk_paned_pack2(GTK_PANED(hpaned), scroll_window_par_values, TRUE, FALSE);
    //////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////
    // Paned horizontal
    // left : Material tree view
    // right: Material values view

    hpaned = gtk_hpaned_new();
    gtk_box_pack_start(GTK_BOX(page_mat_tree), hpaned, TRUE, TRUE, 5);
    gtk_paned_set_position(GTK_PANED(hpaned), 260);

    gtk_paned_pack1(GTK_PANED(hpaned), scroll_window_mat_tree, TRUE, FALSE);
    gtk_paned_pack2(GTK_PANED(hpaned), scroll_window_mat_values, TRUE, FALSE);
    //////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////
    // Paned horizontal
    // left : Graph
    // right: Output

    hpaned = gtk_hpaned_new();
    gtk_box_pack_start(GTK_BOX(hbox2), hpaned, TRUE, TRUE, 5);

    gtk_paned_pack1(GTK_PANED(hpaned), vbox2_graph, TRUE, FALSE);
    gtk_paned_pack2(GTK_PANED(hpaned), frame, TRUE, FALSE);
    //////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////
    // Paned vertical
    // top : Parameter view
    // bottom: Material object view

    vpaned = gtk_vpaned_new();
    //gtk_paned_set_position(GTK_PANED(vpaned), 0);
    gtk_box_pack_start(GTK_BOX(vbox2_left), vpaned, TRUE, TRUE, 5);
    gtk_paned_set_position(GTK_PANED(vpaned), 300);

    gtk_paned_pack1(GTK_PANED(vpaned), vbox_par, TRUE, FALSE);
    gtk_paned_pack2(GTK_PANED(vpaned), vbox_mat, TRUE, FALSE);

    //gint iii = gtk_paned_get_position(GTK_PANED(hpaned));
    //gtk_paned_set_position(GTK_PANED(vpaned), 500);

    //////////////////////////////////////////////////////////////////////////



    /**********   status bar *********************************/

    xgc.statusbar = gtk_statusbar_new();
    xgc.statusbar_context_id = gtk_statusbar_get_context_id(GTK_STATUSBAR(xgc.statusbar), "Statusbar");
    gtk_box_pack_start(GTK_BOX(vbox), xgc.statusbar, FALSE, FALSE, 0);
    gtk_statusbar_push(GTK_STATUSBAR(xgc.statusbar), xgc.statusbar_context_id, "Ready");


    /*********************************************************/



    g_signal_connect(G_OBJECT(window), "destroy", G_CALLBACK(quit_cb), &xgc);


    gtk_widget_show_all(window);
    gtk_window_present(GTK_WINDOW(window));     // NOTE: fixes losing focus when XSvit starts on Windows

    // TODO: remember and set main window size and paned windows position
    gtk_window_maximize(GTK_WINDOW(window));


    set_all_controls_invisible(&xgc);
    set_computation_controls_sensitive(&xgc);

    set_all_mat_controls_invisible(&xgc);

    g_signal_connect(xgc.tb_par, "changed", G_CALLBACK(par_buffer_changed_cb), &xgc);
    g_signal_connect(xgc.tb_mat, "changed", G_CALLBACK(mat_buffer_changed_cb), &xgc);
    g_signal_connect_swapped(xgc.param_notebook, "switch-page", G_CALLBACK(pars_page_switched_cb), &xgc);
    g_signal_connect_swapped(xgc.mat_notebook, "switch-page", G_CALLBACK(material_page_switched_cb), &xgc);

    clear_settings(&(xgc.data.set), FALSE);
    clear_settings_mat(&(xgc.data.set_mat), FALSE);
    if (xgc.loadfile == NULL) {
        file_new_cb(NULL, &xgc);
    } else {
        xgc.data.parfilename = g_strdup(xgc.loadfile);

        file_load(&xgc, xgc.loadfile, NULL);
    }

    //g_signal_connect(xgc.tb_par, "changed", G_CALLBACK(par_buffer_changed_cb), &xgc);
    //g_signal_connect(xgc.tb_mat, "changed", G_CALLBACK(mat_buffer_changed_cb), &xgc);
    //g_signal_connect_swapped(xgc.param_notebook, "switch-page", G_CALLBACK(pars_page_switched_cb), &xgc);
    
    gtk_main();

    return 0;
}

void par_controls_changed_cb(XGControls *xgc)
{
    par_controls_changed(xgc);
}

void mat_controls_changed_cb(XGControls *xgc)
{
    mat_controls_changed(xgc);
}

void
ss_choose_source_file_cb(GtkWidget *widget, XGControls *xgc)
{
    GtkWidget *dialog;
    GtkFileFilter *filter = gtk_file_filter_new();
    gchar  *location = NULL;
    gint id = 0;

    dialog = gtk_file_chooser_dialog_new("Choose source file:",
                                         GTK_WINDOW(xgc->toplevel),
                                         GTK_FILE_CHOOSER_ACTION_OPEN,
                                         GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                         GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                                         NULL);
    gtk_file_filter_set_name(filter, "Source file (*.txt)");
    gtk_file_filter_add_pattern(filter, "*.txt");
    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filter);
    if (xgc->curdir)
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), xgc->curdir);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        location = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
        if (xgc->par_selid >= SET_PSOURCE && xgc->par_selid < SET_POUT) {
            id = xgc->par_selid - SET_PSOURCE;
            g_free(xgc->data.set.ss.pnts[id].source_filename);
            if (xgc->curdir && g_str_has_prefix(location, xgc->curdir))
                xgc->data.set.ss.pnts[id].source_filename = g_strdup(location + g_utf8_strlen(xgc->curdir, -1) + 1);
            else
                xgc->data.set.ss.pnts[id].source_filename = g_strdup(location);

            gtk_entry_set_text(GTK_ENTRY(xgc->cs.sc.source_filename), xgc->data.set.ss.pnts[id].source_filename);
        } else if (xgc->par_selid == SET_SF) {
            g_free(xgc->data.set.ss.sf.source_filename);
            if (xgc->curdir && g_str_has_prefix(location, xgc->curdir))
                xgc->data.set.ss.sf.source_filename = g_strdup(location + g_utf8_strlen(xgc->curdir, -1) + 1);
            else
                xgc->data.set.ss.sf.source_filename = g_strdup(location);

            gtk_entry_set_text(GTK_ENTRY(xgc->cs.sc.source_filename), xgc->data.set.ss.sf.source_filename);
        } else if (xgc->par_selid == SET_TSF) {
            g_free(xgc->data.set.ss.tsf.source_filename);
            if (xgc->curdir && g_str_has_prefix(location, xgc->curdir))
                xgc->data.set.ss.tsf.source_filename = g_strdup(location + g_utf8_strlen(xgc->curdir, -1) + 1);
            else
                xgc->data.set.ss.tsf.source_filename = g_strdup(location);

            gtk_entry_set_text(GTK_ENTRY(xgc->cs.sc.source_filename), xgc->data.set.ss.tsf.source_filename);
        } else if (xgc->par_selid == SET_TSFF) {
            g_free(xgc->data.set.ss.tsff.source_filename);
            if (xgc->curdir && g_str_has_prefix(location, xgc->curdir))
                xgc->data.set.ss.tsff.source_filename = g_strdup(location + g_utf8_strlen(xgc->curdir, -1) + 1);
            else
                xgc->data.set.ss.tsff.source_filename = g_strdup(location);

            gtk_entry_set_text(GTK_ENTRY(xgc->cs.sc.source_filename), xgc->data.set.ss.tsff.source_filename);
        } else if (xgc->par_selid == SET_LTSF) {
            g_free(xgc->data.set.ss.ltsf.source_filename);
            if (xgc->curdir && g_str_has_prefix(location, xgc->curdir))
                xgc->data.set.ss.ltsf.source_filename = g_strdup(location + g_utf8_strlen(xgc->curdir, -1) + 1);
            else
                xgc->data.set.ss.ltsf.source_filename = g_strdup(location);

            gtk_entry_set_text(GTK_ENTRY(xgc->cs.sc.source_filename), xgc->data.set.ss.ltsf.source_filename);
        } else if (xgc->par_selid == SET_LTSFF) {
            g_free(xgc->data.set.ss.ltsff.source_filename);
            if (xgc->curdir && g_str_has_prefix(location, xgc->curdir))
                xgc->data.set.ss.ltsff.source_filename = g_strdup(location + g_utf8_strlen(xgc->curdir, -1) + 1);
            else
                xgc->data.set.ss.ltsff.source_filename = g_strdup(location);

            gtk_entry_set_text(GTK_ENTRY(xgc->cs.sc.source_filename), xgc->data.set.ss.ltsff.source_filename);
        }
        g_free(location);
    }
    gtk_widget_destroy(dialog);

    par_controls_changed_cb(xgc);
}

void
so_choose_output_file_cb(GtkWidget *widget, XGControls *xgc)
{
    GtkWidget *dialog;
    GtkFileFilter *filter = gtk_file_filter_new();
    gchar  *filename, *basename;
    gint id;


    dialog = gtk_file_chooser_dialog_new("Choose output file:",
                                         GTK_WINDOW(xgc->toplevel),
                                         GTK_FILE_CHOOSER_ACTION_SAVE,
                                         GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                         GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                                         NULL);
    gtk_file_filter_set_name(filter, "ASCII output");
    gtk_file_filter_add_pattern(filter, "*.*");
    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filter);
    if (xgc->curdir)
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), xgc->curdir);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        //filename = g_strdup(gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog)));
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

        if (xgc->par_selid == SET_OUT) {
            /* VI. 1. General output properties*/
            g_free(xgc->data.set.so.outfile);
            if (xgc->curdir && g_str_has_prefix(filename, xgc->curdir)) {
                basename = g_path_get_basename(filename);
                xgc->data.set.so.outfile = g_strdup(basename);
                g_free(basename);
            } else
                xgc->data.set.so.outfile = g_strdup(filename);

            gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.general_output_filename), xgc->data.set.so.outfile);
            gtk_widget_set_tooltip_text(xgc->cs.oc.general_output_filename, filename);
        } else if (xgc->par_selid >= SET_POUT && xgc->par_selid < SET_IOUT) {
            /* VI. 2. Point output properties*/
            id = xgc->par_selid - SET_POUT;
            if (xgc->curdir && g_str_has_prefix(filename, xgc->curdir)) {
                basename = g_path_get_basename(filename);
                g_snprintf(xgc->data.set.so.pnts[id].filebase, sizeof(xgc->data.set.so.pnts[id].filebase), "%s", basename);
                g_free(basename);
            } else
                g_snprintf(xgc->data.set.so.pnts[id].filebase, sizeof(xgc->data.set.so.pnts[id].filebase), "%s", filename);

            gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.point_output_filename), xgc->data.set.so.pnts[id].filebase);
            gtk_widget_set_tooltip_text(xgc->cs.oc.point_output_filename, filename);
        } else if (xgc->par_selid >= SET_IIOUT && xgc->par_selid < SET_COUT) {
            /* VI. 3. Image output properties */
        } else if (xgc->par_selid >= SET_IOUT && xgc->par_selid < SET_IIOUT) {
            /* VI. 4. Plane output properties */
            id = xgc->par_selid - SET_IOUT;
            if (xgc->curdir && g_str_has_prefix(filename, xgc->curdir)) {
                basename = g_path_get_basename(filename);
                g_snprintf(xgc->data.set.so.plns[id].filebase, sizeof(xgc->data.set.so.plns[id].filebase), "%s", basename);
                g_free(basename);
            } else
                g_snprintf(xgc->data.set.so.plns[id].filebase, sizeof(xgc->data.set.so.plns[id].filebase), "%s", filename);

            gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.plane_output_filename), xgc->data.set.so.plns[id].filebase);
            gtk_widget_set_tooltip_text(xgc->cs.oc.plane_output_filename, filename);
        } else if (xgc->par_selid >= SET_COUT && xgc->par_selid < SET_SOUT) {
            /* VI. 5. Volume output properties */
            id = xgc->par_selid - SET_COUT;
            if (xgc->curdir && g_str_has_prefix(filename, xgc->curdir)) {
                basename = g_path_get_basename(filename);
                g_snprintf(xgc->data.set.so.cubs[id].filebase, sizeof(xgc->data.set.so.cubs[id].filebase), "%s", basename);
                g_free(basename);
            } else
                g_snprintf(xgc->data.set.so.cubs[id].filebase, sizeof(xgc->data.set.so.cubs[id].filebase), "%s", filename);

            gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.volume_output_filename), xgc->data.set.so.cubs[id].filebase);
            gtk_widget_set_tooltip_text(xgc->cs.oc.volume_output_filename, filename);
        } else if (xgc->par_selid >= SET_SOUT && xgc->par_selid < SET_FOUT) {
            /* VI. 5. Sum output properties */
            id = xgc->par_selid - SET_SOUT;
            if (xgc->curdir && g_str_has_prefix(filename, xgc->curdir)) {
                basename = g_path_get_basename(filename);
                g_snprintf(xgc->data.set.so.sums[id].filename, sizeof(xgc->data.set.so.sums[id].filename), "%s", basename);
                g_free(basename);
            } else
                g_snprintf(xgc->data.set.so.sums[id].filename, sizeof(xgc->data.set.so.sums[id].filename), "%s", filename);

            gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.sum_output_filename), xgc->data.set.so.sums[id].filename);
            gtk_widget_set_tooltip_text(xgc->cs.oc.sum_output_filename, filename);
        } else if (xgc->par_selid >= SET_FOUT && xgc->par_selid < SET_GROW) {
            /* VI. 6. Force output properties */
            id = xgc->par_selid - SET_FOUT;
            if (xgc->curdir && g_str_has_prefix(filename, xgc->curdir)) {
                basename = g_path_get_basename(filename);
                g_snprintf(xgc->data.set.so.forces[id].filename, sizeof(xgc->data.set.so.forces[id].filename), "%s", basename);
                g_free(basename);
            } else
                g_snprintf(xgc->data.set.so.forces[id].filename, sizeof(xgc->data.set.so.forces[id].filename), "%s", filename);

            gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.force_output_filename), xgc->data.set.so.forces[id].filename);
            gtk_widget_set_tooltip_text(xgc->cs.oc.force_output_filename, filename);
        }        

        g_free(filename);
    }
    gtk_widget_destroy(dialog);


    par_controls_changed(xgc);
    if(write_parfile(xgc)) {
        pt_create_tree(xgc);
    }
} /* so_choose_output_file_cb */

void
nfffp_choose_output_file_cb(GtkWidget *widget, XGControls *xgc)
{
    GtkWidget *dialog;
    GtkFileFilter *filter = gtk_file_filter_new();
    gchar  *filename;
    gint id;


    dialog = gtk_file_chooser_dialog_new("Choose output file:",
                                         GTK_WINDOW(xgc->toplevel),
                                         GTK_FILE_CHOOSER_ACTION_SAVE,
                                         GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                         GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                                         NULL);
    gtk_file_filter_set_name(filter, "ASCII output");
    gtk_file_filter_add_pattern(filter, "*.*");
    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filter);
    if (xgc->curdir)
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), xgc->curdir);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

        if (xgc->par_selid >= SET_NFFFP && xgc->par_selid < SET_PNFAREA) {
            /* VII. 2. Near field to far field point */

            id = xgc->par_selid - SET_NFFFP;

            if (xgc->curdir && g_str_has_prefix(filename, xgc->curdir))
                xgc->data.set.sf.source_filename[id] = g_path_get_basename(filename);
            else
                xgc->data.set.sf.source_filename[id] = g_strdup(filename);

            gtk_entry_set_text(GTK_ENTRY(xgc->cs.fc.nfffp_filename), xgc->data.set.sf.source_filename[id]);
        } else if (xgc->par_selid >= SET_PNFFFP) {
            /* VIII. 2. Periodic near field to far field point */

            id = xgc->par_selid - SET_PNFFFP;

            if (xgc->curdir && g_str_has_prefix(filename, xgc->curdir))
                xgc->data.set.spf.source_filename[id] = g_path_get_basename(filename);
            else
                xgc->data.set.spf.source_filename[id] = g_strdup(filename);

            gtk_entry_set_text(GTK_ENTRY(xgc->cs.fc.nfffp_filename), xgc->data.set.spf.source_filename[id]);
        }

        g_free(filename);
    }
    gtk_widget_destroy(dialog);


    par_controls_changed_cb(xgc);
}

void
material_choose_voxel_file_cb(GtkWidget *widget, XGControls *xgc)
{
    GtkWidget *dialog;
    gchar  *filename;
    //GtkFileFilter *filter = gtk_file_filter_new();
    //gchar  *location;

    dialog = gtk_file_chooser_dialog_new("Choose voxel-by-voxel material file:",
                                         GTK_WINDOW(xgc->toplevel),
                                         GTK_FILE_CHOOSER_ACTION_OPEN,
                                         GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                         GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                                         NULL);
    if (xgc->curdir)
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), xgc->curdir);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {        
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

        if (xgc->curdir && g_str_has_prefix(filename, xgc->curdir))
            xgc->data.set.sm.in_voxel_filename = g_path_get_basename(filename);
        else
            xgc->data.set.sm.in_voxel_filename = g_strdup(filename);

        gtk_entry_set_text(GTK_ENTRY(xgc->cs.mc.material_voxel_filename), xgc->data.set.sm.in_voxel_filename);

        g_free(filename);
    }
    gtk_widget_destroy(dialog);

    par_controls_changed_cb(xgc);
}

void
material_choose_vector_file_cb(GtkWidget *widget, XGControls *xgc)
{
    GtkWidget *dialog;
    gchar  *filename;
    //GtkFileFilter *filter = gtk_file_filter_new();
    //gchar  *location;

    dialog = gtk_file_chooser_dialog_new("Choose vector material file:",
                                         GTK_WINDOW(xgc->toplevel),
                                         GTK_FILE_CHOOSER_ACTION_OPEN,
                                         GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                         GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                                         NULL);
    if (xgc->curdir)
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), xgc->curdir);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

        if (xgc->curdir && g_str_has_prefix(filename, xgc->curdir))
            xgc->data.set.sm.in_vector_filename = g_path_get_basename(filename);
        else
            xgc->data.set.sm.in_vector_filename = g_strdup(filename);

        gtk_entry_set_text(GTK_ENTRY(xgc->cs.mc.material_vector_filename), xgc->data.set.sm.in_vector_filename);

        g_free(filename);
    }
    gtk_widget_destroy(dialog);

    par_controls_changed_cb(xgc);
}





/*
void
ss_choose_output_file_cb(GtkWidget *widget, XGOutData *xo)
{
    GtkWidget *dialog;
    GtkFileFilter *filter = gtk_file_filter_new();
    gchar *location;

    dialog = gtk_file_chooser_dialog_new ("Choose general output Gwyddion file:",
        GTK_WINDOW(xo->toplevel),
        GTK_FILE_CHOOSER_ACTION_SAVE,
        GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
        GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
        NULL);
    gtk_file_filter_set_name(filter, "Gwyddion file");
    gtk_file_filter_add_pattern(filter, "*.gwy");
    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filter);
    if (xo->curdir) gtk_file_chooser_set_current_folder (GTK_FILE_CHOOSER (dialog), xo->curdir);

    if (gtk_dialog_run (GTK_DIALOG (dialog)) == GTK_RESPONSE_ACCEPT)
    {
        location = g_strdup(gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog)));
        if (xo->curdir && g_str_has_prefix(location, xo->curdir)) xo->outloc = location + g_utf8_strlen(xo->curdir, 100) + 1;
        else xo->outloc = location;
        gtk_label_set_text(GTK_LABEL(xo->outlab), xo->outloc);
    }
    gtk_widget_destroy (dialog);

    ss_controls_changed(xgc);
}
*/

/*
void mouse(int button, int state, int x, int y)
{
if (button == GLUT_LEFT_BUTTON &amp;&amp; state == GLUT_DOWN)
{
mouseDown = true;

xdiff = x - yrot;
ydiff = -y + xrot;
}
else
mouseDown = false;
}

void mouseMotion(int x, int y)
{
if (mouseDown)
{
yrot = x - xdiff;
xrot = y + ydiff;

glutPostRedisplay();
}
}
*/


/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

