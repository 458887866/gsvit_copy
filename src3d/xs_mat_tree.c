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


/* mat_tree.c :
*  Material file group tree functions
*/

#include "xs_mat_tree.h"
#include "write.h"
#include "messages.h"


void
mt_tree_assemble_row_text(XGControls    *xgc,
                           gchar         *buff,          /* output text */
                           gint          buff_size,      /* output text buffer size */
                           gint          node_id,        /* node group id (expandable node) */
                           gint          index)       /* position in node */
{    
    if (node_id >= SET_MAT_SPHERE && node_id < SET_MAT_VOXEL) {
        SvSphere sphere;
        sphere = g_array_index(xgc->data.set_mat.spheres, SvSphere, index);
        g_snprintf(buff, buff_size, "Sphere radius %g at %g %g %g", sphere.radius, sphere.pnt1[0], sphere.pnt1[1], sphere.pnt1[2]);
    } else if (node_id >= SET_MAT_VOXEL && node_id < SET_MAT_CYLINDER) {
        SvVoxel voxel;
        voxel = g_array_index(xgc->data.set_mat.voxels, SvVoxel, index);
        g_snprintf(buff, buff_size, "Box at %g %g %g ... %g %g %g", voxel.pnt1[0], voxel.pnt1[1], voxel.pnt1[2], voxel.pnt2[0], voxel.pnt2[1], voxel.pnt2[2]);
    } else if (node_id >= SET_MAT_CYLINDER && node_id < SET_MAT_CONE) {
        SvCylinder cyl;
        cyl = g_array_index(xgc->data.set_mat.cylinders, SvCylinder, index);
        g_snprintf(buff, buff_size, "Cylinder radius %g at %g %g %g ... %g %g %g", cyl.radius, cyl.pnt1[0], cyl.pnt1[1], cyl.pnt1[2], cyl.pnt2[0], cyl.pnt2[1], cyl.pnt2[2]);
    } else if (node_id >= SET_MAT_CONE && node_id < SET_MAT_RCONE) {
        SvCone cone;
        cone = g_array_index(xgc->data.set_mat.cones, SvCone, index);
        g_snprintf(buff, buff_size, "Cone radius %g at %g %g %g ... %g %g %g", cone.radius, cone.pnt1[0], cone.pnt1[1], cone.pnt1[2], cone.pnt2[0], cone.pnt2[1], cone.pnt2[2]);
    } else if (node_id >= SET_MAT_RCONE && node_id < SET_MAT_GWYDD) {
        SvRCone rcone;
        rcone = g_array_index(xgc->data.set_mat.rcones, SvRCone, index);
        g_snprintf(buff, buff_size, "Cut cone %g ... %g at %g %g %g ... %g %g %g", rcone.radius1, rcone.radius2, rcone.pnt1[0], rcone.pnt1[1], rcone.pnt1[2], rcone.pnt2[0], rcone.pnt2[1], rcone.pnt2[2]);
    } else if (node_id >= SET_MAT_GWYDD && node_id < SET_MAT_MESH) {
        SvGwydd gwydd = g_array_index(xgc->data.set_mat.gwydds, SvGwydd, index);

        if (xgc->data.set_mat.gwydd_nvx[xgc->data.set_mat.gwydds->len - 1] == 0)
            g_snprintf(buff, buff_size, "Gwyddion field missed: pos %d %d %d span %d %d %d ... %d %d %d", gwydd.i, gwydd.j, gwydd.k,
                       (gint)xgc->data.set_mat.gwydd_xmin[node_id], (gint)xgc->data.set_mat.gwydd_ymin[node_id], (gint)xgc->data.set_mat.gwydd_zmin[node_id],
                       (gint)xgc->data.set_mat.gwydd_xmax[node_id], (gint)xgc->data.set_mat.gwydd_ymax[node_id], (gint)xgc->data.set_mat.gwydd_zmax[node_id]);
        else
            g_snprintf(buff, buff_size, "Gwyddion field at %d %d %d span %d %d %d ... %d %d %d", gwydd.i, gwydd.j, gwydd.k,
                       (gint)xgc->data.set_mat.gwydd_xmin[node_id], (gint)xgc->data.set_mat.gwydd_ymin[node_id], (gint)xgc->data.set_mat.gwydd_zmin[node_id],
                       (gint)xgc->data.set_mat.gwydd_xmax[node_id], (gint)xgc->data.set_mat.gwydd_ymax[node_id], (gint)xgc->data.set_mat.gwydd_zmax[node_id]);
    } else if (node_id >= SET_MAT_MESH && node_id < SET_MAT_TETRAHEDRON) {
        g_snprintf(buff, buff_size, "mesh (%s) at %g %g %g", xgc->data.set_mat.tetgen_filebase[node_id], xgc->data.set_mat.tetgen_xshift[node_id], xgc->data.set_mat.tetgen_yshift[node_id], xgc->data.set_mat.tetgen_zshift[node_id]);
    } else if (node_id >= SET_MAT_TETRAHEDRON) {
        SvTetrahedron tthn;
        tthn = g_array_index(xgc->data.set_mat.tetrahedrons, SvTetrahedron, index);
        g_snprintf(buff, buff_size, "Tetrahedron at %g %g %g ... %g %g %g ... %g %g %g ... %g %g %g", tthn.pnt1[0], tthn.pnt1[1], tthn.pnt1[2], tthn.pnt2[0], tthn.pnt2[1], tthn.pnt2[2], tthn.pnt3[0], tthn.pnt3[1], tthn.pnt3[2], tthn.pnt4[0], tthn.pnt4[1], tthn.pnt4[2]);
    }
}

void mt_create_tree(XGControls *xgc)
{
    GtkTextIter     start;
    GtkTextIter     end;
    GtkWidget       *dialog;
    gchar           buff[256];
    GtkTreeIter     iter, child;
    gint            m, i, /*active,*/isobject;
    //gchar           *contents;
    GError          *err = NULL;
    gchar           *textbuf;

    SvSphere sx;
    SvVoxel vx;
    SvCylinder cx;
    SvCone cnx;
    SvRCone rcnx;
    SvTetrahedron tx;
    SvGwydd gx;
    gint ngtot = 0;
    //SvMatProp mat;

    if (NULL == xgc->tmpmatfilename)
        xgc->tmpmatfilename = get_temporary_mat_filename(xgc);

    if (xgc->tmpmatfilename) {
        mt_remeber_expanded_rows_and_selection(xgc);

        gtk_text_buffer_get_start_iter(xgc->tb_mat, &start);
        gtk_text_buffer_get_end_iter(xgc->tb_mat, &end);
        textbuf = gtk_text_buffer_get_text(xgc->tb_mat, &start, &end, TRUE);
        if (!g_file_set_contents(xgc->tmpfilename, textbuf, -1, &err)) {
            dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                GTK_DIALOG_DESTROY_WITH_PARENT,
                GTK_MESSAGE_ERROR,
                GTK_BUTTONS_CLOSE,
                "Error writing to temporary file: %s",
                g_strdup(err->message));
            gtk_dialog_run(GTK_DIALOG(dialog));
            gtk_widget_destroy(dialog);
        }
        //printf("dir: %s\n", g_get_current_dir());

        /*clear material settings = alloc material data*/
        clear_settings_mat(&xgc->data.set_mat, TRUE);
        alloc_set_mat(&xgc->data.set_mat);
        //////////////////////////////////////////////////////////////////////////

        if (parse_settings_mat(xgc->tmpfilename, &(xgc->data.set), &(xgc->data.set_mat), FALSE)) {
            /* no errors in matfile -> create tree*/

            gtk_tree_store_clear(xgc->ts_mat);

            gtk_tree_store_append(xgc->ts_mat, &iter, NULL);
            if (xgc->data.set.sm.in_vector_filename)
                g_snprintf(buff, sizeof(buff), "File: %s", xgc->data.set.sm.in_vector_filename);
            else
                g_snprintf(buff, sizeof(buff), "File: undefined");


            gtk_tree_store_set(xgc->ts_mat, &iter, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, -1, COLUMN_SHOW_TOGGLE, FALSE, -1);

            for (m = 0; m < (gint)MIN(xgc->data.set_mat.spheres->len, TREE_MAXENT); m++) {
                sx = g_array_index(xgc->data.set_mat.spheres, SvSphere, m);

                mt_tree_assemble_row_text(xgc, buff, sizeof(buff), SET_MAT_SPHERE, m);
                gtk_tree_store_append(xgc->ts_mat, &child, &iter);
                //g_snprintf(buff, sizeof(buff), "Sphere radius %g at %g, %g, %g", sx.radius, sx.pnt1[0], sx.pnt1[1], sx.pnt1[2]);
                gtk_tree_store_set(xgc->ts_mat, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MAT_SPHERE + m, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
            for (m = 0; m < (gint)MIN(xgc->data.set_mat.voxels->len, TREE_MAXENT); m++) {
                vx = g_array_index(xgc->data.set_mat.voxels, SvVoxel, m);

                mt_tree_assemble_row_text(xgc, buff, sizeof(buff), SET_MAT_VOXEL, m);
                gtk_tree_store_append(xgc->ts_mat, &child, &iter);
                //g_snprintf(buff, sizeof(buff), "Box at (%g %g %g), (%g %g %g)", vx.pnt1[0], vx.pnt1[1], vx.pnt1[2], vx.pnt2[0], vx.pnt2[1], vx.pnt2[2]);
                gtk_tree_store_set(xgc->ts_mat, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MAT_VOXEL + m, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
            for (m = 0; m < (gint)MIN(xgc->data.set_mat.cylinders->len, TREE_MAXENT); m++) {
                cx = g_array_index(xgc->data.set_mat.cylinders, SvCylinder, m);

                mt_tree_assemble_row_text(xgc, buff, sizeof(buff), SET_MAT_CYLINDER, m);
                gtk_tree_store_append(xgc->ts_mat, &child, &iter);
                //g_snprintf(buff, sizeof(buff), "Cylinder radius %g at  (%g %g %g), (%g %g %g)", cx.radius, cx.pnt1[0], cx.pnt1[1], cx.pnt1[2], cx.pnt2[0], cx.pnt2[1], cx.pnt2[2]);
                gtk_tree_store_set(xgc->ts_mat, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MAT_CYLINDER + m, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
            for (m = 0; m < (gint)MIN(xgc->data.set_mat.cones->len, TREE_MAXENT); m++) {
                cnx = g_array_index(xgc->data.set_mat.cones, SvCone, m);

                mt_tree_assemble_row_text(xgc, buff, sizeof(buff), SET_MAT_CONE, m);
                gtk_tree_store_append(xgc->ts_mat, &child, &iter);
                //g_snprintf(buff, sizeof(buff), "Cone radius %g at (%g %g %g), (%g %g %g)", cnx.radius, cnx.pnt1[0], cnx.pnt1[1], cnx.pnt1[2], cnx.pnt2[0], cnx.pnt2[1], cnx.pnt2[2]);
                gtk_tree_store_set(xgc->ts_mat, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MAT_CONE + m, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
            for (m = 0; m < (gint)MIN(xgc->data.set_mat.rcones->len, TREE_MAXENT); m++) {
                rcnx = g_array_index(xgc->data.set_mat.rcones, SvRCone, m);

                mt_tree_assemble_row_text(xgc, buff, sizeof(buff), SET_MAT_RCONE, m);
                gtk_tree_store_append(xgc->ts_mat, &child, &iter);
                //g_snprintf(buff, sizeof(buff), "Cut cone %g...%g at (%g %g %g), (%g %g %g)", rcnx.radius1, rcnx.radius2, rcnx.pnt1[0], rcnx.pnt1[1], rcnx.pnt1[2], rcnx.pnt2[0], rcnx.pnt2[1], rcnx.pnt2[2]);
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
                    mt_tree_assemble_row_text(xgc, buff, sizeof(buff), SET_MAT_TETRAHEDRON, m);
                    gtk_tree_store_append(xgc->ts_mat, &child, &iter);
                    //g_snprintf(buff, sizeof(buff), "Tetrahedron starting at %g, %g, %g", tx.pnt1[0], tx.pnt1[1], tx.pnt1[2]);
                    gtk_tree_store_set(xgc->ts_mat, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MAT_TETRAHEDRON + m, COLUMN_SHOW_TOGGLE, TRUE, -1);
                }
            }
            for (m = 0; m < MIN(xgc->data.set_mat.ntetgens, TREE_MAXENT); m++) {
                mt_tree_assemble_row_text(xgc, buff, sizeof(buff), SET_MAT_MESH, m);
                gtk_tree_store_append(xgc->ts_mat, &child, &iter);
                //g_snprintf(buff, sizeof(buff), "mesh (%s) at %g, %g, %g", xgc->data.set_mat.tetgen_filebase[m], xgc->data.set_mat.tetgen_xshift[m], xgc->data.set_mat.tetgen_yshift[m], xgc->data.set_mat.tetgen_zshift[m]);
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
        } else {
            g_snprintf(buff, sizeof(buff), MSG_SB_ERROR_PARSING_PARFILE);
            gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);

            gtk_widget_queue_draw(xgc->view_scene);

            xgc->par_tree_success = FALSE;
            xgc->par_file_success = FALSE;

            gtk_widget_set_sensitive(xgc->page_par_tree, xgc->par_tree_success);
            gtk_widget_set_sensitive(xgc->page_mat_tree, xgc->mat_tree_success);

            dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                GTK_DIALOG_DESTROY_WITH_PARENT,
                GTK_MESSAGE_WARNING,
                GTK_BUTTONS_OK,
                MSG_MB_ERROR_PARSING_PARFILE);
            gtk_dialog_run(GTK_DIALOG(dialog));
            gtk_widget_destroy(dialog);

            gtk_notebook_set_current_page(GTK_NOTEBOOK(xgc->mat_notebook), 1);
        }
        //g_free(filename);

        mt_restore_expanded_rows_and_selection(xgc);

        g_snprintf(buff, sizeof(buff), "Ready");
        gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);

        gtk_widget_queue_draw(xgc->view_scene);

        xgc->mat_tree_success = TRUE;
        xgc->mat_file_success = TRUE;

        gtk_widget_set_sensitive(xgc->page_mat_tree, xgc->mat_tree_success);
    }
}

void mt_remeber_expanded_rows_and_selection(XGControls *xgc)
{
    /* remember expanded rows and selection */

    GtkTreePath *path_expand = NULL, *path_select = NULL;
    gint i;

    for (i = 0; i < TREE_VIEW_MAT_ROOT_ROWS; i++) {
        path_expand = gtk_tree_path_new_from_indices(i, -1);
        xgc->tv_mat_expanded[i] = gtk_tree_view_row_expanded(GTK_TREE_VIEW(xgc->tv_mat), path_expand);
        gtk_tree_path_free(path_expand);
    }

    if (xgc->mat_selid != -1) {
        path_select = gtk_tree_model_get_path(GTK_TREE_MODEL(xgc->ts_mat), &(xgc->mat_seliter));
        if (NULL != xgc->tv_mat_path_select_string)
            g_free(xgc->tv_mat_path_select_string);
        xgc->tv_mat_path_select_string = gtk_tree_path_to_string(path_select);
    }

    xgc->page_mat_hscroll_val = gtk_adjustment_get_value(gtk_scrolled_window_get_hadjustment(GTK_SCROLLED_WINDOW(xgc->cs_mat.sw_values)));
    xgc->page_mat_vscroll_val = gtk_adjustment_get_value(gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(xgc->cs_mat.sw_values)));
}

void mt_restore_expanded_rows_and_selection(XGControls *xgc)
{
    /* restore expanded rows and selection */

    gint i;
    GtkTreePath *path_expand = NULL, *path_select = NULL;
    GtkTreeSelection *selection = NULL;

    for (i = 0; i < TREE_VIEW_MAT_ROOT_ROWS; i++) {
        if (xgc->tv_mat_expanded[i]) {
            path_expand = gtk_tree_path_new_from_indices(i, -1);
            gtk_tree_view_expand_row(GTK_TREE_VIEW(xgc->tv_mat), path_expand, TRUE);
            gtk_tree_path_free(path_expand);
        }
    }

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(xgc->tv_mat));
    if (xgc->tv_mat_path_select_string != NULL) {
        path_select = gtk_tree_path_new_from_string(xgc->tv_mat_path_select_string);
        if (path_select != NULL)
            gtk_tree_selection_select_path(selection, path_select);
    }

    gtk_adjustment_set_value(gtk_scrolled_window_get_hadjustment(GTK_SCROLLED_WINDOW(xgc->cs.sw_values)), xgc->page_hscroll_val);
    gtk_adjustment_set_value(gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(xgc->cs.sw_values)), xgc->page_vscroll_val);
}