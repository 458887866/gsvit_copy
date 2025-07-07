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


/* par_tree.c :
*  Parameter controls functions
*/

#include "xs_par_controls.h"
#include "xs_par_tree.h"
#include "xs_layer_treeview.h"
#include "write.h"
#include "string.h"

static gboolean set_values_general(XGControls* xgc, SvSet *set, SvType type, gint id);
static gboolean set_values_sf(XGControls* xgc, SvSet *set, SvType type, gint id);
static gboolean set_values_output(XGControls* xgc, SvSet *set, SvType type, gint id);
static gboolean set_values_material(XGControls* xgc, SvSet *set, SvType type, gint id);
static gboolean set_values_nfff(XGControls* xgc, SvSet *set, SvType type, gint id);

static void set_boundary_controls_sensitive(XGControls *xgc, gint i, gboolean val);



void
par_controls_changed(XGControls *xgc)
{
    if (xgc->par_lock_controls_update) {
        return;
    }

    gint i;
    guint id, ugval;
    gchar buff[256];
    GtkWidget *focus;
    gchar size_um_buff[256] = {0};

    // remember focus at the beginning and restore at the end of ss_controls_changed()
    // Note: this works only for GtkEntry not for GtkSpinButton
    focus = gtk_window_get_focus(GTK_WINDOW(xgc->toplevel));

    if (xgc->par_selid == SET_POOL) {
        /* I. 1. Computational domain */
        xgc->data.set.sp.xres = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.cc.comp_domain_size_x_spin));
        xgc->data.set.sp.yres = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.cc.comp_domain_size_y_spin));
        xgc->data.set.sp.zres = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.cc.comp_domain_size_z_spin));

        xgc->data.set.sp.dx = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.cc.comp_domain_spacing_x_spin))*1e-6;
        xgc->data.set.sp.dy = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.cc.comp_domain_spacing_y_spin))*1e-6;
        xgc->data.set.sp.dz = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.cc.comp_domain_spacing_z_spin))*1e-6;

        g_snprintf(size_um_buff, sizeof(size_um_buff), "%g", xgc->data.set.sp.xres*xgc->data.set.sp.dx*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.cc.comp_domain_size_x_um), size_um_buff);
        g_snprintf(size_um_buff, sizeof(size_um_buff), "%g", xgc->data.set.sp.yres*xgc->data.set.sp.dy*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.cc.comp_domain_size_y_um), size_um_buff);
        g_snprintf(size_um_buff, sizeof(size_um_buff), "%g", xgc->data.set.sp.zres*xgc->data.set.sp.dz*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.cc.comp_domain_size_z_um), size_um_buff);
    } else if (xgc->par_selid == SET_BASIC) {
        /* II. 1. Basic parameters */
        xgc->data.set.sc.nsteps = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.bac.basic_nsteps_spin));
        xgc->data.set.sc.verbose = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.bac.basic_verbose_spin));
        xgc->data.set.sc.nthreads = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.bac.basic_nthreads_spin));

        xgc->data.set.sc.usegpu = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.bac.basic_usegpu));
        gtk_widget_set_sensitive(xgc->cs.bac.basic_ugpu_index_spin, xgc->data.set.sc.usegpu);

        ugval = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.bac.basic_ugpu_index_spin));
        xgc->data.set.sc.ugpu[0] = xgc->data.set.sc.ugpu[1] = xgc->data.set.sc.ugpu[2] = xgc->data.set.sc.ugpu[3] = 0;
        if (ugval == 0)
            xgc->data.set.sc.ugpu[0] = 1;
        else if (ugval == 1)
            xgc->data.set.sc.ugpu[1] = 1;
        else if (ugval == 2)
            xgc->data.set.sc.ugpu[2] = 1;
        else if (ugval == 3)
            xgc->data.set.sc.ugpu[3] = 1;
    } else if (xgc->par_selid >= SET_PSOURCE && xgc->par_selid < SET_POUT) {
        id = xgc->par_selid - SET_PSOURCE;

        /* III. 1. Point source properties */
        xgc->data.set.ss.pnts[id].point_origin_position_i = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.point_origin_position_i_spin));
        xgc->data.set.ss.pnts[id].point_origin_position_j = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.point_origin_position_j_spin));
        xgc->data.set.ss.pnts[id].point_origin_position_k = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.point_origin_position_k_spin));
        xgc->data.set.ss.pnts[id].point_origin_theta = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.point_origin_theta_spin)) / 180 * G_PI;
        xgc->data.set.ss.pnts[id].point_origin_phi = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.point_origin_phi_spin)) / 180 * G_PI;

        /* III. B. Source group */
        xgc->data.set.ss.pnts[id].source_mode = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.sc.source_mode));

        g_free(xgc->data.set.ss.pnts[id].source_filename);
        xgc->data.set.ss.pnts[id].source_filename = g_strdup(gtk_entry_get_text(GTK_ENTRY(xgc->cs.sc.source_filename)));
        /*if (xgc->data.set.ss.pnts[id].source_mode == 0 && NULL != xgc->data.sfb.ps_source_filename) {
        g_free(xgc->data.set.ss.pnts[id].source_filename);
        xgc->data.set.ss.pnts[id].source_filename = g_strdup(xgc->data.sfb.ps_source_filename);
        }*/

        xgc->data.set.ss.pnts[id].source_amplitude = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_amplitude_spin));
        xgc->data.set.ss.pnts[id].source_wl = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_wavelength_spin))*1e-6;
        xgc->data.set.ss.pnts[id].source_pulsewidth = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_pulsewidth_spin));

        //get_ps_store_text(xgc, id, buff, sizeof(buff));
        //gtk_tree_store_set(xgc->ts_par, &xgc->par_seliter, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, id-SET_PSOURCE, COLUMN_SHOW_TOGGLE, TRUE, -1);
    } else if (xgc->par_selid == SET_SF) {
        /* III. B. Source group */
        xgc->data.set.ss.sf.source_mode = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.sc.source_mode));
        g_free(xgc->data.set.ss.sf.source_filename);
        xgc->data.set.ss.sf.source_filename = g_strdup(gtk_entry_get_text(GTK_ENTRY(xgc->cs.sc.source_filename)));
        /*if (xgc->data.set.ss.sf.source_mode == 0 && NULL != xgc->data.sfb.sf_source_filename) {
        g_free(xgc->data.set.ss.sf.source_filename);
        xgc->data.set.ss.sf.source_filename = g_strdup(xgc->data.sfb.sf_source_filename);
        }*/

        xgc->data.set.ss.sf.source_amplitude = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_amplitude_spin));
        xgc->data.set.ss.sf.source_wl = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_wavelength_spin))*1e-6;
        xgc->data.set.ss.sf.source_pulsewidth = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_pulsewidth_spin));

        /* III. C. Incident angle group */
        xgc->data.set.ss.sf.ia_theta = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.ia_theta_spin)) / 180 * G_PI;
        xgc->data.set.ss.sf.ia_phi = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.ia_phi_spin)) / 180 * G_PI;
        xgc->data.set.ss.sf.ia_psi = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.ia_psi_spin)) / 180 * G_PI;

        //get_sf_store_text(xgc, buff, sizeof(buff));
        //gtk_tree_store_set(xgc->ts_par, &xgc->par_seliter, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_SF, COLUMN_SHOW_TOGGLE, TRUE, -1);
    } else if (xgc->par_selid == SET_TSF) {
        /* III. A. Source box group */
        xgc->data.set.ss.tsf.box_i0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_i0_spin));
        xgc->data.set.ss.tsf.box_j0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_j0_spin));
        xgc->data.set.ss.tsf.box_k0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_k0_spin));
        xgc->data.set.ss.tsf.box_in = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_in_spin));
        xgc->data.set.ss.tsf.box_jn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_jn_spin));
        xgc->data.set.ss.tsf.box_kn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_kn_spin));

        xgc->data.set.ss.tsf.box_boundary_skipi0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipi0));
        xgc->data.set.ss.tsf.box_boundary_skipj0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipj0));
        xgc->data.set.ss.tsf.box_boundary_skipk0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipk0));
        xgc->data.set.ss.tsf.box_boundary_skipin = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipin));
        xgc->data.set.ss.tsf.box_boundary_skipjn = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipjn));
        xgc->data.set.ss.tsf.box_boundary_skipkn = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipkn));

        if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_enable_skipdepth)) == TRUE) {
            xgc->data.set.ss.tsf.box_boundary_skipdepth = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_skipdepth_spin));
            if (xgc->data.set.ss.tsf.box_boundary_skipdepth == -1)
                xgc->data.set.ss.tsf.box_boundary_skipdepth = 0;
        } else
            xgc->data.set.ss.tsf.box_boundary_skipdepth = -1;   /* -1 for non-conformal, none skip depth */

                                                                /* III. B. Source group */
        xgc->data.set.ss.tsf.source_mode = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.sc.source_mode));

        g_free(xgc->data.set.ss.tsf.source_filename);
        xgc->data.set.ss.tsf.source_filename = g_strdup(gtk_entry_get_text(GTK_ENTRY(xgc->cs.sc.source_filename)));
        /*if (xgc->data.set.ss.tsf.source_mode == 0 && NULL != xgc->data.sfb.tsf_source_filename) {
        g_free(xgc->data.set.ss.tsf.source_filename);
        xgc->data.set.ss.tsf.source_filename = g_strdup(xgc->data.sfb.tsf_source_filename);
        }*/

        xgc->data.set.ss.tsf.source_amplitude = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_amplitude_spin));
        xgc->data.set.ss.tsf.source_wl = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_wavelength_spin))*1e-6;
        xgc->data.set.ss.tsf.source_pulsewidth = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_pulsewidth_spin));

        /* III. C. Incident angle group */
        xgc->data.set.ss.tsf.ia_theta = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.ia_theta_spin)) / 180 * G_PI;
        xgc->data.set.ss.tsf.ia_phi = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.ia_phi_spin)) / 180 * G_PI;
        xgc->data.set.ss.tsf.ia_psi = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.ia_psi_spin)) / 180 * G_PI;

        /* III. F. Z multiplier group */
        xgc->data.set.ss.tsf.gaussian = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.gaussian_mult_enable));
        xgc->data.set.ss.tsf.gaussian_fxpos = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_center_pos_i_spin));
        xgc->data.set.ss.tsf.gaussian_fypos = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_center_pos_j_spin));
        xgc->data.set.ss.tsf.gaussian_rx = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_radius_i_spin));
        xgc->data.set.ss.tsf.gaussian_ry = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_radius_j_spin));

        xgc->data.set.ss.tsf.radial = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.radial_mult_enable));
        xgc->data.set.ss.tsf.radial_fxpos = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.radial_mult_center_pos_i_spin));
        xgc->data.set.ss.tsf.radial_fypos = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.radial_mult_center_pos_j_spin));
        xgc->data.set.ss.tsf.radial_rx = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.radial_mult_radius_i_spin));
        xgc->data.set.ss.tsf.radial_ry = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.radial_mult_radius_j_spin));

        xgc->data.set.ss.tsf.fiber = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.fiber_mult_enable));
        xgc->data.set.ss.tsf.fiber_fxpos = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_center_pos_i_spin));
        xgc->data.set.ss.tsf.fiber_fypos = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_center_pos_j_spin));
        xgc->data.set.ss.tsf.fiber_radius = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_radius_spin));
        xgc->data.set.ss.tsf.fiber_cutoff = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_cutoff_spin));
        xgc->data.set.ss.tsf.fiber_epsilon_core = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_epsilon_core_spin));
        xgc->data.set.ss.tsf.fiber_epsilon_cladding = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_epsilon_cladding_spin));

        //get_tsf_store_text(xgc, buff, sizeof(buff));        
        //gtk_tree_store_set(xgc->ts_par, &xgc->par_seliter, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_TSF, COLUMN_SHOW_TOGGLE, TRUE, -1);
    } else if (xgc->par_selid == SET_TSFF) {
        /* III. A. Source box group */
        xgc->data.set.ss.tsff.box_i0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_i0_spin));
        xgc->data.set.ss.tsff.box_j0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_j0_spin));
        xgc->data.set.ss.tsff.box_k0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_k0_spin));
        xgc->data.set.ss.tsff.box_in = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_in_spin));
        xgc->data.set.ss.tsff.box_jn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_jn_spin));
        xgc->data.set.ss.tsff.box_kn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_kn_spin));

        xgc->data.set.ss.tsff.box_boundary_skipi0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipi0));
        xgc->data.set.ss.tsff.box_boundary_skipj0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipj0));
        xgc->data.set.ss.tsff.box_boundary_skipk0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipk0));
        xgc->data.set.ss.tsff.box_boundary_skipin = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipin));
        xgc->data.set.ss.tsff.box_boundary_skipjn = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipjn));
        xgc->data.set.ss.tsff.box_boundary_skipkn = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipkn));

        if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_enable_skipdepth)) == TRUE) {
            xgc->data.set.ss.tsff.box_boundary_skipdepth = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_skipdepth_spin));
            if (xgc->data.set.ss.tsff.box_boundary_skipdepth == -1)
                xgc->data.set.ss.tsff.box_boundary_skipdepth = 0;
        } else
            xgc->data.set.ss.tsff.box_boundary_skipdepth = -1;   /* -1 for non-conformal, none skip depth */

                                                                 /* III. B. Source group */
        xgc->data.set.ss.tsff.source_mode = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.sc.source_mode));

        g_free(xgc->data.set.ss.tsff.source_filename);
        xgc->data.set.ss.tsff.source_filename = g_strdup(gtk_entry_get_text(GTK_ENTRY(xgc->cs.sc.source_filename)));
        /*if (xgc->data.set.ss.tsff.source_mode == 0 && NULL != xgc->data.sfb.tsff_source_filename) {
        g_free(xgc->data.set.ss.tsff.source_filename);
        xgc->data.set.ss.tsff.source_filename = g_strdup(xgc->data.sfb.tsff_source_filename);
        }*/

        xgc->data.set.ss.tsff.source_amplitude = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_amplitude_spin));
        xgc->data.set.ss.tsff.source_wl = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_wavelength_spin))*1e-6;
        xgc->data.set.ss.tsff.source_pulsewidth = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_pulsewidth_spin));

        /* III. D. Focused source group */
        xgc->data.set.ss.tsff.focused_thetamax = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fs_thetamax_deg_spin)) / 180 * G_PI;
        xgc->data.set.ss.tsff.focused_fdist = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fs_fdist_spin));
        xgc->data.set.ss.tsff.focused_pol = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fs_pol_deg_spin)) / 180 * G_PI;
        xgc->data.set.ss.tsff.focused_nip = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fs_nip_spin));
        xgc->data.set.ss.tsff.focused_mip = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fs_mip_spin));

        //get_tsff_store_text(xgc, buff, sizeof(buff));
        //gtk_tree_store_set(xgc->ts_par, &xgc->par_seliter, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_TSFF, COLUMN_SHOW_TOGGLE, TRUE, -1);
    } else if (xgc->par_selid == SET_LTSF) {
        /* III. A. Source box group */
        xgc->data.set.ss.ltsf.box_i0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_i0_spin));
        xgc->data.set.ss.ltsf.box_j0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_j0_spin));
        xgc->data.set.ss.ltsf.box_k0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_k0_spin));
        xgc->data.set.ss.ltsf.box_in = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_in_spin));
        xgc->data.set.ss.ltsf.box_jn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_jn_spin));
        xgc->data.set.ss.ltsf.box_kn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_kn_spin));

        xgc->data.set.ss.ltsf.box_boundary_skipi0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipi0));
        xgc->data.set.ss.ltsf.box_boundary_skipj0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipj0));
        xgc->data.set.ss.ltsf.box_boundary_skipk0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipk0));
        xgc->data.set.ss.ltsf.box_boundary_skipin = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipin));
        xgc->data.set.ss.ltsf.box_boundary_skipjn = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipjn));
        xgc->data.set.ss.ltsf.box_boundary_skipkn = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipkn));

        if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_enable_skipdepth)) == TRUE) {
            xgc->data.set.ss.ltsf.box_boundary_skipdepth = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_skipdepth_spin));
            if (xgc->data.set.ss.ltsf.box_boundary_skipdepth == -1)
                xgc->data.set.ss.ltsf.box_boundary_skipdepth = 0;
        } else
            xgc->data.set.ss.ltsf.box_boundary_skipdepth = -1;   /* -1 for non-conformal, none skip depth */

                                                                 /* III. B. Source group */
        xgc->data.set.ss.ltsf.source_mode = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.sc.source_mode));

        g_free(xgc->data.set.ss.ltsf.source_filename);
        xgc->data.set.ss.ltsf.source_filename = g_strdup(gtk_entry_get_text(GTK_ENTRY(xgc->cs.sc.source_filename)));
        /*if (xgc->data.set.ss.ltsf.source_mode == 0 && NULL != xgc->data.sfb.ltsf_source_filename) {
        g_free(xgc->data.set.ss.ltsf.source_filename);
        xgc->data.set.ss.ltsf.source_filename = g_strdup(xgc->data.sfb.ltsf_source_filename);
        }*/

        xgc->data.set.ss.ltsf.source_amplitude = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_amplitude_spin));
        xgc->data.set.ss.ltsf.source_wl = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_wavelength_spin))*1e-6;
        xgc->data.set.ss.ltsf.source_pulsewidth = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_pulsewidth_spin));

        /* III. C. Incident angle group */
        xgc->data.set.ss.ltsf.ia_theta = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.ia_theta_spin)) / 180 * G_PI;
        xgc->data.set.ss.ltsf.ia_phi = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.ia_phi_spin)) / 180 * G_PI;
        xgc->data.set.ss.ltsf.ia_psi = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.ia_psi_spin)) / 180 * G_PI;

        /* III. E. Layered source group */
        for (i = 0; i < xgc->data.set.ss.ltsf.layered_count; i++)
            g_free(xgc->data.set.ss.ltsf.layered_material[i]);

        LTV_Data data;
        //get_ltv_items(GTK_TREE_VIEW(xgc->cs.sc.ls_prop_wrapper), (gpointer)&data);
        get_ltv_items((TreeViewWrapper*)xgc->cs.sc.ls_prop_wrapper, (gpointer)&data);
        xgc->data.set.ss.ltsf.layered_count = data.nlayers;        
        for (i = 0; i < xgc->data.set.ss.ltsf.layered_count; i++) {
            xgc->data.set.ss.ltsf.layered_zpos[i] = data.zpos[i];
            xgc->data.set.ss.ltsf.layered_epsilon[i] = data.epsilon[i];
            xgc->data.set.ss.ltsf.layered_sigma[i] = data.sigma[i];
            xgc->data.set.ss.ltsf.layered_mu[i] = data.mu[i];
            xgc->data.set.ss.ltsf.layered_sigast[i] = data.sigast[i];
            xgc->data.set.ss.ltsf.layered_material[i] = g_strdup(data.material[i]);
        }

        /* III. F. Z multiplier group */
        xgc->data.set.ss.ltsf.gaussian = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.gaussian_mult_enable));
        xgc->data.set.ss.ltsf.gaussian_fxpos = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_center_pos_i_spin));
        xgc->data.set.ss.ltsf.gaussian_fypos = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_center_pos_j_spin));
        xgc->data.set.ss.ltsf.gaussian_rx = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_radius_i_spin));
        xgc->data.set.ss.ltsf.gaussian_ry = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_radius_j_spin));

        xgc->data.set.ss.ltsf.radial = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.radial_mult_enable));
        xgc->data.set.ss.ltsf.radial_fxpos = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.radial_mult_center_pos_i_spin));
        xgc->data.set.ss.ltsf.radial_fypos = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.radial_mult_center_pos_j_spin));
        xgc->data.set.ss.ltsf.radial_rx = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.radial_mult_radius_i_spin));
        xgc->data.set.ss.ltsf.radial_ry = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.radial_mult_radius_j_spin));

        xgc->data.set.ss.ltsf.fiber = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.fiber_mult_enable));
        xgc->data.set.ss.ltsf.fiber_fxpos = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_center_pos_i_spin));
        xgc->data.set.ss.ltsf.fiber_fypos = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_center_pos_j_spin));
        xgc->data.set.ss.ltsf.fiber_radius = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_radius_spin));
        xgc->data.set.ss.ltsf.fiber_cutoff = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_cutoff_spin));
        xgc->data.set.ss.ltsf.fiber_epsilon_core = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_epsilon_core_spin));
        xgc->data.set.ss.ltsf.fiber_epsilon_cladding = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_epsilon_cladding_spin));

        //get_ltsf_store_text(xgc, buff, sizeof(buff));
        //gtk_tree_store_set(xgc->ts_par, &xgc->par_seliter, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_LTSF, COLUMN_SHOW_TOGGLE, TRUE, -1);
    } else if (xgc->par_selid == SET_LTSFF) {
        /* III. A. Source box group */
        xgc->data.set.ss.ltsff.box_i0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_i0_spin));
        xgc->data.set.ss.ltsff.box_j0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_j0_spin));
        xgc->data.set.ss.ltsff.box_k0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_k0_spin));
        xgc->data.set.ss.ltsff.box_in = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_in_spin));
        xgc->data.set.ss.ltsff.box_jn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_jn_spin));
        xgc->data.set.ss.ltsff.box_kn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_kn_spin));

        xgc->data.set.ss.ltsff.box_boundary_skipi0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipi0));
        xgc->data.set.ss.ltsff.box_boundary_skipj0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipj0));
        xgc->data.set.ss.ltsff.box_boundary_skipk0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipk0));
        xgc->data.set.ss.ltsff.box_boundary_skipin = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipin));
        xgc->data.set.ss.ltsff.box_boundary_skipjn = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipjn));
        xgc->data.set.ss.ltsff.box_boundary_skipkn = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipkn));

        if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_enable_skipdepth)) == TRUE) {
            xgc->data.set.ss.ltsff.box_boundary_skipdepth = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_skipdepth_spin));
            if (xgc->data.set.ss.ltsff.box_boundary_skipdepth == -1)
                xgc->data.set.ss.ltsff.box_boundary_skipdepth = 0;
        } else
            xgc->data.set.ss.ltsff.box_boundary_skipdepth = -1;   /* -1 for non-conformal, none skip depth */

                                                                  /* III. B. Source group */
        xgc->data.set.ss.ltsff.source_mode = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.sc.source_mode));

        g_free(xgc->data.set.ss.ltsff.source_filename);
        xgc->data.set.ss.ltsff.source_filename = g_strdup(gtk_entry_get_text(GTK_ENTRY(xgc->cs.sc.source_filename)));
        /*if (xgc->data.set.ss.ltsff.source_mode == 0 && NULL != xgc->data.sfb.ltsff_source_filename) {
        g_free(xgc->data.set.ss.ltsff.source_filename);
        xgc->data.set.ss.ltsff.source_filename = g_strdup(xgc->data.sfb.ltsff_source_filename);
        }*/

        xgc->data.set.ss.ltsff.source_amplitude = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_amplitude_spin));
        xgc->data.set.ss.ltsff.source_wl = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_wavelength_spin))*1e-6;
        xgc->data.set.ss.ltsff.source_pulsewidth = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_pulsewidth_spin));

        /* III. D. Focused source group */
        xgc->data.set.ss.ltsff.focused_thetamax = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fs_thetamax_deg_spin)) / 180 * G_PI;
        xgc->data.set.ss.ltsff.focused_fdist = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fs_fdist_spin));
        xgc->data.set.ss.ltsff.focused_pol = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fs_pol_deg_spin)) / 180 * G_PI;
        xgc->data.set.ss.ltsff.focused_nip = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fs_nip_spin));
        xgc->data.set.ss.ltsff.focused_mip = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.sc.fs_mip_spin));

        /* III. E. Layered source group */
        for (i = 0; i < xgc->data.set.ss.ltsff.layered_count; i++)
            g_free(xgc->data.set.ss.ltsff.layered_material[i]);

        LTV_Data data;
        //get_ltv_items(GTK_TREE_VIEW(xgc->cs.sc.ls_prop_wrapper), (gpointer)&data);
        get_ltv_items((TreeViewWrapper*)xgc->cs.sc.ls_prop_wrapper, (gpointer)&data);
        xgc->data.set.ss.ltsff.layered_count = data.nlayers;
        for (i = 0; i < xgc->data.set.ss.ltsff.layered_count; i++) {
            xgc->data.set.ss.ltsff.layered_zpos[i] = data.zpos[i];
            xgc->data.set.ss.ltsff.layered_epsilon[i] = data.epsilon[i];
            xgc->data.set.ss.ltsff.layered_sigma[i] = data.sigma[i];
            xgc->data.set.ss.ltsff.layered_mu[i] = data.mu[i];
            xgc->data.set.ss.ltsff.layered_sigast[i] = data.sigast[i];
            xgc->data.set.ss.ltsff.layered_material[i] = g_strdup(data.material[i]);
        }

        //get_ltsff_store_text(xgc, buff, sizeof(buff));
        //gtk_tree_store_set(xgc->ts_par, &xgc->par_seliter, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_LTSFF, COLUMN_SHOW_TOGGLE, TRUE, -1);
    } else if (xgc->par_selid == SET_BND) {
        /* IV. 1. Boundary conditions */
        xgc->data.set.sb.bx0 = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.boc.btype[0]));
        xgc->data.set.sb.depth_bx0 = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bdepth[0]));
        xgc->data.set.sb.m_bx0 = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bm[0]));
        xgc->data.set.sb.sigma_bx0 = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bsigma[0]));
        xgc->data.set.sb.a_bx0 = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.ba[0]));
        xgc->data.set.sb.kappa_bx0 = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bkappa[0]));

        xgc->data.set.sb.bxn = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.boc.btype[1]));
        xgc->data.set.sb.depth_bxn = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bdepth[1]));
        xgc->data.set.sb.m_bxn = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bm[1]));
        xgc->data.set.sb.sigma_bxn = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bsigma[1]));
        xgc->data.set.sb.a_bxn = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.ba[1]));
        xgc->data.set.sb.kappa_bxn = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bkappa[1]));

        xgc->data.set.sb.by0 = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.boc.btype[2]));
        xgc->data.set.sb.depth_by0 = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bdepth[2]));
        xgc->data.set.sb.m_by0 = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bm[2]));
        xgc->data.set.sb.sigma_by0 = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bsigma[2]));
        xgc->data.set.sb.a_by0 = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.ba[2]));
        xgc->data.set.sb.kappa_by0 = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bkappa[2]));

        xgc->data.set.sb.byn = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.boc.btype[3]));
        xgc->data.set.sb.depth_byn = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bdepth[3]));
        xgc->data.set.sb.m_byn = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bm[3]));
        xgc->data.set.sb.sigma_byn = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bsigma[3]));
        xgc->data.set.sb.a_byn = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.ba[3]));
        xgc->data.set.sb.kappa_byn = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bkappa[3]));

        xgc->data.set.sb.bz0 = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.boc.btype[4]));
        xgc->data.set.sb.depth_bz0 = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bdepth[4]));
        xgc->data.set.sb.m_bz0 = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bm[4]));
        xgc->data.set.sb.sigma_bz0 = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bsigma[4]));
        xgc->data.set.sb.a_bz0 = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.ba[4]));
        xgc->data.set.sb.kappa_bz0 = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bkappa[4]));

        xgc->data.set.sb.bzn = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.boc.btype[5]));
        xgc->data.set.sb.depth_bzn = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bdepth[5]));
        xgc->data.set.sb.m_bzn = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bm[5]));
        xgc->data.set.sb.sigma_bzn = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bsigma[5]));
        xgc->data.set.sb.a_bzn = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.ba[5]));
        xgc->data.set.sb.kappa_bzn = gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.bkappa[5]));

        xgc->data.set.smb.bx0 = (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.boc.mb[0]))) ? SV_BOUNDARY_PERIODIC : SV_BOUNDARY_NONE;
        xgc->data.set.smb.bx0pos = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.mbpos[0]));

        xgc->data.set.smb.bxn = (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.boc.mb[1]))) ? SV_BOUNDARY_PERIODIC : SV_BOUNDARY_NONE;
        xgc->data.set.smb.bxnpos = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.mbpos[1]));

        xgc->data.set.smb.by0 = (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.boc.mb[2]))) ? SV_BOUNDARY_PERIODIC : SV_BOUNDARY_NONE;
        xgc->data.set.smb.by0pos = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.mbpos[2]));

        xgc->data.set.smb.byn = (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.boc.mb[3]))) ? SV_BOUNDARY_PERIODIC : SV_BOUNDARY_NONE;
        xgc->data.set.smb.bynpos = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.mbpos[3]));

        xgc->data.set.smb.bz0 = (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.boc.mb[4]))) ? SV_BOUNDARY_PERIODIC : SV_BOUNDARY_NONE;
        xgc->data.set.smb.bz0pos = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.mbpos[4]));

        xgc->data.set.smb.bzn = (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.boc.mb[5]))) ? SV_BOUNDARY_PERIODIC : SV_BOUNDARY_NONE;
        xgc->data.set.smb.bznpos = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(xgc->cs.boc.mbpos[5]));
    } else if (xgc->par_selid == SET_MEDIUM) {
        /* V. 1. Media */

        /* V. 1. Material properties */
        xgc->data.set.sm.in_voxel_filename = g_strdup(gtk_entry_get_text(GTK_ENTRY(xgc->cs.mc.material_voxel_filename)));
        xgc->data.set.sm.in_vector_filename = g_strdup(gtk_entry_get_text(GTK_ENTRY(xgc->cs.mc.material_vector_filename)));
        //xgc->data.set.sm.in_vtk_data = g_strdup(gtk_entry_get_text(GTK_ENTRY(xgc->cs.mc.material_vtk_filename)));
        xgc->data.set.sm.matmode_check = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.mc.material_modecheck));
        xgc->data.set.sm.smooth = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.material_smoothsteps_spin));
    } else if (xgc->par_selid >= SET_GROW && xgc->par_selid < SET_ROUGHNESS) {
        /* V. 2. Add growth modifier */
        id = xgc->par_selid - SET_GROW;

        xgc->data.set.sm.grow_i0[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_i0_spin));
        xgc->data.set.sm.grow_j0[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_j0_spin));
        xgc->data.set.sm.grow_k0[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_k0_spin));
        xgc->data.set.sm.grow_in[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_in_spin));
        xgc->data.set.sm.grow_jn[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_jn_spin));
        xgc->data.set.sm.grow_kn[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_kn_spin));

        xgc->data.set.sm.grow_skipi0[id] = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.mc.grow_skipi0));
        xgc->data.set.sm.grow_skipj0[id] = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.mc.grow_skipj0));
        xgc->data.set.sm.grow_skipk0[id] = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.mc.grow_skipk0));
        xgc->data.set.sm.grow_skipin[id] = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.mc.grow_skipin));
        xgc->data.set.sm.grow_skipjn[id] = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.mc.grow_skipjn));
        xgc->data.set.sm.grow_skipkn[id] = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.mc.grow_skipkn));

        xgc->data.set.sm.grow_addindex[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_addindex_spin));
        xgc->data.set.sm.grow_attachindex[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_attachindex_spin));
        xgc->data.set.sm.grow_subsampling[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_subsampling_spin));
        xgc->data.set.sm.grow_seed[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_seed_spin));
        xgc->data.set.sm.grow_nsteps[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_nsteps_spin));
        xgc->data.set.sm.grow_mobility[id] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_mobility_spin));
        xgc->data.set.sm.grow_probability[id] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_probability_spin));
    } else if (xgc->par_selid >= SET_ROUGHNESS && xgc->par_selid < SET_SPECTRAL) {
        /* V. 3. Add roughness modifier */
        id = xgc->par_selid - SET_ROUGHNESS;

        xgc->data.set.sm.rough_matindex[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.rough_matindex_spin));
        xgc->data.set.sm.rough_voidindex[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.rough_voidindex_spin));
        xgc->data.set.sm.rough_iterations[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.rough_iterations_spin));
        xgc->data.set.sm.rough_seed[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.rough_seed_spin));
        xgc->data.set.sm.rough_radius_peak[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.rough_radiuspeak_spin));
        xgc->data.set.sm.rough_radius_span[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.rough_radiusspan_spin));
        xgc->data.set.sm.rough_probability[id] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.rough_probability_spin));
    } else if (xgc->par_selid >= SET_SPECTRAL && xgc->par_selid < SET_EXPRESSION) {
        /* V. 4. Add spectral modifier */
        id = xgc->par_selid - SET_SPECTRAL;

        xgc->data.set.sm.spectral_sigma[id] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.spectral_sigma_spin));
        xgc->data.set.sm.spectral_t[id] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.spectral_t_spin));
        xgc->data.set.sm.spectral_matindex[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.spectral_matindex_spin));
        xgc->data.set.sm.spectral_seed[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.spectral_seed_spin));
    } else if (xgc->par_selid >= SET_EXPRESSION && xgc->par_selid < SET_NFAREA) {
        /* V. 5. Add expression modifier */
        id = xgc->par_selid - SET_EXPRESSION;

        xgc->data.set.sm.expr_i0[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_i0_spin));
        xgc->data.set.sm.expr_j0[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_j0_spin));
        xgc->data.set.sm.expr_k0[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_k0_spin));
        xgc->data.set.sm.expr_in[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_in_spin));
        xgc->data.set.sm.expr_jn[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_jn_spin));
        xgc->data.set.sm.expr_kn[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_kn_spin));
        xgc->data.set.sm.expr_matindex[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_matindex_spin));
        xgc->data.set.sm.expr_voidindex[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_voidindex_spin));
        xgc->data.set.sm.expr_maxdist[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_maxdist_spin));
        xgc->data.set.sm.expr_distmode[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_distmode_spin));        
        g_snprintf(xgc->data.set.sm.expr_expr[id], sizeof(xgc->data.set.sm.expr_expr), "%s", gtk_entry_get_text(GTK_ENTRY(xgc->cs.mc.expr_expr)));
    } else if (xgc->par_selid >= SET_POUT && xgc->par_selid < SET_IOUT) {
        /* VI. 2. Point output properties */
        id = xgc->par_selid - SET_POUT;
        xgc->data.set.so.pnts[id].i = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.point_output_origin_position_i_spin));
        xgc->data.set.so.pnts[id].j = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.point_output_origin_position_j_spin));
        xgc->data.set.so.pnts[id].k = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.point_output_origin_position_k_spin));
        xgc->data.set.so.pnts[id].component = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.oc.point_output_component));
        xgc->data.set.so.pnts[id].step = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.point_output_step_spin));
        g_snprintf(xgc->data.set.so.pnts[id].filebase, sizeof(xgc->data.set.so.pnts[id].filebase), "%s", gtk_entry_get_text(GTK_ENTRY(xgc->cs.oc.point_output_filename)));
    } else if (xgc->par_selid >= SET_IIOUT && xgc->par_selid < SET_COUT) {
        /* VI. 3. Image output properties */
        id = xgc->par_selid - SET_IIOUT;
        gint i0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.image_output_pos_spin));

        gint dirval = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.oc.image_output_plane));
        if (dirval == 0) {
            xgc->data.set.so.imgs[id].i = i0;
            xgc->data.set.so.imgs[id].j = -1;
            xgc->data.set.so.imgs[id].k = -1;
        } else if (dirval == 1) {
            xgc->data.set.so.imgs[id].j = i0;
            xgc->data.set.so.imgs[id].i = -1;
            xgc->data.set.so.imgs[id].k = -1;
        } else if (dirval == 2) {
            xgc->data.set.so.imgs[id].k = i0;
            xgc->data.set.so.imgs[id].i = -1;
            xgc->data.set.so.imgs[id].j = -1;
        }
        xgc->data.set.so.imgs[id].component = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.oc.image_output_component));
        xgc->data.set.so.imgs[id].step = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.image_output_step_spin));
        g_snprintf(xgc->data.set.so.imgs[id].filebase, sizeof(xgc->data.set.so.imgs[id].filebase), "%s", gtk_entry_get_text(GTK_ENTRY(xgc->cs.oc.image_output_label)));
    } else if (xgc->par_selid >= SET_IOUT && xgc->par_selid < SET_IIOUT) {
        /* VI. 4. Plane output properties */
        id = xgc->par_selid - SET_IOUT;
        gint i0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.plane_output_pos_spin));

        gint dirval = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.oc.plane_output_plane));
        if (dirval == 0) {
            xgc->data.set.so.plns[id].i = i0;
            xgc->data.set.so.plns[id].j = -1;
            xgc->data.set.so.plns[id].k = -1;
        } else if (dirval == 1) {
            xgc->data.set.so.plns[id].j = i0;
            xgc->data.set.so.plns[id].i = -1;
            xgc->data.set.so.plns[id].k = -1;
        } else if (dirval == 2) {
            xgc->data.set.so.plns[id].k = i0;
            xgc->data.set.so.plns[id].i = -1;
            xgc->data.set.so.plns[id].j = -1;
        }
        xgc->data.set.so.plns[id].component = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.oc.plane_output_component));
        xgc->data.set.so.plns[id].step = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.plane_output_step_spin));
        g_snprintf(xgc->data.set.so.plns[id].filebase, sizeof(xgc->data.set.so.plns[id].filebase), "%s", gtk_entry_get_text(GTK_ENTRY(xgc->cs.oc.plane_output_filename)));
        xgc->data.set.so.plns[id].format = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.oc.plane_output_format));
        xgc->data.set.so.plns[id].start = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.plane_output_start_spin));
        xgc->data.set.so.plns[id].stop = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.plane_output_stop_spin));
    } else if (xgc->par_selid >= SET_COUT && xgc->par_selid < SET_SOUT) {
        /* VI. 5. Volume output properties */
        id = xgc->par_selid - SET_COUT;

        xgc->data.set.so.cubs[id].component = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.oc.volume_output_component));
        xgc->data.set.so.cubs[id].step = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.volume_output_step_spin));
        g_snprintf(xgc->data.set.so.cubs[id].filebase, sizeof(xgc->data.set.so.cubs[id].filebase), "%s", gtk_entry_get_text(GTK_ENTRY(xgc->cs.oc.volume_output_filename)));
        xgc->data.set.so.cubs[id].format = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.oc.volume_output_format));
        xgc->data.set.so.cubs[id].start = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.volume_output_start_spin));
        xgc->data.set.so.cubs[id].stop = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.volume_output_stop_spin));
    } else if (xgc->par_selid >= SET_SOUT && xgc->par_selid < SET_FOUT) {
        /* VI. 5. Sum output properties */
        id = xgc->par_selid - SET_SOUT;

        xgc->data.set.so.sums[id].box_i0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_i0_spin));
        xgc->data.set.so.sums[id].box_j0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_j0_spin));
        xgc->data.set.so.sums[id].box_k0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_k0_spin));
        xgc->data.set.so.sums[id].box_in = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_in_spin));
        xgc->data.set.so.sums[id].box_jn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_jn_spin));
        xgc->data.set.so.sums[id].box_kn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_kn_spin));

        xgc->data.set.so.sums[id].component = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->cs.oc.sum_output_component));
        xgc->data.set.so.sums[id].step = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.sum_output_step_spin));
        g_snprintf(xgc->data.set.so.sums[id].filename, sizeof(xgc->data.set.so.sums[id].filename), "%s", gtk_entry_get_text(GTK_ENTRY(xgc->cs.oc.sum_output_filename)));

        xgc->data.set.so.sums[id].layered_epsilon = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.sum_output_epsilon_spin));
        xgc->data.set.so.sums[id].layered_mu = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.sum_output_mu_spin));
        xgc->data.set.so.sums[id].layered_sigma = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.sum_output_sigma_spin));
        xgc->data.set.so.sums[id].layered_sigast = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.sum_output_sigast_spin));

        //set->so.sums[set->so.nsums-1].stringbased = 0;
        //memset(set->so.sums[set->so.nsums-1].string, 0, sizeof(set->so.sums[set->so.nsums-1].string));        
    } else if (xgc->par_selid >= SET_FOUT && xgc->par_selid < SET_GROW) {
        /* VI. 6. Force output properties */
        id = xgc->par_selid - SET_FOUT;

        xgc->data.set.so.forces[id].box_i0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_i0_spin));
        xgc->data.set.so.forces[id].box_j0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_j0_spin));
        xgc->data.set.so.forces[id].box_k0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_k0_spin));
        xgc->data.set.so.forces[id].box_in = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_in_spin));
        xgc->data.set.so.forces[id].box_jn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_jn_spin));
        xgc->data.set.so.forces[id].box_kn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_kn_spin));

        xgc->data.set.so.forces[id].step = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.oc.force_output_step_spin));
        g_snprintf(xgc->data.set.so.forces[id].filename, sizeof(xgc->data.set.so.forces[id].filename), "%s", gtk_entry_get_text(GTK_ENTRY(xgc->cs.oc.force_output_filename)));
    } else if (xgc->par_selid == SET_NFFF) {
        /* VII. 1. Near field to far field transform box */

        xgc->data.set.sf.box_i0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_box_i0_spin));
        xgc->data.set.sf.box_j0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_box_j0_spin));
        xgc->data.set.sf.box_k0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_box_k0_spin));
        xgc->data.set.sf.box_in = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_box_in_spin));
        xgc->data.set.sf.box_jn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_box_jn_spin));
        xgc->data.set.sf.box_kn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_box_kn_spin));

        xgc->data.set.sf.box_boundary_skipi0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.nfff_box_boundary_skipi0));
        xgc->data.set.sf.skipi0_jmin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipi0_jmin_spin));
        xgc->data.set.sf.skipi0_kmin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipi0_kmin_spin));
        xgc->data.set.sf.skipi0_jmax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipi0_jmax_spin));
        xgc->data.set.sf.skipi0_kmax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipi0_kmax_spin));

        xgc->data.set.sf.box_boundary_skipin = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.nfff_box_boundary_skipin));
        xgc->data.set.sf.skipin_jmin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipin_jmin_spin));
        xgc->data.set.sf.skipin_kmin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipin_kmin_spin));
        xgc->data.set.sf.skipin_jmax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipin_jmax_spin));
        xgc->data.set.sf.skipin_kmax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipin_kmax_spin));

        xgc->data.set.sf.box_boundary_skipj0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.nfff_box_boundary_skipj0));
        xgc->data.set.sf.skipj0_imin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipj0_imin_spin));
        xgc->data.set.sf.skipj0_kmin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipj0_kmin_spin));
        xgc->data.set.sf.skipj0_imax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipj0_imax_spin));
        xgc->data.set.sf.skipj0_kmax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipj0_kmax_spin));

        xgc->data.set.sf.box_boundary_skipjn = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.nfff_box_boundary_skipjn));
        xgc->data.set.sf.skipjn_imin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipjn_imin_spin));
        xgc->data.set.sf.skipjn_kmin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipjn_kmin_spin));
        xgc->data.set.sf.skipjn_imax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipjn_imax_spin));
        xgc->data.set.sf.skipjn_kmax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipjn_kmax_spin));

        xgc->data.set.sf.box_boundary_skipk0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.nfff_box_boundary_skipk0));
        xgc->data.set.sf.skipk0_imin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipk0_imin_spin));
        xgc->data.set.sf.skipk0_jmin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipk0_jmin_spin));
        xgc->data.set.sf.skipk0_imax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipk0_imax_spin));
        xgc->data.set.sf.skipk0_jmax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipk0_jmax_spin));

        xgc->data.set.sf.box_boundary_skipkn = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.nfff_box_boundary_skipkn));
        xgc->data.set.sf.skipkn_imin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipkn_imin_spin));
        xgc->data.set.sf.skipkn_jmin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipkn_jmin_spin));
        xgc->data.set.sf.skipkn_imax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipkn_imax_spin));
        xgc->data.set.sf.skipkn_jmax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipkn_jmax_spin));
    } else if (xgc->par_selid >= SET_NFFFP && xgc->par_selid < SET_PNFAREA) {
        /* VII. 2. Near field to far field point */
        id = xgc->par_selid - SET_NFFFP;

        xgc->data.set.sf.ri[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffp_pos_i_spin));
        xgc->data.set.sf.rj[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffp_pos_j_spin));
        xgc->data.set.sf.rk[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffp_pos_k_spin));
        xgc->data.set.sf.source_filename[id] = g_strdup(gtk_entry_get_text(GTK_ENTRY(xgc->cs.fc.nfffp_filename)));
        xgc->data.set.sf.individual[id] = 0;
    } else if (xgc->par_selid >= SET_NFAREA && xgc->par_selid < SET_NFFFP) {
        /* VII. 3. Near field to far field area */
        id = xgc->par_selid - SET_NFAREA;

        xgc->data.set.sf.area_thetares[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_thetares_spin));
        xgc->data.set.sf.area_phires[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_phires_spin));
        xgc->data.set.sf.area_radius[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_radius_spin));
        xgc->data.set.sf.area_thetafrom[id] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_thetafrom_spin))*G_PI / 180;
        xgc->data.set.sf.area_thetato[id] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_thetato_spin))*G_PI / 180;
        xgc->data.set.sf.area_phifrom[id] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_phifrom_spin))*G_PI / 180;
        xgc->data.set.sf.area_phito[id] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_phito_spin))*G_PI / 180;
        xgc->data.set.sf.area_savefile[id] = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.nfffa_savefile));
    } else if (xgc->par_selid == SET_PNFFF) {
        /* VIII. 1. Periodic near field to far field transform box */

        xgc->data.set.spf.box_i0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_box_i0_spin));
        xgc->data.set.spf.box_j0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_box_j0_spin));
        xgc->data.set.spf.box_k0 = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_box_k0_spin));
        xgc->data.set.spf.box_in = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_box_in_spin));
        xgc->data.set.spf.box_jn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_box_jn_spin));
        xgc->data.set.spf.box_kn = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_box_kn_spin));

        xgc->data.set.spf.pimin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_integration_xmin_spin));
        xgc->data.set.spf.pimax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_integration_xmax_spin));
        xgc->data.set.spf.pjmin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_integration_ymin_spin));
        xgc->data.set.spf.pjmax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_integration_ymax_spin));

        xgc->data.set.spf.box_boundary_skipk0 = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.pnfff_box_boundary_skipk0));
        xgc->data.set.spf.skipk0_imin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipk0_imin_spin));
        xgc->data.set.spf.skipk0_jmin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipk0_jmin_spin));
        xgc->data.set.spf.skipk0_imax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipk0_imax_spin));
        xgc->data.set.spf.skipk0_jmax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipk0_jmax_spin));

        xgc->data.set.spf.box_boundary_skipkn = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.pnfff_box_boundary_skipkn));
        xgc->data.set.spf.skipkn_imin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipkn_imin_spin));
        xgc->data.set.spf.skipkn_jmin = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipkn_jmin_spin));
        xgc->data.set.spf.skipkn_imax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipkn_imax_spin));
        xgc->data.set.spf.skipkn_jmax = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipkn_jmax_spin));
    } else if (xgc->par_selid >= SET_PNFAREA && xgc->par_selid < SET_PNFFFP) {
        /* VIII. 3. Periodic near field to far field area */

        id = xgc->par_selid - SET_PNFAREA;

        xgc->data.set.spf.area_thetares[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_thetares_spin));
        xgc->data.set.spf.area_phires[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_phires_spin));
        xgc->data.set.spf.area_radius[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_radius_spin));
        xgc->data.set.spf.area_thetafrom[id] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_thetafrom_spin))*G_PI / 180;
        xgc->data.set.spf.area_thetato[id] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_thetato_spin))*G_PI / 180;
        xgc->data.set.spf.area_phifrom[id] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_phifrom_spin))*G_PI / 180;
        xgc->data.set.spf.area_phito[id] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_phito_spin))*G_PI / 180;
        xgc->data.set.spf.area_savefile[id] = (gint)gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.nfffa_savefile));
    } else if (xgc->par_selid >= SET_PNFFFP) {
        /* VIII. 2. Periodic near field to far field point */

        id = xgc->par_selid - SET_PNFFFP;

        xgc->data.set.spf.ri[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffp_pos_i_spin));
        xgc->data.set.spf.rj[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffp_pos_j_spin));
        xgc->data.set.spf.rk[id] = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffp_pos_k_spin));
        xgc->data.set.spf.source_filename[id] = g_strdup(gtk_entry_get_text(GTK_ENTRY(xgc->cs.fc.nfffp_filename)));
        xgc->data.set.spf.individual[id] = 0;
    }

    //pt_remeber_expanded_rows_and_selection(xgc);
    //write_parfile(xgc);
    //pt_restore_expanded_rows_and_selection(xgc);


    par_controls_sensitive(xgc, xgc->par_selid);

    GtkTreeIter iter_parent, iter;
    GtkTreePath *path;
    gint leaf_pos = 0;

    if (xgc->par_selid == SET_POOL) {
        /* I. 1. Computational domain */
        path = gtk_tree_path_new_from_indices(ROW_POOL_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_POOL, 0);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 0);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_POOL, 1);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 1);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        //gtk_tree_store_set(xgc->ts_par, &xgc->par_seliter, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_pool, COLUMN_ID, SET_POOL, COLUMN_SHOW_TOGGLE, TRUE, -1);
    } else if (xgc->par_selid == SET_BASIC) {
        /* II. 1. Basic parameters */
        path = gtk_tree_path_new_from_indices(ROW_BASIC_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BASIC, 0);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 0);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BASIC, 1);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 1);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BASIC, 2);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 2);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BASIC, 3);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 3);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BASIC, 4);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 4);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid >= SET_PSOURCE && xgc->par_selid < SET_POUT) {
        id = xgc->par_selid - SET_PSOURCE;

        /* III. 1. Point source properties */
        path = gtk_tree_path_new_from_indices(ROW_SOURCES_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_PSOURCE, id);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, id);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid == SET_SF) {
        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_SF, -1);
        gtk_tree_store_set(xgc->ts_par, &xgc->par_seliter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid == SET_TSF) {
        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_TSF, -1);
        gtk_tree_store_set(xgc->ts_par, &xgc->par_seliter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid == SET_TSFF) {
        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_TSFF, -1);
        gtk_tree_store_set(xgc->ts_par, &xgc->par_seliter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid == SET_LTSF) {
        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_LTSF, -1);
        gtk_tree_store_set(xgc->ts_par, &xgc->par_seliter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid == SET_LTSFF) {
        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_LTSFF, -1);
        gtk_tree_store_set(xgc->ts_par, &xgc->par_seliter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid == SET_BND) {
        /* IV. 1. Boundary conditions */
        path = gtk_tree_path_new_from_indices(ROW_BOUNDARY_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 0);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 0);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 1);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 1);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 2);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 2);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 3);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 3);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 4);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 4);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 5);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 5);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        /* optional periodic boundaries */
        if (xgc->data.set.smb.bx0 == SV_BOUNDARY_PERIODIC) {
            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 6);
            if (gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 6))
                gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
            else {
                gtk_tree_store_append(xgc->ts_par, &iter, &iter_parent);
                gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
        }

        if (xgc->data.set.smb.bxn == SV_BOUNDARY_PERIODIC) {
            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 7);
            if (gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 7))
                gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
            else {
                gtk_tree_store_append(xgc->ts_par, &iter, &iter_parent);
                gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
        }

        if (xgc->data.set.smb.by0 == SV_BOUNDARY_PERIODIC) {
            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 8);
            if (gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 8))
                gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
            else {
                gtk_tree_store_append(xgc->ts_par, &iter, &iter_parent);
                gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
        }

        if (xgc->data.set.smb.byn == SV_BOUNDARY_PERIODIC) {
            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 9);
            if (gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 9))
                gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
            else {
                gtk_tree_store_append(xgc->ts_par, &iter, &iter_parent);
                gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
        }

        if (xgc->data.set.smb.bz0 == SV_BOUNDARY_PERIODIC) {
            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 10);
            if (gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 10))
                gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
            else {
                gtk_tree_store_append(xgc->ts_par, &iter, &iter_parent);
                gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
        }

        if (xgc->data.set.smb.bzn == SV_BOUNDARY_PERIODIC) {
            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 11);
            if (gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 11))
                gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
            else {
                gtk_tree_store_append(xgc->ts_par, &iter, &iter_parent);
                gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
        }

        //assemble_row_text(xgc, buff, sizeof(buff), SET_BOUNDS, id);
        //gtk_tree_store_set(xgc->ts_par, &xgc->par_seliter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid == SET_MEDIUM) {
        /* V. 1. Media */

        /* V. 1. Material properties */
        path = gtk_tree_path_new_from_indices(ROW_MEDIA_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_MEDIUM, 0);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 0);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_MEDIUM, 1);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 1);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_MEDIUM, 2);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 2);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_MEDIUM, 3);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 3);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid >= SET_GROW && xgc->par_selid < SET_ROUGHNESS) {
        /* V. 2. Add growth modifier */
        id = xgc->par_selid - SET_GROW;
        leaf_pos = id + NOL_MEDIUM_MAT_PROPS;   /* voxel, vector, material mode and medium smoothing */

        path = gtk_tree_path_new_from_indices(ROW_MEDIA_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_GROW, leaf_pos);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, leaf_pos);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid >= SET_ROUGHNESS && xgc->par_selid < SET_SPECTRAL) {
        /* V. 3. Add roughness modifier */
        id = xgc->par_selid - SET_ROUGHNESS;
        leaf_pos = id + NOL_MEDIUM_MAT_PROPS + xgc->data.set.sm.ngrowths;   /* voxel, vector, material mode and medium smoothing */

        path = gtk_tree_path_new_from_indices(ROW_MEDIA_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_ROUGHNESS, leaf_pos);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, leaf_pos);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid >= SET_SPECTRAL && xgc->par_selid < SET_EXPRESSION) {
        /* V. 4. Add spectral modifier */
    } else if (xgc->par_selid >= SET_EXPRESSION && xgc->par_selid < SET_NFAREA) {
        /* V. 5. Add expression modifier */
    } else if (xgc->par_selid == SET_OUT) {
        /* VI. 1. General output */
        path = gtk_tree_path_new_from_indices(ROW_OUTPUTS_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_OUT, 0);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 0);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid >= SET_POUT && xgc->par_selid < SET_IOUT) {
        /* VI. 2. Point output properties */
        id = xgc->par_selid - SET_POUT;
        leaf_pos = id + NOL_OUTPUT_PROPS;

        path = gtk_tree_path_new_from_indices(ROW_OUTPUTS_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_POUT, leaf_pos);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, leaf_pos);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        /*if (xgc->nfilestoshow < 99) {
            xgc->filestoshow[xgc->nfilestoshow] = g_strdup(xgc->data.set.so.pnts[id].filebase);
            xgc->formattoshow[xgc->nfilestoshow] = xgc->data.set.so.pnts[id].component;
            xgc->nfilestoshow += 1;
        }*/
    } else if (xgc->par_selid >= SET_IIOUT && xgc->par_selid < SET_COUT) {
        /* VI. 3. Image output properties */
        id = xgc->par_selid - SET_IIOUT;
        leaf_pos = id + NOL_OUTPUT_PROPS + xgc->data.set.so.npnts;

        path = gtk_tree_path_new_from_indices(ROW_OUTPUTS_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_IIOUT, leaf_pos);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, leaf_pos);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        //////////////////////////////////////////////////////////////////////////
        /*gint active = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->image_combo));
        if (active >= 0) {
            g_snprintf(buff, sizeof(buff), "%s component, %s, %s",
                       cpstring(xgc->data.set.so.imgs[i].component),
                       planestring(xgc->data.set.so.imgs[i].i, xgc->data.set.so.imgs[i].j, xgc->data.set.so.imgs[i].k), xgc->data.set.so.imgs[i].filebase);
            xgc->imagestoshow[xgc->nimagestoshow] = g_strdup(buff);

            gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(xgc->image_combo), active);
            gtk_combo_box_text_insert_text(GTK_COMBO_BOX_TEXT(xgc->image_combo), active, xgc->imagestoshow[i]);
        }*/

        xgc->nimagestoshow = 0;
        memset(xgc->imagestoshow, 0, IMAGES_TO_SHOW * sizeof(gchar*));
        for (i = 0; i < MIN(xgc->data.set.so.nimgs, TREE_MAXENT); i++) {
            if (xgc->nimagestoshow < 99) {
                g_snprintf(buff, sizeof(buff), "%s component, %s, %s",
                           cpstring(xgc->data.set.so.imgs[i].component),
                           planestring(xgc->data.set.so.imgs[i].i, xgc->data.set.so.imgs[i].j, xgc->data.set.so.imgs[i].k), xgc->data.set.so.imgs[i].filebase);
                xgc->imagestoshow[xgc->nimagestoshow] = g_strdup(buff);
                xgc->nimagestoshow += 1;
            }
        }
        
        for (i = 0; i < 100; i++)
            //gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(xgc->image_combo), 0);
            gtk_combo_box_remove_text(GTK_COMBO_BOX(xgc->image_combo), 0);
        for (i = 0; i < xgc->nimagestoshow; i++) {
            if (xgc->imagestoshow[i] != NULL)
                //gtk_combo_box_text_append_text(GTK_COMBO_BOX(xgc->image_combo), xgc->imagestoshow[i]);
                gtk_combo_box_append_text(GTK_COMBO_BOX(xgc->image_combo), xgc->imagestoshow[i]);
        }

        gint active = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->image_combo));
        if (active >= 0)
            gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->image_combo), active);
        else if (xgc->nimagestoshow > 0)
            gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->image_combo), 0);
        //////////////////////////////////////////////////////////////////////////
    } else if (xgc->par_selid >= SET_IOUT && xgc->par_selid < SET_IIOUT) {
        /* VI. 4. Plane output properties */
        id = xgc->par_selid - SET_IOUT;
        leaf_pos = id + NOL_OUTPUT_PROPS + xgc->data.set.so.npnts + xgc->data.set.so.nimgs;

        path = gtk_tree_path_new_from_indices(ROW_OUTPUTS_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_IOUT, leaf_pos);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, leaf_pos);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid >= SET_COUT && xgc->par_selid < SET_SOUT) {
        /* VI. 5. Volume output properties */
        id = xgc->par_selid - SET_COUT;
        leaf_pos = id + NOL_OUTPUT_PROPS + xgc->data.set.so.npnts + xgc->data.set.so.nimgs + xgc->data.set.so.nplns;

        path = gtk_tree_path_new_from_indices(ROW_OUTPUTS_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_COUT, leaf_pos);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, leaf_pos);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid >= SET_SOUT && xgc->par_selid < SET_FOUT) {
        /* VI. 5. Sum output properties */
        id = xgc->par_selid - SET_SOUT;
        leaf_pos = id + NOL_OUTPUT_PROPS + xgc->data.set.so.npnts + xgc->data.set.so.nimgs + xgc->data.set.so.nplns + xgc->data.set.so.ncubs;

        path = gtk_tree_path_new_from_indices(ROW_OUTPUTS_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_SOUT, leaf_pos);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, leaf_pos);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid >= SET_FOUT && xgc->par_selid < SET_GROW) {
        /* VI. 6. Force output properties */
        id = xgc->par_selid - SET_FOUT;
        leaf_pos = id + NOL_OUTPUT_PROPS + xgc->data.set.so.npnts + xgc->data.set.so.nimgs + xgc->data.set.so.nplns + xgc->data.set.so.ncubs + xgc->data.set.so.nsums;

        path = gtk_tree_path_new_from_indices(ROW_OUTPUTS_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_FOUT, leaf_pos);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, leaf_pos);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid == SET_NFFF) {
        /* VII. 1. Near field to far field transform box */

        path = gtk_tree_path_new_from_indices(ROW_NFFF_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFF, 0);    /* NFFF box */
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 0);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFF, 1);    /* NFFF skips */
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 1);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFF, 2);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 2);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFF, 3);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 3);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFF, 4);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 4);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFF, 5);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 5);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFF, 6);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 6);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid >= SET_NFFFP && xgc->par_selid < SET_PNFAREA) {
        /* VII. 2. Near field to far field point */
        id = xgc->par_selid - SET_NFFFP;
        leaf_pos = id + NOL_NFFF_PROPS;

        path = gtk_tree_path_new_from_indices(ROW_NFFF_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFFP, leaf_pos);    /* Far field point */
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, leaf_pos);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid >= SET_NFAREA && xgc->par_selid < SET_NFFFP) {
        /* VII. 3. Near field to far field area */
        id = xgc->par_selid - SET_NFAREA;
        leaf_pos = id + NOL_NFFF_PROPS + xgc->data.set.sf.nrs;

        path = gtk_tree_path_new_from_indices(ROW_NFFF_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFAREA, leaf_pos);   /* Far field area */
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, leaf_pos);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid == SET_PNFFF) {
        /* VIII. 1. Periodic near field to far field transform box */

        path = gtk_tree_path_new_from_indices(ROW_PNFFF_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_PNFFF, 0);    /* PNFFF box */
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 0);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_PNFFF, 1);    /* PNFFF skips */
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 1);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_PNFFF, 2);
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, 2);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid >= SET_PNFFFP) {
        /* VIII. 2. Periodic near field to far field point */

        id = xgc->par_selid - SET_PNFFFP;
        leaf_pos = id + NOL_PNFFF_PROPS;

        path = gtk_tree_path_new_from_indices(ROW_PNFFF_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_PNFFFP, leaf_pos);    /* Periodic Far field point */
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, leaf_pos);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    } else if (xgc->par_selid >= SET_PNFAREA && xgc->par_selid < SET_PNFFFP) {
        /* VIII. 3. Periodic near field to far field area */

        id = xgc->par_selid - SET_PNFAREA;
        leaf_pos = id + NOL_PNFFF_PROPS + xgc->data.set.spf.nrs;

        path = gtk_tree_path_new_from_indices(ROW_PNFFF_INDEX, -1);
        gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_par), &iter_parent, path);

        pt_assemble_row_text(xgc, buff, sizeof(buff), SET_PNFAREA, leaf_pos);   /* Periodic Far field area */
        gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_par), &iter, &iter_parent, leaf_pos);
        gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, buff, -1);
    }

    gtk_widget_queue_draw(xgc->view_scene);

    // remember focus at the beginning and restore at the end of ss_controls_changed()
    // Note: this works only for GtkEntry not for GtkSpinButton
    if (NULL != focus)
        gtk_window_set_focus(GTK_WINDOW(xgc->toplevel), focus);
} /* par_controls_changed() */

static gboolean
set_values_general(XGControls* xgc, SvSet *set, SvType type, gint id)
{
    gint comp_domain_size_xval = COMP_DOMAIN_SIZE, comp_domain_size_yval = COMP_DOMAIN_SIZE, comp_domain_size_zval = COMP_DOMAIN_SIZE;
    gdouble comp_domain_spacing_xval = COMP_DOMAIN_SPACING, comp_domain_spacing_yval = COMP_DOMAIN_SPACING, comp_domain_spacing_zval = COMP_DOMAIN_SPACING;
    gint basic_nstepsval = BASIC_STEPS, basic_verboseval = BASIC_VERBOSE, basic_nthreadsval = BASIC_NTHREADS, basic_ugval = BASIC_USEGPU;
    gboolean basic_usegpu = BASIC_USEGPU;

    /* computational volume boundary properties */
    gint cvboundary_typeval[6] = {CVBOUNDARY_TYPE, CVBOUNDARY_TYPE, CVBOUNDARY_TYPE, CVBOUNDARY_TYPE, CVBOUNDARY_TYPE, CVBOUNDARY_TYPE};
    gint cvboundary_depthval[6] = {CVBOUNDARY_DEPTH, CVBOUNDARY_DEPTH, CVBOUNDARY_DEPTH, CVBOUNDARY_DEPTH, CVBOUNDARY_DEPTH, CVBOUNDARY_DEPTH};
    gint cvboundary_mval[6] = {CVBOUNDARY_M, CVBOUNDARY_M, CVBOUNDARY_M, CVBOUNDARY_M, CVBOUNDARY_M, CVBOUNDARY_M};
    gdouble cvboundary_sigmaval[6] = {CVBOUNDARY_SIGMA, CVBOUNDARY_SIGMA, CVBOUNDARY_SIGMA, CVBOUNDARY_SIGMA, CVBOUNDARY_SIGMA, CVBOUNDARY_SIGMA};
    gdouble cvboundary_aval[6] = {CVBOUNDARY_A, CVBOUNDARY_A, CVBOUNDARY_A, CVBOUNDARY_A, CVBOUNDARY_A, CVBOUNDARY_A};
    gdouble cvboundary_kappaval[6] = {CVBOUNDARY_KAPPA, CVBOUNDARY_KAPPA, CVBOUNDARY_KAPPA, CVBOUNDARY_KAPPA, CVBOUNDARY_KAPPA, CVBOUNDARY_KAPPA};
    gboolean pboundary_mbval[6] = {CVBOUNDARY_MB, CVBOUNDARY_MB, CVBOUNDARY_MB, CVBOUNDARY_MB, CVBOUNDARY_MB, CVBOUNDARY_MB};
    gint pboundary_mbposval[6] = {CVBOUNDARY_MBPOSX0, CVBOUNDARY_MBPOSXN, CVBOUNDARY_MBPOSY0, CVBOUNDARY_MBPOSYN, CVBOUNDARY_MBPOSZ0, CVBOUNDARY_MBPOSZN};

    gint i;
    gchar size_um_buff[256] = {0};

    if (type == SV_TYPE_POOL) {
        /* I. 1. Computational domain */
        comp_domain_size_xval = set->sp.xres;
        comp_domain_size_yval = set->sp.yres;
        comp_domain_size_zval = set->sp.zres;
        comp_domain_spacing_xval = set->sp.dx;
        comp_domain_spacing_yval = set->sp.dy;
        comp_domain_spacing_zval = set->sp.dz;
    } else if (type == SV_TYPE_BASIC) {
        /* II. 1. Basic parameters */
        basic_nstepsval = set->sc.nsteps;
        basic_verboseval = set->sc.verbose;
        basic_nthreadsval = set->sc.nthreads;
        basic_usegpu = set->sc.usegpu;

        basic_ugval = 0;
        if (set->sc.ugpu[1] == 1)
            basic_ugval = 1;
        else if (set->sc.ugpu[2] == 1)
            basic_ugval = 2;
        else if (set->sc.ugpu[3] == 1)
            basic_ugval = 3;
    } else if (type == SV_TYPE_BND) {
        cvboundary_typeval[0] = set->sb.bx0;
        cvboundary_typeval[1] = set->sb.bxn;
        cvboundary_typeval[2] = set->sb.by0;
        cvboundary_typeval[3] = set->sb.byn;
        cvboundary_typeval[4] = set->sb.bz0;
        cvboundary_typeval[5] = set->sb.bzn;

        cvboundary_depthval[0] = set->sb.depth_bx0;
        cvboundary_depthval[1] = set->sb.depth_bxn;
        cvboundary_depthval[2] = set->sb.depth_by0;
        cvboundary_depthval[3] = set->sb.depth_byn;
        cvboundary_depthval[4] = set->sb.depth_bz0;
        cvboundary_depthval[5] = set->sb.depth_bzn;

        cvboundary_mval[0] = set->sb.m_bx0;
        cvboundary_mval[1] = set->sb.m_bxn;
        cvboundary_mval[2] = set->sb.m_by0;
        cvboundary_mval[3] = set->sb.m_byn;
        cvboundary_mval[4] = set->sb.m_bz0;
        cvboundary_mval[5] = set->sb.m_bzn;

        cvboundary_sigmaval[0] = set->sb.sigma_bx0;
        cvboundary_sigmaval[1] = set->sb.sigma_bxn;
        cvboundary_sigmaval[2] = set->sb.sigma_by0;
        cvboundary_sigmaval[3] = set->sb.sigma_byn;
        cvboundary_sigmaval[4] = set->sb.sigma_bz0;
        cvboundary_sigmaval[5] = set->sb.sigma_bzn;

        cvboundary_aval[0] = set->sb.a_bx0;
        cvboundary_aval[1] = set->sb.a_bxn;
        cvboundary_aval[2] = set->sb.a_by0;
        cvboundary_aval[3] = set->sb.a_byn;
        cvboundary_aval[4] = set->sb.a_bz0;
        cvboundary_aval[5] = set->sb.a_bzn;

        cvboundary_kappaval[0] = set->sb.kappa_bx0;
        cvboundary_kappaval[1] = set->sb.kappa_bxn;
        cvboundary_kappaval[2] = set->sb.kappa_by0;
        cvboundary_kappaval[3] = set->sb.kappa_byn;
        cvboundary_kappaval[4] = set->sb.kappa_bz0;
        cvboundary_kappaval[5] = set->sb.kappa_bzn;


        pboundary_mbval[0] = (set->smb.bx0 == SV_BOUNDARY_PERIODIC) ? TRUE : FALSE;
        pboundary_mbval[1] = (set->smb.bxn == SV_BOUNDARY_PERIODIC) ? TRUE : FALSE;
        pboundary_mbval[2] = (set->smb.by0 == SV_BOUNDARY_PERIODIC) ? TRUE : FALSE;
        pboundary_mbval[3] = (set->smb.byn == SV_BOUNDARY_PERIODIC) ? TRUE : FALSE;
        pboundary_mbval[4] = (set->smb.bz0 == SV_BOUNDARY_PERIODIC) ? TRUE : FALSE;
        pboundary_mbval[5] = (set->smb.bzn == SV_BOUNDARY_PERIODIC) ? TRUE : FALSE;

        pboundary_mbposval[0] = set->smb.bx0pos;
        pboundary_mbposval[1] = set->smb.bxnpos;
        pboundary_mbposval[2] = set->smb.by0pos;
        pboundary_mbposval[3] = set->smb.bynpos;
        pboundary_mbposval[4] = set->smb.bz0pos;
        pboundary_mbposval[5] = set->smb.bznpos;
    }


    xgc->par_lock_controls_update = TRUE;


    if (type == SV_TYPE_POOL) {
        /* I. 1. Computational domain */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.cc.comp_domain_size_x_spin), comp_domain_size_xval);              /* virtual size in voxels */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.cc.comp_domain_size_y_spin), comp_domain_size_yval);              /* virtual size in voxels */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.cc.comp_domain_size_z_spin), comp_domain_size_zval);              /* virtual size in voxels */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.cc.comp_domain_spacing_x_spin), comp_domain_spacing_xval*1e6);    /* virtual/real ratio in micrometers */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.cc.comp_domain_spacing_y_spin), comp_domain_spacing_yval*1e6);    /* virtual/real ratio in micrometers */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.cc.comp_domain_spacing_z_spin), comp_domain_spacing_zval*1e6);    /* virtual/real ratio in micrometers */

                                                                                                                            // real size in micrometers
        g_snprintf(size_um_buff, sizeof(size_um_buff), "%g", comp_domain_size_xval*comp_domain_spacing_xval*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.cc.comp_domain_size_x_um), size_um_buff);
        g_snprintf(size_um_buff, sizeof(size_um_buff), "%g", comp_domain_size_yval*comp_domain_spacing_yval*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.cc.comp_domain_size_y_um), size_um_buff);
        g_snprintf(size_um_buff, sizeof(size_um_buff), "%g", comp_domain_size_zval*comp_domain_spacing_zval*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.cc.comp_domain_size_z_um), size_um_buff);
    } else if (type == SV_TYPE_BASIC) {
        /* II. 1. Basic parameters */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.bac.basic_nsteps_spin), basic_nstepsval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.bac.basic_verbose_spin), basic_verboseval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.bac.basic_nthreads_spin), basic_nthreadsval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.bac.basic_ugpu_index_spin), basic_ugval);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.bac.basic_usegpu), basic_usegpu);
        //gtk_widget_set_sensitive(xgc->cs.bac.basic_ugpu_index_spin, xgc->data.set.sc.usegpu);
    } else if (type == SV_TYPE_BND) {
        /* IV. Boundary conditions */
        for (i = 0; i < 6; i++) {
            gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs.boc.btype[i]), cvboundary_typeval[i]);
            gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.boc.bdepthspin[i]), cvboundary_depthval[i]);
            gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.boc.bmspin[i]), cvboundary_mval[i]);
            gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.boc.bsigmaspin[i]), cvboundary_sigmaval[i]);
            gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.boc.baspin[i]), cvboundary_aval[i]);
            gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.boc.bkappaspin[i]), cvboundary_kappaval[i]);

            /* optional periodic boundaries */
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.boc.mb[i]), pboundary_mbval[i]);
            gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.boc.mbposspin[i]), pboundary_mbposval[i]);
        }
    }

    xgc->par_lock_controls_update = FALSE;

    return FALSE;
} /* set_values_general */

static gboolean
set_values_sf(XGControls* xgc, SvSet *set, SvType type, gint id)
{
    gchar *source_filenameval = NULL;

    gint box_i0val = SOURCE_BOX_VERTEX0, box_j0val = SOURCE_BOX_VERTEX0, box_k0val = SOURCE_BOX_VERTEX0;
    gint box_inval = SOURCE_BOX_VERTEXN, box_jnval = SOURCE_BOX_VERTEXN, box_knval = SOURCE_BOX_VERTEXN;
    gint source_point_origin_position_ival = SOURCE_POINT_ORIGIN_I, source_point_origin_position_jval = SOURCE_POINT_ORIGIN_J, source_point_origin_position_kval = SOURCE_POINT_ORIGIN_K;
    gdouble source_point_origin_thetaval = SOURCE_POINT_ORIGIN_THETA, source_point_origin_phival = SOURCE_POINT_ORIGIN_PHI;
    gdouble source_ia_thetadegval = SOURCE_IA_THETADEG, source_ia_phidegval = SOURCE_IA_PHIDEG, source_ia_psidegval = SOURCE_IA_PSIDEG, source_pulsewidthval = SOURCE_PULSE_WIDTH, source_wlval = SOURCE_WAVELENGTH, source_amplitudeval = SOURCE_AMPLITUDE;
    gdouble fs_thetamaxdegval = SOURCE_FS_THETAMAXDEG, fs_fdistval = SOURCE_FS_FDIST, fs_polarisationdegval = SOURCE_FS_POLARISATIONDEG;
    gint fs_nipval = SOURCE_FS_NIP, fs_mipval = SOURCE_FS_MIP;
    gint box_boundary_skipi0val = SOURCE_BOUNDARY_SKIP, box_boundary_skipinval = SOURCE_BOUNDARY_SKIP, box_boundary_skipj0val = SOURCE_BOUNDARY_SKIP, box_boundary_skipjnval = SOURCE_BOUNDARY_SKIP, box_boundary_skipk0val = SOURCE_BOUNDARY_SKIP, box_boundary_skipknval = SOURCE_BOUNDARY_SKIP;
    gint box_boundary_skipdepthval = SOURCE_BOUNDARY_SKIP_DEPTH;
    gint layered_zposval[NUM_LSMP];
    gdouble layered_epsval[NUM_LSMP], layered_sigmaval[NUM_LSMP], layered_muval[NUM_LSMP], layered_sigastval[NUM_LSMP];
    gchar* layered_materialval[NUM_LSMP];
    gint source_modeval = SOURCE_MODE;
    gint source_gaussian_mult_enable = SOURCE_GAUSSIAN_MULT_ENABLE;
    gdouble source_gaussian_mult_center_pos_i = xgc->data.set.sp.xres / 2.0f, source_gaussian_mult_center_pos_j = xgc->data.set.sp.yres / 2.0f;
    gdouble source_gaussian_mult_radius_i = SOURCE_XXX_MULT_RADIUS, source_gaussian_mult_radius_j = SOURCE_XXX_MULT_RADIUS;
    gint source_radial_mult_enable = SOURCE_RADIAL_MULT_ENABLE;
    gdouble source_radial_mult_center_pos_i = xgc->data.set.sp.xres / 2.0f, source_radial_mult_center_pos_j = xgc->data.set.sp.yres / 2.0f;
    gdouble source_radial_mult_radius_i = SOURCE_XXX_MULT_RADIUS, source_radial_mult_radius_j = SOURCE_XXX_MULT_RADIUS;
    gint source_fiber_mult_enable = SOURCE_FIBER_MULT_ENABLE;
    gdouble source_fiber_mult_center_pos_i = xgc->data.set.sp.xres / 2.0f, source_fiber_mult_center_pos_j = xgc->data.set.sp.yres / 2.0f;
    gdouble source_fiber_mult_radius = SOURCE_XXX_MULT_RADIUS;
    gdouble source_fiber_mult_cutoff = SOURCE_FIBER_MULT_CUTOFF;
    gdouble source_fiber_mult_epsilon_core = SOURCE_FIBER_MULT_EPSILON_CORE;
    gdouble source_fiber_mult_epsilon_cladding = SOURCE_FIBER_MULT_EPSILON_CLADDING;

    gint ilayers, i;

    gchar buff[256] = {0};

    memset(layered_zposval, SOURCE_LAYERED_ZPOS, NUM_LSMP * sizeof(gint));
    memset(layered_epsval, SOURCE_LAYERED_EPS, NUM_LSMP * sizeof(gdouble));
    memset(layered_sigmaval, SOURCE_LAYERED_SIGMA, NUM_LSMP * sizeof(gdouble));
    memset(layered_muval, SOURCE_LAYERED_MU, NUM_LSMP * sizeof(gdouble));
    memset(layered_sigastval, SOURCE_LAYERED_SIGAST, NUM_LSMP * sizeof(gdouble));
    for (i = 0; i < NUM_LSMP; i++)
        layered_materialval[i] = g_strdup(SOURCE_LAYERED_MATERIAL);

    if (type == SV_SRCTYPE_POINT) {
        //name = g_strdup("Point source");

        /* 0b. Point source origin properties */
        source_point_origin_position_ival = set->ss.pnts[id].point_origin_position_i;
        source_point_origin_position_jval = set->ss.pnts[id].point_origin_position_j;
        source_point_origin_position_kval = set->ss.pnts[id].point_origin_position_k;
        source_point_origin_thetaval = set->ss.pnts[id].point_origin_theta * 180 / G_PI;
        source_point_origin_phival = set->ss.pnts[id].point_origin_phi * 180 / G_PI;

        /* 1. Source properties */
        if (set->ss.pnts[id].source_filename)
            source_filenameval = g_strdup(set->ss.pnts[id].source_filename);
        //source_filenameval = g_path_get_basename(set->ss.pnts[id].source_filename);
        else
            source_filenameval = g_strdup(GENERAL_FILENAME);
        source_modeval = CLAMP(set->ss.pnts[id].source_mode, 0, 2);
        source_amplitudeval = set->ss.pnts[id].source_amplitude;
        source_wlval = set->ss.pnts[id].source_wl*1e6;
        source_pulsewidthval = set->ss.pnts[id].source_pulsewidth;
    } else if (type == SV_SRCTYPE_SF) {
        //name = g_strdup("Scattered field source (SF)");

        /* 1. Source properties */
        if (set->ss.sf.source_filename)
            source_filenameval = g_strdup(set->ss.sf.source_filename);
        else
            source_filenameval = g_strdup(GENERAL_FILENAME);
        source_modeval = CLAMP(set->ss.sf.source_mode, 0, 2);
        source_amplitudeval = set->ss.sf.source_amplitude;
        source_wlval = set->ss.sf.source_wl*1e6;
        source_pulsewidthval = set->ss.sf.source_pulsewidth;

        /* III. C. Incident angle group */
        source_ia_thetadegval = set->ss.sf.ia_theta * 180 / G_PI;
        source_ia_phidegval = set->ss.sf.ia_phi * 180 / G_PI;
        source_ia_psidegval = set->ss.sf.ia_psi * 180 / G_PI;
    } else if (type == SV_SRCTYPE_TSF) {
        //name = g_strdup("Total/Scattered field source (TSF)");

        /* 0a. Box properties */
        box_i0val = set->ss.tsf.box_i0;
        box_j0val = set->ss.tsf.box_j0;
        box_k0val = set->ss.tsf.box_k0;
        box_inval = set->ss.tsf.box_in;
        box_jnval = set->ss.tsf.box_jn;
        box_knval = set->ss.tsf.box_kn;
        box_boundary_skipi0val = set->ss.tsf.box_boundary_skipi0;
        box_boundary_skipinval = set->ss.tsf.box_boundary_skipin;
        box_boundary_skipj0val = set->ss.tsf.box_boundary_skipj0;
        box_boundary_skipjnval = set->ss.tsf.box_boundary_skipjn;
        box_boundary_skipk0val = set->ss.tsf.box_boundary_skipk0;
        box_boundary_skipknval = set->ss.tsf.box_boundary_skipkn;
        box_boundary_skipdepthval = set->ss.tsf.box_boundary_skipdepth;

        /* 1. Source properties */
        if (set->ss.tsf.source_filename)
            source_filenameval = g_strdup(set->ss.tsf.source_filename);
        else
            source_filenameval = g_strdup(GENERAL_FILENAME);
        source_modeval = CLAMP(set->ss.tsf.source_mode, 0, 2);
        source_amplitudeval = set->ss.tsf.source_amplitude;
        source_wlval = set->ss.tsf.source_wl*1e6;
        source_pulsewidthval = set->ss.tsf.source_pulsewidth;

        /* III. C. Incident angle group */
        source_ia_thetadegval = set->ss.tsf.ia_theta * 180 / G_PI;
        source_ia_phidegval = set->ss.tsf.ia_phi * 180 / G_PI;
        source_ia_psidegval = set->ss.tsf.ia_psi * 180 / G_PI;

        /* III. F. Z multiplier group */
        source_gaussian_mult_enable = set->ss.tsf.gaussian;
        source_gaussian_mult_center_pos_i = set->ss.tsf.gaussian_fxpos;
        source_gaussian_mult_center_pos_j = set->ss.tsf.gaussian_fypos;
        source_gaussian_mult_radius_i = set->ss.tsf.gaussian_rx;
        source_gaussian_mult_radius_i = set->ss.tsf.gaussian_ry;

        source_radial_mult_enable = set->ss.tsf.radial;
        source_radial_mult_center_pos_i = set->ss.tsf.radial_fxpos;
        source_radial_mult_center_pos_j = set->ss.tsf.radial_fypos;
        source_radial_mult_radius_i = set->ss.tsf.radial_rx;
        source_radial_mult_radius_i = set->ss.tsf.radial_ry;

        source_fiber_mult_enable = set->ss.tsf.fiber;
        source_fiber_mult_center_pos_i = set->ss.tsf.fiber_fxpos;
        source_fiber_mult_center_pos_j = set->ss.tsf.fiber_fypos;
        source_fiber_mult_radius = set->ss.tsf.fiber_radius;
        source_fiber_mult_cutoff = set->ss.tsf.fiber_cutoff;
        source_fiber_mult_epsilon_core = set->ss.tsf.fiber_epsilon_core;
        source_fiber_mult_epsilon_cladding = set->ss.tsf.fiber_epsilon_cladding;
    } else if (type == SV_SRCTYPE_TSFF) {
        //name = g_strdup("Focused Total/Scattered field source (TSFF)");

        if (set->ss.tsff.source_filename)
            source_filenameval = g_strdup(set->ss.tsff.source_filename);
        else
            source_filenameval = g_strdup(GENERAL_FILENAME);
        fs_thetamaxdegval = set->ss.tsff.focused_thetamax * 180 / G_PI;
        fs_fdistval = set->ss.tsff.focused_fdist;
        fs_nipval = set->ss.tsff.focused_nip;
        fs_mipval = set->ss.tsff.focused_mip;
        source_modeval = CLAMP(set->ss.tsff.source_mode, 0, 2);
        source_amplitudeval = set->ss.tsff.source_amplitude;
        source_wlval = set->ss.tsff.source_wl*1e6;
        source_pulsewidthval = set->ss.tsff.source_pulsewidth;
        box_i0val = set->ss.tsff.box_i0;
        box_j0val = set->ss.tsff.box_j0;
        box_k0val = set->ss.tsff.box_k0;
        box_inval = set->ss.tsff.box_in;
        box_jnval = set->ss.tsff.box_jn;
        box_knval = set->ss.tsff.box_kn;
        box_boundary_skipi0val = set->ss.tsff.box_boundary_skipi0;
        box_boundary_skipinval = set->ss.tsff.box_boundary_skipin;
        box_boundary_skipj0val = set->ss.tsff.box_boundary_skipj0;
        box_boundary_skipjnval = set->ss.tsff.box_boundary_skipjn;
        box_boundary_skipk0val = set->ss.tsff.box_boundary_skipk0;
        box_boundary_skipknval = set->ss.tsff.box_boundary_skipkn;
        box_boundary_skipdepthval = set->ss.tsff.box_boundary_skipdepth;

        fs_polarisationdegval = set->ss.tsff.focused_pol * 180 / G_PI;
    } else if (type == SV_SRCTYPE_LTSF) {
        //name = g_strdup("Layered Total/Scattered field source (LTSF)"); 

        if (set->ss.ltsf.source_filename)
            source_filenameval = g_strdup(set->ss.ltsf.source_filename);
        else
            source_filenameval = g_strdup(GENERAL_FILENAME);
        source_ia_thetadegval = set->ss.ltsf.ia_theta * 180 / G_PI;
        source_ia_phidegval = set->ss.ltsf.ia_phi * 180 / G_PI;
        source_ia_psidegval = set->ss.ltsf.ia_psi * 180 / G_PI;
        source_modeval = CLAMP(set->ss.ltsf.source_mode, 0, 2);
        source_amplitudeval = set->ss.ltsf.source_amplitude;
        source_wlval = set->ss.ltsf.source_wl*1e6;
        source_pulsewidthval = set->ss.ltsf.source_pulsewidth;
        box_i0val = set->ss.ltsf.box_i0;
        box_j0val = set->ss.ltsf.box_j0;
        box_k0val = set->ss.ltsf.box_k0;
        box_inval = set->ss.ltsf.box_in;
        box_jnval = set->ss.ltsf.box_jn;
        box_knval = set->ss.ltsf.box_kn;
        box_boundary_skipi0val = set->ss.ltsf.box_boundary_skipi0;
        box_boundary_skipinval = set->ss.ltsf.box_boundary_skipin;
        box_boundary_skipj0val = set->ss.ltsf.box_boundary_skipj0;
        box_boundary_skipjnval = set->ss.ltsf.box_boundary_skipjn;
        box_boundary_skipk0val = set->ss.ltsf.box_boundary_skipk0;
        box_boundary_skipknval = set->ss.ltsf.box_boundary_skipkn;
        box_boundary_skipdepthval = set->ss.ltsf.box_boundary_skipdepth;

        //fs_polarisationval = set->ss.tsff.pol*180/G_PI;

        if (set->ss.ltsf.layered_count > 0) {
            for (ilayers = 0; ilayers < set->ss.ltsf.layered_count; ilayers++) {
                layered_zposval[ilayers] = set->ss.ltsf.layered_zpos[ilayers];
                layered_epsval[ilayers] = set->ss.ltsf.layered_epsilon[ilayers];
                layered_muval[ilayers] = set->ss.ltsf.layered_mu[ilayers];
                layered_sigmaval[ilayers] = set->ss.ltsf.layered_sigma[ilayers];
                layered_sigastval[ilayers] = set->ss.ltsf.layered_sigast[ilayers];
                if (NULL != set->ss.ltsf.layered_material[ilayers])
                    layered_materialval[ilayers] = g_strdup(set->ss.ltsf.layered_material[ilayers]);
            }
        }

        /* III. F. Z multiplier group */
        source_gaussian_mult_enable = set->ss.ltsf.gaussian;
        source_gaussian_mult_center_pos_i = set->ss.ltsf.gaussian_fxpos;
        source_gaussian_mult_center_pos_j = set->ss.ltsf.gaussian_fypos;
        source_gaussian_mult_radius_i = set->ss.ltsf.gaussian_rx;
        source_gaussian_mult_radius_i = set->ss.ltsf.gaussian_ry;

        source_radial_mult_enable = set->ss.ltsf.radial;
        source_radial_mult_center_pos_i = set->ss.ltsf.radial_fxpos;
        source_radial_mult_center_pos_j = set->ss.ltsf.radial_fypos;
        source_radial_mult_radius_i = set->ss.ltsf.radial_rx;
        source_radial_mult_radius_i = set->ss.ltsf.radial_ry;

        source_fiber_mult_enable = set->ss.ltsf.fiber;
        source_fiber_mult_center_pos_i = set->ss.ltsf.fiber_fxpos;
        source_fiber_mult_center_pos_j = set->ss.ltsf.fiber_fypos;
        source_fiber_mult_radius = set->ss.ltsf.fiber_radius;
        source_fiber_mult_cutoff = set->ss.ltsf.fiber_cutoff;
        source_fiber_mult_epsilon_core = set->ss.ltsf.fiber_epsilon_core;
        source_fiber_mult_epsilon_cladding = set->ss.ltsf.fiber_epsilon_cladding;
    } else if (type == SV_SRCTYPE_LTSFF) {
        //name = g_strdup("Layered Focused Total/Scattered field source (LTSFF)");

        if (set->ss.ltsff.source_filename)
            source_filenameval = g_strdup(set->ss.ltsff.source_filename);
        else
            source_filenameval = g_strdup(GENERAL_FILENAME);
        fs_thetamaxdegval = set->ss.ltsff.focused_thetamax * 180 / G_PI;
        fs_fdistval = set->ss.ltsff.focused_fdist;
        fs_nipval = set->ss.ltsff.focused_nip;
        fs_mipval = set->ss.ltsff.focused_mip;
        source_modeval = CLAMP(set->ss.ltsff.source_mode, 0, 2);
        source_amplitudeval = set->ss.ltsff.source_amplitude;
        fs_polarisationdegval = set->ss.ltsff.focused_pol * 180 / G_PI;
        source_wlval = set->ss.ltsff.source_wl*1e6;
        source_pulsewidthval = set->ss.ltsff.source_pulsewidth;
        box_i0val = set->ss.ltsff.box_i0;
        box_j0val = set->ss.ltsff.box_j0;
        box_k0val = set->ss.ltsff.box_k0;
        box_inval = set->ss.ltsff.box_in;
        box_jnval = set->ss.ltsff.box_jn;
        box_knval = set->ss.ltsff.box_kn;
        box_boundary_skipi0val = set->ss.ltsff.box_boundary_skipi0;
        box_boundary_skipinval = set->ss.ltsff.box_boundary_skipin;
        box_boundary_skipj0val = set->ss.ltsff.box_boundary_skipj0;
        box_boundary_skipjnval = set->ss.ltsff.box_boundary_skipjn;
        box_boundary_skipk0val = set->ss.ltsff.box_boundary_skipk0;
        box_boundary_skipknval = set->ss.ltsff.box_boundary_skipkn;
        box_boundary_skipdepthval = set->ss.ltsff.box_boundary_skipdepth;

        if (set->ss.ltsff.layered_count > 0) {
            for (ilayers = 0; ilayers < set->ss.ltsff.layered_count; ilayers++) {
                layered_zposval[ilayers] = set->ss.ltsff.layered_zpos[ilayers];
                layered_epsval[ilayers] = set->ss.ltsff.layered_epsilon[ilayers];
                layered_muval[ilayers] = set->ss.ltsff.layered_mu[ilayers];
                layered_sigmaval[ilayers] = set->ss.ltsff.layered_sigma[ilayers];
                layered_sigastval[ilayers] = set->ss.ltsff.layered_sigast[ilayers];
                if (NULL != set->ss.ltsff.layered_material[ilayers])
                    layered_materialval[ilayers] = g_strdup(set->ss.ltsff.layered_material[ilayers]);
            }
        }
    }

    xgc->par_lock_controls_update = TRUE;

    if (type == SV_SRCTYPE_POINT) {
        /* III. 1. Point source properties */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.point_origin_position_i_spin), source_point_origin_position_ival);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.point_origin_position_j_spin), source_point_origin_position_jval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.point_origin_position_k_spin), source_point_origin_position_kval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.point_origin_theta_spin), source_point_origin_thetaval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.point_origin_phi_spin), source_point_origin_phival);

        /* 1. Source properties */
        gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs.sc.source_mode), source_modeval);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.sc.source_filename), source_filenameval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_amplitude_spin), source_amplitudeval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_wavelength_spin), source_wlval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_pulsewidth_spin), source_pulsewidthval);
    }

    if (type == SV_SRCTYPE_TSF || type == SV_SRCTYPE_TSFF || type == SV_SRCTYPE_LTSF || type == SV_SRCTYPE_LTSFF) {
        /* III. A. Source box group */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_i0_spin), box_i0val);    /* virtual size in voxels */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_j0_spin), box_j0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_k0_spin), box_k0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_in_spin), box_inval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_jn_spin), box_jnval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_kn_spin), box_knval);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_enable_skipdepth), (box_boundary_skipdepthval != -1));
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.box_skipdepth_spin), box_boundary_skipdepthval);

        // real size in micrometers
        g_snprintf(buff, sizeof(buff), "%g", ABS(box_inval - box_i0val)*xgc->data.set.sp.dx*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.sc.box_size_x_um), buff);
        g_snprintf(buff, sizeof(buff), "%g", ABS(box_jnval - box_j0val)*xgc->data.set.sp.dy*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.sc.box_size_y_um), buff);
        g_snprintf(buff, sizeof(buff), "%g", ABS(box_knval - box_k0val)*xgc->data.set.sp.dz*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.sc.box_size_z_um), buff);

        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipi0), box_boundary_skipi0val);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipj0), box_boundary_skipj0val);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipk0), box_boundary_skipk0val);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipin), box_boundary_skipinval);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipjn), box_boundary_skipjnval);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.box_boundary_skipkn), box_boundary_skipknval);
    }

    if (type == SV_SRCTYPE_SF || type == SV_SRCTYPE_TSF || type == SV_SRCTYPE_TSFF || type == SV_SRCTYPE_LTSF || type == SV_SRCTYPE_LTSFF) {
        /* III. B. Source group */
        gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs.sc.source_mode), source_modeval);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.sc.source_filename), source_filenameval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_amplitude_spin), source_amplitudeval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_wavelength_spin), source_wlval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.source_pulsewidth_spin), source_pulsewidthval);
    }

    if (type == SV_SRCTYPE_SF || type == SV_SRCTYPE_TSF || type == SV_SRCTYPE_LTSF) {
        /* III. C. Incident angle group */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.ia_theta_spin), source_ia_thetadegval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.ia_phi_spin), source_ia_phidegval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.ia_psi_spin), source_ia_psidegval);

        /* III. F. Z multiplier group */
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.gaussian_mult_enable), source_gaussian_mult_enable);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_center_pos_i_spin), source_gaussian_mult_center_pos_i);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_center_pos_j_spin), source_gaussian_mult_center_pos_j);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_radius_i_spin), source_gaussian_mult_radius_i);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_radius_j_spin), source_gaussian_mult_radius_j);

        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.radial_mult_enable), source_radial_mult_enable);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_center_pos_i_spin), source_radial_mult_center_pos_i);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_center_pos_j_spin), source_radial_mult_center_pos_j);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_radius_i_spin), source_radial_mult_radius_i);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.gaussian_mult_radius_j_spin), source_radial_mult_radius_j);

        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.sc.fiber_mult_enable), source_fiber_mult_enable);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_center_pos_i_spin), source_fiber_mult_center_pos_i);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_center_pos_j_spin), source_fiber_mult_center_pos_j);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_radius_spin), source_fiber_mult_radius);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_cutoff_spin), source_fiber_mult_cutoff);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_epsilon_core_spin), source_fiber_mult_epsilon_core);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.fiber_mult_epsilon_cladding_spin), source_fiber_mult_epsilon_cladding);
    }

    if (type == SV_SRCTYPE_TSFF || type == SV_SRCTYPE_LTSFF) {
        /* III. D. Focused source group */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.fs_thetamax_deg_spin), fs_thetamaxdegval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.fs_fdist_spin), fs_fdistval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.fs_pol_deg_spin), fs_polarisationdegval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.fs_nip_spin), fs_nipval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.sc.fs_mip_spin), fs_mipval);
    }

    if (type == SV_SRCTYPE_LTSF) {
        /* III. E. Layered source group */
        LTV_Data ltvd;
        memset(&ltvd, 0, sizeof(LTV_Data));

        ltvd.nlayers = set->ss.ltsf.layered_count;
        if (set->ss.ltsf.layered_count > 0) {
            for (ilayers = 0; ilayers < set->ss.ltsf.layered_count; ilayers++) {
                ltvd.zpos[ilayers] = layered_zposval[ilayers];
                ltvd.epsilon[ilayers] = layered_epsval[ilayers];
                ltvd.mu[ilayers] = layered_muval[ilayers];
                ltvd.sigma[ilayers] = layered_sigmaval[ilayers];
                ltvd.sigast[ilayers] = layered_sigastval[ilayers];
                if(layered_materialval[ilayers] != NULL)
                    ltvd.material[ilayers] = g_strdup(layered_materialval[ilayers]);
            }
        }

        //set_ltv_items(GTK_TREE_VIEW(xgc->cs.sc.ls_prop_wrapper), (gpointer)&ltvd);
        set_ltv_items((TreeViewWrapper*)xgc->cs.sc.ls_prop_wrapper, (gpointer)&ltvd);
    }

    if (type == SV_SRCTYPE_LTSFF) {
        /* III. E. Layered source group */
        LTV_Data ltvd;
        memset(&ltvd, 0, sizeof(LTV_Data));

        ltvd.nlayers = set->ss.ltsff.layered_count;
        if (set->ss.ltsff.layered_count > 0) {
            for (ilayers = 0; ilayers < set->ss.ltsff.layered_count; ilayers++) {
                ltvd.zpos[ilayers] = layered_zposval[ilayers];
                ltvd.epsilon[ilayers] = layered_epsval[ilayers];
                ltvd.mu[ilayers] = layered_muval[ilayers];
                ltvd.sigma[ilayers] = layered_sigmaval[ilayers];
                ltvd.sigast[ilayers] = layered_sigastval[ilayers];
                ltvd.material[ilayers] = g_strdup(layered_materialval[ilayers]);
            }
        }

        //set_ltv_items(GTK_TREE_VIEW(xgc->cs.sc.ls_prop_wrapper), (gpointer)&ltvd);
        set_ltv_items((TreeViewWrapper *)xgc->cs.sc.ls_prop_wrapper, (gpointer)&ltvd);
    }

    xgc->par_lock_controls_update = FALSE;

    for (i = 0; i < NUM_LSMP; i++)
        g_free(layered_materialval[i]);

    g_free(source_filenameval);
    //g_free(material_vtx_filenameval);


    return FALSE;
} /* set_values_sf */

static gboolean
set_values_output(XGControls* xgc, SvSet *set, SvType type, gint id)
{
    gchar *filenameval = NULL;
    gint point_ival = POINT_OUTPUT_ORIGIN_I, point_jval = POINT_OUTPUT_ORIGIN_J, point_kval = POINT_OUTPUT_ORIGIN_K;
    gint image_ival = IMAGE_OUTPUT_ORIGIN_I;
    gint plane_ival = PLANE_OUTPUT_ORIGIN_I;
    gint output_dirval = OUTPUT_PLANE;
    gint box_i0val = OUTPUT_BOX_I0, box_j0val = OUTPUT_BOX_J0, box_k0val = OUTPUT_BOX_K0;
    gint box_inval = OUTPUT_BOX_IN, box_jnval = OUTPUT_BOX_JN, box_knval = OUTPUT_BOX_KN;
    gint output_componentval = OUTPUT_COMPONENT, output_stepval = OUTPUT_STEP, output_startval = OUTPUT_START, output_stopval = OUTPUT_STOP, output_formatval = OUTPUT_FORMAT;
    gdouble sum_output_epsilonval = OUTPUT_SUM_EPS, sum_output_muval = OUTPUT_SUM_MU, sum_output_sigmaval = OUTPUT_SUM_SIGMA, sum_output_sigastval = OUTPUT_SUM_SIGAST;

    gchar buff[256] = {0};

    if (type == SV_OUTTYPE_GENERAL) {
        if (set->so.outfile)
            filenameval = g_strdup(set->so.outfile);
        else
            filenameval = g_strdup(GENERAL_FILENAME);
    } else if (type == SV_OUTTYPE_POINT) {
        point_ival = set->so.pnts[id].i;
        point_jval = set->so.pnts[id].j;
        point_kval = set->so.pnts[id].k;

        output_componentval = set->so.pnts[id].component;
        output_stepval = set->so.pnts[id].step;
        output_startval = set->so.pnts[id].start;
        output_stopval = set->so.pnts[id].stop;
        output_formatval = set->so.pnts[id].format;

        if (set->so.pnts[id].filebase)
            filenameval = g_strdup(set->so.pnts[id].filebase);
        else
            filenameval = g_strdup(GENERAL_FILENAME);
    } else if (type == SV_OUTTYPE_IMAGE) {
        if (set->so.imgs[id].i != -1) {
            image_ival = set->so.imgs[id].i;
            output_dirval = 0;
        } else if (set->so.imgs[id].j != -1) {
            image_ival = set->so.imgs[id].j;
            output_dirval = 1;
        } else if (set->so.imgs[id].k != -1) {
            image_ival = set->so.imgs[id].k;
            output_dirval = 2;
        }

        output_componentval = set->so.imgs[id].component;
        output_stepval = set->so.imgs[id].step;

        if (set->so.imgs[id].filebase)
            filenameval = g_strdup(set->so.imgs[id].filebase);
        else
            filenameval = g_strdup(GENERAL_FILENAME);
    } else if (type == SV_OUTTYPE_PLANE) {
        if (set->so.plns[id].i != -1) {
            plane_ival = set->so.plns[id].i;
            output_dirval = 0;
        } else if (set->so.plns[id].j != -1) {
            plane_ival = set->so.plns[id].j;
            output_dirval = 1;
        } else if (set->so.plns[id].k != -1) {
            plane_ival = set->so.plns[id].k;
            output_dirval = 2;
        }

        output_componentval = set->so.plns[id].component;
        output_stepval = set->so.plns[id].step;

        if (set->so.plns[id].filebase)
            filenameval = g_strdup(set->so.plns[id].filebase);
        else
            filenameval = g_strdup(GENERAL_FILENAME);

        output_formatval = set->so.plns[id].format;
        output_startval = set->so.plns[id].start;
        output_stopval = set->so.plns[id].stop;
    } else if (type == SV_OUTTYPE_VOLUME) {
        output_componentval = set->so.cubs[id].component;
        output_stepval = set->so.cubs[id].step;

        if (set->so.cubs[id].filebase)
            filenameval = g_strdup(set->so.cubs[id].filebase);
        else
            filenameval = g_strdup(GENERAL_FILENAME);

        output_formatval = set->so.cubs[id].format;
        output_startval = set->so.cubs[id].start;
        output_stopval = set->so.cubs[id].stop;
    } else if (type == SV_OUTTYPE_SUM || type == SV_OUTTYPE_SUMTAB) {
        box_i0val = set->so.sums[id].box_i0;
        box_j0val = set->so.sums[id].box_j0;
        box_k0val = set->so.sums[id].box_k0;
        box_inval = set->so.sums[id].box_in;
        box_jnval = set->so.sums[id].box_jn;
        box_knval = set->so.sums[id].box_kn;

        output_componentval = set->so.sums[id].component;
        output_stepval = set->so.sums[id].step;

        if (set->so.sums[id].filename)
            filenameval = g_strdup(set->so.sums[id].filename);
        else
            filenameval = g_strdup(GENERAL_FILENAME);

        sum_output_epsilonval = set->so.sums[id].layered_epsilon;
        sum_output_muval = set->so.sums[id].layered_mu;
        sum_output_sigmaval = set->so.sums[id].layered_sigma;
        sum_output_sigastval = set->so.sums[id].layered_sigast;
    } else if (type == SV_OUTTYPE_FORCE) {
        box_i0val = set->so.forces[id].box_i0;
        box_j0val = set->so.forces[id].box_j0;
        box_k0val = set->so.forces[id].box_k0;
        box_inval = set->so.forces[id].box_in;
        box_jnval = set->so.forces[id].box_jn;
        box_knval = set->so.forces[id].box_kn;

        output_stepval = set->so.forces[id].step;

        if (set->so.forces[id].filename)
            filenameval = g_strdup(set->so.forces[id].filename);
        else
            filenameval = g_strdup(GENERAL_FILENAME);
    }

    xgc->par_lock_controls_update = TRUE;

    if (type == SV_OUTTYPE_GENERAL) {
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.general_output_filename), filenameval);
    } else if (type == SV_OUTTYPE_POINT) {
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.point_output_origin_position_i_spin), point_ival);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.point_output_origin_position_j_spin), point_jval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.point_output_origin_position_k_spin), point_kval);

        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.point_output_step_spin), output_stepval);
        gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs.oc.point_output_component), output_componentval);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.point_output_filename), filenameval);
    } else if (type == SV_OUTTYPE_IMAGE) {
        gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs.oc.image_output_plane), output_dirval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.image_output_pos_spin), image_ival);
        gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs.oc.image_output_component), output_componentval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.image_output_step_spin), output_stepval);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.image_output_label), filenameval);
    } else if (type == SV_OUTTYPE_PLANE) {
        gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs.oc.plane_output_plane), output_dirval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.plane_output_pos_spin), plane_ival);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.plane_output_step_spin), output_stepval);
        gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs.oc.plane_output_component), output_componentval);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.plane_output_filename), filenameval);
        gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs.oc.plane_output_format), output_formatval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.plane_output_start_spin), output_startval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.plane_output_stop_spin), output_stopval);
    } else if (type == SV_OUTTYPE_VOLUME) {
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.volume_output_step_spin), output_stepval);
        gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs.oc.volume_output_component), output_componentval);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.volume_output_filename), filenameval);

        gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs.oc.volume_output_format), output_formatval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.volume_output_start_spin), output_startval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.volume_output_stop_spin), output_stopval);
    } else if (type == SV_OUTTYPE_SUM || type == SV_OUTTYPE_SUMTAB) {
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_i0_spin), box_i0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_j0_spin), box_j0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_k0_spin), box_k0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_in_spin), box_inval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_jn_spin), box_jnval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_kn_spin), box_knval);

        // real size in micrometers
        g_snprintf(buff, sizeof(buff), "%g", ABS(box_inval - box_i0val)*xgc->data.set.sp.dx*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.box_output_size_x_um), buff);
        g_snprintf(buff, sizeof(buff), "%g", ABS(box_jnval - box_j0val)*xgc->data.set.sp.dy*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.box_output_size_y_um), buff);
        g_snprintf(buff, sizeof(buff), "%g", ABS(box_knval - box_k0val)*xgc->data.set.sp.dz*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.box_output_size_z_um), buff);

        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.sum_output_step_spin), output_stepval);
        gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs.oc.sum_output_component), output_componentval);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.sum_output_filename), filenameval);

        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.sum_output_epsilon_spin), sum_output_epsilonval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.sum_output_mu_spin), sum_output_muval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.sum_output_sigma_spin), sum_output_sigmaval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.sum_output_sigast_spin), sum_output_sigastval);
    } else if (type == SV_OUTTYPE_FORCE) {
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_i0_spin), box_i0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_j0_spin), box_j0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_k0_spin), box_k0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_in_spin), box_inval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_jn_spin), box_jnval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.box_output_kn_spin), box_knval);

        // real size in micrometers
        g_snprintf(buff, sizeof(buff), "%g", ABS(box_inval - box_i0val)*xgc->data.set.sp.dx*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.box_output_size_x_um), buff);
        g_snprintf(buff, sizeof(buff), "%g", ABS(box_jnval - box_j0val)*xgc->data.set.sp.dy*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.box_output_size_y_um), buff);
        g_snprintf(buff, sizeof(buff), "%g", ABS(box_knval - box_k0val)*xgc->data.set.sp.dz*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.box_output_size_z_um), buff);

        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.oc.force_output_step_spin), output_stepval);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.oc.force_output_filename), filenameval);
    }

    xgc->par_lock_controls_update = FALSE;

    g_free(filenameval);


    return FALSE;
} /* set_values_output() */

static gboolean
set_values_material(XGControls* xgc, SvSet *set, SvType type, gint id)
{
    gboolean material_mode_checkval = MATERIAL_MODECHECK;
    gint material_smooth_stepsval = MATERIAL_SMOOTHSTEPS;

    gchar *material_voxel_filenameval = NULL, *material_vector_filenameval = NULL/*, *material_vtx_filenameval = NULL*/;
    gint material_grow_i0val = MATERIAL_GROW_I0, material_grow_j0val = MATERIAL_GROW_J0, material_grow_k0val = MATERIAL_GROW_K0;
    gint material_grow_inval = set->sp.xres - MATERIAL_GROW_INBORDER;
    gint material_grow_jnval = set->sp.yres - MATERIAL_GROW_JNBORDER;
    gint material_grow_knval = set->sp.zres - MATERIAL_GROW_KNBORDER;
    gint material_grow_skipi0val = MATERIAL_GROW_SKIP, material_grow_skipj0val = MATERIAL_GROW_SKIP, material_grow_skipk0val = MATERIAL_GROW_SKIP, material_grow_skipinval = MATERIAL_GROW_SKIP, material_grow_skipjnval = MATERIAL_GROW_SKIP, material_grow_skipknval = MATERIAL_GROW_SKIP;
    gint material_grow_attachindexval = MATERIAL_GROW_ATTACHINDEX, material_grow_addindexval = MATERIAL_GROW_ADDINDEX, material_grow_subsamplingval = MATERIAL_GROW_SUBSAMPLING, material_grow_seedval = MATERIAL_GROW_SEED, material_grow_nstepsval = MATERIAL_GROW_NSTEPS;
    gdouble material_grow_mobilityval = MATERIAL_GROW_MOBILITY, material_grow_probabilityval = MATERIAL_GROW_PROBABILITY;

    gint material_rough_matindexval = MATERIAL_ROUGHNESS_MATINDEX, material_rough_voidindexval = MATERIAL_ROUGHNESS_VOIDINDEX, material_rough_iterationsval = MATERIAL_ROUGHNESS_ITERATIONS, material_rough_seedval = MATERIAL_ROUGHNESS_SEED, material_rough_radiuspeakval = MATERIAL_ROUGHNESS_RADIUSPEAK, material_rough_radiusspanval = MATERIAL_ROUGHNESS_RADIUSSPAN;
    gdouble material_rough_probabilityval = MATERIAL_ROUGHNESS_PROBABILITY;
    gint material_spectral_matindexval = MATERIAL_SPECTRAL_MATINDEX, material_spectral_seedval = MATERIAL_SPECTRAL_SEED;
    gdouble material_spectral_sigmaval = MATERIAL_SPECTRAL_SIGMA, material_spectral_tval = MATERIAL_SPECTRAL_T;    
    gint material_expr_i0val = MATERIAL_EXPR_I0, material_expr_j0val = MATERIAL_EXPR_J0, material_expr_k0val = MATERIAL_EXPR_K0;
    gint material_expr_inval = set->sp.xres - MATERIAL_EXPR_INBORDER;
    gint material_expr_jnval = set->sp.yres - MATERIAL_EXPR_JNBORDER;
    gint material_expr_knval = set->sp.zres - MATERIAL_EXPR_KNBORDER;
    gint material_expr_matindexval = set->sm.expr_matindex[id];
    gint material_expr_voidindexval = MATERIAL_EXPR_VOIDINDEX, material_expr_maxdistval = MATERIAL_EXPR_MAXDIST, material_expr_distmodeval = MATERIAL_EXPR_DISTMODE;
    gchar material_expr_expr[MATERIAL_EXPR_EXPR_CHARS] = MATERIAL_EXPR_EXPR;

    /* V. Media */

    if (type == SV_TYPE_MATERIAL_PROP) {
        if (set->sm.in_voxel_filename)
            material_voxel_filenameval = g_strdup(set->sm.in_voxel_filename);
        else
            material_voxel_filenameval = g_strdup(GENERAL_FILENAME);
        if (set->sm.in_vector_filename)
            material_vector_filenameval = g_strdup(set->sm.in_vector_filename);
        else
            material_vector_filenameval = g_strdup(GENERAL_FILENAME);

        material_mode_checkval = set->sm.matmode_check;
        material_smooth_stepsval = set->sm.smooth;
    } else if (type == SV_TYPE_MATERIAL_GROW) {
        /* V. 2. Add growth modifier */

        material_grow_i0val = set->sm.grow_i0[id];
        material_grow_j0val = set->sm.grow_j0[id];
        material_grow_k0val = set->sm.grow_k0[id];
        material_grow_inval = set->sm.grow_in[id];
        material_grow_jnval = set->sm.grow_jn[id];
        material_grow_knval = set->sm.grow_kn[id];

        material_grow_skipi0val = set->sm.grow_skipi0[id];
        material_grow_skipj0val = set->sm.grow_skipj0[id];
        material_grow_skipk0val = set->sm.grow_skipk0[id];
        material_grow_skipinval = set->sm.grow_skipin[id];
        material_grow_skipjnval = set->sm.grow_skipjn[id];
        material_grow_skipknval = set->sm.grow_skipkn[id];

        material_grow_attachindexval = set->sm.grow_attachindex[id];
        material_grow_addindexval = set->sm.grow_addindex[id];        
        material_grow_subsamplingval = set->sm.grow_subsampling[id];
        material_grow_seedval = set->sm.grow_seed[id];
        material_grow_nstepsval = set->sm.grow_nsteps[id];
        material_grow_mobilityval = set->sm.grow_mobility[id];
        material_grow_probabilityval = set->sm.grow_probability[id];
    } else if (type == SV_TYPE_MATERIAL_ROUGHNESS) {
        /* V. 3. Add roughness modifier */

        material_rough_matindexval = set->sm.rough_matindex[id];
        material_rough_voidindexval = set->sm.rough_voidindex[id];
        material_rough_iterationsval = set->sm.rough_iterations[id];
        material_rough_seedval = set->sm.rough_seed[id];
        material_rough_radiuspeakval = set->sm.rough_radius_peak[id];
        material_rough_radiusspanval = set->sm.rough_radius_span[id];
        material_rough_probabilityval = set->sm.rough_probability[id];
    } else if (type == SV_TYPE_MATERIAL_SPECTRAL) {
        /* V. 4. Add spectral modifier */

        material_spectral_matindexval = set->sm.spectral_matindex[id];
        material_spectral_sigmaval = set->sm.spectral_sigma[id];
        material_spectral_tval = set->sm.spectral_t[id];
        material_spectral_seedval = set->sm.spectral_seed[id];
    } else if (type == SV_TYPE_MATERIAL_EXPRESSION) {
        /* V. 5. Add expression modifier */

        material_expr_i0val = set->sm.expr_i0[id];
        material_expr_j0val = set->sm.expr_j0[id];
        material_expr_k0val = set->sm.expr_k0[id];
        material_expr_inval = set->sm.expr_in[id];
        material_expr_jnval = set->sm.expr_jn[id];
        material_expr_knval = set->sm.expr_kn[id];
        material_expr_matindexval = set->sm.expr_matindex[id];
        material_expr_voidindexval = set->sm.expr_voidindex[id];
        material_expr_maxdistval = set->sm.expr_maxdist[id];
        material_expr_distmodeval = set->sm.expr_distmode[id];
        g_snprintf(material_expr_expr, sizeof(gchar)*MATERIAL_EXPR_EXPR_CHARS, set->sm.expr_expr[id]);        
    }

    xgc->par_lock_controls_update = TRUE;

    if (type == SV_TYPE_MATERIAL_PROP) {
        /* V. 1. Material properties */

        gtk_entry_set_text(GTK_ENTRY(xgc->cs.mc.material_voxel_filename), material_voxel_filenameval);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.mc.material_vector_filename), material_vector_filenameval);
        //gtk_entry_set_text(GTK_ENTRY(xgc->cs.mc.material_vtx_filename), material_vtx_filenameval);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.mc.material_modecheck), material_mode_checkval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.material_smoothsteps_spin), material_smooth_stepsval);
    } else if (type == SV_TYPE_MATERIAL_GROW) {
        /* V. 2. Add growth modifier */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_i0_spin), material_grow_i0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_j0_spin), material_grow_j0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_k0_spin), material_grow_k0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_in_spin), material_grow_inval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_jn_spin), material_grow_jnval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_kn_spin), material_grow_knval);

        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.mc.grow_skipi0), material_grow_skipi0val);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.mc.grow_skipj0), material_grow_skipj0val);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.mc.grow_skipk0), material_grow_skipk0val);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.mc.grow_skipin), material_grow_skipinval);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.mc.grow_skipjn), material_grow_skipjnval);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.mc.grow_skipkn), material_grow_skipknval);

        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_attachindex_spin), material_grow_attachindexval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_addindex_spin), material_grow_addindexval);        
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_subsampling_spin), material_grow_subsamplingval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_seed_spin), material_grow_seedval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_nsteps_spin), material_grow_nstepsval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_mobility_spin), material_grow_mobilityval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.grow_probability_spin), material_grow_probabilityval);
    } else if (type == SV_TYPE_MATERIAL_ROUGHNESS) {
        /* V. 3. Add roughness modifier */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.rough_matindex_spin), material_rough_matindexval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.rough_voidindex_spin), material_rough_voidindexval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.rough_iterations_spin), material_rough_iterationsval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.rough_seed_spin), material_rough_seedval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.rough_radiuspeak_spin), material_rough_radiuspeakval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.rough_radiusspan_spin), material_rough_radiusspanval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.rough_probability_spin), material_rough_probabilityval);
    } else if (type == SV_TYPE_MATERIAL_SPECTRAL) {
        /* V. 4. Add spectral modifier */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.spectral_sigma_spin), material_spectral_sigmaval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.spectral_t_spin), material_spectral_tval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.spectral_matindex_spin), material_spectral_matindexval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.spectral_matindex_spin), material_spectral_seedval);
    } else if (type == SV_TYPE_MATERIAL_EXPRESSION) {
        /* V. 5. Add expression modifier */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_i0_spin), material_expr_i0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_j0_spin), material_expr_j0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_k0_spin), material_expr_k0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_in_spin), material_expr_inval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_jn_spin), material_expr_jnval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_kn_spin), material_expr_knval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_matindex_spin), material_expr_matindexval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_voidindex_spin), material_expr_voidindexval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_maxdist_spin), material_expr_maxdistval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.mc.expr_distmode_spin), material_expr_distmodeval);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.mc.expr_expr), material_expr_expr);
    }

    xgc->par_lock_controls_update = FALSE;

    g_free(material_voxel_filenameval);
    g_free(material_vector_filenameval);


    return FALSE;
} /* set_values_material() */

static gboolean
set_values_nfff(XGControls* xgc, SvSet *set, SvType type, gint id)
{
    gint nfff_box_i0val = NFFF_BOX_VERTEX0, nfff_box_j0val = NFFF_BOX_VERTEX0, nfff_box_k0val = NFFF_BOX_VERTEX0, nfff_box_inval = NFFF_BOX_VERTEXN, nfff_box_jnval = NFFF_BOX_VERTEXN, nfff_box_knval = NFFF_BOX_VERTEXN;
    gint nfffp_pos_ival = NFFFP_POS_I, nfffp_pos_jval = NFFFP_POS_J, nfffp_pos_kval = NFFFP_POS_K;
    gchar *filenameval = NULL;
    gint nfffa_thetaresval = NFFFA_THETARES, nfffa_phiresval = NFFFA_PHIRES, nfffa_savefileval = NFFFA_SAVEFILE, nfffa_radiusval = NFFFA_RADIUS;
    gdouble nfffa_thetafromval = NFFFA_THETAFROMDEG, nfffa_phifromval = NFFFA_PHIFROMDEG, nfffa_thetatoval = NFFFA_THETATODEG, nfffa_phitoval = NFFFA_PHITODEG;
    gint nfff_box_boundary_skipi0 = NFFF_BOX_BOUNDARY_SKIP, nfff_box_boundary_skipin = NFFF_BOX_BOUNDARY_SKIP, nfff_box_boundary_skipj0 = NFFF_BOX_BOUNDARY_SKIP, nfff_box_boundary_skipjn = NFFF_BOX_BOUNDARY_SKIP, nfff_box_boundary_skipk0 = NFFF_BOX_BOUNDARY_SKIP, nfff_box_boundary_skipkn = NFFF_BOX_BOUNDARY_SKIP;
    gint nfff_skipi0_jmin = NFFF_SKIPI0, nfff_skipi0_kmin = NFFF_SKIPI0, nfff_skipi0_jmax = NFFF_SKIPI0, nfff_skipi0_kmax = NFFF_SKIPI0, nfff_skipin_jmin = NFFF_SKIPIN, nfff_skipin_kmin = NFFF_SKIPIN, nfff_skipin_jmax = NFFF_SKIPIN, nfff_skipin_kmax = NFFF_SKIPIN;
    gint nfff_skipj0_imin = NFFF_SKIPJ0, nfff_skipj0_kmin = NFFF_SKIPJ0, nfff_skipj0_imax = NFFF_SKIPJ0, nfff_skipj0_kmax = NFFF_SKIPJ0, nfff_skipjn_imin = NFFF_SKIPJN, nfff_skipjn_kmin = NFFF_SKIPJN, nfff_skipjn_imax = NFFF_SKIPJN, nfff_skipjn_kmax = NFFF_SKIPJN;
    gint nfff_skipk0_imin = NFFF_SKIPK0, nfff_skipk0_jmin = NFFF_SKIPK0, nfff_skipk0_imax = NFFF_SKIPK0, nfff_skipk0_jmax = NFFF_SKIPK0, nfff_skipkn_imin = NFFF_SKIPKN, nfff_skipkn_jmin = NFFF_SKIPKN, nfff_skipkn_imax = NFFF_SKIPKN, nfff_skipkn_jmax = NFFF_SKIPKN;
    gint pnfff_box_i0val = PNFFF_BOX_VERTEX0, pnfff_box_j0val = PNFFF_BOX_VERTEX0, pnfff_box_k0val = PNFFF_BOX_VERTEX0, pnfff_box_inval = PNFFF_BOX_VERTEXN, pnfff_box_jnval = PNFFF_BOX_VERTEXN, pnfff_box_knval = PNFFF_BOX_VERTEXN;
    gint pnfff_box_boundary_skipk0 = PNFFF_BOX_BOUNDARY_SKIP, pnfff_box_boundary_skipkn = PNFFF_BOX_BOUNDARY_SKIP;
    gint pnfff_skipk0_imin = PNFFF_SKIPK0, pnfff_skipk0_jmin = PNFFF_SKIPK0, pnfff_skipk0_imax = PNFFF_SKIPK0, pnfff_skipk0_jmax = PNFFF_SKIPK0, pnfff_skipkn_imin = PNFFF_SKIPKN, pnfff_skipkn_jmin = PNFFF_SKIPKN, pnfff_skipkn_imax = PNFFF_SKIPKN, pnfff_skipkn_jmax = PNFFF_SKIPKN;
    gint pnfff_integration_xmin = PNFFF_INTEGRATION_XMIN, pnfff_integration_xmax = PNFFF_INTEGRATION_XMAX, pnfff_integration_ymin = PNFFF_INTEGRATION_YMIN, pnfff_integration_ymax = PNFFF_INTEGRATION_YMAX;

    gchar buff[256] = {0};

    if (type == SV_TYPE_NFFF_BOX) {
        /* VII. 1. Near field to far field transform */

        nfff_box_i0val = set->sf.box_i0;
        nfff_box_j0val = set->sf.box_j0;
        nfff_box_k0val = set->sf.box_k0;
        nfff_box_inval = set->sf.box_in;
        nfff_box_jnval = set->sf.box_jn;
        nfff_box_knval = set->sf.box_kn;

        nfff_box_boundary_skipi0 = set->sf.box_boundary_skipi0;
        nfff_skipi0_jmin = set->sf.skipi0_jmin;
        nfff_skipi0_kmin = set->sf.skipi0_kmin;
        nfff_skipi0_jmax = set->sf.skipi0_jmax;
        nfff_skipi0_kmax = set->sf.skipi0_kmax;

        nfff_box_boundary_skipin = set->sf.box_boundary_skipin;
        nfff_skipin_jmin = set->sf.skipin_jmin;
        nfff_skipin_kmin = set->sf.skipin_kmin;
        nfff_skipin_jmax = set->sf.skipin_jmax;
        nfff_skipin_kmax = set->sf.skipin_kmax;

        nfff_box_boundary_skipj0 = set->sf.box_boundary_skipj0;
        nfff_skipj0_imin = set->sf.skipj0_imin;
        nfff_skipj0_kmin = set->sf.skipj0_kmin;
        nfff_skipj0_imax = set->sf.skipj0_imax;
        nfff_skipj0_kmax = set->sf.skipj0_kmax;

        nfff_box_boundary_skipjn = set->sf.box_boundary_skipjn;
        nfff_skipjn_imin = set->sf.skipjn_imin;
        nfff_skipjn_kmin = set->sf.skipjn_kmin;
        nfff_skipjn_imax = set->sf.skipjn_imax;
        nfff_skipjn_kmax = set->sf.skipjn_kmax;

        nfff_box_boundary_skipk0 = set->sf.box_boundary_skipk0;
        nfff_skipk0_imin = set->sf.skipk0_imin;
        nfff_skipk0_jmin = set->sf.skipk0_jmin;
        nfff_skipk0_imax = set->sf.skipk0_imax;
        nfff_skipk0_jmax = set->sf.skipk0_jmax;

        nfff_box_boundary_skipkn = set->sf.box_boundary_skipkn;
        nfff_skipkn_imin = set->sf.skipkn_imin;
        nfff_skipkn_jmin = set->sf.skipkn_jmin;
        nfff_skipkn_imax = set->sf.skipkn_imax;
        nfff_skipkn_jmax = set->sf.skipkn_jmax;
    } else if (type == SV_TYPE_NFFF_POINT) {
        /* VII. 2. Near field to far field point */

        nfffp_pos_ival = set->sf.ri[id];
        nfffp_pos_jval = set->sf.rj[id];
        nfffp_pos_kval = set->sf.rk[id];
        if (set->sf.source_filename[id])
            filenameval = g_strdup(set->sf.source_filename[id]);
        else
            filenameval = g_strdup(GENERAL_FILENAME);
    } else if (type == SV_TYPE_PNFFF_POINT) {
        /* VIII. 2. Periodic near field to far field point */

        nfffp_pos_ival = set->spf.ri[id];
        nfffp_pos_jval = set->spf.rj[id];
        nfffp_pos_kval = set->spf.rk[id];
        if (set->spf.source_filename[id])
            filenameval = g_strdup(set->spf.source_filename[id]);
        else
            filenameval = g_strdup(GENERAL_FILENAME);
    } else if (type == SV_TYPE_NFFF_AREA) {
        /* VII. 3. Near field to far field area */

        nfffa_thetaresval = set->sf.area_thetares[id];
        nfffa_phiresval = set->sf.area_phires[id];
        nfffa_radiusval = set->sf.area_radius[id];
        nfffa_thetafromval = set->sf.area_thetafrom[id] * 180 / G_PI;
        nfffa_thetatoval = set->sf.area_thetato[id] * 180 / G_PI;
        nfffa_phifromval = set->sf.area_phifrom[id] * 180 / G_PI;
        nfffa_phitoval = set->sf.area_phito[id] * 180 / G_PI;
        nfffa_savefileval = set->sf.area_savefile[id];
    } else if (type == SV_TYPE_PNFFF_AREA) {
        /* VIII. 3. Periodic near field to far field area */

        nfffa_thetaresval = set->spf.area_thetares[id];
        nfffa_phiresval = set->spf.area_phires[id];
        nfffa_radiusval = set->spf.area_radius[id];
        nfffa_thetafromval = set->spf.area_thetafrom[id] * 180 / G_PI;
        nfffa_thetatoval = set->spf.area_thetato[id] * 180 / G_PI;
        nfffa_phifromval = set->spf.area_phifrom[id] * 180 / G_PI;
        nfffa_phitoval = set->spf.area_phito[id] * 180 / G_PI;
        nfffa_savefileval = set->spf.area_savefile[id];
    } else if (type == SV_TYPE_PNFFF_BOX) {
        /* VIII. 1. Periodic near field to far field transform */

        pnfff_box_i0val = set->spf.box_i0;
        pnfff_box_j0val = set->spf.box_j0;
        pnfff_box_k0val = set->spf.box_k0;
        pnfff_box_inval = set->spf.box_in;
        pnfff_box_jnval = set->spf.box_jn;
        pnfff_box_knval = set->spf.box_kn;

        pnfff_integration_xmin = set->spf.pimin;
        pnfff_integration_xmax = set->spf.pimax;
        pnfff_integration_ymin = set->spf.pjmin;
        pnfff_integration_ymax = set->spf.pjmax;

        pnfff_box_boundary_skipk0 = set->spf.box_boundary_skipk0;
        pnfff_skipk0_imin = set->spf.skipk0_imin;
        pnfff_skipk0_jmin = set->spf.skipk0_jmin;
        pnfff_skipk0_imax = set->spf.skipk0_imax;
        pnfff_skipk0_jmax = set->spf.skipk0_jmax;

        pnfff_box_boundary_skipkn = set->spf.box_boundary_skipkn;
        pnfff_skipkn_imin = set->spf.skipkn_imin;
        pnfff_skipkn_jmin = set->spf.skipkn_jmin;
        pnfff_skipkn_imax = set->spf.skipkn_imax;
        pnfff_skipkn_jmax = set->spf.skipkn_jmax;
    }

    xgc->par_lock_controls_update = TRUE;

    if (type == SV_TYPE_NFFF_BOX) {
        /* VII. 1. Near field to far field transform box */

        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_box_i0_spin), nfff_box_i0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_box_j0_spin), nfff_box_j0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_box_k0_spin), nfff_box_k0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_box_in_spin), nfff_box_inval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_box_jn_spin), nfff_box_jnval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_box_kn_spin), nfff_box_knval);

        // real size in micrometers
        g_snprintf(buff, sizeof(buff), "%g", ABS(nfff_box_inval - nfff_box_i0val)*xgc->data.set.sp.dx*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.fc.nfff_box_size_x_um), buff);
        g_snprintf(buff, sizeof(buff), "%g", ABS(nfff_box_jnval - nfff_box_j0val)*xgc->data.set.sp.dy*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.fc.nfff_box_size_y_um), buff);
        g_snprintf(buff, sizeof(buff), "%g", ABS(nfff_box_knval - nfff_box_k0val)*xgc->data.set.sp.dz*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.fc.nfff_box_size_z_um), buff);

        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.nfff_box_boundary_skipi0), nfff_box_boundary_skipi0);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipi0_jmin_spin), nfff_skipi0_jmin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipi0_kmin_spin), nfff_skipi0_kmin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipi0_jmax_spin), nfff_skipi0_jmax);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipi0_kmax_spin), nfff_skipi0_kmax);

        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.nfff_box_boundary_skipin), nfff_box_boundary_skipin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipin_jmin_spin), nfff_skipin_jmin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipin_kmin_spin), nfff_skipin_kmin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipin_jmax_spin), nfff_skipin_jmax);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipin_kmax_spin), nfff_skipin_kmax);

        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.nfff_box_boundary_skipj0), nfff_box_boundary_skipj0);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipj0_imin_spin), nfff_skipj0_imin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipj0_kmin_spin), nfff_skipj0_kmin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipj0_imax_spin), nfff_skipj0_imax);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipj0_kmax_spin), nfff_skipj0_kmax);

        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.nfff_box_boundary_skipjn), nfff_box_boundary_skipjn);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipjn_imin_spin), nfff_skipjn_imin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipjn_kmin_spin), nfff_skipjn_kmin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipjn_imax_spin), nfff_skipjn_imax);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipjn_kmax_spin), nfff_skipjn_kmax);

        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.nfff_box_boundary_skipk0), nfff_box_boundary_skipk0);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipk0_imin_spin), nfff_skipk0_imin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipk0_jmin_spin), nfff_skipk0_jmin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipk0_imax_spin), nfff_skipk0_imax);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipk0_jmax_spin), nfff_skipk0_jmax);

        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.nfff_box_boundary_skipkn), nfff_box_boundary_skipkn);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipkn_imin_spin), nfff_skipkn_imin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipkn_jmin_spin), nfff_skipkn_jmin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipkn_imax_spin), nfff_skipkn_imax);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfff_skipkn_jmax_spin), nfff_skipkn_jmax);
    } else if ((type == SV_TYPE_NFFF_POINT) || (type == SV_TYPE_PNFFF_POINT)) {
        /* VII. 2. Near field to far field point */
        /* VIII. 2. Periodic near field to far field point */

        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffp_pos_i_spin), nfffp_pos_ival);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffp_pos_j_spin), nfffp_pos_jval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffp_pos_k_spin), nfffp_pos_kval);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.fc.nfffp_filename), filenameval);
    } else if ((type == SV_TYPE_NFFF_AREA) || (type == SV_TYPE_PNFFF_AREA)) {
        /* VII. 3. Near field to far field area */
        /* VIII. 3. Periodic near field to far field area */

        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_thetares_spin), nfffa_thetaresval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_phires_spin), nfffa_phiresval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_radius_spin), nfffa_radiusval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_thetafrom_spin), nfffa_thetafromval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_thetato_spin), nfffa_thetatoval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_phifrom_spin), nfffa_phifromval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.nfffa_phito_spin), nfffa_phitoval);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.nfffa_savefile), nfffa_savefileval);
    } else if (type == SV_TYPE_PNFFF_BOX) {
        /* VIII. 1. Periodic near field to far field transform box */

        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_box_i0_spin), pnfff_box_i0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_box_j0_spin), pnfff_box_j0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_box_k0_spin), pnfff_box_k0val);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_box_in_spin), pnfff_box_inval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_box_jn_spin), pnfff_box_jnval);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_box_kn_spin), pnfff_box_knval);

        // real size in micrometers
        g_snprintf(buff, sizeof(buff), "%g", ABS(pnfff_box_inval - pnfff_box_i0val)*xgc->data.set.sp.dx*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.fc.pnfff_box_size_x_um), buff);
        g_snprintf(buff, sizeof(buff), "%g", ABS(pnfff_box_jnval - pnfff_box_j0val)*xgc->data.set.sp.dy*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.fc.pnfff_box_size_y_um), buff);
        g_snprintf(buff, sizeof(buff), "%g", ABS(pnfff_box_knval - pnfff_box_k0val)*xgc->data.set.sp.dz*1e6);
        gtk_entry_set_text(GTK_ENTRY(xgc->cs.fc.pnfff_box_size_z_um), buff);

        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_integration_xmin_spin), pnfff_integration_xmin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_integration_xmax_spin), pnfff_integration_xmax);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_integration_ymin_spin), pnfff_integration_ymin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_integration_ymax_spin), pnfff_integration_ymax);

        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.pnfff_box_boundary_skipk0), pnfff_box_boundary_skipk0);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipk0_imin_spin), pnfff_skipk0_imin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipk0_jmin_spin), pnfff_skipk0_jmin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipk0_imax_spin), pnfff_skipk0_imax);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipk0_jmax_spin), pnfff_skipk0_jmax);

        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs.fc.pnfff_box_boundary_skipkn), pnfff_box_boundary_skipkn);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipkn_imin_spin), pnfff_skipkn_imin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipkn_jmin_spin), pnfff_skipkn_jmin);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipkn_imax_spin), pnfff_skipkn_imax);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs.fc.pnfff_skipkn_jmax_spin), pnfff_skipkn_jmax);
    }

    xgc->par_lock_controls_update = FALSE;

    g_free(filenameval);

    return FALSE;
} /* set_values_nfff */

void
par_controls_row_select_process(XGControls *xgc, gint id, GtkTreePath *path)
{
    gboolean writeme = FALSE, matreload = FALSE;
    gint index;

    /*use id to call appropriate callback*/
    if (id == SET_POOL) {
        //writeme = matreload = xsv_get_pool(xgc);
        //writeme = matreload = xsv_set_values_sf(xgc, &(xgc->data.set), SV_TYPE_POOL, id);
        writeme = matreload = set_values_general(xgc, &(xgc->data.set), SV_TYPE_POOL, id);
    } else if (id == SET_BASIC) {
        //writeme = matreload = xsv_get_basic(xgc);
        //writeme = matreload = xsv_set_values_sf(xgc, &(xgc->data.set), SV_TYPE_BASIC, id);
        writeme = matreload = set_values_general(xgc, &(xgc->data.set), SV_TYPE_BASIC, id);
    } else if (id == SET_SF)
        writeme = set_values_sf(xgc, &(xgc->data.set), SV_SRCTYPE_SF, id);
    else if (id == SET_TSF)
        writeme = set_values_sf(xgc, &(xgc->data.set), SV_SRCTYPE_TSF, id);
    else if (id == SET_TSFF)
        writeme = set_values_sf(xgc, &(xgc->data.set), SV_SRCTYPE_TSFF, id);
    else if (id == SET_LTSF)
        writeme = set_values_sf(xgc, &(xgc->data.set), SV_SRCTYPE_LTSF, id);
    else if (id == SET_LTSFF)
        writeme = set_values_sf(xgc, &(xgc->data.set), SV_SRCTYPE_LTSFF, id);
    else if (id >= SET_PSOURCE && id < SET_POUT)
        writeme = set_values_sf(xgc, &(xgc->data.set), SV_SRCTYPE_POINT, id - SET_PSOURCE);
    else if (id == SET_BND)
        //writeme = xsv_get_bnd(xgc);
        //writeme = xsv_set_values_sf(xgc, &(xgc->data.set), SV_TYPE_BND, id);
        writeme = set_values_general(xgc, &(xgc->data.set), SV_TYPE_BND, id);
    else if (id == SET_MEDIUM)
        //writeme = xsv_get_mat(xgc);
        writeme = set_values_material(xgc, &(xgc->data.set), SV_TYPE_MATERIAL_PROP, id);
    else if (id == SET_OUT)                         /* general output */
        writeme = set_values_output(xgc, &(xgc->data.set), SV_OUTTYPE_GENERAL, id);
    //writeme = xsv_get_out(xgc);
    else if (id >= SET_POUT && id < SET_IOUT)       /* point output */
        writeme = set_values_output(xgc, &(xgc->data.set), SV_OUTTYPE_POINT, id - SET_POUT);
    //writeme = xsv_get_pout(xgc, id-SET_POUT, 0);
    else if (id >= SET_IIOUT && id < SET_COUT)      /* image output */
        writeme = set_values_output(xgc, &(xgc->data.set), SV_OUTTYPE_IMAGE, id - SET_IIOUT);
    //writeme = xsv_get_imageout(xgc, id-SET_IIOUT, 1);        
    else if (id >= SET_IOUT && id < SET_IIOUT)      /* plane output */
        writeme = set_values_output(xgc, &(xgc->data.set), SV_OUTTYPE_PLANE, id - SET_IOUT);
    //writeme = xsv_get_imageout(xgc, id-SET_IOUT, 0);
    else if (id >= SET_COUT && id < SET_SOUT)       /* volume output */
        writeme = set_values_output(xgc, &(xgc->data.set), SV_OUTTYPE_VOLUME, id - SET_COUT);
    //writeme = xsv_get_vout(xgc, id-SET_COUT);
    else if (id >= SET_SOUT && id < SET_FOUT)       /* sum output */
        writeme = set_values_output(xgc, &(xgc->data.set), SV_OUTTYPE_SUM, id - SET_SOUT);
    //writeme = xsv_get_pout(xgc, id-SET_SOUT, 1);
    else if (id >= SET_FOUT && id < SET_GROW)       /* force output */
        writeme = set_values_output(xgc, &(xgc->data.set), SV_OUTTYPE_FORCE, id - SET_FOUT);
    //writeme = xsv_get_pout(xgc, id-SET_FOUT, 2);
    else if (id >= SET_GROW && id < SET_ROUGHNESS)        /* V. 2. Add growth modifier */
        writeme = set_values_material(xgc, &(xgc->data.set), SV_TYPE_MATERIAL_GROW, id - SET_GROW);
    else if (id >= SET_ROUGHNESS && id < SET_SPECTRAL)    /* V. 3. Add roughness modifier */
        writeme = set_values_material(xgc, &(xgc->data.set), SV_TYPE_MATERIAL_ROUGHNESS, id - SET_ROUGHNESS);
    else if (id >= SET_SPECTRAL && id < SET_EXPRESSION) /* V. 4. Add spectral modifier */
        writeme = set_values_material(xgc, &(xgc->data.set), SV_TYPE_MATERIAL_SPECTRAL, id - SET_SPECTRAL);
    else if (id >= SET_EXPRESSION && id < SET_NFAREA)    /* V. 5. Add expression modifier */
        writeme = set_values_material(xgc, &(xgc->data.set), SV_TYPE_MATERIAL_EXPRESSION, id - SET_EXPRESSION);
    else if (id == SET_NFFF)                        /* VII. 1. Near field to far field transform */
        //writeme = xsv_get_nfff(xgc, 0);
        writeme = set_values_nfff(xgc, &(xgc->data.set), SV_TYPE_NFFF_BOX, id);
    else if (id == SET_PNFFF)                       /* VIII. 1. Periodic near field to far field transform */
        //writeme = xsv_get_nfff(xgc, 1);
        writeme = set_values_nfff(xgc, &(xgc->data.set), SV_TYPE_PNFFF_BOX, id);
    else if (id >= SET_NFFFP && id < SET_PNFAREA)   /* VII. 2. Near field to far field point */
        //writeme = xsv_get_nfffp(xgc, id-SET_NFFFP, 0);
        writeme = set_values_nfff(xgc, &(xgc->data.set), SV_TYPE_NFFF_POINT, id - SET_NFFFP);
    else if (id >= SET_NFAREA && id < SET_NFFFP)    /* VII. 3. Near field to far field area */
        //writeme = xsv_get_nfarea(xgc, id-SET_NFAREA, 0);
        writeme = set_values_nfff(xgc, &(xgc->data.set), SV_TYPE_NFFF_AREA, id - SET_NFAREA);
    else if (id >= SET_PNFAREA && id < SET_PNFFFP)  /* VIII. 3. Periodic near field to far field area */
        //writeme = xsv_get_nfarea(xgc, id-SET_PNFAREA, 1);
        writeme = set_values_nfff(xgc, &(xgc->data.set), SV_TYPE_PNFFF_AREA, id - SET_PNFAREA);
    else if (id >= SET_PNFFFP)                      /* VIII. 2. Periodic near field to far field point */
        //writeme = xsv_get_nfffp(xgc, id-SET_PNFFFP, 1);
        writeme = set_values_nfff(xgc, &(xgc->data.set), SV_TYPE_PNFFF_POINT, id - SET_PNFFFP);

    //ss_set_controls_visible(xgc, id);


    /*rewrite the temporary parfile and then reload it and interpret completely in order to assure consistency between its text and gui interpretation*/
    if (writeme) {
        if (write_parfile(xgc)) {
            index = gtk_tree_path_get_indices(path)[0];
            path = gtk_tree_path_new_from_indices(index, -1);
            gtk_tree_view_expand_row(GTK_TREE_VIEW(xgc->tv_par), path, TRUE);
        }
    }
    if (matreload) { //in some cases, like complete change of pool geometry, it is safer to reinitialize all the materials as well.
        if (write_matfile(xgc)) {
            index = gtk_tree_path_get_indices(path)[0];
            path = gtk_tree_path_new_from_indices(0, -1);
            gtk_tree_view_expand_row(GTK_TREE_VIEW(xgc->tv_mat), path, TRUE);
        }
    }
}

/* boundary sensitive functions */
static void 
set_boundary_controls_sensitive(XGControls *xgc, gint i, gboolean val)
{
    gtk_widget_set_sensitive(xgc->cs.boc.bdepthspin[i], val);
    gtk_widget_set_sensitive(xgc->cs.boc.bmspin[i], val);
    gtk_widget_set_sensitive(xgc->cs.boc.bkappaspin[i], val);
    gtk_widget_set_sensitive(xgc->cs.boc.baspin[i], val);
    gtk_widget_set_sensitive(xgc->cs.boc.bsigmaspin[i], val);
}

void
par_controls_sensitive(XGControls *xgc, gint id)
{
    if (id >= SET_PSOURCE && id < SET_POUT) {
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.point_origin_theta_spin), (xgc->data.set.ss.pnts[id - SET_PSOURCE].source_mode != 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.point_origin_phi_spin), (xgc->data.set.ss.pnts[id - SET_PSOURCE].source_mode != 0));

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_filename), (xgc->data.set.ss.pnts[id - SET_PSOURCE].source_mode == 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_button_change), (xgc->data.set.ss.pnts[id - SET_PSOURCE].source_mode == 0));

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_amplitude_spin), (xgc->data.set.ss.pnts[id - SET_PSOURCE].source_mode != 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_wavelength_spin), (xgc->data.set.ss.pnts[id - SET_PSOURCE].source_mode != 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_pulsewidth_spin), (xgc->data.set.ss.pnts[id - SET_PSOURCE].source_mode == 2));
    } else if (id == SET_SF) {
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_filename), (xgc->data.set.ss.sf.source_mode == 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_button_change), (xgc->data.set.ss.sf.source_mode == 0));

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_amplitude_spin), (xgc->data.set.ss.sf.source_mode != 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_wavelength_spin), (xgc->data.set.ss.sf.source_mode != 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_pulsewidth_spin), (xgc->data.set.ss.sf.source_mode == 2));
    } else if (id == SET_TSF) {
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.label_box_skipdepth), (xgc->data.set.ss.tsf.box_boundary_skipdepth != -1));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.box_skipdepth_spin), (xgc->data.set.ss.tsf.box_boundary_skipdepth != -1));

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_filename), (xgc->data.set.ss.tsf.source_mode == 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_button_change), (xgc->data.set.ss.tsf.source_mode == 0));

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_amplitude_spin), (xgc->data.set.ss.tsf.source_mode != 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_wavelength_spin), (xgc->data.set.ss.tsf.source_mode != 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_pulsewidth_spin), (xgc->data.set.ss.tsf.source_mode == 2));

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.gaussian_mult_center_pos_i_spin), xgc->data.set.ss.tsf.gaussian);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.gaussian_mult_center_pos_j_spin), xgc->data.set.ss.tsf.gaussian);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.gaussian_mult_radius_i_spin), xgc->data.set.ss.tsf.gaussian);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.gaussian_mult_radius_j_spin), xgc->data.set.ss.tsf.gaussian);

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.radial_mult_center_pos_i_spin), xgc->data.set.ss.tsf.radial);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.radial_mult_center_pos_j_spin), xgc->data.set.ss.tsf.radial);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.radial_mult_radius_i_spin), xgc->data.set.ss.tsf.radial);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.radial_mult_radius_j_spin), xgc->data.set.ss.tsf.radial);

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.fiber_mult_center_pos_i_spin), xgc->data.set.ss.tsf.fiber);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.fiber_mult_center_pos_j_spin), xgc->data.set.ss.tsf.fiber);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.fiber_mult_radius_spin), xgc->data.set.ss.tsf.fiber);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.fiber_mult_cutoff_spin), xgc->data.set.ss.tsf.fiber);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.fiber_mult_epsilon_core_spin), xgc->data.set.ss.tsf.fiber);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.fiber_mult_epsilon_cladding_spin), xgc->data.set.ss.tsf.fiber);
    } else if (id == SET_TSFF) {
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.label_box_skipdepth), (xgc->data.set.ss.tsff.box_boundary_skipdepth != -1));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.box_skipdepth_spin), (xgc->data.set.ss.tsff.box_boundary_skipdepth != -1));

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_filename), (xgc->data.set.ss.tsff.source_mode == 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_button_change), (xgc->data.set.ss.tsff.source_mode == 0));

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_amplitude_spin), (xgc->data.set.ss.tsff.source_mode != 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_wavelength_spin), (xgc->data.set.ss.tsff.source_mode != 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_pulsewidth_spin), (xgc->data.set.ss.tsff.source_mode == 2));
    } else if (id == SET_LTSF) {
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.label_box_skipdepth), (xgc->data.set.ss.ltsf.box_boundary_skipdepth != -1));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.box_skipdepth_spin), (xgc->data.set.ss.ltsf.box_boundary_skipdepth != -1));

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_filename), (xgc->data.set.ss.ltsf.source_mode == 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_button_change), (xgc->data.set.ss.ltsf.source_mode == 0));

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_amplitude_spin), (xgc->data.set.ss.ltsf.source_mode != 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_wavelength_spin), (xgc->data.set.ss.ltsf.source_mode != 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_pulsewidth_spin), (xgc->data.set.ss.ltsf.source_mode == 2));

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.gaussian_mult_center_pos_i_spin), xgc->data.set.ss.ltsf.gaussian);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.gaussian_mult_center_pos_j_spin), xgc->data.set.ss.ltsf.gaussian);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.gaussian_mult_radius_i_spin), xgc->data.set.ss.ltsf.gaussian);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.gaussian_mult_radius_j_spin), xgc->data.set.ss.ltsf.gaussian);

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.radial_mult_center_pos_i_spin), xgc->data.set.ss.ltsf.radial);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.radial_mult_center_pos_j_spin), xgc->data.set.ss.ltsf.radial);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.radial_mult_radius_i_spin), xgc->data.set.ss.ltsf.radial);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.radial_mult_radius_j_spin), xgc->data.set.ss.ltsf.radial);

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.fiber_mult_center_pos_i_spin), xgc->data.set.ss.ltsf.fiber);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.fiber_mult_center_pos_j_spin), xgc->data.set.ss.ltsf.fiber);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.fiber_mult_radius_spin), xgc->data.set.ss.ltsf.fiber);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.fiber_mult_cutoff_spin), xgc->data.set.ss.ltsf.fiber);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.fiber_mult_epsilon_core_spin), xgc->data.set.ss.ltsf.fiber);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.fiber_mult_epsilon_cladding_spin), xgc->data.set.ss.ltsf.fiber);
    } else if (id == SET_LTSFF) {
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.label_box_skipdepth), (xgc->data.set.ss.ltsff.box_boundary_skipdepth != -1));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.box_skipdepth_spin), (xgc->data.set.ss.ltsff.box_boundary_skipdepth != -1));

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_filename), (xgc->data.set.ss.ltsff.source_mode == 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_button_change), (xgc->data.set.ss.ltsff.source_mode == 0));

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_amplitude_spin), (xgc->data.set.ss.ltsff.source_mode != 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_wavelength_spin), (xgc->data.set.ss.ltsff.source_mode != 0));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.sc.source_pulsewidth_spin), (xgc->data.set.ss.ltsff.source_mode == 2));
    } else if (id >= SET_POUT && id < SET_IOUT) {
    } else if (id == SET_BASIC) {
        gtk_widget_set_sensitive(xgc->cs.bac.basic_ugpu_index_spin, xgc->data.set.sc.usegpu);
    } else if (id == SET_BND) {
        set_boundary_controls_sensitive(xgc, 0, (xgc->data.set.sb.bx0 == 3));
        set_boundary_controls_sensitive(xgc, 1, (xgc->data.set.sb.bxn == 3));
        set_boundary_controls_sensitive(xgc, 2, (xgc->data.set.sb.by0 == 3));
        set_boundary_controls_sensitive(xgc, 3, (xgc->data.set.sb.byn == 3));
        set_boundary_controls_sensitive(xgc, 4, (xgc->data.set.sb.bz0 == 3));
        set_boundary_controls_sensitive(xgc, 5, (xgc->data.set.sb.bzn == 3));

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.boc.mbposspin[0]), (xgc->data.set.smb.bx0 == SV_BOUNDARY_PERIODIC));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.boc.mbposspin[1]), (xgc->data.set.smb.bxn == SV_BOUNDARY_PERIODIC));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.boc.mbposspin[2]), (xgc->data.set.smb.by0 == SV_BOUNDARY_PERIODIC));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.boc.mbposspin[3]), (xgc->data.set.smb.byn == SV_BOUNDARY_PERIODIC));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.boc.mbposspin[4]), (xgc->data.set.smb.bz0 == SV_BOUNDARY_PERIODIC));
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.boc.mbposspin[5]), (xgc->data.set.smb.bzn == SV_BOUNDARY_PERIODIC));
    } else if (id == SET_MEDIUM) {
    } else if (id == SET_NFFF) {                        /* VII. 1. Near field to far field transform box */
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipi0_jmin_spin), xgc->data.set.sf.box_boundary_skipi0);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipi0_kmin_spin), xgc->data.set.sf.box_boundary_skipi0);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipi0_jmax_spin), xgc->data.set.sf.box_boundary_skipi0);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipi0_kmax_spin), xgc->data.set.sf.box_boundary_skipi0);

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipin_jmin_spin), xgc->data.set.sf.box_boundary_skipin);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipin_kmin_spin), xgc->data.set.sf.box_boundary_skipin);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipin_jmax_spin), xgc->data.set.sf.box_boundary_skipin);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipin_kmax_spin), xgc->data.set.sf.box_boundary_skipin);

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipj0_imin_spin), xgc->data.set.sf.box_boundary_skipj0);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipj0_kmin_spin), xgc->data.set.sf.box_boundary_skipj0);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipj0_imax_spin), xgc->data.set.sf.box_boundary_skipj0);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipj0_kmax_spin), xgc->data.set.sf.box_boundary_skipj0);

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipjn_imin_spin), xgc->data.set.sf.box_boundary_skipjn);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipjn_kmin_spin), xgc->data.set.sf.box_boundary_skipjn);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipjn_imax_spin), xgc->data.set.sf.box_boundary_skipjn);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipjn_kmax_spin), xgc->data.set.sf.box_boundary_skipjn);

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipk0_imin_spin), xgc->data.set.sf.box_boundary_skipk0);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipk0_jmin_spin), xgc->data.set.sf.box_boundary_skipk0);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipk0_imax_spin), xgc->data.set.sf.box_boundary_skipk0);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipk0_jmax_spin), xgc->data.set.sf.box_boundary_skipk0);

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipkn_imin_spin), xgc->data.set.sf.box_boundary_skipkn);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipkn_jmin_spin), xgc->data.set.sf.box_boundary_skipkn);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipkn_imax_spin), xgc->data.set.sf.box_boundary_skipkn);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.nfff_skipkn_jmax_spin), xgc->data.set.sf.box_boundary_skipkn);
    } else if (id == SET_PNFFF) {                       /* VIII. 1. Periodic near field to far field transform box */
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.pnfff_skipk0_imin_spin), xgc->data.set.spf.box_boundary_skipk0);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.pnfff_skipk0_jmin_spin), xgc->data.set.spf.box_boundary_skipk0);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.pnfff_skipk0_imax_spin), xgc->data.set.spf.box_boundary_skipk0);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.pnfff_skipk0_jmax_spin), xgc->data.set.spf.box_boundary_skipk0);

        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.pnfff_skipkn_imin_spin), xgc->data.set.spf.box_boundary_skipkn);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.pnfff_skipkn_jmin_spin), xgc->data.set.spf.box_boundary_skipkn);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.pnfff_skipkn_imax_spin), xgc->data.set.spf.box_boundary_skipkn);
        gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs.fc.pnfff_skipkn_jmax_spin), xgc->data.set.spf.box_boundary_skipkn);
    } else if (id >= SET_NFFFP && id < SET_PNFAREA) {   /* VII. 2. Near field to far field point */
    } else if (id >= SET_NFAREA && id < SET_NFFFP) {    /* VII. 3. Near field to far field area */
    } else if (id >= SET_PNFAREA && id < SET_PNFFFP) {  /* VIII. 3. Periodic near field to far field area */
    } else if (id >= SET_PNFFFP) {                      /* VIII. 2. Periodic near field to far field point */
    }
}