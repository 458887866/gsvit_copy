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


/* write.c :
*  Parameter and material file write functions
*/

#include "write.h"

gchar* 
get_temporary_par_filename(XGControls *xgc)
{
    gchar *filename = NULL;

    /*write parfile and read it exactly same way as in gsvit*/

    if (xgc->curdir)
        filename = g_build_filename(xgc->curdir, "xsvit_temporary.par", NULL);
    else if (xgc->tmpfilename)
        filename = g_strdup(xgc->tmpfilename);
    else {
#ifdef G_OS_WIN32
        filename = g_build_filename(g_get_home_dir(), ".xsvit", "xsvit_temporary.par", NULL);
#else
        gint ft;
        GError *glib_error = NULL;
        GtkWidget *dialog;

        if (!(ft = g_file_open_tmp("xsvit.parfile.XXXXXX", &filename, &glib_error))) {
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

    return filename;
}

gchar*
get_temporary_mat_filename(XGControls *xgc)
{
    gchar *filename = NULL;

    if (xgc->curdir)
        filename = g_build_filename(xgc->curdir, "xsvit_temporary.mat", NULL);
    /*else if (xgc->tmpfilename)
        filename = g_strdup(xgc->tmpfilename);*/
    else {
#ifdef G_OS_WIN32
        filename = g_build_filename(g_get_home_dir(), ".xsvit", "xsvit_temporary.mat", NULL);
#else
        gint ft;
        GError *glib_error = NULL;
        GtkWidget *dialog;

        if (!(ft = g_file_open_tmp("xsvit.parfile.XXXXXX", &filename, &glib_error))) {
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

    return filename;
}

gchar*
cpstring(SvCompType component)
{
    if (component == SV_COMP_EX)
        return g_strdup("Ex");
    if (component == SV_COMP_EY)
        return g_strdup("Ey");
    if (component == SV_COMP_EZ)
        return g_strdup("Ez");
    if (component == SV_COMP_HX)
        return g_strdup("Hx");
    if (component == SV_COMP_HY)
        return g_strdup("Hy");
    if (component == SV_COMP_HZ)
        return g_strdup("Hz");
    if (component == SV_COMP_ALL)
        return g_strdup("All");
    if (component == SV_COMP_CUR)
        return g_strdup("Cur");
    if (component == SV_COMP_EPSILON)
        return g_strdup("Epsilon");
    if (component == SV_COMP_MU)
        return g_strdup("Mu");
    if (component == SV_COMP_SIGMA)
        return g_strdup("Sigma");
    if (component == SV_COMP_SIGAST)
        return g_strdup("Sigast");
    if (component == SV_COMP_MAT)
        return g_strdup("Material");
    else
        return g_strdup("Unknown");
}

gchar*
vcpstring(SvOutputVolumeType component)
{
    if (component == SV_OVOLUME_EX)
        return g_strdup("Ex");
    if (component == SV_OVOLUME_EY)
        return g_strdup("Ey");
    if (component == SV_OVOLUME_EZ)
        return g_strdup("Ez");
    if (component == SV_OVOLUME_HX)
        return g_strdup("Hx");
    if (component == SV_OVOLUME_HY)
        return g_strdup("Hy");
    if (component == SV_OVOLUME_HZ)
        return g_strdup("Hz");
    if (component == SV_OVOLUME_ALL)
        return g_strdup("All");
    if (component == SV_OVOLUME_EPSILON)
        return g_strdup("Epsilon");
    if (component == SV_OVOLUME_MU)
        return g_strdup("Mu");
    if (component == SV_OVOLUME_SIGMA)
        return g_strdup("Sigma");
    if (component == SV_OVOLUME_SIGAST)
        return g_strdup("Sigast");
    if (component == SV_OVOLUME_MAT)
        return g_strdup("Material");
    if (component == SV_OVOLUME_MATTYPE)
        return g_strdup("Mattype");
    if (component == SV_OVOLUME_ABS)
        return g_strdup("Abs");
    else
        return g_strdup("Unknown");
}

gchar*
scpstring(SvSumType component)
{
    if (component == SV_SUM_EX)
        return g_strdup("Ex");
    if (component == SV_SUM_EY)
        return g_strdup("Ey");
    if (component == SV_SUM_EZ)
        return g_strdup("Ez");
    if (component == SV_SUM_ALL)
        return g_strdup("All");
    if (component == SV_SUM_ABS)
        return g_strdup("Abs");
    else
        return g_strdup("Unknown");
}

gchar*
planestring(gint i, gint j, gint k)
{
    if (j == -1 && k == -1)
        return g_strdup_printf("yz plane at x=%d", i);
    if (i == -1 && k == -1)
        return g_strdup_printf("xz plane at y=%d", j);
    if (i == -1 && j == -1)
        return g_strdup_printf("xy plane at z=%d", k);

    return g_strdup("wrong plane");
}

/* write data to parfile */
gboolean
write_parfile(XGControls *xgc)
{
    FILE *fw;
    GtkWidget *dialog;
    //gchar *filename = NULL;    
    GError* glib_error = NULL;
    gchar *contents;
    gint i;


    /*if (xgc->curdir)
    filename = g_build_filename(xgc->curdir, "xsvit_temporary.par", NULL);
    else {
    #ifdef G_OS_WIN32
    filename = g_build_filename(g_get_home_dir(), ".xsvit", "xsvit_temporary.par", NULL);
    #else
    gint ft;
    if (!(ft = g_file_open_tmp("xsvit.parfile.XXXXXX", &filename, &glib_error)))
    {
    dialog = gtk_message_dialog_new (GTK_WINDOW(xgc->toplevel),
    GTK_DIALOG_DESTROY_WITH_PARENT,
    GTK_MESSAGE_ERROR,
    GTK_BUTTONS_CLOSE,
    "Error opening temporary file: %s",
    g_strdup(glib_error->message));
    gtk_dialog_run (GTK_DIALOG (dialog));
    gtk_widget_destroy (dialog);
    }
    close(ft);
    #endif
    }
    xgc->tmpfilename = g_strdup(filename);*/

    if (xgc->par_file_success == FALSE) {
        printf("Error parsing parameter file. Parameter file will not be saved.\n");
    }

    if (NULL == xgc->tmpfilename)
        xgc->tmpfilename = get_temporary_par_filename(xgc);

    // printf("write parfile will be writing to %s assigning it to temporary location %s\n", filename, xgc->tmpfilename);

    if (xgc->tmpfilename) {
        fw = fopen(xgc->tmpfilename, "w");
        if (!fw) {
            dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                                            GTK_DIALOG_DESTROY_WITH_PARENT,
                                            GTK_MESSAGE_ERROR,
                                            GTK_BUTTONS_CLOSE,
                                            "Error opening temporary file: %s",
                                            g_strdup(glib_error->message));
            gtk_dialog_run(GTK_DIALOG(dialog));
            gtk_widget_destroy(dialog);
        } else {
            fprintf(fw, "POOL\n%d %d %d %g %g %g\n\n", xgc->data.set.sp.xres, xgc->data.set.sp.yres, xgc->data.set.sp.zres, xgc->data.set.sp.dx, xgc->data.set.sp.dy, xgc->data.set.sp.dz);
            fprintf(fw, "COMP\n%d\n\n", xgc->data.set.sc.nsteps);
            fprintf(fw, "VERBOSE\n%d\n\n", xgc->data.set.sc.verbose);
            fprintf(fw, "THREADS\n%d\n\n", xgc->data.set.sc.nthreads);
            fprintf(fw, "GPU\n%d\n\n", xgc->data.set.sc.usegpu);
            if (xgc->data.set.sc.ugpu[0] == 1)
                fprintf(fw, "UGPU\n0\n\n");
            else if (xgc->data.set.sc.ugpu[1] == 1)
                fprintf(fw, "UGPU\n1\n\n");
            else if (xgc->data.set.sc.ugpu[2] == 1)
                fprintf(fw, "UGPU\n2\n\n");
            else if (xgc->data.set.sc.ugpu[3] == 1)
                fprintf(fw, "UGPU\n3\n\n");

            for (i = 0; i < xgc->data.set.ss.npnts; i++) {
                if (xgc->data.set.ss.pnts[i].source_mode == 0)
                    fprintf(fw, "SOURCE_POINT\n%d %d %d 0 %s\n\n",
                            xgc->data.set.ss.pnts[i].point_origin_position_i, xgc->data.set.ss.pnts[i].point_origin_position_j, xgc->data.set.ss.pnts[i].point_origin_position_k, xgc->data.set.ss.pnts[i].source_filename);
                else if (xgc->data.set.ss.pnts[i].source_mode == 1)
                    fprintf(fw, "SOURCE_POINT\n%d %d %d 1 %g %g %g %g\n\n",
                            xgc->data.set.ss.pnts[i].point_origin_position_i, xgc->data.set.ss.pnts[i].point_origin_position_j, xgc->data.set.ss.pnts[i].point_origin_position_k, xgc->data.set.ss.pnts[i].source_wl,
                            xgc->data.set.ss.pnts[i].source_amplitude, xgc->data.set.ss.pnts[i].point_origin_theta, xgc->data.set.ss.pnts[i].point_origin_phi);
                else if (xgc->data.set.ss.pnts[i].source_mode == 2)
                    fprintf(fw, "SOURCE_POINT\n%d %d %d 2 %g %g %g %g %g\n\n",
                            xgc->data.set.ss.pnts[i].point_origin_position_i, xgc->data.set.ss.pnts[i].point_origin_position_j, xgc->data.set.ss.pnts[i].point_origin_position_k, xgc->data.set.ss.pnts[i].source_wl,
                            xgc->data.set.ss.pnts[i].source_pulsewidth, xgc->data.set.ss.pnts[i].source_amplitude,
                            xgc->data.set.ss.pnts[i].point_origin_theta, xgc->data.set.ss.pnts[i].point_origin_phi);
            }
            for (i = 0; i < xgc->data.set.ss.nlocals; i++) {
                fprintf(fw, "SOURCE_LOCAL\n%d %d %d %d %d %d %g %g %g %g %g %d %d\n\n",
                        xgc->data.set.ss.locals[i].box_i0, xgc->data.set.ss.locals[i].box_j0, xgc->data.set.ss.locals[i].box_k0,
                        xgc->data.set.ss.locals[i].box_in, xgc->data.set.ss.locals[i].box_jn, xgc->data.set.ss.locals[i].box_kn,
                        xgc->data.set.ss.locals[i].layered_epsilon, xgc->data.set.ss.locals[i].density, xgc->data.set.ss.locals[i].strength,
                        xgc->data.set.ss.locals[i].lambda_peak, xgc->data.set.ss.locals[i].lambda_width,
                        xgc->data.set.ss.locals[i].source_mode, xgc->data.set.ss.locals[i].startfrom);
            }
            if (xgc->data.set.ss.ext.filebase_ex) {
                fprintf(fw, "SOURCE_EXT\n%d %d %d %d %d %d %d %d %d %d %d %d %s %s %s %s %s %s\n\n",
                        xgc->data.set.ss.ext.i, xgc->data.set.ss.ext.j, xgc->data.set.ss.ext.k,
                        xgc->data.set.ss.ext.ijstart, xgc->data.set.ss.ext.jkstart, xgc->data.set.ss.ext.extxres, xgc->data.set.ss.ext.extyres,
                        xgc->data.set.ss.ext.iextfrom, xgc->data.set.ss.ext.jextfrom, xgc->data.set.ss.ext.iextto, xgc->data.set.ss.ext.jextto,
                        xgc->data.set.ss.ext.shift, xgc->data.set.ss.ext.filebase_ex, xgc->data.set.ss.ext.filebase_ey, xgc->data.set.ss.ext.filebase_ez,
                        xgc->data.set.ss.ext.filebase_hx, xgc->data.set.ss.ext.filebase_hy, xgc->data.set.ss.ext.filebase_hz);
            }

            //if (xgc->data.set.ss.sf.filename) 
            if (xgc->data.set.ss.sf.valid) {
                fprintf(fw, "SOURCE_SF\n%g %g %g %d ", xgc->data.set.ss.sf.ia_theta, xgc->data.set.ss.sf.ia_phi, xgc->data.set.ss.sf.ia_psi, xgc->data.set.ss.sf.source_mode);
                if (xgc->data.set.ss.sf.source_mode == 0)
                    fprintf(fw, "%s\n\n", xgc->data.set.ss.sf.source_filename);
                else if (xgc->data.set.ss.sf.source_mode == 1)
                    fprintf(fw, "%g %g\n\n", xgc->data.set.ss.sf.source_wl, xgc->data.set.ss.sf.source_amplitude);
                else if (xgc->data.set.ss.sf.source_mode == 2)
                    fprintf(fw, "%g %g %g\n\n", xgc->data.set.ss.sf.source_wl, xgc->data.set.ss.sf.source_pulsewidth, xgc->data.set.ss.sf.source_amplitude);
            }

            //if ((xgc->data.set.ss.tsf.i0 + xgc->data.set.ss.tsf.j0 + xgc->data.set.ss.tsf.k0 + xgc->data.set.ss.tsf.i1 + xgc->data.set.ss.tsf.j1 + xgc->data.set.ss.tsf.k1)!=0) {
            if (xgc->data.set.ss.tsf.valid) {
                fprintf(fw, "SOURCE_TSF\n%d %d %d %d %d %d %g %g %g %d ",
                        xgc->data.set.ss.tsf.box_i0, xgc->data.set.ss.tsf.box_j0, xgc->data.set.ss.tsf.box_k0,
                        xgc->data.set.ss.tsf.box_in, xgc->data.set.ss.tsf.box_jn, xgc->data.set.ss.tsf.box_kn,
                        xgc->data.set.ss.tsf.ia_theta, xgc->data.set.ss.tsf.ia_phi, xgc->data.set.ss.tsf.ia_psi, xgc->data.set.ss.tsf.source_mode);
                if (xgc->data.set.ss.tsf.source_mode == 0)
                    fprintf(fw, "%s\n\n", xgc->data.set.ss.tsf.source_filename);
                else if (xgc->data.set.ss.tsf.source_mode == 1)
                    fprintf(fw, "%g %g\n\n", xgc->data.set.ss.tsf.source_wl, xgc->data.set.ss.tsf.source_amplitude);
                else if (xgc->data.set.ss.tsf.source_mode == 2)
                    fprintf(fw, "%g %d %g\n\n", xgc->data.set.ss.tsf.source_wl, xgc->data.set.ss.tsf.source_pulsewidth, xgc->data.set.ss.tsf.source_amplitude);
                if (xgc->data.set.ss.tsf.box_boundary_skipi0 == 1)
                    fprintf(fw, "TSF_SKIP\ni0\n\n");
                if (xgc->data.set.ss.tsf.box_boundary_skipin == 1)
                    fprintf(fw, "TSF_SKIP\nin\n\n");
                if (xgc->data.set.ss.tsf.box_boundary_skipj0 == 1)
                    fprintf(fw, "TSF_SKIP\nj0\n\n");
                if (xgc->data.set.ss.tsf.box_boundary_skipjn == 1)
                    fprintf(fw, "TSF_SKIP\njn\n\n");
                if (xgc->data.set.ss.tsf.box_boundary_skipk0 == 1)
                    fprintf(fw, "TSF_SKIP\nk0\n\n");
                if (xgc->data.set.ss.tsf.box_boundary_skipkn == 1)
                    fprintf(fw, "TSF_SKIP\nkn\n\n");
                if (xgc->data.set.ss.tsf.box_boundary_skipdepth != -1)
                    fprintf(fw, "TSF_SKIP\ndepth %d\n\n", xgc->data.set.ss.tsf.box_boundary_skipdepth);

                if (xgc->data.set.ss.tsf.gaussian == 1) {
                    fprintf(fw, "TSF_GAUSSIAN_Z\n%g %g %g %g\n\n",
                            xgc->data.set.ss.tsf.gaussian_fxpos, xgc->data.set.ss.tsf.gaussian_fypos,
                            xgc->data.set.ss.tsf.gaussian_rx, xgc->data.set.ss.tsf.gaussian_ry);
                }

                if (xgc->data.set.ss.tsf.radial == 1) {
                    fprintf(fw, "TSF_RADIAL_Z\n%g %g %g %g\n\n",
                            xgc->data.set.ss.tsf.radial_fxpos, xgc->data.set.ss.tsf.radial_fypos,
                            xgc->data.set.ss.tsf.radial_rx, xgc->data.set.ss.tsf.radial_ry);
                }

                if (xgc->data.set.ss.tsf.fiber == 1) {
                    fprintf(fw, "TSF_FIBER_Z\n%g %g %g %g %g %g\n\n",
                            xgc->data.set.ss.tsf.fiber_fxpos, xgc->data.set.ss.tsf.fiber_fypos,
                            xgc->data.set.ss.tsf.fiber_radius, xgc->data.set.ss.tsf.fiber_cutoff,
                            xgc->data.set.ss.tsf.fiber_epsilon_core, xgc->data.set.ss.tsf.fiber_epsilon_cladding);
                }

            }

            //if ((xgc->data.set.ss.ltsf.i0 + xgc->data.set.ss.ltsf.j0 + xgc->data.set.ss.ltsf.k0 + xgc->data.set.ss.ltsf.i1 + xgc->data.set.ss.ltsf.j1 + xgc->data.set.ss.ltsf.k1)!=0) {
            if (xgc->data.set.ss.ltsf.valid) {
                fprintf(fw, "SOURCE_LTSF\n%d %d %d %d %d %d %g %g %g %d ",
                        xgc->data.set.ss.ltsf.box_i0, xgc->data.set.ss.ltsf.box_j0, xgc->data.set.ss.ltsf.box_k0,
                        xgc->data.set.ss.ltsf.box_in, xgc->data.set.ss.ltsf.box_jn, xgc->data.set.ss.ltsf.box_kn,
                        xgc->data.set.ss.ltsf.ia_theta, xgc->data.set.ss.ltsf.ia_phi, xgc->data.set.ss.ltsf.ia_psi,
                        xgc->data.set.ss.ltsf.layered_count);

                for (i = 0; i < xgc->data.set.ss.ltsf.layered_count; i++) {
                    if (xgc->data.set.ss.ltsf.layered_material[i] == NULL || strcmp(xgc->data.set.ss.ltsf.layered_material[i], SOURCE_LAYERED_MATERIAL) == 0)
                        fprintf(fw, "%d %g %g %g %g ", xgc->data.set.ss.ltsf.layered_zpos[i], xgc->data.set.ss.ltsf.layered_epsilon[i], xgc->data.set.ss.ltsf.layered_mu[i], xgc->data.set.ss.ltsf.layered_sigma[i], xgc->data.set.ss.ltsf.layered_sigast[i]);
                    else
                        fprintf(fw, "%d 99 %s ", xgc->data.set.ss.ltsf.layered_zpos[i], xgc->data.set.ss.ltsf.layered_material[i]);
                }
                if (xgc->data.set.ss.ltsf.source_mode == 0)
                    fprintf(fw, "%d %s\n\n", xgc->data.set.ss.ltsf.source_mode, xgc->data.set.ss.ltsf.source_filename);
                else if (xgc->data.set.ss.ltsf.source_mode == 1)
                    fprintf(fw, "%d %g %g\n\n", xgc->data.set.ss.ltsf.source_mode, xgc->data.set.ss.ltsf.source_wl, xgc->data.set.ss.ltsf.source_amplitude);
                else if (xgc->data.set.ss.ltsf.source_mode == 2)
                    fprintf(fw, "%d %g %g %g\n\n", xgc->data.set.ss.ltsf.source_mode, xgc->data.set.ss.ltsf.source_wl, xgc->data.set.ss.ltsf.source_pulsewidth, xgc->data.set.ss.ltsf.source_amplitude);

                if (xgc->data.set.ss.ltsf.box_boundary_skipi0 == 1)
                    fprintf(fw, "LTSF_SKIP\ni0\n\n");
                if (xgc->data.set.ss.ltsf.box_boundary_skipin == 1)
                    fprintf(fw, "LTSF_SKIP\nin\n\n");
                if (xgc->data.set.ss.ltsf.box_boundary_skipj0 == 1)
                    fprintf(fw, "LTSF_SKIP\nj0\n\n");
                if (xgc->data.set.ss.ltsf.box_boundary_skipjn == 1)
                    fprintf(fw, "LTSF_SKIP\njn\n\n");
                if (xgc->data.set.ss.ltsf.box_boundary_skipk0 == 1)
                    fprintf(fw, "LTSF_SKIP\nk0\n\n");
                if (xgc->data.set.ss.ltsf.box_boundary_skipkn == 1)
                    fprintf(fw, "LTSF_SKIP\nkn\n\n");
                if (xgc->data.set.ss.ltsf.box_boundary_skipdepth != -1)
                    fprintf(fw, "LTSF_SKIP\ndepth %d\n\n", xgc->data.set.ss.ltsf.box_boundary_skipdepth);

                if (xgc->data.set.ss.ltsf.gaussian == 1) {
                    fprintf(fw, "LTSF_GAUSSIAN\n%g %g %g %g\n\n",
                            xgc->data.set.ss.ltsf.gaussian_fxpos, xgc->data.set.ss.ltsf.gaussian_fypos,
                            xgc->data.set.ss.ltsf.gaussian_rx, xgc->data.set.ss.ltsf.gaussian_ry);
                }

                if (xgc->data.set.ss.ltsf.radial == 1) {
                    fprintf(fw, "LTSF_RADIAL\n%g %g %g %g\n\n",
                            xgc->data.set.ss.ltsf.radial_fxpos, xgc->data.set.ss.ltsf.radial_fypos,
                            xgc->data.set.ss.ltsf.radial_rx, xgc->data.set.ss.ltsf.radial_ry);
                }

                if (xgc->data.set.ss.ltsf.fiber == 1) {
                    fprintf(fw, "LTSF_FIBER\n%g %g %g %g %g %g\n\n",
                            xgc->data.set.ss.ltsf.fiber_fxpos, xgc->data.set.ss.ltsf.fiber_fypos,
                            xgc->data.set.ss.ltsf.fiber_radius, xgc->data.set.ss.ltsf.fiber_cutoff,
                            xgc->data.set.ss.ltsf.fiber_epsilon_core, xgc->data.set.ss.ltsf.fiber_epsilon_cladding);
                }
            }

            //if ((xgc->data.set.ss.tsff.i0 + xgc->data.set.ss.tsff.j0 + xgc->data.set.ss.tsff.k0 + xgc->data.set.ss.tsff.i1 + xgc->data.set.ss.tsff.j1 + xgc->data.set.ss.tsff.k1)!=0) {
            if (xgc->data.set.ss.tsff.valid) {
                fprintf(fw, "SOURCE_TSFF\n%d %d %d %d %d %d %g %g %g %d %d %d ",
                        xgc->data.set.ss.tsff.box_i0, xgc->data.set.ss.tsff.box_j0, xgc->data.set.ss.tsff.box_k0,
                        xgc->data.set.ss.tsff.box_in, xgc->data.set.ss.tsff.box_jn, xgc->data.set.ss.tsff.box_kn,
                        xgc->data.set.ss.tsff.focused_thetamax, xgc->data.set.ss.tsff.focused_fdist, xgc->data.set.ss.tsff.focused_pol, xgc->data.set.ss.tsff.focused_nip, xgc->data.set.ss.tsff.focused_mip,
                        xgc->data.set.ss.tsff.source_mode);
                if (xgc->data.set.ss.tsff.source_mode == 0)
                    fprintf(fw, "%s\n\n", xgc->data.set.ss.tsff.source_filename);
                else if (xgc->data.set.ss.tsff.source_mode == 1)
                    fprintf(fw, "%g %g\n\n", xgc->data.set.ss.tsff.source_wl, xgc->data.set.ss.tsff.source_amplitude);
                else if (xgc->data.set.ss.tsff.source_mode == 2)
                    fprintf(fw, "%g %g %g\n\n", xgc->data.set.ss.tsff.source_wl, xgc->data.set.ss.tsff.source_pulsewidth, xgc->data.set.ss.tsff.source_amplitude);
                if (xgc->data.set.ss.tsff.box_boundary_skipi0 == 1)
                    fprintf(fw, "TSFF_SKIP\ni0\n\n");
                if (xgc->data.set.ss.tsff.box_boundary_skipin == 1)
                    fprintf(fw, "TSFF_SKIP\nin\n\n");
                if (xgc->data.set.ss.tsff.box_boundary_skipj0 == 1)
                    fprintf(fw, "TSFF_SKIP\nj0\n\n");
                if (xgc->data.set.ss.tsff.box_boundary_skipjn == 1)
                    fprintf(fw, "TSFF_SKIP\njn\n\n");
                if (xgc->data.set.ss.tsff.box_boundary_skipk0 == 1)
                    fprintf(fw, "TSFF_SKIP\nk0\n\n");
                if (xgc->data.set.ss.tsff.box_boundary_skipkn == 1)
                    fprintf(fw, "TSFF_SKIP\nkn\n\n");
                if (xgc->data.set.ss.tsff.box_boundary_skipdepth != -1)
                    fprintf(fw, "TSFF_SKIP\ndepth %d\n\n", xgc->data.set.ss.tsff.box_boundary_skipdepth);
            }

            //if ((xgc->data.set.ss.ltsff.i0 + xgc->data.set.ss.ltsff.j0 + xgc->data.set.ss.ltsff.k0 + xgc->data.set.ss.ltsff.i1 + xgc->data.set.ss.ltsff.j1 + xgc->data.set.ss.ltsff.k1)!=0) {
            if (xgc->data.set.ss.ltsff.valid) {
                fprintf(fw, "SOURCE_LTSFF\n%d %d %d %d %d %d %g %g %f %d %d %d ",
                        xgc->data.set.ss.ltsff.box_i0, xgc->data.set.ss.ltsff.box_j0, xgc->data.set.ss.ltsff.box_k0,
                        xgc->data.set.ss.ltsff.box_in, xgc->data.set.ss.ltsff.box_jn, xgc->data.set.ss.ltsff.box_kn,
                        xgc->data.set.ss.ltsff.focused_thetamax, xgc->data.set.ss.ltsff.focused_fdist, xgc->data.set.ss.ltsff.focused_pol, xgc->data.set.ss.ltsff.focused_nip, xgc->data.set.ss.ltsff.focused_mip,
                        xgc->data.set.ss.ltsff.layered_count);

                
                for (i = 0; i < xgc->data.set.ss.ltsff.layered_count; i++) {
                    if (xgc->data.set.ss.ltsff.layered_material[i] == NULL || strcmp(xgc->data.set.ss.ltsff.layered_material[i], SOURCE_LAYERED_MATERIAL) == 0)
                        fprintf(fw, "%d %g %g %g %g ", xgc->data.set.ss.ltsff.layered_zpos[i], xgc->data.set.ss.ltsff.layered_epsilon[i], xgc->data.set.ss.ltsff.layered_mu[i], xgc->data.set.ss.ltsff.layered_sigma[i], xgc->data.set.ss.ltsff.layered_sigast[i]);
                    else
                        fprintf(fw, "%d 99 %s ", xgc->data.set.ss.ltsff.layered_zpos[i], xgc->data.set.ss.ltsff.layered_material[i]);
                }
                if (xgc->data.set.ss.ltsff.source_mode == 0)
                    fprintf(fw, "%d %s\n\n", xgc->data.set.ss.ltsff.source_mode, xgc->data.set.ss.ltsff.source_filename);
                else if (xgc->data.set.ss.ltsff.source_mode == 1)
                    fprintf(fw, "%d %g %g\n\n", xgc->data.set.ss.ltsff.source_mode, xgc->data.set.ss.ltsff.source_wl, xgc->data.set.ss.ltsff.source_amplitude);
                else if (xgc->data.set.ss.ltsff.source_mode == 2)
                    fprintf(fw, "%d %g %g %g\n\n", xgc->data.set.ss.ltsff.source_mode, xgc->data.set.ss.ltsff.source_wl, xgc->data.set.ss.ltsff.source_pulsewidth, xgc->data.set.ss.ltsff.source_amplitude);

                if (xgc->data.set.ss.ltsff.box_boundary_skipi0 == 1)
                    fprintf(fw, "LTSFF_SKIP\ni0\n\n");
                if (xgc->data.set.ss.ltsff.box_boundary_skipin == 1)
                    fprintf(fw, "LTSFF_SKIP\nin\n\n");
                if (xgc->data.set.ss.ltsff.box_boundary_skipj0 == 1)
                    fprintf(fw, "LTSFF_SKIP\nj0\n\n");
                if (xgc->data.set.ss.ltsff.box_boundary_skipjn == 1)
                    fprintf(fw, "LTSFF_SKIP\njn\n\n");
                if (xgc->data.set.ss.ltsff.box_boundary_skipk0 == 1)
                    fprintf(fw, "LTSFF_SKIP\nk0\n\n");
                if (xgc->data.set.ss.ltsff.box_boundary_skipkn == 1)
                    fprintf(fw, "LTSFF_SKIP\nkn\n\n");
                if (xgc->data.set.ss.ltsff.box_boundary_skipdepth != -1)
                    fprintf(fw, "LTSFF_SKIP\ndepth %d\n\n", xgc->data.set.ss.ltsff.box_boundary_skipdepth);
            }
            if ((xgc->data.set.sb.bx0 == xgc->data.set.sb.bxn) && (xgc->data.set.sb.bx0 == xgc->data.set.sb.by0) && (xgc->data.set.sb.bx0 == xgc->data.set.sb.byn) &&
                (xgc->data.set.sb.bx0 == xgc->data.set.sb.bz0) && (xgc->data.set.sb.bx0 == xgc->data.set.sb.bzn) && xgc->data.set.sb.bx0 != 3) {
                if (xgc->data.set.sb.bx0 == 0)
                    fprintf(fw, "BOUNDARY_ALL\nnone\n\n");
                else if (xgc->data.set.sb.bx0 == 1)
                    fprintf(fw, "BOUNDARY_ALL\npec\n\n");
                else if (xgc->data.set.sb.bx0 == 2)
                    fprintf(fw, "BOUNDARY_ALL\nliao\n\n");
                else if (xgc->data.set.sb.bx0 == 3)
                    fprintf(fw, "BOUNDARY_ALL\ncpml %d %d %g %g %g\n\n", xgc->data.set.sb.depth_bx0, xgc->data.set.sb.m_bx0, xgc->data.set.sb.sigma_bx0, xgc->data.set.sb.a_bx0, xgc->data.set.sb.kappa_bx0);

            } else {
                if (xgc->data.set.sb.bx0 == 0)
                    fprintf(fw, "BOUNDARY_X0\nnone\n\n");
                else if (xgc->data.set.sb.bx0 == 1)
                    fprintf(fw, "BOUNDARY_X0\npec\n\n");
                else if (xgc->data.set.sb.bx0 == 2)
                    fprintf(fw, "BOUNDARY_X0\nliao\n\n");
                else if (xgc->data.set.sb.bx0 == 3)
                    fprintf(fw, "BOUNDARY_X0\ncpml %d %d %g %g %g\n\n", xgc->data.set.sb.depth_bx0, xgc->data.set.sb.m_bx0, xgc->data.set.sb.sigma_bx0, xgc->data.set.sb.a_bx0, xgc->data.set.sb.kappa_bx0);

                if (xgc->data.set.sb.bxn == 0)
                    fprintf(fw, "BOUNDARY_XN\nnone\n\n");
                else if (xgc->data.set.sb.bxn == 1)
                    fprintf(fw, "BOUNDARY_XN\npec\n\n");
                else if (xgc->data.set.sb.bxn == 2)
                    fprintf(fw, "BOUNDARY_XN\nliao\n\n");
                else if (xgc->data.set.sb.bxn == 3)
                    fprintf(fw, "BOUNDARY_XN\ncpml %d %d %g %g %g\n\n", xgc->data.set.sb.depth_bxn, xgc->data.set.sb.m_bxn, xgc->data.set.sb.sigma_bxn, xgc->data.set.sb.a_bxn, xgc->data.set.sb.kappa_bxn);

                if (xgc->data.set.sb.by0 == 0)
                    fprintf(fw, "BOUNDARY_Y0\nnone\n\n");
                else if (xgc->data.set.sb.by0 == 1)
                    fprintf(fw, "BOUNDARY_Y0\npec\n\n");
                else if (xgc->data.set.sb.by0 == 2)
                    fprintf(fw, "BOUNDARY_Y0\nliao\n\n");
                else if (xgc->data.set.sb.by0 == 3)
                    fprintf(fw, "BOUNDARY_Y0\ncpml %d %d %g %g %g\n\n", xgc->data.set.sb.depth_by0, xgc->data.set.sb.m_by0, xgc->data.set.sb.sigma_by0, xgc->data.set.sb.a_by0, xgc->data.set.sb.kappa_by0);

                if (xgc->data.set.sb.byn == 0)
                    fprintf(fw, "BOUNDARY_YN\nnone\n\n");
                else if (xgc->data.set.sb.byn == 1)
                    fprintf(fw, "BOUNDARY_YN\npec\n\n");
                else if (xgc->data.set.sb.byn == 2)
                    fprintf(fw, "BOUNDARY_YN\nliao\n\n");
                else if (xgc->data.set.sb.byn == 3)
                    fprintf(fw, "BOUNDARY_YN\ncpml %d %d %g %g %g\n\n", xgc->data.set.sb.depth_byn, xgc->data.set.sb.m_byn, xgc->data.set.sb.sigma_byn, xgc->data.set.sb.a_byn, xgc->data.set.sb.kappa_byn);

                if (xgc->data.set.sb.bz0 == 0)
                    fprintf(fw, "BOUNDARY_Z0\nnone\n\n");
                else if (xgc->data.set.sb.bz0 == 1)
                    fprintf(fw, "BOUNDARY_Z0\npec\n\n");
                else if (xgc->data.set.sb.bz0 == 2)
                    fprintf(fw, "BOUNDARY_Z0\nliao\n\n");
                else if (xgc->data.set.sb.bz0 == 3)
                    fprintf(fw, "BOUNDARY_Z0\ncpml %d %d %g %g %g\n\n", xgc->data.set.sb.depth_bz0, xgc->data.set.sb.m_bz0, xgc->data.set.sb.sigma_bz0, xgc->data.set.sb.a_bz0, xgc->data.set.sb.kappa_bz0);

                if (xgc->data.set.sb.bzn == 0)
                    fprintf(fw, "BOUNDARY_ZN\nnone\n\n");
                else if (xgc->data.set.sb.bzn == 1)
                    fprintf(fw, "BOUNDARY_ZN\npec\n\n");
                else if (xgc->data.set.sb.bzn == 2)
                    fprintf(fw, "BOUNDARY_ZN\nliao\n\n");
                else if (xgc->data.set.sb.bzn == 3)
                    fprintf(fw, "BOUNDARY_ZN\ncpml %d %d %g %g %g\n\n", xgc->data.set.sb.depth_bzn, xgc->data.set.sb.m_bzn, xgc->data.set.sb.sigma_bzn, xgc->data.set.sb.a_bzn, xgc->data.set.sb.kappa_bzn);
            }
            if (xgc->data.set.smb.bx0 == SV_BOUNDARY_PERIODIC)
                fprintf(fw, "MBOUNDARY_X0\nperiodic %d\n\n", xgc->data.set.smb.bx0pos);
            if (xgc->data.set.smb.bxn == SV_BOUNDARY_PERIODIC)
                fprintf(fw, "MBOUNDARY_XN\nperiodic %d\n\n", xgc->data.set.smb.bxnpos);

            if (xgc->data.set.smb.by0 == SV_BOUNDARY_PERIODIC)
                fprintf(fw, "MBOUNDARY_Y0\nperiodic %d\n\n", xgc->data.set.smb.by0pos);
            if (xgc->data.set.smb.byn == SV_BOUNDARY_PERIODIC)
                fprintf(fw, "MBOUNDARY_YN\nperiodic %d\n\n", xgc->data.set.smb.bynpos);

            if (xgc->data.set.smb.bz0 == SV_BOUNDARY_PERIODIC)
                fprintf(fw, "MBOUNDARY_Z0\nperiodic %d\n\n", xgc->data.set.smb.bz0pos);
            if (xgc->data.set.smb.bzn == SV_BOUNDARY_PERIODIC)
                fprintf(fw, "MBOUNDARY_ZN\nperiodic %d\n\n", xgc->data.set.smb.bznpos);

            if (xgc->data.set.sm.matmode_check)
                fprintf(fw, "MATMODE_CHECK\n1\n\n");
            if (xgc->data.set.sm.smooth)
                fprintf(fw, "MEDIUM_SMOOTH\n%d\n\n", xgc->data.set.sm.smooth);
            if (xgc->data.set.sm.in_voxel_filename)
                fprintf(fw, "MEDIUM_LINEAR\n%s\n\n", xgc->data.set.sm.in_voxel_filename);
            if (xgc->data.set.sm.in_vector_filename)
                fprintf(fw, "MEDIUM_VECTOR\n%s\n\n", xgc->data.set.sm.in_vector_filename);
            
            for (i = 0; i < xgc->data.set.sm.ngrowths; i++) {
                fprintf(fw, "MEDIUM_GROW\n%d %d %d %d %d %d %d %d %d %d %g %g %d\n\n",
                        xgc->data.set.sm.grow_i0[i],
                        xgc->data.set.sm.grow_j0[i],
                        xgc->data.set.sm.grow_k0[i],
                        xgc->data.set.sm.grow_in[i],
                        xgc->data.set.sm.grow_jn[i],
                        xgc->data.set.sm.grow_kn[i],
                        xgc->data.set.sm.grow_addindex[i],
                        xgc->data.set.sm.grow_attachindex[i],
                        xgc->data.set.sm.grow_subsampling[i],
                        xgc->data.set.sm.grow_nsteps[i],
                        xgc->data.set.sm.grow_mobility[i],
                        xgc->data.set.sm.grow_probability[i],
                        xgc->data.set.sm.grow_seed[i]);

                if (xgc->data.set.sm.grow_skipi0[i] == 1)
                    fprintf(fw, "MEDIUM_GROW_SKIP_FACE\ni0\n\n");
                if (xgc->data.set.sm.grow_skipin[i] == 1)
                    fprintf(fw, "MEDIUM_GROW_SKIP_FACE\nin\n\n");
                if (xgc->data.set.sm.grow_skipj0[i] == 1)
                    fprintf(fw, "MEDIUM_GROW_SKIP_FACE\nj0\n\n");
                if (xgc->data.set.sm.grow_skipjn[i] == 1)
                    fprintf(fw, "MEDIUM_GROW_SKIP_FACE\njn\n\n");
                if (xgc->data.set.sm.grow_skipk0[i] == 1)
                    fprintf(fw, "MEDIUM_GROW_SKIP_FACE\nk0\n\n");
                if (xgc->data.set.sm.grow_skipkn[i] == 1)
                    fprintf(fw, "MEDIUM_GROW_SKIP_FACE\nkn\n\n");
            }

            for (i = 0; i < xgc->data.set.sm.nroughens; i++) {
                fprintf(fw, "MEDIUM_ROUGHEN\n%d %d %d %g %d %d %d\n\n",
                        xgc->data.set.sm.rough_radius_peak[i], xgc->data.set.sm.rough_radius_span[i],
                        xgc->data.set.sm.rough_iterations[i], xgc->data.set.sm.rough_probability[i],
                        xgc->data.set.sm.rough_matindex[i], xgc->data.set.sm.rough_voidindex[i],
                        xgc->data.set.sm.rough_seed[i]);
            }

            for (i = 0; i < xgc->data.set.sm.nspectrals; i++) {
                fprintf(fw, "MEDIUM_SPECTRAL\n%g %g %d %d\n\n",
                        xgc->data.set.sm.spectral_sigma[i], xgc->data.set.sm.spectral_t[i],
                        xgc->data.set.sm.spectral_matindex[i], xgc->data.set.sm.spectral_seed[i]);
            }

            for (i = 0; i < xgc->data.set.sm.nexprs; i++) {
                fprintf(fw, "MEDIUM_EXPRESSION\n%d %d %d %d %d %d %d %d %d %d %s\n\n",
                        xgc->data.set.sm.expr_i0[i], xgc->data.set.sm.expr_j0[i], xgc->data.set.sm.expr_k0[i],
                        xgc->data.set.sm.expr_in[i], xgc->data.set.sm.expr_jn[i], xgc->data.set.sm.expr_kn[i],
                        xgc->data.set.sm.expr_matindex[i], xgc->data.set.sm.expr_voidindex[i],
                        xgc->data.set.sm.expr_maxdist[i], xgc->data.set.sm.expr_distmode[i],
                        xgc->data.set.sm.expr_expr[i]);
            }

            if (xgc->data.set.so.outfile)
                fprintf(fw, "OUT_FILE\n%s\n\n", xgc->data.set.so.outfile);

            for (i = 0; i < xgc->data.set.so.npnts; i++)
                fprintf(fw, "OUT_POINT\n%s %d %d %d %d %s\n\n",
                        cpstring(xgc->data.set.so.pnts[i].component),
                        xgc->data.set.so.pnts[i].step, xgc->data.set.so.pnts[i].i, xgc->data.set.so.pnts[i].j, xgc->data.set.so.pnts[i].k,
                        xgc->data.set.so.pnts[i].filebase);

            for (i = 0; i < xgc->data.set.so.nplns; i++)
                fprintf(fw, "OUT_PLANE\n%s %d %d %d %d %d %d %d %s\n\n",
                        cpstring(xgc->data.set.so.plns[i].component),
                        xgc->data.set.so.plns[i].step, xgc->data.set.so.plns[i].start, xgc->data.set.so.plns[i].stop, xgc->data.set.so.plns[i].format,
                        xgc->data.set.so.plns[i].i, xgc->data.set.so.plns[i].j, xgc->data.set.so.plns[i].k,
                        xgc->data.set.so.plns[i].filebase);

            for (i = 0; i < xgc->data.set.so.nimgs; i++)
                fprintf(fw, "OUT_IMAGE\n%s %d %d %d %d %s\n\n",
                        cpstring(xgc->data.set.so.imgs[i].component),
                        xgc->data.set.so.imgs[i].step, xgc->data.set.so.imgs[i].i, xgc->data.set.so.imgs[i].j, xgc->data.set.so.imgs[i].k,
                        xgc->data.set.so.imgs[i].filebase);

            for (i = 0; i < xgc->data.set.so.nsums; i++)
                if (xgc->data.set.so.sums[i].stringbased)
                    fprintf(fw, "OUT_SUMTAB\n%s %d %d %d %d %d %d %d %s %s\n\n",
                            scpstring(xgc->data.set.so.sums[i].component),
                            xgc->data.set.so.sums[i].step, xgc->data.set.so.sums[i].box_i0, xgc->data.set.so.sums[i].box_j0, xgc->data.set.so.sums[i].box_k0,
                            xgc->data.set.so.sums[i].box_in, xgc->data.set.so.sums[i].box_jn, xgc->data.set.so.sums[i].box_kn,
                            xgc->data.set.so.sums[i].string,
                            xgc->data.set.so.sums[i].filename);
                else
                    fprintf(fw, "OUT_SUM\n%s %d %d %d %d %d %d %d %g %g %g %g %s\n\n",
                            scpstring(xgc->data.set.so.sums[i].component),
                            xgc->data.set.so.sums[i].step, xgc->data.set.so.sums[i].box_i0, xgc->data.set.so.sums[i].box_j0, xgc->data.set.so.sums[i].box_k0,
                            xgc->data.set.so.sums[i].box_in, xgc->data.set.so.sums[i].box_jn, xgc->data.set.so.sums[i].box_kn,
                            xgc->data.set.so.sums[i].layered_epsilon, xgc->data.set.so.sums[i].layered_mu, xgc->data.set.so.sums[i].layered_sigma, xgc->data.set.so.sums[i].layered_sigast,
                            xgc->data.set.so.sums[i].filename);


            for (i = 0; i < xgc->data.set.so.ncubs; i++)
                fprintf(fw, "OUT_VOLUME\n%s %d %d %d %d %s\n\n",
                        vcpstring(xgc->data.set.so.cubs[i].component),
                        xgc->data.set.so.cubs[i].step, xgc->data.set.so.cubs[i].start,
                        xgc->data.set.so.cubs[i].stop,
                        xgc->data.set.so.cubs[i].format, xgc->data.set.so.cubs[i].filebase);


            for (i = 0; i < xgc->data.set.so.nforces; i++)
                fprintf(fw, "OUT_FORCE\n%d %d %d %d %d %d %d %s\n\n",
                        xgc->data.set.so.forces[i].step, xgc->data.set.so.forces[i].box_i0, xgc->data.set.so.forces[i].box_j0, xgc->data.set.so.forces[i].box_k0,
                        xgc->data.set.so.forces[i].box_in, xgc->data.set.so.forces[i].box_jn, xgc->data.set.so.forces[i].box_kn,
                        xgc->data.set.so.forces[i].filename);

            if (xgc->data.set.sf.nrs || xgc->data.set.sf.nareas) {
                fprintf(fw, "NFFF\n%d %d %d %d %d %d 0 100000\n\n", xgc->data.set.sf.box_i0, xgc->data.set.sf.box_j0, xgc->data.set.sf.box_k0, xgc->data.set.sf.box_in, xgc->data.set.sf.box_jn, xgc->data.set.sf.box_kn);

                if (xgc->data.set.sf.box_boundary_skipi0)
                    fprintf(fw, "NFFF_SKIP\ni0 %d %d %d %d\n\n", xgc->data.set.sf.skipi0_jmin, xgc->data.set.sf.skipi0_kmin, xgc->data.set.sf.skipi0_jmax, xgc->data.set.sf.skipi0_kmax);
                if (xgc->data.set.sf.box_boundary_skipin)
                    fprintf(fw, "NFFF_SKIP\nin %d %d %d %d\n\n", xgc->data.set.sf.skipin_jmin, xgc->data.set.sf.skipin_kmin, xgc->data.set.sf.skipin_jmax, xgc->data.set.sf.skipin_kmax);
                if (xgc->data.set.sf.box_boundary_skipj0)
                    fprintf(fw, "NFFF_SKIP\nj0 %d %d %d %d\n\n", xgc->data.set.sf.skipj0_imin, xgc->data.set.sf.skipj0_kmin, xgc->data.set.sf.skipj0_imax, xgc->data.set.sf.skipj0_kmax);
                if (xgc->data.set.sf.box_boundary_skipjn)
                    fprintf(fw, "NFFF_SKIP\njn %d %d %d %d\n\n", xgc->data.set.sf.skipjn_imin, xgc->data.set.sf.skipjn_kmin, xgc->data.set.sf.skipjn_imax, xgc->data.set.sf.skipjn_kmax);
                if (xgc->data.set.sf.box_boundary_skipk0)
                    fprintf(fw, "NFFF_SKIP\nk0 %d %d %d %d\n\n", xgc->data.set.sf.skipk0_imin, xgc->data.set.sf.skipk0_jmin, xgc->data.set.sf.skipk0_imax, xgc->data.set.sf.skipk0_jmax);
                if (xgc->data.set.sf.box_boundary_skipkn)
                    fprintf(fw, "NFFF_SKIP\nkn %d %d %d %d\n\n", xgc->data.set.sf.skipkn_imin, xgc->data.set.sf.skipkn_jmin, xgc->data.set.sf.skipkn_imax, xgc->data.set.sf.skipkn_jmax);

                for (i = 0; i < xgc->data.set.sf.nrs; i++) {
                    if (xgc->data.set.sf.source_filename[i] && !xgc->data.set.sf.individual[i])
                        fprintf(fw, "NFFF_RAMAHI_POINT\n%d %d %d %s\n\n", xgc->data.set.sf.ri[i], xgc->data.set.sf.rj[i], xgc->data.set.sf.rk[i], xgc->data.set.sf.source_filename[i]);
                }

                for (i = 0; i < xgc->data.set.sf.nareas; i++)
                    fprintf(fw, "NFFF_SPHERICAL_AREA\n%d %d %d %g %g %g %g %d\n\n",
                            xgc->data.set.sf.area_thetares[i], xgc->data.set.sf.area_phires[i],
                            xgc->data.set.sf.area_radius[i],
                            xgc->data.set.sf.area_thetafrom[i], xgc->data.set.sf.area_phifrom[i],
                            xgc->data.set.sf.area_thetato[i], xgc->data.set.sf.area_phito[i],
                            xgc->data.set.sf.area_savefile[i]);
            }

            if (xgc->data.set.spf.nrs || xgc->data.set.spf.nareas) {
                fprintf(fw, "PERIODIC_NFFF\n%d %d %d %d %d %d %d %d %d %d\n\n",
                        xgc->data.set.spf.box_i0, xgc->data.set.spf.box_j0, xgc->data.set.spf.box_k0,
                        xgc->data.set.spf.box_in, xgc->data.set.spf.box_jn, xgc->data.set.spf.box_kn,
                        xgc->data.set.spf.pimin, xgc->data.set.spf.pjmin,
                        xgc->data.set.spf.pimax, xgc->data.set.spf.pjmax);

                if (xgc->data.set.spf.box_boundary_skipk0)
                    fprintf(fw, "PERIODIC_NFFF_SKIP\nk0 %d %d %d %d\n\n", xgc->data.set.spf.skipk0_imin, xgc->data.set.spf.skipk0_jmin, xgc->data.set.spf.skipk0_imax, xgc->data.set.spf.skipk0_jmax);
                if (xgc->data.set.spf.box_boundary_skipkn)
                    fprintf(fw, "PERIODIC_NFFF_SKIP\nkn %d %d %d %d\n\n", xgc->data.set.spf.skipkn_imin, xgc->data.set.spf.skipkn_jmin, xgc->data.set.spf.skipkn_imax, xgc->data.set.spf.skipkn_jmax);
                if (xgc->data.set.spf.postprocess)
                    fprintf(fw, "PERIODIC_NFFF_POSTPROCESS\n1 %d\n\n", xgc->data.set.spf.ppstart);

                for (i = 0; i < xgc->data.set.spf.nrs; i++) {
                    if (xgc->data.set.spf.source_filename[i] && !xgc->data.set.spf.individual[i])
                        fprintf(fw, "PERIODIC_NFFF_RAMAHI_POINT\n%d %d %d %s\n\n", xgc->data.set.spf.ri[i], xgc->data.set.spf.rj[i], xgc->data.set.spf.rk[i], xgc->data.set.spf.source_filename[i]);
                }

                for (i = 0; i < xgc->data.set.spf.nareas; i++) {
                    fprintf(fw, "PERIODIC_NFFF_SPHERICAL_AREA\n%d %d %d %g %g %g %g %d\n\n",
                            xgc->data.set.spf.area_thetares[i], xgc->data.set.spf.area_phires[i],
                            xgc->data.set.spf.area_radius[i],
                            xgc->data.set.spf.area_thetafrom[i], xgc->data.set.spf.area_phifrom[i],
                            xgc->data.set.spf.area_thetato[i], xgc->data.set.spf.area_phito[i],
                            xgc->data.set.spf.area_savefile[i]);
                }
            }

            fclose(fw);

            if (g_file_get_contents(xgc->tmpfilename, &contents, NULL, NULL)) {
                gtk_text_buffer_set_text(xgc->tb_par, contents, -1);
                return TRUE;
            }
        }

        //g_free(filename);
    }

    return FALSE;
} /* write_parfile() */

/* support function */
void
writemat(FILE *fw, SvMatProp mat)
{
    if (mat.type == SV_MAT_LINEAR || mat.type == SV_MAT_LINTAB) {
        fprintf(fw, "%d %g %g %g %g\n", mat.type, mat.epsilon, mat.mu, mat.sigma, mat.sigast);
    } else if (mat.type == SV_MAT_PEC) {
        fprintf(fw, "%d\n", mat.type);
    } else if (mat.type == SV_MAT_CP) {
        fprintf(fw, "%d %g %g %g %g %g %g %g %g %g %g %g\n", mat.type, mat.epsilon, mat.drude_omega_p, mat.drude_nu,
                mat.cp3_a[0], mat.cp3_phi[0], mat.cp3_omega[0], mat.cp3_gamma[0], mat.cp3_a[1], mat.cp3_phi[1], mat.cp3_omega[1], mat.cp3_gamma[1]);
    } else if (mat.type == SV_MAT_PLRC) {
        fprintf(fw, "%d %g 0 %g %g %g %g %g %g %g %g %g %g\n", mat.type, mat.epsilon, mat.drude_omega_p, mat.drude_nu,
                mat.cp3_a[0], mat.cp3_phi[0], mat.cp3_omega[0], mat.cp3_gamma[0], mat.cp3_a[1], mat.cp3_phi[1], mat.cp3_omega[1], mat.cp3_gamma[1]);
    } else if (mat.type == SV_MAT_ADE) {
        fprintf(fw, "%d %g 0 %g %g %g %g %g %g %g %g %g %g\n", mat.type, mat.epsilon, mat.drude_omega_p, mat.drude_nu,
                mat.cp3_a[0], mat.cp3_phi[0], mat.cp3_omega[0], mat.cp3_gamma[0], mat.cp3_a[1], mat.cp3_phi[1], mat.cp3_omega[1], mat.cp3_gamma[1]);
    }
} /* writemat() */


/* write data to material file */
gboolean
write_matfile(XGControls *xgc)
{
    FILE *fw;
    GtkWidget *dialog;
    gchar *filename = NULL;
    GError* glib_error = NULL;
    gchar *contents;
    SvSphere sx;
    SvVoxel vx;
    SvCone cnx;
    SvCylinder cx;
    SvRCone rcnx;
    SvTetrahedron tx;
    SvGwydd gx;
    gint i, m, npos, nskips, ngtot = xgc->data.set_mat.ngtot, nhandled;
    gboolean overriden;
    gint tetcor;


    if (xgc->mat_file_success == FALSE) {
        printf("Error parsing material file. Material file will not be saved.\n");
    }

    if (NULL == xgc->tmpfilename)
        xgc->tmpfilename = get_temporary_par_filename(xgc);    

    if (NULL == xgc->tmpmatfilename)
        xgc->tmpmatfilename = get_temporary_mat_filename(xgc);

    // printf("write parfile will be writing to %s assigning it to temporary location %s\n", filename, xgc->tmpfilename);

    if (xgc->tmpmatfilename) {
        fw = fopen(xgc->tmpmatfilename, "w");
        if (!fw) {
            dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                                            GTK_DIALOG_DESTROY_WITH_PARENT,
                                            GTK_MESSAGE_ERROR,
                                            GTK_BUTTONS_CLOSE,
                                            "Error opening temporary file: %s",
                                            g_strdup(glib_error->message));
            gtk_dialog_run(GTK_DIALOG(dialog));
            gtk_widget_destroy(dialog);
        } else {
            npos = nskips = nhandled = 0;

            tetcor = 0; //correction for tetrahedrons contained in tetgen sets
            for (m = 0; m < xgc->data.set_mat.ntetgens; m++)
                tetcor += xgc->data.set_mat.tetgen_end[m] - xgc->data.set_mat.tetgen_start[m];

            ngtot = xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + xgc->data.set_mat.cylinders->len
                    + xgc->data.set_mat.cones->len + xgc->data.set_mat.rcones->len + xgc->data.set_mat.tetrahedrons->len + xgc->data.set_mat.gwydds->len + xgc->data.set_mat.ntetgens - tetcor;

            do {
                for (m = 0; m < (gint)xgc->data.set_mat.spheres->len; m++) {
                    sx = g_array_index(xgc->data.set_mat.spheres, SvSphere, m);

                    if (npos != sx.n) {
                        nskips++;
                        continue;
                    } else {
                        nskips = 0;
                        nhandled++;
                        npos++;
                    }

                    fprintf(fw, "4 %g %g %g %g ", sx.pnt1[0], sx.pnt1[1], sx.pnt1[2], sx.radius);

                    overriden = FALSE;
                    for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                        if (xgc->data.set_mat.overpos[i] == sx.n) {
                            fprintf(fw, "99 %s\n", xgc->data.set_mat.materials[i]);
                            overriden = TRUE;
                        }
                    }
                    if (!overriden)
                        writemat(fw, sx.mat);
                }

                for (m = 0; m < (gint)xgc->data.set_mat.voxels->len; m++) {
                    vx = g_array_index(xgc->data.set_mat.voxels, SvVoxel, m);

                    if (npos != vx.n) {
                        nskips++;
                        continue;
                    } else {
                        nskips = 0;
                        nhandled++;
                        npos++;
                    }

                    fprintf(fw, "8 %g %g %g %g %g %g ", vx.pnt1[0], vx.pnt1[1], vx.pnt1[2], vx.pnt2[0], vx.pnt2[1], vx.pnt2[2]);

                    overriden = FALSE;
                    for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                        if (xgc->data.set_mat.overpos[i] == vx.n) {
                            fprintf(fw, "99 %s\n", xgc->data.set_mat.materials[i]);
                            overriden = TRUE;
                        }
                    }
                    if (!overriden)
                        writemat(fw, vx.mat);
                }

                for (m = 0; m < (gint)xgc->data.set_mat.cylinders->len; m++) {
                    cx = g_array_index(xgc->data.set_mat.cylinders, SvCylinder, m);

                    if (npos != cx.n) {
                        nskips++;
                        continue;
                    } else {
                        nskips = 0;
                        nhandled++;
                        npos++;
                    }

                    fprintf(fw, "7 %g %g %g %g %g %g %g ", cx.pnt1[0], cx.pnt1[1], cx.pnt1[2], cx.pnt2[0], cx.pnt2[1], cx.pnt2[2], cx.radius);

                    overriden = FALSE;
                    for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                        if (xgc->data.set_mat.overpos[i] == cx.n) {
                            fprintf(fw, "99 %s\n", xgc->data.set_mat.materials[i]);
                            overriden = TRUE;
                        }
                    }
                    if (!overriden)
                        writemat(fw, cx.mat);
                }

                for (m = 0; m < (gint)xgc->data.set_mat.cones->len; m++) {
                    cnx = g_array_index(xgc->data.set_mat.cones, SvCone, m);

                    if (npos != cnx.n) {
                        nskips++;
                        continue;
                    } else {
                        nskips = 0;
                        nhandled++;
                        npos++;
                    }

                    fprintf(fw, "107 %g %g %g %g %g %g %g ", cnx.pnt1[0], cnx.pnt1[1], cnx.pnt1[2], cnx.pnt2[0], cnx.pnt2[1], cnx.pnt2[2], cnx.radius);

                    overriden = FALSE;
                    for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                        if (xgc->data.set_mat.overpos[i] == cnx.n) {
                            fprintf(fw, "99 %s\n", xgc->data.set_mat.materials[i]);
                            overriden = TRUE;
                        }
                    }
                    if (!overriden)
                        writemat(fw, cnx.mat);

                }

                for (m = 0; m < (gint)xgc->data.set_mat.rcones->len; m++) {
                    rcnx = g_array_index(xgc->data.set_mat.rcones, SvRCone, m);

                    if (npos != rcnx.n) {
                        nskips++;
                        continue;
                    } else {
                        nskips = 0;
                        nhandled++;
                        npos++;
                    }

                    fprintf(fw, "108 %g %g %g %g %g %g %g %g ", rcnx.pnt1[0], rcnx.pnt1[1], rcnx.pnt1[2], rcnx.pnt2[0], rcnx.pnt2[1], rcnx.pnt2[2], rcnx.radius1, rcnx.radius2);

                    overriden = FALSE;
                    for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                        if (xgc->data.set_mat.overpos[i] == rcnx.n) {
                            fprintf(fw, "99 %s\n", xgc->data.set_mat.materials[i]);
                            overriden = TRUE;
                        }
                    }

                    if (!overriden)
                        writemat(fw, rcnx.mat);
                }

                for (m = 0; m < (gint)xgc->data.set_mat.tetrahedrons->len; m++) {
                    tx = g_array_index(xgc->data.set_mat.tetrahedrons, SvTetrahedron, m);
                    if (tx.setpart != 0)
                        continue;

                    if (npos != tx.n) {
                        nskips++;
                        continue;
                    } else {
                        nskips = 0;
                        nhandled++;
                        npos++;
                    }

                    fprintf(fw, "12 %g %g %g %g %g %g %g %g %g %g %g %g ", tx.pnt1[0], tx.pnt1[1], tx.pnt1[2], tx.pnt2[0], tx.pnt2[1], tx.pnt2[2], tx.pnt3[0], tx.pnt3[1], tx.pnt3[2], tx.pnt4[0], tx.pnt4[1], tx.pnt4[2]);

                    overriden = FALSE;
                    for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                        if (xgc->data.set_mat.overpos[i] == tx.n) {
                            fprintf(fw, "99 %s\n", xgc->data.set_mat.materials[i]);
                            overriden = TRUE;
                        }
                    }
                    if (!overriden)
                        writemat(fw, tx.mat);
                }

                for (m = 0; m < (gint)xgc->data.set_mat.gwydds->len; m++) {
                    gx = g_array_index(xgc->data.set_mat.gwydds, SvGwydd, m);

                    if (npos != gx.n) {
                        nskips++;
                        continue;
                    } else {
                        nskips = 0;
                        nhandled++;
                        npos++;
                    }

                    fprintf(fw, "22 %s %d %d %d %d %d %g %g %g %d ", gx.filebase, gx.channel, gx.mask, gx.i, gx.j, gx.k, gx.xoffset, gx.yoffset, gx.zoffset, gx.depth);
                    overriden = FALSE;
                    for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                        if (xgc->data.set_mat.overpos[i] == gx.n) {
                            fprintf(fw, "99 %s\n", xgc->data.set_mat.materials[i]);
                            overriden = TRUE;
                        }
                    }
                    if (!overriden)
                        writemat(fw, gx.mat);

                }
                for (m = 0; m < xgc->data.set_mat.ntetgens; m++) {
                    if (npos == xgc->data.set_mat.tetgen_start[m]) {
                        fprintf(fw, "21 %s %d %d %g %g %g %g %g %g ", xgc->data.set_mat.tetgen_filebase[m], xgc->data.set_mat.tetgen_attribute_pos[m], xgc->data.set_mat.tetgen_attribute_val[m],
                                xgc->data.set_mat.tetgen_xshift[m], xgc->data.set_mat.tetgen_yshift[m], xgc->data.set_mat.tetgen_zshift[m], xgc->data.set_mat.tetgen_xmult[m], xgc->data.set_mat.tetgen_ymult[m], xgc->data.set_mat.tetgen_zmult[m]);
                        overriden = FALSE;

                        for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                            if (xgc->data.set_mat.overpos[i] == xgc->data.set_mat.tetgen_start[m]) {
                                //printf("writing overriden mat %d %s\n", i, xgc->materials[i]);

                                fprintf(fw, "99 %s\n", xgc->data.set_mat.materials[i]);
                                overriden = TRUE;
                            }
                        }
                        if (!overriden) {
                            tx = g_array_index(xgc->data.set_mat.tetrahedrons, SvTetrahedron, xgc->data.set_mat.tetgen_start[m]);
                            writemat(fw, tx.mat);
                        }

                        npos = xgc->data.set_mat.tetgen_end[m];
                        //printf("next npos will be %d\n", npos);
                        nskips = 0;
                        nhandled++;
                    } else
                        nskips++;
                }

                // printf("nh %d  np %d  ns %d ntot %d\n", nhandled, npos, nskips, ngtot);

                if (nskips > 100) {
                    printf("Cannot find object %d, it was probably deleted\n", npos);
                    npos++;
                    nskips = 0;
                }   //this was probably a deleted item, don't search for it anymore
            } while (nhandled < ngtot && npos < (1e9));  //second condition (to prevent infinite loop at non-contiguous sets) should be modified

            fclose(fw);

            if (g_file_get_contents(xgc->tmpmatfilename, &contents, NULL, NULL)) {
                gtk_text_buffer_set_text(xgc->tb_mat, contents, -1);

                return TRUE;
            }
        }

        g_free(filename);
    }

    return FALSE;
} /* write_matfile() */

/* write text buffer to parfile */
gboolean
write_textbuffer_to_parfile(XGControls *xgc)
{
    GtkTextIter     start;
    GtkTextIter     end;
    GtkWidget       *dialog;
    GError          *err = NULL;
    gchar           *textbuf;

    if (NULL == xgc->tmpfilename)
        xgc->tmpfilename = get_temporary_par_filename(xgc);

    if (xgc->tmpfilename) {
        gtk_text_buffer_get_start_iter(xgc->tb_par, &start);
        gtk_text_buffer_get_end_iter(xgc->tb_par, &end);
        textbuf = gtk_text_buffer_get_text(xgc->tb_par, &start, &end, TRUE);
        if (!g_file_set_contents(xgc->tmpfilename, textbuf, -1, &err)) {
            dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                                            GTK_DIALOG_DESTROY_WITH_PARENT,
                                            GTK_MESSAGE_ERROR,
                                            GTK_BUTTONS_CLOSE,
                                            "Error writing to temporary file: %s",
                                            g_strdup(err->message));
            gtk_dialog_run(GTK_DIALOG(dialog));
            gtk_widget_destroy(dialog);

            return 0;
        }
    }

    return 1;
}

/* write text buffer to matfile */
gboolean
write_textbuffer_to_matfile(XGControls *xgc)
{
    GtkTextIter     start;
    GtkTextIter     end;
    GtkWidget       *dialog;
    GError          *err = NULL;
    gchar           *textbuf;

    if (NULL == xgc->tmpmatfilename)
        xgc->tmpmatfilename = get_temporary_mat_filename(xgc);

    if (xgc->tmpmatfilename) {
        gtk_text_buffer_get_start_iter(xgc->tb_mat, &start);
        gtk_text_buffer_get_end_iter(xgc->tb_mat, &end);
        textbuf = gtk_text_buffer_get_text(xgc->tb_mat, &start, &end, TRUE);
        if (!g_file_set_contents(xgc->tmpmatfilename, textbuf, -1, &err)) {
            dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                GTK_DIALOG_DESTROY_WITH_PARENT,
                GTK_MESSAGE_ERROR,
                GTK_BUTTONS_CLOSE,
                "Error writing to temporary file: %s",
                g_strdup(err->message));
            gtk_dialog_run(GTK_DIALOG(dialog));
            gtk_widget_destroy(dialog);

            return 0;
        }
    }

    return 1;
}