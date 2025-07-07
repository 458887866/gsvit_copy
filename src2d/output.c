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


/*  output.c : 
 *  all the output functions
 */
#include <stdlib.h>
#include "source.h"
#include "pool.h"
#include "constants.h"
#include "settings.h"
#include <math.h>
#include <libgwymodule/gwymodule-file.h>
//#include <app/gwyapp.h>

#define MAGIC "GWYO"
#define MAGIC2 "GWYP"
#define MAGIC_SIZE (sizeof(MAGIC)-1)

void sv_pool_allocate_output(SvPool *mp, SvSet *set);

SvOutput* 
sv_output_new()
{
    SvOutput* out = (SvOutput*)g_malloc(sizeof(SvOutput));
    
    out->data = gwy_container_new();
    out->last = 0;
   
    return out;
}


int is_output(SvPool *mp, SvSet *set)
{
    gint i;
    for (i=0; i<set->so.nimgs; i++)
    {
	    if (set->sc.step_act%set->so.imgs[i].skip == 0) return 1;
    }
    for (i=0; i<set->so.nplns; i++)
    {
	    if (set->sc.step_act%set->so.plns[i].skip == 0 
                && set->sc.step_act>=set->so.plns[i].start 
                && set->sc.step_act<=set->so.plns[i].stop ) return 1;
    }
     for (i=0; i<set->so.npnts; i++)
    {
	    if (set->sc.step_act%set->so.pnts[i].skip == 0) return 1;
    }
    for (i=0; i<set->so.nsums; i++)
    {
	    if (set->sc.step_act%set->so.sums[i].skip == 0) return 1;
    }
    return 0;

}

void sv_pool_allocate_output(SvPool *mp, SvSet *set)
{
    gint i, n;

    mp->out->outsumdata = (gdouble **)g_malloc(set->so.nsums*sizeof(gdouble *));
    for (n=0; n<set->so.nsums; n++)
    {
        mp->out->outsumdata[n] = (gdouble *)g_malloc(set->sc.nsteps*sizeof(gdouble));
        for (i=0; i<set->sc.nsteps; i++)
        {
            mp->out->outsumdata[n][i] = 0;
        }
    }

}

void sv_pool_getsum(SvPool *mp, SvSet *set)
{
    if (set->so.nsums == 0) 
        return;

    if (set->sc.verbose) {
        printf("Collecting sums...   ");
        fflush(stdout);
    }
/*  

    for (n=0; n<set->so.nsums; n++)
    {
        if (set->plan.gpumode == SV_GPUMODE_NONE) {

            for (i=set->so.sums[n].i0; i<set->so.sums[n].i1; i++)
            {
                for (j=set->so.sums[n].j0; j<set->so.sums[n].j1; j++)
                {
                    for (k=set->so.sums[n].k0; k<set->so.sums[n].k1; k++)
                    {
                        if (fabs(set->so.sums[n].epsilon - mp->epsilon->data[i][j][k])<(EPSILON_0/100.0))
                        {
                            val = 0;
                            if (set->so.sums[n].component==SV_SUM_ALL || set->so.sums[n].component==SV_SUM_ABS || set->so.sums[n].component==SV_SUM_EX)
                                val += mp->ex->data[i][j][k]*mp->ex->data[i][j][k];
                            if (set->so.sums[n].component==SV_SUM_ALL || set->so.sums[n].component==SV_SUM_ABS || set->so.sums[n].component==SV_SUM_EY)
                                val += mp->ey->data[i][j][k]*mp->ey->data[i][j][k];
                            if (set->so.sums[n].component==SV_SUM_ALL || set->so.sums[n].component==SV_SUM_ABS || set->so.sums[n].component==SV_SUM_EZ)
                                val += mp->ez->data[i][j][k]*mp->ez->data[i][j][k];

                            if (set->so.sums[n].component==SV_SUM_ABS)
                                mp->out->outsumdata[n][set->sc.step_act] += mp->sigma->data[i][j][k]*val;
                            else
                                mp->out->outsumdata[n][set->sc.step_act] += val;

                        }
                    }
                }
            }
        }
    }
    */
    if (set->sc.verbose) 
        printf("done.\n");

}

void sv_pool_output(SvPool *mp, SvSet *set)
{
    gint i, j, n;
    gboolean outim = FALSE;
    gchar key[32],titlekey[32], description[100], filename[100];
    GByteArray *buffer;
    GwyDataField *df;
    gdouble val;
    FILE *fh;
    gchar *component_string[] = {"Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "All", "Cur", "Epsilon", "Sigma", "Mu", "Sigast"};
    gint xres = set->sp.xres;
    gint yres = set->sp.yres;


    if (set->sc.verbose) {
        printf("Running output...   ");
        fflush(stdout);
    }

//    a = sv_dcube_get_sum(mp->hx);
//    b = sv_dcube_get_sum(mp->hy);
//    c = sv_dcube_get_sum(mp->hz);

//   fprintf(stderr, "Set diff: (%g %g %g)\n", a, b, c); 


//    fprintf(stderr, "Set diff: (%g %g %g)\n", (gdouble)(mp->hx->xres-1)*(mp->hx->yres-1)*(mp->hx->zres-1)*1.0 - a, 
//                               (gdouble)(mp->hx->xres-1)*(mp->hx->yres-1)*(mp->hx->zres-1)*2.0 - b,
//                               (gdouble)(mp->hx->xres-1)*(mp->hx->yres-1)*(mp->hx->zres-1)*3.0 - c);
/*
    for (i=0; i<set->so.npnts; i++) {
        if (set->sc.step_act%set->so.pnts[i].skip != 0) continue;
        

        if (set->sc.step_act==0) {
             fh = fopen(set->so.pnts[i].filebase, "w");
             if (set->so.pnts[i].component == SV_COMP_ALL) fprintf(fh, "# step    E_abs    E_x   E_y   E_z\n");
        }
        else fh = fopen(set->so.pnts[i].filebase, "a");

        switch (set->so.pnts[i].component)
        {
            case SV_COMP_EX:
            val = mp->ex->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k];
            break;

            case SV_COMP_EY:
            val = mp->ey->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k];
            break;

            case SV_COMP_EZ:
            val = mp->ez->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k];
            break;

            case SV_COMP_HX:
            val = mp->hx->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k];
            break;

            case SV_COMP_HY:
            val = mp->hy->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k];
            break;

            case SV_COMP_HZ:
            val = mp->hz->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k];
            break;

            default:
            break;
        }

        if (set->so.pnts[i].component == SV_COMP_ALL) {
            fprintf(fh, "%d %g %g %g %g\n",set->sc.step_act,
                                         sqrt(mp->ex->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k]
                                                     *mp->ex->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k] +
                                              mp->ey->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k]
                                                     *mp->ey->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k] +
                                              mp->ez->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k]
                                                     *mp->ez->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k]),
                                       mp->ex->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k],
                                       mp->ey->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k],
                                       mp->ez->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k]);
        }
        else {
            fprintf(fh,"%d %g\n", set->sc.step_act, val);
        }
        fclose(fh);

    }
    for (n=0; n<set->so.nsums; n++)
    {
        if (set->sc.step_act%set->so.sums[n].skip != 0 && set->sc.step_act<(set->sc.nsteps-1)) continue;

        fh = fopen(set->so.sums[n].filename, "w");
        for (i=0; i<set->sc.nsteps; i++)
        {
            fprintf(fh, "%d %g\n", i, mp->out->outsumdata[n][i]*set->sp.dx*set->sp.dy*set->sp.dz);
        }
        fclose(fh);
    }

*/
    for (i=0; i<set->so.nimgs; i++) {
        if (set->sc.step_act%set->so.imgs[i].skip != 0) 
            continue;
        else 
            outim = TRUE;

	
        df = gwy_data_field_new_alike(mp->ex, 1);
        gwy_data_field_set_si_unit_xy(df, gwy_si_unit_new("pixels"));

        gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("V"));
       

        if (set->so.imgs[i].component == SV_COMP_EX) gwy_data_field_copy(mp->ex, df, FALSE);
        else if (set->so.imgs[i].component == SV_COMP_EY) gwy_data_field_copy(mp->ey, df, FALSE);
        else if (set->so.imgs[i].component == SV_COMP_EZ) gwy_data_field_copy(mp->ez, df, FALSE);
        else if (set->so.imgs[i].component == SV_COMP_HX) gwy_data_field_copy(mp->hx, df, FALSE);
        else if (set->so.imgs[i].component == SV_COMP_HY) gwy_data_field_copy(mp->hy, df, FALSE);
        else if (set->so.imgs[i].component == SV_COMP_HZ) gwy_data_field_copy(mp->hz, df, FALSE);
        else if (set->so.imgs[i].component == SV_COMP_EPSILON || set->so.imgs[i].component == SV_COMP_SIGMA
               || set->so.imgs[i].component == SV_COMP_MU || set->so.imgs[i].component == SV_COMP_SIGAST) {
            //TODO, here interpret also tabulated and Drude materials
            if (set->so.imgs[i].component == SV_COMP_EPSILON && mp->epsilon) {
               gwy_data_field_copy(mp->epsilon, df, FALSE);
               gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("F/m"));
            }
            if (set->so.imgs[i].component == SV_COMP_SIGMA && mp->sigma) {
               gwy_data_field_copy(mp->sigma, df, FALSE);
               gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("S/m"));
            }
            if (set->so.imgs[i].component == SV_COMP_MU && mp->mu) {
               gwy_data_field_copy(mp->mu, df, FALSE);
               gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("H/m"));
            }
            if (set->so.imgs[i].component == SV_COMP_SIGAST && mp->sigast) {
               gwy_data_field_copy(mp->sigast, df, FALSE);
               gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("H"));
            }
        }

        g_snprintf(key, sizeof(key), "/%u/data", mp->out->last);
        g_snprintf(titlekey, sizeof(titlekey), "/%u/data/title", mp->out->last++);
        g_snprintf(description, sizeof(description), "%s step %d, %s", set->so.imgs[i].filebase, set->sc.step_act, 
                   component_string[set->so.imgs[i].component]);

        gwy_container_set_object_by_name(mp->out->data, g_strdup(key), df);
        gwy_container_set_string_by_name(mp->out->data, g_strdup(titlekey), (guchar*)g_strdup(description));

        g_object_unref(df);
    }

    if (outim) {
        /*save container*/
        fh = fopen(set->so.outfile, "wb");
        if (!fh) 
            fprintf(stderr, "Error opening file\n");
        buffer = gwy_serializable_serialize(G_OBJECT(mp->out->data), NULL);

        if (fwrite(MAGIC2, 1, MAGIC_SIZE, fh) != MAGIC_SIZE || fwrite(buffer->data, 1, buffer->len, fh) != buffer->len) {
            fprintf(stderr, "Error saving file\n");
        }

        fclose(fh);         
        g_byte_array_free(buffer, TRUE); 
    }

    for (n=0; n<set->so.nplns; n++) {
        if (set->sc.step_act%set->so.plns[n].skip != 0) 
            continue;
        if (set->sc.step_act<set->so.plns[n].start) 
            continue;
        if (set->sc.step_act>=set->so.plns[n].stop) 
            continue;
         
        df = gwy_data_field_new_alike(mp->ex, 1);

        if (set->so.plns[n].component == SV_COMP_EX) gwy_data_field_copy(mp->ex, df, FALSE);
        else if (set->so.plns[n].component == SV_COMP_EY) gwy_data_field_copy(mp->ey, df, FALSE);
        else if (set->so.plns[n].component == SV_COMP_EZ) gwy_data_field_copy(mp->ez, df, FALSE);
        else if (set->so.plns[n].component == SV_COMP_HX) gwy_data_field_copy(mp->hx, df, FALSE);
        else if (set->so.plns[n].component == SV_COMP_HY) gwy_data_field_copy(mp->hy, df, FALSE);
        else if (set->so.plns[n].component == SV_COMP_HZ) gwy_data_field_copy(mp->hz, df, FALSE);
        else if (set->so.plns[n].component == SV_COMP_EPSILON || set->so.plns[n].component == SV_COMP_SIGMA
               || set->so.plns[n].component == SV_COMP_MU || set->so.plns[n].component == SV_COMP_SIGAST) {
            //TODO, here interpret also tabulated and Drude materials
            if (set->so.plns[n].component == SV_COMP_EPSILON && mp->epsilon) {
               gwy_data_field_copy(mp->epsilon, df, FALSE);
            }
            if (set->so.plns[n].component == SV_COMP_SIGMA && mp->sigma) {
               gwy_data_field_copy(mp->sigma, df, FALSE);
            }
            if (set->so.plns[n].component == SV_COMP_MU && mp->mu) {
               gwy_data_field_copy(mp->mu, df, FALSE);
            }
            if (set->so.plns[n].component == SV_COMP_SIGAST && mp->sigast) {
               gwy_data_field_copy(mp->sigast, df, FALSE);
            }
        }

        if (set->so.plns[n].logscale) { //this means ascii mode here
            g_snprintf(filename, 100, "%s_%.4d", set->so.plns[n].filebase, set->sc.step_act);

            fh = fopen(filename, "w");
            if (fh == NULL) {
                fprintf(stderr, "Error: cannot open file %s for writing\n", filename); 
                continue; 
            }

            for (i = 0; i < xres; i++) {
                for (j = 0; j < yres; j++) {
                    fprintf(fh, "%g ", gwy_data_field_get_val(df, i, j));
                }
                fprintf(fh, "\n");
            }

            fclose(fh);

        } else {
            g_snprintf(filename, 100, "%s_%.4d.raw", set->so.plns[n].filebase, set->sc.step_act);

            fh = fopen(filename, "wb");
            if (fh == NULL) {
                fprintf(stderr, "Error: cannot open file %s for writing\n", filename); 
                continue; 
            }

            for (i = 0; i < xres; i++) {
                for (j = 0; j < yres; j++) {
                    val = gwy_data_field_get_val(df, i, j);
                    fwrite(&val, sizeof(gdouble), 1, fh);
               }
            }

            fclose(fh);
        }



    }


    if (set->sc.verbose) 
        printf("done.\n");
}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

