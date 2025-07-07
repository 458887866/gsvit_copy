
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


/*  farfield.c : 
 *   *  near-to-far-field calculation algorithm based on Kirchoff surface integral
 *    */


#include "farfield.h"
#include "pool.h"
#include "settings.h"
#include "constants.h"
#include <math.h>

SvFarfield* 
sv_farfield_new()
{
    SvFarfield* farfield = (SvFarfield*)g_malloc(sizeof(SvFarfield));
  
    farfield->nfs = (SvNFFFStore*)g_malloc(sizeof(SvNFFFStore)); 
    return farfield;
}


void sv_pool_farfield(SvPool *mp, SvSet *set)
{
    gint i, j;
    gint jposlr, jpostb;
    gint nsteps;
    gdouble dt, dx, dy;
    FILE *fw = NULL;

    if (set->sf.savefile==NULL) return;
    
    if (set->sc.verbose) {
        printf("Running NFFF...   ");
        fflush(stdout);
    }

    if (set->sc.step_act==0) {

        if (!set->sf.save) {

            fw = fopen(set->sf.savefile, "rb");    
            fread(&(set->sf.i0), sizeof(gint), 1, fw);
            fread(&(set->sf.j0), sizeof(gint), 1, fw);
            fread(&(set->sf.i1), sizeof(gint), 1, fw);
            fread(&(set->sf.j1), sizeof(gint), 1, fw);
            fread(&nsteps, sizeof(gint), 1, fw);
            fread(&dt, sizeof(gdouble), 1, fw);
            fread(&dx, sizeof(gdouble), 1, fw);
            fread(&dy, sizeof(gdouble), 1, fw);
            if (nsteps != set->sc.nsteps) fprintf(stderr, "Error: different number of steps in farfield data and pool settings\n");
            if (dt != set->plan.dt) fprintf(stderr, "Error: different time steps in loaded data and pool settings\n");
            if (set->sc.verbose) {printf("\nLoading %dx%d and %dx%d fields for NFFF applications... ", set->sf.j1 - set->sf.j0, set->sc.nsteps,
                                         set->sf.j1 - set->sf.j0, set->sc.nsteps); fflush(stdout);}

        }

        if (set->tmmode) {
            if (set->sc.verbose) {printf("\nAllocating %dx%d and %dx%d fields for NFFF storage... ", set->sf.j1 - set->sf.j0, set->sc.nsteps,
                                         set->sf.j1 - set->sf.j0, set->sc.nsteps); fflush(stdout);}
            mp->farfield->nfs->tm_left_hyp = gwy_data_field_new(set->sf.j1 - set->sf.j0 + 1, set->sc.nsteps, set->sf.j1 - set->sf.j0 + 1, set->sc.nsteps, 1);
            mp->farfield->nfs->tm_left_hym = gwy_data_field_new_alike(mp->farfield->nfs->tm_left_hyp, 0);
            mp->farfield->nfs->tm_left_ez = gwy_data_field_new_alike(mp->farfield->nfs->tm_left_hyp, 0);

            mp->farfield->nfs->tm_right_hyp = gwy_data_field_new_alike(mp->farfield->nfs->tm_left_hyp, 0);
            mp->farfield->nfs->tm_right_hym = gwy_data_field_new_alike(mp->farfield->nfs->tm_left_hyp, 0);
            mp->farfield->nfs->tm_right_ez = gwy_data_field_new_alike(mp->farfield->nfs->tm_left_hyp, 0);

            mp->farfield->nfs->tm_top_hxp = gwy_data_field_new(set->sf.i1 - set->sf.i0 + 1, set->sc.nsteps, set->sf.i1 - set->sf.i0 + 1, set->sc.nsteps, 1);
            mp->farfield->nfs->tm_top_hxm = gwy_data_field_new_alike(mp->farfield->nfs->tm_top_hxp, 0);
            mp->farfield->nfs->tm_top_ez = gwy_data_field_new_alike(mp->farfield->nfs->tm_top_hxp, 0);

            mp->farfield->nfs->tm_bottom_hxp = gwy_data_field_new_alike(mp->farfield->nfs->tm_top_hxp, 0);
            mp->farfield->nfs->tm_bottom_hxm = gwy_data_field_new_alike(mp->farfield->nfs->tm_top_hxp, 0);
            mp->farfield->nfs->tm_bottom_ez = gwy_data_field_new_alike(mp->farfield->nfs->tm_top_hxp, 0);

            if (!set->sf.save) {
                fread(mp->farfield->nfs->tm_left_hyp->data, sizeof(gdouble), (set->sf.j1 - set->sf.j0 + 1)*set->sc.nsteps, fw);
                fread(mp->farfield->nfs->tm_left_hym->data, sizeof(gdouble), (set->sf.j1 - set->sf.j0 + 1)*set->sc.nsteps, fw);
                fread(mp->farfield->nfs->tm_left_ez->data, sizeof(gdouble), (set->sf.j1 - set->sf.j0 + 1)*set->sc.nsteps, fw);
                fread(mp->farfield->nfs->tm_right_hyp->data, sizeof(gdouble), (set->sf.j1 - set->sf.j0 + 1)*set->sc.nsteps, fw);
                fread(mp->farfield->nfs->tm_right_hym->data, sizeof(gdouble), (set->sf.j1 - set->sf.j0 + 1)*set->sc.nsteps, fw);
                fread(mp->farfield->nfs->tm_right_ez->data, sizeof(gdouble), (set->sf.j1 - set->sf.j0 + 1)*set->sc.nsteps, fw);

                fread(mp->farfield->nfs->tm_top_hxp->data, sizeof(gdouble), (set->sf.i1 - set->sf.i0 + 1)*set->sc.nsteps, fw);
                fread(mp->farfield->nfs->tm_top_hxm->data, sizeof(gdouble), (set->sf.i1 - set->sf.i0 + 1)*set->sc.nsteps, fw);
                fread(mp->farfield->nfs->tm_top_ez->data, sizeof(gdouble), (set->sf.i1 - set->sf.i0 + 1)*set->sc.nsteps, fw);
                fread(mp->farfield->nfs->tm_bottom_hxp->data, sizeof(gdouble), (set->sf.i1 - set->sf.i0 + 1)*set->sc.nsteps, fw);
                fread(mp->farfield->nfs->tm_bottom_hxm->data, sizeof(gdouble), (set->sf.i1 - set->sf.i0 + 1)*set->sc.nsteps, fw);
                fread(mp->farfield->nfs->tm_bottom_ez->data, sizeof(gdouble), (set->sf.i1 - set->sf.i0 + 1)*set->sc.nsteps, fw);

                fclose(fw);

            }
        }
    }
    
    jposlr = set->sc.step_act*(set->sf.j1 - set->sf.j0 + 1);
    jpostb = set->sc.step_act*(set->sf.i1 - set->sf.i0 + 1);

    if (set->sf.save) {
        for (j=set->sf.j0; j<=set->sf.j1; j++) {
            if (set->sf.ramahi) {
                mp->farfield->nfs->tm_left_hyp->data[(j - set->sf.j0) + jposlr] = mp->ez->data[set->sf.i0+1 + j*set->sp.xres];
                mp->farfield->nfs->tm_left_hym->data[(j - set->sf.j0) + jposlr] = mp->ez->data[set->sf.i0-1 + j*set->sp.xres];
                mp->farfield->nfs->tm_left_ez->data[(j - set->sf.j0) + jposlr] = mp->ez->data[set->sf.i0 + j*set->sp.xres];

                mp->farfield->nfs->tm_right_hyp->data[(j - set->sf.j0) + jposlr] = mp->ez->data[set->sf.i1+1 + j*set->sp.xres];
                mp->farfield->nfs->tm_right_hym->data[(j - set->sf.j0) + jposlr] = mp->ez->data[set->sf.i1-1 + j*set->sp.xres];
                mp->farfield->nfs->tm_right_ez->data[(j - set->sf.j0) + jposlr] = mp->ez->data[set->sf.i1 + j*set->sp.xres];
            }
            else {
                mp->farfield->nfs->tm_left_hyp->data[(j - set->sf.j0) + jposlr] = mp->hy->data[set->sf.i0 + j*set->sp.xres];
                mp->farfield->nfs->tm_left_hym->data[(j - set->sf.j0) + jposlr] = mp->hy->data[set->sf.i0-1 + j*set->sp.xres];
                mp->farfield->nfs->tm_left_ez->data[(j - set->sf.j0) + jposlr] = mp->ez->data[set->sf.i0 + j*set->sp.xres];

                mp->farfield->nfs->tm_right_hyp->data[(j - set->sf.j0) + jposlr] = mp->hy->data[set->sf.i1 + j*set->sp.xres];
                mp->farfield->nfs->tm_right_hym->data[(j - set->sf.j0) + jposlr] = mp->hy->data[set->sf.i1-1 + j*set->sp.xres];
                mp->farfield->nfs->tm_right_ez->data[(j - set->sf.j0) + jposlr] = mp->ez->data[set->sf.i1 + j*set->sp.xres];
            }
        }

        for (i=set->sf.i0; i<=set->sf.i1; i++) {
            if (set->sf.ramahi) {
                mp->farfield->nfs->tm_top_hxp->data[(i - set->sf.i0) + jpostb] = mp->ez->data[i + (set->sf.j0+1)*set->sp.xres];
                mp->farfield->nfs->tm_top_hxm->data[(i - set->sf.i0) + jpostb] = mp->ez->data[i + (set->sf.j0-1)*set->sp.xres];
                mp->farfield->nfs->tm_top_ez->data[(i - set->sf.i0) + jpostb] = mp->ez->data[i + set->sf.j0*set->sp.xres];

                mp->farfield->nfs->tm_bottom_hxp->data[(i - set->sf.i0) + jpostb] = mp->ez->data[i + (set->sf.j1+1)*set->sp.xres];
                mp->farfield->nfs->tm_bottom_hxm->data[(i - set->sf.i0) + jpostb] = mp->ez->data[i + (set->sf.j1-1)*set->sp.xres];
                mp->farfield->nfs->tm_bottom_ez->data[(i - set->sf.i0) + jpostb] = mp->ez->data[i + set->sf.j1*set->sp.xres];
            }
            else {
                mp->farfield->nfs->tm_top_hxp->data[(i - set->sf.i0) + jpostb] = mp->hx->data[i + set->sf.j0*set->sp.xres];
                mp->farfield->nfs->tm_top_hxm->data[(i - set->sf.i0) + jpostb] = mp->hx->data[i + (set->sf.j0-1)*set->sp.xres];
                mp->farfield->nfs->tm_top_ez->data[(i - set->sf.i0) + jpostb] = mp->ez->data[i + set->sf.j0*set->sp.xres];

                mp->farfield->nfs->tm_bottom_hxp->data[(i - set->sf.i0) + jpostb] = mp->hx->data[i + set->sf.j1*set->sp.xres];
                mp->farfield->nfs->tm_bottom_hxm->data[(i - set->sf.i0) + jpostb] = mp->hx->data[i + (set->sf.j1-1)*set->sp.xres];
                mp->farfield->nfs->tm_bottom_ez->data[(i - set->sf.i0) + jpostb] = mp->ez->data[i + set->sf.j1*set->sp.xres];

            }
        }
    } else {
        for (j=set->sf.j0; j<=set->sf.j1; j++) {
            mp->hy->data[set->sf.i0 + j*set->sp.xres] = mp->farfield->nfs->tm_left_hyp->data[(j - set->sf.j0) + jposlr];
            mp->hy->data[set->sf.i0-1 + j*set->sp.xres] = mp->farfield->nfs->tm_left_hym->data[(j - set->sf.j0) + jposlr];
            mp->ez->data[set->sf.i0 + j*set->sp.xres] = mp->farfield->nfs->tm_left_ez->data[(j - set->sf.j0) + jposlr];

            mp->hy->data[set->sf.i1 + j*set->sp.xres] = mp->farfield->nfs->tm_right_hyp->data[(j - set->sf.j0) + jposlr];
            mp->hy->data[set->sf.i1-1 + j*set->sp.xres] = mp->farfield->nfs->tm_right_hym->data[(j - set->sf.j0) + jposlr];
            mp->ez->data[set->sf.i1 + j*set->sp.xres] = mp->farfield->nfs->tm_right_ez->data[(j - set->sf.j0) + jposlr];
        }

        for (i=set->sf.i0; i<=set->sf.i1; i++) {
            mp->hx->data[i + set->sf.j0*set->sp.xres] = mp->farfield->nfs->tm_top_hxp->data[(i - set->sf.i0) + jpostb];
            mp->hx->data[i + (set->sf.j0-1)*set->sp.xres] = mp->farfield->nfs->tm_top_hxm->data[(i - set->sf.i0) + jpostb];
            mp->ez->data[i + set->sf.j0*set->sp.xres] = mp->farfield->nfs->tm_top_ez->data[(i - set->sf.i0) + jpostb];

            mp->hx->data[i + set->sf.j1*set->sp.xres] = mp->farfield->nfs->tm_bottom_hxp->data[(i - set->sf.i0) + jpostb];
            mp->hx->data[i + (set->sf.j1-1)*set->sp.xres] = mp->farfield->nfs->tm_bottom_hxm->data[(i - set->sf.i0) + jpostb];
            mp->ez->data[i + set->sf.j1*set->sp.xres] = mp->farfield->nfs->tm_bottom_ez->data[(i - set->sf.i0) + jpostb];
        }
 
    }

    if (set->sf.save && set->sc.step_act==(set->sc.nsteps-1)) {
        if (set->sc.verbose) {printf("\nStoring NFF data to %s... ", set->sf.savefile); fflush(stdout);}

        fw = fopen(set->sf.savefile, "wb");    
        fwrite(&(set->sf.i0), sizeof(gint), 1, fw);
        fwrite(&(set->sf.j0), sizeof(gint), 1, fw);
        fwrite(&(set->sf.i1), sizeof(gint), 1, fw);
        fwrite(&(set->sf.j1), sizeof(gint), 1, fw);
        fwrite(&(set->sc.nsteps), sizeof(gint), 1, fw);
        fwrite(&(set->plan.dt), sizeof(gdouble), 1, fw);
        fwrite(&(set->sp.dx), sizeof(gdouble), 1, fw);
        fwrite(&(set->sp.dy), sizeof(gdouble), 1, fw);



        fwrite(mp->farfield->nfs->tm_left_hyp->data, sizeof(gdouble), (set->sf.j1 - set->sf.j0 + 1)*set->sc.nsteps, fw);
        fwrite(mp->farfield->nfs->tm_left_hym->data, sizeof(gdouble), (set->sf.j1 - set->sf.j0 + 1)*set->sc.nsteps, fw);
        fwrite(mp->farfield->nfs->tm_left_ez->data, sizeof(gdouble), (set->sf.j1 - set->sf.j0 + 1)*set->sc.nsteps, fw);
        fwrite(mp->farfield->nfs->tm_right_hyp->data, sizeof(gdouble), (set->sf.j1 - set->sf.j0 + 1)*set->sc.nsteps, fw);
        fwrite(mp->farfield->nfs->tm_right_hym->data, sizeof(gdouble), (set->sf.j1 - set->sf.j0 + 1)*set->sc.nsteps, fw);
        fwrite(mp->farfield->nfs->tm_right_ez->data, sizeof(gdouble), (set->sf.j1 - set->sf.j0 + 1)*set->sc.nsteps, fw);

        fwrite(mp->farfield->nfs->tm_top_hxp->data, sizeof(gdouble), (set->sf.i1 - set->sf.i0 + 1)*set->sc.nsteps, fw);
        fwrite(mp->farfield->nfs->tm_top_hxm->data, sizeof(gdouble), (set->sf.i1 - set->sf.i0 + 1)*set->sc.nsteps, fw);
        fwrite(mp->farfield->nfs->tm_top_ez->data, sizeof(gdouble), (set->sf.i1 - set->sf.i0 + 1)*set->sc.nsteps, fw);
        fwrite(mp->farfield->nfs->tm_bottom_hxp->data, sizeof(gdouble), (set->sf.i1 - set->sf.i0 + 1)*set->sc.nsteps, fw);
        fwrite(mp->farfield->nfs->tm_bottom_hxm->data, sizeof(gdouble), (set->sf.i1 - set->sf.i0 + 1)*set->sc.nsteps, fw);
        fwrite(mp->farfield->nfs->tm_bottom_ez->data, sizeof(gdouble), (set->sf.i1 - set->sf.i0 + 1)*set->sc.nsteps, fw);

        fclose(fw);
    }
 
    if (set->sc.verbose) printf("done.\n");

}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

