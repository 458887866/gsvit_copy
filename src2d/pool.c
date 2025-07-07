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


 /*  pool.c :
  *  main algorithms (Yee solver, including PRLC, PEC, parts of TSF and CPML)
  */

#if defined (_WIN32)
#include "../config_win.h"
#else
#include "../config.h"
#endif
#include <glib.h>
#include <stdlib.h>
#include "pool.h"
#include "constants.h"
#include "source.h"
#include "boundary.h"
#include <math.h>

#ifdef UCUDA
SvGpu* sv_gpu_new(SvPool *mp, SvSet *set);
#endif

SvPool*
sv_pool_new(SvSet *set)
{
    SvPool *mp = (SvPool *)g_malloc(sizeof(SvPool));

    mp->ex = gwy_data_field_new(set->sp.xres, set->sp.yres,
                                set->sp.xres*set->sp.dx,
                                set->sp.yres*set->sp.dy,
                                1);

    mp->ey = gwy_data_field_new_alike(mp->ex, 1);
    mp->ez = gwy_data_field_new_alike(mp->ex, 1);
    mp->hx = gwy_data_field_new_alike(mp->ex, 1);
    mp->hy = gwy_data_field_new_alike(mp->ex, 1);
    mp->hz = gwy_data_field_new_alike(mp->ex, 1);
    mp->epsilon = NULL;
    mp->mu = NULL;
    mp->sigma = NULL;
    mp->sigast = NULL;

    //je-li treba, alokovat integerovou materialovou krychli

    /*if full material is requested, otherwise use matmode from settings*/
    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC) {

        mp->epsilon = gwy_data_field_new(set->sp.xres, set->sp.yres,
                                         set->sp.xres*set->sp.dx,
                                         set->sp.yres*set->sp.dy,
                                         1);
        mp->sigma = gwy_data_field_new_alike(mp->epsilon, 1);
        gwy_data_field_fill(mp->epsilon, EPSILON_0);

    }
    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_MAGNETIC) {
        mp->mu =
            gwy_data_field_new(set->sp.xres, set->sp.yres,
                               set->sp.xres*set->sp.dx,
                               set->sp.yres*set->sp.dy,
                               1);
        mp->sigast = gwy_data_field_new_alike(mp->mu, 1);
        gwy_data_field_fill(mp->mu, MU_0);
    }

    /*alloc other structures*/
    mp->src = sv_source_new();
    mp->out = sv_output_new();
    mp->bnd = sv_boundary_new(set);
    mp->farfield = sv_farfield_new();
    mp->nmat = 0;
    mp->mats = NULL;

    if (set->sc.verbose) printf("Pool allocated\n");

    return mp;
}

void
sv_pool_allocate_gpus(SvPool *mp, SvSet *set)
{
#ifdef UCUDA
    mp->gpu = sv_gpu_new(mp, set);
#endif

}

int
sv_pool_ystep_h(SvPool *mp, SvSet *set)
{
    gint xnres = set->sp.xres - 1;
    gint ynres = set->sp.yres - 1;
    gint xres = set->sp.xres;
    gint xnzero = 0;
    gint ynzero = 0;
    gint i = 0, j = 0;
    gdouble da = 0, db = 0;
    gdouble dx = set->sp.dx;
    gdouble dt = set->plan.dt;

    if (set->sc.verbose) {
        printf("Running Yee H step...   ");
        fflush(stdout);
        //      printf("\n");
        //    fprintf(stderr, "hstep\n");
    }

#pragma omp parallel default(shared) private(i, j, da, db)
#pragma omp for nowait
    for (i = xnzero; i < xnres; i++) {
        for (j = ynzero; j < ynres; j++) {
            if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_MAGNETIC) {
                da = (1 - mp->sigast->data[i + j*xres] * dt / 2 / mp->mu->data[i + j*xres]) / (1 + mp->sigast->data[i + j*xres] * dt / 2 / mp->mu->data[i + j*xres]);
                db = (dt / mp->mu->data[i + j*xres] / dx) / (1 + mp->sigast->data[i + j*xres] * dt / 2 / mp->mu->data[i + j*xres]);
            } else {
                da = 1;
                db = dt / MU_0 / dx;
            }
            if (set->tmmode) //tm mode
            {
                mp->hx->data[i + j*xres] = da*mp->hx->data[i + j*xres] + db*(-(mp->ez->data[i + (j + 1)*xres] - mp->ez->data[i + j*xres]));

                mp->hy->data[i + j*xres] = da*mp->hy->data[i + j*xres] + db*((mp->ez->data[(i + 1) + j*xres] - mp->ez->data[i + j*xres]));

            } else {           //te mode                 
                mp->hz->data[i + j*xres] = da*mp->hz->data[i + j*xres] + db*((mp->ex->data[i + (j + 1)*xres] - mp->ex->data[i + j*xres]) -
                    (mp->ey->data[(i + 1) + j*xres] - mp->ey->data[i + j*xres]));

                /*
                            if ((set->sc.step_act)==61) {


                                if (i==60 && j==60) fprintf(stderr, "%d %d da: %g phz: %g  db: %g diff %g ([%d %d] %g [%d %d] %g [%d %d] %g [%d %d] %g)   = %g\n", i, j, da, phz, db,
                                                                              ((mp->ex->data[i+(j+1)*xres] - mp->ex->data[i+j*xres]) - (mp->ey->data[(i+1)+j*xres] - mp->ey->data[i+j*xres])),
                                                                              i, (j+1), mp->ex->data[i+(j+1)*xres],
                                                                              i, j, mp->ex->data[i+(j)*xres],
                                                                              i+1, j, mp->ey->data[(i+1)+j*xres],
                                                                              i, j, mp->ey->data[i+j*xres],
                                                                              mp->hz->data[i+j*xres]);

                                if (i==60 && j==99) fprintf(stderr, "%d %d da: %g phz: %g  db: %g diff %g ([%d %d] %g [%d %d] %g   [%d %d] %g  [%d %d] %g)   = %g\n", i, j, da, phz, db, ((mp->ex->data[i+(j+1)*xres] - mp->ex->data[i+j*xres]) - (mp->ey->data[(i+1)+j*xres] - mp->ey->data[i+j*xres])), i, j+1, mp->ex->data[i+(j+1)*xres], i, j, mp->ex->data[i+(j)*xres], i+1, j, mp->ey->data[(i+1)+j*xres], i, j, mp->ey->data[i+j*xres], mp->hz->data[i+j*xres]);
                                if (i==99 && j==60) fprintf(stderr, "%d %d da: %g phz: %g  db: %g diff %g ([%d %d] %g [%d %d] %g   [%d %d] %g  [%d %d] %g)   = %g\n", i, j, da, phz, db, ((mp->ex->data[i+(j+1)*xres] - mp->ex->data[i+j*xres]) - (mp->ey->data[(i+1)+j*xres] - mp->ey->data[i+j*xres])), i, j+1, mp->ex->data[i+(j+1)*xres], i, j, mp->ex->data[i+(j)*xres], i+1, j, mp->ey->data[(i+1)+j*xres], i, j, mp->ey->data[i+j*xres], mp->hz->data[i+j*xres]);
                                if (i==99 && j==99) fprintf(stderr, "%d %d da: %g phz: %g  db: %g diff %g ([%d %d] %g [%d %d] %g   [%d %d] %g  [%d %d] %g)   = %g\n", i, j, da, phz, db, ((mp->ex->data[i+(j+1)*xres] - mp->ex->data[i+j*xres]) - (mp->ey->data[(i+1)+j*xres] - mp->ey->data[i+j*xres])), i, j+1, mp->ex->data[i+(j+1)*xres], i, j, mp->ex->data[i+(j)*xres], i+1, j, mp->ey->data[(i+1)+j*xres], i, j, mp->ey->data[i+j*xres], mp->hz->data[i+j*xres]);


                }


                                if (i>=57 && i<63 && j>=57 && j<63) {
                                   printf("%+5.3e  ", mp->hz->data[i+j*xres]);
                                }

                                if (i>=57 && i<63 && j==64) {
                                   printf(" ||   ");
                                }


                                if (i>=57 && i<63 && j>=97 && j<103) {
                                   printf("%+5.3e  ", mp->hz->data[i+j*xres]);
                                }


                                if (i>=97 && i<103 && j>=57 && j<63) {
                                   printf("%+5.3e  ", mp->hz->data[i+j*xres]);
                                }

                                if (i>=97 && i<103 && j==64) {
                                   printf(" ||   ");
                                }


                                if (i>=97 && i<103 && j>=97 && j<103) {
                                   printf("%+5.3e  ", mp->hz->data[i+j*xres]);
                                }
                */
            }

        }
        /*
                  if (i>=57 && i<63) {
                     printf("\n");
                }
                if (i==64) printf("\n");


                if (i>=97 && i<103) {
                     printf("\n");
                }
        */
    }

    if (set->sc.verbose) printf("done.\n");
    return 0;

}


int
sv_pool_ystep_e(SvPool *mp, SvSet *set)
{
    gint xres = set->sp.xres;
    gint yres = set->sp.yres;
    gint matindex = 0;
    gint xnzero = 1;
    gint ynzero = 1;
    gint i = 0, j = 0;
    gdouble ca = 0, cb = 0, cax = 0, cbx = 0, cay = 0, cby = 0;
    gdouble dx = set->sp.dx;
    gdouble dt = set->plan.dt;
    gdouble dtepsdx = dt / EPSILON_0 / dx;

    if (set->sc.verbose) {
        printf("Running Yee E step...   ");
        fflush(stdout);
        // printf("\n");
        // fprintf(stderr, "E step\n");
    }

#pragma omp parallel default(shared) private(i, j, matindex, ca, cb, cax, cbx, cay, cby)    
#pragma omp for nowait    
    for (i = xnzero; i < xres; i++) {
        for (j = ynzero; j < yres; j++) {            

            ca = cb = cax = cbx = cay = cby = 0;

            if (mp->nmat > 0)
                matindex = (gint)mp->mat->data[i + j*xres];
            else
                matindex = -1;

            if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC) {
                if (mp->nmat == 0 || matindex == 0) { //linear local material, no other information
                    ca = (1 - mp->sigma->data[i + j*xres] * dt / 2 / mp->epsilon->data[i + j*xres])
                        / (1 + mp->sigma->data[i + j*xres] * dt / 2 / mp->epsilon->data[i + j*xres]);
                    cb = (dt / mp->epsilon->data[i + j*xres] / dx) / (1 + mp->sigma->data[i + j*xres] * dt / 2 / mp->epsilon->data[i + j*xres]);

                    cax = (1 - 0.5*(mp->sigma->data[i + j*xres] + mp->sigma->data[(i - 1) + j*xres])*dt / 2 / (0.5*(mp->epsilon->data[i + j*xres] + mp->epsilon->data[(i - 1) + j*xres])))
                        / (1 + 0.5*(mp->sigma->data[i + j*xres] + mp->sigma->data[(i - 1) + j*xres])*dt / 2 / (0.5*(mp->epsilon->data[i + j*xres] + mp->epsilon->data[(i - 1) + j*xres])));
                    cbx = (dt / (0.5*(mp->epsilon->data[i + j*xres] + mp->epsilon->data[(i - 1) + j*xres])) / dx) / (1 + 0.5*(mp->sigma->data[i + j*xres] + mp->sigma->data[(i - 1) + j*xres])*dt / 2 / (0.5*(mp->epsilon->data[i + j*xres] + mp->epsilon->data[(i - 1) + j*xres])));

                    cay = (1 - 0.5*(mp->sigma->data[i + j*xres] + mp->sigma->data[i + (j - 1)*xres])*dt / 2 / (0.5*(mp->epsilon->data[i + j*xres] + mp->epsilon->data[i + (j - 1)*xres])))
                        / (1 + 0.5*(mp->sigma->data[i + j*xres] + mp->sigma->data[i + (j - 1)*xres])*dt / 2 / (0.5*(mp->epsilon->data[i + j*xres] + mp->epsilon->data[i + (j - 1)*xres])));
                    cby = (dt / (0.5*(mp->epsilon->data[i + j*xres] + mp->epsilon->data[i + (j - 1)*xres])) / dx) / (1 + 0.5*(mp->sigma->data[i + j*xres] + mp->sigma->data[i + (j - 1)*xres])*dt / 2 / (0.5*(mp->epsilon->data[i + j*xres] + mp->epsilon->data[i + (j - 1)*xres])));


                } else if (mp->mats[matindex].type == SV_MAT_LINTAB) { // some tabulated data inside local material
                    ca = (1 - mp->mats[matindex].sigma*dt / 2 / mp->mats[matindex].epsilon)
                        / (1 + mp->mats[matindex].sigma*dt / 2 / mp->mats[matindex].epsilon);
                    cb = (dt / mp->mats[matindex].epsilon / dx)
                        / (1 + mp->mats[matindex].sigma*dt / 2 / mp->mats[matindex].epsilon);
                }
            } else if (mp->nmat > 0 && mp->mats[matindex].type == SV_MAT_LINTAB) {
                ca = (1 - mp->mats[matindex].sigma*dt / 2 / mp->mats[matindex].epsilon)
                    / (1 + mp->mats[matindex].sigma*dt / 2 / mp->mats[matindex].epsilon);
                cb = (dt / mp->mats[matindex].epsilon / dx)
                    / (1 + mp->mats[matindex].sigma*dt / 2 / mp->mats[matindex].epsilon);
            } else { //vacuum
                ca = 1;
                cb = dtepsdx;
            }

            if (!set->tmmode) { //te mode
                mp->ex->data[i + j*xres] = cay*mp->ex->data[i + j*xres] + cby*((mp->hz->data[i + j*xres] - mp->hz->data[i + (j - 1)*xres]));

                mp->ey->data[i + j*xres] = cax*mp->ey->data[i + j*xres] + cbx*(-(mp->hz->data[i + j*xres] - mp->hz->data[(i - 1) + j*xres]));

                /*
                                if (i>=57 && i<63 && j>=57 && j<63) {
                                   printf("(%d %d %5.3g %+5.3e) ", i, j, cb, mp->ex->data[i+j*xres]);
                                }

                                if (i>=57 && i<63 && j==64) {
                                   printf(" ||   ");
                                }


                                if (i>=57 && i<63 && j>=97 && j<103) {
                                   printf("(%d %d %5.3g %+5.3e) ", i, j, cb, mp->ex->data[i+j*xres]);
                                }


                                if ((set->sc.step_act)==61) {

                                if (i==60 && j==60) fprintf(stderr, "%d %d ca: %g pex: %g pey: %g cb: %g  diff: %g ([%d %d] %g [%d %d] %g [%d %d] %g) = %g, %g\n", i, j, ca, pex, pey, cb, -(mp->hz->data[i+j*xres]-mp->hz->data[(i-1)+j*xres]), i, j, mp->hz->data[i+j*xres], i-1, j, mp->hz->data[(i-1)+j*xres], i, j-1, mp->hz->data[i+(j-1)*xres], mp->ex->data[i+j*xres], mp->ey->data[i+j*xres]);
                                if (i==60 && j==99) fprintf(stderr, "%d %d ca: %g pex: %g pey: %g cb: %g  diff: %g ([%d %d] %g [%d %d] %g [%d %d] %g) = %g, %g\n", i, j, ca, pex, pey, cb, -(mp->hz->data[i+j*xres]-mp->hz->data[(i-1)+j*xres]), i, j, mp->hz->data[i+j*xres], i-1, j, mp->hz->data[(i-1)+j*xres], i, j-1, mp->hz->data[i+(j-1)*xres], mp->ex->data[i+j*xres], mp->ey->data[i+j*xres]);
                                if (i==99 && j==60) fprintf(stderr, "%d %d ca: %g pex: %g pey: %g cb: %g  diff: %g ([%d %d] %g [%d %d] %g [%d %d] %g) = %g, %g\n", i, j, ca, pex, pey, cb, -(mp->hz->data[i+j*xres]-mp->hz->data[(i-1)+j*xres]), i, j, mp->hz->data[i+j*xres], i-1, j, mp->hz->data[(i-1)+j*xres], i, j-1, mp->hz->data[i+(j-1)*xres], mp->ex->data[i+j*xres], mp->ey->data[i+j*xres]);
                                if (i==99 && j==99) fprintf(stderr, "%d %d ca: %g pex: %g pey: %g cb: %g  diff: %g ([%d %d] %g [%d %d] %g [%d %d] %g) = %g, %g\n", i, j, ca, pex, pey, cb, -(mp->hz->data[i+j*xres]-mp->hz->data[(i-1)+j*xres]), i, j, mp->hz->data[i+j*xres], i-1, j, mp->hz->data[(i-1)+j*xres], i, j-1, mp->hz->data[i+(j-1)*xres], mp->ex->data[i+j*xres], mp->ey->data[i+j*xres]);


                }
                                if (i>=97 && i<103 && j>=57 && j<63) {
                                   printf("(%d %d %5.3g %+5.3e) ", i, j, cb, mp->ex->data[i+j*xres]);
                                }

                                if (i>=97 && i<103 && j==64) {
                                   printf(" ||   ");
                                }


                                if (i>=97 && i<103 && j>=97 && j<103) {
                                   printf("(%d %d %5.3g %+5.3e) ", i, j, cb, mp->ex->data[i+j*xres]);
                                }
                */
            } else { //tm mode
                mp->ez->data[i + j*xres] = ca*mp->ez->data[i + j*xres] + cb*((mp->hy->data[i + j*xres] - mp->hy->data[(i - 1) + j*xres]) -
                    (mp->hx->data[i + j*xres] - mp->hx->data[i + (j - 1)*xres]));
            }
        }

        /*
                if (i>=57 && i<63) {
                     printf("\n");
                }
                if (i==64) printf("\n");


                if (i>=97 && i<103) {
                     printf("\n");
                }
        */
    }
    if (set->sc.verbose)
        printf("done.\n");

    return 0;
}




/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */



