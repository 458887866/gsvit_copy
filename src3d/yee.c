
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

#if defined (_WIN32) && defined (_MSC_VER)
#include "../config_win.h"
#else
#include "../config.h"
#endif
#include <glib.h>
#include <stdlib.h>
#include "yee.h"
#include "constants.h"
#include "source.h"
#include "boundary.h"
#include <math.h>
#include <omp.h>



SvYeeData*
sv_yee_data_new(SvSet *set, gdouble dt)
{
    SvYeeData *d = (SvYeeData *)g_malloc(sizeof(SvYeeData));

    /*alloc other structures*/
    d->mat = NULL;
    d->dplrcx = NULL;
    d->dplrcy = NULL;
    d->dplrcz = NULL;
    d->plrcx = NULL;
    d->plrcy = NULL;
    d->plrcz = NULL;
    d->iplrcx = NULL;
    d->iplrcy = NULL;
    d->iplrcz = NULL;

    d->pxp = NULL;
    d->pyp = NULL;
    d->pzp = NULL;
    d->px = NULL;
    d->py = NULL;
    d->pz = NULL;

    d->dpxp = NULL;
    d->dpyp = NULL;
    d->dpzp = NULL;
    d->dpx = NULL;
    d->dpy = NULL;
    d->dpz = NULL;

    d->dt = dt;

    if (set->sc.verbose > 1)
        printf("SvYeeData initialized\n");

    return d;
}

void
sv_yee_data_allocate(SvYeeData *d, SvSet *set, gint xres, gint yres, gint zres, gdouble dx, gdouble dy, gdouble dz)
{
    d->ex = sv_dcube_new(xres, yres, zres, xres*dx, yres*dy, zres*dz, 1);
    d->ey = sv_dcube_new_alike(d->ex, 1);
    d->ez = sv_dcube_new_alike(d->ex, 1);
    d->hx = sv_dcube_new_alike(d->ex, 1);
    d->hy = sv_dcube_new_alike(d->ex, 1);
    d->hz = sv_dcube_new_alike(d->ex, 1);
    d->epsilon = NULL;
    d->mu = NULL;
    d->sigma = NULL;
    d->sigast = NULL;

    d->dx = dx;
    d->dy = dy;
    d->dz = dz;
    d->xres = xres;
    d->yres = yres;
    d->zres = zres;

    /*internal boundaries*/
    d->bnds = sv_icube_new(xres, yres, zres, xres*dx, yres*dy, zres*dz, 1);
    //sv_icube_fill(d->bnds, 0);

    /*if full material is requested, otherwise use matmode from settings*/
    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC) {
        d->epsilon = sv_fcube_new(xres, yres, zres, xres*dx, yres*dy, zres*dz, 1);
        d->sigma = sv_fcube_new_alike(d->epsilon, 1);
        sv_fcube_fill(d->epsilon, (gfloat)EPSILON_0);
    }
    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_MAGNETIC) {
        d->mu = sv_fcube_new(xres, yres, zres, xres*dx, yres*dy, zres*dz, 1);
        d->sigast = sv_fcube_new_alike(d->mu, 1);
        sv_fcube_fill(d->mu, (gfloat)MU_0);
    }

    if (set->sc.verbose > 1) 
        printf("Yee data allocated to %d %d %d\n", d->xres, d->yres, d->zres);

    return;
}

void
sv_yee_data_free(SvYeeData *d, SvSet *set)
{
    sv_dcube_free(d->ex);
    sv_dcube_free(d->ey);
    sv_dcube_free(d->ez);
    sv_dcube_free(d->hx);
    sv_dcube_free(d->hy);
    sv_dcube_free(d->hz);

    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC) {
        sv_fcube_free(d->epsilon);
        sv_fcube_free(d->sigma);
    }
    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_MAGNETIC) {
        sv_fcube_free(d->mu);
        sv_fcube_free(d->sigast);
    }
}

/*Yee step: calculate H from E*/
gint
sv_yee_data_ystep_h(SvYeeData *d, SvSetBoundary *sb, SvBoundary *bnd, SvSet *set, SvMatProp *mats, gint nmat)
{
    gint xnres = d->xres;
    gint matindex;
    gint ynres = d->yres;
    gint znres = d->zres;
    gint xnzero = 0;
    gint ynzero = 0;
    gint znzero = 0;
    gint i, j, k;
    gdouble da, dbx, dby, dbz;
    gdouble kappax, kappay, kappaz;
    gdouble ex, ey, ez;

    gdouble dt = d->dt;

    matindex = 0;
    i = j = k = 0;
    da = dbx = dby = dbz = 0;
    kappax = kappay = kappaz = 0;
    ex = ey = ez = 0;

    if (set->sc.verbose > 1) {
        printf("Running Yee H step...   ");
        fflush(stdout);
    }

#pragma omp parallel default(shared) private(i, j, k, matindex, kappax, kappay, kappaz, ex, ey, ez, da, dbx, dby, dbz)

#pragma omp for nowait

    for (i = xnzero; i < xnres; i++) {
        for (j = ynzero; j < ynres; j++) {
            for (k = znzero; k < znres; k++) {
                if (nmat > 0)
                    matindex = d->mat->data[i][j][k]; /*data are tabulated*/
                else
                    matindex = -1;

                /*CPML stretching*/
                kappax = 1;
                kappay = 1;
                kappaz = 1;

                if (sb) {
                   if (sb->bx0 == SV_BOUNDARY_CPML && i < sb->depth_bx0)
                       kappax = MAX(kappax, bnd->cpml.kappah_x0[i]);//

                   if (sb->bxn == SV_BOUNDARY_CPML && i >= (d->xres - sb->depth_bxn))
                       kappax = MAX(kappax, bnd->cpml.kappah_xn[i - (d->xres - sb->depth_bxn)]);

                   if (sb->by0 == SV_BOUNDARY_CPML && j < sb->depth_by0)
                       kappay = MAX(kappay, bnd->cpml.kappah_y0[j]);//

                   if (sb->byn == SV_BOUNDARY_CPML && j >= (d->yres - sb->depth_byn))
                       kappay = MAX(kappay, bnd->cpml.kappah_yn[j - (d->yres - sb->depth_byn)]);

                   if (sb->bz0 == SV_BOUNDARY_CPML && k < sb->depth_bz0)
                       kappaz = MAX(kappaz, bnd->cpml.kappah_z0[k]);//

                   if (sb->bzn == SV_BOUNDARY_CPML && k >= (d->zres - sb->depth_bzn))
                       kappaz = MAX(kappaz, bnd->cpml.kappah_zn[k - (d->zres - sb->depth_bzn)]);
                }
            

                ex = d->ex->data[i][j][k]; ey = d->ey->data[i][j][k]; ez = d->ez->data[i][j][k];

                /*determine local material parameters: copied from sv_pool_get_da_db() for speedup*/
                if (nmat > 0)
                    matindex = d->mat->data[i][j][k];
                else
                    matindex = -1;

                /*there is a local tabulated material in this point: use it always*/
                if (matindex > 0 && (mats[matindex].type == SV_MAT_LINTAB
                                     || mats[matindex].type == SV_MAT_DRUDE
                                     || mats[matindex].type == SV_MAT_PLRC
                                     || mats[matindex].type == SV_MAT_ADE
                                     || mats[matindex].type == SV_MAT_CP
                                     || mats[matindex].type == SV_MAT_CP3)) {
                    da = (1 - mats[matindex].sigast*dt / 2 / mats[matindex].mu) / (1 + mats[matindex].sigast*dt / 2 / mats[matindex].mu);
                    dbx = (dt / mats[matindex].mu / d->dx) / (1 + mats[matindex].sigast*dt / 2 / mats[matindex].mu);
                    dby = (dt / mats[matindex].mu / d->dy) / (1 + mats[matindex].sigast*dt / 2 / mats[matindex].mu);
                    dbz = (dt / mats[matindex].mu / d->dz) / (1 + mats[matindex].sigast*dt / 2 / mats[matindex].mu);
                }
                /*there is no local material, but there are some voxel-by-voxel data*/
                else if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_MAGNETIC) {

                    da = (1 - d->sigast->data[i][j][k] * dt / 2 / d->mu->data[i][j][k]) / (1 + d->sigast->data[i][j][k] * dt / 2 / d->mu->data[i][j][k]);
                    dbx = (dt / d->mu->data[i][j][k] / d->dx) / (1 + d->sigast->data[i][j][k] * dt / 2 / d->mu->data[i][j][k]);
                    dby = (dt / d->mu->data[i][j][k] / d->dy) / (1 + d->sigast->data[i][j][k] * dt / 2 / d->mu->data[i][j][k]);
                    dbz = (dt / d->mu->data[i][j][k] / d->dz) / (1 + d->sigast->data[i][j][k] * dt / 2 / d->mu->data[i][j][k]);
                } else {
                    da = 1;
                    dbx = dt / MU_0 / d->dx;
                    dby = dt / MU_0 / d->dy;
                    dbz = dt / MU_0 / d->dz;
                }

                /*update H*/
                if (!(j == (ynres - 1) || k == (znres - 1)))
                    d->hx->data[i][j][k] = da*d->hx->data[i][j][k] + dbx*((d->ey->data[i][j][k + 1] - ey) / kappaz -
                    (d->ez->data[i][j + 1][k] - ez) / kappay);

                if (!(i == (xnres - 1) || k == (znres - 1)))
                    d->hy->data[i][j][k] = da*d->hy->data[i][j][k] + dby*((d->ez->data[i + 1][j][k] - ez) / kappax -
                    (d->ex->data[i][j][k + 1] - ex) / kappaz);

                if (!(i == (xnres - 1) || j == (ynres - 1)))
                    d->hz->data[i][j][k] = da*d->hz->data[i][j][k] + dbz*((d->ex->data[i][j + 1][k] - ex) / kappay -
                    (d->ey->data[i + 1][j][k] - ey) / kappax);
            } // k
        } // j
    } // i

    if (set->sc.verbose > 1)
        printf("done.\n");

    return 0;
}

/*set of functions for PLRC*/
static inline gdouble
dch0(gdouble omegadiv, gdouble dt)
{
    return -omegadiv*omegadiv*(1.0 - exp(-LIGHT_SPEED*dt)) * (1.0 - exp(-LIGHT_SPEED*dt));

}

static inline gdouble
ch0(gdouble omegadiv, gdouble dt)
{
    return omegadiv*omegadiv*LIGHT_SPEED*dt - omegadiv*omegadiv*(1.0 - exp(-LIGHT_SPEED*dt));

}

static inline gdouble
dch1(gdouble omegadiv, gdouble dt)
{
    return -omegadiv*omegadiv * exp(-LIGHT_SPEED*dt) * (1.0 - exp(-LIGHT_SPEED*dt)) * (1.0 - exp(-LIGHT_SPEED*dt));
}

static void
cmult(gdouble a, gdouble b, gdouble c, gdouble d, gdouble *ra, gdouble *rb)
{
    *ra = a*c - b*d;
    *rb = a*d + b*c;
}

// RS - can this be replaced by cmult ?
static void
mcmult(gdouble a, gdouble b, gdouble c, gdouble d, gdouble *ra, gdouble *rb)
{
    gdouble pra, prb;

    pra = a*c - b*d;
    prb = a*d + b*c;

    *ra = pra;
    *rb = prb;
}

static void
cdiv(gdouble a, gdouble b, gdouble c, gdouble d, gdouble *ra, gdouble *rb)
{
    gdouble sub = c*c + d*d;
    //printf("%g %g %g %g %g\n", a, b, c, d, sub);
    *ra = (a*c + b*d) / sub;
    *rb = (b*c - a*d) / sub;
}

static inline void
fract(gdouble a, gdouble gamma, gdouble omega, gdouble phi, gdouble *re, gdouble *im)
{
    gdouble r_sub, i_sub;
    gdouble r_sup, i_sup;

    r_sub = gamma;
    i_sub = -omega;
    r_sup = -2 * a*omega*sin(phi);//was plus
    i_sup = -2 * a*omega*cos(phi);

    cdiv(r_sup*1e-16, i_sup*1e-16, r_sub*1e-16, i_sub*1e-16, re, im); //trying to prevent too large numbers division
}

//fdtd_metal_t.pdf
static void
update_plrc_cp3(gdouble *plrc, gdouble *iplrc, gdouble e, SvMatProp *prop, gdouble dt, gint n)
{
    gdouble r_efact, i_efact;
    gdouble r_sup, i_sup;

    r_efact = 1 - exp(-prop->cp3_gamma[n] * dt)*cos(prop->cp3_omega[n] * dt);
    i_efact = -exp(-prop->cp3_gamma[n] * dt)*sin(prop->cp3_omega[n] * dt);

    mcmult(r_efact, i_efact, r_efact, i_efact, &r_efact, &i_efact);//konecna zavorka na druhou

    fract(prop->cp3_a[n], prop->cp3_gamma[n], prop->cp3_omega[n], prop->cp3_phi[n], &r_sup, &i_sup);
    mcmult(r_sup, i_sup, r_efact, i_efact, &r_sup, &i_sup);

    mcmult(r_sup, i_sup, *plrc, *iplrc, &r_sup, &i_sup); //old multiplied

    r_efact = exp(-prop->cp3_gamma[n] * dt)*cos(prop->cp3_omega[n] * dt);
    i_efact = exp(-prop->cp3_gamma[n] * dt)*sin(prop->cp3_omega[n] * dt);

    *plrc = r_sup + r_efact*e;
    *iplrc = i_sup + i_efact*e;

    //printf("returning %g %g (%g %g %g)\n", *plrc, *iplrc, r_efact, i_efact, e);
    //
    //printf("n=%d, a %g, gamma %g omega %g phi %g\n", n, prop->cp3_a[n], prop->cp3_gamma[n], prop->cp3_omega[n], prop->cp3_phi[n]);
}

gdouble
xi(SvMatProp *prop, gdouble dt)
{
    gint n;
    gdouble r_efact, i_efact;
    gdouble r_sup, i_sup;
    gdouble sum_sup = 0;

    for (n = 0; n < 3; n++) {
        r_efact = 1 - exp(-prop->cp3_gamma[n] * dt)*cos(prop->cp3_omega[n] * dt);
        i_efact = -exp(-prop->cp3_gamma[n] * dt)*sin(prop->cp3_omega[n] * dt);

        fract(prop->cp3_a[n], prop->cp3_gamma[n], prop->cp3_omega[n], prop->cp3_phi[n], &r_sup, &i_sup);
        mcmult(r_sup, i_sup, r_efact, i_efact, &r_sup, &i_sup);

        sum_sup += r_sup;
    }

    //    printf("%g %g    %g\n", -prop->sigma*dt, sum_sup, -prop->sigma*dt + sum_sup);
    return -prop->sigma*dt + sum_sup;
}

void
dxxi(SvMatProp *prop, gdouble dt, gdouble *re, gdouble *im, gint n)
{
    gdouble r_efact, i_efact;
    gdouble r_sup, i_sup;

    r_efact = 1 - exp(-prop->cp3_gamma[n] * dt)*cos(prop->cp3_omega[n] * dt);
    i_efact = -exp(-prop->cp3_gamma[n] * dt)*sin(prop->cp3_omega[n] * dt);

    fract(prop->cp3_a[n], prop->cp3_gamma[n], prop->cp3_omega[n], prop->cp3_phi[n], &r_sup, &i_sup);
    cmult(r_efact, i_efact, r_efact, i_efact, &r_efact, &i_efact);
    cmult(r_sup, i_sup, r_efact, i_efact, &r_sup, &i_sup);

    *re = r_sup;
    *im = i_sup;
}

gdouble
xxi0(SvMatProp *prop, gdouble dt)
{
    gint n;
    gdouble r_efact, i_efact;
    gdouble omega, gamma;
    gdouble r_sup, i_sup;
    gdouble sum_sup = 0;

    for (n = 0; n < 2; n++) {
        r_efact = 1 - exp(-prop->cp3_gamma[n] * dt)*cos(prop->cp3_omega[n] * dt);
        i_efact = -exp(-prop->cp3_gamma[n] * dt)*sin(prop->cp3_omega[n] * dt);

        //        printf("data for material, point %d   a %g gamma %g omega %g phi %g\n", n, prop->cp3_a[n], prop->cp3_gamma[n], prop->cp3_omega[n], prop->cp3_phi[n]);
        //        printf("r/f %g %g\n", r_efact, i_efact);

        fract(prop->cp3_a[n], prop->cp3_gamma[n], prop->cp3_omega[n], prop->cp3_phi[n], &r_sup, &i_sup);

        //        printf("ff %g %g\n", r_sup, i_sup);

        cmult(r_sup, i_sup, r_efact, i_efact, &r_sup, &i_sup);

        //      printf("cres %g %g, I will be adding %g to %g\n", r_sup, i_sup, r_sup, sum_sup);

        sum_sup += r_sup;
    }

    omega = prop->drude_omega_p;
    gamma = prop->drude_nu;

    //    printf("xxi0: omega %g   gamma %g     %g %g %g,   result %g\n", omega, gamma, -(omega*omega/gamma/gamma)*(1-exp(-gamma*dt)), dt*omega*omega/gamma, sum_sup, -(omega*omega/gamma/gamma)*(1-exp(-gamma*dt)) + dt*omega*omega/gamma + sum_sup);
    return -(omega*omega / gamma / gamma)*(1 - exp(-gamma*dt)) + dt*omega*omega / gamma + (sum_sup);
}

static void
update_plrc_cp(gdouble *plrc, gdouble *iplrc, gdouble e, SvMatProp *prop, gdouble dt, gint n)
{
    gdouble r_efact, i_efact;
    gdouble r_exp, i_exp;
    gdouble r_res, i_res;

    r_exp = exp(-prop->cp3_gamma[n] * dt)*cos(prop->cp3_omega[n] * dt);
    i_exp = exp(-prop->cp3_gamma[n] * dt)*sin(prop->cp3_omega[n] * dt);
    cmult(r_exp, i_exp, *plrc, *iplrc, &r_res, &i_res);

    dxxi(prop, dt, &r_efact, &i_efact, n);

    *plrc = r_res + r_efact*e;
    *iplrc = i_res + i_efact*e;

    //printf("returning %g %g (%g %g %g)\n", *plrc, *iplrc, r_efact, i_efact, e);
    //
    //printf("n=%d, a %g, gamma %g omega %g phi %g\n", n, prop->cp3_a[n], prop->cp3_gamma[n], prop->cp3_omega[n], prop->cp3_phi[n]);
}

gint
sv_yee_data_get_tbc(SvYeeData *d, SvMatProp *mats, gint nmat, SvSet *set, gint i, gint j, gint k, gint *xtbc, gint *ytbc, gint *ztbc)
{
    gint matindex;
    (*xtbc) = (*ytbc) = (*ztbc) = 1; //ex, ey, ez component still to be computed after dispersive material treatment


 
    if (nmat > 0)
        matindex = d->mat->data[i][j][k];
    else {
        matindex = -1;
        return matindex;
    }

    /*decide which components to compute*/
    if (matindex > 0 &&                                            //do this if we are at dispersive material point
        (mats[matindex].type == SV_MAT_DRUDE || mats[matindex].type == SV_MAT_CP3 || mats[matindex].type == SV_MAT_CP ||
         mats[matindex].type == SV_MAT_ADE || mats[matindex].type == SV_MAT_PLRC)) {

        (*xtbc) = (*ytbc) = (*ztbc) = 0;

    } else {
        if (i > 0 && d->mat->data[i - 1][j][k] > 0 &&                       //or if we are at plus boundary with a dispersive material point
            (mats[d->mat->data[i - 1][j][k]].type == SV_MAT_DRUDE || mats[d->mat->data[i - 1][j][k]].type == SV_MAT_CP3 ||
             mats[d->mat->data[i - 1][j][k]].type == SV_MAT_CP ||
             mats[d->mat->data[i - 1][j][k]].type == SV_MAT_ADE || mats[d->mat->data[i - 1][j][k]].type == SV_MAT_PLRC)) {
            matindex = d->mat->data[i - 1][j][k];
            (*ytbc) = (*ztbc) = 0;
        }
        if (j > 0 && d->mat->data[i][j - 1][k] > 0 &&                       //or if we are at plus boundary with a dispersive material point
            (mats[d->mat->data[i][j - 1][k]].type == SV_MAT_DRUDE || mats[d->mat->data[i][j - 1][k]].type == SV_MAT_CP3 ||
             mats[d->mat->data[i][j - 1][k]].type == SV_MAT_CP ||
             mats[d->mat->data[i][j - 1][k]].type == SV_MAT_ADE || mats[d->mat->data[i][j - 1][k]].type == SV_MAT_PLRC)) {
            matindex = d->mat->data[i][j - 1][k];
            (*xtbc) = (*ztbc) = 0;
        }

        if (k > 0 && d->mat->data[i][j][k - 1] > 0 &&                       //or if we are at plus boundary with a dispersive material point
            (mats[d->mat->data[i][j][k - 1]].type == SV_MAT_DRUDE || mats[d->mat->data[i][j][k - 1]].type == SV_MAT_CP3 ||
             mats[d->mat->data[i][j][k - 1]].type == SV_MAT_CP ||
             mats[d->mat->data[i][j][k - 1]].type == SV_MAT_ADE || mats[d->mat->data[i][j][k - 1]].type == SV_MAT_PLRC)) {
            matindex = d->mat->data[i][j][k - 1];
            (*xtbc) = (*ytbc) = 0;
        }

        if (i > 0 && j > 0 && d->mat->data[i - 1][j - 1][k] > 0 &&                       //or if we are at plus boundary with a dispersive material point
            (mats[d->mat->data[i - 1][j - 1][k]].type == SV_MAT_DRUDE || mats[d->mat->data[i - 1][j - 1][k]].type == SV_MAT_CP3 ||
             mats[d->mat->data[i - 1][j - 1][k]].type == SV_MAT_CP ||
             mats[d->mat->data[i - 1][j - 1][k]].type == SV_MAT_ADE || mats[d->mat->data[i - 1][j - 1][k]].type == SV_MAT_PLRC)) {
            matindex = d->mat->data[i - 1][j - 1][k];
            *ztbc = 0;

        }
        if (j > 0 && k > 0 && d->mat->data[i][j - 1][k - 1] > 0 &&                       //or if we are at plus boundary with a dispersive material point
            (mats[d->mat->data[i][j - 1][k - 1]].type == SV_MAT_DRUDE || mats[d->mat->data[i][j - 1][k - 1]].type == SV_MAT_CP3 ||
             mats[d->mat->data[i][j - 1][k - 1]].type == SV_MAT_CP ||
             mats[d->mat->data[i][j - 1][k - 1]].type == SV_MAT_ADE || mats[d->mat->data[i][j - 1][k - 1]].type == SV_MAT_PLRC)) {
            matindex = d->mat->data[i][j - 1][k - 1];
            *xtbc = 0;
        }

        if (i > 0 && k > 0 && d->mat->data[i - 1][j][k - 1] > 0 &&                       //or if we are at plus boundary with a dispersive material point
            (mats[d->mat->data[i - 1][j][k - 1]].type == SV_MAT_DRUDE || mats[d->mat->data[i - 1][j][k - 1]].type == SV_MAT_CP3 ||
             mats[d->mat->data[i - 1][j][k - 1]].type == SV_MAT_CP ||
             mats[d->mat->data[i - 1][j][k - 1]].type == SV_MAT_ADE || mats[d->mat->data[i - 1][j][k - 1]].type == SV_MAT_PLRC)) {
            matindex = d->mat->data[i - 1][j][k - 1];
            *ytbc = 0;
        }
    }

    return matindex;

}

/*Yee step: calculate E from H*/
gint
sv_yee_data_ystep_e(SvYeeData *d, SvSetBoundary *sb, SvBoundary *bnd, SvSource *src, SvSet *set, SvMatProp *mats, gint nmat)
{
    gint xnres = d->xres;
    gint ynres = d->yres;
    gint znres = d->zres;
    gint matindex = 0;
    gint xnzero = 0;
    gint ynzero = 0;
    gint znzero = 0;
    gint i = 0, j = 0, k = 0, n = 0, ir = 0, jr = 0, kr = 0;
    gdouble ca = 1, cb = 0, cc = 0, dd = 0, ve = 0;
    gdouble cax = 1, cbx = 0, cay = 1, cby = 0, caz = 1, cbz = 0;
    gdouble dx = d->dx;
    gdouble dy = d->dy;
    gdouble dz = d->dz;
    gdouble pex = 0, pey = 0, pez = 0;
    gdouble kappax = 0, kappay = 0, kappaz = 0;
    gdouble hx = 0, hy = 0, hz = 0;
    gdouble vexp = 0, veyp = 0, vezp = 0;
    gdouble vdpx = 0, vdpy = 0, vdpz = 0;
    gdouble vpx[2] = {0, 0}, vpy[2] = {0, 0}, vpz[2] = {0, 0};
    gdouble sumxi0 = 0, sumchi0 = 0, sumpsix = 0, sumpsiy = 0, sumpsiz = 0;
    gdouble rc = 0, ic = 0;
    gdouble be = 0, ce = 0;
    gint xtbc = 0, ytbc = 0, ztbc = 0;
    gdouble vch0, vdch0, vdch1, epsilon, omega, gamma;

    gdouble dt = d->dt;

    if (set->sc.verbose > 1) {
        printf("Running Yee E step...   ");
        fflush(stdout);
    }    


#pragma omp parallel default(shared) private(i, j, k, n, matindex, kappax, kappay, kappaz, pex, pey, pez, ca, cb, cc, dd, ve, vexp, veyp, vezp, vdpx, vdpy, vdpz, vpx, vpy, vpz, sumxi0, sumchi0, sumpsix, sumpsiy, sumpsiz, rc, ic, be, ce, ir, jr, kr, xtbc, ytbc, ztbc, cax, cay, caz, cbx, cby, cbz, hx, hy, hz)

#pragma omp for nowait

    for (i = xnzero; i < xnres; i++) {
        for (j = ynzero; j < ynres; j++) {
            for (k = znzero; k < znres; k++) {
                /*CPML stretching*/
                kappax = 1;
                kappay = 1;
                kappaz = 1;

                cb = 0;

                if (sb) {
                   if (sb->bx0 == SV_BOUNDARY_CPML && i < sb->depth_bx0)
                       kappax = MAX(kappax, bnd->cpml.kappae_x0[i]);

                   if (sb->bxn == SV_BOUNDARY_CPML && i >= (d->xres - sb->depth_bxn))
                       kappax = MAX(kappax, bnd->cpml.kappae_xn[i - (d->xres - sb->depth_bxn)]);//

                   if (sb->by0 == SV_BOUNDARY_CPML && j < sb->depth_by0)
                       kappay = MAX(kappay, bnd->cpml.kappae_y0[j]);

                   if (sb->byn == SV_BOUNDARY_CPML && j >= (d->yres - sb->depth_byn))
                       kappay = MAX(kappay, bnd->cpml.kappae_yn[j - (d->yres - sb->depth_byn)]);//

                   if (sb->bz0 == SV_BOUNDARY_CPML && k < sb->depth_bz0)
                       kappaz = MAX(kappaz, bnd->cpml.kappae_z0[k]);

                   if (sb->bzn == SV_BOUNDARY_CPML && k >= (d->zres - sb->depth_bzn))
                       kappaz = MAX(kappaz, bnd->cpml.kappae_zn[k - (d->zres - sb->depth_bzn)]);
                }

                hx = d->hx->data[i][j][k];
                hy = d->hy->data[i][j][k];
                hz = d->hz->data[i][j][k];

                matindex = sv_yee_data_get_tbc(d, mats, nmat, set, i, j, k, &xtbc, &ytbc, &ztbc);

                if (nmat > 0 && (d->plrcx || d->dplrcx || d->pxp)) {      //do this if any dispersive material exists
                    if (xtbc*ytbc*ztbc == 0) { //not everything left for future, so compute something for this pixel, at least for some component
                        pex = d->ex->data[i][j][k];
                        pey = d->ey->data[i][j][k];
                        pez = d->ez->data[i][j][k];

                        if (mats[matindex].type == SV_MAT_DRUDE) {

                            //Enhanced optical transmission through a silver plate with a slit array
                            //Yong Fu, Kang Li , Fanmin Kong
                            //Optik 121 (2010) 259–262

                            vch0 = ch0(mats[matindex].drude_omega_p/LIGHT_SPEED, dt);
                            vdch1 = dch1(mats[matindex].drude_omega_p/LIGHT_SPEED, dt);
                            vdch0 = dch0(mats[matindex].drude_omega_p/LIGHT_SPEED, dt);
                            epsilon = mats[matindex].epsilon/EPSILON_0;

                            if (j > 0 && k > 0 && xtbc == 0) d->ex->data[i][j][k] = (epsilon+vdch0)/(epsilon+vch0)*pex
                                + 1.0/(epsilon+vch0)*d->dplrcx[0]->data[i][j][k]
                                + dt/(epsilon+vch0)/EPSILON_0/dx
                                *((d->hz->data[i][j][k] - d->hz->data[i][j-1][k]) - (d->hy->data[i][j][k] - d->hy->data[i][j][k-1]));

                            if (i > 0 && k > 0 && ytbc == 0) d->ey->data[i][j][k] = (epsilon+vdch0)/(epsilon+vch0)*pey
                                + 1.0/(epsilon+vch0)*d->dplrcy[0]->data[i][j][k]
                                + dt/(epsilon+vch0)/EPSILON_0/dx
                                *((d->hx->data[i][j][k] - d->hx->data[i][j][k-1]) - (d->hz->data[i][j][k] - d->hz->data[i-1][j][k]));

                            if (i > 0 && j > 0 && ztbc == 0) d->ez->data[i][j][k] = (epsilon+vdch0)/(epsilon+vch0)*pez
                                + 1.0/(epsilon+vch0)*d->dplrcz[0]->data[i][j][k]
                                + dt/(epsilon+vch0)/EPSILON_0/dx
                                *((d->hy->data[i][j][k] - d->hy->data[i-1][j][k]) - (d->hx->data[i][j][k] - d->hx->data[i][j-1][k]));

                            d->dplrcx[0]->data[i][j][k] = exp(-mats[matindex].drude_nu*dt) * d->dplrcx[0]->data[i][j][k]
                                + vdch1*pex;

                            d->dplrcy[0]->data[i][j][k] = exp(-mats[matindex].drude_nu*dt) * d->dplrcy[0]->data[i][j][k]
                                + vdch1*pey;

                            d->dplrcz[0]->data[i][j][k] = exp(-mats[matindex].drude_nu*dt) * d->dplrcz[0]->data[i][j][k]
                                + vdch1*pez;

                        } 
                        else if (mats[matindex].type == SV_MAT_CP) {
                            // A. Vial, J. Opt. A: Pure Appl. Opt. 9 (2007) 745–748

                            omega = mats[matindex].drude_omega_p;
                            gamma = mats[matindex].drude_nu;
                            epsilon = mats[matindex].epsilon/EPSILON_0;

                            d->dplrcx[0]->data[i][j][k] = exp(-gamma*dt) * d->dplrcx[0]->data[i][j][k]
                                - (omega*omega/gamma/gamma) * (1-exp(-gamma*dt)) * (1-exp(-gamma*dt))*pex;

                            d->dplrcy[0]->data[i][j][k] = exp(-gamma*dt) * d->dplrcy[0]->data[i][j][k]
                                - (omega*omega/gamma/gamma) * (1-exp(-gamma*dt)) * (1-exp(-gamma*dt))*pey;

                            d->dplrcz[0]->data[i][j][k] = exp(-gamma*dt) * d->dplrcz[0]->data[i][j][k]
                                - (omega*omega/gamma/gamma) * (1-exp(-gamma*dt)) * (1-exp(-gamma*dt))*pez;

                            for (n = 0; n < 2; n++) {
                                update_plrc_cp(&(d->plrcx[n]->data[i][j][k]), &(d->iplrcx[n]->data[i][j][k]), d->ex->data[i][j][k],
                                               &(mats[matindex]), dt, n);

                                update_plrc_cp(&(d->plrcy[n]->data[i][j][k]), &(d->iplrcy[n]->data[i][j][k]), d->ey->data[i][j][k],
                                               &(mats[matindex]), dt, n);

                                update_plrc_cp(&(d->plrcz[n]->data[i][j][k]), &(d->iplrcz[n]->data[i][j][k]), d->ez->data[i][j][k],
                                               &(mats[matindex]), dt, n);
                            }

                            ca = epsilon / (epsilon + xxi0(&(mats[matindex]), dt));
                            cb = dt / (dx*EPSILON_0 * (epsilon + xxi0(&(mats[matindex]), dt)));
                            cc = 1.0 / (epsilon + xxi0(&(mats[matindex]), dt));

                            if (j > 0 && k > 0 && xtbc == 0) d->ex->data[i][j][k] = ca*pex
                                + cc*(d->dplrcx[0]->data[i][j][k] + d->plrcx[0]->data[i][j][k]+d->plrcx[1]->data[i][j][k])
                                + cb*((d->hz->data[i][j][k] - d->hz->data[i][j-1][k]) - (d->hy->data[i][j][k] - d->hy->data[i][j][k-1]));

                            if (i > 0 && k > 0 && ytbc == 0) d->ey->data[i][j][k] = ca*pey
                                + cc*(d->dplrcy[0]->data[i][j][k] + d->plrcy[0]->data[i][j][k]+d->plrcy[1]->data[i][j][k])
                                + cb*((d->hx->data[i][j][k] - d->hx->data[i][j][k-1]) - (d->hz->data[i][j][k] - d->hz->data[i-1][j][k]));

                            if (i > 0 && j > 0 && ztbc == 0) d->ez->data[i][j][k] = ca*pez
                                + cc*(d->dplrcz[0]->data[i][j][k] + d->plrcz[0]->data[i][j][k]+d->plrcz[1]->data[i][j][k])
                                + cb*((d->hy->data[i][j][k] - d->hy->data[i-1][j][k]) - (d->hx->data[i][j][k] - d->hx->data[i][j-1][k]));
                        } // else if SV_MAT_CP


                        else if (mats[matindex].type == SV_MAT_CP3) {
                            //J. Y. Lu, Y. H. Chang, Superlattices and Microstructures, 47 (2010) 60

                            epsilon = mats[matindex].epsilon/EPSILON_0;

                            for (n = 0; n < 3; n++) {
                                update_plrc_cp3(&(d->plrcx[n]->data[i][j][k]), &(d->iplrcx[n]->data[i][j][k]),
                                                d->ex->data[i][j][k], &(mats[matindex]), dt, n);

                                update_plrc_cp3(&(d->plrcy[n]->data[i][j][k]), &(d->iplrcy[n]->data[i][j][k]),
                                                d->ey->data[i][j][k], &(mats[matindex]), dt, n);

                                update_plrc_cp3(&(d->plrcz[n]->data[i][j][k]), &(d->iplrcz[n]->data[i][j][k]),
                                                d->ez->data[i][j][k], &(mats[matindex]), dt, n);
                            }

                            ca = epsilon / (epsilon + xi(&(mats[matindex]), dt));
                            cb = dt / (dx*EPSILON_0 * (epsilon + xi(&(mats[matindex]), dt)));
                            cc = 1.0 / (epsilon + xi(&(mats[matindex]), dt));

                            if (j > 0 && k > 0 && xtbc == 0) d->ex->data[i][j][k] = ca*pex
                                + cc*(d->plrcx[0]->data[i][j][k]+d->plrcx[1]->data[i][j][k]+d->plrcx[2]->data[i][j][k])
                                + cb*((d->hz->data[i][j][k] - d->hz->data[i][j-1][k]) - (d->hy->data[i][j][k] - d->hy->data[i][j][k-1]));

                            if (i > 0 && k > 0 && ytbc == 0) d->ey->data[i][j][k] = ca*pey
                                + cc*(d->plrcy[0]->data[i][j][k]+d->plrcy[1]->data[i][j][k]+d->plrcy[2]->data[i][j][k])
                                + cb*((d->hx->data[i][j][k] - d->hx->data[i][j][k-1]) - (d->hz->data[i][j][k] -d->hz->data[i-1][j][k]));

                            if (i > 0 && j > 0 && ztbc == 0) d->ez->data[i][j][k] = ca*pez
                                + cc*(d->plrcz[0]->data[i][j][k]+d->plrcz[1]->data[i][j][k]+d->plrcz[2]->data[i][j][k])
                                + cb*((d->hy->data[i][j][k] - d->hy->data[i-1][j][k]) - (d->hx->data[i][j][k] - d->hx->data[i][j-1][k]));

                        } // else if SV_MAT_CP3

                        //K. Chun, H. Kim, H. Kim, Y Chung, PIERS 135, 373-390, 2013
                        else if (mats[matindex].type == SV_MAT_PLRC) {
                            /*sum plrc terms*/
                            sumxi0 = mats[matindex].plrc_d_xi + mats[matindex].plrc_p_xi[0] + mats[matindex].plrc_p_xi[1];
                            sumchi0 = mats[matindex].plrc_d_chi + mats[matindex].plrc_p_chi[0] + mats[matindex].plrc_p_chi[1];
                            sumpsix = d->dplrcx[0]->data[i][j][k] + d->plrcx[0]->data[i][j][k] + d->plrcx[1]->data[i][j][k];
                            sumpsiy = d->dplrcy[0]->data[i][j][k] + d->plrcy[0]->data[i][j][k] + d->plrcy[1]->data[i][j][k];
                            sumpsiz = d->dplrcz[0]->data[i][j][k] + d->plrcz[0]->data[i][j][k] + d->plrcz[1]->data[i][j][k];

                            ca = (2 * EPSILON_0 * (mats[matindex].epsilon / EPSILON_0 - sumxi0) - mats[matindex].sigma * dt) /
                                (2 * EPSILON_0 * (mats[matindex].epsilon / EPSILON_0 - sumxi0 + sumchi0) + mats[matindex].sigma * dt);
                            cb = 2 * EPSILON_0 / (2 * EPSILON_0 * (mats[matindex].epsilon / EPSILON_0 - sumxi0 + sumchi0) + mats[matindex].sigma * dt);
                            cc = 2 * dt / (2 * EPSILON_0 * (mats[matindex].epsilon / EPSILON_0 - sumxi0 + sumchi0) + mats[matindex].sigma * dt);

                            /*update electric field*/
                            if (j > 0 && k > 0 && xtbc == 0)
                                d->ex->data[i][j][k] = ca*pex
                                + cb*sumpsix
                                + (cc / dx)*((d->hz->data[i][j][k] - d->hz->data[i][j - 1][k]) / kappay - (d->hy->data[i][j][k] - d->hy->data[i][j][k - 1]) / kappaz);

                            if (i > 0 && k > 0 && ytbc == 0)
                                d->ey->data[i][j][k] = ca*pey
                                + cb*sumpsiy
                                + (cc / dy)*((d->hx->data[i][j][k] - d->hx->data[i][j][k - 1]) / kappaz - (d->hz->data[i][j][k] - d->hz->data[i - 1][j][k]) / kappax);

                            if (i > 0 && j > 0 && ztbc == 0)
                                d->ez->data[i][j][k] = ca*pez
                                + cb*sumpsiz
                                + (cc / dz)*((d->hy->data[i][j][k] - d->hy->data[i - 1][j][k]) / kappax - (d->hx->data[i][j][k] - d->hx->data[i][j - 1][k]) / kappay);


                            /*try to do cpml update here*/
                            if (sb) {
                               if (sb->bx0 == SV_BOUNDARY_CPML && i > 1 && i < sb->depth_bx0) {
                                  cby = (dt / mats[matindex].epsilon / dy) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  cbz = (dt / mats[matindex].epsilon / dz) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  be = bnd->cpml.be_x0[i];
                                  ce = bnd->cpml.ce_x0[i];

                                  bnd->cpml.peyx_x0->data[i][j][k] = be * bnd->cpml.peyx_x0->data[i][j][k]
                                      + ce * ((d->hz->data[i][j][k] - d->hz->data[i - 1][j][k])) / dy;
                                  bnd->cpml.pezx_x0->data[i][j][k] = be * bnd->cpml.pezx_x0->data[i][j][k]
                                      + ce * ((d->hy->data[i][j][k] - d->hy->data[i - 1][j][k])) / dz;

                                  if (!(j == 0 || k == 0) && ytbc == 0)
                                      d->ey->data[i][j][k] -= dy * cby * bnd->cpml.peyx_x0->data[i][j][k];
                                  if (!(j == 0 || k == 0) && ztbc == 0)
                                      d->ez->data[i][j][k] += dz * cbz * bnd->cpml.pezx_x0->data[i][j][k];
                              }

                              if (sb->by0 == SV_BOUNDARY_CPML && j > 1 && j < sb->depth_by0) {
                                  cbx = (dt / mats[matindex].epsilon / dx) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  cbz = (dt / mats[matindex].epsilon / dz) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  be = bnd->cpml.be_y0[j];
                                  ce = bnd->cpml.ce_y0[j];

                                  bnd->cpml.pexy_y0->data[i][j][k] = be*bnd->cpml.pexy_y0->data[i][j][k]
                                      + ce * ((d->hz->data[i][j][k] - d->hz->data[i][j - 1][k])) / dx;
                                  bnd->cpml.pezy_y0->data[i][j][k] = be*bnd->cpml.pezy_y0->data[i][j][k]
                                      + ce * ((d->hx->data[i][j][k] - d->hx->data[i][j - 1][k])) / dz;

                                  if (!(i == 0 || k == 0) && xtbc == 0)
                                      d->ex->data[i][j][k] += dx * cbx * bnd->cpml.pexy_y0->data[i][j][k];
                                  if (!(i == 0 || k == 0) && ztbc == 0)
                                      d->ez->data[i][j][k] -= dz * cbz * bnd->cpml.pezy_y0->data[i][j][k];

                               }

                               if (sb->bz0 == SV_BOUNDARY_CPML && k > 1 && k < sb->depth_bz0) {
                                  cbx = (dt / mats[matindex].epsilon / dx) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  cby = (dt / mats[matindex].epsilon / dy) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  be = bnd->cpml.be_z0[k];
                                  ce = bnd->cpml.ce_z0[k];

                                  bnd->cpml.peyz_z0->data[i][j][k] = be*bnd->cpml.peyz_z0->data[i][j][k]
                                      + ce * ((d->hx->data[i][j][k] - d->hx->data[i][j][k - 1])) / dy;
                                  bnd->cpml.pexz_z0->data[i][j][k] = be*bnd->cpml.pexz_z0->data[i][j][k]
                                      + ce * ((d->hy->data[i][j][k] - d->hy->data[i][j][k - 1])) / dx;

                                  if (!(i == 0 || j == 0) && ytbc == 0)
                                      d->ey->data[i][j][k] += dy * cby * bnd->cpml.peyz_z0->data[i][j][k];
                                  if (!(i == 0 || j == 0) && xtbc == 0)
                                      d->ex->data[i][j][k] -= dx * cbx * bnd->cpml.pexz_z0->data[i][j][k];
                               }

                               if (sb->bxn == SV_BOUNDARY_CPML && (i >= (set->sp.xres - sb->depth_bxn + 1)) && i < set->sp.xres) {
                                  ir = i - (set->sp.xres - sb->depth_bxn);
                                  cby = (dt / mats[matindex].epsilon / dy) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  cbz = (dt / mats[matindex].epsilon / dz) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  be = bnd->cpml.be_xn[ir];
                                  ce = bnd->cpml.ce_xn[ir];


                                  bnd->cpml.peyx_xn->data[ir][j][k] = be*bnd->cpml.peyx_xn->data[ir][j][k]
                                      + ce * ((d->hz->data[i][j][k] - d->hz->data[i - 1][j][k])) / dy;
                                  bnd->cpml.pezx_xn->data[ir][j][k] = be*bnd->cpml.pezx_xn->data[ir][j][k]
                                      + ce * ((d->hy->data[i][j][k] - d->hy->data[i - 1][j][k])) / dz;

                                  if (!(j == 0 || k == 0) && ytbc == 0)
                                      d->ey->data[i][j][k] -= dy * cby * bnd->cpml.peyx_xn->data[ir][j][k];
                                  if (!(j == 0 || k == 0) && ztbc == 0)
                                      d->ez->data[i][j][k] += dz * cbz * bnd->cpml.pezx_xn->data[ir][j][k];
                               }

                               if (sb->byn == SV_BOUNDARY_CPML && (j >= (set->sp.yres - sb->depth_byn + 1)) && j < set->sp.yres) {
                                  jr = j - (set->sp.yres - sb->depth_byn);
                                  cbx = (dt / mats[matindex].epsilon / dx) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  cbz = (dt / mats[matindex].epsilon / dz) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  be = bnd->cpml.be_yn[jr];
                                  ce = bnd->cpml.ce_yn[jr];

                                  bnd->cpml.pexy_yn->data[i][jr][k] = be*bnd->cpml.pexy_yn->data[i][jr][k]
                                      + ce * ((d->hz->data[i][j][k] - d->hz->data[i][j - 1][k])) / dx;
                                  bnd->cpml.pezy_yn->data[i][jr][k] = be*bnd->cpml.pezy_yn->data[i][jr][k]
                                      + ce * ((d->hx->data[i][j][k] - d->hx->data[i][j - 1][k])) / dz;

                                  if (!(i == 0 || k == 0) && xtbc == 0)
                                      d->ex->data[i][j][k] += dx * cbx * bnd->cpml.pexy_yn->data[i][jr][k];
                                  if (!(i == 0 || k == 0) && ztbc == 0)
                                      d->ez->data[i][j][k] -= dz * cbz * bnd->cpml.pezy_yn->data[i][jr][k];
                               }
  
                              if (sb->bzn == SV_BOUNDARY_CPML && (k >= (set->sp.zres - sb->depth_bzn + 1)) && k < set->sp.zres) {
                                  kr = k - (set->sp.zres - sb->depth_bzn);
                                  cbx = (dt / mats[matindex].epsilon / dx) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  cby = (dt / mats[matindex].epsilon / dy) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);

                                  be = bnd->cpml.be_zn[kr];
                                  ce = bnd->cpml.ce_zn[kr];

                                  bnd->cpml.peyz_zn->data[i][j][kr] = be*bnd->cpml.peyz_zn->data[i][j][kr]
                                      + ce * ((d->hx->data[i][j][k] - d->hx->data[i][j][k - 1])) / dy;
                                  bnd->cpml.pexz_zn->data[i][j][kr] = be*bnd->cpml.pexz_zn->data[i][j][kr]
                                      + ce * ((d->hy->data[i][j][k] - d->hy->data[i][j][k - 1])) / dx;

                                  if (!(i == 0 || j == 0) && ytbc == 0)
                                      d->ey->data[i][j][k] += dy * cby * bnd->cpml.peyz_zn->data[i][j][kr];
                                  if (!(i == 0 || j == 0) && xtbc == 0)
                                      d->ex->data[i][j][k] -= dx * cbx * bnd->cpml.pexz_zn->data[i][j][kr];
                              }
                            }
                            /*end of cpml update*/


                            /*update plrc terms*/

                            d->dplrcx[0]->data[i][j][k] = (mats[matindex].plrc_d_dchi - mats[matindex].plrc_d_dxi) * d->ex->data[i][j][k]
                                + mats[matindex].plrc_d_dxi * pex + exp(-mats[matindex].drude_nu*dt) * d->dplrcx[0]->data[i][j][k];

                            d->dplrcy[0]->data[i][j][k] = (mats[matindex].plrc_d_dchi - mats[matindex].plrc_d_dxi) * d->ey->data[i][j][k]
                                + mats[matindex].plrc_d_dxi * pey + exp(-mats[matindex].drude_nu*dt) * d->dplrcy[0]->data[i][j][k];

                            d->dplrcz[0]->data[i][j][k] = (mats[matindex].plrc_d_dchi - mats[matindex].plrc_d_dxi) * d->ez->data[i][j][k]
                                + mats[matindex].plrc_d_dxi * pez + exp(-mats[matindex].drude_nu*dt) * d->dplrcz[0]->data[i][j][k];

                            //d->dplrcx[0]->data[i][j][k] = d->dplrcy[0]->data[i][j][k] = d->dplrcz[0]->data[i][j][k] = 0;

                            for (n = 0; n < 2; n++) {
                                cmult(exp(-dt*mats[matindex].cp3_gamma[n]) * cos(dt*mats[matindex].cp3_omega[n]),
                                      -exp(-dt*mats[matindex].cp3_gamma[n]) * sin(dt*mats[matindex].cp3_omega[n]),
                                      d->plrcx[n]->data[i][j][k], d->iplrcx[n]->data[i][j][k], &rc, &ic);

                                d->plrcx[n]->data[i][j][k] = (mats[matindex].plrc_p_dchir[n] - mats[matindex].plrc_p_dxir[n]) * d->ex->data[i][j][k]
                                    + mats[matindex].plrc_p_dxir[n] * pex + rc;
                                d->iplrcx[n]->data[i][j][k] = (mats[matindex].plrc_p_dchii[n] - mats[matindex].plrc_p_dxii[n]) * d->ex->data[i][j][k]
                                    + mats[matindex].plrc_p_dxii[n] * pex + ic;

                                cmult(exp(-dt*mats[matindex].cp3_gamma[n]) * cos(dt*mats[matindex].cp3_omega[n]),
                                      -exp(-dt*mats[matindex].cp3_gamma[n]) * sin(dt*mats[matindex].cp3_omega[n]),
                                      d->plrcy[n]->data[i][j][k], d->iplrcy[n]->data[i][j][k], &rc, &ic);

                                d->plrcy[n]->data[i][j][k] = (mats[matindex].plrc_p_dchir[n] - mats[matindex].plrc_p_dxir[n]) * d->ey->data[i][j][k]
                                    + mats[matindex].plrc_p_dxir[n] * pey + rc;
                                d->iplrcy[n]->data[i][j][k] = (mats[matindex].plrc_p_dchii[n] - mats[matindex].plrc_p_dxii[n]) * d->ey->data[i][j][k]
                                    + mats[matindex].plrc_p_dxii[n] * pey + ic;

                                cmult(exp(-dt*mats[matindex].cp3_gamma[n]) * cos(dt*mats[matindex].cp3_omega[n]),
                                      -exp(-dt*mats[matindex].cp3_gamma[n]) * sin(dt*mats[matindex].cp3_omega[n]),
                                      d->plrcz[n]->data[i][j][k], d->iplrcz[n]->data[i][j][k], &rc, &ic);

                                d->plrcz[n]->data[i][j][k] = (mats[matindex].plrc_p_dchir[n] - mats[matindex].plrc_p_dxir[n]) * d->ez->data[i][j][k]
                                    + mats[matindex].plrc_p_dxir[n] * pez + rc;
                                d->iplrcz[n]->data[i][j][k] = (mats[matindex].plrc_p_dchii[n] - mats[matindex].plrc_p_dxii[n]) * d->ez->data[i][j][k]
                                    + mats[matindex].plrc_p_dxii[n] * pez + ic;

                                // d->plrcx[n]->data[i][j][k] = d->iplrcx[n]->data[i][j][k] = d->plrcy[n]->data[i][j][k] = d->iplrcy[n]->data[i][j][k] = d->plrcz[n]->data[i][j][k] = d->iplrcz[n]->data[i][j][k] = 0;
                            } // n
                        } // else if SV_MAT_PLRC			
                        else if (mats[matindex].type == SV_MAT_ADE) {
                            //K. Chun, H. Kim, H. Kim, Y Chung, PIERS 135, 373-390, 2013

                            /*store temporarily e n-1*/
                            vexp = d->ex->data[i][j][k];
                            veyp = d->ey->data[i][j][k];
                            vezp = d->ez->data[i][j][k];

                            /*update E*/
                            if (j > 0 && k > 0 && xtbc == 0) {
                                d->ex->data[i][j][k] = mats[matindex].ade_c0 *
                                    ((d->hz->data[i][j][k] - d->hz->data[i][j - 1][k]) / kappay - (d->hy->data[i][j][k] - d->hy->data[i][j][k - 1]) / kappaz) / dx
                                    + mats[matindex].ade_c1*(-mats[matindex].ade_a0*d->dpxp[0]->data[i][j][k] + (1 - mats[matindex].ade_a1)*d->dpx[0]->data[i][j][k])
                                    + mats[matindex].ade_c1*(-mats[matindex].ade_bp0[0] * d->pxp[0]->data[i][j][k] + (1 - mats[matindex].ade_bp1[0])*d->px[0]->data[i][j][k])
                                    + mats[matindex].ade_c1*(-mats[matindex].ade_bp0[1] * d->pxp[1]->data[i][j][k] + (1 - mats[matindex].ade_bp1[1])*d->px[1]->data[i][j][k])
                                    + mats[matindex].ade_c2*d->exp->data[i][j][k] + mats[matindex].ade_c3*d->ex->data[i][j][k];
                            }

                            if (i > 0 && k > 0 && ytbc == 0) {
                                d->ey->data[i][j][k] = mats[matindex].ade_c0 *
                                    ((d->hx->data[i][j][k] - d->hx->data[i][j][k - 1]) / kappaz - (d->hz->data[i][j][k] - d->hz->data[i - 1][j][k]) / kappax) / dy
                                    + mats[matindex].ade_c1*(-mats[matindex].ade_a0*d->dpyp[0]->data[i][j][k] + (1 - mats[matindex].ade_a1)*d->dpy[0]->data[i][j][k])
                                    + mats[matindex].ade_c1*(-mats[matindex].ade_bp0[0] * d->pyp[0]->data[i][j][k] + (1 - mats[matindex].ade_bp1[0])*d->py[0]->data[i][j][k])
                                    + mats[matindex].ade_c1*(-mats[matindex].ade_bp0[1] * d->pyp[1]->data[i][j][k] + (1 - mats[matindex].ade_bp1[1])*d->py[1]->data[i][j][k])
                                    + mats[matindex].ade_c2*d->eyp->data[i][j][k] + mats[matindex].ade_c3*d->ey->data[i][j][k];
                            }

                            if (i > 0 && j > 0 && ztbc == 0) {
                                d->ez->data[i][j][k] = mats[matindex].ade_c0 *
                                    ((d->hy->data[i][j][k] - d->hy->data[i - 1][j][k]) / kappax - (d->hx->data[i][j][k] - d->hx->data[i][j - 1][k]) / kappay) / dz
                                    + mats[matindex].ade_c1*(-mats[matindex].ade_a0*d->dpzp[0]->data[i][j][k] + (1 - mats[matindex].ade_a1)*d->dpz[0]->data[i][j][k])
                                    + mats[matindex].ade_c1*(-mats[matindex].ade_bp0[0] * d->pzp[0]->data[i][j][k] + (1 - mats[matindex].ade_bp1[0])*d->pz[0]->data[i][j][k])
                                    + mats[matindex].ade_c1*(-mats[matindex].ade_bp0[1] * d->pzp[1]->data[i][j][k] + (1 - mats[matindex].ade_bp1[1])*d->pz[1]->data[i][j][k])
                                    + mats[matindex].ade_c2*d->ezp->data[i][j][k] + mats[matindex].ade_c3*d->ez->data[i][j][k];
                            }

                            /*try to do cpml update here*/
                            if (sb) {
                              if (sb->bx0 == SV_BOUNDARY_CPML && i > 1 && i < sb->depth_bx0) {
                                  cby = (dt / mats[matindex].epsilon / dy) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  cbz = (dt / mats[matindex].epsilon / dz) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  be = bnd->cpml.be_x0[i];
                                  ce = bnd->cpml.ce_x0[i];

                                  bnd->cpml.peyx_x0->data[i][j][k] = be * bnd->cpml.peyx_x0->data[i][j][k]
                                      + ce * ((d->hz->data[i][j][k] - d->hz->data[i - 1][j][k])) / dy;
                                  bnd->cpml.pezx_x0->data[i][j][k] = be * bnd->cpml.pezx_x0->data[i][j][k]
                                      + ce * ((d->hy->data[i][j][k] - d->hy->data[i - 1][j][k])) / dz;

                                  if (!(j == 0 || k == 0) && ytbc == 0)
                                      d->ey->data[i][j][k] -= dy * cb * bnd->cpml.peyx_x0->data[i][j][k];
                                  if (!(j == 0 || k == 0) && ztbc == 0)
                                      d->ez->data[i][j][k] += dz * cb * bnd->cpml.pezx_x0->data[i][j][k];
                              }

                              if (sb->by0 == SV_BOUNDARY_CPML && j > 1 && j < sb->depth_by0) {
                                  cbx = (dt / mats[matindex].epsilon / dx) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  cbz = (dt / mats[matindex].epsilon / dz) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  be = bnd->cpml.be_y0[j];
                                  ce = bnd->cpml.ce_y0[j];

                                  bnd->cpml.pexy_y0->data[i][j][k] = be*bnd->cpml.pexy_y0->data[i][j][k]
                                      + ce * ((d->hz->data[i][j][k] - d->hz->data[i][j - 1][k])) / dz;
                                  bnd->cpml.pezy_y0->data[i][j][k] = be*bnd->cpml.pezy_y0->data[i][j][k]
                                      + ce * ((d->hx->data[i][j][k] - d->hx->data[i][j - 1][k])) / dx;

                                  if (!(i == 0 || k == 0) && xtbc == 0)
                                      d->ex->data[i][j][k] += dx * cbx * bnd->cpml.pexy_y0->data[i][j][k];
                                  if (!(i == 0 || k == 0) && ztbc == 0)
                                      d->ez->data[i][j][k] -= dz * cbz * bnd->cpml.pezy_y0->data[i][j][k];
                              }

                              if (sb->bz0 == SV_BOUNDARY_CPML && k > 1 && k < sb->depth_bz0) {
                                  cbx = (dt / mats[matindex].epsilon / dx) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  cby = (dt / mats[matindex].epsilon / dy) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  be = bnd->cpml.be_z0[k];
                                  ce = bnd->cpml.ce_z0[k];

                                  bnd->cpml.peyz_z0->data[i][j][k] = be*bnd->cpml.peyz_z0->data[i][j][k]
                                      + ce * ((d->hx->data[i][j][k] - d->hx->data[i][j][k - 1])) / dy;
                                  bnd->cpml.pexz_z0->data[i][j][k] = be*bnd->cpml.pexz_z0->data[i][j][k]
                                      + ce * ((d->hy->data[i][j][k] - d->hy->data[i][j][k - 1])) / dx;

                                  if (!(i == 0 || j == 0) && ytbc == 0)
                                      d->ey->data[i][j][k] += dy * cby * bnd->cpml.peyz_z0->data[i][j][k];
                                  if (!(i == 0 || j == 0) && xtbc == 0)
                                      d->ex->data[i][j][k] -= dx * cbx * bnd->cpml.pexz_z0->data[i][j][k];
                              }

                              if (sb->bxn == SV_BOUNDARY_CPML && (i >= (set->sp.xres - sb->depth_bxn + 1)) && i < set->sp.xres) {
                                  ir = i - (set->sp.xres - sb->depth_bxn);
                                  cby = (dt / mats[matindex].epsilon / dy) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  cbz = (dt / mats[matindex].epsilon / dz) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  be = bnd->cpml.be_xn[ir];
                                  ce = bnd->cpml.ce_xn[ir];

                                  bnd->cpml.peyx_xn->data[ir][j][k] = be*bnd->cpml.peyx_xn->data[ir][j][k]
                                      + ce * ((d->hz->data[i][j][k] - d->hz->data[i - 1][j][k])) / dy;
                                  bnd->cpml.pezx_xn->data[ir][j][k] = be*bnd->cpml.pezx_xn->data[ir][j][k]
                                      + ce * ((d->hy->data[i][j][k] - d->hy->data[i - 1][j][k])) / dz;

                                  if (!(j == 0 || k == 0) && ytbc == 0)
                                      d->ey->data[i][j][k] -= dy * cby * bnd->cpml.peyx_xn->data[ir][j][k];
                                  if (!(j == 0 || k == 0) && ztbc == 0)
                                      d->ez->data[i][j][k] += dz * cbz * bnd->cpml.pezx_xn->data[ir][j][k];
                              }

                              if (sb->byn == SV_BOUNDARY_CPML && (j >= (set->sp.yres - sb->depth_byn + 1)) && j < set->sp.yres) {
                                  jr = j - (set->sp.yres - sb->depth_byn);
                                  cbx = (dt / mats[matindex].epsilon / dx) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  cbz = (dt / mats[matindex].epsilon / dz) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  be = bnd->cpml.be_yn[jr];
                                  ce = bnd->cpml.ce_yn[jr];

                                  bnd->cpml.pexy_yn->data[i][jr][k] = be*bnd->cpml.pexy_yn->data[i][jr][k]
                                      + ce * ((d->hz->data[i][j][k] - d->hz->data[i][j - 1][k])) / dx;
                                  bnd->cpml.pezy_yn->data[i][jr][k] = be*bnd->cpml.pezy_yn->data[i][jr][k]
                                      + ce * ((d->hx->data[i][j][k] - d->hx->data[i][j - 1][k])) / dz;
  
                                  if (!(i == 0 || k == 0) && xtbc == 0)
                                    d->ex->data[i][j][k] += dx * cbx * bnd->cpml.pexy_yn->data[i][jr][k];
                                  if (!(i == 0 || k == 0) && ztbc == 0)
                                    d->ez->data[i][j][k] -= dz * cbz * bnd->cpml.pezy_yn->data[i][jr][k];
                              }

                              if (sb->bzn == SV_BOUNDARY_CPML && (k >= (set->sp.zres - sb->depth_bzn + 1)) && k < set->sp.zres) {
                                  kr = k - (set->sp.zres - sb->depth_bzn);
                                  cbx = (dt / mats[matindex].epsilon / dx) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  cby = (dt / mats[matindex].epsilon / dy) / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
                                  be = bnd->cpml.be_zn[kr];
                                  ce = bnd->cpml.ce_zn[kr];

                                  bnd->cpml.peyz_zn->data[i][j][kr] = be*bnd->cpml.peyz_zn->data[i][j][kr]
                                      + ce * ((d->hx->data[i][j][k] - d->hx->data[i][j][k - 1])) / dy;
                                  bnd->cpml.pexz_zn->data[i][j][kr] = be*bnd->cpml.pexz_zn->data[i][j][kr]
                                      + ce * ((d->hy->data[i][j][k] - d->hy->data[i][j][k - 1])) / dx;

                                  if (!(i == 0 || j == 0) && ytbc == 0)
                                      d->ey->data[i][j][k] += dy * cby * bnd->cpml.peyz_zn->data[i][j][kr];
                                  if (!(i == 0 || j == 0) && xtbc == 0)
                                      d->ex->data[i][j][k] -= dx * cbx * bnd->cpml.pexz_zn->data[i][j][kr];
                              }
                            }
                            /*end of cpml update*/

                            /*store previous p terms*/
                            vdpx = d->dpx[0]->data[i][j][k];
                            vdpy = d->dpy[0]->data[i][j][k];
                            vdpz = d->dpz[0]->data[i][j][k];
                            for (n = 0; n < 2; n++) {
                                vpx[n] = d->px[n]->data[i][j][k];
                                vpy[n] = d->py[n]->data[i][j][k];
                                vpz[n] = d->pz[n]->data[i][j][k];
                            }

                            /*update all p terms*/
                            d->dpx[0]->data[i][j][k] = mats[matindex].ade_a0 * d->dpxp[0]->data[i][j][k] + mats[matindex].ade_a1 * d->dpx[0]->data[i][j][k]
                                + mats[matindex].ade_a2 * (d->exp->data[i][j][k] + 2 * vexp + d->ex->data[i][j][k]);
                            d->dpy[0]->data[i][j][k] = mats[matindex].ade_a0 * d->dpyp[0]->data[i][j][k] + mats[matindex].ade_a1 * d->dpy[0]->data[i][j][k]
                                + mats[matindex].ade_a2 * (d->eyp->data[i][j][k] + 2 * veyp + d->ey->data[i][j][k]);
                            d->dpz[0]->data[i][j][k] = mats[matindex].ade_a0 * d->dpzp[0]->data[i][j][k] + mats[matindex].ade_a1 * d->dpz[0]->data[i][j][k]
                                + mats[matindex].ade_a2 * (d->ezp->data[i][j][k] + 2 * vezp + d->ez->data[i][j][k]);

                            //d->dpx[0]->data[i][j][k] = d->dpy[0]->data[i][j][k] = d->dpz[0]->data[i][j][k] = 0;
                            for (n = 0; n < 2; n++) {
                                d->px[n]->data[i][j][k] = mats[matindex].ade_bp0[n] * d->pxp[n]->data[i][j][k] + mats[matindex].ade_bp1[n] * d->px[n]->data[i][j][k] +
                                    mats[matindex].ade_bp2[n] * d->exp->data[i][j][k] + mats[matindex].ade_bp3[n] * vexp +
                                    mats[matindex].ade_bp4[n] * d->ex->data[i][j][k];

                                d->py[n]->data[i][j][k] = mats[matindex].ade_bp0[n] * d->pyp[n]->data[i][j][k] + mats[matindex].ade_bp1[n] * d->py[n]->data[i][j][k] +
                                    mats[matindex].ade_bp2[n] * d->eyp->data[i][j][k] + mats[matindex].ade_bp3[n] * veyp +
                                    mats[matindex].ade_bp4[n] * d->ey->data[i][j][k];

                                d->pz[n]->data[i][j][k] = mats[matindex].ade_bp0[n] * d->pzp[n]->data[i][j][k] + mats[matindex].ade_bp1[n] * d->pz[n]->data[i][j][k] +
                                    mats[matindex].ade_bp2[n] * d->ezp->data[i][j][k] + mats[matindex].ade_bp3[n] * vezp +
                                    mats[matindex].ade_bp4[n] * d->ez->data[i][j][k];

                                // d->px[n]->data[i][j][k] = d->py[n]->data[i][j][k] = d->pz[n]->data[i][j][k] = 0;  
                            }

                            /*save e n-1*/
                            d->exp->data[i][j][k] = vexp;
                            d->eyp->data[i][j][k] = veyp;
                            d->ezp->data[i][j][k] = vezp;

                            /*save previous p terms*/
                            d->dpxp[0]->data[i][j][k] = vdpx;
                            d->dpyp[0]->data[i][j][k] = vdpy;
                            d->dpzp[0]->data[i][j][k] = vdpz;
                            for (n = 0; n < 2; n++) {
                                d->pxp[n]->data[i][j][k] = vpx[n];
                                d->pyp[n]->data[i][j][k] = vpy[n];
                                d->pzp[n]->data[i][j][k] = vpz[n];
                            }
                        } // else if SV_MAT_ADE
                    } //if tbc
                } //if (d->nmat > 0 && (d->plrcx || d->dplrcx || d->pxp))


                /*normal state or compute what rested from before*/
                if ((xtbc + ytbc + ztbc) != 0) {

                    sv_yee_data_get_cabs(d, set, mats, nmat, dt, i, j, k, &cax, &cay, &caz, &cbx, &cby, &cbz);

                    if (j > 0 && k > 0 && xtbc)
                        d->ex->data[i][j][k] = cax*d->ex->data[i][j][k] + cbx*((hz - d->hz->data[i][j - 1][k]) / kappay -
                        (hy - d->hy->data[i][j][k - 1]) / kappaz);

                    if (i > 0 && k > 0 && ytbc)
                        d->ey->data[i][j][k] = cay*d->ey->data[i][j][k] + cby*((hx - d->hx->data[i][j][k - 1]) / kappaz -
                        (hz - d->hz->data[i - 1][j][k]) / kappax);

                    if (i > 0 && j > 0 && ztbc)
                        d->ez->data[i][j][k] = caz*d->ez->data[i][j][k] + cbz*((hy - d->hy->data[i - 1][j][k]) / kappax -
                        (hx - d->hx->data[i][j - 1][k]) / kappay);
                        


                    


                    //                  printf("%g %g %g %g %g %g\n", cax, cay, caz, cbx, cby, cbz);
                    //                    if (i==58 && j==22 && k==50) printf("%d    %.10f %.10f %.10f    %.10f %.10f %.10f:  %g %g %g   ds %g %g %g\n", set->sc.step_act, cax, cay, caz, cbx, cby, cbz, d->ex->data[i][j][k], d->ey->data[i][j][k], d->ez->data[i][j][k], d->dx, d->dy, d->dz);
                }

                /*eventually apply also PEC and scattered field source on it*/
                if (matindex >= 0 && i > 0 && j > 0 && k > 0 &&
                    (mats[matindex].type == SV_MAT_PEC || mats[d->mat->data[i - 1][j][k]].type == SV_MAT_PEC ||
                     mats[d->mat->data[i][j - 1][k]].type == SV_MAT_PEC || mats[d->mat->data[i][j][k - 1]].type == SV_MAT_PEC ||
                     mats[d->mat->data[i - 1][j - 1][k]].type == SV_MAT_PEC || mats[d->mat->data[i][j - 1][k - 1]].type == SV_MAT_PEC ||
                     mats[d->mat->data[i - 1][j][k - 1]].type == SV_MAT_PEC || mats[d->mat->data[i - 1][j - 1][k - 1]].type == SV_MAT_PEC)) {
                    if ((mats[d->mat->data[i - 1][j][k]].type != SV_MAT_PEC && mats[d->mat->data[i][j][k]].type == SV_MAT_PEC) ||
                        (mats[d->mat->data[i - 1][j][k]].type == SV_MAT_PEC && mats[d->mat->data[i][j][k]].type != SV_MAT_PEC)) {
                        if (src && src->sf) {
                            dd = dcomp(i, j, k,
                                      d->xres, d->yres, d->zres,
                                      src->sf->ia_theta, src->sf->ia_phi,
                                      0, d->xres, 0, d->yres, 0, d->zres);
                            ve = sv_dline_get_dval(src->sf->jp->e, dd);

                            d->ey->data[i][j][k] = -gey(ve, src->sf->ia_theta, src->sf->ia_phi, src->sf->ia_psi);
                            d->ez->data[i][j][k] = -gez(ve, src->sf->ia_theta, src->sf->ia_phi, src->sf->ia_psi);
                        } else {
                            d->ey->data[i][j][k] = 0;
                            d->ez->data[i][j][k] = 0;
                        }
                    }

                    if ((mats[d->mat->data[i][j - 1][k]].type != SV_MAT_PEC && mats[d->mat->data[i][j][k]].type == SV_MAT_PEC) ||
                        (mats[d->mat->data[i][j - 1][k]].type == SV_MAT_PEC && mats[d->mat->data[i][j][k]].type != SV_MAT_PEC)) {
                        if (src && src->sf) {
                            dd = dcomp(i, j, k,
                                      d->xres, d->yres, d->zres,
                                      src->sf->ia_theta, src->sf->ia_phi,
                                      0, d->xres, 0, d->yres, 0, d->zres);
                            ve = sv_dline_get_dval(src->sf->jp->e, dd);

                            d->ex->data[i][j][k] = -gex(ve, src->sf->ia_theta, src->sf->ia_phi, src->sf->ia_psi);
                            d->ez->data[i][j][k] = -gez(ve, src->sf->ia_theta, src->sf->ia_phi, src->sf->ia_psi);
                        } else {
                            d->ex->data[i][j][k] = 0;
                            d->ez->data[i][j][k] = 0;
                        }
                    }

                    if ((mats[d->mat->data[i][j][k - 1]].type != SV_MAT_PEC && mats[d->mat->data[i][j][k]].type == SV_MAT_PEC) ||
                        (mats[d->mat->data[i][j][k - 1]].type == SV_MAT_PEC && mats[d->mat->data[i][j][k]].type != SV_MAT_PEC)) {
                        if (src && src->sf) {
                            dd = dcomp(i, j, k,
                                      d->xres, d->yres, d->zres,
                                      src->sf->ia_theta, src->sf->ia_phi,
                                      0, d->xres, 0, d->yres, 0, d->zres);
                            ve = sv_dline_get_dval(src->sf->jp->e, dd);

                            d->ex->data[i][j][k] = -gex(ve, src->sf->ia_theta, src->sf->ia_phi, src->sf->ia_psi);
                            d->ey->data[i][j][k] = -gey(ve, src->sf->ia_theta, src->sf->ia_phi, src->sf->ia_psi);
                        } else {
                            d->ex->data[i][j][k] = 0;
                            d->ey->data[i][j][k] = 0;
                        }
                    }

                    if ((mats[d->mat->data[i - 1][j - 1][k]].type == SV_MAT_PEC && mats[d->mat->data[i][j][k]].type != SV_MAT_PEC)) {
                        if (src && src->sf) {
                            dd = dcomp(i, j, k,
                                      d->xres, d->yres, d->zres,
                                      src->sf->ia_theta, src->sf->ia_phi,
                                      0, d->xres, 0, d->yres, 0, d->zres);
                            ve = sv_dline_get_dval(src->sf->jp->e, dd);

                            d->ez->data[i][j][k] = -gez(ve, src->sf->ia_theta, src->sf->ia_phi, src->sf->ia_psi);
                        } else {
                            d->ez->data[i][j][k] = 0;
                        }
                    }

                    if ((mats[d->mat->data[i - 1][j][k - 1]].type == SV_MAT_PEC && mats[d->mat->data[i][j][k]].type != SV_MAT_PEC)) {
                        if (src && src->sf) {
                            dd = dcomp(i, j, k,
                                      d->xres, d->yres, d->zres,
                                      src->sf->ia_theta, src->sf->ia_phi,
                                      0, d->xres, 0, d->yres, 0, d->zres);
                            ve = sv_dline_get_dval(src->sf->jp->e, dd);

                            d->ey->data[i][j][k] = -gey(ve, src->sf->ia_theta, src->sf->ia_phi, src->sf->ia_psi);
                        } else {
                            d->ey->data[i][j][k] = 0;
                        }
                    }

                    if ((mats[d->mat->data[i][j - 1][k - 1]].type == SV_MAT_PEC && mats[d->mat->data[i][j][k]].type != SV_MAT_PEC)) {
                        if (src && src->sf) {
                            dd = dcomp(i, j, k,
                                      d->xres, d->yres, d->zres,
                                      src->sf->ia_theta, src->sf->ia_phi,
                                      0, d->xres, 0, d->yres, 0, d->zres);
                            ve = sv_dline_get_dval(src->sf->jp->e, dd);

                            d->ex->data[i][j][k] = -gex(ve, src->sf->ia_theta, src->sf->ia_phi, src->sf->ia_psi);
                        } else {
                            d->ex->data[i][j][k] = 0;
                        }
                    }
                } //pec

            } // k
        } // j
    } // i

    if (set->sc.verbose > 1)
        printf("done.\n");

    return 0;
}

gboolean
sv_yee_data_is_pec(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k)
{
    gint matindex;

    if (nmat > 0)
        matindex = d->mat->data[i][j][k];
    else
        matindex = -1;

    /*there is a local tabulated material in this point: use it always*/
    if (matindex > 0 && (mats[matindex].type == SV_MAT_PEC)) return 1;

    return 0;
}

gdouble
sv_yee_data_get_epsilon(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k)
{
    gint matindex;

    if (nmat > 0)
        matindex = d->mat->data[i][j][k];
    else
        matindex = -1;

    /*there is a local tabulated material in this point: use it always*/
    if (matindex > 0
        && (mats[matindex].type == SV_MAT_LINTAB
            || mats[matindex].type == SV_MAT_DRUDE            
            || mats[matindex].type == SV_MAT_PLRC
            || mats[matindex].type == SV_MAT_ADE
            || mats[matindex].type == SV_MAT_CP
            || mats[matindex].type == SV_MAT_CP3))
        return mats[matindex].epsilon;
    /*there is no local material, but there are some voxel-by-voxel data*/
    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC) {
        return d->epsilon->data[i][j][k];
    } else
        return EPSILON_0;
}

gdouble sv_yee_data_get_sigma(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k)
{
    gint matindex;

    if (nmat > 0)
        matindex = d->mat->data[i][j][k];
    else
        matindex = -1;

    /*there is a local tabulated material in this point: use it always*/
    if (matindex > 0
        && (mats[matindex].type == SV_MAT_LINTAB
            || mats[matindex].type == SV_MAT_DRUDE
            || mats[matindex].type == SV_MAT_PLRC
            || mats[matindex].type == SV_MAT_ADE
            || mats[matindex].type == SV_MAT_CP
            || mats[matindex].type == SV_MAT_CP3))
        return mats[matindex].sigma;
    /*there is no local material, but there are some voxel-by-voxel data*/
    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC)
        return d->sigma->data[i][j][k];
    else
        return 0;
}


gint
sv_yee_data_get_epsilon_sigma(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k, gdouble *epsval, gdouble *sigmaval)
{
    gint matindex;
    gdouble epsilon = EPSILON_0;
    gdouble sigma = 0;

    if (nmat > 0)
        matindex = d->mat->data[i][j][k];
    else
        matindex = -1;

    /*there is a local tabulated material in this point: use it always*/
    if (matindex > 0
        && (mats[matindex].type == SV_MAT_LINTAB
            || mats[matindex].type == SV_MAT_DRUDE
            || mats[matindex].type == SV_MAT_PLRC
            || mats[matindex].type == SV_MAT_ADE
            || mats[matindex].type == SV_MAT_CP
            || mats[matindex].type == SV_MAT_CP3)) {
        epsilon = mats[matindex].epsilon;
        sigma = mats[matindex].sigma;
    }
    /*there is no local material, but there are some voxel-by-voxel data*/
    else if (d->epsilon != NULL && d->sigma != NULL) {
        epsilon = d->epsilon->data[i][j][k];
        sigma = d->sigma->data[i][j][k];
    }

    *epsval = epsilon;
    *sigmaval = sigma;

    return matindex;
}

gdouble
sv_yee_data_get_cb(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gdouble dt, gint i, gint j, gint k)
{
    gint matindex;
    gdouble sumxi0, sumchi0, vch0, epsilon;

    if (nmat > 0)
        matindex = d->mat->data[i][j][k];
    else
        matindex = -1;

    /*there is a local tabulated material in this point: use it always*/
    if (matindex > 0
        && (mats[matindex].type == SV_MAT_LINTAB
            || mats[matindex].type == SV_MAT_CP
            || mats[matindex].type == SV_MAT_CP3)) {
        return (dt / mats[matindex].epsilon / d->dx) / (1 + mats[matindex].sigma * dt / 2 / mats[matindex].epsilon);
    }
    if (matindex > 0 && (mats[matindex].type == SV_MAT_DRUDE)) {
        vch0 = ch0(mats[matindex].drude_omega_p, dt);
        epsilon = mats[matindex].epsilon;

        return dt/(epsilon+vch0)/EPSILON_0/d->dx;
    }

    if (matindex > 0 && (mats[matindex].type == SV_MAT_ADE)) {
        return mats[matindex].ade_c0 / d->dx;
    }

    if (matindex > 0 && (mats[matindex].type == SV_MAT_PLRC)) {
        sumxi0 = mats[matindex].plrc_d_xi + mats[matindex].plrc_p_xi[0] + mats[matindex].plrc_p_xi[1];
        sumchi0 = mats[matindex].plrc_d_chi + mats[matindex].plrc_p_chi[0] + mats[matindex].plrc_p_chi[1];

        return 2 * dt / (2 * EPSILON_0 * (mats[matindex].epsilon / EPSILON_0 - sumxi0 + sumchi0) + mats[matindex].sigma * dt) / d->dx;
    }

    /*there is no local material, but there are some voxel-by-voxel data*/
    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC) {
        return (dt / d->epsilon->data[i][j][k] / d->dx) / (1 + d->sigma->data[i][j][k] * dt / 2 / d->epsilon->data[i][j][k]);
    }

    return dt / EPSILON_0 / d->dx;
}

void
sv_yee_data_get_cabs(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gdouble dt, gint i, gint j, gint k, gdouble *cax, gdouble *cay, gdouble *caz, gdouble *cbx, gdouble *cby, gdouble *cbz)
{
    gdouble epsilonx, epsilony, epsilonz, epsilon;
    gdouble sigmax, sigmay, sigmaz, sigma;
    gint nx, ny, nz;

    sv_yee_data_get_epsilon_sigma(d, set, mats, nmat, i, j, k, &epsilon, &sigma);
    epsilonx = epsilony = epsilonz = epsilon;
    sigmax = sigmay = sigmaz = sigma;

    /*do spatial interpolation only if we are close to materials boundary*/

    
    if (d->bnds->data[i][j][k] > 0) {
        nx = ny = nz = 1;

        if (i > 0/* && !(nmat>0 && (d->mat->data[i-1][j][k]==SV_MAT_PLRC || d->mat->data[i-1][j][k]==SV_MAT_ADE))*/) {  //do not include ade/plrc epsilon or sigma unless you know it from dispersion formula
            sv_yee_data_get_epsilon_sigma(d, set, mats, nmat, i - 1, j, k, &epsilon, &sigma);
            epsilony += epsilon;
            sigmay += sigma;
            epsilonz += epsilon;
            sigmaz += sigma;
            ny++;
            nz++;
        }


        if (j > 0/* && !(nmat>0 && (d->mat->data[i][j-1][k]==SV_MAT_PLRC || d->mat->data[i][j-1][k]==SV_MAT_ADE))*/) { //do not include ade/plrc epsilon or sigma unless you know it from dispersion formula
            sv_yee_data_get_epsilon_sigma(d, set, mats, nmat, i, j - 1, k, &epsilon, &sigma);
            epsilonx += epsilon;
            sigmax += sigma;
            epsilonz += epsilon;
            sigmaz += sigma;
            nx++;
            nz++;
        }

        if (k > 0/* && !(nmat>0 && (d->mat->data[i][j][k-1]==SV_MAT_PLRC || d->mat->data[i][j][k-1]==SV_MAT_ADE))*/) { //do not include ade/plrc epsilon or sigma unless you know it from dispersion formula
            sv_yee_data_get_epsilon_sigma(d, set, mats, nmat, i, j, k - 1, &epsilon, &sigma);
            epsilonx += epsilon;
            sigmax += sigma;
            epsilony += epsilon;
            sigmay += sigma;
            nx++;
            ny++;
        }

        if (j > 0 && k > 0 /*&&
                !(nmat>0 && (d->mat->data[i][j-1][k]==SV_MAT_PLRC || d->mat->data[i][j-1][k]==SV_MAT_ADE || d->mat->data[i][j][k-1]==SV_MAT_PLRC || d->mat->data[i][j][k-1]==SV_MAT_ADE))*/) {
            sv_yee_data_get_epsilon_sigma(d, set, mats, nmat, i, j - 1, k - 1, &epsilon, &sigma);
            epsilonx += epsilon;
            sigmax += sigma;
            nx++;
        }
        if (i > 0 && k > 0 /*&&
                    !(nmat>0 && (d->mat->data[i-1][j][k]==SV_MAT_PLRC || d->mat->data[i-1][j][k]==SV_MAT_ADE || d->mat->data[i][j][k-1]==SV_MAT_PLRC || d->mat->data[i][j][k-1]==SV_MAT_ADE))*/) {
            sv_yee_data_get_epsilon_sigma(d, set, mats, nmat, i - 1, j, k - 1, &epsilon, &sigma);
            epsilony += epsilon;
            sigmay += sigma;
            ny++;
        }
        if (i > 0 && j > 0 /*&&
                !(nmat>0 && (d->mat->data[i-1][j][k]==SV_MAT_PLRC || d->mat->data[i-1][j][k]==SV_MAT_ADE || d->mat->data[i][j-1][k]==SV_MAT_PLRC || d->mat->data[i][j-1][k]==SV_MAT_ADE))*/) {
            sv_yee_data_get_epsilon_sigma(d, set, mats, nmat, i - 1, j - 1, k, &epsilon, &sigma);
            epsilonz += epsilon;
            sigmaz += sigma;
            nz++;
        }


        epsilonx /= nx;
        sigmax /= nx;
        epsilony /= ny;
        sigmay /= ny;
        epsilonz /= nz;
        sigmaz /= nz;
    }

    *cax = (1 - sigmax*dt / 2 / epsilonx)
        / (1 + sigmax*dt / 2 / epsilonx);
    *cbx = (dt / epsilonx / d->dx) / (1 + sigmax*dt / 2 / epsilonx);


    if (d->bnds->data[i][j][k] > 0 || (d->dx != d->dy || d->dx != d->dz)) {
        *cay = (1 - sigmay*dt / 2 / epsilony)
            / (1 + sigmay*dt / 2 / epsilony);
        *cby = (dt / epsilony / d->dy) / (1 + sigmay*dt / 2 / epsilony);

        *caz = (1 - sigmaz*dt / 2 / epsilonz)
            / (1 + sigmaz*dt / 2 / epsilonz);
        *cbz = (dt / epsilonz / d->dz) / (1 + sigmaz*dt / 2 / epsilonz);
    } else {
        *cay = *caz = *cax;
        *cby = *cbz = *cbx;
    }
}


gint
sv_yee_data_get_ca_cb(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gdouble dt, gint i, gint j, gint k, gdouble *caval, gdouble *cbval)
{
    gint matindex;
    gdouble ca = 1;
    gdouble cb = dt / EPSILON_0 / d->dx;

    if (nmat > 0)
        matindex = d->mat->data[i][j][k];
    else
        matindex = -1;

    /*there is a local tabulated material in this point: use it always*/
    if (matindex > 0
        && (mats[matindex].type == SV_MAT_LINTAB
            || mats[matindex].type == SV_MAT_DRUDE
            || mats[matindex].type == SV_MAT_PLRC
            || mats[matindex].type == SV_MAT_ADE
            || mats[matindex].type == SV_MAT_CP
            || mats[matindex].type == SV_MAT_CP3)) {
        ca = (1 - mats[matindex].sigma*dt / 2 / mats[matindex].epsilon)
            / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
        cb = (dt / mats[matindex].epsilon / d->dx)
            / (1 + mats[matindex].sigma*dt / 2 / mats[matindex].epsilon);
    }
    /*there is no local material, but there are some voxel-by-voxel data*/
    else if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC) {
        ca = (1 - d->sigma->data[i][j][k] * dt / 2 / d->epsilon->data[i][j][k])
            / (1 + d->sigma->data[i][j][k] * dt / 2 / d->epsilon->data[i][j][k]);
        cb = (dt / d->epsilon->data[i][j][k] / d->dx) / (1 + d->sigma->data[i][j][k] * dt / 2 / d->epsilon->data[i][j][k]);
    }

    *caval = ca;
    *cbval = cb;

    return matindex;
}


gdouble
sv_yee_data_get_mu(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k)
{
    gint matindex;

    if (nmat > 0)
        matindex = d->mat->data[i][j][k];
    else
        matindex = -1;

    /*there is a local tabulated material in this point: use it always*/
    if (matindex > 0
        && (mats[matindex].type == SV_MAT_LINTAB
            || mats[matindex].type == SV_MAT_DRUDE
            || mats[matindex].type == SV_MAT_PLRC
            || mats[matindex].type == SV_MAT_ADE
            || mats[matindex].type == SV_MAT_CP
            || mats[matindex].type == SV_MAT_CP3))
        return mats[matindex].mu;
    /*there is no local material, but there are some voxel-by-voxel data*/
    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_MAGNETIC)
        return d->mu->data[i][j][k];
    else
        return MU_0;
}

gdouble sv_yee_data_get_sigast(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k)
{
    gint matindex;

    if (nmat > 0)
        matindex = d->mat->data[i][j][k];
    else
        matindex = -1;

    /*there is a local tabulated material in this point: use it always*/
    if (matindex > 0
        && (mats[matindex].type == SV_MAT_LINTAB
            || mats[matindex].type == SV_MAT_DRUDE
            || mats[matindex].type == SV_MAT_PLRC
            || mats[matindex].type == SV_MAT_ADE
            || mats[matindex].type == SV_MAT_CP
            || mats[matindex].type == SV_MAT_CP3))
        return mats[matindex].sigast;
    /*there is no local material, but there are some voxel-by-voxel data*/
    else if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_MAGNETIC)
        return d->sigast->data[i][j][k];
    else
        return 0;
}


gint
sv_yee_data_get_mu_sigast(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k, gdouble *muval, gdouble *sigastval)
{
    gint matindex;
    gdouble mu = MU_0;
    gdouble sigast = 0;

    if (nmat > 0)
        matindex = d->mat->data[i][j][k];
    else
        matindex = -1;

    /*there is a local tabulated material in this point: use it always*/
    if (matindex > 0
        && (mats[matindex].type == SV_MAT_LINTAB
            || mats[matindex].type == SV_MAT_DRUDE
            || mats[matindex].type == SV_MAT_PLRC
            || mats[matindex].type == SV_MAT_ADE
            || mats[matindex].type == SV_MAT_CP
            || mats[matindex].type == SV_MAT_CP3)) {
        mu = mats[matindex].mu;
        sigast = mats[matindex].sigast;
    }
    /*there is no local material, but there are some voxel-by-voxel data*/
    else if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_MAGNETIC) {
        mu = d->mu->data[i][j][k];
        sigast = d->sigast->data[i][j][k];
    }

    *muval = mu;
    *sigastval = sigast;

    return matindex;
}
gdouble sv_yee_data_get_db(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gdouble dt, gint i, gint j, gint k)
{
    gint matindex;

    if (nmat > 0)
        matindex = d->mat->data[i][j][k];
    else
        matindex = -1;

    /*there is a local tabulated material in this point: use it always*/
    if (matindex > 0
        && (mats[matindex].type == SV_MAT_LINTAB
            || mats[matindex].type == SV_MAT_DRUDE
            || mats[matindex].type == SV_MAT_PLRC
            || mats[matindex].type == SV_MAT_ADE
            || mats[matindex].type == SV_MAT_CP
            || mats[matindex].type == SV_MAT_CP3)) {
        return (dt / mats[matindex].mu / d->dx) / (1 + mats[matindex].sigast*dt / 2 / mats[matindex].mu);
    }


    /*there is no local material, but there are some voxel-by-voxel data*/

    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_MAGNETIC)
        return (dt / d->mu->data[i][j][k] / d->dx) / (1 + d->sigast->data[i][j][k] * dt / 2 / d->mu->data[i][j][k]);
    else
        return dt / MU_0 / d->dx;
}


gint
sv_yee_data_get_da_db(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gdouble dt, gint i, gint j, gint k, gdouble *daval, gdouble *dbval)
{
    gint matindex;
    gdouble da = 1;
    gdouble db = dt / MU_0 / d->dx;

    if (nmat > 0)
        matindex = d->mat->data[i][j][k];
    else
        matindex = -1;

    /*there is a local tabulated material in this point: use it always*/
    if (matindex > 0 && (mats[matindex].type == SV_MAT_LINTAB
                         || mats[matindex].type == SV_MAT_DRUDE
                         || mats[matindex].type == SV_MAT_PLRC
                         || mats[matindex].type == SV_MAT_ADE
                         || mats[matindex].type == SV_MAT_CP
                         || mats[matindex].type == SV_MAT_CP3)) {
        da = (1 - mats[matindex].sigast*dt / 2 / mats[matindex].mu) / (1 + mats[matindex].sigast*dt / 2 / mats[matindex].mu);
        db = (dt / mats[matindex].mu / d->dx) / (1 + mats[matindex].sigast*dt / 2 / mats[matindex].mu);
    }
    /*there is no local material, but there are some voxel-by-voxel data*/
    else if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_MAGNETIC) {
        da = (1 - d->sigast->data[i][j][k] * dt / 2 / d->mu->data[i][j][k]) / (1 + d->sigast->data[i][j][k] * dt / 2 / d->mu->data[i][j][k]);
        db = (dt / d->mu->data[i][j][k] / d->dx) / (1 + d->sigast->data[i][j][k] * dt / 2 / d->mu->data[i][j][k]);
    }

    *daval = da;
    *dbval = db;

    return matindex;
}

/*get mat type information in form relevant to TSF source treatment - now accepting only vacuum*/
SvMatType
sv_yee_data_get_tsf_mattype(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k)
{
    gint matindex;

    if (nmat > 0)
        matindex = d->mat->data[i][j][k];
    else
        matindex = -1;

    if (matindex == 0 || matindex == -1) {
        //printf("%d %g %g\n", matindex, d->epsilon->data[i][j][k], EPSILON_0);
        if ((set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC) && fabs(1.0 - d->epsilon->data[i][j][k] / EPSILON_0) < 0.001)
            return 0;
        if (set->plan.matmode == SV_MATMODE_NONE)
            return 0;
    }

    if (matindex > 0 && mats[matindex].type == SV_MAT_LINTAB && fabs(1.0 - mats[matindex].epsilon / EPSILON_0) < 0.001)
        return 0;

    return 1;
}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
