
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


/*  pfarfield.c : 
 *  periodic near-to-far field transformation
 */

#define MAX_STORAGE 1

#include "pfarfield.h"
#include "pool.h"
#include "settings.h"
#include "constants.h"
#include <math.h>
#include <omp.h>

SvPFarfield*
sv_pfarfield_new(SvPool *mp, SvSet *set)
{
    gint n;
    gdouble div;
    SvPFarfield* ff = (SvPFarfield*)g_malloc(sizeof(SvPFarfield));

    ff->ndata = 1000 + 2*set->sc.nsteps + (gint)(sqrt(set->sp.xres*set->sp.xres + set->sp.yres*set->sp.yres + set->sp.zres*set->sp.zres) / (LIGHT_SPEED*set->plan.dt/set->sp.dx));
    ff->prpoints = (SvPRPoint *) g_malloc(set->spf.nrs*sizeof(SvPRPoint));

    if (set->sc.verbose > 1) {
        printf("Initializing periodic far field positions data (%d points)... ", set->spf.nrs);
        fflush(stdout);
    }

    div = (LIGHT_SPEED * set->plan.dt/set->sp.dx);
    for (n = 0; n < set->spf.nrs; n++) {
        ff->prpoints[n].ex = gwy_data_line_new(ff->ndata, ff->ndata, TRUE);
        ff->prpoints[n].ey = gwy_data_line_new(ff->ndata, ff->ndata, TRUE);
        ff->prpoints[n].ez = gwy_data_line_new(ff->ndata, ff->ndata, TRUE);
        ff->prpoints[n].istart = (gint)MAX(0, sqrt((set->spf.ri[n]/div)*(set->spf.ri[n]/div) + (set->spf.rj[n]/div)*(set->spf.rj[n]/div) + (set->spf.rk[n]/div)*(set->spf.rk[n]/div)) - 500);
    }

    if (set->spf.postprocess) {
	if (set->sc.verbose>1) {
	    printf("Allocating postprocess fields (%dx%dx%d)...", set->spf.box_in-set->spf.box_i0, set->spf.box_jn-set->spf.box_j0, set->sc.nsteps - set->spf.ppstart);
	    fflush(stdout);
	}
 
	ff->ex_k0 =  sv_dcube_new(set->spf.box_in-set->spf.box_i0, set->spf.box_jn-set->spf.box_j0, set->sc.nsteps - set->spf.ppstart, 1, 1, 1, 1);
	ff->dex_k0 = sv_dcube_new(set->spf.box_in-set->spf.box_i0, set->spf.box_jn-set->spf.box_j0, set->sc.nsteps - set->spf.ppstart, 1, 1, 1, 1);

	ff->ey_k0 =  sv_dcube_new(set->spf.box_in-set->spf.box_i0, set->spf.box_jn-set->spf.box_j0, set->sc.nsteps - set->spf.ppstart, 1, 1, 1, 1);
	ff->dey_k0 = sv_dcube_new(set->spf.box_in-set->spf.box_i0, set->spf.box_jn-set->spf.box_j0, set->sc.nsteps - set->spf.ppstart, 1, 1, 1, 1);

	ff->ez_k0 =  sv_dcube_new(set->spf.box_in-set->spf.box_i0, set->spf.box_jn-set->spf.box_j0, set->sc.nsteps - set->spf.ppstart, 1, 1, 1, 1);
	ff->dez_k0 = sv_dcube_new(set->spf.box_in-set->spf.box_i0, set->spf.box_jn-set->spf.box_j0, set->sc.nsteps - set->spf.ppstart, 1, 1, 1, 1);
 
	ff->ex_kn =  sv_dcube_new(set->spf.box_in-set->spf.box_i0, set->spf.box_jn-set->spf.box_j0, set->sc.nsteps - set->spf.ppstart, 1, 1, 1, 1);
	ff->dex_kn = sv_dcube_new(set->spf.box_in-set->spf.box_i0, set->spf.box_jn-set->spf.box_j0, set->sc.nsteps - set->spf.ppstart, 1, 1, 1, 1);

	ff->ey_kn =  sv_dcube_new(set->spf.box_in-set->spf.box_i0, set->spf.box_jn-set->spf.box_j0, set->sc.nsteps - set->spf.ppstart, 1, 1, 1, 1);
	ff->dey_kn = sv_dcube_new(set->spf.box_in-set->spf.box_i0, set->spf.box_jn-set->spf.box_j0, set->sc.nsteps - set->spf.ppstart, 1, 1, 1, 1);

	ff->ez_kn =  sv_dcube_new(set->spf.box_in-set->spf.box_i0, set->spf.box_jn-set->spf.box_j0, set->sc.nsteps - set->spf.ppstart, 1, 1, 1, 1);
	ff->dez_kn = sv_dcube_new(set->spf.box_in-set->spf.box_i0, set->spf.box_jn-set->spf.box_j0, set->sc.nsteps - set->spf.ppstart, 1, 1, 1, 1);
    }

    if (set->sc.verbose > 1) printf("done.\n");

    return ff;
}

/*
  int
  sv_pool_farfield_allocate_storage(SvPool *mp, SvSet *set)
  {
  gint ix, iy, iz, i, j, k, n, x, y, z;
  gint i0, i1, j0, j1, k0, k1;
  SvFarfield* ff = mp->farfield;

  i0 = set->sf.i0;
  i1 = set->sf.i1;
  j0 = set->sf.j0;
  j1 = set->sf.j1;
  k0 = set->sf.k0;
  k1 = set->sf.k1;

  ix = i1 - i0;
  iy = j1 - j0;
  iz = k1 - k0;

  for (n=0; n<set->sf.nrs; n++) {
  ff->rpoints[n].r_front       = gwy_data_field_new(ix, iy, ix, iy, TRUE);
  ff->rpoints[n].ctheta_front  = gwy_data_field_new(ix, iy, ix, iy, TRUE);

  x = set->sf.ri[n];
  y = set->sf.rj[n];
  z = set->sf.rk[n];

  }
  }
*/

int
sv_pool_pfarfield(SvPool *mp, SvSet *set)
{
    gint n, i, j, tit;
    gint i0, i1, j0, j1, k0, k1, ired = 0, jred = 0;
    gdouble dx = set->sp.dx;
    gdouble dy = set->sp.dy;
    gdouble dt = set->plan.dt;
    gdouble size = 0, r = 0, a = 0, b = 0, c = 0, ctheta0 = 0, nred = 0;
    gint nn = 0;
    gint pi = 0, pj = 0;
    gint pri = 0, prj = 0;
    gdouble mult = 0;
    gdouble cvala = 0, cvalb = 0, aval = 0, ex = 0, ey = 0, ez = 0;

    if (set->spf.nrs == 0) 
        return 0;

    i0 = set->spf.box_i0;
    i1 = set->spf.box_in;
    j0 = set->spf.box_j0;
    j1 = set->spf.box_jn;
    k0 = set->spf.box_k0;
    k1 = set->spf.box_kn;

    if (mp->set->sc.verbose > 1) {
        printf("Running periodic farfield...  ");
        fflush(stdout);
    }

    mult = 1.0/(4.0*G_PI);
    size = dx*dy;

    if (set->spf.postprocess) {
        if (set->sc.step_act >= set->spf.ppstart) {
            if (mp->set->sc.verbose > 1) printf("storing values... "); 

            for (i = i0; i < i1; i++) {
                for (j = j0; j < j1; j++) {
                    mp->pfarfield->ex_k0->data[i-i0][j-j0][set->sc.step_act - set->spf.ppstart] = mp->d->ex->data[i][j][k0];
                    mp->pfarfield->dex_k0->data[i-i0][j-j0][set->sc.step_act - set->spf.ppstart] = mp->d->ex->data[i][j][k0+1]-mp->d->ex->data[i][j][k0-1];

                    mp->pfarfield->ey_k0->data[i-i0][j-j0][set->sc.step_act - set->spf.ppstart] =  mp->d->ey->data[i][j][k0];
                    mp->pfarfield->dey_k0->data[i-i0][j-j0][set->sc.step_act - set->spf.ppstart] = mp->d->ey->data[i][j][k0+1]-mp->d->ey->data[i][j][k0-1];

                    mp->pfarfield->ez_k0->data[i-i0][j-j0][set->sc.step_act - set->spf.ppstart] =  mp->d->ez->data[i][j][k0];
                    mp->pfarfield->dez_k0->data[i-i0][j-j0][set->sc.step_act - set->spf.ppstart] = mp->d->ez->data[i][j][k0+1]-mp->d->ez->data[i][j][k0-1];

                    mp->pfarfield->ex_kn->data[i-i0][j-j0][set->sc.step_act - set->spf.ppstart] =  mp->d->ex->data[i][j][k1];
                    mp->pfarfield->dex_kn->data[i-i0][j-j0][set->sc.step_act - set->spf.ppstart] = mp->d->ex->data[i][j][k1+1]-mp->d->ex->data[i][j][k1-1];

                    mp->pfarfield->ey_kn->data[i-i0][j-j0][set->sc.step_act - set->spf.ppstart] =  mp->d->ey->data[i][j][k1];
                    mp->pfarfield->dey_kn->data[i-i0][j-j0][set->sc.step_act - set->spf.ppstart] = mp->d->ey->data[i][j][k1+1]-mp->d->ey->data[i][j][k1-1];

                    mp->pfarfield->ez_kn->data[i-i0][j-j0][set->sc.step_act - set->spf.ppstart] =  mp->d->ez->data[i][j][k1];
                    mp->pfarfield->dez_kn->data[i-i0][j-j0][set->sc.step_act - set->spf.ppstart] = mp->d->ez->data[i][j][k1+1]-mp->d->ez->data[i][j][k1-1];
                } // j
            } // i
        } // if
        if (set->sc.step_act == (set->sc.nsteps-1)) {
            if (mp->set->sc.verbose > 1) {
                printf("postprocessing...  ");
                fflush(stdout);
            }
            
#pragma omp parallel default(shared) private(n, i, j, pi, pj, pri, prj, ired, jred, r, ctheta0, nn, a, b, c, ex, ey, ez, cvala, aval, cvalb, nred)
#pragma omp for nowait
	    
            for (n = 0; n < set->spf.nrs; n++) {
                for (i = i0; i < i1; i++) {
                    for (j = j0; j < j1; j++) {
                        for (pi = set->spf.pimin; pi < set->spf.pimax; pi++) {
                            for (pj = set->spf.pjmin; pj < set->spf.pjmax; pj++) {
                                pri = (i1-i0)*pi;
                                prj = (j1-j0)*pj;
                                ired = i-i0;
                                jred = j-j0;

                                /*z0*/
                                if (!(set->spf.box_boundary_skipk0 && i >= set->spf.skipk0_imin && j >= set->spf.skipk0_jmin 
                                      && i <= set->spf.skipk0_imax && j <= set->spf.skipk0_jmax)) {
				    
                                    r = sqrt(((float)set->spf.ri[n] - (i+pri)) * ((float)set->spf.ri[n] - (i+pri))
                                             + ((float)set->spf.rj[n] - (j+prj)) * ((float)set->spf.rj[n] - (j+prj))
                                             + ((float)set->spf.rk[n] - k0) * ((float)set->spf.rk[n] - k0));
                                    ctheta0 = ((gdouble)(k0 - set->spf.rk[n])/r);
                                    r *= dx;
                                    a = mult/r;                                   //1.0/(4.0*G_PI*r);
                                    b = -mult*ctheta0/(r*r)*size;              //-1.0/(4*G_PI)*ctheta0/(r*r);
                                    c = -mult*ctheta0/(LIGHT_SPEED*r);    //-1.0/(4*G_PI)*ctheta0/(LIGHT_SPEED*r);

                                    cvala = c/2.0/dt*size;
                                    aval = a/2.0/dx*size;
                                    nred = set->spf.ppstart + 1 + r/(LIGHT_SPEED*dt);

                                    for (tit = 0; tit < (set->sc.nsteps - set->spf.ppstart); tit++) {
                                        nn = (gint)(tit + nred) - mp->pfarfield->prpoints[n].istart;

                                        if (nn < 3 || nn >= mp->pfarfield->ndata) {
                                            //printf("Warning: farfield ndata %d, trying access %d\n", mp->pfarfield->ndata, nn);
                                        }
					else {
                                            ex = mp->pfarfield->ex_k0->data[ired][jred][tit];
                                            cvalb = cvala*ex;
                                            mp->pfarfield->prpoints[n].ex->data[nn] += -cvalb;
                                            mp->pfarfield->prpoints[n].ex->data[nn-1] += (-aval*mp->pfarfield->dex_k0->data[ired][jred][tit] + b*ex);
                                            mp->pfarfield->prpoints[n].ex->data[nn-2] += cvalb;

                                            ey = mp->pfarfield->ey_k0->data[ired][jred][tit];
                                            cvalb = cvala*ey;
                                            mp->pfarfield->prpoints[n].ey->data[nn] += -cvalb;
                                            mp->pfarfield->prpoints[n].ey->data[nn-1] += (-aval*mp->pfarfield->dey_k0->data[ired][jred][tit] + b*ey);
                                            mp->pfarfield->prpoints[n].ey->data[nn-2] += cvalb;

                                            ez = mp->pfarfield->ez_k0->data[ired][jred][tit];
                                            cvalb = cvala*ez;
                                            mp->pfarfield->prpoints[n].ez->data[nn] += -cvalb;
                                            mp->pfarfield->prpoints[n].ez->data[nn-1] += (-aval*mp->pfarfield->dez_k0->data[ired][jred][tit] + b*ez);
                                            mp->pfarfield->prpoints[n].ez->data[nn-2] += cvalb;
                                        }
                                    }
                                }
                                /*z1*/
                                if (!(set->spf.box_boundary_skipkn && i >= set->spf.skipkn_imin && j >= set->spf.skipkn_jmin
                                      && i <= set->spf.skipkn_imax && j <= set->spf.skipkn_jmax)) {
                                    r = sqrt(((float)set->spf.ri[n] - (i+pri)) * ((float)set->spf.ri[n] - (i+pri))
                                             + ((float)set->spf.rj[n] - (j+prj)) * ((float)set->spf.rj[n] - (j+prj))
                                             + ((float)set->spf.rk[n] - k1) * ((float)set->spf.rk[n] - k1));
                                    ctheta0 = ((gdouble)(-(k1 - set->spf.rk[n]))/r);
                                    r *= dx;
                                    a = mult/r;                                   //1.0/(4.0*G_PI*r);
                                    b = -mult*ctheta0/(r*r)*size;              //-1.0/(4*G_PI)*ctheta0/(r*r);
                                    c = -mult*ctheta0/(LIGHT_SPEED*r);    //-1.0/(4*G_PI)*ctheta0/(LIGHT_SPEED*r);

                                    cvala = c/2.0/dt*size;
                                    aval = a/2.0/dx*size;
                                    nred = set->spf.ppstart + 1 + r/(LIGHT_SPEED*dt);

                                    for (tit = 0; tit < (set->sc.nsteps - set->spf.ppstart); tit++) {
                                        nn = (gint)(tit + nred) - mp->pfarfield->prpoints[n].istart;

                                        if (nn < 3 || nn >= mp->pfarfield->ndata) {
                                            //printf("Warning: farfield ndata %d, trying access %d\n", mp->pfarfield->ndata, nn);
                                        }
					else {
                                            ex = mp->pfarfield->ex_kn->data[ired][jred][tit];
                                            cvalb = cvala*ex;
                                            mp->pfarfield->prpoints[n].ex->data[nn] += -cvalb;
                                            mp->pfarfield->prpoints[n].ex->data[nn-1] += (aval*mp->pfarfield->dex_kn->data[ired][jred][tit] + b*ex);
                                            mp->pfarfield->prpoints[n].ex->data[nn-2] += cvalb;

                                            ey = mp->pfarfield->ey_kn->data[ired][jred][tit];
                                            cvalb = cvala*ey;
                                            mp->pfarfield->prpoints[n].ey->data[nn] += -cvalb;
                                            mp->pfarfield->prpoints[n].ey->data[nn-1] += (aval*mp->pfarfield->dey_kn->data[ired][jred][tit] + b*ey);
                                            mp->pfarfield->prpoints[n].ey->data[nn-2] += cvalb;

                                            ez = mp->pfarfield->ez_kn->data[ired][jred][tit];
                                            cvalb = cvala*ez;
                                            mp->pfarfield->prpoints[n].ez->data[nn] += -cvalb;
                                            mp->pfarfield->prpoints[n].ez->data[nn-1] += (aval*mp->pfarfield->dez_kn->data[ired][jred][tit] + b*ez);
                                            mp->pfarfield->prpoints[n].ez->data[nn-2] += cvalb;
                                        } // else
                                    } // tit
                                } // if (!(set->spf.skipkn ...
                            } // pj
                        } // pi
                    } // j
                } // i
            } // n
        } // if (set->sc.step_act == (set->sc.nsteps-1))
    } //  if (set->spf.postprocess)
    else {

#pragma omp parallel default(shared) private(n, i, j, pi, pj, pri, prj, r, ctheta0, nn, a, b, c, ex, ey, ez, cvala, aval, cvalb)
#pragma omp for nowait
	
        for (n = 0; n < set->spf.nrs; n++) {
	    for (i = i0; i < i1; i++) {
                for (j = j0; j < j1; j++) {
                    for (pi = set->spf.pimin; pi < set->spf.pimax; pi++) {
                        for (pj = set->spf.pjmin; pj < set->spf.pjmax; pj++) {
                            pri = (i1-i0)*pi;
                            prj = (j1-j0)*pj;

                            /*z0*/
                            if (!(set->spf.box_boundary_skipk0 && i >= set->spf.skipk0_imin && j >= set->spf.skipk0_jmin 
                                  && i <= set->spf.skipk0_imax && j <= set->spf.skipk0_jmax)) {
				
                                r = sqrt(((float)set->spf.ri[n] - (i+pri)) * ((float)set->spf.ri[n] - (i+pri))
                                         + ((float)set->spf.rj[n] - (j+prj)) * ((float)set->spf.rj[n] - (j+prj))
                                         + ((float)set->spf.rk[n] - k0) * ((float)set->spf.rk[n] - k0));
                                ctheta0 = ((gdouble)(k0 - set->spf.rk[n])/r);
                                r *= dx;

                                nn = (gint)(set->sc.step_act + 1 + r/(LIGHT_SPEED*dt)) - mp->pfarfield->prpoints[n].istart;

                                if (nn < 3 || nn >= mp->pfarfield->ndata) {
                                    //printf("Warning: farfield ndata %d, trying access %d\n", mp->pfarfield->ndata, nn);
                                }
				else {
                                    a = mult/r;                                   //1.0/(4.0*G_PI*r);
                                    b = -mult*ctheta0/(r*r);              //-1.0/(4*G_PI)*ctheta0/(r*r);
                                    c = -mult*ctheta0/(LIGHT_SPEED*r);    //-1.0/(4*G_PI)*ctheta0/(LIGHT_SPEED*r);

                                    cvala = c/2.0/dt*size;
                                    aval = a/2.0/dx;

                                    ex = mp->d->ex->data[i][j][k0];
                                    cvalb = cvala*ex;
                                    mp->pfarfield->prpoints[n].ex->data[nn] += -cvalb;
                                    mp->pfarfield->prpoints[n].ex->data[nn-1] += (-aval*(mp->d->ex->data[i][j][k0+1]-mp->d->ex->data[i][j][k0-1]) + b*ex)*size;
                                    mp->pfarfield->prpoints[n].ex->data[nn-2] += cvalb;

                                    ey = mp->d->ey->data[i][j][k0];
                                    cvalb = cvala*ey;
                                    mp->pfarfield->prpoints[n].ey->data[nn] += -cvalb;
                                    mp->pfarfield->prpoints[n].ey->data[nn-1] += (-aval*(mp->d->ey->data[i][j][k0+1]-mp->d->ey->data[i][j][k0-1]) + b*ey)*size;
                                    mp->pfarfield->prpoints[n].ey->data[nn-2] += cvalb;

                                    ez = mp->d->ez->data[i][j][k0];
                                    cvalb = cvala*ez;
                                    mp->pfarfield->prpoints[n].ez->data[nn] += -cvalb;
                                    mp->pfarfield->prpoints[n].ez->data[nn-1] += (-aval*(mp->d->ez->data[i][j][k0+1]-mp->d->ez->data[i][j][k0-1]) + b*ez)*size;
                                    mp->pfarfield->prpoints[n].ez->data[nn-2] += cvalb;
                                } // else
                            } // if
                            /*z1*/
                            if (!(set->spf.box_boundary_skipkn && i >= set->spf.skipkn_imin && j >= set->spf.skipkn_jmin 
                                  && i <= set->spf.skipkn_imax && j <= set->spf.skipkn_jmax)) {

                                r = sqrt(((float)set->spf.ri[n] - (i+pri)) * ((float)set->spf.ri[n] - (i+pri))
                                         + ((float)set->spf.rj[n] - (j+prj)) * ((float)set->spf.rj[n] - (j+prj))
                                         + ((float)set->spf.rk[n] - k1) * ((float)set->spf.rk[n] - k1));
                                ctheta0 = ((gdouble)(-(k1 - set->spf.rk[n]))/r);
                                r *= dx;

                                nn = (gint)(set->sc.step_act + 1 + r/(LIGHT_SPEED*dt)) - mp->pfarfield->prpoints[n].istart;

                                if (nn < 3 || nn >= mp->pfarfield->ndata) {
                                    //printf("Warning: farfield ndata %d, trying access %d\n", mp->pfarfield->ndata, nn);
                                }
				else {
                                    a = mult/r;                                   //1.0/(4.0*G_PI*r);
                                    b = -mult*ctheta0/(r*r);              //-1.0/(4*G_PI)*ctheta0/(r*r);
                                    c = -mult*ctheta0/(LIGHT_SPEED*r);    //-1.0/(4*G_PI)*ctheta0/(LIGHT_SPEED*r);

                                    cvala = c/2.0/dt*size;
                                    aval = a/2.0/dx;

                                    ex = mp->d->ex->data[i][j][k1];
                                    cvalb = cvala*ex;
                                    mp->pfarfield->prpoints[n].ex->data[nn] += -cvalb;
                                    mp->pfarfield->prpoints[n].ex->data[nn-1] += (aval*(mp->d->ex->data[i][j][k1+1]-mp->d->ex->data[i][j][k1-1]) + b*ex)*size;
                                    mp->pfarfield->prpoints[n].ex->data[nn-2] += cvalb;

                                    ey = mp->d->ey->data[i][j][k1];
                                    cvalb = cvala*ey;
                                    mp->pfarfield->prpoints[n].ey->data[nn] += -cvalb;
                                    mp->pfarfield->prpoints[n].ey->data[nn-1] += (aval*(mp->d->ey->data[i][j][k1+1]-mp->d->ey->data[i][j][k1-1]) + b*ey)*size;
                                    mp->pfarfield->prpoints[n].ey->data[nn-2] += cvalb;

                                    ez = mp->d->ez->data[i][j][k1];
                                    cvalb = cvala*ez;
                                    mp->pfarfield->prpoints[n].ez->data[nn] += -cvalb;
                                    mp->pfarfield->prpoints[n].ez->data[nn-1] += (aval*(mp->d->ez->data[i][j][k1+1]-mp->d->ez->data[i][j][k1-1]) + b*ez)*size;
                                    mp->pfarfield->prpoints[n].ez->data[nn-2] += cvalb;
                                } // else
                            } // if (!(set->spf.skipkn
                        } // pj
                    } // pi
                } // j
            } // i
        } // n
    } // else
    if (mp->set->sc.verbose > 1) printf("done.\n");

    //  printf("\n%g %g %g %d   %d %d %d %g\n", a, b, c, nn, ipr, jpr, kpr, probe);

    return 0;
}    

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
