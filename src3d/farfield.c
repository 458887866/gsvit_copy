
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
 *  near-to-far-field calculation algorithm based on Kirchoff surface integral
 */

#include "farfield.h"
#include "pool.h"
#include "settings.h"
#include "constants.h"
#include <math.h>
#include <omp.h>

static gdouble
dist(gdouble x1, gdouble y1, gdouble z1, gdouble x2, gdouble y2, gdouble z2)
{
    return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
}

static gdouble
angvec(gdouble x1, gdouble y1, gdouble z1, gdouble x2, gdouble y2, gdouble z2)
{
    return ((gdouble)(x1*x2 + y1*y2 + z1*z2) / (sqrt((gdouble)(x1*x1 + y1*y1 + z1*z1) * (x2*x2 + y2*y2 + z2*z2))));
}

static gdouble
angle(gdouble x1, gdouble y1, gdouble z1, gdouble x2, gdouble y2, gdouble z2, gint side) /*0 top, 1 bottom, 2 left, 3 right, 4 front, 5 back*/
{
    //cout << "angle: " << x2-x1 << " " << y2-y1 << " " << z2-z1 << endl;
    if (side == 0) return angvec(x2 - x1, y2 - y1, z2 - z1, 0, -1, 0);
    else if (side == 1) return angvec(x2 - x1, y2 - y1, z2 - z1, 0, 1, 0);
    else if (side == 2) return angvec(x2 - x1, y2 - y1, z2 - z1, -1, 0, 0);
    else if (side == 3) return angvec(x2 - x1, y2 - y1, z2 - z1, 1, 0, 0);
    else if (side == 4) return angvec(x2 - x1, y2 - y1, z2 - z1, 0, 0, -1);
    else return angvec(x2 - x1, y2 - y1, z2 - z1, 0, 0, 1);
}

SvFarfield*
sv_farfield_new(SvPool *mp, SvSet *set)
{
    gint n, actual_individual, last_istart;
    gdouble div;
    SvFarfield* ff = (SvFarfield*)g_malloc(sizeof(SvFarfield));

    ff->ndata = (gint)(2000 + 2*set->sc.nsteps + sqrt(set->sp.xres*set->sp.xres + set->sp.yres*set->sp.yres + set->sp.zres*set->sp.zres) /
		       (LIGHT_SPEED * set->plan.dt/set->sp.dx));
    ff->rpoints = (SvRPoint *) g_malloc(set->sf.nrs*sizeof(SvRPoint));

    if (set->sc.verbose > 1) {
        printf("Initializing far field positions data (%d points)... ", set->sf.nrs);
        fflush(stdout);
    }
    div = (LIGHT_SPEED*set->plan.dt/set->sp.dx);

    actual_individual = last_istart = 0;

    for (n = 0; n < set->sf.nrs; n++) {

        ff->rpoints[n].ex = gwy_data_line_new(ff->ndata, ff->ndata, TRUE);
        ff->rpoints[n].ey = gwy_data_line_new(ff->ndata, ff->ndata, TRUE);
        ff->rpoints[n].ez = gwy_data_line_new(ff->ndata, ff->ndata, TRUE);

        if (set->sf.individual[n] > 0) {
            if (set->sf.individual[n] != actual_individual) { //changing individual set: pick actual istart value for this set

	     	ff->rpoints[n].istart = (gint)(MAX(0, sqrt(((set->sf.ri[n]/div)*(set->sf.ri[n]/div) + (set->sf.rj[n]/div)*(set->sf.rj[n]/div) +
							    (set->sf.rk[n]/div)*(set->sf.rk[n]/div))) - 500));
		    last_istart = ff->rpoints[n].istart;
		    actual_individual = set->sf.individual[n];
        }
	    else
		    ff->rpoints[n].istart = last_istart;
        }
        else { //use independent value for all the independent farfield points
	        ff->rpoints[n].istart = (gint)(MAX(0, sqrt(((set->sf.ri[n]/div)*(set->sf.ri[n]/div) + (set->sf.rj[n]/div)*(set->sf.rj[n]/div) +
							(set->sf.rk[n]/div)*(set->sf.rk[n]/div))) - 500));
        }
        //printf("istart %d\n", ff->rpoints[n].istart);
    }

    if (set->sc.verbose > 1) printf("done.\n");

    return ff;
}

int
sv_pool_farfield_allocate_storage(SvPool *mp, SvSet *set)
{
    gint ix, iy, iz, i, j, k, n, x, y, z;
    gint i0, i1, j0, j1, k0, k1;
    SvFarfield* ff = mp->farfield;

    i0 = set->sf.box_i0;
    i1 = set->sf.box_in;
    j0 = set->sf.box_j0;
    j1 = set->sf.box_jn;
    k0 = set->sf.box_k0;
    k1 = set->sf.box_kn;

    ix = i1 - i0;
    iy = j1 - j0;
    iz = k1 - k0;

    for (n = 0; n < set->sf.nrs; n++) {

        ff->rpoints[n].r_left       = gwy_data_field_new(iy, iz, iy, iz, TRUE);
        ff->rpoints[n].ctheta_left  = gwy_data_field_new(iy, iz, iy, iz, TRUE);
        ff->rpoints[n].r_right      = gwy_data_field_new(iy, iz, iy, iz, TRUE);
        ff->rpoints[n].ctheta_right = gwy_data_field_new(iy, iz, iy, iz, TRUE);

        ff->rpoints[n].r_top         = gwy_data_field_new(ix, iz, ix, iz, TRUE);
        ff->rpoints[n].ctheta_top    = gwy_data_field_new(ix, iz, ix, iz, TRUE);
        ff->rpoints[n].r_bottom      = gwy_data_field_new(ix, iz, ix, iz, TRUE);
        ff->rpoints[n].ctheta_bottom = gwy_data_field_new(ix, iz, ix, iz, TRUE);

        ff->rpoints[n].r_front       = gwy_data_field_new(ix, iy, ix, iy, TRUE);
        ff->rpoints[n].ctheta_front  = gwy_data_field_new(ix, iy, ix, iy, TRUE);
        ff->rpoints[n].r_back        = gwy_data_field_new(ix, iy, ix, iy, TRUE);
        ff->rpoints[n].ctheta_back   = gwy_data_field_new(ix, iy, ix, iy, TRUE);

        x = set->sf.ri[n];
        y = set->sf.rj[n];
        z = set->sf.rk[n];

        for (j = j0; j < j1; j++) {
            for (k = k0; k < k1; k++) {
                gwy_data_field_set_val(ff->rpoints[n].r_left, j-j0, k-k0, dist(x, y, z, i0, j, k));
                gwy_data_field_set_val(ff->rpoints[n].ctheta_left, j-j0, k-k0, angle(i0, j, k, x, y, z, 2));

                gwy_data_field_set_val(ff->rpoints[n].r_right, j-j0, k-k0, dist(x, y, z, i1, j, k));
                gwy_data_field_set_val(ff->rpoints[n].ctheta_right, j-j0, k-k0, angle(i1, j, k, x, y, z, 3));
            }
        }
        for (i = i0; i < i1; i++) {
            for (k = k0; k < k1; k++) {
                gwy_data_field_set_val(ff->rpoints[n].r_top, i-i0, k-k0, dist(x, y, z, i, j0, k));
                gwy_data_field_set_val(ff->rpoints[n].ctheta_top, i-i0, k-k0, angle(i, j0, k, x, y, z, 0));

                gwy_data_field_set_val(ff->rpoints[n].r_bottom, i-i0, k-k0, dist(x, y, z, i, j1, k));
                gwy_data_field_set_val(ff->rpoints[n].ctheta_bottom, i-i0, k-k0, angle(i, j1, k, x, y, z, 1));
            }
        }
        for (i = i0; i < i1; i++) {
            for (j = j0; j < j1; j++) {
                gwy_data_field_set_val(ff->rpoints[n].r_front, i-i0, j-j0, dist(x, y, z, i, j, k0));
                gwy_data_field_set_val(ff->rpoints[n].ctheta_front, i-i0, j-j0, angle(i, j, k0, x, y, z, 4));

                gwy_data_field_set_val(ff->rpoints[n].r_back, i-i0, j-j0, dist(x, y, z, i, j, k1));
                gwy_data_field_set_val(ff->rpoints[n].ctheta_back, i-i0, j-j0, angle(i, j, k1, x, y, z, 5));
            }
        }
    } // n
    return 0;
}

int
sv_pool_farfield(SvPool *mp, SvSet *set)
{
    gint n = 0, i = 0, j = 0, k = 0;
    gint i0 = 0, i1 = 0, j0 = 0, j1 = 0, k0 = 0, k1 = 0;
    gdouble dx = set->sp.dx;
    gdouble dy = set->sp.dy;
    gdouble dz = set->sp.dz;
    gdouble dt = set->plan.dt;
    gdouble size = 0, r = 0, a = 0, b = 0, c = 0;
    gint nn = 0;
    GwyDataField *r0 = NULL, *r1 = NULL, *ctheta0 = NULL, *ctheta1 = NULL;

    if (set->sf.nrs == 0) 
        return 0;

    i0 = set->sf.box_i0;
    i1 = set->sf.box_in;
    j0 = set->sf.box_j0;
    j1 = set->sf.box_jn;
    k0 = set->sf.box_k0;
    k1 = set->sf.box_kn;

    if (set->sc.verbose>1) {
        printf("Running farfield... \n");
        fflush(stdout);
    }

#pragma omp parallel default(shared) private(n, i, j, k, r, r0, r1, ctheta0, ctheta1, nn, a, b, c, size)
#pragma omp for nowait
    
    for (n = 0; n < set->sf.nrs; n++) {

        /*perform the summation*/
        r0 = mp->farfield->rpoints[n].r_left;
        r1 = mp->farfield->rpoints[n].r_right;
        ctheta0 = mp->farfield->rpoints[n].ctheta_left;
        ctheta1 = mp->farfield->rpoints[n].ctheta_right;

        for (j = j0; j < j1; j++) {
            for (k = k0; k < k1; k++) {
                size = dy*dz;
                /*x0*/
                if (!(set->sf.box_boundary_skipi0 && j >= set->sf.skipi0_jmin && k >= set->sf.skipi0_kmin && j <= set->sf.skipi0_jmax && k <= set->sf.skipi0_kmax)) {
                    r = gwy_data_field_get_val(r0, j-j0, k-k0)*dx;
                    nn = (int)(set->sc.step_act + 1 + r/(LIGHT_SPEED*dt)) - mp->farfield->rpoints[n].istart;
                    if (nn < 3 || nn >= mp->farfield->ndata) {
                        //printf("Warning: farfield ndata %d, trying access %d\n", mp->farfield->ndata, nn);
                    }
		    else {
                        a = 1.0/(4.0*G_PI*r);
                        b = -1.0/(4.0*G_PI)*gwy_data_field_get_val(ctheta0, j-j0, k-k0)/(r*r);
                        c = -1.0/(4.0*G_PI)*gwy_data_field_get_val(ctheta0, j-j0, k-k0)/(LIGHT_SPEED*r);

                        mp->farfield->rpoints[n].ex->data[nn] += -c/2.0/dt*mp->d->ex->data[i0][j][k]*size;
                        mp->farfield->rpoints[n].ex->data[nn-1] += (-a*(mp->d->ex->data[i0+1][j][k]-mp->d->ex->data[i0-1][j][k])/2.0/dx + b*mp->d->ex->data[i0][j][k])*size;
                        mp->farfield->rpoints[n].ex->data[nn-2] += c/2.0/dt*mp->d->ex->data[i0][j][k]*size;

                        mp->farfield->rpoints[n].ey->data[nn] += -c/2.0/dt*mp->d->ey->data[i0][j][k]*size;
                        mp->farfield->rpoints[n].ey->data[nn-1] += (-a*(mp->d->ey->data[i0+1][j][k]-mp->d->ey->data[i0-1][j][k])/2.0/dx + b*mp->d->ey->data[i0][j][k])*size;
                        mp->farfield->rpoints[n].ey->data[nn-2] += c/2.0/dt*mp->d->ey->data[i0][j][k]*size;

                        mp->farfield->rpoints[n].ez->data[nn] += -c/2.0/dt*mp->d->ez->data[i0][j][k]*size;
                        mp->farfield->rpoints[n].ez->data[nn-1] += (-a*(mp->d->ez->data[i0+1][j][k]-mp->d->ez->data[i0-1][j][k])/2.0/dx + b*mp->d->ez->data[i0][j][k])*size;
                        mp->farfield->rpoints[n].ez->data[nn-2] += c/2.0/dt*mp->d->ez->data[i0][j][k]*size;
                    }
                }

                /*xn*/
                if (!(set->sf.box_boundary_skipin && j >= set->sf.skipin_jmin && k >= set->sf.skipin_kmin && j <= set->sf.skipin_jmax && k <= set->sf.skipin_kmax)) {

                    r = gwy_data_field_get_val(r1, j-j0, k-k0)*dx;
                    nn = (int)(set->sc.step_act + 1 + r/(LIGHT_SPEED*dt)) - mp->farfield->rpoints[n].istart;
                    if (nn < 3 || nn >= mp->farfield->ndata) {
                        //printf("Warning: farfield ndata %d, trying access %d\n", mp->farfield->ndata, nn);
                    }
		    else {
                        a = 1.0/(4.0*G_PI*r);
                        b = -1.0/(4.0*G_PI)*gwy_data_field_get_val(ctheta1, j-j0, k-k0)/(r*r);
                        c = -1.0/(4.0*G_PI)*gwy_data_field_get_val(ctheta1, j-j0, k-k0)/(LIGHT_SPEED*r);

                        mp->farfield->rpoints[n].ex->data[nn] += -c/2.0/dt*mp->d->ex->data[i1][j][k]*size;
                        mp->farfield->rpoints[n].ex->data[nn-1] += (a*(mp->d->ex->data[i1+1][j][k]-mp->d->ex->data[i1-1][j][k])/2.0/dx + b*mp->d->ex->data[i1][j][k])*size;
                        mp->farfield->rpoints[n].ex->data[nn-2] += c/2.0/dt*mp->d->ex->data[i1][j][k]*size;

                        mp->farfield->rpoints[n].ey->data[nn] += -c/2.0/dt*mp->d->ey->data[i1][j][k]*size;
                        mp->farfield->rpoints[n].ey->data[nn-1] += (a*(mp->d->ey->data[i1+1][j][k]-mp->d->ey->data[i1-1][j][k])/2.0/dx + b*mp->d->ey->data[i1][j][k])*size;
                        mp->farfield->rpoints[n].ey->data[nn-2] += c/2.0/dt*mp->d->ey->data[i1][j][k]*size;

                        mp->farfield->rpoints[n].ez->data[nn] += -c/2.0/dt*mp->d->ez->data[i1][j][k]*size;
                        mp->farfield->rpoints[n].ez->data[nn-1] += (a*(mp->d->ez->data[i1+1][j][k]-mp->d->ez->data[i1-1][j][k])/2.0/dx + b*mp->d->ez->data[i1][j][k])*size;
                        mp->farfield->rpoints[n].ez->data[nn-2] += c/2.0/dt*mp->d->ez->data[i1][j][k]*size;
                    }
                }
            } // k
        } // j
	
        /*y planes*/
        r0 = mp->farfield->rpoints[n].r_top;
        r1 = mp->farfield->rpoints[n].r_bottom;
        ctheta0 = mp->farfield->rpoints[n].ctheta_top;
        ctheta1 = mp->farfield->rpoints[n].ctheta_bottom;

        for (i = i0; i < i1; i++) {
            for (k = k0; k < k1; k++) {
                size = dx*dz;

                /*y0*/
                if (!(set->sf.box_boundary_skipj0 && i >= set->sf.skipj0_imin && k >= set->sf.skipj0_kmin && i <= set->sf.skipj0_imax && k <= set->sf.skipj0_kmax)) {
                    r = gwy_data_field_get_val(r0, i-i0, k-k0)*dx;
                    nn = (int)(set->sc.step_act + 1 + r/(LIGHT_SPEED*dt)) - mp->farfield->rpoints[n].istart;
                    if (nn < 3 || nn >= mp->farfield->ndata) {
                        //printf("Warning: farfield ndata %d, trying access %d\n", mp->farfield->ndata, nn);
                    }
		    else {
                        a = 1.0/(4.0*G_PI*r);
                        b = -1.0/(4*G_PI)*gwy_data_field_get_val(ctheta0, i-i0, k-k0)/(r*r);
                        c = -1.0/(4*G_PI)*gwy_data_field_get_val(ctheta0, i-i0, k-k0)/(LIGHT_SPEED*r);

                        mp->farfield->rpoints[n].ex->data[nn] += -c/2.0/dt*mp->d->ex->data[i][j0][k]*size;
                        mp->farfield->rpoints[n].ex->data[nn-1] += (-a*(mp->d->ex->data[i][j0+1][k]-mp->d->ex->data[i][j0-1][k])/2.0/dx + b*mp->d->ex->data[i][j0][k])*size;
                        mp->farfield->rpoints[n].ex->data[nn-2] += c/2.0/dt*mp->d->ex->data[i][j0][k]*size;

                        mp->farfield->rpoints[n].ey->data[nn] += -c/2.0/dt*mp->d->ey->data[i][j0][k]*size;
                        mp->farfield->rpoints[n].ey->data[nn-1] += (-a*(mp->d->ey->data[i][j0+1][k]-mp->d->ey->data[i][j0-1][k])/2.0/dx + b*mp->d->ey->data[i][j0][k])*size;
                        mp->farfield->rpoints[n].ey->data[nn-2] += c/2.0/dt*mp->d->ey->data[i][j0][k]*size;

                        mp->farfield->rpoints[n].ez->data[nn] += -c/2.0/dt*mp->d->ez->data[i][j0][k]*size;
                        mp->farfield->rpoints[n].ez->data[nn-1] += (-a*(mp->d->ez->data[i][j0+1][k]-mp->d->ez->data[i][j0-1][k])/2.0/dx + b*mp->d->ez->data[i][j0][k])*size;
                        mp->farfield->rpoints[n].ez->data[nn-2] += c/2.0/dt*mp->d->ez->data[i][j0][k]*size;
                    }
                }

                /*yn*/
                if (!(set->sf.box_boundary_skipjn && i >= set->sf.skipjn_imin && k >= set->sf.skipjn_kmin && i <= set->sf.skipjn_imax && k <= set->sf.skipjn_kmax)) {
                    r = gwy_data_field_get_val(r1, i-i0, k-k0)*dx;
                    nn = (int)(set->sc.step_act + 1 + r/(LIGHT_SPEED*dt)) - mp->farfield->rpoints[n].istart;
                    if (nn < 3 || nn >= mp->farfield->ndata) {
                        //printf("Warning: farfield ndata %d, trying access %d\n", mp->farfield->ndata, nn);
                    }
		    else {
                        a = 1.0/(4.0*G_PI*r);
                        b = -1.0/(4*G_PI)*gwy_data_field_get_val(ctheta1, i-i0, k-k0)/(r*r);
                        c = -1.0/(4*G_PI)*gwy_data_field_get_val(ctheta1, i-i0, k-k0)/(LIGHT_SPEED*r);

                        mp->farfield->rpoints[n].ex->data[nn] += -c/2.0/dt*mp->d->ex->data[i][j1][k]*size;
                        mp->farfield->rpoints[n].ex->data[nn-1] += (a*(mp->d->ex->data[i][j1+1][k]-mp->d->ex->data[i][j1-1][k])/2.0/dx + b*mp->d->ex->data[i][j1][k])*size;
                        mp->farfield->rpoints[n].ex->data[nn-2] += c/2.0/dt*mp->d->ex->data[i][j1][k]*size;

                        mp->farfield->rpoints[n].ey->data[nn] += -c/2.0/dt*mp->d->ey->data[i][j1][k]*size;
                        mp->farfield->rpoints[n].ey->data[nn-1] += (a*(mp->d->ey->data[i][j1+1][k]-mp->d->ey->data[i][j1-1][k])/2.0/dx + b*mp->d->ey->data[i][j1][k])*size;
                        mp->farfield->rpoints[n].ey->data[nn-2] += c/2.0/dt*mp->d->ey->data[i][j1][k]*size;

                        mp->farfield->rpoints[n].ez->data[nn] += -c/2.0/dt*mp->d->ez->data[i][j1][k]*size;
                        mp->farfield->rpoints[n].ez->data[nn-1] += (a*(mp->d->ez->data[i][j1+1][k]-mp->d->ez->data[i][j1-1][k])/2.0/dx + b*mp->d->ez->data[i][j1][k])*size;
                        mp->farfield->rpoints[n].ez->data[nn-2] += c/2.0/dt*mp->d->ez->data[i][j1][k]*size;
                    }
                }
            } // k
        } // i
	
        /*z planes*/
        r0 = mp->farfield->rpoints[n].r_front;
        r1 = mp->farfield->rpoints[n].r_back;
        ctheta0 = mp->farfield->rpoints[n].ctheta_front;
        ctheta1 = mp->farfield->rpoints[n].ctheta_back;

        for (i = i0; i < i1; i++) {
            for (j = j0; j < j1; j++) {
                size = dx*dy;

                /*z0*/
                if (!(set->sf.box_boundary_skipk0 && i >= set->sf.skipk0_imin && j >= set->sf.skipk0_jmin && i <= set->sf.skipk0_imax && j <= set->sf.skipk0_jmax)) {
                    r = gwy_data_field_get_val(r0, i-i0, j-j0)*dx;
                    nn = (int)(set->sc.step_act + 1 + r/(LIGHT_SPEED*dt)) - mp->farfield->rpoints[n].istart;
                    if (nn < 3 || nn >= mp->farfield->ndata) {
                        //printf("Warning: farfield ndata %d, trying access %d\n", mp->farfield->ndata, nn);
                    }
		    else {
                        a = 1.0/(4.0*G_PI*r);
                        b = -1.0/(4*G_PI)*gwy_data_field_get_val(ctheta0, i-i0, j-j0)/(r*r);
                        c = -1.0/(4*G_PI)*gwy_data_field_get_val(ctheta0, i-i0, j-j0)/(LIGHT_SPEED*r);

                        mp->farfield->rpoints[n].ex->data[nn] += -c/2.0/dt*mp->d->ex->data[i][j][k0]*size;
                        mp->farfield->rpoints[n].ex->data[nn-1] += (-a*(mp->d->ex->data[i][j][k0+1]-mp->d->ex->data[i][j][k0-1])/2.0/dx + b*mp->d->ex->data[i][j][k0])*size;
                        mp->farfield->rpoints[n].ex->data[nn-2] += c/2.0/dt*mp->d->ex->data[i][j][k0]*size;

                        mp->farfield->rpoints[n].ey->data[nn] += -c/2.0/dt*mp->d->ey->data[i][j][k0]*size;
                        mp->farfield->rpoints[n].ey->data[nn-1] += (-a*(mp->d->ey->data[i][j][k0+1]-mp->d->ey->data[i][j][k0-1])/2.0/dx + b*mp->d->ey->data[i][j][k0])*size;
                        mp->farfield->rpoints[n].ey->data[nn-2] += c/2.0/dt*mp->d->ey->data[i][j][k0]*size;

                        mp->farfield->rpoints[n].ez->data[nn] += -c/2.0/dt*mp->d->ez->data[i][j][k0]*size;
                        mp->farfield->rpoints[n].ez->data[nn-1] += (-a*(mp->d->ez->data[i][j][k0+1]-mp->d->ez->data[i][j][k0-1])/2.0/dx + b*mp->d->ez->data[i][j][k0])*size;
                        mp->farfield->rpoints[n].ez->data[nn-2] += c/2.0/dt*mp->d->ez->data[i][j][k0]*size;
                    }
                }

                /*zn*/
                if (!(set->sf.box_boundary_skipkn && i >= set->sf.skipkn_imin && j >= set->sf.skipkn_jmin && i <= set->sf.skipkn_imax && j <= set->sf.skipkn_jmax)) {
                    r = gwy_data_field_get_val(r1, i-i0, j-j0)*dx;
                    nn = (int)(set->sc.step_act + 1 + r/(LIGHT_SPEED*dt)) - mp->farfield->rpoints[n].istart;
                    if (nn < 3 || nn >= mp->farfield->ndata) {
                        //printf("Warning: farfield ndata %d, trying access %d\n", mp->farfield->ndata, nn);
                    }
		    else {
                        a = 1.0/(4.0*G_PI*r);
                        b = -1.0/(4*G_PI)*gwy_data_field_get_val(ctheta1, i-i0, j-j0)/(r*r);
                        c = -1.0/(4*G_PI)*gwy_data_field_get_val(ctheta1, i-i0, j-j0)/(LIGHT_SPEED*r);

                        mp->farfield->rpoints[n].ex->data[nn] += -c/2.0/dt*mp->d->ex->data[i][j][k1]*size;
                        mp->farfield->rpoints[n].ex->data[nn-1] += (a*(mp->d->ex->data[i][j][k1+1]-mp->d->ex->data[i][j][k1-1])/2.0/dx + b*mp->d->ex->data[i][j][k1])*size;
                        mp->farfield->rpoints[n].ex->data[nn-2] += c/2.0/dt*mp->d->ex->data[i][j][k1]*size;

                        mp->farfield->rpoints[n].ey->data[nn] += -c/2.0/dt*mp->d->ey->data[i][j][k1]*size;
                        mp->farfield->rpoints[n].ey->data[nn-1] += (a*(mp->d->ey->data[i][j][k1+1]-mp->d->ey->data[i][j][k1-1])/2.0/dx + b*mp->d->ey->data[i][j][k1])*size;
                        mp->farfield->rpoints[n].ey->data[nn-2] += c/2.0/dt*mp->d->ey->data[i][j][k1]*size;

                        mp->farfield->rpoints[n].ez->data[nn] += -c/2.0/dt*mp->d->ez->data[i][j][k1]*size;
                        mp->farfield->rpoints[n].ez->data[nn-1] += (a*(mp->d->ez->data[i][j][k1+1]-mp->d->ez->data[i][j][k1-1])/2.0/dx + b*mp->d->ez->data[i][j][k1])*size;
                        mp->farfield->rpoints[n].ez->data[nn-2] += c/2.0/dt*mp->d->ez->data[i][j][k1]*size;

			// if (i==180 && j==170) printf("theta pfff %g r %g nn %d ex %g exorig %g abc %g %g %g\n", gwy_data_field_get_val(ctheta1, i-i0, j-j0), r, nn, mp->farfield->rpoints[n].ex->data[nn], mp->d->ex->data[i][j][k1], a, b, c);
                    }
                }
            } // j
        } // i
    } // n

    if (set->sc.verbose > 1) printf("done.\n");

    return 0;
}    

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
