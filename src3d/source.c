
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


 /*  source.c :
  *  algorithms for creating elmag field sources
  */

#include "source.h"
#include "pool.h"
#include "settings.h"
#include "constants.h"
#include <math.h>
#include <omp.h>
#if defined (_WIN32) && defined (_MSC_VER)
#include "../config_win.h"
#else
#include "../config.h"
#include <sys/time.h>
#endif

SvSource*
sv_source_new()
{
    SvSource* src = (SvSource*)g_malloc(sizeof(SvSource));

    src->tsf = NULL;
    src->ltsf = NULL;
    src->tsff = NULL;
    src->ltsff = NULL;
    src->sf = NULL;
    src->sh = NULL;
    return src;
}

gdouble
dcomp(gint i, gint j, gint k, gint xres, gint yres, gint zres, gdouble theta, gdouble phi,
      gint i0, gint i1, gint j0, gint j1, gint k0, gint k1)
{
    gdouble rx, ry, rz;
    gdouble ax, ay, az;

    ax = sin(theta)*cos(phi);
    ay = sin(theta)*sin(phi);
    az = cos(theta);

    rx = 0;

    if (theta >= 0 && theta <= (G_PI / 2.0)) {
        if (phi >= 0 && phi <= (G_PI / 2.0)) {
            rx = i - i0;
            ry = j - j0;
            rz = k - k0;
        } else if (phi > (G_PI / 2.0) && phi <= G_PI) {
            rx = i - i1;
            ry = j - j0;
            rz = k - k0;
        } else if (phi > G_PI && phi <= (3.0*G_PI / 2.0)) {
            rx = i - i1;
            ry = j - j1;
            rz = k - k0;
        } else {
            rx = i - i0;
            ry = j - j1;
            rz = k - k0;
        }
    } else if (theta < G_PI && theta >(G_PI / 2.0)) {
        if (phi >= 0 && phi <= (G_PI / 2.0)) {
            rx = i - i0;
            ry = j - j0;
            rz = k - k1;
        } else if (phi > (G_PI / 2.0) && phi <= G_PI) {
            ry = j - j0;
            rz = k - k1;
        } else if (phi > G_PI && phi <= (3.0*G_PI / 2.0)) {
            rx = i - i1;
            ry = j - j1;
            rz = k - k1;
        } else {
            rx = i - i0;
            ry = j - j1;
            rz = k - k1;
        }
    } else {
        fprintf(stderr, "Theta angles out of interval 0-PI unsupported\n");
        return -1;
    }

    return 10 + (ax*rx + ay*ry + az*rz);
}

gdouble
rdcomp(gint i, gint j, gint k, gint xres, gint yres, gint zres, gdouble theta, gdouble phi,
       gint i0, gint i1, gint j0, gint j1, gint k0, gint k1)
{

    gdouble rx, ry, rz;
    gdouble ax, ay, az;

    ax = sin(theta)*cos(phi);
    ay = sin(theta)*sin(phi);
    az = cos(theta);

    rx = -i + xres / 2;
    ry = -j + yres / 2;
    rz = -k + zres / 2;

    if ((ax*rx + ay*ry + az*rz) >= 200)
        fprintf(stderr, "Error! Negative field position (%g)!\n", (ax*rx + ay*ry + az*rz));
    return 200 - (ax*rx + ay*ry + az*rz);
}

gdouble
gex(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return field*(cos(psi)*sin(phi) - sin(psi)*cos(theta)*cos(phi));
}

gdouble
gey(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return field*(-cos(psi)*cos(phi) - sin(psi)*cos(theta)*sin(phi));
}

gdouble
gez(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return field*(sin(psi)*sin(theta));
}

gdouble
ghx(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return field*(sin(psi)*sin(phi) + cos(psi)*cos(theta)*cos(phi));
}

gdouble
ghy(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return field*(-sin(psi)*cos(phi) + cos(psi)*cos(theta)*sin(phi));
}

gdouble
ghz(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return field*(-cos(psi)*sin(theta));
}

gdouble
fs_amplitude(gdouble freal, gdouble an, gdouble bm, gdouble thetan)
{
    return 1e15*an*bm*(freal / (2 * G_PI*LIGHT_SPEED))*sqrt(cos(thetan))*sin(thetan);
}

gdouble
srdcomp(gint i, gint j, gint k, gint xres, gint yres, gint zres, gdouble stheta, gdouble sphi, gdouble ctheta, gdouble cphi,
        gint i0, gint i1, gint j0, gint j1, gint k0, gint k1, gdouble xshift, gdouble yshift, gdouble zshift)
{

    gdouble rx, ry, rz;
    gdouble ax, ay, az;

    ax = stheta*cphi;
    ay = stheta*sphi;
    az = ctheta;

    rx = -i + xres / 2 + xshift;
    ry = -j + yres / 2 + yshift;
    rz = -k + zres / 2 + zshift;

    if ((ax*rx + ay*ry + az*rz) >= 250)
        fprintf(stderr, "Error! Negative field position (%g)!\n", (ax*rx + ay*ry + az*rz));

    return 250 - (ax*rx + ay*ry + az*rz);
}

gdouble
sgex(gdouble field, gdouble stheta, gdouble sphi, gdouble spsi, gdouble ctheta, gdouble cphi, gdouble cpsi)
{
    return field*(cpsi*sphi - spsi*ctheta*cphi);
}

gdouble
sgey(gdouble field, gdouble stheta, gdouble sphi, gdouble spsi, gdouble ctheta, gdouble cphi, gdouble cpsi)
{
    return field*(-cpsi*cphi - spsi*ctheta*sphi);
}

gdouble
sgez(gdouble field, gdouble stheta, gdouble sphi, gdouble spsi, gdouble ctheta, gdouble cphi, gdouble cpsi)
{
    return field*(spsi*stheta);
}

gdouble
sghx(gdouble field, gdouble stheta, gdouble sphi, gdouble spsi, gdouble ctheta, gdouble cphi, gdouble cpsi)
{
    return field*(spsi*sphi + cpsi*ctheta*cphi);
}

gdouble
sghy(gdouble field, gdouble stheta, gdouble sphi, gdouble spsi, gdouble ctheta, gdouble cphi, gdouble cpsi)
{
    return field*(-spsi*cphi + cpsi*ctheta*sphi);
}

gdouble
sghz(gdouble field, gdouble stheta, gdouble sphi, gdouble spsi, gdouble ctheta, gdouble cphi, gdouble cpsi)
{
    return field*(-cpsi*stheta);
}

gdouble
gaussmult(gint x, gint y, gint z, gdouble xc, gdouble yc, gdouble rx, gdouble ry)
{
    return exp(-((gdouble)(x - xc))*((gdouble)(x - xc)) / rx / rx / 2.0) * exp(-((gdouble)(y - yc))*((gdouble)(y - yc)) / ry / ry / 2.0);
}

gdouble
fibermult(gint x, gint y, gint z, gdouble xc, gdouble yc, gdouble radius, gdouble cutoff,
          gdouble epsilon_core, gdouble epsilon_cladding, gdouble dx, gdouble lambda)
{
    gdouble v, omegazero, radval;

    v = 2 * G_PI / lambda * radius*dx*(epsilon_core - epsilon_cladding);
    omegazero = (0.65 + 1.619 / v / sqrt(v) + 2.879 / pow(v, 6)) * radius*dx;
    radval = ((x - xc)*(x - xc) + (y - yc)*(y - yc));

    if (radval < (cutoff*cutoff))
        return exp(-radval*dx*dx / omegazero / omegazero);

    return 0;
}


gboolean mat_too_close_check(SvPool *mp, SvSet *set, gint ipos, gint jpos, gint kpos)
{
    if (set->ss.tsf.box_boundary_skipdepth == -1 && set->ss.ltsf.box_boundary_skipdepth == -1)
        return TRUE;

    /*use precomputed data*/
    if (mp->d->sourceskip->data[ipos][jpos][kpos])
        return FALSE;

    /*everything is OK, no material around*/
    return TRUE;
}

/*
void
getcoefs(gint i, gint j, gint k, gdouble theta, gdouble phi, gint step, gint *zpos, gdouble *timepos, gdouble dx, gdouble dt)
{
    double alpha = atan2(j, i);
    double radius = sqrt(i*i + j*j);
    radius *= cos(alpha-phi);               //was radius - i originally for phi=0

//    if (i==100 && j==100)
//    printf("rad %g %g %g %g\n", (double)sqrt(i*i + j*j), cos(alpha-phi), radius, (gdouble)step - radius*dx*sin(theta)/(LIGHT_SPEED*dt));

    *zpos = k;
    *timepos = (gdouble)step - radius*dx*sin(theta)/(LIGHT_SPEED*dt);

    if (*timepos<0) *timepos = 0;

    //printf("%d %g/%g = %g, %g  %d\n", step, ((gdouble)i)*dx, (LIGHT_SPEED*dt),  ((gdouble)i)*dx*sin(theta)/(LIGHT_SPEED*dt), *timepos, i);
}*/

void
getcoefs(gint i, gint j, gint k, gdouble theta, gdouble phi, gint step, gint *zpos, gint *tp, gdouble *atp, gdouble dx, gdouble dt,
         gint i0, gint j0, gint k0, gint i1, gint j1, gint k1)
{
    gint rx = 0, ry = 0;
    gdouble alpha;
    gdouble radius;
    gdouble timepos;

    if (theta >= 0 && theta <= (G_PI / 2.0)) {
        if (phi >= 0 && phi <= (G_PI / 2.0)) {
            rx = i - i0;
            ry = j - j0;
        } else if (phi > (G_PI / 2.0) && phi <= G_PI) {
            rx = i - i1;
            ry = j - j0;
        } else if (phi > G_PI && phi <= (3.0*G_PI / 2.0)) {
            rx = i - i1;
            ry = j - j1;
        } else {
            rx = i - i0;
            ry = j - j1;
        }
    } else
        fprintf(stderr, "Error: unsupported LTSF angle!\n");

    radius = sqrt(rx*rx + ry*ry);
    alpha = atan2(ry, rx);
    radius *= cos(alpha - phi);

    *zpos = k; //rz;
    timepos = 3.0*((gdouble)step - radius*dx*sin(theta) / (LIGHT_SPEED*dt));

    *tp = (gint)timepos;
    *atp = timepos - (gdouble)(*tp);
    if (*tp < 0)
        *tp = 0;
}

void
sgetcoefs(gint i, gint j, gint k, gdouble theta, gdouble phi, gint step, gint *zpos, gint *tp, gdouble *atp, gdouble dx, gdouble dt,
          gint i0, gint j0, gint k0, gint i1, gint j1, gint k1, gint xres, gint yres, gint zres)
{
    gint rx = 0, ry = 0;
    gdouble alpha;
    gdouble radius;
    gdouble timepos;

    if (theta >= 0 && theta <= (G_PI / 2.0)) {
        if (phi >= 0 && phi <= (G_PI / 2.0)) {
            rx = i - i0;
            ry = j - j0;
        } else if (phi > (G_PI / 2.0) && phi <= G_PI) {
            rx = i - i1;
            ry = j - j0;
        } else if (phi > G_PI && phi <= (3.0*G_PI / 2.0)) {
            rx = i - i1;
            ry = j - j1;
        } else {
            rx = i - i0;
            ry = j - j1;
        }
    } else
        fprintf(stderr, "Error: unsupported LTSF angle!\n");

    radius = sqrt(rx*rx + ry*ry);
    alpha = atan2(ry, rx);
    radius *= cos(alpha - phi);

    *zpos = k;
    timepos = 3.0*((gdouble)step - radius*dx*sin(theta) / (LIGHT_SPEED*dt));

    *tp = (gint)timepos;
    *atp = timepos - (gdouble)(*tp);
    if ((*tp) < 0)
        *tp = 0;
}

static gdouble
get_gaussian_rnd(GRand *rnd, gdouble center, gdouble width)
{
    double x1, x2, w;

    do {
        x1 = 2.0 * g_rand_double_range(rnd, 0, 1) - 1.0;
        x2 = 2.0 * g_rand_double_range(rnd, 0, 1) - 1.0;
        w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);

    return center + x1*w*width;
}

void
sv_pool_apply_source_estep(SvPool *mp, SvSet *set)
{
    gint i = 0, j = 0, k = 0, m = 0, pos = 0, im = 0, i1 = 0, zpos = 0, tp = 0, sm = 0;
    gdouble d = 0, vh = 0, ecor = 0, theta = 0, phi = 0, psi = 0, corr = 0, acorr = 0, an = 0, bm = 0, epsilon = 0, sigma = 0;
    gdouble stheta, sphi, spsi, ctheta, cphi, cpsi, val;
    gdouble hx = 0, hy = 0, hz = 0, atp = 0, vhh = 0, vhe = 0, gecor = 0;
    gchar filename[100];
    FILE *fhx, *fhy, *fhz;
    gfloat vhx, vhy, vhz;
    GwyDataField *hxs, *hys, *hzs;
    gdouble *hxsdata, *hysdata, *hzsdata;
    gint ipos, jpos;
    gdouble xdir, ydir, zdir, mag;
    gdouble rphi = 0, xradmult = 0, yradmult = 0;

    GRand *rnd;
#ifndef G_OS_WIN32
    struct timeval time;
#endif

    if (set->sc.verbose > 1) {
        printf("Running source estep...   ");
        fflush(stdout);
    }

    for (i = 0; i < set->ss.npnts; i++) {
        /*if we are lucky, data are not shifted by user and can be used directly*/
        pos = -1;
        if (set->sc.step_act <= mp->src->sp[i].sdata.ndata && mp->src->sp[i].sdata.layered_zpos[set->sc.step_act] == set->sc.step_act)
            pos = set->sc.step_act;
        else { /*otherwise we need to search for right position*/
            for (j = 0; j < mp->src->sp[i].sdata.ndata; j++) {
                if (mp->src->sp[i].sdata.layered_zpos[j] == set->sc.step_act) {
                    pos = j;
                    break;
                }
            } // j
        }
        if (pos >= 0 && pos < mp->src->sp[i].sdata.ndata) {


            if (mp->src->sp[i].sdata.ex[pos] != 0)
                mp->d->ex->data[mp->src->sp[i].i][mp->src->sp[i].j][mp->src->sp[i].k] = mp->src->sp[i].sdata.ex[pos];
            if (mp->src->sp[i].sdata.ey[pos] != 0)
                mp->d->ey->data[mp->src->sp[i].i][mp->src->sp[i].j][mp->src->sp[i].k] = mp->src->sp[i].sdata.ey[pos];
            if (mp->src->sp[i].sdata.ez[pos] != 0)
                mp->d->ez->data[mp->src->sp[i].i][mp->src->sp[i].j][mp->src->sp[i].k] = mp->src->sp[i].sdata.ez[pos];


//            printf("%g %g %g\n", mp->src->sp[i].sdata.ex[pos], mp->src->sp[i].sdata.ey[pos], mp->src->sp[i].sdata.ez[pos]);
//            if (fabs(mp->src->sp[i].sdata.ex[pos]) > 1e-10)
//                mp->d->ex->data[mp->src->sp[i].i][mp->src->sp[i].j][mp->src->sp[i].k] += mp->src->sp[i].sdata.ex[pos];
//            if (fabs(mp->src->sp[i].sdata.ey[pos]) > 1e-10)
//                mp->d->ey->data[mp->src->sp[i].i][mp->src->sp[i].j][mp->src->sp[i].k] += mp->src->sp[i].sdata.ey[pos];
//            if (fabs(mp->src->sp[i].sdata.ez[pos]) > 1e-10)
//                mp->d->ez->data[mp->src->sp[i].i][mp->src->sp[i].j][mp->src->sp[i].k] += mp->src->sp[i].sdata.ez[pos];
        }
    } // i

    if (set->ss.nlocals > 0) {
        if (set->sc.step_act == 0) {
            if (set->sc.verbose > 1) {
                printf("generating local sources...   ");
                fflush(stdout);
            }
            mp->src->sl = (SvSourceLocal *)g_malloc(set->ss.nlocals * sizeof(SvSourceLocal));
        }

#ifdef G_OS_WIN32
        rnd = g_rand_new();
        g_rand_set_seed(rnd, g_random_int() & 0x7fffffff);
#else
        gettimeofday(&time, NULL);
        rnd = g_rand_new_with_seed((time.tv_sec * 1000) + (time.tv_usec / 1000));
#endif

        for (m = 0; m < set->ss.nlocals; m++) {
            if (set->sc.step_act == 0) {  //generate positions at zero step
                mp->src->sl[m].xpos = (gint *)g_malloc(100000 * sizeof(gint)); //fixed maximum value
                mp->src->sl[m].ypos = (gint *)g_malloc(100000 * sizeof(gint)); //fixed maximum value
                mp->src->sl[m].zpos = (gint *)g_malloc(100000 * sizeof(gint)); //fixed maximum value
                mp->src->sl[m].xdir = (gdouble *)g_malloc(100000 * sizeof(gdouble)); //fixed maximum value
                mp->src->sl[m].ydir = (gdouble *)g_malloc(100000 * sizeof(gdouble)); //fixed maximum value
                mp->src->sl[m].zdir = (gdouble *)g_malloc(100000 * sizeof(gdouble)); //fixed maximum value
                mp->src->sl[m].source_amplitude = (gdouble *)g_malloc(100000 * sizeof(gdouble)); //fixed maximum value
                mp->src->sl[m].phase = (gdouble *)g_malloc(100000 * sizeof(gdouble)); //fixed maximum value
                mp->src->sl[m].lambda = (gdouble *)g_malloc(100000 * sizeof(gdouble)); //fixed maximum value
                mp->src->sl[m].avgfield = 0;

                mp->src->sl[m].npos = 0;
                //determine positions
                for (i = set->ss.locals[m].box_i0; i < set->ss.locals[m].box_in; i++) {
                    for (j = set->ss.locals[m].box_j0; j < set->ss.locals[m].box_jn; j++) {
                        for (k = set->ss.locals[m].box_k0; k < set->ss.locals[m].box_kn; k++) {
                            epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, j, k);

                            if (((set->ss.locals[m].layered_epsilon / EPSILON_0) == -1 || fabs(set->ss.locals[m].layered_epsilon - epsilon) < (EPSILON_0 / 100.0)) &&
                                g_rand_double_range(rnd, 0, 1) < set->ss.locals[m].density && mp->src->sl[m].npos < 100000) {
                                mp->src->sl[m].xpos[mp->src->sl[m].npos] = i;
                                mp->src->sl[m].ypos[mp->src->sl[m].npos] = j;
                                mp->src->sl[m].zpos[mp->src->sl[m].npos] = k;

                                xdir = g_rand_double_range(rnd, -1, 1);
                                ydir = g_rand_double_range(rnd, -1, 1);
                                zdir = g_rand_double_range(rnd, -1, 1);
                                mag = sqrt(xdir*xdir + ydir*ydir + zdir*zdir);

                                mp->src->sl[m].xdir[mp->src->sl[m].npos] = xdir / mag;
                                mp->src->sl[m].ydir[mp->src->sl[m].npos] = ydir / mag;
                                mp->src->sl[m].zdir[mp->src->sl[m].npos] = zdir / mag;

                                mp->src->sl[m].source_amplitude[mp->src->sl[m].npos] = 1;
                                mp->src->sl[m].phase[mp->src->sl[m].npos] = g_rand_double_range(rnd, 0, 50); //randomly chosen limit, fix this to wavelength 
                                mp->src->sl[m].lambda[mp->src->sl[m].npos] = get_gaussian_rnd(rnd, set->ss.locals[m].lambda_peak, set->ss.locals[m].lambda_width);

                                /* printf("source No %d at %d %d %d, (%g %g %g), amplitude %g phase %g, lambda %g\n", mp->src->sl[m].npos,
                                    mp->src->sl[m].xpos[mp->src->sl[m].npos], mp->src->sl[m].ypos[mp->src->sl[m].npos], mp->src->sl[m].zpos[mp->src->sl[m].npos],
                                    mp->src->sl[m].xdir[mp->src->sl[m].npos], mp->src->sl[m].ydir[mp->src->sl[m].npos], mp->src->sl[m].zdir[mp->src->sl[m].npos],
                                    mp->src->sl[m].amplitude[mp->src->sl[m].npos], mp->src->sl[m].phase[mp->src->sl[m].npos],
                                    mp->src->sl[m].lambda[mp->src->sl[m].npos]); */

                                mp->src->sl[m].npos++;
                            } // if
                        } // k
                    } // j
                } // i

                if (set->sc.verbose > 1) {
                    printf("%d point sources added...   ", mp->src->sl[m].npos);
                    fflush(stdout);
                }
            } // if (set->sc.step_act==0)

            for (sm = 0; sm < mp->src->sl[m].npos; sm++) {
                val = mp->src->sl[m].source_amplitude[sm] * sin(LIGHT_SPEED * 2 * G_PI*set->plan.dt*(gdouble)(set->sc.step_act + mp->src->sl[m].phase[sm]) / mp->src->sl[m].lambda[sm]);
                if (set->sc.step_act < 100)
                    val *= (0.01)*((gdouble)(set->sc.step_act));

                i = mp->src->sl[m].xpos[sm];
                j = mp->src->sl[m].ypos[sm];
                k = mp->src->sl[m].zpos[sm];

                if (set->ss.locals[m].source_mode == 0)
                    val *= set->ss.locals[m].strength;
                else if (set->ss.locals[m].source_mode == 1 || set->ss.locals[m].source_mode == 2) {
                    if (set->sc.step_act < set->ss.locals[m].startfrom) { //prevent feedback effects: evaluate field only prior to adding sources
                        mp->src->sl[m].avgfield = mp->src->sl[m].avgfield + (mp->d->ex->data[i][j][k] * mp->d->ex->data[i][j][k] + mp->d->ey->data[i][j][k] * mp->d->ey->data[i][j][k] + mp->d->ez->data[i][j][k] * mp->d->ez->data[i][j][k]) - mp->src->sl[m].avgfield / 100.0;
                    }

                    val *= (mp->src->sl[m].avgfield / 100.0) * set->ss.locals[m].strength;

                    if (set->ss.locals[m].source_mode == 2) {
                        sigma = sv_yee_data_get_sigma(mp->d, set, mp->mats, mp->nmat, i, j, k);
                        val *= sigma;
                    }
                    // if (sm==0) printf("\navav: %g  val %g\n", mp->src->sl[m].avgfield/100.0, (mp->src->sl[m].avgfield/100.0)*set->ss.locals[m].strength);                           
                } // else if

                if (set->sc.step_act > set->ss.locals[m].startfrom) {
                    mp->d->ex->data[i][j][k] += mp->src->sl[m].xdir[sm] * val;
                    mp->d->ey->data[i][j][k] += mp->src->sl[m].ydir[sm] * val;
                    mp->d->ez->data[i][j][k] += mp->src->sl[m].zdir[sm] * val;
                }
            } // sm                
        } // if (set->sc.step_act==0)
    } // m

    if (mp->src->tsf) {
        ecor = set->plan.dt / set->sp.dx;

        /*i0, i1*/
#pragma omp parallel default(shared) private(j, k, d, vh, gecor, epsilon, rphi, xradmult, yradmult)
#pragma omp for nowait

        for (j = mp->src->tsf->box_j0; j <= mp->src->tsf->box_jn; j++) {
            for (k = mp->src->tsf->box_k0; k <= mp->src->tsf->box_kn; k++) {
                //i0
                if (!set->ss.tsf.box_boundary_skipi0 && mat_too_close_check(mp, set, mp->src->tsf->box_i0, j, k)) {
                    if (set->ss.tsf.gaussian)
                        gecor = ecor * gaussmult(mp->src->tsf->box_i0 - 1, j, k, set->ss.tsf.gaussian_fxpos, set->ss.tsf.gaussian_fypos,
                                                 set->ss.tsf.gaussian_rx, set->ss.tsf.gaussian_ry);
                    else if (set->ss.tsf.fiber)
                        gecor = ecor * fibermult(mp->src->tsf->box_i0 - 1, j, k, set->ss.tsf.fiber_fxpos, set->ss.tsf.fiber_fypos,
                                                 set->ss.tsf.fiber_radius, set->ss.tsf.fiber_cutoff,
                                                 set->ss.tsf.fiber_epsilon_core, set->ss.tsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        gecor = ecor;

                    d = mp->src->tsf->corr * dcomp(mp->src->tsf->box_i0 - 1, j, k, set->sp.xres, set->sp.yres, set->sp.zres,
                                                   mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres);
                    vh = sv_dline_get_dval(mp->src->tsf->jp->h, d);
                    epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, mp->src->tsf->box_i0, j, k);

                    xradmult = yradmult = 1;
                    if (set->ss.tsf.radial == 1) { //z radial polarisation test
                        rphi = atan2(i - set->ss.tsf.radial_fxpos, j - set->ss.tsf.radial_fypos);
                        xradmult = sin(rphi);
                        yradmult = cos(rphi);
                        if (((i - set->ss.tsf.radial_fxpos)*(i - set->ss.tsf.radial_fxpos) + (j - set->ss.tsf.radial_fypos)*(j - set->ss.tsf.radial_fypos)) < 10)
                            xradmult = yradmult = 0;
                    }


                    if (j < mp->src->tsf->box_jn)
                        mp->d->ey->data[mp->src->tsf->box_i0][j][k] += yradmult * gecor / epsilon * ghz(vh, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                    if (k < mp->src->tsf->box_kn)
                        mp->d->ez->data[mp->src->tsf->box_i0][j][k] -= gecor / epsilon * ghy(vh, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                } // if (!set->ss.tsf.skipi0)

                //i1
                if (!set->ss.tsf.box_boundary_skipin && mat_too_close_check(mp, set, mp->src->tsf->box_in, j, k)) {
                    if (set->ss.tsf.gaussian)
                        gecor = ecor * gaussmult(mp->src->tsf->box_in, j, k, set->ss.tsf.gaussian_fxpos, set->ss.tsf.gaussian_fypos, set->ss.tsf.gaussian_rx, set->ss.tsf.gaussian_ry);
                    else if (set->ss.tsf.fiber)
                        gecor = ecor * fibermult(mp->src->tsf->box_in, j, k, set->ss.tsf.fiber_fxpos, set->ss.tsf.fiber_fypos,
                                                 set->ss.tsf.fiber_radius, set->ss.tsf.fiber_cutoff,
                                                 set->ss.tsf.fiber_epsilon_core, set->ss.tsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        gecor = ecor;

                    d = mp->src->tsf->corr * dcomp(mp->src->tsf->box_in, j, k, set->sp.xres, set->sp.yres, set->sp.zres,
                                                   mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres);
                    vh = sv_dline_get_dval(mp->src->tsf->jp->h, d);
                    epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, mp->src->tsf->box_in, j, k);

                    xradmult = yradmult = 1;
                    if (set->ss.tsf.radial == 1) { //z radial polarisation test
                        rphi = atan2(i - set->ss.tsf.radial_fxpos, j - set->ss.tsf.radial_fypos);
                        xradmult = sin(rphi);
                        yradmult = cos(rphi);
                        if (((i - set->ss.tsf.radial_fxpos)*(i - set->ss.tsf.radial_fxpos) + (j - set->ss.tsf.radial_fypos)*(j - set->ss.tsf.radial_fypos)) < 10)
                            xradmult = yradmult = 0;
                    }


                    if (j < mp->src->tsf->box_jn)
                        mp->d->ey->data[mp->src->tsf->box_in][j][k] -= yradmult * gecor / epsilon * ghz(vh, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                    if (k < mp->src->tsf->box_kn)
                        mp->d->ez->data[mp->src->tsf->box_in][j][k] += gecor / epsilon * ghy(vh, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                } // (!set->ss.tsf.skipin)
            } // k
        } // j

    /*j0, j1*/
#pragma omp parallel default(shared) private(i, k, d, vh, gecor, epsilon, rphi, xradmult, yradmult)
#pragma omp for nowait

        for (i = mp->src->tsf->box_i0; i <= mp->src->tsf->box_in; i++) {
            for (k = mp->src->tsf->box_k0; k <= mp->src->tsf->box_kn; k++) {
                //j0
                if (!set->ss.tsf.box_boundary_skipj0 && mat_too_close_check(mp, set, i, mp->src->tsf->box_j0, k)) {
                    if (set->ss.tsf.gaussian)
                        gecor = ecor * gaussmult(i, mp->src->tsf->box_j0 - 1, k, set->ss.tsf.gaussian_fxpos, set->ss.tsf.gaussian_fypos,
                                                 set->ss.tsf.gaussian_rx, set->ss.tsf.gaussian_ry);
                    else if (set->ss.tsf.fiber)
                        gecor = ecor * fibermult(i, mp->src->tsf->box_j0 - 1, k, set->ss.tsf.fiber_fxpos, set->ss.tsf.fiber_fypos,
                                                 set->ss.tsf.fiber_radius, set->ss.tsf.fiber_cutoff,
                                                 set->ss.tsf.fiber_epsilon_core, set->ss.tsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        gecor = ecor;

                    d = mp->src->tsf->corr * dcomp(i, mp->src->tsf->box_j0 - 1, k, set->sp.xres, set->sp.yres, set->sp.zres,
                                                   mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres);
                    vh = sv_dline_get_dval(mp->src->tsf->jp->h, d);
                    epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, mp->src->tsf->box_j0, k);

                    xradmult = yradmult = 1;
                    if (set->ss.tsf.radial == 1) { //z radial polarisation test
                        rphi = atan2(i - set->ss.tsf.radial_fxpos, j - set->ss.tsf.radial_fypos);
                        xradmult = sin(rphi);
                        yradmult = cos(rphi);
                        if (((i - set->ss.tsf.radial_fxpos)*(i - set->ss.tsf.radial_fxpos) + (j - set->ss.tsf.radial_fypos)*(j - set->ss.tsf.radial_fypos)) < 10)
                            xradmult = yradmult = 0;
                    }


                    // horni vec
                    if (i < mp->src->tsf->box_in)
                        mp->d->ex->data[i][mp->src->tsf->box_j0][k] -= xradmult * gecor / epsilon * ghz(vh, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                    if (k < mp->src->tsf->box_kn)
                        mp->d->ez->data[i][mp->src->tsf->box_j0][k] += gecor / epsilon * ghx(vh, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                }

                //j1
                if (!set->ss.tsf.box_boundary_skipjn && mat_too_close_check(mp, set, i, mp->src->tsf->box_jn, k)) {
                    if (set->ss.tsf.gaussian)
                        gecor = ecor*gaussmult(i, mp->src->tsf->box_jn, k, set->ss.tsf.gaussian_fxpos, set->ss.tsf.gaussian_fypos,
                                               set->ss.tsf.gaussian_rx, set->ss.tsf.gaussian_ry);
                    else if (set->ss.tsf.fiber)
                        gecor = ecor*fibermult(i, mp->src->tsf->box_jn, k, set->ss.tsf.fiber_fxpos, set->ss.tsf.fiber_fypos,
                                               set->ss.tsf.fiber_radius, set->ss.tsf.fiber_cutoff,
                                               set->ss.tsf.fiber_epsilon_core, set->ss.tsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        gecor = ecor;

                    d = mp->src->tsf->corr * dcomp(i, mp->src->tsf->box_jn, k, set->sp.xres, set->sp.yres, set->sp.zres,
                                                   mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres);
                    vh = sv_dline_get_dval(mp->src->tsf->jp->h, d);
                    epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, mp->src->tsf->box_jn, k);

                    xradmult = yradmult = 1;
                    if (set->ss.tsf.radial == 1) { //z radial polarisation test
                        rphi = atan2(i - set->ss.tsf.radial_fxpos, j - set->ss.tsf.radial_fypos);
                        xradmult = sin(rphi);
                        yradmult = cos(rphi);
                        if (((i - set->ss.tsf.radial_fxpos)*(i - set->ss.tsf.radial_fxpos) + (j - set->ss.tsf.radial_fypos)*(j - set->ss.tsf.radial_fypos)) < 10)
                            xradmult = yradmult = 0;
                    }

                    if (i < mp->src->tsf->box_in)
                        mp->d->ex->data[i][mp->src->tsf->box_jn][k] += xradmult * gecor / epsilon * ghz(vh, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                    if (k < mp->src->tsf->box_kn)
                        mp->d->ez->data[i][mp->src->tsf->box_jn][k] -= gecor / epsilon * ghx(vh, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                } // if (!set->ss.tsf.skipjn)
            } // k
        } // i

        /*k0, k1*/
#pragma omp parallel default(shared) private(i, j, d, vh, gecor, epsilon, rphi, xradmult, yradmult)
#pragma omp for nowait

        for (i = (mp->src->tsf->box_i0); i <= mp->src->tsf->box_in; i++) {
            for (j = (mp->src->tsf->box_j0); j <= mp->src->tsf->box_jn; j++) {
                //k0
                if (!set->ss.tsf.box_boundary_skipk0 && mat_too_close_check(mp, set, i, j, mp->src->tsf->box_k0)) {
                    if (set->ss.tsf.gaussian)
                        gecor = ecor * gaussmult(i, j, mp->src->tsf->box_k0 - 1, set->ss.tsf.gaussian_fxpos, set->ss.tsf.gaussian_fypos,
                                                 set->ss.tsf.gaussian_rx, set->ss.tsf.gaussian_ry);
                    else if (set->ss.tsf.fiber)
                        gecor = ecor * fibermult(i, j, mp->src->tsf->box_k0 - 1, set->ss.tsf.fiber_fxpos, set->ss.tsf.fiber_fypos,
                                                 set->ss.tsf.fiber_radius, set->ss.tsf.fiber_cutoff,
                                                 set->ss.tsf.fiber_epsilon_core, set->ss.tsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        gecor = ecor;

                    d = mp->src->tsf->corr * dcomp(i, j, mp->src->tsf->box_k0 - 1, set->sp.xres, set->sp.yres, set->sp.zres,
                                                   mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres);
                    vh = sv_dline_get_dval(mp->src->tsf->jp->h, d);
                    epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, j, mp->src->tsf->box_k0);

                    //nahore a dole 

                    //AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa

                    xradmult = yradmult = 1;
                    if (set->ss.tsf.radial == 1) { //z radial polarisation test
                        rphi = atan2(i - set->ss.tsf.radial_fxpos, j - set->ss.tsf.radial_fypos);
                        xradmult = sin(rphi);
                        yradmult = cos(rphi);
                        if (((i - set->ss.tsf.radial_fxpos)*(i - set->ss.tsf.radial_fxpos) + (j - set->ss.tsf.radial_fypos)*(j - set->ss.tsf.radial_fypos)) < 10)
                            xradmult = yradmult = 0;
                    }

                    if (i < mp->src->tsf->box_in && !(sv_yee_data_is_pec(mp->d, set, mp->mats, mp->nmat, i, j, mp->src->tsf->box_k0)))
                        mp->d->ex->data[i][j][mp->src->tsf->box_k0] += xradmult*gecor / epsilon * ghy(vh, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                    if (j < mp->src->tsf->box_jn && !(sv_yee_data_is_pec(mp->d, set, mp->mats, mp->nmat, i, j, mp->src->tsf->box_k0)))
                        mp->d->ey->data[i][j][mp->src->tsf->box_k0] -= yradmult*gecor / epsilon * ghx(vh, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                }
                //k1
                if (!set->ss.tsf.box_boundary_skipkn && mat_too_close_check(mp, set, i, j, mp->src->tsf->box_kn)) {
                    if (set->ss.tsf.gaussian)
                        gecor = ecor * gaussmult(i, j, mp->src->tsf->box_kn, set->ss.tsf.gaussian_fxpos, set->ss.tsf.gaussian_fypos,
                                                 set->ss.tsf.gaussian_rx, set->ss.tsf.gaussian_ry);
                    else if (set->ss.tsf.fiber)
                        gecor = ecor * fibermult(i, j, mp->src->tsf->box_kn, set->ss.tsf.fiber_fxpos, set->ss.tsf.fiber_fypos,
                                                 set->ss.tsf.fiber_radius, set->ss.tsf.fiber_cutoff,
                                                 set->ss.tsf.fiber_epsilon_core, set->ss.tsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        gecor = ecor;

                    d = mp->src->tsf->corr * dcomp(i, j, mp->src->tsf->box_kn, set->sp.xres, set->sp.yres, set->sp.zres,
                                                   mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres);
                    vh = sv_dline_get_dval(mp->src->tsf->jp->h, d);
                    epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, j, mp->src->tsf->box_kn);

                    xradmult = yradmult = 1;
                    if (set->ss.tsf.radial == 1) { //z radial polarisation test
                        rphi = atan2(i - set->ss.tsf.radial_fxpos, j - set->ss.tsf.radial_fypos);
                        xradmult = sin(rphi);
                        yradmult = cos(rphi);
                        if (((i - set->ss.tsf.radial_fxpos)*(i - set->ss.tsf.radial_fxpos) + (j - set->ss.tsf.radial_fypos)*(j - set->ss.tsf.radial_fypos)) < 10)
                            xradmult = yradmult = 0;
                    }

                    if (i < mp->src->tsf->box_in)
                        mp->d->ex->data[i][j][mp->src->tsf->box_kn] -= xradmult*gecor / epsilon * ghy(vh, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                    if (j < mp->src->tsf->box_jn)
                        mp->d->ey->data[i][j][mp->src->tsf->box_kn] += yradmult*gecor / epsilon * ghx(vh, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                } // if (!set->ss.tsf.skipkn)
            } // j
        } // i
    } // if (mp->src->tsf)

    if (mp->src->tsff) {
        ecor = set->plan.dt / set->sp.dx;

        for (i1 = 0; i1 < mp->src->tsff->focused_nip; i1++) {
            for (im = 0; im < mp->src->tsff->focused_mip; im++) {
                //   if (i1!=0 || im!=0) continue;
                //if (!(i1==0 && im==0)) continue;

                theta = mp->src->tsff->ia_theta[i1][im];
                phi = mp->src->tsff->ia_phi[i1][im];
                psi = mp->src->tsff->ia_psi[i1][im];
                corr = mp->src->tsff->corr[i1][im];
                an = mp->src->tsff->an[i1][im];
                bm = mp->src->tsff->bm[i1][im];

                ctheta = cos(theta); stheta = sin(theta);
                cphi = cos(phi); sphi = sin(phi);
                cpsi = cos(psi); spsi = sin(psi);

                acorr = fs_amplitude(mp->src->tsff->focused_fdist*set->sp.dx, an, bm, theta);

                //printf("src run: %g %g %g\n", theta, phi, acorr);

                /*i0, i1*/
#pragma omp parallel default(shared) private(j, k, d, vh, epsilon)
#pragma omp for nowait

                for (j = mp->src->tsff->box_j0; j <= mp->src->tsff->box_jn; j++) {
                    for (k = mp->src->tsff->box_k0; k <= mp->src->tsff->box_kn; k++) {
                        /*i0*/
                        if (!set->ss.tsff.box_boundary_skipi0) {
                            d = corr * srdcomp(mp->src->tsff->box_i0 - 1, j, k, set->sp.xres, set->sp.yres, set->sp.zres,
                                               stheta, sphi, ctheta, cphi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres,
                                               set->ss.tsff.xshift, set->ss.tsff.yshift, set->ss.tsff.zshift);

                            vh = acorr * sv_dline_get_dval(mp->src->tsff->jp[i1][im]->h, d);
                            epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, mp->src->tsff->box_i0, j, k);

                            if (j < mp->src->tsff->box_jn)
                                mp->d->ey->data[mp->src->tsff->box_i0][j][k] += ecor / epsilon * vh * (-cpsi*stheta);
                            // *sghz(vh, stheta, sphi, spsi, ctheta, cphi, cpsi);
                            if (k < mp->src->tsff->box_kn)
                                mp->d->ez->data[mp->src->tsff->box_i0][j][k] -= ecor / epsilon * vh * (-spsi*cphi + cpsi*ctheta*sphi);
                            // *sghy(vh, stheta, sphi, spsi, ctheta, cphi, cpsi);
                        }

                        //i1
                        if (!set->ss.tsff.box_boundary_skipin) {
                            d = corr * srdcomp(mp->src->tsff->box_in, j, k, set->sp.xres, set->sp.yres, set->sp.zres,
                                               stheta, sphi, ctheta, cphi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres,
                                               set->ss.tsff.xshift, set->ss.tsff.yshift, set->ss.tsff.zshift);
                            vh = acorr*sv_dline_get_dval(mp->src->tsff->jp[i1][im]->h, d);
                            epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, mp->src->tsff->box_in, j, k);

                            if (j < mp->src->tsff->box_jn)
                                mp->d->ey->data[mp->src->tsff->box_in][j][k] -= ecor / epsilon * vh * (-cpsi*stheta);
                            //*sghz(vh, stheta, sphi, spsi, ctheta, cphi, cpsi);
                            if (k < mp->src->tsff->box_kn)
                                mp->d->ez->data[mp->src->tsff->box_in][j][k] += ecor / epsilon * vh * (-spsi*cphi + cpsi*ctheta*sphi);
                            //*sghy(vh, stheta, sphi, spsi, ctheta, cphi, cpsi);
                        }
                    } // k
                } // i

                /*j0, j1*/
#pragma omp parallel default(shared) private(i, k, d, vh, epsilon)
#pragma omp for nowait

                for (i = (mp->src->tsff->box_i0); i <= mp->src->tsff->box_in; i++) {
                    for (k = (mp->src->tsff->box_k0); k <= mp->src->tsff->box_kn; k++) {
                        /*j0*/
                        if (!set->ss.tsff.box_boundary_skipj0) {
                            d = corr * srdcomp(i, mp->src->tsff->box_j0 - 1, k, set->sp.xres, set->sp.yres, set->sp.zres,
                                               stheta, sphi, ctheta, cphi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres,
                                               set->ss.tsff.xshift, set->ss.tsff.yshift, set->ss.tsff.zshift);
                            vh = acorr * sv_dline_get_dval(mp->src->tsff->jp[i1][im]->h, d);
                            epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, mp->src->tsff->box_j0, k);

                            //horni vec
                            if (i < mp->src->tsff->box_in)
                                mp->d->ex->data[i][mp->src->tsff->box_j0][k] -= ecor / epsilon * vh * (-cpsi*stheta);
                            //sghz(vh, stheta, sphi, spsi, ctheta, cphi, cpsi);
                            if (k < mp->src->tsff->box_kn)
                                mp->d->ez->data[i][mp->src->tsff->box_j0][k] += ecor / epsilon * vh * (spsi*sphi + cpsi*ctheta*cphi);
                            //sghx(vh, stheta, sphi, spsi, ctheta, cphi, cpsi);
                        }

                        //j1
                        if (!set->ss.tsff.box_boundary_skipjn) {
                            d = corr * srdcomp(i, mp->src->tsff->box_jn, k, set->sp.xres, set->sp.yres, set->sp.zres,
                                               stheta, sphi, ctheta, cphi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres,
                                               set->ss.tsff.xshift, set->ss.tsff.yshift, set->ss.tsff.zshift);
                            vh = acorr * sv_dline_get_dval(mp->src->tsff->jp[i1][im]->h, d);
                            epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, mp->src->tsff->box_jn, k);

                            if (i < mp->src->tsff->box_in)
                                mp->d->ex->data[i][mp->src->tsff->box_jn][k] += ecor / epsilon * vh * (-cpsi*stheta);
                            //sghz(vh, stheta, sphi, spsi, ctheta, cphi, cpsi);
                            if (k < mp->src->tsff->box_kn)
                                mp->d->ez->data[i][mp->src->tsff->box_jn][k] -= ecor / epsilon * vh * (spsi*sphi + cpsi*ctheta*cphi);
                            //sghx(vh, stheta, sphi, spsi, ctheta, cphi, cpsi);
                        }
                    } // k
                } // i

                /*k0, k1*/
#pragma omp parallel default(shared) private(i, j, d, vh, epsilon)
#pragma omp for nowait

                for (i = mp->src->tsff->box_i0; i <= mp->src->tsff->box_in; i++) {
                    for (j = mp->src->tsff->box_j0; j <= mp->src->tsff->box_jn; j++) {
                        /*k0*/
                        if (!set->ss.tsff.box_boundary_skipk0) {
                            d = corr * srdcomp(i, j, mp->src->tsff->box_k0 - 1, set->sp.xres, set->sp.yres, set->sp.zres,
                                               stheta, sphi, ctheta, cphi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres,
                                               set->ss.tsff.xshift, set->ss.tsff.yshift, set->ss.tsff.zshift);
                            vh = acorr * sv_dline_get_dval(mp->src->tsff->jp[i1][im]->h, d);
                            epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, j, mp->src->tsff->box_k0);

                            //nahore a dole 
                            if (i < mp->src->tsff->box_in)
                                mp->d->ex->data[i][j][mp->src->tsff->box_k0] += ecor / epsilon * vh * (-spsi*cphi + cpsi*ctheta*sphi);
                            //*sghy(vh, stheta, sphi, spsi, ctheta, cphi, cpsi);
                            if (j < mp->src->tsff->box_jn)
                                mp->d->ey->data[i][j][mp->src->tsff->box_k0] -= ecor / epsilon * vh * (spsi*sphi + cpsi*ctheta*cphi);
                            //sghx(vh, stheta, sphi, spsi, ctheta, cphi, cpsi);
                        }
                        //k1
                        if (!set->ss.tsff.box_boundary_skipkn) {
                            d = corr * srdcomp(i, j, mp->src->tsff->box_kn, set->sp.xres, set->sp.yres, set->sp.zres,
                                               stheta, sphi, ctheta, cphi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres,
                                               set->ss.tsff.xshift, set->ss.tsff.yshift, set->ss.tsff.zshift);
                            vh = acorr * sv_dline_get_dval(mp->src->tsff->jp[i1][im]->h, d);
                            epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, j, mp->src->tsff->box_kn);

                            if (i < mp->src->tsff->box_in)
                                mp->d->ex->data[i][j][mp->src->tsff->box_kn] -= ecor / epsilon * vh * (-spsi*cphi + cpsi*ctheta*sphi);
                            //*sghy(vh, stheta, sphi, spsi, ctheta, cphi, cpsi);
                            if (j < mp->src->tsff->box_jn)
                                mp->d->ey->data[i][j][mp->src->tsff->box_kn] += ecor / epsilon * vh * (spsi*sphi + cpsi*ctheta*cphi);
                            //*sghx(vh, stheta, sphi, spsi, ctheta, cphi, cpsi);
                        }
                    } // j
                } // i
            } // im
        } // in
    } //  if (mp->src->tsf)

    if (mp->src->ltsf) {
        ecor = set->plan.dt / set->sp.dx;

        /*i0, i1*/
#pragma omp parallel default(shared) private(j, k, zpos, tp, atp, hy, hz, gecor, vhh, vhe, epsilon)
#pragma omp for nowait

        for (j = mp->src->ltsf->box_j0; j <= mp->src->ltsf->box_jn; j++) {
            for (k = mp->src->ltsf->box_k0; k <= mp->src->ltsf->box_kn; k++) {
                /*i0*/
                if (!set->ss.ltsf.box_boundary_skipi0) {
                    if (set->ss.ltsf.gaussian)
                        gecor = ecor * gaussmult(mp->src->ltsf->box_i0 - 1, j, k, set->ss.ltsf.gaussian_fxpos, set->ss.ltsf.gaussian_fypos,
                                                 set->ss.ltsf.gaussian_rx, set->ss.ltsf.gaussian_ry);
                    else if (set->ss.ltsf.fiber)
                        gecor = ecor * fibermult(mp->src->ltsf->box_i0 - 1, j, k, set->ss.ltsf.fiber_fxpos, set->ss.ltsf.fiber_fypos,
                                                 set->ss.ltsf.fiber_radius, set->ss.ltsf.fiber_cutoff,
                                                 set->ss.ltsf.fiber_epsilon_core, set->ss.ltsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        gecor = ecor;

                    getcoefs(mp->src->ltsf->box_i0 - 1, j, k, mp->src->ltsf->ia_theta, mp->src->ltsf->ia_phi,
                             set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                             mp->src->ltsf->box_i0, mp->src->ltsf->box_j0, mp->src->ltsf->box_k0,
                             mp->src->ltsf->box_in, mp->src->ltsf->box_jn, mp->src->ltsf->box_kn);

                    vhh = (1.0 - atp) * mp->src->ltsf->jp->dhh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhh->data[(tp + 1)*set->sp.zres + zpos];
                    vhe = (1.0 - atp) * mp->src->ltsf->jp->dhe->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhe->data[(tp + 1)*set->sp.zres + zpos];
                    hz = (1.0 - atp) * mp->src->ltsf->jp->dhz->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhz->data[(tp + 1)*set->sp.zres + zpos];

                    hy = vhh*sin(mp->src->ltsf->ia_phi) - vhe*cos(mp->src->ltsf->ia_phi);
                    hy = -hy;

                    epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, mp->src->ltsf->box_i0, j, k);

                    if (j < mp->src->ltsf->box_jn)
                        mp->d->ey->data[mp->src->ltsf->box_i0][j][k] += gecor / epsilon * hz;
                    if (k < mp->src->ltsf->box_kn)
                        mp->d->ez->data[mp->src->ltsf->box_i0][j][k] -= gecor / epsilon * hy;

                } // if (!set->ss.ltsf.skipi0)

                        //i1
                if (!set->ss.ltsf.box_boundary_skipin) {
                    if (set->ss.ltsf.gaussian)
                        gecor = ecor * gaussmult(mp->src->ltsf->box_in, j, k, set->ss.ltsf.gaussian_fxpos, set->ss.ltsf.gaussian_fypos,
                                                 set->ss.ltsf.gaussian_rx, set->ss.ltsf.gaussian_ry);
                    else if (set->ss.ltsf.fiber)
                        gecor = ecor*fibermult(mp->src->ltsf->box_in, j, k, set->ss.ltsf.fiber_fxpos, set->ss.ltsf.fiber_fypos,
                                               set->ss.ltsf.fiber_radius, set->ss.ltsf.fiber_cutoff,
                                               set->ss.ltsf.fiber_epsilon_core, set->ss.ltsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        gecor = ecor;

                    getcoefs(mp->src->ltsf->box_in, j, k, mp->src->ltsf->ia_theta, mp->src->ltsf->ia_phi,
                             set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                             mp->src->ltsf->box_i0, mp->src->ltsf->box_j0, mp->src->ltsf->box_k0,
                             mp->src->ltsf->box_in, mp->src->ltsf->box_jn, mp->src->ltsf->box_kn);

                    vhh = (1.0 - atp) * mp->src->ltsf->jp->dhh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhh->data[(tp + 1)*set->sp.zres + zpos];
                    vhe = (1.0 - atp) * mp->src->ltsf->jp->dhe->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhe->data[(tp + 1)*set->sp.zres + zpos];
                    hz = (1.0 - atp) * mp->src->ltsf->jp->dhz->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhz->data[(tp + 1)*set->sp.zres + zpos];

                    hy = vhh*sin(mp->src->ltsf->ia_phi) - vhe*cos(mp->src->ltsf->ia_phi);
                    hy = -hy;
                    epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, mp->src->ltsf->box_in, j, k);

                    if (j < mp->src->ltsf->box_jn)
                        mp->d->ey->data[mp->src->ltsf->box_in][j][k] -= gecor / epsilon * hz;
                    if (k < mp->src->ltsf->box_kn)
                        mp->d->ez->data[mp->src->ltsf->box_in][j][k] += gecor / epsilon * hy;
                } //  if (!set->ss.ltsf.skipin)
            } // k
        } // j

        /*j0, j1*/
#pragma omp parallel default(shared) private(i, k, zpos, tp, atp, hx, hz, gecor, vhh, vhe, epsilon)
#pragma omp for nowait

        for (i = mp->src->ltsf->box_i0; i <= mp->src->ltsf->box_in; i++) {
            for (k = mp->src->ltsf->box_k0; k <= mp->src->ltsf->box_kn; k++) {
                /*j0*/
                if (!set->ss.ltsf.box_boundary_skipj0) {
                    if (set->ss.ltsf.gaussian)
                        gecor = ecor*gaussmult(i, mp->src->ltsf->box_j0 - 1, k, set->ss.ltsf.gaussian_fxpos, set->ss.ltsf.gaussian_fypos,
                                               set->ss.ltsf.gaussian_rx, set->ss.ltsf.gaussian_ry);
                    else if (set->ss.ltsf.fiber)
                        gecor = ecor*fibermult(i, mp->src->ltsf->box_j0 - 1, k, set->ss.ltsf.fiber_fxpos, set->ss.ltsf.fiber_fypos,
                                               set->ss.ltsf.fiber_radius, set->ss.ltsf.fiber_cutoff,
                                               set->ss.ltsf.fiber_epsilon_core, set->ss.ltsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        gecor = ecor;

                    getcoefs(i, mp->src->ltsf->box_j0 - 1, k, mp->src->ltsf->ia_theta, mp->src->ltsf->ia_phi,
                             set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                             mp->src->ltsf->box_i0, mp->src->ltsf->box_j0, mp->src->ltsf->box_k0,
                             mp->src->ltsf->box_in, mp->src->ltsf->box_jn, mp->src->ltsf->box_kn);

                    vhh = (1.0 - atp) * mp->src->ltsf->jp->dhh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhh->data[(tp + 1)*set->sp.zres + zpos];
                    vhe = (1.0 - atp) * mp->src->ltsf->jp->dhe->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhe->data[(tp + 1)*set->sp.zres + zpos];
                    hz = (1.0 - atp) * mp->src->ltsf->jp->dhz->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhz->data[(tp + 1)*set->sp.zres + zpos];

                    hx = vhh*cos(mp->src->ltsf->ia_phi) + vhe*sin(mp->src->ltsf->ia_phi);
                    hx = -hx;
                    epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, mp->src->ltsf->box_j0, k);

                    if (i < mp->src->ltsf->box_in)
                        mp->d->ex->data[i][mp->src->ltsf->box_j0][k] -= gecor / epsilon * hz;
                    if (k < mp->src->ltsf->box_kn)
                        mp->d->ez->data[i][mp->src->ltsf->box_j0][k] += gecor / epsilon * hx;
                } //  if (!set->ss.ltsf.skipj0)

                //j1
                if (!set->ss.ltsf.box_boundary_skipjn) {
                    if (set->ss.ltsf.gaussian)
                        gecor = ecor * gaussmult(i, mp->src->ltsf->box_jn, k, set->ss.ltsf.gaussian_fxpos, set->ss.ltsf.gaussian_fypos,
                                                 set->ss.ltsf.gaussian_rx, set->ss.ltsf.gaussian_ry);
                    else if (set->ss.ltsf.fiber)
                        gecor = ecor*fibermult(i, mp->src->ltsf->box_jn, k, set->ss.ltsf.fiber_fxpos, set->ss.ltsf.fiber_fypos,
                                               set->ss.ltsf.fiber_radius, set->ss.ltsf.fiber_cutoff,
                                               set->ss.ltsf.fiber_epsilon_core, set->ss.ltsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        gecor = ecor;

                    getcoefs(i, mp->src->ltsf->box_jn, k, mp->src->ltsf->ia_theta, mp->src->ltsf->ia_phi,
                             set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                             mp->src->ltsf->box_i0, mp->src->ltsf->box_j0, mp->src->ltsf->box_k0,
                             mp->src->ltsf->box_in, mp->src->ltsf->box_jn, mp->src->ltsf->box_kn);

                    vhh = (1.0 - atp) * mp->src->ltsf->jp->dhh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhh->data[(tp + 1)*set->sp.zres + zpos];
                    vhe = (1.0 - atp) * mp->src->ltsf->jp->dhe->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhe->data[(tp + 1)*set->sp.zres + zpos];
                    hz = (1.0 - atp) * mp->src->ltsf->jp->dhz->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhz->data[(tp + 1)*set->sp.zres + zpos];

                    hx = vhh*cos(mp->src->ltsf->ia_phi) + vhe*sin(mp->src->ltsf->ia_phi);
                    hx = -hx;
                    epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, mp->src->ltsf->box_jn, k);

                    if (i < mp->src->ltsf->box_in)
                        mp->d->ex->data[i][mp->src->ltsf->box_jn][k] += gecor / epsilon * hz;
                    if (k < mp->src->ltsf->box_kn)
                        mp->d->ez->data[i][mp->src->ltsf->box_jn][k] -= gecor / epsilon * hx;
                } //  if (!set->ss.ltsf.skipjn)
            } // k
        } // i

        /*k0, k1*/
#pragma omp parallel default(shared) private(i, j, zpos, tp, atp, hx, hy, gecor, vhh, vhe, epsilon)
#pragma omp for nowait

        for (i = mp->src->ltsf->box_i0; i <= mp->src->ltsf->box_in; i++) {
            for (j = mp->src->ltsf->box_j0; j <= mp->src->ltsf->box_jn; j++) {
                /*k0*/
                if (!set->ss.ltsf.box_boundary_skipk0  && mat_too_close_check(mp, set, i, j, mp->src->ltsf->box_k0)) {
                    if (set->ss.ltsf.gaussian)
                        gecor = ecor * gaussmult(i, j, mp->src->ltsf->box_k0 - 1, set->ss.ltsf.gaussian_fxpos, set->ss.ltsf.gaussian_fypos,
                                                 set->ss.ltsf.gaussian_rx, set->ss.ltsf.gaussian_ry);
                    else if (set->ss.ltsf.fiber)
                        gecor = ecor * fibermult(i, j, mp->src->ltsf->box_k0 - 1, set->ss.ltsf.fiber_fxpos, set->ss.ltsf.fiber_fypos,
                                                 set->ss.ltsf.fiber_radius, set->ss.ltsf.fiber_cutoff,
                                                 set->ss.ltsf.fiber_epsilon_core, set->ss.ltsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        gecor = ecor;

                    getcoefs(i, j, mp->src->ltsf->box_k0 - 1, mp->src->ltsf->ia_theta, mp->src->ltsf->ia_phi,
                             set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                             mp->src->ltsf->box_i0, mp->src->ltsf->box_j0, mp->src->ltsf->box_k0,
                             mp->src->ltsf->box_in, mp->src->ltsf->box_jn, mp->src->ltsf->box_kn);

                    vhh = (1.0 - atp) * mp->src->ltsf->jp->dhh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhh->data[(tp + 1)*set->sp.zres + zpos];
                    vhe = (1.0 - atp) * mp->src->ltsf->jp->dhe->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhe->data[(tp + 1)*set->sp.zres + zpos];

                    hx = vhh*cos(mp->src->ltsf->ia_phi) + vhe*sin(mp->src->ltsf->ia_phi);
                    hy = vhh*sin(mp->src->ltsf->ia_phi) - vhe*cos(mp->src->ltsf->ia_phi);
                    hx = -hx;
                    hy = -hy;
                    epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, j, mp->src->ltsf->box_k0);

                    if (i < mp->src->ltsf->box_in)
                        mp->d->ex->data[i][j][mp->src->ltsf->box_k0] += gecor / epsilon * hy;
                    if (j < mp->src->ltsf->box_jn)
                        mp->d->ey->data[i][j][mp->src->ltsf->box_k0] -= gecor / epsilon * hx;
                } //  if (!set->ss.ltsf.skipk0)
                //k1
                //

                if (!set->ss.ltsf.box_boundary_skipkn) {
                    if (set->ss.ltsf.gaussian)
                        gecor = ecor * gaussmult(i, j, mp->src->ltsf->box_kn, set->ss.ltsf.gaussian_fxpos, set->ss.ltsf.gaussian_fypos,
                                                 set->ss.ltsf.gaussian_rx, set->ss.ltsf.gaussian_ry);
                    else if (set->ss.ltsf.fiber)
                        gecor = ecor * fibermult(i, j, mp->src->ltsf->box_kn, set->ss.ltsf.fiber_fxpos, set->ss.ltsf.fiber_fypos,
                                                 set->ss.ltsf.fiber_radius, set->ss.ltsf.fiber_cutoff,
                                                 set->ss.ltsf.fiber_epsilon_core, set->ss.ltsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        gecor = ecor;

                    getcoefs(i, j, mp->src->ltsf->box_kn, mp->src->ltsf->ia_theta, mp->src->ltsf->ia_phi,
                             set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                             mp->src->ltsf->box_i0, mp->src->ltsf->box_j0, mp->src->ltsf->box_k0,
                             mp->src->ltsf->box_in, mp->src->ltsf->box_jn, mp->src->ltsf->box_kn);

                    vhh = (1.0 - atp) * mp->src->ltsf->jp->dhh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhh->data[(tp + 1)*set->sp.zres + zpos];
                    vhe = (1.0 - atp) * mp->src->ltsf->jp->dhe->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dhe->data[(tp + 1)*set->sp.zres + zpos];

                    hx = vhh*cos(mp->src->ltsf->ia_phi) + vhe*sin(mp->src->ltsf->ia_phi);
                    hy = vhh*sin(mp->src->ltsf->ia_phi) - vhe*cos(mp->src->ltsf->ia_phi);
                    hx = -hx;
                    hy = -hy;
                    epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, j, mp->src->ltsf->box_kn);

                    if (i < mp->src->ltsf->box_in)
                        mp->d->ex->data[i][j][mp->src->ltsf->box_kn] -= gecor / epsilon * hy;
                    if (j < mp->src->ltsf->box_jn)
                        mp->d->ey->data[i][j][mp->src->ltsf->box_kn] += gecor / epsilon * hx;
                } //  if (!set->ss.ltsf.skipkn)
            } // j
        } // i
    } // if (mp->src->ltsf)

    if (mp->src->ltsff) {
        ecor = set->plan.dt / set->sp.dx;

        for (i1 = 0; i1 < mp->src->ltsff->focused_nip; i1++) {
            for (im = 0; im < mp->src->ltsff->focused_mip; im++) {
                theta = mp->src->ltsff->ia_theta[i1][im];
                phi = mp->src->ltsff->ia_phi[i1][im];
                psi = mp->src->ltsff->ia_psi[i1][im];
                corr = mp->src->ltsff->corr[i1][im];
                an = mp->src->ltsff->an[i1][im];
                bm = mp->src->ltsff->bm[i1][im];

                ctheta = cos(theta); stheta = sin(theta);
                cphi = cos(phi); sphi = sin(phi);
                cpsi = cos(psi); spsi = sin(psi);

                acorr = fs_amplitude(mp->src->ltsff->focused_fdist * set->sp.dx, an, bm, theta);

                /*i0, i1*/
#pragma omp parallel default(shared) private(j, k, zpos, tp, atp, hy, hz, vhh, vhe, epsilon)
#pragma omp for nowait
                for (j = mp->src->ltsff->box_j0; j <= mp->src->ltsff->box_jn; j++) {
                    for (k = mp->src->ltsff->box_k0; k <= mp->src->ltsff->box_kn; k++) {
                        /*i0*/
                        if (!set->ss.ltsff.box_boundary_skipi0) {
                            sgetcoefs(mp->src->ltsff->box_i0 - 1, j, k, mp->src->ltsff->ia_theta[i1][im], mp->src->ltsff->ia_phi[i1][im],
                                      set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                                      mp->src->ltsff->box_i0, mp->src->ltsff->box_j0, mp->src->ltsff->box_k0, mp->src->ltsff->box_in,
                                      mp->src->ltsff->box_jn, mp->src->ltsff->box_kn,
                                      set->sp.xres, set->sp.yres, set->sp.zres);

                            vhh = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhh->data[(tp + 1)*set->sp.zres + zpos];
                            vhe = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhe->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhe->data[(tp + 1)*set->sp.zres + zpos];
                            hz = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhz->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhz->data[(tp + 1)*set->sp.zres + zpos];

                            hy = vhh*sin(mp->src->ltsff->ia_phi[i1][im]) - vhe*cos(mp->src->ltsff->ia_phi[i1][im]);
                            hy = -hy;
                            epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, mp->src->ltsff->box_i0, j, k);

                            if (j < mp->src->ltsff->box_jn)
                                mp->d->ey->data[mp->src->ltsff->box_i0][j][k] += acorr * ecor / epsilon * hz;
                            if (k < mp->src->ltsff->box_kn)
                                mp->d->ez->data[mp->src->ltsff->box_i0][j][k] -= acorr * ecor / epsilon * hy;
                        } // (!set->ss.ltsff.skipi0)

                        //i1
                        if (!set->ss.ltsff.box_boundary_skipin) {
                            sgetcoefs(mp->src->ltsff->box_in, j, k, mp->src->ltsff->ia_theta[i1][im], mp->src->ltsff->ia_phi[i1][im],
                                      set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                                      mp->src->ltsff->box_i0, mp->src->ltsff->box_j0, mp->src->ltsff->box_k0, mp->src->ltsff->box_in,
                                      mp->src->ltsff->box_jn, mp->src->ltsff->box_kn,
                                      set->sp.xres, set->sp.yres, set->sp.zres);

                            vhh = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhh->data[(tp + 1)*set->sp.zres + zpos];
                            vhe = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhe->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhe->data[(tp + 1)*set->sp.zres + zpos];
                            hz = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhz->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhz->data[(tp + 1)*set->sp.zres + zpos];

                            hy = vhh*sin(mp->src->ltsff->ia_phi[i1][im]) - vhe*cos(mp->src->ltsff->ia_phi[i1][im]);
                            hy = -hy;
                            epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, mp->src->ltsff->box_in, j, k);

                            if (j < mp->src->ltsff->box_jn)
                                mp->d->ey->data[mp->src->ltsff->box_in][j][k] -= acorr * ecor / epsilon * hz;
                            if (k < mp->src->ltsff->box_kn)
                                mp->d->ez->data[mp->src->ltsff->box_in][j][k] += acorr * ecor / epsilon * hy;
                        } // if (!set->ss.ltsff.skipin)
                    } // k
                } // j

                /*j0, j1*/
#pragma omp parallel default(shared) private(i, k, zpos, tp, atp, hx, hz, vhh, vhe, epsilon)
#pragma omp for nowait

                for (i = mp->src->ltsff->box_i0; i <= mp->src->ltsff->box_in; i++) {
                    for (k = mp->src->ltsff->box_k0; k <= mp->src->ltsff->box_kn; k++) {
                        /*j0*/
                        if (!set->ss.ltsff.box_boundary_skipj0) {
                            sgetcoefs(i, mp->src->ltsff->box_j0 - 1, k, mp->src->ltsff->ia_theta[i1][im], mp->src->ltsff->ia_phi[i1][im],
                                      set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                                      mp->src->ltsff->box_i0, mp->src->ltsff->box_j0, mp->src->ltsff->box_k0,
                                      mp->src->ltsff->box_in, mp->src->ltsff->box_jn, mp->src->ltsff->box_kn,
                                      set->sp.xres, set->sp.yres, set->sp.zres);

                            vhh = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhh->data[(tp + 1)*set->sp.zres + zpos];
                            vhe = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhe->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhe->data[(tp + 1)*set->sp.zres + zpos];
                            hz = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhz->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhz->data[(tp + 1)*set->sp.zres + zpos];

                            hx = vhh*cos(mp->src->ltsff->ia_phi[i1][im]) + vhe*sin(mp->src->ltsff->ia_phi[i1][im]);
                            hx = -hx;
                            epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, mp->src->ltsff->box_j0, k);

                            if (i < mp->src->ltsff->box_in)
                                mp->d->ex->data[i][mp->src->ltsff->box_j0][k] -= acorr * ecor / epsilon * hz;
                            if (k < mp->src->ltsff->box_kn)
                                mp->d->ez->data[i][mp->src->ltsff->box_j0][k] += acorr * ecor / epsilon * hx;
                        } // if (!set->ss.ltsff.skipj0)

                        //j1
                        if (!set->ss.ltsff.box_boundary_skipjn) {
                            sgetcoefs(i, mp->src->ltsff->box_jn, k, mp->src->ltsff->ia_theta[i1][im], mp->src->ltsff->ia_phi[i1][im],
                                      set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                                      mp->src->ltsff->box_i0, mp->src->ltsff->box_j0, mp->src->ltsff->box_k0,
                                      mp->src->ltsff->box_in, mp->src->ltsff->box_jn, mp->src->ltsff->box_kn,
                                      set->sp.xres, set->sp.yres, set->sp.zres);

                            vhh = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhh->data[(tp + 1)*set->sp.zres + zpos];
                            vhe = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhe->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhe->data[(tp + 1)*set->sp.zres + zpos];
                            hz = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhz->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhz->data[(tp + 1)*set->sp.zres + zpos];

                            hx = vhh*cos(mp->src->ltsff->ia_phi[i1][im]) + vhe*sin(mp->src->ltsff->ia_phi[i1][im]);
                            hx = -hx;
                            epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, mp->src->ltsff->box_jn, k);

                            if (i < mp->src->ltsff->box_in)
                                mp->d->ex->data[i][mp->src->ltsff->box_jn][k] += acorr * ecor / epsilon * hz;
                            if (k < mp->src->ltsff->box_kn)
                                mp->d->ez->data[i][mp->src->ltsff->box_jn][k] -= acorr * ecor / epsilon * hx;
                        } // if (!set->ss.ltsff.skipjn)
                    } // k
                } // i

                /*k0, k1*/
#pragma omp parallel default(shared) private(i, j, zpos, tp, atp, hx, hy, vhh, vhe, epsilon)
#pragma omp for nowait

                for (i = mp->src->ltsff->box_i0; i <= mp->src->ltsff->box_in; i++) {
                    for (j = mp->src->ltsff->box_j0; j <= mp->src->ltsff->box_jn; j++) {
                        /*k0*/
                        if (!set->ss.ltsff.box_boundary_skipk0) {
                            sgetcoefs(i, j, mp->src->ltsff->box_k0 - 1, mp->src->ltsff->ia_theta[i1][im], mp->src->ltsff->ia_phi[i1][im],
                                      set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                                      mp->src->ltsff->box_i0, mp->src->ltsff->box_j0, mp->src->ltsff->box_k0,
                                      mp->src->ltsff->box_in, mp->src->ltsff->box_jn, mp->src->ltsff->box_kn,
                                      set->sp.xres, set->sp.yres, set->sp.zres);

                            vhh = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhh->data[(tp + 1)*set->sp.zres + zpos];
                            vhe = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhe->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhe->data[(tp + 1)*set->sp.zres + zpos];

                            hx = vhh*cos(mp->src->ltsff->ia_phi[i1][im]) + vhe*sin(mp->src->ltsff->ia_phi[i1][im]);
                            hy = vhh*sin(mp->src->ltsff->ia_phi[i1][im]) - vhe*cos(mp->src->ltsff->ia_phi[i1][im]);
                            hx = -hx;
                            hy = -hy;
                            epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, j, mp->src->ltsff->box_k0);

                            if (i < mp->src->ltsff->box_in)
                                mp->d->ex->data[i][j][mp->src->ltsff->box_k0] += acorr * ecor / epsilon * hy;
                            if (j < mp->src->ltsff->box_jn)
                                mp->d->ey->data[i][j][mp->src->ltsff->box_k0] -= acorr * ecor / epsilon * hx;
                        } // if (!set->ss.ltsff.skipk0)
                        //k1
                        //

                        if (!set->ss.ltsff.box_boundary_skipkn) {
                            sgetcoefs(i, j, mp->src->ltsff->box_kn, mp->src->ltsff->ia_theta[i1][im], mp->src->ltsff->ia_phi[i1][im],
                                      set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                                      mp->src->ltsff->box_i0, mp->src->ltsff->box_j0, mp->src->ltsff->box_k0,
                                      mp->src->ltsff->box_in, mp->src->ltsff->box_jn, mp->src->ltsff->box_kn,
                                      set->sp.xres, set->sp.yres, set->sp.zres);

                            vhh = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhh->data[(tp + 1)*set->sp.zres + zpos];
                            vhe = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dhe->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dhe->data[(tp + 1)*set->sp.zres + zpos];

                            hx = vhh*cos(mp->src->ltsff->ia_phi[i1][im]) + vhe*sin(mp->src->ltsff->ia_phi[i1][im]);
                            hy = vhh*sin(mp->src->ltsff->ia_phi[i1][im]) - vhe*cos(mp->src->ltsff->ia_phi[i1][im]);
                            hx = -hx;
                            hy = -hy;
                            epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, j, mp->src->ltsff->box_kn);

                            if (i < mp->src->ltsff->box_in)
                                mp->d->ex->data[i][j][mp->src->ltsff->box_kn] -= acorr * ecor / epsilon * hy;
                            if (j < mp->src->ltsff->box_jn)
                                mp->d->ey->data[i][j][mp->src->ltsff->box_kn] += acorr * ecor / epsilon * hx;
                        } //  if (!set->ss.ltsff.skipkn)
                    } // j
                } // i
            } // im
        } // in
    } //  if (mp->src->ltsff)

    if (mp->src->sf) {
        ecor = set->plan.dt / set->sp.dx;

        /*update 1d field*/
        sv_1dpool_ystep_e(mp->src->sf->jp);
        sv_1dpool_boundary(mp->src->sf->jp);
        sv_1dpool_boundary_copy(mp->src->sf->jp);
        sv_1dpool_apply_source(mp->src->sf->jp, mp->src->sf->e[set->sc.step_act]);
        sv_1dpool_ystep_h(mp->src->sf->jp);

        if (set->sc.step_act == 0) { /*one step more to synchronize with TSF source*/
            sv_1dpool_ystep_e(mp->src->sf->jp);
            sv_1dpool_boundary(mp->src->sf->jp);
            sv_1dpool_boundary_copy(mp->src->sf->jp);
            sv_1dpool_apply_source(mp->src->sf->jp, mp->src->sf->e[set->sc.step_act]);
            sv_1dpool_ystep_h(mp->src->sf->jp);
        }

        /*tangential field values are applied in pool_ystep_e*/
    } // if (mp->src->sf)

    if (set->ss.ext.filebase_ex) { /*external source*/
        g_snprintf(filename, 100, "%s_%.4d", set->ss.ext.filebase_hx, set->sc.step_act - set->ss.ext.shift);
        fhx = fopen(filename, "r");
        g_snprintf(filename, 100, "%s_%.4d", set->ss.ext.filebase_hy, set->sc.step_act - set->ss.ext.shift);
        fhy = fopen(filename, "r");
        g_snprintf(filename, 100, "%s_%.4d", set->ss.ext.filebase_hz, set->sc.step_act - set->ss.ext.shift);
        fhz = fopen(filename, "r");

        if (fhx == NULL || fhy == NULL || fhz == NULL) {
            fprintf(stderr, "Error: cannot open external source files (e.g. %s shifted step %d), skipping.\n", filename, set->sc.step_act - set->ss.ext.shift);
        } else {
            /*load entire file*/
            hxs = gwy_data_field_new(set->ss.ext.extxres, set->ss.ext.extyres, set->ss.ext.extxres, set->ss.ext.extyres, FALSE);
            hys = gwy_data_field_new(set->ss.ext.extxres, set->ss.ext.extyres, set->ss.ext.extxres, set->ss.ext.extyres, FALSE);
            hzs = gwy_data_field_new(set->ss.ext.extxres, set->ss.ext.extyres, set->ss.ext.extxres, set->ss.ext.extyres, FALSE);

            hxsdata = gwy_data_field_get_data(hxs);
            hysdata = gwy_data_field_get_data(hys);
            hzsdata = gwy_data_field_get_data(hzs);

            for (j = 0; j < set->ss.ext.extyres; j++) {
                for (i = 0; i < set->ss.ext.extxres; i++) {
                    fscanf(fhx, "%f", &vhx);
                    hxsdata[j*set->ss.ext.extxres + i] = vhx;

                    fscanf(fhy, "%f", &vhy);
                    hysdata[j*set->ss.ext.extxres + i] = vhy;

                    fscanf(fhz, "%f", &vhz);
                    hzsdata[j*set->ss.ext.extxres + i] = vhz;
                }
            }
            fclose(fhx);
            fclose(fhy);
            fclose(fhz);

            /* now use requested part of the plane to feed it into main pool*/
            if (set->ss.ext.i == -1 && set->ss.ext.j == -1) {
                for (i = set->ss.ext.iextfrom; i < set->ss.ext.iextto; i++) {
                    for (j = set->ss.ext.jextfrom; j < set->ss.ext.jextto; j++) {
                        ipos = i - set->ss.ext.iextfrom + set->ss.ext.ijstart;
                        jpos = j - set->ss.ext.jextfrom + set->ss.ext.jkstart;
                        if (ipos < 0 && ipos >= set->sp.xres)
                            continue;
                        if (jpos < 0 && jpos >= set->sp.yres)
                            continue;

                        ecor = set->plan.dt / set->sp.dx / sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, ipos, jpos, set->ss.ext.k);

                        mp->d->ex->data[ipos][jpos][set->ss.ext.k] += ecor * hysdata[j*set->ss.ext.extxres + i];
                        mp->d->ey->data[ipos][jpos][set->ss.ext.k] -= ecor * hxsdata[j*set->ss.ext.extxres + i];
                    } // j
                } // i
            } // if

            g_object_unref(hxs);
            g_object_unref(hys);
            g_object_unref(hzs);
        } // else
    } //  if (set->ss.ext.filebase_ex)

    if (mp->src->sh) {
        //printf("sh!\n");
    }

    if (set->sc.verbose > 1)
        printf("done.\n");
}

static gdouble
get_shifted_val(gdouble *e, gint i, gdouble timeshift, gint n)
{
    gint its = (gint)timeshift;
    gint ipos = i - its;
    gdouble diff = timeshift - (gdouble)its;
    if ((ipos > 0) && (ipos < (n - 1)))
        return e[ipos] + diff * (e[ipos + 1] - e[ipos]);
    else return 0;
}

void
sv_pool_apply_source_hstep(SvPool *mp, SvSet *set)
{
    gint i = 0, j = 0, k = 0, i1 = 0, im = 0, pos = 0, zpos = 0, tp = 0;
    gdouble d = 0, ve = 0, hcor = 0, theta = 0, phi = 0, psi = 0, corr = 0, acorr = 0, an = 0, bm = 0;
    gdouble stheta = 0, sphi = 0, spsi = 0, ctheta = 0, cphi = 0, cpsi = 0;
    gdouble ex = 0, ey = 0, ez = 0, atp = 0, veh = 0, vee = 0, ghcor = 0;
    gdouble rphi = 0, xradmult = 0, yradmult = 0;

    gchar filename[100];
    FILE *fex, *fey, *fez;
    gfloat vex, vey, vez;
    GwyDataField *exs, *eys, *ezs;
    gdouble *exsdata, *eysdata, *ezsdata;

    if (set->sc.verbose > 1) {
        printf("Running source hstep...   ");
        fflush(stdout);
    }

    for (i = 0; i < set->ss.npnts; i++) {
        /*if we are lucky, data are not shifted by user and can be used directly*/
        pos = -1;
        if (set->sc.step_act <= mp->src->sp[i].sdata.ndata && mp->src->sp[i].sdata.layered_zpos[set->sc.step_act] == set->sc.step_act)
            pos = set->sc.step_act;
        else { /*otherwise we need to search for right position*/
            for (j = 0; j < mp->src->sp[i].sdata.ndata; j++) {
                if (mp->src->sp[i].sdata.layered_zpos[j] == set->sc.step_act) {
                    pos = j;
                    break;
                }
            }
        }
        if (pos >= 0 && pos < mp->src->sp[i].sdata.ndata) {
            if (mp->src->sp[i].sdata.hx[pos] != 0)
                mp->d->hx->data[mp->src->sp[i].i][mp->src->sp[i].j][mp->src->sp[i].k] = mp->src->sp[i].sdata.hx[pos];
            if (mp->src->sp[i].sdata.hy[pos] != 0)
                mp->d->hy->data[mp->src->sp[i].i][mp->src->sp[i].j][mp->src->sp[i].k] = mp->src->sp[i].sdata.hy[pos];
            if (mp->src->sp[i].sdata.hz[pos] != 0)
                mp->d->hz->data[mp->src->sp[i].i][mp->src->sp[i].j][mp->src->sp[i].k] = mp->src->sp[i].sdata.hz[pos];
                
        }
    } // i

    if (mp->src->tsf) {
        hcor = set->plan.dt / set->sp.dx / MU_0;

        if (mp->src->tsf->jp == NULL) fprintf(stderr, "Error, tsf jp not initialized\n"); //FIXME

        sv_1dpool_ystep_e(mp->src->tsf->jp);
        sv_1dpool_boundary(mp->src->tsf->jp);
        sv_1dpool_boundary_copy(mp->src->tsf->jp);
        sv_1dpool_apply_source(mp->src->tsf->jp, mp->src->tsf->e[set->sc.step_act]);
        sv_1dpool_ystep_h(mp->src->tsf->jp);

        /*i0, i1*/
#pragma omp parallel default(shared) private(j, k, d, ve, ghcor, rphi, xradmult, yradmult)
#pragma omp for nowait

        for (j = mp->src->tsf->box_j0; j <= mp->src->tsf->box_jn; j++) {
            for (k = mp->src->tsf->box_k0; k <= mp->src->tsf->box_kn; k++) {
                if (!set->ss.tsf.box_boundary_skipi0 && mat_too_close_check(mp, set, mp->src->tsf->box_i0, j, k)) {
                    if (set->ss.tsf.gaussian)
                        ghcor = hcor * gaussmult(mp->src->tsf->box_i0, j, k, set->ss.tsf.gaussian_fxpos, set->ss.tsf.gaussian_fypos,
                                                 set->ss.tsf.gaussian_rx, set->ss.tsf.gaussian_ry);
                    else if (set->ss.tsf.fiber)
                        ghcor = hcor*fibermult(mp->src->tsf->box_i0, j, k, set->ss.tsf.fiber_fxpos, set->ss.tsf.fiber_fypos,
                                               set->ss.tsf.fiber_radius, set->ss.tsf.fiber_cutoff,
                                               set->ss.tsf.fiber_epsilon_core, set->ss.tsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        ghcor = hcor;

                    d = mp->src->tsf->corr * dcomp(mp->src->tsf->box_i0, j, k, set->sp.xres, set->sp.yres, set->sp.zres, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres);
                    ve = sv_dline_get_dval(mp->src->tsf->jp->e, d);

                    xradmult = yradmult = 1;
                    if (set->ss.tsf.radial == 1) { //z radial polarisation test
                        rphi = atan2(i - set->ss.tsf.radial_fxpos, j - set->ss.tsf.radial_fypos);
                        xradmult = sin(rphi);
                        yradmult = cos(rphi);
                        if (((i - set->ss.tsf.radial_fxpos)*(i - set->ss.tsf.radial_fxpos) + (j - set->ss.tsf.radial_fypos)*(j - set->ss.tsf.radial_fypos)) < 10) xradmult = yradmult = 0;
                    }


                    if (k < mp->src->tsf->box_kn)
                        mp->d->hy->data[mp->src->tsf->box_i0 - 1][j][k] -= xradmult * ghcor * gez(ve, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                    if (j < mp->src->tsf->box_jn)
                        mp->d->hz->data[mp->src->tsf->box_i0 - 1][j][k] += ghcor * gey(ve, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                } //  if (!set->ss.tsf.skipi0)

                //i1
                if (!set->ss.tsf.box_boundary_skipin && mat_too_close_check(mp, set, mp->src->tsf->box_in, j, k)) {
                    if (set->ss.tsf.gaussian)
                        ghcor = hcor * gaussmult(mp->src->tsf->box_in, j, k, set->ss.tsf.gaussian_fxpos, set->ss.tsf.gaussian_fypos,
                                                 set->ss.tsf.gaussian_rx, set->ss.tsf.gaussian_ry);
                    else if (set->ss.tsf.fiber)
                        ghcor = hcor * fibermult(mp->src->tsf->box_in, j, k, set->ss.tsf.fiber_fxpos, set->ss.tsf.fiber_fypos,
                                                 set->ss.tsf.fiber_radius, set->ss.tsf.fiber_cutoff,
                                                 set->ss.tsf.fiber_epsilon_core, set->ss.tsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        ghcor = hcor;

                    d = mp->src->tsf->corr * dcomp(mp->src->tsf->box_in, j, k, set->sp.xres, set->sp.yres, set->sp.zres,
                                                   mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres);
                    ve = sv_dline_get_dval(mp->src->tsf->jp->e, d);

                    xradmult = yradmult = 1;
                    if (set->ss.tsf.radial == 1) { //z radial polarisation test
                        rphi = atan2(i - set->ss.tsf.radial_fxpos, j - set->ss.tsf.radial_fypos);
                        xradmult = sin(rphi);
                        yradmult = cos(rphi);
                        if (((i - set->ss.tsf.radial_fxpos)*(i - set->ss.tsf.radial_fxpos) + (j - set->ss.tsf.radial_fypos)*(j - set->ss.tsf.radial_fypos)) < 10)
                            xradmult = yradmult = 0;
                    }

                    if (k < mp->src->tsf->box_kn)
                        mp->d->hy->data[mp->src->tsf->box_in][j][k] += xradmult * ghcor * gez(ve, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                    if (j < mp->src->tsf->box_jn)
                        mp->d->hz->data[mp->src->tsf->box_in][j][k] -= ghcor * gey(ve, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                } // if (!set->ss.tsf.skipin)
            } // k
        } // j

        /*j0, j1*/
#pragma omp parallel default(shared) private(i, k, d, ve, ghcor, rphi, xradmult, yradmult)
#pragma omp for nowait

        for (i = mp->src->tsf->box_i0; i <= mp->src->tsf->box_in; i++) {
            for (k = mp->src->tsf->box_k0; k <= mp->src->tsf->box_kn; k++) {
                //j0
                if (!set->ss.tsf.box_boundary_skipj0 && mat_too_close_check(mp, set, i, mp->src->tsf->box_j0, k)) {
                    if (set->ss.tsf.gaussian)
                        ghcor = hcor * gaussmult(i, mp->src->tsf->box_j0, k, set->ss.tsf.gaussian_fxpos, set->ss.tsf.gaussian_fypos, set->ss.tsf.gaussian_rx, set->ss.tsf.gaussian_ry);
                    else if (set->ss.tsf.fiber)
                        ghcor = hcor * fibermult(i, mp->src->tsf->box_j0, k, set->ss.tsf.fiber_fxpos, set->ss.tsf.fiber_fypos, set->ss.tsf.fiber_radius, set->ss.tsf.fiber_cutoff, set->ss.tsf.fiber_epsilon_core, set->ss.tsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        ghcor = hcor;

                    d = mp->src->tsf->corr * dcomp(i, mp->src->tsf->box_j0, k, set->sp.xres, set->sp.yres, set->sp.zres, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres);
                    ve = sv_dline_get_dval(mp->src->tsf->jp->e, d);

                    xradmult = yradmult = 1;
                    if (set->ss.tsf.radial == 1) { //z radial polarisation test                      
                        rphi = atan2(i - set->ss.tsf.radial_fxpos, j - set->ss.tsf.radial_fypos);
                        xradmult = sin(rphi);
                        yradmult = cos(rphi);
                        if (((i - set->ss.tsf.radial_fxpos)*(i - set->ss.tsf.radial_fxpos) + (j - set->ss.tsf.radial_fypos)*(j - set->ss.tsf.radial_fypos)) < 10)
                            xradmult = yradmult = 0;
                    }


                    //horni vec
                    if (k < mp->src->tsf->box_kn)
                        mp->d->hx->data[i][mp->src->tsf->box_j0 - 1][k] += yradmult * ghcor * gez(ve, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                    if (i < mp->src->tsf->box_in)
                        mp->d->hz->data[i][mp->src->tsf->box_j0 - 1][k] -= ghcor * gex(ve, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                } //  if (!set->ss.tsf.skipj0)

                //j1
                if (!set->ss.tsf.box_boundary_skipjn && mat_too_close_check(mp, set, i, mp->src->tsf->box_jn, k)) {
                    if (set->ss.tsf.gaussian)
                        ghcor = hcor * gaussmult(i, mp->src->tsf->box_jn, k, set->ss.tsf.gaussian_fxpos, set->ss.tsf.gaussian_fypos, set->ss.tsf.gaussian_rx, set->ss.tsf.gaussian_ry);
                    else if (set->ss.tsf.fiber)
                        ghcor = hcor * fibermult(i, mp->src->tsf->box_jn, k, set->ss.tsf.fiber_fxpos, set->ss.tsf.fiber_fypos, set->ss.tsf.fiber_radius, set->ss.tsf.fiber_cutoff, set->ss.tsf.fiber_epsilon_core, set->ss.tsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        ghcor = hcor;

                    d = mp->src->tsf->corr * dcomp(i, mp->src->tsf->box_jn, k, set->sp.xres, set->sp.yres, set->sp.zres, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres);
                    ve = sv_dline_get_dval(mp->src->tsf->jp->e, d);

                    xradmult = yradmult = 1;
                    if (set->ss.tsf.radial == 1) { //z radial polarisation test
                        rphi = atan2(i - set->ss.tsf.radial_fxpos, j - set->ss.tsf.radial_fypos);
                        xradmult = sin(rphi);
                        yradmult = cos(rphi);
                        if (((i - set->ss.tsf.radial_fxpos)*(i - set->ss.tsf.radial_fxpos) + (j - set->ss.tsf.radial_fypos)*(j - set->ss.tsf.radial_fypos)) < 10)
                            xradmult = yradmult = 0;
                    }


                    if (k < mp->src->tsf->box_kn)
                        mp->d->hx->data[i][mp->src->tsf->box_jn][k] -= yradmult * ghcor * gez(ve, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                    if (i < mp->src->tsf->box_in)
                        mp->d->hz->data[i][mp->src->tsf->box_jn][k] += ghcor * gex(ve, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                } // if (!set->ss.tsf.skipjn)
            } // k
        } // i

        /*k0, k1*/
#pragma omp parallel default(shared) private(i, j, d, ve, ghcor, xradmult, yradmult, rphi)
#pragma omp for nowait

        for (i = mp->src->tsf->box_i0; i <= mp->src->tsf->box_in; i++) {
            for (j = mp->src->tsf->box_j0; j <= mp->src->tsf->box_jn; j++) {
                //k0
                if (!set->ss.tsf.box_boundary_skipk0 && mat_too_close_check(mp, set, i, j, mp->src->tsf->box_k0)) {
                    if (set->ss.tsf.gaussian)
                        ghcor = hcor * gaussmult(i, j, mp->src->tsf->box_k0, set->ss.tsf.gaussian_fxpos, set->ss.tsf.gaussian_fypos, set->ss.tsf.gaussian_rx, set->ss.tsf.gaussian_ry);
                    else if (set->ss.tsf.fiber)
                        ghcor = hcor*fibermult(i, j, mp->src->tsf->box_k0, set->ss.tsf.fiber_fxpos, set->ss.tsf.fiber_fypos, set->ss.tsf.fiber_radius, set->ss.tsf.fiber_cutoff, set->ss.tsf.fiber_epsilon_core, set->ss.tsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        ghcor = hcor;


                    d = mp->src->tsf->corr * dcomp(i, j, mp->src->tsf->box_k0, set->sp.xres, set->sp.yres, set->sp.zres, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres);
                    ve = sv_dline_get_dval(mp->src->tsf->jp->e, d);

                    xradmult = yradmult = 1;
                    if (set->ss.tsf.radial == 1) { //z radial polarisation test
                        rphi = atan2(i - set->ss.tsf.radial_fxpos, j - set->ss.tsf.radial_fypos);
                        xradmult = sin(rphi);
                        yradmult = cos(rphi);

                        if (((i - set->ss.tsf.radial_fxpos)*(i - set->ss.tsf.radial_fxpos) + (j - set->ss.tsf.radial_fypos)*(j - set->ss.tsf.radial_fypos)) < 10)
                            xradmult = yradmult = 0;
                    }


                    if (j < mp->src->tsf->box_jn && !(sv_yee_data_is_pec(mp->d, set, mp->mats, mp->nmat, i, j, mp->src->tsf->box_k0)))
                        mp->d->hx->data[i][j][mp->src->tsf->box_k0 - 1] -= yradmult * ghcor * gey(ve, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                    if (i < mp->src->tsf->box_in && !(sv_yee_data_is_pec(mp->d, set, mp->mats, mp->nmat, i, j, mp->src->tsf->box_k0)))
                        mp->d->hy->data[i][j][mp->src->tsf->box_k0 - 1] += xradmult * ghcor * gex(ve, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                } // if (!set->ss.tsf.skipk0)

                //k1
                if (!set->ss.tsf.box_boundary_skipkn && mat_too_close_check(mp, set, i, j, mp->src->tsf->box_kn)) {
                    if (set->ss.tsf.gaussian)
                        ghcor = hcor * gaussmult(i, j, mp->src->tsf->box_kn, set->ss.tsf.gaussian_fxpos, set->ss.tsf.gaussian_fypos, set->ss.tsf.gaussian_rx, set->ss.tsf.gaussian_ry);
                    else if (set->ss.tsf.fiber)
                        ghcor = hcor * fibermult(i, j, mp->src->tsf->box_kn, set->ss.tsf.fiber_fxpos, set->ss.tsf.fiber_fypos, set->ss.tsf.fiber_radius, set->ss.tsf.fiber_cutoff, set->ss.tsf.fiber_epsilon_core, set->ss.tsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        ghcor = hcor;


                    d = mp->src->tsf->corr * dcomp(i, j, mp->src->tsf->box_kn, set->sp.xres, set->sp.yres, set->sp.zres, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres);
                    ve = sv_dline_get_dval(mp->src->tsf->jp->e, d);

                    xradmult = yradmult = 1;
                    if (set->ss.tsf.radial == 1) { //z radial polarisation test
                        rphi = atan2(i - set->ss.tsf.radial_fxpos, j - set->ss.tsf.radial_fypos);
                        xradmult = sin(rphi);
                        yradmult = cos(rphi);
                        if (((i - set->ss.tsf.radial_fxpos)*(i - set->ss.tsf.radial_fxpos) + (j - set->ss.tsf.radial_fypos)*(j - set->ss.tsf.radial_fypos)) < 10) xradmult = yradmult = 0;
                    }

                    if (j < mp->src->tsf->box_jn)
                        mp->d->hx->data[i][j][mp->src->tsf->box_kn] += yradmult * ghcor * gey(ve, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                    if (i < mp->src->tsf->box_in)
                        mp->d->hy->data[i][j][mp->src->tsf->box_kn] -= xradmult * ghcor * gex(ve, mp->src->tsf->ia_theta, mp->src->tsf->ia_phi, mp->src->tsf->ia_psi);
                } //  if (!set->ss.tsf.skipkn)
            } // j
        } // i
    } //  if (mp->src->tsf)

    if (mp->src->ltsf) {
        hcor = set->plan.dt / set->sp.dx / MU_0;

        if (set->sc.step_act == 0) {
            if (set->sc.verbose > 1) {
                printf("(precomputing LTSF data...   ");
                fflush(stdout);
            }
            for (i = 0; i < (3 * set->sc.nsteps); i++) { //we use dt/3 for the auxiliary LTSF data
                sv_zpool_boundary_copy(mp->src->ltsf->jp);
                sv_zpool_ystep_e(mp->src->ltsf->jp);
                sv_zpool_boundary(mp->src->ltsf->jp);
                sv_zpool_apply_source(mp->src->ltsf->jp, get_shifted_val(mp->src->ltsf->e, i, mp->src->ltsf->timeshift, mp->src->ltsf->ndata));
                sv_zpool_ystep_h(mp->src->ltsf->jp);
                sv_zpool_store(mp->src->ltsf->jp, i);

                /*
                                sv_zpool_boundary_copy(mp->src->ltsf->jp);
                                sv_zpool_ystep_e(mp->src->ltsf->jp);
                                sv_zpool_boundary(mp->src->ltsf->jp);
                                sv_zpool_apply_source(mp->src->ltsf->jp, get_shifted_val(mp->src->ltsf->e, i, mp->src->ltsf->timeshift, mp->src->ltsf->ndata));
                                sv_zpool_ystep_h(mp->src->ltsf->jp);
                                sv_zpool_store(mp->src->ltsf->jp, i);
                */
            }
        }

        /*i0, i1*/
#pragma omp parallel default(shared) private(j, k, zpos, tp, atp, ey, ez, ghcor, veh, vee)
#pragma omp for nowait

        for (j = mp->src->ltsf->box_j0; j <= mp->src->ltsf->box_jn; j++) {
            for (k = mp->src->ltsf->box_k0; k <= mp->src->ltsf->box_kn; k++) {
                if (!set->ss.ltsf.box_boundary_skipi0) {
                    if (set->ss.ltsf.gaussian)
                        ghcor = hcor * gaussmult(mp->src->ltsf->box_i0, j, k, set->ss.ltsf.gaussian_fxpos, set->ss.ltsf.gaussian_fypos, set->ss.ltsf.gaussian_rx, set->ss.ltsf.gaussian_ry);
                    else if (set->ss.ltsf.fiber)
                        ghcor = hcor*fibermult(mp->src->ltsf->box_i0, j, k, set->ss.ltsf.fiber_fxpos, set->ss.ltsf.fiber_fypos, set->ss.ltsf.fiber_radius, set->ss.ltsf.fiber_cutoff, set->ss.ltsf.fiber_epsilon_core, set->ss.ltsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        ghcor = hcor;

                    getcoefs(mp->src->ltsf->box_i0, j, k, mp->src->ltsf->ia_theta, mp->src->ltsf->ia_phi,
                             set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                             mp->src->ltsf->box_i0, mp->src->ltsf->box_j0, mp->src->ltsf->box_k0,
                             mp->src->ltsf->box_in, mp->src->ltsf->box_jn, mp->src->ltsf->box_kn);

                    veh = (1.0 - atp) * mp->src->ltsf->jp->deh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->deh->data[(tp + 1)*set->sp.zres + zpos];
                    vee = (1.0 - atp) * mp->src->ltsf->jp->dee->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dee->data[(tp + 1)*set->sp.zres + zpos];
                    ez = (1.0 - atp) * mp->src->ltsf->jp->dez->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dez->data[(tp + 1)*set->sp.zres + zpos];

                    ey = vee*sin(mp->src->ltsf->ia_phi) + veh*cos(mp->src->ltsf->ia_phi);

                    if (k < mp->src->ltsf->box_kn)
                        mp->d->hy->data[mp->src->ltsf->box_i0 - 1][j][k] -= ghcor*ez;
                    if (j < mp->src->ltsf->box_jn)
                        mp->d->hz->data[mp->src->ltsf->box_i0 - 1][j][k] += ghcor*ey;
                } //  if (!set->ss.ltsf.skipi0)

                //i1
                if (!set->ss.ltsf.box_boundary_skipin) {
                    if (set->ss.ltsf.gaussian)
                        ghcor = hcor * gaussmult(mp->src->ltsf->box_in, j, k, set->ss.ltsf.gaussian_fxpos, set->ss.ltsf.gaussian_fypos, set->ss.ltsf.gaussian_rx, set->ss.ltsf.gaussian_ry);
                    if (set->ss.ltsf.fiber)
                        ghcor = hcor * fibermult(mp->src->ltsf->box_in, j, k, set->ss.ltsf.fiber_fxpos, set->ss.ltsf.fiber_fypos,
                                                 set->ss.ltsf.fiber_radius, set->ss.ltsf.fiber_cutoff,
                                                 set->ss.ltsf.fiber_epsilon_core, set->ss.ltsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        ghcor = hcor;

                    getcoefs(mp->src->ltsf->box_in, j, k, mp->src->ltsf->ia_theta, mp->src->ltsf->ia_phi,
                             set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                             mp->src->ltsf->box_i0, mp->src->ltsf->box_j0, mp->src->ltsf->box_k0,
                             mp->src->ltsf->box_in, mp->src->ltsf->box_jn, mp->src->ltsf->box_kn);

                    veh = (1.0 - atp) * mp->src->ltsf->jp->deh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->deh->data[(tp + 1)*set->sp.zres + zpos];
                    vee = (1.0 - atp) * mp->src->ltsf->jp->dee->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dee->data[(tp + 1)*set->sp.zres + zpos];
                    ez = (1.0 - atp) * mp->src->ltsf->jp->dez->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dez->data[(tp + 1)*set->sp.zres + zpos];

                    ey = vee*sin(mp->src->ltsf->ia_phi) + veh*cos(mp->src->ltsf->ia_phi);

                    if (k < mp->src->ltsf->box_kn)
                        mp->d->hy->data[mp->src->ltsf->box_in][j][k] += ghcor*ez;
                    if (j < mp->src->ltsf->box_jn)
                        mp->d->hz->data[mp->src->ltsf->box_in][j][k] -= ghcor*ey;
                } // (!set->ss.ltsf.skipin)
            } // k
        } // j

        /*j0, j1*/
#pragma omp parallel default(shared) private(i, k, zpos, tp, atp, ex, ez, ghcor, veh, vee)
#pragma omp for nowait

        for (i = mp->src->ltsf->box_i0; i <= mp->src->ltsf->box_in; i++) {
            for (k = mp->src->ltsf->box_k0; k <= mp->src->ltsf->box_kn; k++) {
                /*j0*/
                if (!set->ss.ltsf.box_boundary_skipj0) {
                    if (set->ss.ltsf.gaussian)
                        ghcor = hcor * gaussmult(i, mp->src->ltsf->box_j0, k, set->ss.ltsf.gaussian_fxpos, set->ss.ltsf.gaussian_fypos, set->ss.ltsf.gaussian_rx, set->ss.ltsf.gaussian_ry);
                    else if (set->ss.ltsf.fiber)
                        ghcor = hcor*fibermult(i, mp->src->ltsf->box_j0, k, set->ss.ltsf.fiber_fxpos, set->ss.ltsf.fiber_fypos,
                                               set->ss.ltsf.fiber_radius, set->ss.ltsf.fiber_cutoff,
                                               set->ss.ltsf.fiber_epsilon_core, set->ss.ltsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        ghcor = hcor;

                    getcoefs(i, mp->src->ltsf->box_j0, k, mp->src->ltsf->ia_theta, mp->src->ltsf->ia_phi,
                             set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                             mp->src->ltsf->box_i0, mp->src->ltsf->box_j0, mp->src->ltsf->box_k0,
                             mp->src->ltsf->box_in, mp->src->ltsf->box_jn, mp->src->ltsf->box_kn);

                    veh = (1.0 - atp) * mp->src->ltsf->jp->deh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->deh->data[(tp + 1)*set->sp.zres + zpos];
                    vee = (1.0 - atp) * mp->src->ltsf->jp->dee->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dee->data[(tp + 1)*set->sp.zres + zpos];
                    ez = (1.0 - atp) * mp->src->ltsf->jp->dez->data[tp*set->sp.zres + zpos] + atp *mp->src->ltsf->jp->dez->data[(tp + 1)*set->sp.zres + zpos];

                    ex = vee*cos(mp->src->ltsf->ia_phi) - veh*sin(mp->src->ltsf->ia_phi);

                    if (k < mp->src->ltsf->box_kn)
                        mp->d->hx->data[i][mp->src->ltsf->box_j0 - 1][k] += ghcor*ez;
                    if (i < mp->src->ltsf->box_in)
                        mp->d->hz->data[i][mp->src->ltsf->box_j0 - 1][k] -= ghcor*ex;
                } //  if (!set->ss.ltsf.skipj0)

                //j1
                if (!set->ss.ltsf.box_boundary_skipjn) {
                    if (set->ss.ltsf.gaussian)
                        ghcor = hcor * gaussmult(i, mp->src->ltsf->box_jn, k, set->ss.ltsf.gaussian_fxpos, set->ss.ltsf.gaussian_fypos, set->ss.ltsf.gaussian_rx, set->ss.ltsf.gaussian_ry);
                    else if (set->ss.ltsf.fiber)
                        ghcor = hcor * fibermult(i, mp->src->ltsf->box_jn, k, set->ss.ltsf.fiber_fxpos, set->ss.ltsf.fiber_fypos,
                                                 set->ss.ltsf.fiber_radius, set->ss.ltsf.fiber_cutoff,
                                                 set->ss.ltsf.fiber_epsilon_core, set->ss.ltsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        ghcor = hcor;

                    getcoefs(i, mp->src->ltsf->box_jn, k, mp->src->ltsf->ia_theta, mp->src->ltsf->ia_phi,
                             set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                             mp->src->ltsf->box_i0, mp->src->ltsf->box_j0, mp->src->ltsf->box_k0,
                             mp->src->ltsf->box_in, mp->src->ltsf->box_jn, mp->src->ltsf->box_kn);

                    veh = (1.0 - atp) * mp->src->ltsf->jp->deh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->deh->data[(tp + 1)*set->sp.zres + zpos];
                    vee = (1.0 - atp) * mp->src->ltsf->jp->dee->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dee->data[(tp + 1)*set->sp.zres + zpos];
                    ez = (1.0 - atp) * mp->src->ltsf->jp->dez->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dez->data[(tp + 1)*set->sp.zres + zpos];

                    ex = vee*cos(mp->src->ltsf->ia_phi) - veh*sin(mp->src->ltsf->ia_phi);

                    if (k < mp->src->ltsf->box_kn)
                        mp->d->hx->data[i][mp->src->ltsf->box_jn][k] -= ghcor*ez;
                    if (i < mp->src->ltsf->box_in)
                        mp->d->hz->data[i][mp->src->ltsf->box_jn][k] += ghcor*ex;
                } // if (!set->ss.ltsf.skipjn)
            } // k
        } // i

        /*k0, k1*/
#pragma omp parallel default(shared) private(i, j, zpos, tp, atp, ex, ey, ghcor, veh, vee)
#pragma omp for nowait

        for (i = mp->src->ltsf->box_i0; i <= mp->src->ltsf->box_in; i++) {
            for (j = mp->src->ltsf->box_j0; j <= mp->src->ltsf->box_jn; j++) {
                /*k0*/
                if (!set->ss.ltsf.box_boundary_skipk0  && mat_too_close_check(mp, set, i, j, mp->src->ltsf->box_k0)) {
                    if (set->ss.ltsf.gaussian)
                        ghcor = hcor * gaussmult(i, j, mp->src->ltsf->box_k0, set->ss.ltsf.gaussian_fxpos, set->ss.ltsf.gaussian_fypos, set->ss.ltsf.gaussian_rx, set->ss.ltsf.gaussian_ry);
                    else if (set->ss.ltsf.fiber)
                        ghcor = hcor*fibermult(i, j, mp->src->ltsf->box_k0, set->ss.ltsf.fiber_fxpos, set->ss.ltsf.fiber_fypos,
                                               set->ss.ltsf.fiber_radius, set->ss.ltsf.fiber_cutoff,
                                               set->ss.ltsf.fiber_epsilon_core, set->ss.ltsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        ghcor = hcor;

                    getcoefs(i, j, mp->src->ltsf->box_k0, mp->src->ltsf->ia_theta, mp->src->ltsf->ia_phi,
                             set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                             mp->src->ltsf->box_i0, mp->src->ltsf->box_j0, mp->src->ltsf->box_k0,
                             mp->src->ltsf->box_in, mp->src->ltsf->box_jn, mp->src->ltsf->box_kn);

                    veh = (1.0 - atp) * mp->src->ltsf->jp->deh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->deh->data[(tp + 1)*set->sp.zres + zpos];
                    vee = (1.0 - atp) * mp->src->ltsf->jp->dee->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dee->data[(tp + 1)*set->sp.zres + zpos];

                    ex = vee*cos(mp->src->ltsf->ia_phi) - veh*sin(mp->src->ltsf->ia_phi);
                    ey = vee*sin(mp->src->ltsf->ia_phi) + veh*cos(mp->src->ltsf->ia_phi);

                    if (j < mp->src->ltsf->box_jn)
                        mp->d->hx->data[i][j][mp->src->ltsf->box_k0 - 1] -= ghcor*ey;
                    if (i < mp->src->ltsf->box_in)
                        mp->d->hy->data[i][j][mp->src->ltsf->box_k0 - 1] += ghcor*ex;
                } //  if (!set->ss.ltsf.skipk0)

                //k1
                if (!set->ss.ltsf.box_boundary_skipkn) {
                    if (set->ss.ltsf.gaussian)
                        ghcor = hcor * gaussmult(i, j, mp->src->ltsf->box_kn, set->ss.ltsf.gaussian_fxpos, set->ss.ltsf.gaussian_fypos,
                                                 set->ss.ltsf.gaussian_rx, set->ss.ltsf.gaussian_ry);
                    else if (set->ss.ltsf.fiber)
                        ghcor = hcor * fibermult(i, j, mp->src->ltsf->box_kn, set->ss.ltsf.fiber_fxpos, set->ss.ltsf.fiber_fypos,
                                                 set->ss.ltsf.fiber_radius, set->ss.ltsf.fiber_cutoff,
                                                 set->ss.ltsf.fiber_epsilon_core, set->ss.ltsf.fiber_epsilon_cladding, set->sp.dx, set->ss.lambda_center);
                    else
                        ghcor = hcor;

                    getcoefs(i, j, mp->src->ltsf->box_kn, mp->src->ltsf->ia_theta, mp->src->ltsf->ia_phi,
                             set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                             mp->src->ltsf->box_i0, mp->src->ltsf->box_j0, mp->src->ltsf->box_k0,
                             mp->src->ltsf->box_in, mp->src->ltsf->box_jn, mp->src->ltsf->box_kn);

                    veh = (1.0 - atp) * mp->src->ltsf->jp->deh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->deh->data[(tp + 1)*set->sp.zres + zpos];
                    vee = (1.0 - atp) * mp->src->ltsf->jp->dee->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsf->jp->dee->data[(tp + 1)*set->sp.zres + zpos];

                    ex = vee*cos(mp->src->ltsf->ia_phi) - veh*sin(mp->src->ltsf->ia_phi);
                    ey = vee*sin(mp->src->ltsf->ia_phi) + veh*cos(mp->src->ltsf->ia_phi);

                    if (j < mp->src->ltsf->box_jn)
                        mp->d->hx->data[i][j][mp->src->ltsf->box_kn] += ghcor*ey;
                    if (i < mp->src->ltsf->box_in)
                        mp->d->hy->data[i][j][mp->src->ltsf->box_kn] -= ghcor*ex;
                } // if (!set->ss.ltsf.skipkn)
            } // j
        } // i
    } //  if (mp->src->ltsf)

    if (mp->src->tsff) {
        hcor = set->plan.dt / set->sp.dx / MU_0;

        for (i1 = 0; i1 < mp->src->tsff->focused_nip; i1++) {
            for (im = 0; im < mp->src->tsff->focused_mip; im++) {
                //  if (in!=(mp->src->tsff->nint-1) || im!=(mp->src->tsff->mint-1)) continue;

              ///if (!(i1==0 && im==0)) continue;


                sv_1dpool_ystep_e(mp->src->tsff->jp[i1][im]);
                sv_1dpool_boundary(mp->src->tsff->jp[i1][im]);
                sv_1dpool_boundary_copy(mp->src->tsff->jp[i1][im]);
                sv_1dpool_apply_source(mp->src->tsff->jp[i1][im], mp->src->tsff->e[set->sc.step_act]);
                sv_1dpool_ystep_h(mp->src->tsff->jp[i1][im]);

                theta = mp->src->tsff->ia_theta[i1][im];
                phi = mp->src->tsff->ia_phi[i1][im];
                psi = mp->src->tsff->ia_psi[i1][im];
                corr = mp->src->tsff->corr[i1][im];
                an = mp->src->tsff->an[i1][im];
                bm = mp->src->tsff->bm[i1][im];

                ctheta = cos(theta); stheta = sin(theta);
                cphi = cos(phi); sphi = sin(phi);
                cpsi = cos(psi); spsi = sin(psi);

                acorr = fs_amplitude(mp->src->tsff->focused_fdist*set->sp.dx, an, bm, theta);

                /*i0, i1*/
#pragma omp parallel default(shared) private(j, k, d, ve)
#pragma omp for nowait

                for (j = mp->src->tsff->box_j0; j <= mp->src->tsff->box_jn; j++) {
                    for (k = mp->src->tsff->box_k0; k <= mp->src->tsff->box_kn; k++) {
                        if (!set->ss.tsff.box_boundary_skipi0) {
                            d = corr * srdcomp(mp->src->tsff->box_i0, j, k, set->sp.xres, set->sp.yres, set->sp.zres,
                                               stheta, sphi, ctheta, cphi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres,
                                               set->ss.tsff.xshift, set->ss.tsff.yshift, set->ss.tsff.zshift);
                            ve = acorr * sv_dline_get_dval(mp->src->tsff->jp[i1][im]->e, d);

                            if (k < mp->src->tsff->box_kn)
                                mp->d->hy->data[mp->src->tsff->box_i0 - 1][j][k] -= hcor * ve * (spsi*stheta); //*sgez(ve, stheta, sphi, spsi, ctheta, cphi, cpsi);
                            if (j < mp->src->tsff->box_jn)
                                mp->d->hz->data[mp->src->tsff->box_i0 - 1][j][k] += hcor * ve * (-cpsi*cphi - spsi*ctheta*sphi); //*sgey(ve, stheta, sphi, spsi, ctheta, cphi, cpsi);

                        }

                        //i1
                        if (!set->ss.tsff.box_boundary_skipin) {
                            d = corr * srdcomp(mp->src->tsff->box_in, j, k, set->sp.xres, set->sp.yres, set->sp.zres,
                                               stheta, sphi, ctheta, cphi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres,
                                               set->ss.tsff.xshift, set->ss.tsff.yshift, set->ss.tsff.zshift);
                            ve = acorr * sv_dline_get_dval(mp->src->tsff->jp[i1][im]->e, d);

                            if (k < mp->src->tsff->box_kn)
                                mp->d->hy->data[mp->src->tsff->box_in][j][k] += hcor * ve * (spsi*stheta); //*sgez(ve, stheta, sphi, spsi, ctheta, cphi, cpsi);
                            if (j < mp->src->tsff->box_jn)
                                mp->d->hz->data[mp->src->tsff->box_in][j][k] -= hcor * ve * (-cpsi*cphi - spsi*ctheta*sphi); //*sgey(ve, stheta, sphi, spsi, ctheta, cphi, cpsi);
                        }
                    } // k
                } // j

                /*j0, j1*/
#pragma omp parallel default(shared) private(i, k, d, ve)
#pragma omp for nowait

                for (i = mp->src->tsff->box_i0; i <= mp->src->tsff->box_in; i++) {
                    for (k = mp->src->tsff->box_k0; k <= mp->src->tsff->box_kn; k++) {
                        /*j0*/
                        if (!set->ss.tsff.box_boundary_skipj0) {
                            d = corr * srdcomp(i, mp->src->tsff->box_j0, k, set->sp.xres, set->sp.yres, set->sp.zres, stheta, sphi, ctheta, cphi,
                                               0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres, set->ss.tsff.xshift, set->ss.tsff.yshift, set->ss.tsff.zshift);
                            ve = acorr * sv_dline_get_dval(mp->src->tsff->jp[i1][im]->e, d);

                            //horni vec
                            if (k < mp->src->tsff->box_kn)
                                mp->d->hx->data[i][mp->src->tsff->box_j0 - 1][k] += hcor * ve * (spsi*stheta); //*sgez(ve, stheta, sphi, spsi, ctheta, cphi, cpsi);
                            if (i < mp->src->tsff->box_in)
                                mp->d->hz->data[i][mp->src->tsff->box_j0 - 1][k] -= hcor * ve * (cpsi*sphi - spsi*ctheta*cphi); //*sgex(ve, stheta, sphi, spsi, ctheta, cphi, cpsi);
                        }

                        //j1
                        if (!set->ss.tsff.box_boundary_skipjn) {
                            d = corr * srdcomp(i, mp->src->tsff->box_jn, k, set->sp.xres, set->sp.yres, set->sp.zres,
                                               stheta, sphi, ctheta, cphi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres,
                                               set->ss.tsff.xshift, set->ss.tsff.yshift, set->ss.tsff.zshift);
                            ve = acorr * sv_dline_get_dval(mp->src->tsff->jp[i1][im]->e, d);

                            if (k < mp->src->tsff->box_kn)
                                mp->d->hx->data[i][mp->src->tsff->box_jn][k] -= hcor * ve * (spsi*stheta); //*sgez(ve, stheta, sphi, spsi, ctheta, cphi, cpsi);
                            if (i < mp->src->tsff->box_in)
                                mp->d->hz->data[i][mp->src->tsff->box_jn][k] += hcor * ve * (cpsi*sphi - spsi*ctheta*cphi); //*sgex(ve, stheta, sphi, spsi, ctheta, cphi, cpsi);
                        }
                    } // k
                } // i

                /*k0, k1*/
#pragma omp parallel default(shared) private(i, j, d, ve)
#pragma omp for nowait

                for (i = mp->src->tsff->box_i0; i <= mp->src->tsff->box_in; i++) {
                    for (j = mp->src->tsff->box_j0; j <= mp->src->tsff->box_jn; j++) {
                        /*k0*/
                        if (!set->ss.tsff.box_boundary_skipk0) {
                            d = corr * srdcomp(i, j, mp->src->tsff->box_k0, set->sp.xres, set->sp.yres, set->sp.zres,
                                               stheta, sphi, ctheta, cphi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres,
                                               set->ss.tsff.xshift, set->ss.tsff.yshift, set->ss.tsff.zshift);
                            ve = acorr * sv_dline_get_dval(mp->src->tsff->jp[i1][im]->e, d);

                            if (j < mp->src->tsff->box_jn)
                                mp->d->hx->data[i][j][mp->src->tsff->box_k0 - 1] -= hcor*ve*(-cpsi*cphi - spsi*ctheta*sphi); //*sgey(ve, stheta, sphi, spsi, ctheta, cphi, cpsi);
                            if (i < mp->src->tsff->box_in)
                                mp->d->hy->data[i][j][mp->src->tsff->box_k0 - 1] += hcor*ve*(cpsi*sphi - spsi*ctheta*cphi);  //*sgex(ve, stheta, sphi, spsi, ctheta, cphi, cpsi);
                        }

                        //k1
                        if (!set->ss.tsff.box_boundary_skipkn) {
                            d = corr * srdcomp(i, j, mp->src->tsff->box_kn, set->sp.xres, set->sp.yres, set->sp.zres,
                                               stheta, sphi, ctheta, cphi, 0, set->sp.xres, 0, set->sp.yres, 0, set->sp.zres,
                                               set->ss.tsff.xshift, set->ss.tsff.yshift, set->ss.tsff.zshift);
                            ve = acorr * sv_dline_get_dval(mp->src->tsff->jp[i1][im]->e, d);

                            if (j < mp->src->tsff->box_jn)
                                mp->d->hx->data[i][j][mp->src->tsff->box_kn] += hcor * ve * (-cpsi*cphi - spsi*ctheta*sphi); //*sgey(ve, stheta, sphi, spsi, ctheta, cphi, cpsi);
                            if (i < mp->src->tsff->box_in)
                                mp->d->hy->data[i][j][mp->src->tsff->box_kn] -= hcor * ve * (cpsi*sphi - spsi*ctheta*cphi);  //*sgex(ve, stheta, sphi, spsi, ctheta, cphi, cpsi);
                        }
                    } // j
                } // i
            } // im
        } // in
    } // if (mp->src->tsff)

    if (mp->src->ltsff) {
        hcor = set->plan.dt / set->sp.dx / MU_0;

        for (i1 = 0; i1 < mp->src->ltsff->focused_nip; i1++) {
            for (im = 0; im < mp->src->ltsff->focused_mip; im++) {
                if (set->sc.step_act == 0) {
                    if (set->sc.verbose > 1) {
                        printf("(precomputing LTSF data...   ");
                        fflush(stdout);
                    }

                    for (i = 0; i < (3 * set->sc.nsteps); i++) { //we use dt/3 for the auxiliary LTSF data  //add timeshift here to source
                        sv_zpool_ystep_e(mp->src->ltsff->jp[i1][im]);
                        sv_zpool_boundary(mp->src->ltsff->jp[i1][im]);
                        sv_zpool_boundary_copy(mp->src->ltsff->jp[i1][im]);
                        sv_zpool_apply_source(mp->src->ltsff->jp[i1][im], get_shifted_val(mp->src->ltsff->e, i, mp->src->ltsff->timeshift[i1][im], mp->src->ltsff->ndata));
                        sv_zpool_ystep_h(mp->src->ltsff->jp[i1][im]);
                        sv_zpool_store(mp->src->ltsff->jp[i1][im], i);
                    }
                }

                theta = mp->src->ltsff->ia_theta[i1][im];
                phi = mp->src->ltsff->ia_phi[i1][im];
                psi = mp->src->ltsff->ia_psi[i1][im];
                corr = mp->src->ltsff->corr[i1][im];
                an = mp->src->ltsff->an[i1][im];
                bm = mp->src->ltsff->bm[i1][im];

                ctheta = cos(theta); stheta = sin(theta);
                cphi = cos(phi); sphi = sin(phi);
                cpsi = cos(psi); spsi = sin(psi);

                acorr = fs_amplitude(mp->src->ltsff->focused_fdist*set->sp.dx, an, bm, theta);

                /*i0, i1*/
#pragma omp parallel default(shared) private(j, k, zpos, tp, atp, ey, ez, veh, vee)
#pragma omp for nowait

                for (j = mp->src->ltsff->box_j0; j <= mp->src->ltsff->box_jn; j++) {
                    for (k = mp->src->ltsff->box_k0; k <= mp->src->ltsff->box_kn; k++) {
                        if (!set->ss.ltsff.box_boundary_skipi0) {
                            sgetcoefs(mp->src->ltsff->box_i0, j, k, mp->src->ltsff->ia_theta[i1][im], mp->src->ltsff->ia_phi[i1][im],
                                      set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                                      mp->src->ltsff->box_i0, mp->src->ltsff->box_j0, mp->src->ltsff->box_k0,
                                      mp->src->ltsff->box_in, mp->src->ltsff->box_jn, mp->src->ltsff->box_kn,
                                      set->sp.xres, set->sp.yres, set->sp.zres);

                            veh = (1.0 - atp) *mp->src->ltsff->jp[i1][im]->deh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->deh->data[(tp + 1)*set->sp.zres + zpos];
                            vee = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dee->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dee->data[(tp + 1)*set->sp.zres + zpos];
                            ez = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dez->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dez->data[(tp + 1)*set->sp.zres + zpos];

                            ey = vee*sin(mp->src->ltsff->ia_phi[i1][im]) + veh*cos(mp->src->ltsff->ia_phi[i1][im]);

                            if (k < mp->src->ltsff->box_kn)
                                mp->d->hy->data[mp->src->ltsff->box_i0 - 1][j][k] -= acorr*hcor*ez;
                            if (j < mp->src->ltsff->box_jn)
                                mp->d->hz->data[mp->src->ltsff->box_i0 - 1][j][k] += acorr*hcor*ey;
                        }

                        //i1
                        if (!set->ss.ltsff.box_boundary_skipin) {
                            sgetcoefs(mp->src->ltsff->box_in, j, k, mp->src->ltsff->ia_theta[i1][im], mp->src->ltsff->ia_phi[i1][im],
                                      set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                                      mp->src->ltsff->box_i0, mp->src->ltsff->box_j0, mp->src->ltsff->box_k0,
                                      mp->src->ltsff->box_in, mp->src->ltsff->box_jn, mp->src->ltsff->box_kn,
                                      set->sp.xres, set->sp.yres, set->sp.zres);

                            veh = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->deh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->deh->data[(tp + 1)*set->sp.zres + zpos];
                            vee = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dee->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dee->data[(tp + 1)*set->sp.zres + zpos];
                            ez = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dez->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dez->data[(tp + 1)*set->sp.zres + zpos];

                            ey = vee*sin(mp->src->ltsff->ia_phi[i1][im]) + veh*cos(mp->src->ltsff->ia_phi[i1][im]);

                            if (k < mp->src->ltsff->box_kn)
                                mp->d->hy->data[mp->src->ltsff->box_in][j][k] += acorr*hcor*ez;
                            if (j < mp->src->ltsff->box_jn)
                                mp->d->hz->data[mp->src->ltsff->box_in][j][k] -= acorr*hcor*ey;
                        }
                    } // k
                } // i

                /*j0, j1*/
#pragma omp parallel default(shared) private(i, k, zpos, tp, atp, ex, ez, veh, vee)
#pragma omp for nowait

                for (i = mp->src->ltsff->box_i0; i <= mp->src->ltsff->box_in; i++) {
                    for (k = mp->src->ltsff->box_k0; k <= mp->src->ltsff->box_kn; k++) {
                        /*j0*/
                        if (!set->ss.ltsff.box_boundary_skipj0) {
                            sgetcoefs(i, mp->src->ltsff->box_j0, k, mp->src->ltsff->ia_theta[i1][im], mp->src->ltsff->ia_phi[i1][im],
                                      set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                                      mp->src->ltsff->box_i0, mp->src->ltsff->box_j0, mp->src->ltsff->box_k0,
                                      mp->src->ltsff->box_in, mp->src->ltsff->box_jn, mp->src->ltsff->box_kn,
                                      set->sp.xres, set->sp.yres, set->sp.zres);

                            veh = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->deh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->deh->data[(tp + 1)*set->sp.zres + zpos];
                            vee = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dee->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dee->data[(tp + 1)*set->sp.zres + zpos];
                            ez = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dez->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dez->data[(tp + 1)*set->sp.zres + zpos];

                            ex = vee*cos(mp->src->ltsff->ia_phi[i1][im]) - veh*sin(mp->src->ltsff->ia_phi[i1][im]);

                            if (k < mp->src->ltsff->box_kn)
                                mp->d->hx->data[i][mp->src->ltsff->box_j0 - 1][k] += acorr*hcor*ez;
                            if (i < mp->src->ltsff->box_in)
                                mp->d->hz->data[i][mp->src->ltsff->box_j0 - 1][k] -= acorr*hcor*ex;
                        }

                        //j1
                        if (!set->ss.ltsff.box_boundary_skipjn) {
                            sgetcoefs(i, mp->src->ltsff->box_jn, k, mp->src->ltsff->ia_theta[i1][im], mp->src->ltsff->ia_phi[i1][im],
                                      set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                                      mp->src->ltsff->box_i0, mp->src->ltsff->box_j0, mp->src->ltsff->box_k0,
                                      mp->src->ltsff->box_in, mp->src->ltsff->box_jn, mp->src->ltsff->box_kn,
                                      set->sp.xres, set->sp.yres, set->sp.zres);

                            veh = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->deh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->deh->data[(tp + 1)*set->sp.zres + zpos];
                            vee = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dee->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dee->data[(tp + 1)*set->sp.zres + zpos];
                            ez = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dez->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dez->data[(tp + 1)*set->sp.zres + zpos];

                            ex = vee*cos(mp->src->ltsff->ia_phi[i1][im]) - veh*sin(mp->src->ltsff->ia_phi[i1][im]);

                            if (k < mp->src->ltsff->box_kn)
                                mp->d->hx->data[i][mp->src->ltsff->box_jn][k] -= acorr*hcor*ez;
                            if (i < mp->src->ltsff->box_in)
                                mp->d->hz->data[i][mp->src->ltsff->box_jn][k] += acorr*hcor*ex;
                        }
                    } // k
                } // i

                /*k0, k1*/
#pragma omp parallel default(shared) private(i, j, zpos, tp, atp, ex, ey, veh, vee)
#pragma omp for nowait

                for (i = mp->src->ltsff->box_i0; i <= mp->src->ltsff->box_in; i++) {
                    for (j = mp->src->ltsff->box_j0; j <= mp->src->ltsff->box_jn; j++) {
                        /*k0*/
                        if (!set->ss.ltsff.box_boundary_skipk0) {
                            sgetcoefs(i, j, mp->src->ltsff->box_k0, mp->src->ltsff->ia_theta[i1][im], mp->src->ltsff->ia_phi[i1][im],
                                      set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                                      mp->src->ltsff->box_i0, mp->src->ltsff->box_j0, mp->src->ltsff->box_k0,
                                      mp->src->ltsff->box_in, mp->src->ltsff->box_jn, mp->src->ltsff->box_kn,
                                      set->sp.xres, set->sp.yres, set->sp.zres);

                            veh = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->deh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->deh->data[(tp + 1)*set->sp.zres + zpos];
                            vee = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dee->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dee->data[(tp + 1)*set->sp.zres + zpos];

                            ex = vee*cos(mp->src->ltsff->ia_phi[i1][im]) - veh*sin(mp->src->ltsff->ia_phi[i1][im]);
                            ey = vee*sin(mp->src->ltsff->ia_phi[i1][im]) + veh*cos(mp->src->ltsff->ia_phi[i1][im]);

                            if (j < mp->src->ltsff->box_jn)
                                mp->d->hx->data[i][j][mp->src->ltsff->box_k0 - 1] -= acorr*hcor*ey;
                            if (i < mp->src->ltsff->box_in)
                                mp->d->hy->data[i][j][mp->src->ltsff->box_k0 - 1] += acorr*hcor*ex;
                        }

                        //k1
                        if (!set->ss.ltsff.box_boundary_skipkn) {
                            sgetcoefs(i, j, mp->src->ltsff->box_kn, mp->src->ltsff->ia_theta[i1][im], mp->src->ltsff->ia_phi[i1][im],
                                      set->sc.step_act, &zpos, &tp, &atp, set->sp.dx, set->plan.dt,
                                      mp->src->ltsff->box_i0, mp->src->ltsff->box_j0, mp->src->ltsff->box_k0,
                                      mp->src->ltsff->box_in, mp->src->ltsff->box_jn, mp->src->ltsff->box_kn,
                                      set->sp.xres, set->sp.yres, set->sp.zres);

                            veh = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->deh->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->deh->data[(tp + 1)*set->sp.zres + zpos];
                            vee = (1.0 - atp) * mp->src->ltsff->jp[i1][im]->dee->data[tp*set->sp.zres + zpos] + atp * mp->src->ltsff->jp[i1][im]->dee->data[(tp + 1)*set->sp.zres + zpos];

                            ex = vee*cos(mp->src->ltsff->ia_phi[i1][im]) - veh*sin(mp->src->ltsff->ia_phi[i1][im]);
                            ey = vee*sin(mp->src->ltsff->ia_phi[i1][im]) + veh*cos(mp->src->ltsff->ia_phi[i1][im]);

                            if (j < mp->src->ltsff->box_jn)
                                mp->d->hx->data[i][j][mp->src->ltsff->box_kn] += acorr*hcor*ey;
                            if (i < mp->src->ltsff->box_in)
                                mp->d->hy->data[i][j][mp->src->ltsff->box_kn] -= acorr*hcor*ex;
                        }
                    } // j
                } // i
            } // im
        } // in
    } //  if (mp->src->tsff)

    if (set->ss.ext.filebase_ex) { /*external source*/
        g_snprintf(filename, 100, "%s_%.4d", set->ss.ext.filebase_ex, set->sc.step_act - set->ss.ext.shift);
        fex = fopen(filename, "r");
        g_snprintf(filename, 100, "%s_%.4d", set->ss.ext.filebase_ey, set->sc.step_act - set->ss.ext.shift);
        fey = fopen(filename, "r");
        g_snprintf(filename, 100, "%s_%.4d", set->ss.ext.filebase_ez, set->sc.step_act - set->ss.ext.shift);
        fez = fopen(filename, "r");

        if (fex == NULL || fey == NULL || fez == NULL) {
            fprintf(stderr, "Error: cannot open external source files (e.g. %s, shifted step %d), skipping.\n", filename, set->sc.step_act - set->ss.ext.shift);
        } else {
            /*load entire file*/
            exs = gwy_data_field_new(set->ss.ext.extxres, set->ss.ext.extyres, set->ss.ext.extxres, set->ss.ext.extyres, FALSE);
            eys = gwy_data_field_new(set->ss.ext.extxres, set->ss.ext.extyres, set->ss.ext.extxres, set->ss.ext.extyres, FALSE);
            ezs = gwy_data_field_new(set->ss.ext.extxres, set->ss.ext.extyres, set->ss.ext.extxres, set->ss.ext.extyres, FALSE);

            exsdata = gwy_data_field_get_data(exs);
            eysdata = gwy_data_field_get_data(eys);
            ezsdata = gwy_data_field_get_data(ezs);

            for (j = 0; j < set->ss.ext.extyres; j++) {
                for (i = 0; i < set->ss.ext.extxres; i++) {

                    fscanf(fex, "%f", &vex);
                    exsdata[j*set->ss.ext.extxres + i] = vex;

                    fscanf(fey, "%f", &vey);
                    eysdata[j*set->ss.ext.extxres + i] = vey;

                    fscanf(fez, "%f", &vez);
                    ezsdata[j*set->ss.ext.extxres + i] = vez;
                }
            }

            fclose(fex);
            fclose(fey);
            fclose(fez);

            hcor = set->plan.dt / set->sp.dx / MU_0;

            /*now use requested part of the plane to feed it into main pool*/
            if (set->ss.ext.i == -1 && set->ss.ext.j == -1) {
                for (i = set->ss.ext.iextfrom; i < set->ss.ext.iextto; i++) {
                    for (j = set->ss.ext.jextfrom; j < set->ss.ext.jextto; j++) {
                        if ((i - set->ss.ext.iextfrom + set->ss.ext.ijstart) < 0 && (i - set->ss.ext.iextfrom + set->ss.ext.ijstart) >= set->sp.xres)
                            continue;
                        if ((j - set->ss.ext.jextfrom + set->ss.ext.jkstart) < 0 && (j - set->ss.ext.jextfrom + set->ss.ext.jkstart) >= set->sp.yres)
                            continue;

                        mp->d->hx->data[i - set->ss.ext.iextfrom + set->ss.ext.ijstart][j - set->ss.ext.jextfrom + set->ss.ext.jkstart][set->ss.ext.k - 1]
                            -= hcor*eysdata[j*set->ss.ext.extxres + i];
                        mp->d->hy->data[i - set->ss.ext.iextfrom + set->ss.ext.ijstart][j - set->ss.ext.jextfrom + set->ss.ext.jkstart][set->ss.ext.k - 1]
                            += hcor*exsdata[j*set->ss.ext.extxres + i];
                    } // j
                } // i
            } // if

            g_object_unref(exs);
            g_object_unref(eys);
            g_object_unref(ezs);
        } // else
    } //  if (set->ss.ext.filebase_ex)

    if (set->sc.verbose > 1)
        printf("done.\n");
}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
