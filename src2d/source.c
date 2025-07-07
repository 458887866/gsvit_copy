
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

SvSource* 
sv_source_new()
{
    SvSource* src = (SvSource*)g_malloc(sizeof(SvSource));
  
    src->tsf = NULL; 
    return src;
}

gdouble
dcomp(gint i, gint j, gint k, gint xres, gint yres, gint zres, gdouble theta, gdouble phi,
      gint i0, gint i1, gint j0, gint j1, gint k0, gint k1)
{
        gdouble rx, ry, rz;
        gdouble ax, ay, az;
        
        //ax = sin(theta)*cos(phi);
        //ay = sin(theta)*sin(phi);
        //az = cos(theta);
        ax = cos(phi);
        ay = sin(phi);
        az = 0;


        if (phi >= 0 && phi <= (G_PI/2.0))
        {
            rx = i - i0;
            ry = j - j0;
            rz = k - k0;
        }
        else if (phi > (G_PI/2.0) && phi <= G_PI)
        {
            rx = i - i1;
            ry = j - j0;
            rz = k - k0;
        }
        else if (phi > G_PI && phi <= (3.0*G_PI/2.0))
        {
            rx = i - i1;
            ry = j - j1;
            rz = k - k0;
        }
        else
        {
            rx = i - i0;
            ry = j - j1;
            rz = k - k0;
        }

        return 10 + (ax*rx + ay*ry + az*rz);

}


gdouble gex(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return field*(cos(psi)*sin(phi));
}
gdouble gey(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return field*(-cos(psi)*cos(phi));
}
gdouble gez(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return field*(sin(psi));
}
gdouble ghx(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return field*(sin(psi)*sin(phi));
}
gdouble ghy(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return field*(-sin(psi)*cos(phi));
}
gdouble ghz(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return field*(-cos(psi));
}


/*
gdouble gex(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return 0;
}
gdouble gey(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return 0;
}
gdouble gez(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return field;
}
gdouble ghx(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return field*sin(phi);
}
gdouble ghy(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return field*cos(phi);
}
gdouble ghz(gdouble field, gdouble theta, gdouble phi, gdouble psi)
{
    return 0;
}
*/




void sv_pool_apply_source_estep(SvPool *mp, SvSet *set)
{
    gint i, j, pos;
    gdouble d, vh, ecor;
    
    if (set->sc.verbose) {
        printf("Running source...   ");
        fflush(stdout);
    }

    for (i=0; i<set->ss.npnts; i++)
    {
        /*if we are lucky, data are not shifted by user and can be used directly*/
        pos = -1;
        if (set->sc.step_act <= mp->src->sp[i].sdata.ndata && 
                 mp->src->sp[i].sdata.pos[set->sc.step_act] == set->sc.step_act) 
            pos = set->sc.step_act;
        else {  /*otherwise we need to search for right position*/
            for (j=0; j<mp->src->sp[i].sdata.ndata; j++) {
                if (mp->src->sp[i].sdata.pos[j] == set->sc.step_act)
                {
                    pos = j; 
                    break;
                }
            }
        }
        if (pos>=0 && pos<mp->src->sp[i].sdata.ndata)
        {
             if (mp->src->sp[i].sdata.ex[pos] != 0) 
                 mp->ex->data[mp->src->sp[i].j*set->sp.xres+mp->src->sp[i].i] = mp->src->sp[i].sdata.ex[pos];
             if (mp->src->sp[i].sdata.ey[pos] != 0) 
                 mp->ey->data[mp->src->sp[i].j*set->sp.xres+mp->src->sp[i].i] = mp->src->sp[i].sdata.ey[pos];
             if (mp->src->sp[i].sdata.ez[pos] != 0) 
                 mp->ez->data[mp->src->sp[i].j*set->sp.xres+mp->src->sp[i].i] = mp->src->sp[i].sdata.ez[pos];

        }

    }

    if (mp->src->tsf) {

        //i0,i1
        ecor = set->plan.dt/set->sp.dx;
        
        for (j=mp->src->tsf->j0; j<=mp->src->tsf->j1; j++)
        {
                 //i0, left
                 d = mp->src->tsf->corr*dcomp(mp->src->tsf->i0-1, j, 0,
                           set->sp.xres, set->sp.yres, 1,
                           mp->src->tsf->theta, mp->src->tsf->phi,
                           mp->src->tsf->i0, mp->src->tsf->i1,
                           mp->src->tsf->j0, mp->src->tsf->j1,
                           0, 0);

                 vh = sv_dline_get_dval(mp->src->tsf->jp->h, d);

                 if (!set->tmmode) { //te mode
                     if (j<mp->src->tsf->j1) mp->ey->data[mp->src->tsf->i0+j*set->sp.xres]
                         += ecor/mp->epsilon->data[mp->src->tsf->i0+j*set->sp.xres]*ghz(vh, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 } else { //tm mode
                     mp->ez->data[mp->src->tsf->i0+j*set->sp.xres]
                         -= ecor/mp->epsilon->data[mp->src->tsf->i0+j*set->sp.xres]*ghy(vh, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 }


                 //i1, right
                 d = mp->src->tsf->corr*dcomp(mp->src->tsf->i1, j, 0,
                           set->sp.xres, set->sp.yres, 1,
                           mp->src->tsf->theta, mp->src->tsf->phi,
                           mp->src->tsf->i0, mp->src->tsf->i1,
                           mp->src->tsf->j0, mp->src->tsf->j1,
                           0, 0);
                 vh = sv_dline_get_dval(mp->src->tsf->jp->h, d);


                 if (!set->tmmode) { //te mode
                     if (j<mp->src->tsf->j1) mp->ey->data[mp->src->tsf->i1+j*set->sp.xres]
                         -= ecor/mp->epsilon->data[mp->src->tsf->i1+j*set->sp.xres]*ghz(vh, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 } else { //tm mode
                     mp->ez->data[mp->src->tsf->i1+j*set->sp.xres]
                         += ecor/mp->epsilon->data[mp->src->tsf->i1+j*set->sp.xres]*ghy(vh, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 }

        }
        //j0, j1
        for (i=(mp->src->tsf->i0); i<=mp->src->tsf->i1; i++)
        {
                 //j0 top
                 d = mp->src->tsf->corr*dcomp(i, mp->src->tsf->j0-1, 0,
                           set->sp.xres, set->sp.yres, 0,
                           mp->src->tsf->theta, mp->src->tsf->phi,
                           mp->src->tsf->i0, mp->src->tsf->i1,
                           mp->src->tsf->j0, mp->src->tsf->j1,
                           0, 0);
                 vh = sv_dline_get_dval(mp->src->tsf->jp->h, d);

                 if (!set->tmmode) { //te mode
                     if (i<mp->src->tsf->i1) mp->ex->data[i + mp->src->tsf->j0*set->sp.xres]
                         -= ecor/mp->epsilon->data[i + mp->src->tsf->j0*set->sp.xres]*ghz(vh, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 }
                 else { //tm mode
                     mp->ez->data[i+mp->src->tsf->j0*set->sp.xres]
                         += ecor/mp->epsilon->data[i + mp->src->tsf->j0*set->sp.xres]*ghx(vh, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 }


                 //j1 bottom
                 d = mp->src->tsf->corr*dcomp(i, mp->src->tsf->j1, 0,
                           set->sp.xres, set->sp.yres, 0,
                           mp->src->tsf->theta, mp->src->tsf->phi,
                           mp->src->tsf->i0, mp->src->tsf->i1,
                           mp->src->tsf->j0, mp->src->tsf->j1,
                           0, 0);
                 vh = sv_dline_get_dval(mp->src->tsf->jp->h, d);

                 if (!set->tmmode) { //te mode
                     if (i<mp->src->tsf->i1) mp->ex->data[i + mp->src->tsf->j1*set->sp.xres]
                         += ecor/mp->epsilon->data[i + mp->src->tsf->j1*set->sp.xres]*ghz(vh, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 }
                 else { //tm mode
                     mp->ez->data[i + mp->src->tsf->j1*set->sp.xres]
                     -= ecor/mp->epsilon->data[i + mp->src->tsf->j1*set->sp.xres]*ghx(vh, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 }

        }
    }

    if (set->sc.verbose) printf("done.\n");

}



void sv_pool_apply_source_hstep(SvPool *mp, SvSet *set)
{
    gint i, j, pos;
    gdouble d, ve, hcor;
    
    if (set->sc.verbose) {
        printf("Running source...   ");
        fflush(stdout);
    }

    for (i=0; i<set->ss.npnts; i++)
    {
        /*if we are lucky, data are not shifted by user and can be used directly*/
        pos = -1;
        if (set->sc.step_act <= mp->src->sp[i].sdata.ndata && 
                 mp->src->sp[i].sdata.pos[set->sc.step_act] == set->sc.step_act) 
            pos = set->sc.step_act;
        else {  /*otherwise we need to search for right position*/
            for (j=0; j<mp->src->sp[i].sdata.ndata; j++) {
                if (mp->src->sp[i].sdata.pos[j] == set->sc.step_act)
                {
                    pos = j; 
                    break;
                }
            }
        }
        if (pos>=0 && pos<mp->src->sp[i].sdata.ndata)
        {
             if (mp->src->sp[i].sdata.hx[pos] != 0) 
                 mp->hx->data[mp->src->sp[i].j*set->sp.xres+mp->src->sp[i].i] = mp->src->sp[i].sdata.hx[pos];
             if (mp->src->sp[i].sdata.hy[pos] != 0) 
                 mp->hy->data[mp->src->sp[i].j*set->sp.xres+mp->src->sp[i].i] = mp->src->sp[i].sdata.hy[pos];
             if (mp->src->sp[i].sdata.hz[pos] != 0) 
                 mp->hz->data[mp->src->sp[i].j*set->sp.xres+mp->src->sp[i].i] = mp->src->sp[i].sdata.hz[pos];

        }

    }

    if (mp->src->tsf) {

        hcor = set->plan.dt/set->sp.dx/MU_0;


        /*update 1d field*/
        sv_1dpool_ystep_e(mp->src->tsf->jp);
        sv_1dpool_boundary(mp->src->tsf->jp);
        sv_1dpool_boundary_copy(mp->src->tsf->jp);
        sv_1dpool_apply_source(mp->src->tsf->jp, mp->src->tsf->e[set->sc.step_act]);
        sv_1dpool_ystep_h(mp->src->tsf->jp);


        for (j=mp->src->tsf->j0; j<=mp->src->tsf->j1; j++)
        {
                 //i0, left
                 d = mp->src->tsf->corr*dcomp(mp->src->tsf->i0, j, 0,
                           set->sp.xres, set->sp.yres, 1,
                           mp->src->tsf->theta, mp->src->tsf->phi,
                           mp->src->tsf->i0, mp->src->tsf->i1,
                           mp->src->tsf->j0, mp->src->tsf->j1,
                           0, 0);

                 ve = sv_dline_get_dval(mp->src->tsf->jp->e, d);

                 if (!set->tmmode) { //te mode
                     if (j<mp->src->tsf->j1) mp->hz->data[mp->src->tsf->i0-1+j*set->sp.xres] 
                         += hcor*gey(ve, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 } else { //tm mode
                     mp->hy->data[mp->src->tsf->i0-1+j*set->sp.xres]
                         -= hcor*gez(ve, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 }


                 //i1, right
                 d = mp->src->tsf->corr*dcomp(mp->src->tsf->i1, j, 0,
                           set->sp.xres, set->sp.yres, 1,
                           mp->src->tsf->theta, mp->src->tsf->phi,
                           mp->src->tsf->i0, mp->src->tsf->i1,
                           mp->src->tsf->j0, mp->src->tsf->j1,
                           0, 0);
                 ve = sv_dline_get_dval(mp->src->tsf->jp->e, d);


                 if (!set->tmmode) { //te mode
                     if (j<mp->src->tsf->j1) mp->hz->data[mp->src->tsf->i1+j*set->sp.xres]
                         -= hcor*gey(ve, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 } else { //tm mode
                     mp->hy->data[mp->src->tsf->i1+j*set->sp.xres]
                         += hcor*gez(ve, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 }

        }
        //j0, j1
        for (i=(mp->src->tsf->i0); i<=mp->src->tsf->i1; i++)
        {
                 //j0 top
                 d = mp->src->tsf->corr*dcomp(i, mp->src->tsf->j0, 0,
                           set->sp.xres, set->sp.yres, 0,
                           mp->src->tsf->theta, mp->src->tsf->phi,
                           mp->src->tsf->i0, mp->src->tsf->i1,
                           mp->src->tsf->j0, mp->src->tsf->j1,
                           0, 0);
                 ve = sv_dline_get_dval(mp->src->tsf->jp->e, d);

                 if (!set->tmmode) { //te mode
                     if (i<mp->src->tsf->i1) mp->hz->data[i + (mp->src->tsf->j0-1)*set->sp.xres]
                         -= hcor*gex(ve, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 }
                 else { //tm mode
                     mp->hx->data[i+(mp->src->tsf->j0-1)*set->sp.xres]
                         += hcor*gez(ve, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 }


                 //j1 bottom
                 d = mp->src->tsf->corr*dcomp(i, mp->src->tsf->j1, 0,
                           set->sp.xres, set->sp.yres, 0,
                           mp->src->tsf->theta, mp->src->tsf->phi,
                           mp->src->tsf->i0, mp->src->tsf->i1,
                           mp->src->tsf->j0, mp->src->tsf->j1,
                           0, 0);
                 ve = sv_dline_get_dval(mp->src->tsf->jp->e, d);


                 if (!set->tmmode) { //te mode
                     if (i<mp->src->tsf->i1) mp->hz->data[i + mp->src->tsf->j1*set->sp.xres]
                         += hcor*gex(ve, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 }
                 else { //tm mode
                     mp->hx->data[i + mp->src->tsf->j1*set->sp.xres]
                     -= hcor*gez(ve, mp->src->tsf->theta, mp->src->tsf->phi, mp->src->tsf->psi);
                 }

        }

    }

    if (set->sc.verbose) printf("done.\n");

}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

