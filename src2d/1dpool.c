
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


/*  1dpool.c : 
 *  one-dimensional FDTD for use in TSF and SF sources
 */

#if defined (_WIN32) && defined (_MSC_VER)
#include "../config_win.h"
#else
#include "../config.h"
#endif
#include <glib.h>
#include <stdlib.h>
#include "1dpool.h"
#include "constants.h"


Sv1DPool* 
sv_1dpool_new(gint res, gdouble dx, gdouble dt)
{
    Sv1DPool *jp = (Sv1DPool *)g_malloc(sizeof(Sv1DPool));
    
    jp->e = sv_dline_new(res, res*dx, 1); 
    jp->h = sv_dline_new_alike(jp->e, 1);
    jp->epsilon = sv_dline_new_alike(jp->e, 1);
    sv_dline_fill(jp->epsilon, 1);
    jp->mu = sv_dline_new_alike(jp->e, 1);
    sv_dline_fill(jp->mu, 1);
    jp->sigma = sv_dline_new_alike(jp->e, 1);
    jp->sigast = sv_dline_new_alike(jp->e, 1);
    jp->source = sv_dline_new_alike(jp->e, 1);

    jp->step = 0;
    jp->dx = dx;
    jp->dt = dt;
    jp->res = res;

    jp->fbnx0_e_0=0;
    jp->fbnx0_e_1=0;
    jp->fbnx0_h_0=0;
    jp->fbnx0_h_1=0;
    jp->fbnxn_e_0=0;
    jp->fbnxn_e_1=0;
    jp->fbnxn_h_0=0;
    jp->fbnxn_h_1=0;

    return jp;
}

void 
sv_1dpool_set_source_pos(Sv1DPool *jp, gint i)
{
    jp->source->data[i] = 1;
}

void
sv_1dpool_set_material(Sv1DPool *jp, gdouble epsilon, gdouble mu)
{
    gint i;
    for (i=0; i<jp->res; i++)
    {
        jp->mu->data[i] = mu;
        jp->epsilon->data[i] = epsilon;
    }
}

void 
sv_1dpool_ystep_e(Sv1DPool *jp)
{
    gint i;

    gdouble angle_mult=1; //for dAngle=0 or 90 FIXME for others
    gdouble multe=jp->dt/(EPSILON_0*jp->dx)/angle_mult;

    for (i=0; i<jp->res; i++) //x
    {

        /*update E*/
        if (i>0)
        {
            jp->e->data[i]+= 1.0/jp->epsilon->data[i]*multe*(jp->h->data[i-1] - jp->h->data[i]);
        }

    }


    jp->step++;
}

void 
sv_1dpool_ystep_h(Sv1DPool *jp)
{
    gint i;
//    gchar buff[50];
//    FILE *fr;

    gdouble angle_mult=1;
    gdouble multh=jp->dt/(MU_0*jp->dx)/angle_mult;
 
//    sprintf(buff, "tsrc%.3d.txt", jp->step); 
//    fr = fopen(buff, "w");


    for (i=0; i<jp->res; i++) //x
    {

        /*update H*/
        if (i<(jp->res-1))
        {
            jp->h->data[i]+= 1.0/jp->mu->data[i]*multh*(jp->e->data[i] - jp->e->data[i+1]);

        }
    }

    
//    for (i=0; i<jp->res; i++)
//    {
//        fprintf(fr, "%d %g %g\n", i, jp->e->data[i], jp->h->data[i]);
//    }
//    fclose(fr);
    

}

void 
sv_1dpool_boundary_copy(Sv1DPool *jp)
{
    jp->fbnx0_e_0=jp->e->data[0];
    jp->fbnx0_e_1=jp->e->data[1];
    jp->fbnx0_h_0=jp->h->data[0];
    jp->fbnx0_h_1=jp->h->data[1];
    jp->fbnxn_e_0=jp->e->data[jp->res-1];
    jp->fbnxn_e_1=jp->e->data[jp->res-2];
    jp->fbnxn_h_0=jp->h->data[jp->res-1];
    jp->fbnxn_h_1=jp->h->data[jp->res-2];
}

void 
sv_1dpool_boundary(Sv1DPool *jp)
{
    jp->e->data[0]=jp->fbnx0_e_1 +
        (jp->dt*LIGHT_SPEED-jp->dx)/(jp->dt*LIGHT_SPEED+jp->dx)*
        (jp->e->data[1] - jp->fbnx0_e_0);

    jp->h->data[0]=jp->fbnx0_h_1 +
        (jp->dt*LIGHT_SPEED-jp->dx)/(jp->dt*LIGHT_SPEED+jp->dx)*
        (jp->h->data[1] - jp->fbnx0_h_0);

    jp->e->data[jp->res-1]=jp->fbnxn_e_1 +
        (jp->dt*LIGHT_SPEED-jp->dx)/(jp->dt*LIGHT_SPEED+jp->dx)*
        (jp->e->data[jp->res-2] - jp->fbnxn_e_0);

    jp->h->data[jp->res-1]=jp->fbnxn_h_1 +
        (jp->dt*LIGHT_SPEED-jp->dx)/(jp->dt*LIGHT_SPEED+jp->dx)*
        (jp->h->data[jp->res-2] - jp->fbnxn_h_0);

}

void 
sv_1dpool_apply_source(Sv1DPool *jp, gdouble eval)
{
    gint i;

    /*for each point marked as the source apply source*/
    for (i=0; i<jp->res; i++)
    {
        if (jp->source->data[i]==1)
        {
            jp->e->data[i] = eval;
        }
    }

}


/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */



