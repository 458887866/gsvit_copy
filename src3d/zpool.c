
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


/*  zpool.c : 
 *  one-dimensional FDTD for use in TSF and SF sources
 */

#if defined (_WIN32) && defined (_MSC_VER)
#include "../config_win.h"
#else
#include "../config.h"
#endif
#include <glib.h>
#include <stdlib.h>
#include "zpool.h"
#include "constants.h"
#include <math.h>

SvZPool*
sv_zpool_new(gint res, gdouble dx, gdouble dt, gint nsteps)
{
    SvZPool *jp = (SvZPool *)g_malloc(sizeof(SvZPool));
    
    jp->dt = dt/3;
    jp->res = res;

    jp->ee = sv_dline_new(res, jp->res*dx, 1); 
    jp->eh = sv_dline_new_alike(jp->ee, 1); 
    jp->ez = sv_dline_new_alike(jp->ee, 1); 
    jp->he = sv_dline_new_alike(jp->ee, 1);
    jp->hh = sv_dline_new_alike(jp->ee, 1);
    jp->hz = sv_dline_new_alike(jp->ee, 1);
    jp->e = sv_dline_new(res, res*dx, 1);
    jp->h = sv_dline_new_alike(jp->e, 1);

    jp->nsteps = 3*nsteps;
    jp->dee = gwy_data_field_new(res, jp->nsteps, res, jp->nsteps, 1);
    jp->deh = gwy_data_field_new(res, jp->nsteps, res, jp->nsteps, 1);
    jp->dez = gwy_data_field_new(res, jp->nsteps, res, jp->nsteps, 1);
    jp->dhe = gwy_data_field_new(res, jp->nsteps, res, jp->nsteps, 1);
    jp->dhh = gwy_data_field_new(res, jp->nsteps, res, jp->nsteps, 1);
    jp->dhz = gwy_data_field_new(res, jp->nsteps, res, jp->nsteps, 1);

    jp->epsilon = sv_dline_new_alike(jp->ee, 1);
    sv_dline_fill(jp->epsilon, 1);
    jp->mu = sv_dline_new_alike(jp->ee, 1);
    sv_dline_fill(jp->mu, 1);
    
    jp->sigma = sv_dline_new_alike(jp->ee, 1);
    jp->sigast = sv_dline_new_alike(jp->ee, 1);
    jp->source = sv_dline_new_alike(jp->ee, 1);

    jp->step = 0;
    jp->dx = dx;
    jp->corr = 1;
    
    jp->fbnx0_e_0=0;
    jp->fbnx0_e_1=0;
    jp->fbnx0_h_0=0;
    jp->fbnx0_h_1=0;
    jp->fbnxn_e_0=0;
    jp->fbnxn_e_1=0;
    jp->fbnxn_h_0=0;
    jp->fbnxn_h_1=0;
  
    jp->fbnx0_ee_0=0;
    jp->fbnx0_ee_1=0;
    jp->fbnx0_he_0=0;
    jp->fbnx0_he_1=0;
    jp->fbnxn_ee_0=0;
    jp->fbnxn_ee_1=0;
    jp->fbnxn_he_0=0;
    jp->fbnxn_he_1=0;
  
    jp->fbnx0_eh_0=0;
    jp->fbnx0_eh_1=0;
    jp->fbnx0_hh_0=0;
    jp->fbnx0_hh_1=0;
    jp->fbnxn_eh_0=0;
    jp->fbnxn_eh_1=0;
    jp->fbnxn_hh_0=0;
    jp->fbnxn_hh_1=0;
     
    return jp;
}

void
sv_zpool_set_source_pos(SvZPool *jp, gint i)
{
    jp->source->data[i] = 1;
}

void
sv_zpool_set_correction(SvZPool *jp, gdouble corr)
{
    jp->corr = corr;
}

void
sv_zpool_set_material(SvZPool *jp, gint nlayers, gint *lpos, gdouble *epsilon, gdouble *mu, gdouble *sigma, gdouble *sigast)
{
    gint i, k;

    for (k = 0; k < nlayers; k++) {
        for (i = lpos[k]; i < jp->res; i++) {
            jp->mu->data[i] = mu[k];
            jp->epsilon->data[i] = epsilon[k];
            jp->sigma->data[i] = sigma[k];
            jp->sigast->data[i] = sigast[k];
	}
    }
}

void
sv_zpool_set_angle(SvZPool *jp, gdouble theta, gdouble psi)
{
    jp->theta = theta;
    jp->psi = psi;
}

void
sv_zpool_store(SvZPool *jp, gint pos)
{
    gint i;
    gdouble *eedata = gwy_data_field_get_data(jp->dee);
    gdouble *ehdata = gwy_data_field_get_data(jp->deh);
    gdouble *ezdata = gwy_data_field_get_data(jp->dez);
    gdouble *hedata = gwy_data_field_get_data(jp->dhe);
    gdouble *hhdata = gwy_data_field_get_data(jp->dhh);
    gdouble *hzdata = gwy_data_field_get_data(jp->dhz);

    for (i = 0; i < jp->res; i++) {
        eedata[pos*jp->res + i] = jp->ee->data[i];
        ehdata[pos*jp->res + i] = jp->eh->data[i];
        ezdata[pos*jp->res + i] = jp->ez->data[i];
        hedata[pos*jp->res + i] = jp->he->data[i];
        hhdata[pos*jp->res + i] = jp->hh->data[i];
        hzdata[pos*jp->res + i] = jp->hz->data[i];


    }      
     /*  sprintf(buff, "src%.3d.txt", pos); 
       fr = fopen(buff, "w");

       for (i=0; i<jp->res; i++)
       {
       fprintf(fr, "%d %g    %g %g %g   %g %g %g\n", i, jp->theta, eedata[pos*jp->res + i], ehdata[pos*jp->res + i], ezdata[pos*jp->res + i], hedata[pos*jp->res + i], hhdata[pos*jp->res + i], hzdata[pos*jp->res + i]);
       }
       fclose(fr);
    */
}

void
sv_zpool_ystep_e(SvZPool *jp)
{
    gint i;
    
    gdouble angle_mult = jp->corr; //for dAngle=0 or 90 FIXME for others
    gdouble multe = jp->dt / (EPSILON_0*jp->dx) / angle_mult;

   for (i = 0; i < jp->res; i++) { //x
        /*update source for 1d tsf, vacuum*/
        jp->e->data[i]+= 1.0 * multe * (jp->h->data[i-1] - jp->h->data[i]);
    }

 

    for (i = 0; i < jp->res; i++) { //x
        /*update splitted fields*/
        if (i > 0) {
            jp->ee->data[i] = jp->ee->data[i] + 1.0/jp->epsilon->data[i] * multe * (jp->he->data[i-1] - jp->he->data[i]);
            jp->eh->data[i] = jp->eh->data[i] + 1.0/(jp->epsilon->data[i] - sin(jp->theta)*sin(jp->theta)) * multe * (jp->hh->data[i-1] - jp->hh->data[i]);

            jp->ez->data[i] = -jp->ez->data[i] - 2*sqrt(MU_0/EPSILON_0)*sin(jp->theta) / jp->epsilon->data[i] * jp->he->data[i];
        }
    }
    jp->ee->data[10] -= jp->dt / jp->dx / (EPSILON_0*jp->epsilon->data[10]) * jp->h->data[9] * sin(jp->psi); //inverted to match tsf
    jp->eh->data[10] -= jp->dt / jp->dx / (EPSILON_0*jp->epsilon->data[10]) * jp->h->data[9] * cos(jp->psi); //inverted to match tsf

    //printf("------ %g %g\n", jp->dt/jp->dx/EPSILON_0*jp->h->data[9]*sin(jp->psi), jp->dt/jp->dx/EPSILON_0*jp->h->data[9]*cos(jp->psi));
    jp->step++;
}

void
sv_zpool_ystep_h(SvZPool *jp)
{
    gint i;

    gdouble angle_mult = jp->corr;
    gdouble multh = jp->dt / (MU_0*jp->dx) / angle_mult;

    for (i=0; i<jp->res; i++) { //x
        /*update source for 1d tsf, vacuum*/
        jp->h->data[i]+= 1.0*multh * (jp->e->data[i] - jp->e->data[i+1]);
    }

 
    for (i=0; i<jp->res; i++) { //x
        /*update splited fields*/
        if (i < (jp->res - 1)) {
            jp->he->data[i] = jp->he->data[i] + 1.0/jp->mu->data[i] * jp->epsilon->data[i]/(jp->epsilon->data[i]-sin(jp->theta)*sin(jp->theta)) *
		              multh * (jp->ee->data[i] - jp->ee->data[i+1]);
            jp->hh->data[i] = jp->hh->data[i] + 1.0/jp->mu->data[i] * multh * (jp->eh->data[i] - jp->eh->data[i+1]);
            
            jp->hz->data[i] = -jp->hz->data[i] + 2*sqrt(EPSILON_0/MU_0) * sin(jp->theta) * jp->eh->data[i];
        }
    }
    jp->he->data[9] -= jp->dt/jp->dx/MU_0*jp->e->data[10]*sin(jp->psi);
    jp->hh->data[9] -= jp->dt/jp->dx/MU_0*jp->e->data[10]*cos(jp->psi); 

/*
    if (jp->step == 1) {
        fm = fopen("material.txt", "w");
        for (i = 0; i < jp->res; i++)
	    fprintf(fm, "%d %g %g %g %g\n", i, jp->epsilon->data[i], jp->sigma->data[i], jp->mu->data[i], jp->sigast->data[i]);
        fclose(fm);
    }
*/
}

void
sv_zpool_boundary_copy(SvZPool *jp)
{
    jp->fbnx0_e_0 = jp->e->data[0];
    jp->fbnx0_e_1 = jp->e->data[1];
    jp->fbnx0_h_0 = jp->h->data[0];
    jp->fbnx0_h_1 = jp->h->data[1];
    jp->fbnxn_e_0 = jp->e->data[jp->res-1];
    jp->fbnxn_e_1 = jp->e->data[jp->res-2];
    jp->fbnxn_h_0 = jp->h->data[jp->res-1];
    jp->fbnxn_h_1 = jp->h->data[jp->res-2];
 
    jp->fbnx0_ee_0 = jp->ee->data[0];
    jp->fbnx0_ee_1 = jp->ee->data[1];
    jp->fbnx0_he_0 = jp->he->data[0];
    jp->fbnx0_he_1 = jp->he->data[1];
    jp->fbnxn_ee_0 = jp->ee->data[jp->res-1];
    jp->fbnxn_ee_1 = jp->ee->data[jp->res-2];
    jp->fbnxn_he_0 = jp->he->data[jp->res-1];
    jp->fbnxn_he_1 = jp->he->data[jp->res-2];
  
    jp->fbnx0_eh_0 = jp->eh->data[0];
    jp->fbnx0_eh_1 = jp->eh->data[1];
    jp->fbnx0_hh_0 = jp->hh->data[0];
    jp->fbnx0_hh_1 = jp->hh->data[1];
    jp->fbnxn_eh_0 = jp->eh->data[jp->res-1];
    jp->fbnxn_eh_1 = jp->eh->data[jp->res-2];
    jp->fbnxn_hh_0 = jp->hh->data[jp->res-1];
    jp->fbnxn_hh_1 = jp->hh->data[jp->res-2];
}

void
sv_zpool_boundary(SvZPool *jp)
{

    gdouble ind0 = sqrt(jp->epsilon->data[0]);
    gdouble ind1 = sqrt(jp->epsilon->data[jp->res-1]);

    jp->e->data[0] = jp->fbnx0_e_1 + (jp->dt*LIGHT_SPEED/ind0-jp->dx)/(jp->dt*LIGHT_SPEED/ind0+jp->dx) * (jp->e->data[1] - jp->fbnx0_e_0);
//    jp->h->data[0] = jp->fbnx0_h_1 + (jp->dt*LIGHT_SPEED/ind0-jp->dx)/(jp->dt*LIGHT_SPEED/ind0+jp->dx) * (jp->h->data[1] - jp->fbnx0_h_0);

    jp->e->data[jp->res-1] = jp->fbnxn_e_1 + (jp->dt*LIGHT_SPEED/ind1-jp->dx)/(jp->dt*LIGHT_SPEED/ind1+jp->dx) * (jp->e->data[jp->res-2] - jp->fbnxn_e_0);
//    jp->h->data[jp->res-1] = jp->fbnxn_h_1 + (jp->dt*LIGHT_SPEED/ind1-jp->dx)/(jp->dt*LIGHT_SPEED/ind1+jp->dx) * (jp->h->data[jp->res-2] - jp->fbnxn_h_0);
 
    jp->ee->data[0] = jp->fbnx0_ee_1 + (jp->dt*LIGHT_SPEED/ind0-jp->dx)/(jp->dt*LIGHT_SPEED/ind0+jp->dx) * (jp->ee->data[1] - jp->fbnx0_ee_0);
 //   jp->he->data[0] = jp->fbnx0_he_1 + (jp->dt*LIGHT_SPEED/ind0-jp->dx)/(jp->dt*LIGHT_SPEED/ind0+jp->dx) * (jp->he->data[1] - jp->fbnx0_he_0);

    jp->ee->data[jp->res-1] = jp->fbnxn_ee_1 + (jp->dt*LIGHT_SPEED/ind1-jp->dx)/(jp->dt*LIGHT_SPEED/ind1+jp->dx)*(jp->ee->data[jp->res-2] - jp->fbnxn_ee_0);
 //   jp->he->data[jp->res-1] = jp->fbnxn_he_1 + (jp->dt*LIGHT_SPEED/ind1-jp->dx)/(jp->dt*LIGHT_SPEED/ind1+jp->dx)*(jp->he->data[jp->res-2] - jp->fbnxn_he_0);

    jp->eh->data[0] = jp->fbnx0_eh_1 + (jp->dt*LIGHT_SPEED/ind0-jp->dx)/(jp->dt*LIGHT_SPEED/ind0+jp->dx) * (jp->eh->data[1] - jp->fbnx0_eh_0);
 //   jp->hh->data[0] = jp->fbnx0_hh_1 + (jp->dt*LIGHT_SPEED/ind0-jp->dx)/(jp->dt*LIGHT_SPEED/ind0+jp->dx) * (jp->hh->data[1] - jp->fbnx0_hh_0);

    jp->eh->data[jp->res-1] = jp->fbnxn_eh_1 + (jp->dt*LIGHT_SPEED/ind1-jp->dx)/(jp->dt*LIGHT_SPEED/ind1+jp->dx)*(jp->eh->data[jp->res-2] - jp->fbnxn_eh_0);
 //   jp->hh->data[jp->res-1] = jp->fbnxn_hh_1 + (jp->dt*LIGHT_SPEED/ind1-jp->dx)/(jp->dt*LIGHT_SPEED/ind1+jp->dx)*(jp->hh->data[jp->res-2] - jp->fbnxn_hh_0);
}

void
sv_zpool_boundary_estep(SvZPool *jp)
{

    gdouble ind0 = sqrt(jp->epsilon->data[0]);
    gdouble ind1 = sqrt(jp->epsilon->data[jp->res-1]);

    jp->e->data[0] = jp->fbnx0_e_1 + (jp->dt*LIGHT_SPEED/ind0-jp->dx)/(jp->dt*LIGHT_SPEED/ind0+jp->dx) * (jp->e->data[1] - jp->fbnx0_e_0);
    jp->e->data[jp->res-1] = jp->fbnxn_e_1 + (jp->dt*LIGHT_SPEED/ind1-jp->dx)/(jp->dt*LIGHT_SPEED/ind1+jp->dx) * (jp->e->data[jp->res-2] - jp->fbnxn_e_0);
    jp->ee->data[0] = jp->fbnx0_ee_1 + (jp->dt*LIGHT_SPEED/ind0-jp->dx)/(jp->dt*LIGHT_SPEED/ind0+jp->dx) * (jp->ee->data[1] - jp->fbnx0_ee_0);
    jp->ee->data[jp->res-1] = jp->fbnxn_ee_1 + (jp->dt*LIGHT_SPEED/ind1-jp->dx)/(jp->dt*LIGHT_SPEED/ind1+jp->dx)*(jp->ee->data[jp->res-2] - jp->fbnxn_ee_0);
    jp->eh->data[0] = jp->fbnx0_eh_1 + (jp->dt*LIGHT_SPEED/ind0-jp->dx)/(jp->dt*LIGHT_SPEED/ind0+jp->dx) * (jp->eh->data[1] - jp->fbnx0_eh_0);
    jp->eh->data[jp->res-1] = jp->fbnxn_eh_1 + (jp->dt*LIGHT_SPEED/ind1-jp->dx)/(jp->dt*LIGHT_SPEED/ind1+jp->dx)*(jp->eh->data[jp->res-2] - jp->fbnxn_eh_0);
}

void
sv_zpool_boundary_hstep(SvZPool *jp)
{

    gdouble ind0 = sqrt(jp->epsilon->data[0]);
    gdouble ind1 = sqrt(jp->epsilon->data[jp->res-1]);

    jp->h->data[0] = jp->fbnx0_h_1 + (jp->dt*LIGHT_SPEED/ind0-jp->dx)/(jp->dt*LIGHT_SPEED/ind0+jp->dx) * (jp->h->data[1] - jp->fbnx0_h_0);
    jp->h->data[jp->res-1] = jp->fbnxn_h_1 + (jp->dt*LIGHT_SPEED/ind1-jp->dx)/(jp->dt*LIGHT_SPEED/ind1+jp->dx) * (jp->h->data[jp->res-2] - jp->fbnxn_h_0);
    jp->he->data[0] = jp->fbnx0_he_1 + (jp->dt*LIGHT_SPEED/ind0-jp->dx)/(jp->dt*LIGHT_SPEED/ind0+jp->dx) * (jp->he->data[1] - jp->fbnx0_he_0);
    jp->he->data[jp->res-1] = jp->fbnxn_he_1 + (jp->dt*LIGHT_SPEED/ind1-jp->dx)/(jp->dt*LIGHT_SPEED/ind1+jp->dx)*(jp->he->data[jp->res-2] - jp->fbnxn_he_0);
    jp->hh->data[0] = jp->fbnx0_hh_1 + (jp->dt*LIGHT_SPEED/ind0-jp->dx)/(jp->dt*LIGHT_SPEED/ind0+jp->dx) * (jp->hh->data[1] - jp->fbnx0_hh_0);
    jp->hh->data[jp->res-1] = jp->fbnxn_hh_1 + (jp->dt*LIGHT_SPEED/ind1-jp->dx)/(jp->dt*LIGHT_SPEED/ind1+jp->dx)*(jp->hh->data[jp->res-2] - jp->fbnxn_hh_0);
}

void
sv_zpool_apply_source(SvZPool *jp, gdouble eval)
{
    gint i;
    
    /*for each point marked as the source apply source*/
    for (i = 0; i < jp->res; i++) {
        if (jp->source->data[i] == 1)
            jp->e->data[i] = eval;
    }
}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
