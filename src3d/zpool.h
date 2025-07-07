
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


/*  zpool.h : 
 *  one-dimensional FDTD for use in TSF and SF sources
 */

#ifndef ZPOOL_H
#define ZPOOL_H

#include "dline.h"
#include <libprocess/gwyprocess.h>
#include "settings.h"


typedef struct _SvZPool {
   SvDLine *hh;        /*basic data arrays*/
   SvDLine *he;
   SvDLine *eh;
   SvDLine *ee;
   SvDLine *ez;
   SvDLine *hz;

   GwyDataField  *deh; /*all fields history*/
   GwyDataField  *dee;
   GwyDataField  *dez;
   GwyDataField  *dhh;
   GwyDataField  *dhe;
   GwyDataField  *dhz;

   SvDLine *e;
   SvDLine *h;

   SvDLine *epsilon;   /*detailed material properties*/
   SvDLine *mu;
   SvDLine *sigma;
   SvDLine *sigast;

   SvDLine *source;

   gdouble fbnx0_ee_0;  //previous value at the boundary
   gdouble fbnx0_he_0;
   gdouble fbnx0_ee_1;  //previous value 1 step before boundary
   gdouble fbnx0_he_1;
   gdouble fbnxn_ee_0;
   gdouble fbnxn_he_0;
   gdouble fbnxn_ee_1;
   gdouble fbnxn_he_1;

   gdouble fbnx0_eh_0;  //previous value at the boundary
   gdouble fbnx0_hh_0;
   gdouble fbnx0_eh_1;  //previous value 1 step before boundary
   gdouble fbnx0_hh_1;
   gdouble fbnxn_eh_0;
   gdouble fbnxn_hh_0;
   gdouble fbnxn_eh_1;
   gdouble fbnxn_hh_1;

   gdouble fbnx0_e_0;  //previous value at the boundary
   gdouble fbnx0_h_0;
   gdouble fbnxn_e_0;
   gdouble fbnxn_h_0;
   gdouble fbnx0_e_1;  //previous value 1 step before boundary
   gdouble fbnx0_h_1;
   gdouble fbnxn_e_1;
   gdouble fbnxn_h_1;


   /*boundary parameters*/
   gint bx0; //status of the x=0 plane boundary
   gint bxn; //status of the x=max plane boundary

   gdouble dx;
   gint res;
   gdouble dt; 
   gint step;
   gint nsteps;
   gdouble theta;
   gdouble psi;
   gdouble corr;
   
} SvZPool;

SvZPool * sv_zpool_new(gint res, gdouble dx, gdouble dt, gint nsteps);
void sv_zpool_set_source_pos(SvZPool *jp, gint i);

void sv_zpool_ystep_e(SvZPool *jp); 
void sv_zpool_ystep_h(SvZPool *jp);
void sv_zpool_apply_source(SvZPool *jp, gdouble eval);
void sv_zpool_boundary(SvZPool *jp);
void sv_zpool_boundary_estep(SvZPool *jp);
void sv_zpool_boundary_hstep(SvZPool *jp);
void sv_zpool_boundary_copy(SvZPool *jp);
void sv_zpool_set_material(SvZPool *jp, gint nlayers, gint *lpos, gdouble *epsilon, gdouble *mu, gdouble *sigma, gdouble *sigast);
void sv_zpool_set_angle(SvZPool *jp, gdouble theta, gdouble psi);
void sv_zpool_set_correction(SvZPool *jp, gdouble corr);
void sv_zpool_store(SvZPool *jp, gint pos);

#endif /* ZPOOL_H */

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

