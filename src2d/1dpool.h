
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


/*  1dpool.h : 
 *  one-dimensional FDTD for use in TSF and SF sources
 */

#ifndef SV_JPOOL
#define SV_JPOOL

#include "dline.h"
#include <libprocess/gwyprocess.h>
#include "settings.h"


typedef struct _Sv1DPool {
   SvDLine *e;        /*basic data arrays*/
   SvDLine *h;

   SvDLine *epsilon;   /*detailed material properties*/
   SvDLine *mu;
   SvDLine *sigma;
   SvDLine *sigast;

   SvDLine *source;

   gdouble fbnx0_e_0;  //previous value at the boundary
   gdouble fbnx0_h_0;
   gdouble fbnx0_e_1;  //previous value 1 step before boundary
   gdouble fbnx0_h_1;
   gdouble fbnxn_e_0;
   gdouble fbnxn_h_0;
   gdouble fbnxn_e_1;
   gdouble fbnxn_h_1;

   /*boundary parameters*/
   gint bx0; //status of the x=0 plane boundary
   gint bxn; //status of the x=max plane boundary

   gdouble dx;
   gint res;
   gdouble dt; 
   gint step;
   
} Sv1DPool;

Sv1DPool * sv_1dpool_new(gint res, gdouble dx, gdouble dt);
void sv_1dpool_set_source_pos(Sv1DPool *jp, gint i);

void sv_1dpool_ystep_e(Sv1DPool *jp); 
void sv_1dpool_ystep_h(Sv1DPool *jp);
void sv_1dpool_apply_source(Sv1DPool *jp, gdouble eval);
void sv_1dpool_boundary(Sv1DPool *jp);
void sv_1dpool_boundary_copy(Sv1DPool *jp);
void sv_1dpool_set_material(Sv1DPool *jp, gdouble epsilon, gdouble mu);


#endif
/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

