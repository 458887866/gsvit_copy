
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


/*  boundary.h : 
 *  boundary conditions implementation
 */

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <glib.h>
#include "dcube.h"

/*positions in stored boundary condition data*/
#define EX00 0
#define EY00 1
#define EZ00 2
#define HX00 3
#define HY00 4
#define HZ00 5
#define EXN0 6
#define EYN0 7
#define EZN0 8
#define HXN0 9
#define HYN0 10
#define HZN0 11

#define EX0P 12
#define EY0P 13
#define EZ0P 14
#define HX0P 15
#define HY0P 16
#define HZ0P 17
#define EXNP 18
#define EYNP 19
#define EZNP 20
#define HXNP 21
#define HYNP 22
#define HZNP 23

#define NPLANES 24


/*Liao BC data*/
typedef struct{
    SvDCube *fbnx;
    SvDCube *fbny;
    SvDCube *fbnz;
} SvLiaoBc;

/*CPML data*/
typedef struct{
    
    SvDCube *peyx_x0;
    SvDCube *pezx_x0;

    SvDCube *peyx_xn;
    SvDCube *pezx_xn;

    SvDCube *phzx_x0;
    SvDCube *phyx_x0;

    SvDCube *phzx_xn;
    SvDCube *phyx_xn;


    gdouble *be_x0;
    gdouble *ce_x0;
    gdouble *bh_x0;
    gdouble *ch_x0;
    gdouble *kappae_x0;
    gdouble *kappah_x0;

    gdouble *be_xn;
    gdouble *ce_xn;
    gdouble *bh_xn;
    gdouble *ch_xn;
    gdouble *kappae_xn;
    gdouble *kappah_xn;


    SvDCube *pexy_y0;
    SvDCube *pezy_y0;

    SvDCube *pexy_yn;
    SvDCube *pezy_yn;

    SvDCube *phxy_y0;
    SvDCube *phzy_y0;

    SvDCube *phxy_yn;
    SvDCube *phzy_yn;


    gdouble *be_y0;
    gdouble *ce_y0;
    gdouble *bh_y0;
    gdouble *ch_y0;
    gdouble *kappae_y0;
    gdouble *kappah_y0;

    gdouble *be_yn;
    gdouble *ce_yn;
    gdouble *bh_yn;
    gdouble *ch_yn;
    gdouble *kappae_yn;
    gdouble *kappah_yn;

    SvDCube *peyz_z0;
    SvDCube *pexz_z0;

    SvDCube *peyz_zn;
    SvDCube *pexz_zn;

    SvDCube *phyz_z0;
    SvDCube *phxz_z0;

    SvDCube *phyz_zn;
    SvDCube *phxz_zn;

    gdouble *be_z0;
    gdouble *ce_z0;
    gdouble *bh_z0;
    gdouble *ch_z0;
    gdouble *kappae_z0;
    gdouble *kappah_z0;

    gdouble *be_zn;
    gdouble *ce_zn;
    gdouble *bh_zn;
    gdouble *ch_zn;
    gdouble *kappae_zn;
    gdouble *kappah_zn;
} SvCpmlBc;

typedef struct
{
   SvLiaoBc liao;
   SvCpmlBc cpml;

} SvBoundary;




#endif /* BOUNDARY_H */
/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
