
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


/*  source.h : 
 *  data for computation of elmag field sources
 */


#ifndef SV_SOURCE
#define SV_SOURCE

#include <glib.h>
#include "1dpool.h"

typedef struct
{
    gint ndata;
    gint *pos;
    gdouble* ex;
    gdouble* ey;
    gdouble* ez;
    gdouble* hx;
    gdouble* hy;
    gdouble* hz;
} SvSrcData;

typedef struct
{
    gint i;
    gint j;
    SvSrcData sdata;
} SvSourcePoint;

typedef struct
{
    Sv1DPool *jp;
    gint i0;
    gint j0;
    gint i1;
    gint j1;
    gdouble theta;
    gdouble phi;
    gdouble psi;
    gdouble *e;
    gint *pos;
    gdouble epsilon;
    gdouble sigma;
    gdouble mu;
    gdouble sigast;
    gdouble corr;
    int ndata;

} SvSourceTSF;

typedef struct
{
    SvSourcePoint *sp;
    SvSourceTSF *tsf;
} SvSource;

SvSource* sv_source_new();


#endif
/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
