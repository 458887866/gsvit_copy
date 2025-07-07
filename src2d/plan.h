
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


/*  plan.h : 
 *  create plan for computation, load data files
 */


#ifndef SV_PLAN
#define SV_PLAN

#include "settings.h"
#include "pool.h"


typedef struct {
    double pnt1[2];
    double pnt2[2];
    double pnt3[2];
    SvMatProp mat;
    int n;
} SvTriangle;

typedef struct {
    double pnt1[2];
    double pnt2[2];
    SvMatProp mat;
    int n;
} SvRectangle;

typedef struct {
    double pnt1[2];
    double radius;
    SvMatProp mat;
    int n;
} SvCircle;

SvPool* init_and_make_plan(SvSet *set);


#endif
/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
