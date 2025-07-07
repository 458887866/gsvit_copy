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

#ifndef SV_BND
#define SV_BND

#include <glib.h>
#include "dcube.h"

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

#include "settings.h"

typedef struct{
    SvDCube *fbnx;
    SvDCube *fbny;
} SvLiaoBc;

typedef struct
{
   SvLiaoBc liao;
} SvBoundary;

SvBoundary* sv_boundary_new(SvSet *set);



#endif
/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
