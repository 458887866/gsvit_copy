
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


/*  dcube.h : 
 *  Simple 3D data representation.
 *  Should be replaced by GwyBrick in future.
 */



#ifndef SV_DCUBE
#define SV_DCUBE

#include <glib.h>
#include <libprocess/gwyprocess.h>

typedef struct {
   gint xres;
   gint yres;
   gint zres;
   gdouble xreal;
   gdouble yreal;
   gdouble zreal;
   gdouble ***data;
} SvDCube;

SvDCube *sv_dcube_new(gint xres, gint yres, gint zres,
                      gdouble xreal, gdouble yreal, gdouble zreal,
                      gboolean nullme);
SvDCube *sv_dcube_new_alike   (SvDCube *dcube, gboolean nullme);
void sv_dcube_free            (SvDCube *dcube);
void sv_dcube_fill            (SvDCube *dcube, gdouble value);
gdouble sv_dcube_get_sum      (SvDCube *dcube);
void sv_dcube_add             (SvDCube *dcube, gdouble value);
void sv_dcube_multiply        (SvDCube *dcube, gdouble value);
void sv_dcube_get_datafield   (SvDCube *dcube, GwyDataField *dfield, gint i, gint j, gint k);

#endif

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

