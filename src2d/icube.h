
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


/*  icube.h : 
 *  3D volume of integer values
 */

#ifndef SV_ICUBE
#define SV_ICUBE

#include <glib.h>
#include <libprocess/gwyprocess.h>

typedef struct {
   gint xres;
   gint yres;
   gint zres;
   gdouble xreal;
   gdouble yreal;
   gdouble zreal;
   gint ***data;
} SvICube;

SvICube *sv_icube_new(gint xres, gint yres, gint zres,
                      gdouble xreal, gdouble yreal, gdouble zreal,
                      gboolean nullme);

void        sv_icube_free       (SvICube *icube);
void        sv_icube_fill       (SvICube *icube, gint value);
void        sv_icube_add        (SvICube *icube, gint value);
void        sv_icube_multiply   (SvICube *icube, gint value);
void        sv_icube_get_datafield (SvICube *icube, GwyDataField *dfield, gint ipos, gint jpos, gint kpos);

#endif
/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

