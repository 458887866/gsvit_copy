
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


/*  dline.h : 
 *  simple 1D data representation
 *  FIXME: Should be replaced by GwyDataLine
 */

#ifndef DLINE_H
#define DLINE_H

#include <glib.h>
#include <libprocess/gwyprocess.h>

typedef struct {
   gint res;
   gdouble real;
   gdouble *data;
} SvDLine;

SvDLine *sv_dline_new(gint res, 
                      gdouble real,
                      gboolean nullme);
SvDLine *sv_dline_new_alike   (SvDLine *dline, gboolean nullme);
void sv_dline_free            (SvDLine *dline);
void sv_dline_fill            (SvDLine *dline, gdouble value);
void sv_dline_add             (SvDLine *dline, gdouble value);
void sv_dline_multiply        (SvDLine *dline, gdouble value);
void sv_dline_get_dataline    (SvDLine *dline, GwyDataLine *gwydline);
gdouble sv_dline_get_dval     (SvDLine *dline, gdouble pos);

#endif /* DLINE_H */

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

