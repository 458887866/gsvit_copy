
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


/*  modifiers.h : 
 *  header file for material modifiers
 */

#ifndef MODIFIERS_H
#define MODIFIERS_H

#include <glib.h>
#include "dcube.h"
#include "pool.h"

/*alters boundaries position by random shifts in 3d*/
void modifier_fftshift(SvICube *mat, gint material, gdouble xsigma, gdouble xt, gdouble ysigma, gdouble yt, gdouble zsigma, gdouble zt, gint seed);

/*alters boundaries position by random shifts in 3d, experimental antialiased version*/
SvMatProp * modifier_fftshift_antialiased(SvICube *mat, SvMatProp *mats, gint *nmat, gint material, gdouble xsigma, gdouble xt, gdouble ysigma, gdouble yt, gdouble zsigma, gdouble zt, gint seed, gint aa);


#endif /* MODIFIERS_H */
/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
