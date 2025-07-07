
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


/*  pfarfield.h :
 *  periodic near-to-far field transformation
 */

#ifndef PFARFIELD_H
#define PFARFIELD_H

#include <glib.h>
#include <libprocess/gwyprocess.h>
#include <libgwyddion/gwycontainer.h>
#include "dcube.h"


typedef struct
{
    GwyDataLine *ex;
    GwyDataLine *ey;
    GwyDataLine *ez;
    gint istart;

} SvPRPoint;

typedef struct
{
    SvPRPoint *prpoints;
    gint ndata;

    SvDCube *ex_k0;
    SvDCube *ey_k0;
    SvDCube *ez_k0;

    SvDCube *ex_kn;
    SvDCube *ey_kn;
    SvDCube *ez_kn;

    SvDCube *dex_k0;
    SvDCube *dey_k0;
    SvDCube *dez_k0;

    SvDCube *dex_kn;
    SvDCube *dey_kn;
    SvDCube *dez_kn;



} SvPFarfield;


#endif /* PFARFIELD_H */

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
