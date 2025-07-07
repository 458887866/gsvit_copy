
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


/*  farfield.h : 
 *  near-to-far-field calculation algorithm based on Kirchoff surface integral
 */


#ifndef SV_NFFF
#define SV_NFFF

#include <glib.h>
#include <libprocess/gwyprocess.h>
#include "settings.h"
#include "1dpool.h"

typedef struct
{
    GwyDataField *tm_left_hyp;
    GwyDataField *tm_left_hym;
    GwyDataField *tm_left_ez;

    GwyDataField *tm_right_hyp;
    GwyDataField *tm_right_hym;
    GwyDataField *tm_right_ez;

    GwyDataField *tm_top_hxp;
    GwyDataField *tm_top_hxm;
    GwyDataField *tm_top_ez;

    GwyDataField *tm_bottom_hxp;
    GwyDataField *tm_bottom_hxm;
    GwyDataField *tm_bottom_ez;

} SvNFFFStore;


typedef struct
{
    SvNFFFStore *nfs;
} SvFarfield;

SvFarfield* sv_farfield_new();


#endif
/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
