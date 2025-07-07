
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


/*  subgrid.h : 
 *  subgrid implementation
 */

#ifndef SUBGRID_H
#define SUBGRID_H

#define MAXSUBGRID 10

#include <glib.h>
#include "yee.h"

/*Subgrid data*/
typedef struct{
    
    gint ifrom;
    gint jfrom;
    gint kfrom;
    gint ito;
    gint jto;
    gint kto;
    gint division;

    /*all the local computational domain arrays*/
    SvYeeData *d;

    /*previous h field values*/
    gint phxres;
    gint phyres;
    gint phzres;

    GwyDataField *px0hy; //step before
    GwyDataField *px0hz;
    GwyDataField *py0hx;
    GwyDataField *py0hz;
    GwyDataField *pz0hx;
    GwyDataField *pz0hy;

    GwyDataField *ppx0hy; //two steps before
    GwyDataField *ppx0hz;
    GwyDataField *ppy0hx;
    GwyDataField *ppy0hz;
    GwyDataField *ppz0hx;
    GwyDataField *ppz0hy;

    GwyDataField *ax0hy; //actual value
    GwyDataField *ax0hz;
    GwyDataField *ay0hx;
    GwyDataField *ay0hz;
    GwyDataField *az0hx;
    GwyDataField *az0hy;

    GwyDataField *nx0hy; //instantenous updates
    GwyDataField *nx0hz;
    GwyDataField *ny0hx;
    GwyDataField *ny0hz;
    GwyDataField *nz0hx;
    GwyDataField *nz0hy;


    GwyDataField *pxnhy; //step before
    GwyDataField *pxnhz;
    GwyDataField *pynhx;
    GwyDataField *pynhz;
    GwyDataField *pznhx;
    GwyDataField *pznhy;

    GwyDataField *ppxnhy; //two steps before
    GwyDataField *ppxnhz;
    GwyDataField *ppynhx;
    GwyDataField *ppynhz;
    GwyDataField *ppznhx;
    GwyDataField *ppznhy;

    GwyDataField *axnhy; //actual value
    GwyDataField *axnhz;
    GwyDataField *aynhx;
    GwyDataField *aynhz;
    GwyDataField *aznhx;
    GwyDataField *aznhy;


    GwyDataField *nxnhy; //instantenous updates
    GwyDataField *nxnhz;
    GwyDataField *nynhx;
    GwyDataField *nynhz;
    GwyDataField *nznhx;
    GwyDataField *nznhy;


} SvSg;

/*all the subgrids*/
typedef struct
{
   SvSg **sg;
   int nsubgrid;

} SvSubgrid;




#endif /* BOUNDARY_H */
/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
