
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


/*  output.h : 
 *  all the output functions
 */

#ifndef OUTPUT_H
#define OUTPUT_H

#include <glib.h>
#include <libprocess/gwyprocess.h>
#include <libgwyddion/gwycontainer.h>
#include "pool.h"


typedef struct
{
    gdouble **outsumdata;
    gdouble **outforcedata;
    gdouble **outpointdata;

    SvDCube **outabsdata;

    gint last;                 /*last output datafield key*/
    gint glast;                /*last output dataline key*/
    GwyContainer *data;
} SvOutput;

SvOutput* sv_output_new();


void output_ivtk(int ***data, int xres, int yres, int zres, char *filename, char *description);


#endif /* OUTPUT_H */
/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
