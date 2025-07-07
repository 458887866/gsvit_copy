
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


/*  icube.c : 
 *  3D volume of integers
 */

#include <stdio.h>
#include "icube.h"


SvICube *sv_icube_new(gint xres, gint yres, gint zres,
                      gdouble xreal, gdouble yreal, gdouble zreal,
                      gboolean nullme)
{
    gint i, j;
    SvICube *ic = (SvICube *)g_malloc(sizeof(SvICube));

    ic->xres = xres;
    ic->yres = yres; 
    ic->zres = zres;
    ic->xreal = xreal;
    ic->yreal = yreal;
    ic->zreal = zreal;
    ic->data= (gint ***) g_malloc(xres*sizeof(gint**));
    for (i=0; i<xres; i++)
    {
        ic->data[i]= (gint **) g_malloc(yres*sizeof(gint*));
        for (j=0; j<yres; j++)
        {
            ic->data[i][j]= (gint *) g_malloc(zres*sizeof(gint));
        }
    }
    if (nullme) sv_icube_fill(ic, 0);
    return ic; 
}

void sv_icube_free(SvICube *icube)
{
    gint i, j;
    for (i=0; i<icube->xres; i++)
    {
        for (j=0; j<icube->yres; j++)
        {
            g_free((void **) icube->data[i][j]);
        }
        g_free((void *) icube->data[i]);
    }

    g_free((void *)icube->data);
    icube->data=NULL;

}

void sv_icube_fill(SvICube *icube, gint value)
{
    gint i, j, k;
//#pragma omp parallel default(shared) private(j, k)
    for (i=0; i<icube->xres; i++)
    {
        for (j=0; j<icube->yres; j++)
        {
            for (k=0; k<icube->zres; k++)
                icube->data[i][j][k] = value;
        }
    }
}

void sv_icube_add(SvICube *icube, gint value)
{
    gint i, j, k;
//#pragma omp parallel default(shared) private(j, k)
    for (i=0; i<icube->xres; i++)
    {
        for (j=0; j<icube->yres; j++)
        {
            for (k=0; k<icube->zres; k++)
                icube->data[i][j][k] += value;
        }
    }
}

void sv_icube_multiply(SvICube *icube, gint value)
{
    gint i, j, k;
//#pragma omp parallel default(shared) private(j, k)
    for (i=0; i<icube->xres; i++)
    {
        for (j=0; j<icube->yres; j++)
        {
            for (k=0; k<icube->zres; k++)
                icube->data[i][j][k] *= value;
        }
    }
}

void sv_icube_get_datafield   (SvICube *icube, GwyDataField *dfield, gint ipos, gint jpos, gint kpos)
{
    gint i, j, k;
    if (ipos==-1 && jpos==-1)
    {
        gwy_data_field_resample(dfield, icube->xres, icube->yres, GWY_INTERPOLATION_NONE);
        for (i=0; i<icube->xres; i++)
        {
            for (j=0; j<icube->yres; j++)
                gwy_data_field_set_val(dfield, i, j, icube->data[i][j][kpos]);
        }
    }
    else if (ipos==-1 && kpos==-1)
    {
        gwy_data_field_resample(dfield, icube->xres, icube->zres, GWY_INTERPOLATION_NONE);
        for (i=0; i<icube->xres; i++)
        {
            for (k=0; k<icube->zres; k++)
                gwy_data_field_set_val(dfield, i, k, icube->data[i][jpos][k]);
        }

    }
    else if (jpos==-1 && kpos==-1)
    {
        gwy_data_field_resample(dfield, icube->yres, icube->zres, GWY_INTERPOLATION_NONE);
        for (j=0; j<icube->yres; j++)
        {
            for (k=0; k<icube->zres; k++)
                gwy_data_field_set_val(dfield, j, k, icube->data[ipos][j][k]);
        }
    }
}



/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

