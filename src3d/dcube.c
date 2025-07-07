 
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

/*  dcube.c : 
 *  Simple 3D data representation.
 *  TODO: should be replaced by GwyBrick in future
 */

#include <stdio.h>
#include <omp.h>
#include "dcube.h"

SvDCube*
sv_dcube_new(gint xres, gint yres, gint zres, gdouble xreal, gdouble yreal, gdouble zreal, gboolean nullme)
{
    gint i, j;
    SvDCube *dc = (SvDCube *)g_malloc(sizeof(SvDCube));

    dc->xres = xres;
    dc->yres = yres; 
    dc->zres = zres;
    dc->xreal = xreal;
    dc->yreal = yreal;
    dc->zreal = zreal;
    dc->data = (gdouble ***) g_malloc(xres*sizeof(gdouble**));

    for (i = 0; i < xres; i++) {
        dc->data[i] = (gdouble **) g_malloc(yres*sizeof(gdouble*));
        for (j = 0; j < yres; j++)
            dc->data[i][j] = (gdouble *) g_malloc(zres*sizeof(gdouble));
    }
    if (nullme) 
        sv_dcube_fill(dc, 0);

    return dc; 
}

SvDCube*
sv_dcube_new_alike (SvDCube *dcube, gboolean nullme)
{
    return sv_dcube_new(dcube->xres, dcube->yres, dcube->zres, dcube->xreal, dcube->yreal, dcube->zreal, nullme);
}

void sv_dcube_free(SvDCube *dcube)
{
    gint i, j;
    for (i = 0; i < dcube->xres; i++) {
        for (j = 0; j < dcube->yres; j++)
            g_free((void **) dcube->data[i][j]);
        g_free((void *) dcube->data[i]);
    }

    g_free((void *)dcube->data);
    dcube->data=NULL;
}

//! paralelizovt, az tohle dela alokaci dat!!!
void
sv_dcube_fill(SvDCube *dcube, gdouble value)
{
    gint i, j, k;
    
#ifndef G_OS_WIN32
#pragma omp parallel default(shared) private(j, k)
#pragma omp for nowait
#endif
    
    for (i = 0; i < dcube->xres; i++) {
        for (j = 0; j < dcube->yres; j++) {
            for (k = 0; k < dcube->zres; k++)
                dcube->data[i][j][k] = value;
        }
    }
}

void
sv_dcube_add(SvDCube *dcube, gdouble value)
{
    gint i, j, k;
    
#ifndef G_OS_WIN32
#pragma omp parallel default(shared) private(j, k)
#pragma omp for nowait
#endif
    
    for (i = 0; i < dcube->xres; i++) {
        for (j = 0; j < dcube->yres; j++) {
            for (k = 0; k < dcube->zres; k++)
                dcube->data[i][j][k] += value;
        }
    }
}

gdouble
sv_dcube_get_sum(SvDCube *dcube)
{
    gint i, j, k;
    gdouble val = 0;
    
#ifndef G_OS_WIN32
#pragma omp parallel default(shared) private(j, k)
#pragma omp for nowait
#endif
 
    for (i = 0; i < dcube->xres; i++) {
        for (j = 0; j < dcube->yres; j++) {
            for (k = 0; k < dcube->zres; k++)
                val += dcube->data[i][j][k];
        }
    }
    return val;
}

void
sv_dcube_multiply(SvDCube *dcube, gdouble value)
{
    gint i = 0, j = 0, k = 0;
    
#ifndef G_OS_WIN32
#pragma omp parallel default(shared) private(j, k)
#pragma omp for nowait
#endif

    for (i = 0; i < dcube->xres; i++) {
        for (j = 0; j < dcube->yres; j++) {
            for (k = 0; k < dcube->zres; k++)
                dcube->data[i][j][k] *= value;
        }
    }
}

void
sv_dcube_get_datafield(SvDCube *dcube, GwyDataField *dfield, gint ipos, gint jpos, gint kpos)
{
    gint i, j, k;
    if (ipos == -1 && jpos == -1) {
        gwy_data_field_resample(dfield, dcube->xres, dcube->yres, GWY_INTERPOLATION_NONE);
        for (i = 0; i < dcube->xres; i++) {
            for (j = 0; j < dcube->yres; j++)
                gwy_data_field_set_val(dfield, i, j, dcube->data[i][j][kpos]);
        }
    }
    else if (ipos == -1 && kpos == -1) {
        gwy_data_field_resample(dfield, dcube->xres, dcube->zres, GWY_INTERPOLATION_NONE);
        for (i = 0; i < dcube->xres; i++) {
            for (k = 0; k < dcube->zres; k++)
                gwy_data_field_set_val(dfield, i, k, dcube->data[i][jpos][k]);
        }
    }
    else if (jpos == -1 && kpos == -1) {
        gwy_data_field_resample(dfield, dcube->yres, dcube->zres, GWY_INTERPOLATION_NONE);
        for (j = 0; j < dcube->yres; j++) {
            for (k = 0; k < dcube->zres; k++)  
                gwy_data_field_set_val(dfield, j, k, dcube->data[ipos][j][k]);
        }
    }
}

void
sv_dcube_get_datafield_area(SvDCube *dcube, GwyDataField *dfield, gint ipos, gint jpos, gint kpos, gint ioffset, gint joffset, gint koffset)
{
    gint i, j, k;
    if (ipos == -1 && jpos == -1) {
        for (i = 0; i < gwy_data_field_get_xres(dfield); i++) {
            for (j = 0; j < gwy_data_field_get_yres(dfield); j++){
//                printf("get_dfield_area ij:  dcube %d %d %d   to field %d %d\n", i+ioffset, j+joffset, kpos, i, j);
                gwy_data_field_set_val(dfield, i, j, dcube->data[i+ioffset][j+joffset][kpos]);
            }
        }
    }
    else if (ipos == -1 && kpos == -1) {
        for (i = 0; i < gwy_data_field_get_xres(dfield); i++) {
            for (k = 0; k < gwy_data_field_get_yres(dfield); k++) {
//                printf("get_dfield_area ik:  dcube %d %d %d   to field %d %d\n", i+ioffset, jpos, k+koffset, i, k);
                 gwy_data_field_set_val(dfield, i, k, dcube->data[i+ioffset][jpos][k+koffset]);
            }
        }
    }
    else if (jpos == -1 && kpos == -1) {
        for (j = 0; j < gwy_data_field_get_xres(dfield); j++) {
            for (k = 0; k < gwy_data_field_get_yres(dfield); k++) {
//                printf("get_dfield_area jk:  dcube %d %d %d   to field %d %d\n", ipos, j+joffset, k+koffset, j, k);
                 gwy_data_field_set_val(dfield, j, k, dcube->data[ipos][j+joffset][k+koffset]);
            }
        }
    }
}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
