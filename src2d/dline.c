
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


/*  dline.c : 
 *  simple 1D data representation
 *  FIXME: Should be replaced by GwyDataLine
 */

#include <stdio.h>
#include <math.h>
#include "dline.h"


SvDLine *
sv_dline_new(gint res, 
                      gdouble real,
                      gboolean nullme)
{
    SvDLine *dl = (SvDLine *)g_malloc(sizeof(SvDLine));

    dl->res = res;
    dl->real = real;
    dl->data= (gdouble *) g_malloc(res*sizeof(gdouble));
    if (nullme) sv_dline_fill(dl, 0);
    return dl; 
}

SvDLine*
sv_dline_new_alike (SvDLine *dline, gboolean nullme)
{
    return sv_dline_new(dline->res,
                        dline->real, nullme);
}

void 
sv_dline_free(SvDLine *dline)
{
    g_free((void *)dline->data);
    dline->data=NULL;
}

void 
sv_dline_fill(SvDLine *dline, gdouble value)
{
    gint i;
    for (i=0; i<dline->res; i++)
        dline->data[i] = value;
}

void 
sv_dline_add(SvDLine *dline, gdouble value)
{
    gint i;
    for (i=0; i<dline->res; i++)
        dline->data[i] += value;
}

void 
sv_dline_multiply(SvDLine *dline, gdouble value)
{
    gint i;
    for (i=0; i<dline->res; i++)
        dline->data[i] *= value;
}

void 
sv_dline_get_dataline   (SvDLine *dline, GwyDataLine *gwydline)
{
    gint i;
    gwy_data_line_resample(gwydline, dline->res, GWY_INTERPOLATION_NONE);
    for (i=0; i<dline->res; i++)
        gwy_data_line_set_val(gwydline, i, dline->data[i]);

}

gdouble 
sv_dline_get_dval(SvDLine *dline, gdouble x)
{
    gdouble w1, w2, w3, w4;
    gint l = (gint)floor(x);
    gdouble a = x-(gdouble)l;

    if (x>=1 && x<(dline->res-1))
    {
        w1=a+1; w2=a; w3=1-a; w4=2-a;
        w1=4-8*w1+5*w1*w1-w1*w1*w1;
        w2=1-2*w2*w2+w2*w2*w2;
        w3=1-2*w3*w3+w3*w3*w3;
        w4=4-8*w4+5*w4*w4-w4*w4*w4;
        return w1*dline->data[l-1]+w2*dline->data[l]+w3*dline->data[l+1]+w4*dline->data[l+2];
    } else if ((x<1 && x>=0) || (x>=(dline->res-1) && x<(dline->res)))
    {
        return (1-a)*dline->data[l]+a*dline->data[l+1];
    }
    else return 0;
}


/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

