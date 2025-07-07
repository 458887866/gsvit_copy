
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


/*  source.h : 
 *  data for computation of elmag field sources
 */

#ifndef SOURCE_H
#define SOURCE_H

#include <glib.h>
#include "1dpool.h"
#include "zpool.h"

typedef struct
{
    gint ndata;
    gint *layered_zpos;
    gdouble* ex;
    gdouble* ey;
    gdouble* ez;
    gdouble* hx;
    gdouble* hy;
    gdouble* hz;
} SvSrcData;

typedef struct
{
    gint i;
    gint j;
    gint k;
    SvSrcData sdata;
} SvSourcePoint;

typedef struct
{
    gint npos;
    gint *xpos;
    gint *ypos;
    gint *zpos;
    gdouble *xdir;
    gdouble *ydir;
    gdouble *zdir;
    gdouble *source_amplitude;
    gdouble *phase; 
    gdouble *lambda;
    gdouble avgfield; 
} SvSourceLocal;

typedef struct
{
    Sv1DPool *jp;
    gint box_i0;
    gint box_j0;
    gint box_k0;
    gint box_in;
    gint box_jn;
    gint box_kn;
    gint box_boundary_skipi0;
    gint box_boundary_skipin;
    gint box_boundary_skipj0;
    gint box_boundary_skipjn;
    gint box_boundary_skipk0;
    gint box_boundary_skipkn;
    gdouble ia_theta;   //angle
    gdouble ia_phi;     //angle
    gdouble ia_psi;     //polarisation
    gdouble *e;
    gint *layered_zpos;
    gdouble layered_epsilon;
    gdouble layered_sigma;
    gdouble layered_mu;
    gdouble layered_sigast;
    gdouble corr;
    gint ndata;

} SvSourceTSF;

typedef struct
{
    SvZPool *jp;
    gint box_i0;
    gint box_j0;
    gint box_k0;
    gint box_in;
    gint box_jn;
    gint box_kn;
    gint box_boundary_skipi0;
    gint box_boundary_skipin;
    gint box_boundary_skipj0;
    gint box_boundary_skipjn;
    gint box_boundary_skipk0;
    gint box_boundary_skipkn;
    gdouble ia_theta;   //angle
    gdouble ia_phi;     //angle
    gdouble ia_psi;     //polarisation
    gdouble *e;
    gint *layered_zpos;
    gint lpos[50];
    gdouble layered_epsilon[50];
    gdouble layered_sigma[50];
    gdouble layered_mu[50];
    gdouble layered_sigast[50];
    gdouble corr;
    gint ndata;
    gint layered_count;
    gdouble timeshift;

} SvSourceLTSF;

typedef struct
{
    Sv1DPool ***jp;
    gint box_i0;
    gint box_j0;
    gint box_k0;
    gint box_in;
    gint box_jn;
    gint box_kn;
    gint box_boundary_skipi0;
    gint box_boundary_skipin;
    gint box_boundary_skipj0;
    gint box_boundary_skipjn;
    gint box_boundary_skipk0;
    gint box_boundary_skipkn;
    gdouble focused_thetamax; //aperture max angle
    gdouble focused_fdist;     //focal length
    gint focused_nip;     //N integration points
    gint focused_mip;     //M integration points
    gdouble *e;    //incident field (to be recalculated)
    gint *layered_zpos;
    gdouble layered_epsilon;
    gdouble layered_sigma;
    gdouble layered_mu;
    gdouble layered_sigast;
    gint ndata;
    gdouble **ia_theta;
    gdouble **ia_phi;
    gdouble **ia_psi;
    gdouble **corr;
    gdouble **an;
    gdouble **bm;
} SvSourceTSFF;

typedef struct
{
    SvZPool ***jp;
    gint box_i0;
    gint box_j0;
    gint box_k0;
    gint box_in;
    gint box_jn;
    gint box_kn;
    gint box_boundary_skipi0;
    gint box_boundary_skipin;
    gint box_boundary_skipj0;
    gint box_boundary_skipjn;
    gint box_boundary_skipk0;
    gint box_boundary_skipkn;
    gdouble focused_thetamax; //aperture max angle
    gdouble focused_fdist;     //focal length
    gint focused_nip;     //N integration points
    gint focused_mip;     //M integration points
    gdouble *e;
    gint *layered_zpos;
    gint lpos[50];
    gdouble layered_epsilon[50];
    gdouble layered_sigma[50];
    gdouble layered_mu[50];
    gdouble layered_sigast[50];
    gint layered_count;
    gint ndata;
    gdouble **ia_theta;
    gdouble **ia_phi;
    gdouble **ia_psi;
    gdouble **corr;
    gdouble **an;
    gdouble **bm;
    gdouble **timeshift;
} SvSourceLTSFF;

typedef struct
{
    Sv1DPool *jp;
    gdouble ia_theta;  //angle
    gdouble ia_phi;    //angle
    gdouble ia_psi;    //polarisation
    gdouble *e;
    gint *layered_zpos;
    gdouble layered_epsilon;
    gdouble layered_sigma;
    gdouble layered_mu;
    gdouble layered_sigast;
    int ndata;

} SvSourceSF;

typedef struct
{
    gdouble **phi;
    gdouble **exs;
    gdouble **eys;
    gdouble **hxs;
    gdouble **hys;
} SvSourceSheet;

typedef struct
{
    SvSourcePoint *sp;
    SvSourceLocal *sl;
    SvSourceTSF *tsf;
    SvSourceLTSF *ltsf;
    SvSourceTSFF *tsff;
    SvSourceLTSFF *ltsff;
    SvSourceSF *sf;
    SvSourceSheet *sh;
} SvSource;

SvSource* sv_source_new();

gdouble dcomp(gint i, gint j, gint k, gint xres, gint yres, gint zres, gdouble theta, gdouble phi,
            gint i0, gint i1, gint j0, gint j1, gint k0, gint k1);
gdouble gex(gdouble field, gdouble theta, gdouble phi, gdouble psi);
gdouble gey(gdouble field, gdouble theta, gdouble phi, gdouble psi);
gdouble gez(gdouble field, gdouble theta, gdouble phi, gdouble psi);
gdouble ghx(gdouble field, gdouble theta, gdouble phi, gdouble psi);
gdouble ghy(gdouble field, gdouble theta, gdouble phi, gdouble psi);
gdouble ghz(gdouble field, gdouble theta, gdouble phi, gdouble psi);


#endif /* SOURCE_H */
/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
