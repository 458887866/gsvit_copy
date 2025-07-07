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


/*  yee.h : 
 *  computation domain data and core algorithms
 */


#ifndef YEE_H
#define YEE_H

#if defined (_WIN32) && defined (_MSC_VER)
#include "../config_win.h"
#else
#include "../config.h"
#endif

#include <libprocess/gwyprocess.h>
#include "dcube.h"
#include "fcube.h"
#include "icube.h"
#include "settings.h"
#include "boundary.h"
#include "source.h"
#include "settings.h"


typedef struct _SvYeeData {
   SvDCube *ex;        /*basic data arrays*/
   SvDCube *ey;
   SvDCube *ez;
   SvDCube *hx;
   SvDCube *hy;
   SvDCube *hz;

   SvDCube **dplrcx;   /*plrc calculation*/
   SvDCube **dplrcy;
   SvDCube **dplrcz;
   SvDCube **plrcx;
   SvDCube **plrcy;
   SvDCube **plrcz;
   SvDCube **iplrcx;
   SvDCube **iplrcy;
   SvDCube **iplrcz;

   SvDCube *exp;        /*plrc storage*/
   SvDCube *eyp;
   SvDCube *ezp;

   SvDCube **pxp;        /*ADE storage*/
   SvDCube **pyp;
   SvDCube **pzp;
   SvDCube **px;       
   SvDCube **py;
   SvDCube **pz;
   SvDCube **dpxp;
   SvDCube **dpyp;
   SvDCube **dpzp;
   SvDCube **dpx;       
   SvDCube **dpy;
   SvDCube **dpz;

   SvFCube *epsilon;   /*detailed material properties*/
   SvFCube *mu;
   SvFCube *sigma;
   SvFCube *sigast;
   
   SvICube *mat;       /*listed material properties*/
   SvICube *bnds;      /*local boundaries*/

   SvICube *sourceskip;  /*materials that should be skipped in TSF source, including boundaries*/

   gdouble dx;
   gdouble dy;
   gdouble dz;
   gdouble dt;

   gint xres;
   gint yres;
   gint zres;



} SvYeeData;


SvYeeData* sv_yee_data_new(SvSet *set, gdouble dt);
void sv_yee_data_allocate(SvYeeData *d, SvSet *set, gint xres, gint yres, gint zres, gdouble dx, gdouble dy, gdouble dz);
void sv_yee_data_free(SvYeeData *d, SvSet *set);

gint sv_yee_data_ystep_h(SvYeeData *d, SvSetBoundary *sb, SvBoundary *bnd, SvSet *set, SvMatProp *mats, gint nmat);
gint sv_yee_data_ystep_e(SvYeeData *d, SvSetBoundary *sb, SvBoundary *bnd, SvSource *src, SvSet *set, SvMatProp *mats, gint nmat);

gdouble sv_yee_data_get_epsilon(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k);
gdouble sv_yee_data_get_sigma(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k);
gint sv_yee_data_get_epsilon_sigma(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k, gdouble *epsval, gdouble *sigmaval);

gdouble sv_yee_data_get_cb(SvYeeData *mp, SvSet *set, SvMatProp *mats, gint nmat, gdouble dt, gint i, gint j, gint k);
void sv_yee_data_get_cabs(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gdouble dt, gint i, gint j, gint k, gdouble *cax, gdouble *cay, gdouble *caz, gdouble *cbx, gdouble *cby, gdouble *cbz);
gint sv_yee_data_get_ca_cb(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gdouble dt, gint i, gint j, gint k, gdouble *caval, gdouble *cbval);

gdouble sv_yee_data_get_mu(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k);
gdouble sv_yee_data_get_sigast(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k);
gint sv_yee_data_get_mu_sigast(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k, gdouble *muval, gdouble *sigastval);

gdouble sv_yee_data_get_db(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gdouble dt, gint i, gint j, gint k);
gint sv_yee_data_get_da_db(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gdouble dt, gint i, gint j, gint k, gdouble *daval, gdouble *dbval);

gboolean sv_yee_data_is_pec(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k);


SvMatType sv_yee_data_get_tsf_mattype(SvYeeData *d, SvSet *set, SvMatProp *mats, gint nmat, gint i, gint j, gint k);

gint sv_yee_data_get_tbc(SvYeeData *d, SvMatProp *mats, gint nmat, SvSet *set, gint i, gint j, gint k, gint *xtbc, gint *ytbc, gint *ztbc);


#endif /* YEE_H */

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

