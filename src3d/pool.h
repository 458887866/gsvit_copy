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


/*  pool.h : 
 *  most of the data and/or pointers to them 
 */


#ifndef POOL_H
#define POOL_H

#if defined (_WIN32) && defined (_MSC_VER)
#include "../config_win.h"
#else
#include "../config.h"
#endif
#include <libprocess/gwyprocess.h>
#include "settings.h"
#include "dcube.h"
#include "fcube.h"
#include "icube.h"
#include "source.h"
#include "yee.h"
#include "output.h"
#include "boundary.h"
#include "gpu.h"
#include "farfield.h"
#include "pfarfield.h"
#include "subgrid.h"

/*synchronisation between cpu and gpu*/
typedef enum {
	SV_SYNC_ALL = 0, /*all the fields including material parameters*/
	SV_SYNC_EH = 1,  /*only electric and magnetic fields*/
	SV_SYNC_E = 2,   /*only electric fields*/
	SV_SYNC_H = 3    /*only magnetic fields*/
} SvGpuSync;

/*main data structures*/
typedef struct _SvPool {

   SvYeeData *d;  /*all the computational domain data*/

   SvMatProp *mats;
   gint nmat;

   /*other data structures (note that settings are in settings.h*/
   SvSource *src;            /*sources data*/
   SvOutput *out;            /*outputs data*/
   SvBoundary *bnd;          /*external and internal boundary conditions*/
   SvFarfield *farfield;     /*far field calculation*/
   SvPFarfield *pfarfield;   /*periodic far field calculation*/
   SvSubgrid *sg;       /*data for local mesh refinement*/

#ifdef UCUDA
   SvGpu* gpu;
#endif

   SvSet *set;

} SvPool;



/*all the basic functions. Mostly called in main() in svit.c*/

SvPool * sv_pool_new(SvSet *set);
void sv_pool_allocate_gpus(SvPool *mp, SvSet *set);
int sv_pool_ystep_h(SvPool *mp, SvSet *set); 
int sv_pool_ystep_e(SvPool *mp, SvSet *set); 

void sv_pool_apply_source(SvPool *mp, SvSet *set);
void sv_pool_boundary_hstep(SvPool *mp, SvSet *set);
void sv_pool_boundary_estep(SvPool *mp, SvSet *set);
void sv_pool_boundary_copy(SvPool *mp, SvSet *set);
void sv_pool_getsum(SvPool *mp, SvSet *set);
void sv_pool_output_allocate(SvPool *mp, SvSet *set);
void sv_pool_output(SvPool *mp, SvSet *set);
void sv_pool_getforce(SvPool *mp, SvSet *set);

int sv_pool_gpu_init(SvPool *mp, SvSet *set);
int sv_pool_gpu_apply_source_hstep(SvPool *mp, SvSet *set);
int sv_pool_gpu_apply_source_estep(SvPool *mp, SvSet *set);

int sv_pool_gpu_ystep_h(SvPool *mp, SvSet *set); 
int sv_pool_gpu_ystep_e(SvPool *mp, SvSet *set); 
int sv_pool_gpu_copyto(SvPool *mp, SvSet *set, SvGpuSync what); 
int sv_pool_gpu_copyfrom(SvPool *mp, SvSet *set, SvGpuSync what);
int sv_pool_gpu_getsum(SvPool *mp, SvSet *set);
int sv_pool_gpu_boundary_copy(SvPool *mp, SvSet *set);
int sv_pool_gpu_boundary_hstep(SvPool *mp, SvSet *set);
int sv_pool_gpu_boundary_estep(SvPool *mp, SvSet *set);

void sv_pool_apply_source_hstep(SvPool *mp, SvSet *set);
void sv_pool_apply_source_estep(SvPool *mp, SvSet *set);
int sv_pool_gpu_apply_farfield(SvPool *mp, SvSet *set);
int sv_pool_gpu_getpoints(SvPool *mp, SvSet *set);
int sv_pool_gpu_getforce(SvPool *mp, SvSet *set);

void sv_pool_allocate(SvSet *set, SvPool *mp);
void sv_pool_allocate_output(SvPool *mp, SvSet *set);

int sv_pool_farfield(SvPool *mp, SvSet *set);
int sv_pool_farfield_allocate_storage(SvPool *mp, SvSet *set);

int sv_pool_pfarfield(SvPool *mp, SvSet *set);
int sv_pool_gpu_apply_pfarfield(SvPool *mp, SvSet *set);

int is_output(SvPool *mp, SvSet *set);

void sv_pool_getpoints(SvPool *mp, SvSet *set);


int sv_pool_subgrid_copy(SvPool *mp, SvSet *set);
int sv_pool_subgrid_step(SvPool *mp, SvSet *set);


SvBoundary* sv_boundary_new(SvPool *mp, SvSet *set);
SvFarfield* sv_farfield_new(SvPool *mp, SvSet *set);
SvPFarfield* sv_pfarfield_new(SvPool *mp, SvSet *set);
SvSubgrid *sv_subgrid_new(SvPool *mp, SvSet *set);


/*functions to get local material values*/
gdouble sv_pool_get_epsilon(SvPool *mp, SvSet *set, gint i, gint j, gint k);
gdouble sv_pool_get_sigma(SvPool *mp, SvSet *set, gint i, gint j, gint k);
gint    sv_pool_get_epsilon_sigma(SvPool *mp, SvSet *set, gint i, gint j, gint k, gdouble *epsilonval, gdouble *sigmaval);
gdouble sv_pool_get_cb(SvPool *mp, SvSet *set, gint i, gint j, gint k);
gint    sv_pool_get_ca_cb(SvPool *mp, SvSet *set, gint i, gint j, gint k, gdouble *caval, gdouble *cbval);

void
sv_pool_get_cabs(SvPool *mp, SvSet *set, gint i, gint j, gint k, gdouble *cax, gdouble *cay, gdouble *caz, gdouble *cbx, gdouble *cby, gdouble *cbz);

gdouble sv_pool_get_mu(SvPool *mp, SvSet *set, gint i, gint j, gint k);
gdouble sv_pool_get_sigast(SvPool *mp, SvSet *set, gint i, gint j, gint k);
gint    sv_pool_get_mu_sigast(SvPool *mp, SvSet *set, gint i, gint j, gint k, gdouble *muval, gdouble *sigastval);
gdouble sv_pool_get_db(SvPool *mp, SvSet *set, gint i, gint j, gint k);
gint    sv_pool_get_da_db(SvPool *mp, SvSet *set, gint i, gint j, gint k, gdouble *daval, gdouble *dbval);

SvMatType sv_pool_get_tsf_mattype(SvPool *mp, SvSet *set, gint i, gint j, gint k);

gint
get_tbc(SvPool *mp, SvSet *set, gint i, gint j, gint k, gint *xtbc, gint *ytbc, gint *ztbc);



void sv_pool_free(SvPool *mp, SvSet *set);

#ifdef UCUDA
SvGpu* sv_gpu_new(SvPool *mp, SvSet *set);
#endif


#endif /* POOL_H */

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

