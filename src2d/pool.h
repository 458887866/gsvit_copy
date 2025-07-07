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


#ifndef SV_POOL
#define SV_POOL

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
#include "farfield.h"
#include "output.h"
#include "boundary.h"
#include "gpu.h"

typedef enum {
	SV_SYNC_ALL = 0, /*all the fields including material parameters*/
	SV_SYNC_EH = 1,  /*only electric and magnetic fields*/
	SV_SYNC_E = 2,   /*only electric fields*/
	SV_SYNC_H = 3    /*only magnetic fields*/
} SvGpuSync;

/*material properties*/
typedef struct {
   gint type;
   gdouble epsilon;
   gdouble mu;
   gdouble sigma;
   gdouble sigast;
   gdouble drude_omega_p;
   gdouble drude_nu;
   gdouble cp3_a[3];
   gdouble cp3_phi[3];
   gdouble cp3_omega[3];
   gdouble cp3_gamma[3];
   gint pos;
} SvMatProp;

/*main data structure*/
typedef struct _SvPool {
   GwyDataField *ex;        /*basic data arrays*/
   GwyDataField *ey;
   GwyDataField *ez;
   GwyDataField *hx;
   GwyDataField *hy;
   GwyDataField *hz;

   //SvDCube **plrcx;
   //SvDCube **plrcy;
   //SvDCube **plrcz;
   //SvDCube **iplrcx;
   //SvDCube **iplrcy;
   //SvDCube **iplrcz;

   GwyDataField *epsilon;   /*detailed material properties*/
   GwyDataField *mu;
   GwyDataField *sigma;
   GwyDataField *sigast;
   
   GwyDataField *mat;       /*listed material properties*/
   SvMatProp *mats;
   gint nmat;

   /*other data structures*/
   SvSource *src;
   SvOutput *out;
   SvBoundary *bnd;
   SvFarfield *farfield;
#ifdef UCUDA
   SvGpu* gpu;
#endif

   SvSet *set;

} SvPool;

SvPool *sv_pool_new(SvSet *set);
void sv_pool_allocate_gpus(SvPool *mp, SvSet *set);
int sv_pool_ystep_h(SvPool *mp, SvSet *set); 
int sv_pool_ystep_e(SvPool *mp, SvSet *set); 

void sv_pool_add_material(SvPool *mp, SvMatProp *mat);

void sv_pool_apply_source_estep(SvPool *mp, SvSet *set);
void sv_pool_apply_source_hstep(SvPool *mp, SvSet *set);

void sv_pool_boundary_hstep(SvPool *mp, SvSet *set);
void sv_pool_boundary_estep(SvPool *mp, SvSet *set);
void sv_pool_boundary_copy(SvPool *mp, SvSet *set);
void sv_pool_getsum(SvPool *mp, SvSet *set);
void sv_pool_output_allocate(SvPool *mp, SvSet *set);
void sv_pool_output(SvPool *mp, SvSet *set);

void sv_pool_farfield(SvPool *mp, SvSet *set);

void sv_pool_gpu_apply_source_estep(SvPool *mp, SvSet *set);
int sv_pool_gpu_apply_source_hstep(SvPool *mp, SvSet *set);

int sv_pool_gpu_init(SvPool *mp, SvSet *set);

int sv_pool_gpu_ystep_h(SvPool *mp, SvSet *set); 
int sv_pool_gpu_ystep_e(SvPool *mp, SvSet *set); 
int sv_pool_gpu_copyto(SvPool *mp, SvSet *set, SvGpuSync what); 
int sv_pool_gpu_copyfrom(SvPool *mp, SvSet *set, SvGpuSync what);
int sv_pool_gpu_getsum(SvPool *mp, SvSet *set);
int sv_pool_gpu_boundary_copy(SvPool *mp, SvSet *set);
int sv_pool_gpu_boundary_hstep(SvPool *mp, SvSet *set);
int sv_pool_gpu_boundary_estep(SvPool *mp, SvSet *set);

int is_output(SvPool *mp, SvSet *set);


#endif
/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

