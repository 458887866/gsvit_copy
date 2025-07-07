
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

/*  gpu.c : 
 *  all the GPU related functions that are performed on both CPU and GPU (preparation of data, etc).
 */

#if defined (_WIN32) && defined (_MSC_VER)
#include "../config_win.h"
#else
#include "../config.h"
#endif
#include "pool.h"
#include "settings.h"
#include "gpu.h"

#ifdef UCUDA

#include <cuda_runtime.h>
#define AT(i, j, k) ((i)*plan->zres*plan->yres + (j)*plan->zres + (k))

#define GPU_PERIODIC_HLPS_LIMIT 120


void
sv_pool_sync(SvPool *mp, SvSet *set);

int
sv_pool_gpu_init(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose>1) printf("pool_gpu_init\n");
    for (n=0; n<gpu->ngset; n++)
       g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_INIT));
    sv_pool_sync(mp, set);
    return 0;
}

int
sv_pool_gpu_ystep_h(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running Yee H GPU step...           ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++)
       g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_YSTEPH));
    sv_pool_sync(mp, set);
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_ystep_e(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running Yee E GPU step...           ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++)
       g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_YSTEPE));
    sv_pool_sync(mp, set);
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;

}

int
sv_pool_gpu_liao_cpy(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running copy boundary GPU step...   ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++) {
       g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_LIAO_CPY));
    }
    sv_pool_sync(mp, set);
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_cpml_hstep(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running CPML boundary H GPU step...   ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++) {
       g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_CPML_HSTEP));
    }
    sv_pool_sync(mp, set);
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_cpml_estep(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running CPML boundary E GPU step...   ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++) {
       g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_CPML_ESTEP));
    }
    sv_pool_sync(mp, set);
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_liao_bnd(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running Liao boundary GPU step...   ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++)
       g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_LIAO_BND));
    sv_pool_sync(mp, set);
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_getsum(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpuPlan *plan;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (!set->so.nsums) return 0;

    if (set->sc.verbose > 1) {
	    printf("Running summing GPU step...   ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++) {
        plan = gpu->gset[n].aplan->plan;
        plan->step = set->sc.step_act;
    }

    for (n=0; n<gpu->ngset; n++) {
       g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_SUM));
    }
    sv_pool_sync(mp, set);

    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_getforce(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpuPlan *plan;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (!set->so.nforces) return 0;

    if (set->sc.verbose > 1) {
	    printf("Running summing force GPU step...   ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++) {
        plan = gpu->gset[n].aplan->plan;
        plan->step = set->sc.step_act;
    }

    for (n=0; n<gpu->ngset; n++) {
       g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_FORCE));
    }
    sv_pool_sync(mp, set);

    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_getpoints(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpuPlan *plan;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (!set->so.npnts) return 0;

    if (set->sc.verbose > 1) {
	    printf("Running output points GPU step...   ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++) {
        plan = gpu->gset[n].aplan->plan;
        plan->step = set->sc.step_act;
    }

    for (n=0; n<gpu->ngset; n++) {
       g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_OPNT));
    }
    sv_pool_sync(mp, set);

    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_mbnd_e(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running internal boundary GPU step(s)...   ");
	    fflush(stdout);
    }

    if (set->smb.bx0 == 4 || set->smb.bxn == 4)
    {
        for (n=0; n<gpu->ngset; n++)
            g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_MBNDX_E));
        sv_pool_sync(mp, set);
    }
    if (set->smb.by0 == 4 || set->smb.byn == 4)
    {
        for (n=0; n<gpu->ngset; n++)
            g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_MBNDY_E));
        sv_pool_sync(mp, set);
    }
    if (set->smb.bz0 == 4 || set->smb.bzn == 4)
    {
        for (n=0; n<gpu->ngset; n++)
            g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_MBNDZ_E));
        sv_pool_sync(mp, set);
    }
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_mbnd_h(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running internal boundary GPU step(s)...   ");
	    fflush(stdout);
    }

    if (set->smb.bx0 == 4 || set->smb.bxn == 4)
    {
        for (n=0; n<gpu->ngset; n++)
            g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_MBNDX_H));
        sv_pool_sync(mp, set);
    }
    if (set->smb.by0 == 4 || set->smb.byn == 4)
    {
        for (n=0; n<gpu->ngset; n++)
            g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_MBNDY_H));
        sv_pool_sync(mp, set);
    }
    if (set->smb.bz0 == 4 || set->smb.bzn == 4)
    {
        for (n=0; n<gpu->ngset; n++)
            g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_MBNDZ_H));
        sv_pool_sync(mp, set);
    }
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_point_srce(SvPool *mp, SvSet *set)
{
    gint i, j, n, pos;
    SvGpuPlan *plan;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running source GPU step...          ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++)
    {
	    plan = gpu->gset[n].aplan->plan;

	    for (i=0; i<set->ss.npnts; i++)
	    {
		    /*if we are lucky, data are not shifted by user and can be used directly*/
		    pos = -1;
		    if (set->sc.step_act <= mp->src->sp[i].sdata.ndata &&
				    mp->src->sp[i].sdata.layered_zpos[set->sc.step_act] == set->sc.step_act)
			    pos = set->sc.step_act;
		    else {  /*otherwise we need to search for right position*/
			    for (j=0; j<mp->src->sp[i].sdata.ndata; j++) {
				    if (mp->src->sp[i].sdata.layered_zpos[j] == set->sc.step_act)
				    {
					    pos = j;
					    break;
				    }
			    }
		    }
		    if (pos>=0)
		    {
			    if (mp->src->sp[i].i < plan->xfrom || mp->src->sp[i].i >= plan->xto) continue;

			    plan->source_i = mp->src->sp[i].i - plan->xfrom;
			    plan->source_j = mp->src->sp[i].j;
			    plan->source_k = mp->src->sp[i].k;
			    plan->source_ex = (float)mp->src->sp[i].sdata.ex[pos];
			    plan->source_ey = (float)mp->src->sp[i].sdata.ey[pos];
			    plan->source_ez = (float)mp->src->sp[i].sdata.ez[pos];
	
			    g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_SRCE_POINT));
			    sv_pool_sync(mp, set);
		    }
	    }
    }
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_tsf_estep(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running TSF source GPU estep...          ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++)
    {
        g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_TSF_ESTEP));
        sv_pool_sync(mp, set);
    }
    
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_tsff_estep(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running TSFF source GPU estep...          ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++)
    {
        g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_TSFF_ESTEP));
        sv_pool_sync(mp, set);
    }
    
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_apply_farfield(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;
    SvGpuPlan *plan;

    if (set->sf.nrs == 0) return 0;

    if (set->sc.verbose > 1) {
	    printf("Running farfield GPU estep...       ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++)
    {
        plan = gpu->gset[n].aplan->plan;
        if (plan[n].nhlps==0) g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_FF));
        else g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_FFA));
        sv_pool_sync(mp, set);
    }
    
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_apply_pfarfield(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->spf.nrs == 0) return 0;
    if (set->sc.verbose > 1) {
	    printf("Running periodic farfield GPU estep...       ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++)
    {
        g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_PFF));
        sv_pool_sync(mp, set);
    }
    
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_tsf_hstep(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpuPlan *plan;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running TSF 1dpool step...          ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++)
    {
	plan = gpu->gset[n].aplan->plan;
        plan->source_ex = (float)mp->src->tsf->e[set->sc.step_act];
        g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_TSF_JSTEP));
        sv_pool_sync(mp, set);
    }
    if (set->sc.verbose > 1) printf("done.\n");


    if (set->sc.verbose > 1) {
	    printf("Running TSF source GPU hstep...          ");
	    fflush(stdout);
    }

    
    for (n=0; n<gpu->ngset; n++)
    {
        g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_TSF_HSTEP));
        sv_pool_sync(mp, set);
    }
    
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_tsff_hstep(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpuPlan *plan;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running TSFF 1dpool step...          ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++)
    {
	plan = gpu->gset[n].aplan->plan;
        plan->source_ex = (float)mp->src->tsff->e[set->sc.step_act];
        g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_TSFF_JSTEP));
        sv_pool_sync(mp, set);
    }
    if (set->sc.verbose > 1) printf("done.\n");


    if (set->sc.verbose > 1) {
	    printf("Running TSFF source GPU hstep...          ");
	    fflush(stdout);
    }

    
    for (n=0; n<gpu->ngset; n++)
    {
        g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_TSFF_HSTEP));
        sv_pool_sync(mp, set);
    }
    
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_sf_hstep(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpuPlan *plan;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running SF 1dpool step...          ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++)
    {
	plan = gpu->gset[n].aplan->plan;
        plan->source_ex = (float)mp->src->sf->e[MAX(0, set->sc.step_act-1)];
        g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_SF_JSTEP));
        sv_pool_sync(mp, set);
    }
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
}

int
sv_pool_gpu_boundary_copy(SvPool *mp, SvSet *set)
{
//    if (set->sb.bx0 == SV_BOUNDARY_LIAO || set->sb.bxn == SV_BOUNDARY_LIAO ||
//        set->sb.by0 == SV_BOUNDARY_LIAO || set->sb.byn == SV_BOUNDARY_LIAO ||
//        set->sb.bz0 == SV_BOUNDARY_LIAO || set->sb.bzn == SV_BOUNDARY_LIAO)
   
    //run this always to increase step, if there is no liao boundary it does nothing else.
    sv_pool_gpu_liao_cpy(mp, set);
    return 0;
}

int
sv_pool_gpu_boundary_hstep(SvPool *mp, SvSet *set)
{    
    if (set->smb.bx0 == 4 || set->smb.bxn == 4 
        || set->smb.by0 == 4 || set->smb.byn == 4 
        || set->smb.bz0 == 4 || set->smb.bzn == 4) 
        sv_pool_gpu_mbnd_h(mp, set);    

    if (set->sb.bx0 == SV_BOUNDARY_CPML || set->sb.bxn == SV_BOUNDARY_CPML ||
        set->sb.by0 == SV_BOUNDARY_CPML || set->sb.byn == SV_BOUNDARY_CPML ||
        set->sb.bz0 == SV_BOUNDARY_CPML || set->sb.bzn == SV_BOUNDARY_CPML)
    sv_pool_gpu_cpml_hstep(mp, set);

    return 0;
}

int
sv_pool_gpu_boundary_estep(SvPool *mp, SvSet *set)
{    
    if (set->sb.bx0 == SV_BOUNDARY_LIAO || set->sb.bxn == SV_BOUNDARY_LIAO ||
        set->sb.by0 == SV_BOUNDARY_LIAO || set->sb.byn == SV_BOUNDARY_LIAO ||
        set->sb.bz0 == SV_BOUNDARY_LIAO || set->sb.bzn == SV_BOUNDARY_LIAO)
    sv_pool_gpu_liao_bnd(mp, set);

    if (set->sb.bx0 == SV_BOUNDARY_CPML || set->sb.bxn == SV_BOUNDARY_CPML ||
        set->sb.by0 == SV_BOUNDARY_CPML || set->sb.byn == SV_BOUNDARY_CPML ||
        set->sb.bz0 == SV_BOUNDARY_CPML || set->sb.bzn == SV_BOUNDARY_CPML)
    sv_pool_gpu_cpml_estep(mp, set);


    if (set->smb.bx0 == 4 || set->smb.bxn == 4 
        || set->smb.by0 == 4 || set->smb.byn == 4 
        || set->smb.bz0 == 4 || set->smb.bzn == 4) 
        sv_pool_gpu_mbnd_e(mp, set);
    return 0;
}

int
sv_pool_gpu_apply_source_estep(SvPool *mp, SvSet *set)
{
    sv_pool_gpu_point_srce(mp, set);
    if (mp->src->tsf) sv_pool_gpu_tsf_estep(mp, set);
    if (mp->src->tsff) sv_pool_gpu_tsff_estep(mp, set);
    return 0;
}

int
sv_pool_gpu_apply_source_hstep(SvPool *mp, SvSet *set)
{
    if (mp->src->tsf) sv_pool_gpu_tsf_hstep(mp, set);
    if (mp->src->tsff) sv_pool_gpu_tsff_hstep(mp, set);
    if (mp->src->sf) sv_pool_gpu_sf_hstep(mp, set);
    return 0;
}

int
sv_pool_gpu_copyto(SvPool *mp, SvSet *set, SvGpuSync what)
{
    gint n, i, j, k;
    SvGpu* gpu = (SvGpu *)mp->gpu;
    SvGpuPlan *plan;
    long int pos, ip;

    if (set->sc.verbose>1) {
        printf("Copying data to GPU...");
        fflush(stdout);
    }
    for (n=0; n<gpu->ngset; n++)
    {
	    plan = gpu->gset[n].aplan->plan;
	    for (i=0; i<plan->xres; i++)
		    for (j=0; j<plan->yres; j++)
			    for (k=0; k<plan->zres; k++)
			    {
				    pos = AT(i, j, k);
				    ip = i+plan->xfrom;
				    if (what == SV_SYNC_ALL) {
					    if (set->plan.matmode==SV_MATMODE_FULL || set->plan.matmode==SV_MATMODE_ELECTRIC) 
					    {
						    plan->h_epsilon[pos] = mp->d->epsilon->data[ip][j][k];
						    plan->h_sigma[pos] = mp->d->sigma->data[ip][j][k];
					    }
					    if (set->plan.matmode==SV_MATMODE_FULL || set->plan.matmode==SV_MATMODE_MAGNETIC) 
					    {
						    plan->h_mu[pos] = mp->d->mu->data[ip][j][k];
						    plan->h_sigast[pos] = mp->d->sigast->data[ip][j][k];
					    }
				    }
				    if (what == SV_SYNC_ALL || what == SV_SYNC_EH || what == SV_SYNC_E) 
				    {
					    plan->h_ex[pos] = (float)mp->d->ex->data[ip][j][k];
					    plan->h_ey[pos] = (float)mp->d->ey->data[ip][j][k];
					    plan->h_ez[pos] = (float)mp->d->ez->data[ip][j][k];
				    }
				    if (what == SV_SYNC_ALL || what == SV_SYNC_EH || what == SV_SYNC_H)
				    {
					    plan->h_hx[pos] = (float)mp->d->hx->data[ip][j][k];
					    plan->h_hy[pos] = (float)mp->d->hy->data[ip][j][k];
					    plan->h_hz[pos] = (float)mp->d->hz->data[ip][j][k];
				    }
				    if (mp->nmat)
				    {
					    plan->h_mat[pos] = mp->d->mat->data[ip][j][k];
				    }
			    }
	    g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_COPYTO));
    }
    sv_pool_sync(mp, set);
    if (set->sc.verbose > 1) printf("done.\n");

    return 0;
}

int
sv_pool_gpu_copyfrom(SvPool *mp, SvSet *set, SvGpuSync what)
{
    gint n, i, j, k;
    SvGpu* gpu = (SvGpu *)mp->gpu;
    SvGpuPlan *plan;
    long int pos, ip;

    if (set->sc.verbose>1) {
        printf("Copying data from GPU...");
        fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++)
       g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_COPYFROM));

    sv_pool_sync(mp, set);

    for (n=0; n<gpu->ngset; n++)
    {
	    plan = gpu->gset[n].aplan->plan;
	    for (i=0; i<plan->xres; i++)
		    for (j=0; j<plan->yres; j++)
			    for (k=0; k<plan->zres; k++)
			    {
				    pos = AT(i, j, k);
				    ip = i+plan->xfrom;
				    if (what == SV_SYNC_ALL || what == SV_SYNC_EH || what == SV_SYNC_E)
				    {
					    mp->d->ex->data[ip][j][k] = plan->h_ex[pos];
					    mp->d->ey->data[ip][j][k] = plan->h_ey[pos];
					    mp->d->ez->data[ip][j][k] = plan->h_ez[pos];
				    }
				    if (what == SV_SYNC_ALL || what == SV_SYNC_EH || what == SV_SYNC_H)
				    {
					    mp->d->hx->data[ip][j][k] = plan->h_hx[pos];
					    mp->d->hy->data[ip][j][k] = plan->h_hy[pos];
					    mp->d->hz->data[ip][j][k] = plan->h_hz[pos];
				    }
			    }

	    if (set->sf.nrs) {
		    j = 0;
		    for (k=0; k<set->sf.nrs; k++) {
			    for (i=0; i<mp->farfield->ndata; i++)
			    {
				   mp->farfield->rpoints[k].ex->data[i] = plan->h_ff_ex[j];
                                   mp->farfield->rpoints[k].ey->data[i] = plan->h_ff_ey[j];
                                   mp->farfield->rpoints[k].ez->data[i] = plan->h_ff_ez[j];
                                   j++;
			    }
		    }
	    }
	    if (set->spf.nrs) {
                    printf("copy back %d x %d = %d data\n", set->spf.nrs, mp->pfarfield->ndata, set->spf.nrs*mp->pfarfield->ndata);
		    j = 0;
		    for (k=0; k<set->spf.nrs; k++) {
			    for (i=0; i<mp->pfarfield->ndata; i++)
			    {
				   mp->pfarfield->prpoints[k].ex->data[i] = plan->h_pff_ex[j];
                                   mp->pfarfield->prpoints[k].ey->data[i] = plan->h_pff_ey[j];
                                   mp->pfarfield->prpoints[k].ez->data[i] = plan->h_pff_ez[j];
                           //        printf("%g %g %g\n", mp->pfarfield->prpoints[k].ex->data[i],  mp->pfarfield->prpoints[k].ey->data[i], mp->pfarfield->prpoints[k].ez->data[i]);
                                   j++;
			    }
		    }
	    }


	    for (n=0; n<set->so.nsums; n++)
	    {
		    for (i=0; i<set->sc.nsteps; i++) {
			    mp->out->outsumdata[n][i] = plan->h_sums[set->sc.nsteps*n + i];
		    }
	    
            }
       	    for (n=0; n<set->so.nforces; n++)
	    {
		    for (i=0; i<(3*set->sc.nsteps); i++) {
			    mp->out->outforcedata[n][i] = plan->h_forces[3*set->sc.nsteps*n + i];
		    }
	    
            }
             //for (i=0; i<(6*set->sc.nsteps*set->so.npnts); i++)
           /// {
                   
           //         printf("%d %g\n", i, plan->h_outpointdata[i]);
           // }

	    for (n=0; n<set->so.npnts; n++)
	    {
		    for (i=0; i<(6*set->sc.nsteps); i++) {
			    mp->out->outpointdata[n][i] = plan->h_outpointdata[6*set->sc.nsteps*n + i];
		    }
	    }
    }
    if (set->sc.verbose > 1) printf("done.\n");
    return 0;
 }

void
sv_pool_sync(SvPool *mp, SvSet *set)
{
    gint n;
    gint all;
    SvGpu* gpu = (SvGpu *)mp->gpu;
    //if (set->sc.verbose>1) printf("pool_gpu_sync...\n");
    do {
    //   usleep(100);
       all = 0;
       for (n=0; n<gpu->ngset; n++)
       {
	   all += (int)g_queue_is_empty(gpu->gset[n].aplan->todos);
       }
    } while (all<gpu->ngset);
    //if (set->sc.verbose>1) printf("sync done.\n");
}

void
gpu_init(SvGpuPlan *plan)
{
    //int dev;
    cudaError_t err;

    err = cudaSetDevice(plan->device);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    //err = cudaGetDevice(&dev);
    //if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));


    err = cudaMalloc((void **)((void *)&(plan->d_ex)), plan->xres*plan->yres*plan->zres*sizeof(float));
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    err = cudaMalloc((void **)((void *)&(plan->d_ey)), plan->xres*plan->yres*plan->zres*sizeof(float));
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    err = cudaMalloc((void **)((void *)&(plan->d_ez)), plan->xres*plan->yres*plan->zres*sizeof(float));
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    err = cudaMalloc((void **)((void *)&(plan->d_hx)), plan->xres*plan->yres*plan->zres*sizeof(float));
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    err = cudaMalloc((void **)((void *)&(plan->d_hy)), plan->xres*plan->yres*plan->zres*sizeof(float));
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    err = cudaMalloc((void **)((void *)&(plan->d_hz)), plan->xres*plan->yres*plan->zres*sizeof(float));
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

    if (plan->matmode==SV_MATMODE_FULL || plan->matmode==SV_MATMODE_ELECTRIC)
    { 
        cudaMalloc((void **)((void *)&(plan->d_epsilon)), plan->xres*plan->yres*plan->zres*sizeof(float));
        cudaMalloc((void **)((void *)&(plan->d_sigma)), plan->xres*plan->yres*plan->zres*sizeof(float));
    } 
    if (plan->matmode==SV_MATMODE_FULL || plan->matmode==SV_MATMODE_MAGNETIC) 
    {
        cudaMalloc((void **)((void *)&(plan->d_mu)), plan->xres*plan->yres*plan->zres*sizeof(float));
        cudaMalloc((void **)((void *)&(plan->d_sigast)), plan->xres*plan->yres*plan->zres*sizeof(float));
    }
    if (plan->nmat) {
        cudaMalloc((void **)((void *)&(plan->d_mat)), plan->xres*plan->yres*plan->zres*sizeof(int));
        cudaMalloc((void **)((void *)&(plan->d_mattype)), plan->nmat*sizeof(int));
        cudaMalloc((void **)((void *)&(plan->d_mattab)), SV_GM_N*plan->nmat*sizeof(float));
    }

    if (plan->h_isplrc) {
	    err = cudaMalloc((void **)((void *)&(plan->d_plrcx)), 7*plan->xres*plan->yres*plan->zres*sizeof(float));
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMalloc((void **)((void *)&(plan->d_plrcy)), 7*plan->xres*plan->yres*plan->zres*sizeof(float));
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMalloc((void **)((void *)&(plan->d_plrcz)), 7*plan->xres*plan->yres*plan->zres*sizeof(float));
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->h_isade) {
	    err = cudaMalloc((void **)((void *)&(plan->d_dpx)), 3*plan->xres*plan->yres*plan->zres*sizeof(float));
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMalloc((void **)((void *)&(plan->d_dpxp)), 3*plan->xres*plan->yres*plan->zres*sizeof(float));
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMalloc((void **)((void *)&(plan->d_px)), 6*plan->xres*plan->yres*plan->zres*sizeof(float));
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMalloc((void **)((void *)&(plan->d_pxp)), 6*plan->xres*plan->yres*plan->zres*sizeof(float));
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

	    err = cudaMalloc((void **)((void *)&(plan->d_exp)), plan->xres*plan->yres*plan->zres*sizeof(float));
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMalloc((void **)((void *)&(plan->d_eyp)), plan->xres*plan->yres*plan->zres*sizeof(float));
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMalloc((void **)((void *)&(plan->d_ezp)), plan->xres*plan->yres*plan->zres*sizeof(float));
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }

    
    if (plan->bndx0==SV_BOUNDARY_LIAO) {
        err = cudaMalloc((void **)((void *)&(plan->d_liao_x0)),12*plan->yres*plan->zres*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->bndxn==SV_BOUNDARY_LIAO) {
        err = cudaMalloc((void **)((void *)&(plan->d_liao_xn)),12*plan->yres*plan->zres*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->bndy0==SV_BOUNDARY_LIAO) {
        err = cudaMalloc((void **)((void *)&(plan->d_liao_y0)),12*plan->xres*plan->zres*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->bndyn==SV_BOUNDARY_LIAO) {
        err = cudaMalloc((void **)((void *)&(plan->d_liao_yn)),12*plan->xres*plan->zres*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->bndz0==SV_BOUNDARY_LIAO) {
        err = cudaMalloc((void **)((void *)&(plan->d_liao_z0)),12*plan->xres*plan->yres*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->bndzn==SV_BOUNDARY_LIAO) {
        err = cudaMalloc((void **)((void *)&(plan->d_liao_zn)),12*plan->xres*plan->yres*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }

    err = cudaMalloc((void **)((void *)&(plan->d_bnds)), 6*sizeof(int));
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

    if (plan->bndx0==SV_BOUNDARY_CPML || plan->bndxn==SV_BOUNDARY_CPML || plan->bndy0==SV_BOUNDARY_CPML || plan->bndyn==SV_BOUNDARY_CPML
            || plan->bndz0==SV_BOUNDARY_CPML || plan->bndzn==SV_BOUNDARY_CPML) {

        err = cudaMalloc((void **)((void *)&(plan->d_depths)), 6*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }


    if (plan->bndx0==SV_BOUNDARY_CPML) {

        err = cudaMalloc((void **)((void *)&(plan->d_cpml_be_x0)), plan->cpml_depth_bx0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_ce_x0)), plan->cpml_depth_bx0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_bh_x0)), plan->cpml_depth_bx0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_ch_x0)), plan->cpml_depth_bx0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_kappae_x0)), plan->cpml_depth_bx0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_kappah_x0)), plan->cpml_depth_bx0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }

    if (plan->bndxn==SV_BOUNDARY_CPML) {
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_be_xn)), plan->cpml_depth_bxn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_ce_xn)), plan->cpml_depth_bxn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_bh_xn)), plan->cpml_depth_bxn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_ch_xn)), plan->cpml_depth_bxn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_kappae_xn)), plan->cpml_depth_bxn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_kappah_xn)), plan->cpml_depth_bxn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }

    if (plan->bndy0==SV_BOUNDARY_CPML) {
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_be_y0)), plan->cpml_depth_by0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_ce_y0)), plan->cpml_depth_by0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_bh_y0)), plan->cpml_depth_by0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_ch_y0)), plan->cpml_depth_by0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_kappae_y0)), plan->cpml_depth_by0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_kappah_y0)), plan->cpml_depth_by0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }

    if (plan->bndyn==SV_BOUNDARY_CPML) {
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_be_yn)), plan->cpml_depth_byn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_ce_yn)), plan->cpml_depth_byn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_bh_yn)), plan->cpml_depth_byn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_ch_yn)), plan->cpml_depth_byn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_kappae_yn)), plan->cpml_depth_byn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_kappah_yn)), plan->cpml_depth_byn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }

    if (plan->bndz0==SV_BOUNDARY_CPML) {
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_be_z0)), plan->cpml_depth_bz0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_ce_z0)), plan->cpml_depth_bz0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_bh_z0)), plan->cpml_depth_bz0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_ch_z0)), plan->cpml_depth_bz0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_kappae_z0)), plan->cpml_depth_bz0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_kappah_z0)), plan->cpml_depth_bz0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }

    if (plan->bndzn==SV_BOUNDARY_CPML) {
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_be_zn)), plan->cpml_depth_bzn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_ce_zn)), plan->cpml_depth_bzn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_bh_zn)), plan->cpml_depth_bzn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_ch_zn)), plan->cpml_depth_bzn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_kappae_zn)), plan->cpml_depth_bzn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_kappah_zn)), plan->cpml_depth_bzn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }

    if (plan->bndx0==SV_BOUNDARY_CPML) {
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_p_x0)), 4*plan->yres*plan->zres*plan->cpml_depth_bx0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->bndxn==SV_BOUNDARY_CPML) {
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_p_xn)), 4*plan->yres*plan->zres*plan->cpml_depth_bxn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->bndy0==SV_BOUNDARY_CPML) {    
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_p_y0)), 4*plan->xres*plan->zres*plan->cpml_depth_by0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->bndyn==SV_BOUNDARY_CPML) {
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_p_yn)), 4*plan->xres*plan->zres*plan->cpml_depth_byn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->bndz0==SV_BOUNDARY_CPML) {
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_p_z0)), 4*plan->xres*plan->yres*plan->cpml_depth_bz0*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->bndzn==SV_BOUNDARY_CPML) {
        err = cudaMalloc((void **)((void *)&(plan->d_cpml_p_zn)), 4*plan->xres*plan->yres*plan->cpml_depth_bzn*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }


    if (plan->h_tsfset) {
        err = cudaMalloc((void **)((void *)&(plan->d_tsf_jpool_e)),plan->tsf_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsf_jpool_h)),plan->tsf_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsf_jpool_epsilon)),plan->tsf_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsf_jpool_mu)),plan->tsf_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsf_jpool_sigma)),plan->tsf_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsf_jpool_sigast)),plan->tsf_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsfset)),29*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsf_jpvals)),8*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->h_tsffset) {
        err = cudaMalloc((void **)((void *)&(plan->d_tsff_jpool_e)),plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsff_jpool_h)),plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsff_jpool_epsilon)),plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsff_jpool_mu)),plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsff_jpool_sigma)),plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsff_jpool_sigast)),plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsffset)),16*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsff_jpvals)),8*plan->tsff_nxm*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
 }
     if (plan->h_sfset) {
        err = cudaMalloc((void **)((void *)&(plan->d_sf_jpool_e)),plan->sf_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sf_jpool_h)),plan->sf_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sf_jpool_epsilon)),plan->sf_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sf_jpool_mu)),plan->sf_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sf_jpool_sigma)),plan->sf_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sf_jpool_sigast)),plan->sf_jpool_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sfset)),10*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sf_jpvals)),8*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
     if (plan->ff_size!=0) {
        err = cudaMalloc((void **)((void *)&(plan->d_iset)), plan->iset_size*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_ff_ex)), plan->ff_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_ff_ey)), plan->ff_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_ff_ez)), plan->ff_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_ff_hlp_ex)), plan->ff_size*plan->nhlps*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_ff_hlp_ey)), plan->ff_size*plan->nhlps*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_ff_hlp_ez)), plan->ff_size*plan->nhlps*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
     if (plan->pff_size!=0) {
        err = cudaMalloc((void **)((void *)&(plan->d_piset)), plan->piset_size*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_pff_ex)), plan->pff_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_pff_ey)), plan->pff_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_pff_ez)), plan->pff_size*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	if (plan->h_piset[10]<GPU_PERIODIC_HLPS_LIMIT) { printf("cpt\n");
		err = cudaMalloc((void **)((void *)&(plan->d_pff_hlp_ex)), plan->pff_size*(plan->h_piset[8]-plan->h_piset[6])*(plan->h_piset[9]-plan->h_piset[7])*sizeof(float));
		if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
		err = cudaMalloc((void **)((void *)&(plan->d_pff_hlp_ey)), plan->pff_size*(plan->h_piset[8]-plan->h_piset[6])*(plan->h_piset[9]-plan->h_piset[7])*sizeof(float));
		if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
		err = cudaMalloc((void **)((void *)&(plan->d_pff_hlp_ez)), plan->pff_size*(plan->h_piset[8]-plan->h_piset[6])*(plan->h_piset[9]-plan->h_piset[7])*sizeof(float));
		if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	}
     }

    if (plan->nsums>0)
    {
        err = cudaMalloc((void **)((void *)&(plan->d_sums)),plan->nsums*plan->nsteps*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_epsilon)),plan->nsums*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_mu)),plan->nsums*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_sigma)),plan->nsums*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_sigast)),plan->nsums*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_mode)),plan->nsums*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_i0)),plan->nsums*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_i1)),plan->nsums*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_j0)),plan->nsums*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_j1)),plan->nsums*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_k0)),plan->nsums*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_k1)),plan->nsums*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

        if (plan->dir==2) {
            err = cudaMalloc((void **)((void *)&(plan->d_sum_accumulator)),plan->xres*plan->yres*plan->nsums*sizeof(float));
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

        }
    }
    if (plan->nforces>0)
    {
        err = cudaMalloc((void **)((void *)&(plan->d_forces)),3*plan->nforces*plan->nsteps*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_force_i0)),plan->nforces*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_force_i1)),plan->nforces*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_force_j0)),plan->nforces*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_force_j1)),plan->nforces*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_force_k0)),plan->nforces*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_force_k1)),plan->nforces*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
     if (plan->noutpoints>0) {
        err = cudaMalloc((void **)((void *)&(plan->d_outpointdata)),6*plan->noutpoints*plan->nsteps*sizeof(float));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_outpoint_pos)),3*plan->noutpoints*sizeof(int));
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }

}

void
gpu_copyto(SvGpuPlan *plan, SvGpuSync what)
{
    
    cudaError_t err;

//    int dev;
//    err = cudaGetDevice(&dev);
//    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    //printf("thread reported that device %d is set\n", dev);


    if (what == SV_SYNC_ALL || what == SV_SYNC_EH || what == SV_SYNC_E) {
	    err = cudaMemcpy(plan->d_ex, plan->h_ex, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_ey, plan->h_ey, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_ez, plan->h_ez, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (what == SV_SYNC_ALL || what == SV_SYNC_EH || what == SV_SYNC_H) {
	    err = cudaMemcpy(plan->d_hx, plan->h_hx, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_hy, plan->h_hy, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_hz, plan->h_hz, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (what == SV_SYNC_ALL && (plan->matmode==SV_MATMODE_FULL || plan->matmode==SV_MATMODE_ELECTRIC))
    {
	    err = cudaMemcpy(plan->d_epsilon, plan->h_epsilon, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_sigma, plan->h_sigma, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    } 
    if (what == SV_SYNC_ALL && (plan->matmode==SV_MATMODE_FULL || plan->matmode==SV_MATMODE_MAGNETIC)) 
    {
	    err = cudaMemcpy(plan->d_mu, plan->h_mu, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_sigast, plan->h_sigast, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }

    if (what == SV_SYNC_ALL && plan->noutpoints>0) {
  	    err = cudaMemcpy(plan->d_outpointdata, plan->h_outpointdata, 6*plan->noutpoints*plan->nsteps*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("sakrble, %s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_outpoint_pos, plan->h_outpoint_pos, 3*plan->noutpoints*sizeof(int), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }

    if (what == SV_SYNC_ALL && plan->h_tsfset!=NULL) {
	    err = cudaMemcpy(plan->d_tsf_jpool_e, plan->h_tsf_jpool_e, plan->tsf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
   	    err = cudaMemcpy(plan->d_tsf_jpool_h, plan->h_tsf_jpool_h, plan->tsf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
   	    err = cudaMemcpy(plan->d_tsf_jpool_epsilon, plan->h_tsf_jpool_epsilon, plan->tsf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
   	    err = cudaMemcpy(plan->d_tsf_jpool_sigma, plan->h_tsf_jpool_sigma, plan->tsf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_tsf_jpool_mu, plan->h_tsf_jpool_mu, plan->tsf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    	err = cudaMemcpy(plan->d_tsf_jpool_sigast, plan->h_tsf_jpool_sigast, plan->tsf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    	err = cudaMemcpy(plan->d_tsfset, plan->h_tsfset, 29*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    	err = cudaMemcpy(plan->d_tsf_jpvals, plan->h_tsf_jpvals, 8*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    
    }
    if (what == SV_SYNC_ALL && plan->h_tsffset!=NULL) {
	    err = cudaMemcpy(plan->d_tsff_jpool_e, plan->h_tsff_jpool_e, plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_tsff_jpool_h, plan->h_tsff_jpool_h, plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_tsff_jpool_epsilon, plan->h_tsff_jpool_epsilon, plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_tsff_jpool_sigma, plan->h_tsff_jpool_sigma, plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_tsff_jpool_mu, plan->h_tsff_jpool_mu, plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_tsff_jpool_sigast, plan->h_tsff_jpool_sigast, plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_tsffset, plan->h_tsffset, 16*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_tsff_jpvals, plan->h_tsff_jpvals, plan->tsff_nxm*8*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
     if (what == SV_SYNC_ALL && plan->h_sfset!=NULL) {
	    err = cudaMemcpy(plan->d_sf_jpool_e, plan->h_sf_jpool_e, plan->sf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
   	    err = cudaMemcpy(plan->d_sf_jpool_h, plan->h_sf_jpool_h, plan->sf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
   	    err = cudaMemcpy(plan->d_sf_jpool_epsilon, plan->h_sf_jpool_epsilon, plan->sf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
   	    err = cudaMemcpy(plan->d_sf_jpool_sigma, plan->h_sf_jpool_sigma, plan->sf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sf_jpool_mu, plan->h_sf_jpool_mu, plan->sf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    	err = cudaMemcpy(plan->d_sf_jpool_sigast, plan->h_sf_jpool_sigast, plan->sf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    	err = cudaMemcpy(plan->d_sfset, plan->h_sfset, 10*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    	err = cudaMemcpy(plan->d_sf_jpvals, plan->h_sf_jpvals, 8*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
      
    }
    if (what == SV_SYNC_ALL && plan->ff_size!=0) {
	    err = cudaMemcpy(plan->d_iset, plan->h_iset, plan->iset_size*sizeof(int), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_ff_ex, plan->h_ff_ex, plan->ff_size*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_ff_ey, plan->h_ff_ey, plan->ff_size*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_ff_ez, plan->h_ff_ez, plan->ff_size*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_ff_hlp_ex, plan->h_ff_hlp_ex, plan->ff_size*plan->nhlps*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_ff_hlp_ey, plan->h_ff_hlp_ey, plan->ff_size*plan->nhlps*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_ff_hlp_ez, plan->h_ff_hlp_ez, plan->ff_size*plan->nhlps*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

    }
    if (what == SV_SYNC_ALL && plan->pff_size!=0) {
	    err = cudaMemcpy(plan->d_piset, plan->h_piset, plan->piset_size*sizeof(int), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_pff_ex, plan->h_pff_ex, plan->pff_size*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_pff_ey, plan->h_pff_ey, plan->pff_size*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_pff_ez, plan->h_pff_ez, plan->pff_size*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    if (plan->h_piset[10]<GPU_PERIODIC_HLPS_LIMIT) {
		    err = cudaMemcpy(plan->d_pff_hlp_ex, plan->h_pff_hlp_ex, plan->pff_size*(plan->h_piset[8]-plan->h_piset[6])*(plan->h_piset[9]-plan->h_piset[7])*sizeof(float), cudaMemcpyHostToDevice);
		    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
		    err = cudaMemcpy(plan->d_pff_hlp_ey, plan->h_pff_hlp_ey, plan->pff_size*(plan->h_piset[8]-plan->h_piset[6])*(plan->h_piset[9]-plan->h_piset[7])*sizeof(float), cudaMemcpyHostToDevice);
		    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
		    err = cudaMemcpy(plan->d_pff_hlp_ez, plan->h_pff_hlp_ez, plan->pff_size*(plan->h_piset[8]-plan->h_piset[6])*(plan->h_piset[9]-plan->h_piset[7])*sizeof(float), cudaMemcpyHostToDevice);
		    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    }


    }


    if (what == SV_SYNC_ALL && plan->nsums > 0) {
        err = cudaMemcpy(plan->d_sums, plan->h_sums, plan->nsums*plan->nsteps*sizeof(float), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_epsilon, plan->h_sum_epsilon, plan->nsums*sizeof(float), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_mu, plan->h_sum_mu, plan->nsums*sizeof(float), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_sigma, plan->h_sum_sigma, plan->nsums*sizeof(float), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_sigast, plan->h_sum_sigast, plan->nsums*sizeof(float), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_mode, plan->h_sum_mode, plan->nsums*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_i0, plan->h_sum_i0, plan->nsums*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_i1, plan->h_sum_i1, plan->nsums*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_j0, plan->h_sum_j0, plan->nsums*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_j1, plan->h_sum_j1, plan->nsums*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_k0, plan->h_sum_k0, plan->nsums*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_k1, plan->h_sum_k1, plan->nsums*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (what == SV_SYNC_ALL && plan->nforces > 0) {
        err = cudaMemcpy(plan->d_forces, plan->h_forces, 3*plan->nforces*plan->nsteps*sizeof(float), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_force_i0, plan->h_force_i0, plan->nforces*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_force_i1, plan->h_force_i1, plan->nforces*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_force_j0, plan->h_force_j0, plan->nforces*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_force_j1, plan->h_force_j1, plan->nforces*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_force_k0, plan->h_force_k0, plan->nforces*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_force_k1, plan->h_force_k1, plan->nforces*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
     if (what == SV_SYNC_ALL && plan->nmat>0) {
        err = cudaMemcpy(plan->d_mat, plan->h_mat, plan->xres*plan->yres*plan->zres*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_mattype, plan->h_mattype, plan->nmat*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_mattab, plan->h_mattab, SV_GM_N*plan->nmat*sizeof(float), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

    }
    if (what == SV_SYNC_ALL) {

        err = cudaMemcpy(plan->d_bnds, plan->h_bnds, 6*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
         
        if (plan->bndx0==SV_BOUNDARY_CPML || plan->bndxn==SV_BOUNDARY_CPML || plan->bndy0==SV_BOUNDARY_CPML || plan->bndyn==SV_BOUNDARY_CPML
            || plan->bndz0==SV_BOUNDARY_CPML || plan->bndzn==SV_BOUNDARY_CPML)
	{
            err = cudaMemcpy(plan->d_depths, plan->h_depths, 6*sizeof(int), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

	}

	if (plan->bndx0==SV_BOUNDARY_CPML) {
            err = cudaMemcpy(plan->d_cpml_be_x0, plan->h_cpml_be_x0, plan->cpml_depth_bx0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_ce_x0, plan->h_cpml_ce_x0, plan->cpml_depth_bx0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_bh_x0, plan->h_cpml_bh_x0, plan->cpml_depth_bx0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_ch_x0, plan->h_cpml_ch_x0, plan->cpml_depth_bx0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_kappae_x0, plan->h_cpml_kappae_x0, plan->cpml_depth_bx0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_kappah_x0, plan->h_cpml_kappah_x0, plan->cpml_depth_bx0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

            err = cudaMemcpy(plan->d_cpml_p_x0, plan->h_cpml_p_x0, 4*plan->yres*plan->zres*plan->cpml_depth_bx0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

        }

        if (plan->bndxn==SV_BOUNDARY_CPML) {
            err = cudaMemcpy(plan->d_cpml_be_xn, plan->h_cpml_be_xn, plan->cpml_depth_bxn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_ce_xn, plan->h_cpml_ce_xn, plan->cpml_depth_bxn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_bh_xn, plan->h_cpml_bh_xn, plan->cpml_depth_bxn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_ch_xn, plan->h_cpml_ch_xn, plan->cpml_depth_bxn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_kappae_xn, plan->h_cpml_kappae_xn, plan->cpml_depth_bxn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_kappah_xn, plan->h_cpml_kappah_xn, plan->cpml_depth_bxn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

            err = cudaMemcpy(plan->d_cpml_p_xn, plan->h_cpml_p_xn, 4*plan->yres*plan->zres*plan->cpml_depth_bxn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

        }

        if (plan->bndy0==SV_BOUNDARY_CPML) {
            err = cudaMemcpy(plan->d_cpml_be_y0, plan->h_cpml_be_y0, plan->cpml_depth_by0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_ce_y0, plan->h_cpml_ce_y0, plan->cpml_depth_by0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_bh_y0, plan->h_cpml_bh_y0, plan->cpml_depth_by0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_ch_y0, plan->h_cpml_ch_y0, plan->cpml_depth_by0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_kappae_y0, plan->h_cpml_kappae_y0, plan->cpml_depth_by0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_kappah_y0, plan->h_cpml_kappah_y0, plan->cpml_depth_by0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

            err = cudaMemcpy(plan->d_cpml_p_y0, plan->h_cpml_p_y0, 4*plan->xres*plan->zres*plan->cpml_depth_by0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

        }

        if (plan->bndyn==SV_BOUNDARY_CPML) {
            err = cudaMemcpy(plan->d_cpml_be_yn, plan->h_cpml_be_yn, plan->cpml_depth_byn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_ce_yn, plan->h_cpml_ce_yn, plan->cpml_depth_byn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_bh_yn, plan->h_cpml_bh_yn, plan->cpml_depth_byn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_ch_yn, plan->h_cpml_ch_yn, plan->cpml_depth_byn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_kappae_yn, plan->h_cpml_kappae_yn, plan->cpml_depth_byn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_kappah_yn, plan->h_cpml_kappah_yn, plan->cpml_depth_byn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

            err = cudaMemcpy(plan->d_cpml_p_yn, plan->h_cpml_p_yn, 4*plan->xres*plan->zres*plan->cpml_depth_byn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

        }

        if (plan->bndz0==SV_BOUNDARY_CPML) {
            err = cudaMemcpy(plan->d_cpml_be_z0, plan->h_cpml_be_z0, plan->cpml_depth_bz0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_ce_z0, plan->h_cpml_ce_z0, plan->cpml_depth_bz0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_bh_z0, plan->h_cpml_bh_z0, plan->cpml_depth_bz0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_ch_z0, plan->h_cpml_ch_z0, plan->cpml_depth_bz0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_kappae_z0, plan->h_cpml_kappae_z0, plan->cpml_depth_bz0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_kappah_z0, plan->h_cpml_kappah_z0, plan->cpml_depth_bz0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

            err = cudaMemcpy(plan->d_cpml_p_z0, plan->h_cpml_p_z0, 4*plan->xres*plan->yres*plan->cpml_depth_bz0*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

        }

        if (plan->bndzn==SV_BOUNDARY_CPML) {
            err = cudaMemcpy(plan->d_cpml_be_zn, plan->h_cpml_be_zn, plan->cpml_depth_bzn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_ce_zn, plan->h_cpml_ce_zn, plan->cpml_depth_bzn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_bh_zn, plan->h_cpml_bh_zn, plan->cpml_depth_bzn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_ch_zn, plan->h_cpml_ch_zn, plan->cpml_depth_bzn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_kappae_zn, plan->h_cpml_kappae_zn, plan->cpml_depth_bzn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_cpml_kappah_zn, plan->h_cpml_kappah_zn, plan->cpml_depth_bzn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

            err = cudaMemcpy(plan->d_cpml_p_zn, plan->h_cpml_p_zn, 4*plan->xres*plan->yres*plan->cpml_depth_bzn*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

        }

    }

    err = cudaThreadSynchronize();
    if (err != cudaSuccess) printf("E returned \"%s\"\n", cudaGetErrorString(err));
    //printf("gpu_copyto finished\n");

}

void
gpu_copyfrom(SvGpuPlan *plan)
{
    cudaError_t err;
    int dev;
    err = cudaGetDevice(&dev);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));


    /*collect hlps for farfield, if there are any*/
      if (plan->nhlps) {
	    err = wrap_ffKernel_gethlps(plan);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    } 
    /*we cannot always use hlps for periodic farfield, this eats too much memory*/
    if (plan->pff_size!=0) {
          if (plan->h_piset[10]<GPU_PERIODIC_HLPS_LIMIT) {   //FIXME this should be related to card abilities
printf("collect copyfrom\n");
	         err = wrap_pffKernel_gethlps(plan);
             if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
         }
    }


    err = cudaMemcpy(plan->h_ex, plan->d_ex, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    err = cudaMemcpy(plan->h_ey, plan->d_ey, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    err = cudaMemcpy(plan->h_ez, plan->d_ez, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    err = cudaMemcpy(plan->h_hx, plan->d_hx, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    err = cudaMemcpy(plan->h_hy, plan->d_hy, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    err = cudaMemcpy(plan->h_hz, plan->d_hz, plan->xres*plan->yres*plan->zres*sizeof(float), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

    if (plan->ff_size!=0) {
	    err = cudaMemcpy(plan->h_iset, plan->d_iset, plan->iset_size*sizeof(int), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->h_ff_ex, plan->d_ff_ex, plan->ff_size*sizeof(float), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->h_ff_ey, plan->d_ff_ey, plan->ff_size*sizeof(float), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->h_ff_ez, plan->d_ff_ez, plan->ff_size*sizeof(float), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->pff_size!=0) {
 printf("copy back %d data\n", plan->pff_size);
	    err = cudaMemcpy(plan->h_piset, plan->d_piset, plan->piset_size*sizeof(int), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->h_pff_ex, plan->d_pff_ex, plan->pff_size*sizeof(float), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->h_pff_ey, plan->d_pff_ey, plan->pff_size*sizeof(float), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->h_pff_ez, plan->d_pff_ez, plan->pff_size*sizeof(float), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }

    /*
    err = cudaMemcpy(plan->h_tsff_jpool_e, plan->d_tsff_jpool_e, plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    int i;
    for (i=0; i<plan->tsff_nxm*plan->tsff_jpool_size; i++) printf("%d %g\n", i, plan->h_tsff_jpool_e[i]);

     err = cudaMemcpy(plan->h_tsff_jpvals, plan->d_tsff_jpvals, plan->tsff_nxm*8*sizeof(float), cudaMemcpyDeviceToHost);
     if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
   for (i=0; i<plan->tsff_nxm*8; i++) printf("jp %d %g\n", i, plan->h_tsff_jpvals[i]);
  */ 

    if (plan->noutpoints > 0) {
  	    err = cudaMemcpy(plan->h_outpointdata, plan->d_outpointdata, 6*plan->noutpoints*plan->nsteps*sizeof(float), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->h_outpoint_pos, plan->d_outpoint_pos, 3*plan->noutpoints*sizeof(int), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }

    if (plan->nsums > 0) {
	    err = cudaMemcpy(plan->h_sums, plan->d_sums, plan->nsums*plan->nsteps*sizeof(float), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->nforces > 0) {
	    err = cudaMemcpy(plan->h_forces, plan->d_forces, 3*plan->nforces*plan->nsteps*sizeof(float), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    }
		
    err = cudaThreadSynchronize();
    if (err != cudaSuccess) printf("E returned \"%s\"\n", cudaGetErrorString(err));

}

void
gpu_ystep_e(SvGpuPlan *plan)
{
    cudaError_t err = wrap_eKernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_ystep_h(SvGpuPlan *plan)
{
    cudaError_t err = wrap_hKernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_liao_cpy(SvGpuPlan *plan)
{
    cudaError_t err;
    
    plan->step += 1;
    err = wrap_liao_cpybnd(plan);
    
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_liao_bnd(SvGpuPlan *plan)
{
    cudaError_t err = wrap_liao_applybnd(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_cpml_estep(SvGpuPlan *plan)
{
    cudaError_t err = wrap_cpml_estep(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_cpml_hstep(SvGpuPlan *plan)
{
    cudaError_t err = wrap_cpml_hstep(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_srce(SvGpuPlan *plan)
{
    cudaError_t err = wrap_srceKernel(plan, plan->source_i, plan->source_j, plan->source_k, plan->source_ex, plan->source_ey, plan->source_ez);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_tsf_estep(SvGpuPlan *plan)
{
    cudaError_t err = wrap_tsf_e_Kernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_tsf_hstep(SvGpuPlan *plan)
{
    cudaError_t err = wrap_tsf_h_Kernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_tsff_estep(SvGpuPlan *plan)
{
    cudaError_t err = wrap_tsff_e_Kernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_tsff_hstep(SvGpuPlan *plan)
{
    cudaError_t err = wrap_tsff_h_Kernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_mbndx_e(SvGpuPlan *plan)
{
    cudaError_t err = wrap_mbndx_e_Kernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_mbndy_e(SvGpuPlan *plan)
{
    cudaError_t err = wrap_mbndy_e_Kernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_mbndz_e(SvGpuPlan *plan)
{
    cudaError_t err = wrap_mbndz_e_Kernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_mbndx_h(SvGpuPlan *plan)
{
    cudaError_t err = wrap_mbndx_h_Kernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_mbndy_h(SvGpuPlan *plan)
{
    cudaError_t err = wrap_mbndy_h_Kernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_mbndz_h(SvGpuPlan *plan)
{
    cudaError_t err = wrap_mbndz_h_Kernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_tsf_jstep(SvGpuPlan *plan)
{
    cudaError_t err = wrap_tsf_jstepKernel(plan, plan->source_ex);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_tsff_jstep(SvGpuPlan *plan)
{
    cudaError_t err = wrap_tsff_jstepKernel(plan, plan->source_ex);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_sf_jstep(SvGpuPlan *plan)
{
    //if (plan->step == 0) {   //one step more to sync with tsf source
    //   cudaError_t err = wrap_sf_jstepKernel(plan, plan->source_ex);
    //if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));

    //}
     cudaError_t err = wrap_sf_jstepKernel(plan, plan->source_ex);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_sum(SvGpuPlan *plan)
{
    cudaError_t err = wrap_fastsumKernel(plan);//wrap_fastsumKernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_force(SvGpuPlan *plan)
{
    cudaError_t err = wrap_forceKernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_opnt(SvGpuPlan *plan)
{
    cudaError_t err = wrap_outpointKernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_ff(SvGpuPlan *plan)
{
    cudaError_t err = wrap_ffKernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void
gpu_ffa(SvGpuPlan *plan)
{
    cudaError_t err = wrap_ffKernel_hlps(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

/*in contrast to full farfield, here we always use helpers*/
void
gpu_pff(SvGpuPlan *plan)
{
    cudaError_t err;

    if (plan->h_piset[10]<GPU_PERIODIC_HLPS_LIMIT)
       err = wrap_pffKernel_hlps(plan);
    else 
       err = wrap_pffKernel(plan);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
}

void* gpu_main(void *data)
{
    SvGpuAPlan *aplan = (SvGpuAPlan *)data;
    SvGpuPlan *plan = aplan->plan;
    SvGTType todo;

    while (1)
    {
	/*check whether there is a request for action*/
	if (g_queue_is_empty(aplan->todos)) continue;
	else todo = (SvGTType)GPOINTER_TO_INT(g_queue_peek_tail(aplan->todos));

	switch ((int)todo) {
	    case SV_GT_INIT:
		//printf("Init request on device %d\n", plan->device);
		gpu_init(plan);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_COPYTO:
		//printf("Copyto on device %d\n", plan->device);
		gpu_copyto(plan, SV_SYNC_ALL);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_COPYFROM:
		//printf("Copyfrom on device %d\n", plan->device);
		gpu_copyfrom(plan);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_LIAO_CPY:
		//printf("Liao copy on device %d\n", plan->device);
		gpu_liao_cpy(plan);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_SRCE_POINT:
		//printf("E source on device %d\n", plan->device);
		gpu_srce(plan);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_LIAO_BND:
		//printf("Liao bnd on device %d\n", plan->device);
		gpu_liao_bnd(plan);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_CPML_ESTEP:
		//printf("TSF on device %d\n", plan->device);
		gpu_cpml_estep(plan);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_CPML_HSTEP:
		//printf("TSF on device %d\n", plan->device);
		gpu_cpml_hstep(plan);
		g_queue_pop_tail(aplan->todos);
		break;
		case SV_GT_YSTEPH:
		//printf("H Ystep on device %d\n", plan->device);
		gpu_ystep_h(plan);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_YSTEPE:
		//printf("E Ystep on device %d\n", plan->device);
		gpu_ystep_e(plan);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_TSF_ESTEP:
		//printf("TSF on device %d\n", plan->device);
		gpu_tsf_estep(plan);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_TSF_HSTEP:
		//printf("TSF on device %d\n", plan->device);
		gpu_tsf_hstep(plan);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_TSFF_ESTEP:
		//printf("TSF on device %d\n", plan->device);
		gpu_tsff_estep(plan);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_TSFF_HSTEP:
		//printf("TSF on device %d\n", plan->device);
		gpu_tsff_hstep(plan);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_TSF_JSTEP:
		//printf("J step on device %d\n", plan->device);
		gpu_tsf_jstep(plan);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_TSFF_JSTEP:
		//printf("J step on device %d\n", plan->device);
		gpu_tsff_jstep(plan);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_SF_JSTEP:
		//printf("J step on device %d\n", plan->device);
		gpu_sf_jstep(plan);
		g_queue_pop_tail(aplan->todos);
		break;
	    case SV_GT_MBNDX_E:
		//printf("MBNDX on device %d\n", plan->device);
		gpu_mbndx_e(plan);
		g_queue_pop_tail(aplan->todos);
		break;		
	    case SV_GT_MBNDY_E:
		//printf("MBNDX on device %d\n", plan->device);
		gpu_mbndy_e(plan);
		g_queue_pop_tail(aplan->todos);
		break;		
	    case SV_GT_MBNDZ_E:
		//printf("MBNDX on device %d\n", plan->device);
		gpu_mbndz_e(plan);
		g_queue_pop_tail(aplan->todos);
		break;		
	    case SV_GT_MBNDX_H:
		//printf("MBNDX on device %d\n", plan->device);
		gpu_mbndx_h(plan);
		g_queue_pop_tail(aplan->todos);
		break;		
	    case SV_GT_MBNDY_H:
		//printf("MBNDX on device %d\n", plan->device);
		gpu_mbndy_h(plan);
		g_queue_pop_tail(aplan->todos);
		break;		
	    case SV_GT_MBNDZ_H:
		//printf("MBNDX on device %d\n", plan->device);
		gpu_mbndz_h(plan);
		g_queue_pop_tail(aplan->todos);
		break;		        
	    case SV_GT_OPNT:
		gpu_opnt(plan);
		g_queue_pop_tail(aplan->todos);
		break;		        
            case SV_GT_SUM:
		//printf("Sum on device %d\n", plan->device);
		gpu_sum(plan);
		g_queue_pop_tail(aplan->todos);
		break;		        
       	    case SV_GT_FF:
		gpu_ff(plan);
		g_queue_pop_tail(aplan->todos);
		break;		        
       	    case SV_GT_PFF:
		gpu_pff(plan);
		g_queue_pop_tail(aplan->todos);
		break;		        
       	    case SV_GT_FFA:
		gpu_ffa(plan);
		g_queue_pop_tail(aplan->todos);
		break;		        
       	    case SV_GT_FORCE:
		gpu_force(plan);
		g_queue_pop_tail(aplan->todos);
		break;		        
         case SV_GT_NOTHING:
		//printf("Nothing on device %d\n", plan->device);
//		sleep(1);
		g_queue_pop_tail(aplan->todos);
		break;
            default:
		break;
				    
	}
    }
}


/*allocate structures and start worker threads*/
SvGpu*
sv_gpu_new(SvPool *mp, SvSet *set)
{
    int i, j, k, n, splitsize;
    cudaError_t err;
    GError *gerr = NULL;
    struct cudaDeviceProp prop;
    SvGpuPlan *plan;
    SvGpu* gpu;
    int maxthreads[MAX_GPUS];

    int ngpu = 0, ngset=0, nprescribed, sel_gpu[MAX_GPUS];

    if (set->sc.usegpu == 0) 
        return NULL;

    gpu = (SvGpu*)g_malloc(sizeof(SvGpu));

    if (set->sc.devicequery) {
	    printf("running getdevicecount\n");
	    err = cudaGetDeviceCount(&ngpu);
	    printf("finished\n");
	    if (err != cudaSuccess) {
                fprintf(stderr, "Error: %s\n", cudaGetErrorString(err));
		ngpu = 0;
	    }
	    if (set->sc.verbose) 
                printf("Found %d GPUs\n", ngpu);


	    for (i=0; i<ngpu; i++) {
		    err = cudaGetDeviceProperties(&prop, i);
		    if (set->sc.verbose>1) {
			    printf("The Properties of the Device with ID %d are\n", i);
			    printf("Device Name             : %s\n", prop.name);
			    printf("Device Memory Size      : %lu\n", (gulong) prop.totalGlobalMem);
			    printf("Block Shared memory size: %d\n", (gint)prop.sharedMemPerBlock);
			    printf("Max grid size           : %dx%dx%d\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
			    printf("Max threads dim         : %dx%dx%d\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
			    printf("Revision                : %d.%d\n", prop.major, prop.minor);
			    printf("Multiprocessors         : %d\n", prop.multiProcessorCount);
			    printf("Clock rate              : %d\n", prop.clockRate);
		    }
		    maxthreads[i] = prop.maxThreadsDim[0];
	    }

	    if (set->sc.verbose) 
            printf("Gpu found and will be used\n");
    } else {
        ngpu = set->sc.usegpu;   /*no test for gpus, just believe user, unsafe*/
      	if (set->sc.verbose) 
            printf("Warning: GPU detection skipped by your request (GPU_QUERY=0), assuming having at least %d GPUs\n", ngpu);
    }

    if (ngpu) 
        set->plan.gpumode = SV_GPUMODE_FULL;

    /*control settings*/
    nprescribed = 0;

    for (i=0; i<MAX_GPUS; i++) {
        sel_gpu[i] = 0;
        nprescribed += set->sc.ugpu[i];
    }
    if (set->sc.verbose > 1) 
        printf("prescribed %d GPU positions\n", nprescribed);
    if (nprescribed==0 && set->sc.usegpu>0) { /*just pick first few devices*/
        if (set->sc.usegpu <= ngpu) 
            ngset = set->sc.usegpu;
        for (i=0; i<ngset; i++) 
            sel_gpu[i] = i;
    } else if (nprescribed>0) {
        ngset = 0;
        for (i=0; i<MAX_GPUS; i++) {
            if (set->sc.ugpu[i] && i<ngpu) 
                sel_gpu[ngset++] = i;
            else if (set->sc.ugpu[i] && i>=ngpu)
                fprintf(stderr, "Error: requesting GPU outside physical range (No. %d of %d)\n", i, ngpu);
        }
    }

    if (set->sc.verbose > 1) printf("%d devices used: %d%d%d%d %d%d%d%d %d%d%d%d %d%d%d%d\n", 
		    ngset, sel_gpu[0], sel_gpu[1], sel_gpu[2], sel_gpu[3],
		    sel_gpu[4], sel_gpu[5], sel_gpu[6], sel_gpu[7],
		    sel_gpu[8], sel_gpu[9], sel_gpu[10], sel_gpu[11],
            sel_gpu[12], sel_gpu[13], sel_gpu[14], sel_gpu[15]);

    /*perform allocation*/
    gpu->ngset = ngset;
    gpu->gset = (SvGpuTSet *)g_malloc(ngset*sizeof(SvGpuTSet));
    for (n=0; n<ngset; n++)
    {
        if (set->sc.verbose > 1) 
            printf("Creating plan for card %d\n", sel_gpu[n]);
        gpu->gset[n].aplan = (SvGpuAPlan *)g_malloc(sizeof(SvGpuAPlan));
        gpu->gset[n].aplan->plan = (SvGpuPlan *)g_malloc(sizeof(SvGpuPlan));
        gpu->gset[n].aplan->todos = g_queue_new();
        gpu->gset[n].aplan->doing = SV_GT_NOTHING;
        gpu->gset[n].aplan->mutex = g_mutex_new();

        plan = gpu->gset[n].aplan->plan;
        plan->step = -1;
        plan->matmode = set->plan.matmode;
        plan->device = sel_gpu[n];
        plan->maxthreads = maxthreads[n];

        plan->h_bnds = (int *)g_malloc(6*sizeof(int));
        plan->h_bnds[0] = set->sb.bx0; 
        plan->h_bnds[1] = set->sb.bxn;
        plan->h_bnds[2] = set->sb.by0;
        plan->h_bnds[3] = set->sb.byn;
        plan->h_bnds[4] = set->sb.bz0;
        plan->h_bnds[5] = set->sb.bzn;

        plan->bndx0 = set->sb.bx0;
        plan->bndxn = set->sb.bxn;
        plan->bndy0 = set->sb.by0;
        plan->bndyn = set->sb.byn;
        plan->bndz0 = set->sb.bz0;
        plan->bndzn = set->sb.bzn;

        plan->zres = set->sp.zres;
        plan->yres = set->sp.yres;
        if (ngset==1) {
            plan->xfrom = 0;
            plan->xto = set->sp.xres;
        } else {
            plan->xfrom = MAX(0, (n*set->sp.xres/ngset) - 1);
            plan->xto = MIN(set->sp.xres, ((n+1)*set->sp.xres/ngset) + 1);
        }
        plan->xres = plan->xto - plan->xfrom;
        plan->dir = 2;
        plan->dx = (float)set->sp.dx;
        plan->dy = (float)set->sp.dy;
        plan->dz = (float)set->sp.dz;
        plan->dt = (float)set->plan.dt;
        plan->nsteps = set->sc.nsteps;
        plan->nhlps = 0;

        if (set->sc.verbose > 1) 
            printf("Allocating %dx%dx%d\n", plan->xres, plan->yres, plan->zres);
        /*
           if (ugpu>= 0 || ngpu==1) plan[n].bnd = 3;
           else if (n==0) plan[n].bnd = 1;
           else if (n==(ngpu-1)) plan[n].bnd = 2;
           else plan[n].bnd = 0;
           */
        plan->h_ex = (float *)g_malloc(plan->xres*plan->yres*plan->zres*sizeof(float));
        plan->h_ey = (float *)g_malloc(plan->xres*plan->yres*plan->zres*sizeof(float));
        plan->h_ez = (float *)g_malloc(plan->xres*plan->yres*plan->zres*sizeof(float));
        plan->h_hx = (float *)g_malloc(plan->xres*plan->yres*plan->zres*sizeof(float));
        plan->h_hy = (float *)g_malloc(plan->xres*plan->yres*plan->zres*sizeof(float));
        plan->h_hz = (float *)g_malloc(plan->xres*plan->yres*plan->zres*sizeof(float));

        if (set->plan.matmode==SV_MATMODE_FULL || set->plan.matmode==SV_MATMODE_ELECTRIC) {
            plan->h_epsilon = (float *)g_malloc(plan->xres*plan->yres*plan->zres*sizeof(float));
            plan->h_sigma =   (float *)g_malloc(plan->xres*plan->yres*plan->zres*sizeof(float));
        }
        if (set->plan.matmode==SV_MATMODE_FULL || set->plan.matmode==SV_MATMODE_MAGNETIC) {
            plan->h_mu =      (float *)g_malloc(plan->xres*plan->yres*plan->zres*sizeof(float));
            plan->h_sigast =  (float *)g_malloc(plan->xres*plan->yres*plan->zres*sizeof(float));
        }
        if (mp->nmat) {
            plan->h_mat = (int *)g_malloc(plan->xres*plan->yres*plan->zres*sizeof(int));
            plan->h_mattype =  (int *)g_malloc(mp->nmat*sizeof(int));
            plan->h_mattab =  (float *)g_malloc(SV_GM_N*mp->nmat*sizeof(float));
            plan->nmat = mp->nmat;
            for (i=0; i<mp->nmat; i++) {
                plan->h_mattype[i] = mp->mats[i].type;
                plan->h_mattab[SV_GM_N*i + SV_GM_EPSILON]       = (float)mp->mats[i].epsilon;
                plan->h_mattab[SV_GM_N*i + SV_GM_MU]            = (float)mp->mats[i].mu;
                plan->h_mattab[SV_GM_N*i + SV_GM_SIGMA]         = (float)mp->mats[i].sigma;
                plan->h_mattab[SV_GM_N*i + SV_GM_SIGAST]        = (float)mp->mats[i].sigast;
                plan->h_mattab[SV_GM_N*i + SV_GM_DRUDE_OMEGA_P] = (float)mp->mats[i].drude_omega_p;
                plan->h_mattab[SV_GM_N*i + SV_GM_DRUDE_NU]      = (float)mp->mats[i].drude_nu;
                plan->h_mattab[SV_GM_N*i + SV_GM_CP3_A0]        = (float)mp->mats[i].cp3_a[0];
                plan->h_mattab[SV_GM_N*i + SV_GM_CP3_A1]        = (float)mp->mats[i].cp3_a[1];
                plan->h_mattab[SV_GM_N*i + SV_GM_CP3_A2]        = (float)mp->mats[i].cp3_a[2];
                plan->h_mattab[SV_GM_N*i + SV_GM_CP3_PHI0]      = (float)mp->mats[i].cp3_phi[0];
                plan->h_mattab[SV_GM_N*i + SV_GM_CP3_PHI1]      = (float)mp->mats[i].cp3_phi[1];
                plan->h_mattab[SV_GM_N*i + SV_GM_CP3_PHI2]      = (float)mp->mats[i].cp3_phi[2];
                plan->h_mattab[SV_GM_N*i + SV_GM_CP3_OMEGA0]    = (float)mp->mats[i].cp3_omega[0];
                plan->h_mattab[SV_GM_N*i + SV_GM_CP3_OMEGA1]    = (float)mp->mats[i].cp3_omega[1];
                plan->h_mattab[SV_GM_N*i + SV_GM_CP3_OMEGA2]    = (float)mp->mats[i].cp3_omega[2];
                plan->h_mattab[SV_GM_N*i + SV_GM_CP3_GAMMA0]    = (float)mp->mats[i].cp3_gamma[0];
                plan->h_mattab[SV_GM_N*i + SV_GM_CP3_GAMMA1]    = (float)mp->mats[i].cp3_gamma[1];
                plan->h_mattab[SV_GM_N*i + SV_GM_CP3_GAMMA2]    = (float)mp->mats[i].cp3_gamma[2];
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_A0]        = (float)mp->mats[i].ade_a0;
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_A1]        = (float)mp->mats[i].ade_a1;
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_A2]        = (float)mp->mats[i].ade_a2;
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_BP0]       = (float)mp->mats[i].ade_bp0[0];
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_BP0+1]     = (float)mp->mats[i].ade_bp0[1];
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_BP1]       = (float)mp->mats[i].ade_bp1[0];
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_BP1+1]     = (float)mp->mats[i].ade_bp1[1];
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_BP2]       = (float)mp->mats[i].ade_bp2[0];
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_BP2+1]     = (float)mp->mats[i].ade_bp2[1];
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_BP3]       = (float)mp->mats[i].ade_bp3[0];
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_BP3+1]     = (float)mp->mats[i].ade_bp3[1];
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_BP4]       = (float)mp->mats[i].ade_bp4[0];
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_BP4+1]     = (float)mp->mats[i].ade_bp4[1];
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_C0]        = (float)mp->mats[i].ade_c0;
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_C1]        = (float)mp->mats[i].ade_c1;
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_C2]        = (float)mp->mats[i].ade_c2;
                plan->h_mattab[SV_GM_N*i + SV_GM_ADE_C3]        = (float)mp->mats[i].ade_c3;
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_D_CHI]    = (float)mp->mats[i].plrc_d_chi;
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_D_XI]     = (float)mp->mats[i].plrc_d_xi;
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_D_DCHI]   = (float)mp->mats[i].plrc_d_dchi;
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_D_DXI]    = (float)mp->mats[i].plrc_d_dxi;
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_P_CHI]    = (float)mp->mats[i].plrc_p_chi[0];
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_P_CHI+1]  = (float)mp->mats[i].plrc_p_chi[1];
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_P_XI]     = (float)mp->mats[i].plrc_p_xi[0];
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_P_XI+1]   = (float)mp->mats[i].plrc_p_xi[1];
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_P_DCHIR]  = (float)mp->mats[i].plrc_p_dchir[0];
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_P_DCHIR+1]= (float)mp->mats[i].plrc_p_dchir[1];
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_P_DXIR]   = (float)mp->mats[i].plrc_p_dxir[0];
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_P_DXIR+1] = (float)mp->mats[i].plrc_p_dxir[1];
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_P_DCHII]  = (float)mp->mats[i].plrc_p_dchii[0];
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_P_DCHII+1]= (float)mp->mats[i].plrc_p_dchii[1];
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_P_DXII]   = (float)mp->mats[i].plrc_p_dxii[0];
                plan->h_mattab[SV_GM_N*i + SV_GM_PLRC_P_DXII+1] = (float)mp->mats[i].plrc_p_dxii[1];
            }
        } else {
            plan->h_mat = NULL;
            plan->h_mattype = NULL;
            plan->h_mattab = NULL;
            plan->nmat = 0;
        }

	    plan->h_isplrc = plan->h_isade = FALSE;
	    for (i=0; i<mp->nmat; i++)
	    {
	    	if (mp->mats[i].type==SV_MAT_DRUDE || mp->mats[i].type==SV_MAT_CP3 || mp->mats[i].type==SV_MAT_CP || mp->mats[i].type==SV_MAT_PLRC) 
	    		plan->h_isplrc = TRUE;

	    	if (mp->mats[i].type==SV_MAT_ADE) 
	    		plan->h_isade = TRUE;

	    }


        plan->noutpoints = set->so.npnts;
        plan->h_outpoint_pos = (int *)g_malloc(3*plan->noutpoints*sizeof(int));
        for (i=0; i<plan->noutpoints; i++) {
            plan->h_outpoint_pos[3*i] = set->so.pnts[i].i;
            plan->h_outpoint_pos[3*i + 1] = set->so.pnts[i].j;
            plan->h_outpoint_pos[3*i + 2] = set->so.pnts[i].k;
        }
        plan->h_outpointdata = (float *)g_malloc(6*plan->noutpoints*set->sc.nsteps*sizeof(float));
        for (i=0; i<(6*plan->noutpoints*plan->nsteps); i++) plan->h_outpointdata[i] = 0;

        if (plan->bndx0==SV_BOUNDARY_CPML || plan->bndxn==SV_BOUNDARY_CPML || plan->bndy0==SV_BOUNDARY_CPML || plan->bndyn==SV_BOUNDARY_CPML
            || plan->bndz0==SV_BOUNDARY_CPML || plan->bndzn==SV_BOUNDARY_CPML)
        {
            plan->h_depths = (int *)g_malloc(6*sizeof(int));

            plan->h_depths[0] = plan->cpml_depth_bx0 = set->sb.depth_bx0;
            plan->h_depths[1] = plan->cpml_depth_bxn = set->sb.depth_bxn;
            plan->h_depths[2] = plan->cpml_depth_by0 = set->sb.depth_by0; 
            plan->h_depths[3] = plan->cpml_depth_byn = set->sb.depth_byn;
            plan->h_depths[4] = plan->cpml_depth_bz0 = set->sb.depth_bz0;
            plan->h_depths[5] = plan->cpml_depth_bzn = set->sb.depth_bzn;
        } 

        if (plan->bndx0 == SV_BOUNDARY_CPML) {            
            plan->h_cpml_be_x0 = (float *)g_malloc(plan->cpml_depth_bx0*sizeof(float));
            plan->h_cpml_ce_x0 = (float *)g_malloc(plan->cpml_depth_bx0*sizeof(float));
            plan->h_cpml_bh_x0 = (float *)g_malloc(plan->cpml_depth_bx0*sizeof(float));
            plan->h_cpml_ch_x0 = (float *)g_malloc(plan->cpml_depth_bx0*sizeof(float));
            plan->h_cpml_kappae_x0 = (float *)g_malloc(plan->cpml_depth_bx0*sizeof(float));
            plan->h_cpml_kappah_x0 = (float *)g_malloc(plan->cpml_depth_bx0*sizeof(float));

            plan->h_cpml_p_x0 = (float *)g_malloc(4*plan->yres*plan->zres*plan->cpml_depth_bx0*sizeof(float));

            for (i=0; i<(4*plan->yres*plan->zres*plan->cpml_depth_bx0); i++) 
                plan->h_cpml_p_x0[i] = 0;
            
            for (i=0; i<set->sb.depth_bx0; i++)
            {
                plan->h_cpml_be_x0[i] = (float)mp->bnd->cpml.be_x0[i];
                plan->h_cpml_ce_x0[i] = (float)mp->bnd->cpml.ce_x0[i];
                plan->h_cpml_bh_x0[i] = (float)mp->bnd->cpml.bh_x0[i];
                plan->h_cpml_ch_x0[i] = (float)mp->bnd->cpml.ch_x0[i];
                plan->h_cpml_kappae_x0[i] = (float)mp->bnd->cpml.kappae_x0[i];
                plan->h_cpml_kappah_x0[i] = (float)mp->bnd->cpml.kappah_x0[i];
            }
        }
        if (plan->bndxn == SV_BOUNDARY_CPML) {

            plan->h_cpml_be_xn = (float *)g_malloc(plan->cpml_depth_bxn*sizeof(float));
            plan->h_cpml_ce_xn = (float *)g_malloc(plan->cpml_depth_bxn*sizeof(float));
            plan->h_cpml_bh_xn = (float *)g_malloc(plan->cpml_depth_bxn*sizeof(float));
            plan->h_cpml_ch_xn = (float *)g_malloc(plan->cpml_depth_bxn*sizeof(float));
            plan->h_cpml_kappae_xn = (float *)g_malloc(plan->cpml_depth_bxn*sizeof(float));
            plan->h_cpml_kappah_xn = (float *)g_malloc(plan->cpml_depth_bxn*sizeof(float));

            plan->h_cpml_p_xn = (float *)g_malloc(4*plan->yres*plan->zres*plan->cpml_depth_bxn*sizeof(float));

            for (i=0; i<(4*plan->yres*plan->zres*plan->cpml_depth_bxn); i++) 
                plan->h_cpml_p_xn[i] = 0;


            for (i=0; i<set->sb.depth_bxn; i++) {
                plan->h_cpml_be_xn[i] = (float)mp->bnd->cpml.be_xn[i];
                plan->h_cpml_ce_xn[i] = (float)mp->bnd->cpml.ce_xn[i];
                plan->h_cpml_bh_xn[i] = (float)mp->bnd->cpml.bh_xn[i];
                plan->h_cpml_ch_xn[i] = (float)mp->bnd->cpml.ch_xn[i];
                plan->h_cpml_kappae_xn[i] = (float)mp->bnd->cpml.kappae_xn[i];
                plan->h_cpml_kappah_xn[i] = (float)mp->bnd->cpml.kappah_xn[i];
            }
        }

        if (plan->bndy0 == SV_BOUNDARY_CPML) {

            plan->h_cpml_be_y0 = (float *)g_malloc(plan->cpml_depth_by0*sizeof(float));
            plan->h_cpml_ce_y0 = (float *)g_malloc(plan->cpml_depth_by0*sizeof(float));
            plan->h_cpml_bh_y0 = (float *)g_malloc(plan->cpml_depth_by0*sizeof(float));
            plan->h_cpml_ch_y0 = (float *)g_malloc(plan->cpml_depth_by0*sizeof(float));
            plan->h_cpml_kappae_y0 = (float *)g_malloc(plan->cpml_depth_by0*sizeof(float));
            plan->h_cpml_kappah_y0 = (float *)g_malloc(plan->cpml_depth_by0*sizeof(float));

            plan->h_cpml_p_y0 = (float *)g_malloc(4*plan->xres*plan->zres*plan->cpml_depth_by0*sizeof(float));

            for (i=0; i<(4*plan->xres*plan->zres*plan->cpml_depth_by0); i++) 
                plan->h_cpml_p_y0[i] = 0;

            for (i=0; i<set->sb.depth_by0; i++) {
                plan->h_cpml_be_y0[i] = (float)mp->bnd->cpml.be_y0[i];
                plan->h_cpml_ce_y0[i] = (float)mp->bnd->cpml.ce_y0[i];
                plan->h_cpml_bh_y0[i] = (float)mp->bnd->cpml.bh_y0[i];
                plan->h_cpml_ch_y0[i] = (float)mp->bnd->cpml.ch_y0[i];
                plan->h_cpml_kappae_y0[i] = (float)mp->bnd->cpml.kappae_y0[i];
                plan->h_cpml_kappah_y0[i] = (float)mp->bnd->cpml.kappah_y0[i];
            }
        }

        if (plan->bndyn == SV_BOUNDARY_CPML) {
            plan->h_cpml_be_yn = (float *)g_malloc(plan->cpml_depth_byn*sizeof(float));
            plan->h_cpml_ce_yn = (float *)g_malloc(plan->cpml_depth_byn*sizeof(float));
            plan->h_cpml_bh_yn = (float *)g_malloc(plan->cpml_depth_byn*sizeof(float));
            plan->h_cpml_ch_yn = (float *)g_malloc(plan->cpml_depth_byn*sizeof(float));
            plan->h_cpml_kappae_yn = (float *)g_malloc(plan->cpml_depth_byn*sizeof(float));
            plan->h_cpml_kappah_yn = (float *)g_malloc(plan->cpml_depth_byn*sizeof(float));

            plan->h_cpml_p_yn = (float *)g_malloc(4*plan->xres*plan->zres*plan->cpml_depth_byn*sizeof(float));

            for (i=0; i<(4*plan->xres*plan->zres*plan->cpml_depth_byn); i++) 
                plan->h_cpml_p_yn[i] = 0;

            for (i=0; i<set->sb.depth_byn; i++) {
                plan->h_cpml_be_yn[i] = (float)mp->bnd->cpml.be_yn[i];
                plan->h_cpml_ce_yn[i] = (float)mp->bnd->cpml.ce_yn[i];
                plan->h_cpml_bh_yn[i] = (float)mp->bnd->cpml.bh_yn[i];
                plan->h_cpml_ch_yn[i] = (float)mp->bnd->cpml.ch_yn[i];
                plan->h_cpml_kappae_yn[i] = (float)mp->bnd->cpml.kappae_yn[i];
                plan->h_cpml_kappah_yn[i] = (float)mp->bnd->cpml.kappah_yn[i];
            }
        }

        if (plan->bndz0 == SV_BOUNDARY_CPML) {
            plan->h_cpml_be_z0 = (float *)g_malloc(plan->cpml_depth_bz0*sizeof(float));
            plan->h_cpml_ce_z0 = (float *)g_malloc(plan->cpml_depth_bz0*sizeof(float));
            plan->h_cpml_bh_z0 = (float *)g_malloc(plan->cpml_depth_bz0*sizeof(float));
            plan->h_cpml_ch_z0 = (float *)g_malloc(plan->cpml_depth_bz0*sizeof(float));
            plan->h_cpml_kappae_z0 = (float *)g_malloc(plan->cpml_depth_bz0*sizeof(float));
            plan->h_cpml_kappah_z0 = (float *)g_malloc(plan->cpml_depth_bz0*sizeof(float));

            plan->h_cpml_p_z0 = (float *)g_malloc(4*plan->xres*plan->yres*plan->cpml_depth_bz0*sizeof(float));

            for (i=0; i<(4*plan->xres*plan->yres*plan->cpml_depth_bz0); i++) 
                plan->h_cpml_p_z0[i] = 0;

            for (i=0; i<set->sb.depth_bz0; i++) {
                plan->h_cpml_be_z0[i] = (float)mp->bnd->cpml.be_z0[i];
                plan->h_cpml_ce_z0[i] = (float)mp->bnd->cpml.ce_z0[i];
                plan->h_cpml_bh_z0[i] = (float)mp->bnd->cpml.bh_z0[i];
                plan->h_cpml_ch_z0[i] = (float)mp->bnd->cpml.ch_z0[i];
                plan->h_cpml_kappae_z0[i] = (float)mp->bnd->cpml.kappae_z0[i];
                plan->h_cpml_kappah_z0[i] = (float)mp->bnd->cpml.kappah_z0[i];
            }
        }

        if (plan->bndzn == SV_BOUNDARY_CPML) {
            plan->h_cpml_be_zn = (float *)g_malloc(plan->cpml_depth_bzn*sizeof(float));
            plan->h_cpml_ce_zn = (float *)g_malloc(plan->cpml_depth_bzn*sizeof(float));
            plan->h_cpml_bh_zn = (float *)g_malloc(plan->cpml_depth_bzn*sizeof(float));
            plan->h_cpml_ch_zn = (float *)g_malloc(plan->cpml_depth_bzn*sizeof(float));
            plan->h_cpml_kappae_zn = (float *)g_malloc(plan->cpml_depth_bzn*sizeof(float));
            plan->h_cpml_kappah_zn = (float *)g_malloc(plan->cpml_depth_bzn*sizeof(float));

            plan->h_cpml_p_zn = (float *)g_malloc(4*plan->xres*plan->yres*plan->cpml_depth_bzn*sizeof(float));

            for (i=0; i<(4*plan->xres*plan->yres*plan->cpml_depth_bzn); i++) 
                plan->h_cpml_p_zn[i] = 0;

            for (i=0; i<set->sb.depth_bzn; i++) {
                plan->h_cpml_be_zn[i] = (float)mp->bnd->cpml.be_zn[i];
                plan->h_cpml_ce_zn[i] = (float)mp->bnd->cpml.ce_zn[i];
                plan->h_cpml_bh_zn[i] = (float)mp->bnd->cpml.bh_zn[i];
                plan->h_cpml_ch_zn[i] = (float)mp->bnd->cpml.ch_zn[i];
                plan->h_cpml_kappae_zn[i] = (float)mp->bnd->cpml.kappae_zn[i];
                plan->h_cpml_kappah_zn[i] = (float)mp->bnd->cpml.kappah_zn[i];
            }
        }


        if (mp->src->tsf) {
            plan->h_tsfset = (float *)g_malloc(29*sizeof(float));
            plan->h_tsf_jpool_e = (float *)g_malloc(mp->src->tsf->ndata*sizeof(float));
            plan->h_tsf_jpool_h = (float *)g_malloc(mp->src->tsf->ndata*sizeof(float));
            plan->h_tsf_jpool_epsilon = (float *)g_malloc(mp->src->tsf->ndata*sizeof(float));
            plan->h_tsf_jpool_sigma = (float *)g_malloc(mp->src->tsf->ndata*sizeof(float));
            plan->h_tsf_jpool_mu = (float *)g_malloc(mp->src->tsf->ndata*sizeof(float));
            plan->h_tsf_jpool_sigast = (float *)g_malloc(mp->src->tsf->ndata*sizeof(float));
            plan->h_tsf_jpvals = (float *)g_malloc(8*sizeof(float));
            plan->tsf_jpool_size = mp->src->tsf->ndata;
            for (i=0; i<plan->tsf_jpool_size; i++) {
                plan->h_tsf_jpool_e[i] = 0;
                plan->h_tsf_jpool_h[i] = 0;
                plan->h_tsf_jpool_sigma[i] = 0;   /* natvrdo !! */
                plan->h_tsf_jpool_sigast[i] = 0;
                plan->h_tsf_jpool_mu[i] = (float)mp->src->tsf->layered_mu;
                plan->h_tsf_jpool_epsilon[i] = (float)mp->src->tsf->layered_epsilon;
            }
            for (i=0; i<8; i++) 
                plan->h_tsf_jpvals[i] = 0;
            plan->h_tsfset[0]  = (float)mp->src->tsf->box_i0;
            plan->h_tsfset[1]  = (float)mp->src->tsf->box_j0;
            plan->h_tsfset[2]  = (float)mp->src->tsf->box_k0;
            plan->h_tsfset[3]  = (float)mp->src->tsf->box_in;
            plan->h_tsfset[4]  = (float)mp->src->tsf->box_jn;
            plan->h_tsfset[5]  = (float)mp->src->tsf->box_kn;
            plan->h_tsfset[6]  = (float)mp->src->tsf->ia_theta;
            plan->h_tsfset[7]  = (float)mp->src->tsf->ia_phi;
            plan->h_tsfset[8]  = (float)mp->src->tsf->ia_psi;
            plan->h_tsfset[9]  = (float)mp->src->tsf->corr;
            plan->h_tsfset[10] = (float)mp->src->tsf->box_boundary_skipi0;
            plan->h_tsfset[11] = (float)mp->src->tsf->box_boundary_skipin;
            plan->h_tsfset[12] = (float)mp->src->tsf->box_boundary_skipj0;
            plan->h_tsfset[13] = (float)mp->src->tsf->box_boundary_skipjn;
            plan->h_tsfset[14] = (float)mp->src->tsf->box_boundary_skipk0;
            plan->h_tsfset[15] = (float)mp->src->tsf->box_boundary_skipkn; 
            plan->h_tsfset[16] = (float)set->ss.tsf.gaussian;    
            plan->h_tsfset[17] = (float)set->ss.tsf.gaussian_fxpos;
            plan->h_tsfset[18] = (float)set->ss.tsf.gaussian_fypos;
            plan->h_tsfset[19] = (float)set->ss.tsf.gaussian_rx;
            plan->h_tsfset[20] = (float)set->ss.tsf.gaussian_ry;
            // FIXME: implement RADIAL MULTIPLIER (tsf.radial parameters)
            plan->h_tsfset[21] = (float)set->ss.tsf.fiber;
            plan->h_tsfset[22] = (float)set->ss.tsf.fiber_fxpos;
            plan->h_tsfset[23] = (float)set->ss.tsf.fiber_fypos;
            plan->h_tsfset[24] = (float)set->ss.tsf.fiber_radius;
            plan->h_tsfset[25] = (float)set->ss.tsf.fiber_cutoff;
            plan->h_tsfset[26] = (float)set->ss.tsf.fiber_epsilon_core;
            plan->h_tsfset[27] = (float)set->ss.tsf.fiber_epsilon_cladding;
            plan->h_tsfset[28] = (float)set->ss.lambda_center;
        } else 
            plan->h_tsfset = NULL;

        if (mp->src->tsff) {
            plan->h_tsffset = (float *)g_malloc(16*sizeof(float));

            plan->tsff_jpool_size = mp->src->tsff->ndata;
            plan->tsff_nxm = mp->src->tsff->focused_mip*mp->src->tsff->focused_nip;
            plan->tsff_freal = (float)(mp->src->tsff->focused_fdist*set->sp.dx);

            plan->h_tsff_jpool_e = (float *)g_malloc(plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float));
            plan->h_tsff_jpool_h = (float *)g_malloc(plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float));
            plan->h_tsff_jpool_epsilon = (float *)g_malloc(plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float));
            plan->h_tsff_jpool_sigma = (float *)g_malloc(plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float));
            plan->h_tsff_jpool_mu = (float *)g_malloc(plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float));
            plan->h_tsff_jpool_sigast = (float *)g_malloc(plan->tsff_nxm*plan->tsff_jpool_size*sizeof(float));
            plan->h_tsff_jpvals = (float *)g_malloc(plan->tsff_nxm*8*sizeof(float));
            plan->tsff_data = (float *)g_malloc(plan->tsff_nxm*6*sizeof(float));

            for (i=0; i<(plan->tsff_nxm*plan->tsff_jpool_size); i++) {
                plan->h_tsff_jpool_e[i] = 0;
                plan->h_tsff_jpool_h[i] = 0;
                plan->h_tsff_jpool_sigma[i] = 0;  
                plan->h_tsff_jpool_sigast[i] = 0;
                plan->h_tsff_jpool_mu[i] = (float)mp->src->tsff->layered_mu;
                plan->h_tsff_jpool_epsilon[i] = (float)mp->src->tsff->layered_epsilon;
            }
            for (i=0; i<(plan->tsff_nxm*8); i++) 
                plan->h_tsff_jpvals[i] = 0;
            plan->h_tsffset[0] = (float)mp->src->tsff->box_i0;
            plan->h_tsffset[1] = (float)mp->src->tsff->box_j0;
            plan->h_tsffset[2] = (float)mp->src->tsff->box_k0;
            plan->h_tsffset[3] = (float)mp->src->tsff->box_in;
            plan->h_tsffset[4] = (float)mp->src->tsff->box_jn;
            plan->h_tsffset[5] = (float)mp->src->tsff->box_kn;
            plan->h_tsffset[6] = 0.0f;
            plan->h_tsffset[7] = 0.0f;
            plan->h_tsffset[8] = 0.0f;
            plan->h_tsffset[9] = 0.0f;
            plan->h_tsffset[10] = (float)mp->src->tsff->box_boundary_skipi0;
            plan->h_tsffset[11] = (float)mp->src->tsff->box_boundary_skipin;
            plan->h_tsffset[12] = (float)mp->src->tsff->box_boundary_skipj0;
            plan->h_tsffset[13] = (float)mp->src->tsff->box_boundary_skipjn;
            plan->h_tsffset[14] = (float)mp->src->tsff->box_boundary_skipk0;
            plan->h_tsffset[15] = (float)mp->src->tsff->box_boundary_skipkn;
            k = 0;
            for (i=0; i<mp->src->tsff->focused_mip; i++) {
                for (j=0; j<mp->src->tsff->focused_nip; j++) {
                    plan->tsff_data[k + 0] = (float)mp->src->tsff->ia_theta[i][j];
                    plan->tsff_data[k + 1] = (float)mp->src->tsff->ia_phi[i][j];
                    plan->tsff_data[k + 2] = (float)mp->src->tsff->ia_psi[i][j];
                    plan->tsff_data[k + 3] = (float)mp->src->tsff->corr[i][j];
                    plan->tsff_data[k + 4] = (float)mp->src->tsff->an[i][j];
                    plan->tsff_data[k + 5] = (float)mp->src->tsff->bm[i][j];
                    k+=6;
                }
            }

        } else 
            plan->h_tsffset = NULL;


        if (mp->src->sf) {
            plan->h_sfset = (float *)g_malloc(10*sizeof(float));
            plan->h_sf_jpool_e = (float *)g_malloc(mp->src->sf->ndata*sizeof(float));
            plan->h_sf_jpool_h = (float *)g_malloc(mp->src->sf->ndata*sizeof(float));
            plan->h_sf_jpool_epsilon = (float *)g_malloc(mp->src->sf->ndata*sizeof(float));
            plan->h_sf_jpool_sigma = (float *)g_malloc(mp->src->sf->ndata*sizeof(float));
            plan->h_sf_jpool_mu = (float *)g_malloc(mp->src->sf->ndata*sizeof(float));
            plan->h_sf_jpool_sigast = (float *)g_malloc(mp->src->sf->ndata*sizeof(float));
            plan->h_sf_jpvals = (float *)g_malloc(8*sizeof(float));
            plan->sf_jpool_size = mp->src->sf->ndata;
            for (i=0; i<plan->sf_jpool_size; i++) {
                plan->h_sf_jpool_e[i] = 0;
                plan->h_sf_jpool_h[i] = 0;
                plan->h_sf_jpool_sigma[i] = 0;   /* natvrdo !! */
                plan->h_sf_jpool_sigast[i] = 0;
                plan->h_sf_jpool_mu[i] = (float)mp->src->sf->layered_mu;
                plan->h_sf_jpool_epsilon[i] = (float)mp->src->sf->layered_epsilon;
            }
            for (i=0; i<8; i++) 
                plan->h_sf_jpvals[i] = 0;
            /*first five values are reserved for future use*/
            plan->h_sfset[6] = (float)mp->src->sf->ia_theta;
            plan->h_sfset[7] = (float)mp->src->sf->ia_phi;
            plan->h_sfset[8] = (float)mp->src->sf->ia_psi;

        } else 
            plan->h_sfset = NULL;


        if (set->sf.nrs) {
	        plan->ff_size = 0;
            splitsize = MAX(MAX(set->sf.box_in-set->sf.box_i0, set->sf.box_jn-set->sf.box_j0), set->sf.box_kn-set->sf.box_k0);
            
	    if (set->sf.nrs<=256) 
            plan->nhlps = (MIN(splitsize, plan->maxthreads)/set->sf.nrs);
        else 
            plan->nhlps = 0;

        if (set->sc.verbose > 1) 
            printf("We will use %d helpers for NFFF computation.\n", plan->nhlps);


	    plan->iset_size = 37 + 5*set->sf.nrs;
	    plan->h_iset = (int *)g_malloc(plan->iset_size*sizeof(int));

        plan->h_iset[0] = set->sf.box_i0;
        plan->h_iset[1] = set->sf.box_in;
        plan->h_iset[2] = set->sf.box_j0;
        plan->h_iset[3] = set->sf.box_jn;
        plan->h_iset[4] = set->sf.box_k0;
        plan->h_iset[5] = set->sf.box_kn;
        plan->h_iset[6] = set->sf.box_boundary_skipi0;
        plan->h_iset[7] = set->sf.skipi0_jmin;
        plan->h_iset[8] = set->sf.skipi0_kmin;
        plan->h_iset[9] = set->sf.skipi0_jmax;
        plan->h_iset[10] = set->sf.skipi0_kmax;
        plan->h_iset[11] = set->sf.box_boundary_skipin;
        plan->h_iset[12] = set->sf.skipin_jmin;
        plan->h_iset[13] = set->sf.skipin_kmin;
        plan->h_iset[14] = set->sf.skipin_jmax;
        plan->h_iset[15] = set->sf.skipin_kmax;
        plan->h_iset[16] = set->sf.box_boundary_skipj0;
        plan->h_iset[17] = set->sf.skipj0_imin;
        plan->h_iset[18] = set->sf.skipj0_kmin;
        plan->h_iset[19] = set->sf.skipj0_imax;
        plan->h_iset[20] = set->sf.skipj0_kmax;
        plan->h_iset[21] = set->sf.box_boundary_skipjn;
        plan->h_iset[22] = set->sf.skipjn_imin;
        plan->h_iset[23] = set->sf.skipjn_kmin;
        plan->h_iset[24] = set->sf.skipjn_imax;
        plan->h_iset[25] = set->sf.skipjn_kmax;
        plan->h_iset[26] = set->sf.box_boundary_skipk0;
        plan->h_iset[27] = set->sf.skipk0_imin;
        plan->h_iset[28] = set->sf.skipk0_jmin;
        plan->h_iset[29] = set->sf.skipk0_imax;
        plan->h_iset[30] = set->sf.skipk0_jmax;
        plan->h_iset[31] = set->sf.box_boundary_skipkn;
        plan->h_iset[32] = set->sf.skipkn_imin;
        plan->h_iset[33] = set->sf.skipkn_jmin;
        plan->h_iset[34] = set->sf.skipkn_imax;
        plan->h_iset[35] = set->sf.skipkn_jmax;
        plan->h_iset[36] = set->sf.nrs; 
	    for (i = 0; i<set->sf.nrs; i++) {
		    plan->h_iset[5*i+37] = mp->farfield->ndata;
		    plan->ff_size += mp->farfield->ndata;
		    plan->h_iset[5*i+38] = (int)mp->farfield->rpoints[i].istart; //istart
		    plan->h_iset[5*i+39] = (int)set->sf.ri[i];
		    plan->h_iset[5*i+40] = (int)set->sf.rj[i];
		    plan->h_iset[5*i+41] = (int)set->sf.rk[i];
	    }
	    plan->h_ff_ex = (float *)g_malloc(plan->ff_size*sizeof(float));
	    plan->h_ff_ey = (float *)g_malloc(plan->ff_size*sizeof(float));
	    plan->h_ff_ez = (float *)g_malloc(plan->ff_size*sizeof(float));
        for (i=0; i<plan->ff_size; i++) {
                plan->h_ff_ex[i] = plan->h_ff_ey[i] = plan->h_ff_ez[i] = 0;
        }
            
	    plan->h_ff_hlp_ex = (float *)g_malloc(plan->ff_size*plan->nhlps*sizeof(float));
	    plan->h_ff_hlp_ey = (float *)g_malloc(plan->ff_size*plan->nhlps*sizeof(float));
	    plan->h_ff_hlp_ez = (float *)g_malloc(plan->ff_size*plan->nhlps*sizeof(float));
        for (i=0; i<(plan->ff_size*plan->nhlps); i++) {
                plan->h_ff_hlp_ex[i] = plan->h_ff_hlp_ey[i] = plan->h_ff_hlp_ez[i] = 0;
        }
	} else plan->ff_size = 0;  
        if (set->spf.nrs) {
	        plan->pff_size = 0;
	        plan->piset_size = 21 + 5*set->spf.nrs;
	        plan->h_piset = (int *)g_malloc(plan->piset_size*sizeof(int));

            plan->h_piset[0] = set->spf.box_i0;
            plan->h_piset[1] = set->spf.box_in;
            plan->h_piset[2] = set->spf.box_j0;
            plan->h_piset[3] = set->spf.box_jn;
            plan->h_piset[4] = set->spf.box_k0;
            plan->h_piset[5] = set->spf.box_kn;
            plan->h_piset[6] = set->spf.pimin;
            plan->h_piset[7] = set->spf.pjmin;
            plan->h_piset[8] = set->spf.pimax;
            plan->h_piset[9] = set->spf.pjmax;
            plan->h_piset[10] = set->spf.nrs; //
            plan->h_piset[11] = set->spf.box_boundary_skipk0;
            plan->h_piset[12] = set->spf.skipk0_imin;
            plan->h_piset[13] = set->spf.skipk0_jmin;
            plan->h_piset[14] = set->spf.skipk0_imax;
            plan->h_piset[15] = set->spf.skipk0_jmax;
            plan->h_piset[16] = set->spf.box_boundary_skipkn;
            plan->h_piset[17] = set->spf.skipkn_imin;
            plan->h_piset[18] = set->spf.skipkn_jmin;
            plan->h_piset[19] = set->spf.skipkn_imax;
            plan->h_piset[20] = set->spf.skipkn_jmax;
            for (i = 0; i<set->spf.nrs; i++) {
		        plan->h_piset[5*i+21] = mp->pfarfield->ndata;
		        plan->pff_size += mp->pfarfield->ndata;
		        plan->h_piset[5*i+22] = (int)mp->pfarfield->prpoints[i].istart; //istart
		        plan->h_piset[5*i+23] = (int)set->spf.ri[i];
		        plan->h_piset[5*i+24] = (int)set->spf.rj[i];
		        plan->h_piset[5*i+25] = (int)set->spf.rk[i];
	        }
            printf("pisize: %d %d = %d\n", set->spf.nrs, mp->pfarfield->ndata, plan->pff_size);
	        plan->h_pff_ex = (float *)g_malloc(plan->pff_size*sizeof(float));
	        plan->h_pff_ey = (float *)g_malloc(plan->pff_size*sizeof(float));
	        plan->h_pff_ez = (float *)g_malloc(plan->pff_size*sizeof(float));
            for (i=0; i<plan->pff_size; i++) {
                plan->h_pff_ex[i] = plan->h_pff_ey[i] = plan->h_pff_ez[i] = 0;
            }
	    if (plan->h_piset[10]<GPU_PERIODIC_HLPS_LIMIT) {
		    plan->h_pff_hlp_ex = (float *)g_malloc(plan->pff_size*(plan->h_piset[8]-plan->h_piset[6])*(plan->h_piset[9]-plan->h_piset[7])*sizeof(float));
		    plan->h_pff_hlp_ey = (float *)g_malloc(plan->pff_size*(plan->h_piset[8]-plan->h_piset[6])*(plan->h_piset[9]-plan->h_piset[7])*sizeof(float));
		    plan->h_pff_hlp_ez = (float *)g_malloc(plan->pff_size*(plan->h_piset[8]-plan->h_piset[6])*(plan->h_piset[9]-plan->h_piset[7])*sizeof(float));
		    for (i=0; i<(plan->pff_size*(plan->h_piset[8]-plan->h_piset[6])*(plan->h_piset[9]-plan->h_piset[7])); i++) {
			    plan->h_pff_hlp_ex[i] = plan->h_pff_hlp_ey[i] = plan->h_pff_hlp_ez[i] = 0;
		    }
	    }
	} else plan->pff_size = 0;  
        if (set->so.nsums) {
            plan->nsteps = set->sc.nsteps;
            plan->nsums = set->so.nsums;
            plan->h_sums = (float *)g_malloc(set->so.nsums*set->sc.nsteps*sizeof(float));
            plan->h_sum_epsilon = (float *)g_malloc(set->so.nsums*sizeof(float));
            plan->h_sum_sigma = (float *)g_malloc(set->so.nsums*sizeof(float));
            plan->h_sum_mu = (float *)g_malloc(set->so.nsums*sizeof(float));
            plan->h_sum_sigast = (float *)g_malloc(set->so.nsums*sizeof(float));
            plan->h_sum_mode = (int *)g_malloc(set->so.nsums*sizeof(int));
            plan->h_sum_i0 = (int *)g_malloc(set->so.nsums*sizeof(int));
            plan->h_sum_i1 = (int *)g_malloc(set->so.nsums*sizeof(int));
            plan->h_sum_j0 = (int *)g_malloc(set->so.nsums*sizeof(int));
            plan->h_sum_j1 = (int *)g_malloc(set->so.nsums*sizeof(int));
            plan->h_sum_k0 = (int *)g_malloc(set->so.nsums*sizeof(int));
            plan->h_sum_k1 = (int *)g_malloc(set->so.nsums*sizeof(int));
                  
            for (i=0; i<(set->so.nsums*set->sc.nsteps); i++)
                plan->h_sums[i] = 0;

            for (i=0; i<(set->so.nsums); i++) {
                plan->h_sum_epsilon[i] = (float)set->so.sums[i].layered_epsilon;
                printf("sum eps %g\n", plan->h_sum_epsilon[i]);
                plan->h_sum_mu[i] = (float)set->so.sums[i].layered_mu;
                plan->h_sum_sigma[i] = (float)set->so.sums[i].layered_sigma;
                plan->h_sum_sigast[i] = (float)set->so.sums[i].layered_sigast;
                plan->h_sum_mode[i] = set->so.sums[i].component;
                plan->h_sum_i0[i] = set->so.sums[i].box_i0;
                plan->h_sum_i1[i] = set->so.sums[i].box_in;
                plan->h_sum_j0[i] = set->so.sums[i].box_j0;
                plan->h_sum_j1[i] = set->so.sums[i].box_jn;
                plan->h_sum_k0[i] = set->so.sums[i].box_k0;
                plan->h_sum_k1[i] = set->so.sums[i].box_kn;
            }
        } else 
            plan->nsums = 0;

        if (set->so.nforces) {
            plan->nsteps = set->sc.nsteps;
            plan->nforces = set->so.nforces;
            plan->h_forces = (float *)g_malloc(3*set->so.nforces*set->sc.nsteps*sizeof(float));
            plan->h_force_i0 = (int *)g_malloc(set->so.nforces*sizeof(int));
            plan->h_force_i1 = (int *)g_malloc(set->so.nforces*sizeof(int));
            plan->h_force_j0 = (int *)g_malloc(set->so.nforces*sizeof(int));
            plan->h_force_j1 = (int *)g_malloc(set->so.nforces*sizeof(int));
            plan->h_force_k0 = (int *)g_malloc(set->so.nforces*sizeof(int));
            plan->h_force_k1 = (int *)g_malloc(set->so.nforces*sizeof(int));
                  
            for (i=0; i<(3*set->so.nforces*set->sc.nsteps); i++)
                plan->h_forces[i] = 0;

            for (i=0; i<(set->so.nforces); i++) {
                plan->h_force_i0[i] = set->so.forces[i].box_i0;
                plan->h_force_i1[i] = set->so.forces[i].box_in;
                plan->h_force_j0[i] = set->so.forces[i].box_j0;
                plan->h_force_j1[i] = set->so.forces[i].box_jn;
                plan->h_force_k0[i] = set->so.forces[i].box_k0;
                plan->h_force_k1[i] = set->so.forces[i].box_kn;
            }
        } else 
            plan->nforces = 0;
             

        plan->mb_bx0 = set->smb.bx0;
        plan->mb_by0 = set->smb.by0;
        plan->mb_bz0 = set->smb.bz0;
        plan->mb_bxn = set->smb.bxn;
        plan->mb_byn = set->smb.byn;
        plan->mb_bzn = set->smb.bzn;
        plan->mb_bx0pos = set->smb.bx0pos;
        plan->mb_by0pos = set->smb.by0pos;
        plan->mb_bz0pos = set->smb.bz0pos;
        plan->mb_bxnpos = set->smb.bxnpos;
        plan->mb_bynpos = set->smb.bynpos;
        plan->mb_bznpos = set->smb.bznpos;


        /*start thread*/
        g_thread_create((GThreadFunc)gpu_main, (void *)(gpu->gset[n].aplan), TRUE, &gerr);

    }
    return gpu;
}

#endif

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
