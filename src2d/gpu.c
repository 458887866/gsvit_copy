
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
#define AT(i, j) ((j)*plan->xres + (i))

#define GPU_PERIODIC_HLPS_LIMIT 120


void sv_pool_sync(SvPool *mp, SvSet *set);

int sv_pool_gpu_init(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose>1) 
        printf("pool_gpu_init\n");
    for (n=0; n<gpu->ngset; n++)
       g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_INIT));
    sv_pool_sync(mp, set);

    return 0;
}

int sv_pool_gpu_ystep_h(SvPool *mp, SvSet *set)
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

    if (set->sc.verbose > 1) 
        printf("done.\n");

    return 0;
}

int sv_pool_gpu_ystep_e(SvPool *mp, SvSet *set)
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

    if (set->sc.verbose > 1) 
        printf("done.\n");

    return 0;

}

int sv_pool_gpu_liao_cpy(SvPool *mp, SvSet *set)
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
    if (set->sc.verbose > 1) 
        printf("done.\n");

    return 0;
}

int sv_pool_gpu_liao_bnd(SvPool *mp, SvSet *set)
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
    if (set->sc.verbose > 1) 
        printf("done.\n");

    return 0;
}

int sv_pool_gpu_getsum(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpuPlan *plan;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (!set->so.nsums) 
        return 0;

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

    if (set->sc.verbose > 1) 
        printf("done.\n");

    return 0;
}

int sv_pool_gpu_getpoints(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpuPlan *plan;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (!set->so.npnts) 
        return 0;

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

    if (set->sc.verbose > 1) 
        printf("done.\n");

    return 0;
}

int sv_pool_gpu_mbnd_e(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running internal boundary GPU step(s)...   ");
	    fflush(stdout);
    }

    if (set->smb.bx0 == 4 || set->smb.bxn == 4) {
        for (n=0; n<gpu->ngset; n++)
            g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_MBNDX_E));
        sv_pool_sync(mp, set);
    }
    if (set->smb.by0 == 4 || set->smb.byn == 4) {
        for (n=0; n<gpu->ngset; n++)
            g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_MBNDY_E));
        sv_pool_sync(mp, set);
    }
    if (set->sc.verbose > 1) 
        printf("done.\n");

    return 0;
}

int sv_pool_gpu_mbnd_h(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running internal boundary GPU step(s)...   ");
	    fflush(stdout);
    }

    if (set->smb.bx0 == 4 || set->smb.bxn == 4) {
        for (n=0; n<gpu->ngset; n++)
            g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_MBNDX_H));
        sv_pool_sync(mp, set);
    }
    if (set->smb.by0 == 4 || set->smb.byn == 4) {
        for (n=0; n<gpu->ngset; n++)
            g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_MBNDY_H));
        sv_pool_sync(mp, set);
    }
    if (set->sc.verbose > 1) 
        printf("done.\n");

    return 0;
}

int sv_pool_gpu_point_srce(SvPool *mp, SvSet *set)
{
    gint i, j, n, pos;
    SvGpuPlan *plan;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running source GPU step...          ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++) {
	    plan = gpu->gset[n].aplan->plan;
	    for (i=0; i<set->ss.npnts; i++) {
		    /*if we are lucky, data are not shifted by user and can be used directly*/
		    pos = -1;
		    if (set->sc.step_act <= mp->src->sp[i].sdata.ndata && mp->src->sp[i].sdata.pos[set->sc.step_act] == set->sc.step_act)
			    pos = set->sc.step_act;
		    else {  /*otherwise we need to search for right position*/
			    for (j=0; j<mp->src->sp[i].sdata.ndata; j++) {
				    if (mp->src->sp[i].sdata.pos[j] == set->sc.step_act) {
					    pos = j;
					    break;
				    }
			    }
		    }
		    if (pos>=0) {
			    if (mp->src->sp[i].i < plan->xfrom || mp->src->sp[i].i >= plan->xto) 
                    continue;

			    plan->source_i = mp->src->sp[i].i - plan->xfrom;
			    plan->source_j = mp->src->sp[i].j;
			    plan->source_ex = (float)mp->src->sp[i].sdata.ex[pos];
			    plan->source_ey = (float)mp->src->sp[i].sdata.ey[pos];
			    plan->source_ez = (float)mp->src->sp[i].sdata.ez[pos];
	
			    g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_SRCE_POINT));
			    sv_pool_sync(mp, set);
		    }
	    }
    }

    if (set->sc.verbose > 1) 
        printf("done.\n");

    return 0;
}

int sv_pool_gpu_tsf_estep(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running TSF source GPU estep...          ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++) {
        g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_TSF_ESTEP));
        sv_pool_sync(mp, set);
    }
    
    if (set->sc.verbose > 1) 
        printf("done.\n");

    return 0;
}

int sv_pool_gpu_apply_farfield(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpu* gpu = (SvGpu *)mp->gpu;
    SvGpuPlan *plan;

    if (set->sf.nrs == 0) 
        return 0;

    if (set->sc.verbose > 1) {
	    printf("Running farfield GPU estep...       ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++) {
        plan = gpu->gset[n].aplan->plan;
        if (plan[n].nhlps==0) 
            g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_FF));
        else 
            g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_FFA));
        sv_pool_sync(mp, set);
    }
    
    if (set->sc.verbose > 1) 
        printf("done.\n");

    return 0;
}


int sv_pool_gpu_tsf_hstep(SvPool *mp, SvSet *set)
{
    gint n;
    SvGpuPlan *plan;
    SvGpu* gpu = (SvGpu *)mp->gpu;

    if (set->sc.verbose > 1) {
	    printf("Running TSF 1dpool step...          ");
	    fflush(stdout);
    }

    for (n=0; n<gpu->ngset; n++) {
	    plan = gpu->gset[n].aplan->plan;
        plan->source_ex = (float)mp->src->tsf->e[set->sc.step_act];
        g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_TSF_JSTEP));
        sv_pool_sync(mp, set);
    }
    if (set->sc.verbose > 1) 
        printf("done.\n");


    if (set->sc.verbose > 1) {
	    printf("Running TSF source GPU hstep...          ");
	    fflush(stdout);
    }

    
    for (n=0; n<gpu->ngset; n++) {
        g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_TSF_HSTEP));
        sv_pool_sync(mp, set);
    }
    
    if (set->sc.verbose > 1) 
        printf("done.\n");

    return 0;
}


int
sv_pool_gpu_boundary_copy(SvPool *mp, SvSet *set)
{
    //run this always to increase step, if there is no liao boundary it does nothing else.
    sv_pool_gpu_liao_cpy(mp, set);

    return 0;
}

int
sv_pool_gpu_boundary_hstep(SvPool *mp, SvSet *set)
{    
    if (set->smb.bx0 == 4 || set->smb.bxn == 4 || set->smb.by0 == 4 || set->smb.byn == 4) 
        sv_pool_gpu_mbnd_h(mp, set);    

    return 0;
}

int
sv_pool_gpu_boundary_estep(SvPool *mp, SvSet *set)
{    
    if (set->sb.bx0 == SV_BOUNDARY_LIAO || set->sb.bxn == SV_BOUNDARY_LIAO ||
        set->sb.by0 == SV_BOUNDARY_LIAO || set->sb.byn == SV_BOUNDARY_LIAO)
        sv_pool_gpu_liao_bnd(mp, set);

    if (set->smb.bx0 == 4 || set->smb.bxn == 4 
        || set->smb.by0 == 4 || set->smb.byn == 4) 
        sv_pool_gpu_mbnd_e(mp, set);

    return 0;
}




void
sv_pool_gpu_apply_source_estep(SvPool *mp, SvSet *set)
{
    sv_pool_gpu_point_srce(mp, set);
    if (mp->src->tsf) 
        sv_pool_gpu_tsf_estep(mp, set);
}

int
sv_pool_gpu_apply_source_hstep(SvPool *mp, SvSet *set)
{
    if (mp->src->tsf) 
        sv_pool_gpu_tsf_hstep(mp, set);

    return 0;
}



int sv_pool_gpu_copyto(SvPool *mp, SvSet *set, SvGpuSync what)
{
    gint n, i, j;
    SvGpu* gpu = (SvGpu *)mp->gpu;
    SvGpuPlan *plan;
    long int pos, ip;

    if (set->sc.verbose>1) {
        printf("Copying data to GPU...");
        fflush(stdout);
    }
    for (n=0; n<gpu->ngset; n++) {
	    plan = gpu->gset[n].aplan->plan;
	    for (i=0; i<plan->xres; i++)
		    for (j=0; j<plan->yres; j++) {
			    pos = AT(i, j);
			    ip = i+plan->xfrom;
			    if (what == SV_SYNC_ALL) {
				    if (set->plan.matmode==SV_MATMODE_FULL || set->plan.matmode==SV_MATMODE_ELECTRIC) {
					    plan->h_epsilon[pos] = (float)mp->epsilon->data[ip + j*plan->xres];
					    plan->h_sigma[pos] = (float)mp->sigma->data[ip + j*plan->xres];
				    }
				    if (set->plan.matmode==SV_MATMODE_FULL || set->plan.matmode==SV_MATMODE_MAGNETIC) {
					    plan->h_mu[pos] = (float)mp->mu->data[ip + j*plan->xres];
					    plan->h_sigast[pos] = (float)mp->sigast->data[ip + j*plan->xres];
				    }
			    }
			    if (what == SV_SYNC_ALL || what == SV_SYNC_EH || what == SV_SYNC_E) {
				    plan->h_ex[pos] = (float)mp->ex->data[ip + j*plan->xres];
				    plan->h_ey[pos] = (float)mp->ey->data[ip + j*plan->xres];
				    plan->h_ez[pos] = (float)mp->ez->data[ip + j*plan->xres];
			    }
			    if (what == SV_SYNC_ALL || what == SV_SYNC_EH || what == SV_SYNC_H) {
				    plan->h_hx[pos] = (float)mp->hx->data[ip + j*plan->xres];
				    plan->h_hy[pos] = (float)mp->hy->data[ip + j*plan->xres];
				    plan->h_hz[pos] = (float)mp->hz->data[ip + j*plan->xres];
			    }
			    if (mp->nmat) {
				    plan->h_mat[pos] = (int)mp->mat->data[ip + j*plan->xres];
			    }
		    }
	    g_queue_push_tail(gpu->gset[n].aplan->todos, GINT_TO_POINTER(SV_GT_COPYTO));
    }
    sv_pool_sync(mp, set);
    if (set->sc.verbose > 1) 
        printf("done.\n");

    return 0;
}

int sv_pool_gpu_copyfrom(SvPool *mp, SvSet *set, SvGpuSync what)
{
    gint n, i, j;
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

    for (n=0; n<gpu->ngset; n++) {
	    plan = gpu->gset[n].aplan->plan;
	    for (i=0; i<plan->xres; i++)
		    for (j=0; j<plan->yres; j++) {
				    pos = AT(i, j);
				    ip = i+plan->xfrom;
				    if (what == SV_SYNC_ALL || what == SV_SYNC_EH || what == SV_SYNC_E) {
					    mp->ex->data[ip + j*plan->xres] = plan->h_ex[pos];
					    mp->ey->data[ip + j*plan->xres] = plan->h_ey[pos];
					    mp->ez->data[ip + j*plan->xres] = plan->h_ez[pos];
				    }
				    if (what == SV_SYNC_ALL || what == SV_SYNC_EH || what == SV_SYNC_H) {
					    mp->hx->data[ip + j*plan->xres] = plan->h_hx[pos];
					    mp->hy->data[ip + j*plan->xres] = plan->h_hy[pos];
					    mp->hz->data[ip + j*plan->xres] = plan->h_hz[pos];
				    }
			    }

	    if (set->sf.nrs) {
/*		    j = 0;
		    for (k=0; k<set->sf.nrs; k++) {
			    for (i=0; i<mp->farfield->ndata; i++)
			    {
				   mp->farfield->rpoints[k].ex->data[i] = plan->h_ff_ex[j];
                                   mp->farfield->rpoints[k].ey->data[i] = plan->h_ff_ey[j];
                                   mp->farfield->rpoints[k].ez->data[i] = plan->h_ff_ez[j];
                                   j++;
			    }
		    }
*/
	    }
	    for (n=0; n<set->so.nsums; n++) {
		    for (i=0; i<set->sc.nsteps; i++) {
			    mp->out->outsumdata[n][i] = plan->h_sums[set->sc.nsteps*n + i];
		    }
	    
        }

	    for (n=0; n<set->so.npnts; n++) {
		 //   for (i=0; i<(6*set->sc.nsteps); i++) {
		//	    mp->out->outpointdata[n][i] = plan->h_outpointdata[6*set->sc.nsteps*n + i];
		 //   }
	    }
    }

    if (set->sc.verbose > 1) 
        printf("done.\n");

    return 0;
 }

void sv_pool_sync(SvPool *mp, SvSet *set)
{
    gint n;
    gint all;
    SvGpu* gpu = (SvGpu *)mp->gpu;
    //if (set->sc.verbose>1) printf("pool_gpu_sync...\n");
    do {
    //   usleep(100);
        all = 0;
        for (n=0; n<gpu->ngset; n++) {
	        all += (int)g_queue_is_empty(gpu->gset[n].aplan->todos);
        }
    } while (all<gpu->ngset);
    //if (set->sc.verbose>1) printf("sync done.\n");
}

void gpu_init(SvGpuPlan *plan)
{
    //int dev;
    cudaError_t err;

    err = cudaSetDevice(plan->device);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
    //err = cudaGetDevice(&dev);
    //if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));


    err = cudaMalloc((void **)((void *)&(plan->d_ex)), plan->xres*plan->yres*sizeof(float));
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
    err = cudaMalloc((void **)((void *)&(plan->d_ey)), plan->xres*plan->yres*sizeof(float));
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
    err = cudaMalloc((void **)((void *)&(plan->d_ez)), plan->xres*plan->yres*sizeof(float));
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
    err = cudaMalloc((void **)((void *)&(plan->d_hx)), plan->xres*plan->yres*sizeof(float));
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
    err = cudaMalloc((void **)((void *)&(plan->d_hy)), plan->xres*plan->yres*sizeof(float));
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
    err = cudaMalloc((void **)((void *)&(plan->d_hz)), plan->xres*plan->yres*sizeof(float));
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));

    if (plan->matmode==SV_MATMODE_FULL || plan->matmode==SV_MATMODE_ELECTRIC) { 
        cudaMalloc((void **)((void *)&(plan->d_epsilon)), plan->xres*plan->yres*sizeof(float));
        cudaMalloc((void **)((void *)&(plan->d_sigma)), plan->xres*plan->yres*sizeof(float));
    } 
    if (plan->matmode==SV_MATMODE_FULL || plan->matmode==SV_MATMODE_MAGNETIC)  {
        cudaMalloc((void **)((void *)&(plan->d_mu)), plan->xres*plan->yres*sizeof(float));
        cudaMalloc((void **)((void *)&(plan->d_sigast)), plan->xres*plan->yres*sizeof(float));
    }
    if (plan->nmat) {
        cudaMalloc((void **)((void *)&(plan->d_mat)), plan->xres*plan->yres*sizeof(int));
        cudaMalloc((void **)((void *)&(plan->d_mattype)), plan->nmat*sizeof(int));
        cudaMalloc((void **)((void *)&(plan->d_mattab)), 18*plan->nmat*sizeof(float));
    }

    if (plan->h_isplrc) {
	    err = cudaMalloc((void **)((void *)&(plan->d_plrcx)), 7*plan->xres*plan->yres*sizeof(float));
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
	    err = cudaMalloc((void **)((void *)&(plan->d_plrcy)), 7*plan->xres*plan->yres*sizeof(float));
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
	    err = cudaMalloc((void **)((void *)&(plan->d_plrcz)), 7*plan->xres*plan->yres*sizeof(float));
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }


    
    if (plan->bndx0==SV_BOUNDARY_LIAO) {
        err = cudaMalloc((void **)((void *)&(plan->d_liao_x0)),12*plan->yres*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->bndxn==SV_BOUNDARY_LIAO) {
        err = cudaMalloc((void **)((void *)&(plan->d_liao_xn)),12*plan->yres*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->bndy0==SV_BOUNDARY_LIAO) {
        err = cudaMalloc((void **)((void *)&(plan->d_liao_y0)),12*plan->xres*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->bndyn==SV_BOUNDARY_LIAO) {
        err = cudaMalloc((void **)((void *)&(plan->d_liao_yn)),12*plan->xres*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }

    err = cudaMalloc((void **)((void *)&(plan->d_bnds)), 6*sizeof(int));
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));


    if (plan->h_tsfset) {
        err = cudaMalloc((void **)((void *)&(plan->d_tsf_jpool_e)),plan->tsf_jpool_size*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsf_jpool_h)),plan->tsf_jpool_size*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsf_jpool_epsilon)),plan->tsf_jpool_size*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsf_jpool_mu)),plan->tsf_jpool_size*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsf_jpool_sigma)),plan->tsf_jpool_size*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsf_jpool_sigast)),plan->tsf_jpool_size*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsfset)),21*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_tsf_jpvals)),8*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }
     if (plan->ff_size!=0) {
        err = cudaMalloc((void **)((void *)&(plan->d_iset)), plan->iset_size*sizeof(int));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_ff_ex)), plan->ff_size*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_ff_ey)), plan->ff_size*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_ff_ez)), plan->ff_size*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_ff_hlp_ex)), plan->ff_size*plan->nhlps*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_ff_hlp_ey)), plan->ff_size*plan->nhlps*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_ff_hlp_ez)), plan->ff_size*plan->nhlps*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }
    if (plan->nsums>0) {
        err = cudaMalloc((void **)((void *)&(plan->d_sums)),plan->nsums*plan->nsteps*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_epsilon)),plan->nsums*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_mu)),plan->nsums*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_sigma)),plan->nsums*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_sigast)),plan->nsums*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_mode)),plan->nsums*sizeof(int));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_i0)),plan->nsums*sizeof(int));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_i1)),plan->nsums*sizeof(int));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_j0)),plan->nsums*sizeof(int));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_j1)),plan->nsums*sizeof(int));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_k0)),plan->nsums*sizeof(int));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_sum_k1)),plan->nsums*sizeof(int));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));

        if (plan->dir==2) {
            err = cudaMalloc((void **)((void *)&(plan->d_sum_accumulator)),plan->xres*plan->yres*plan->nsums*sizeof(float));
            if (err != cudaSuccess) 
                printf("%s\n", cudaGetErrorString(err));

        }
    }
    if (plan->noutpoints>0) {
        err = cudaMalloc((void **)((void *)&(plan->d_outpointdata)),6*plan->noutpoints*plan->nsteps*sizeof(float));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMalloc((void **)((void *)&(plan->d_outpoint_pos)),3*plan->noutpoints*sizeof(int));
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }

}

void gpu_copyto(SvGpuPlan *plan, SvGpuSync what)
{
    
    cudaError_t err;

//    int dev;
//    err = cudaGetDevice(&dev);
//    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    //printf("thread reported that device %d is set\n", dev);


    if (what == SV_SYNC_ALL || what == SV_SYNC_EH || what == SV_SYNC_E) {
	    err = cudaMemcpy(plan->d_ex, plan->h_ex, plan->xres*plan->yres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_ey, plan->h_ey, plan->xres*plan->yres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_ez, plan->h_ez, plan->xres*plan->yres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }
    if (what == SV_SYNC_ALL || what == SV_SYNC_EH || what == SV_SYNC_H) {
	    err = cudaMemcpy(plan->d_hx, plan->h_hx, plan->xres*plan->yres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_hy, plan->h_hy, plan->xres*plan->yres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_hz, plan->h_hz, plan->xres*plan->yres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }
    if (what == SV_SYNC_ALL && (plan->matmode==SV_MATMODE_FULL || plan->matmode==SV_MATMODE_ELECTRIC))
    {
	    err = cudaMemcpy(plan->d_epsilon, plan->h_epsilon, plan->xres*plan->yres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_sigma, plan->h_sigma, plan->xres*plan->yres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    } 
    if (what == SV_SYNC_ALL && (plan->matmode==SV_MATMODE_FULL || plan->matmode==SV_MATMODE_MAGNETIC)) 
    {
	    err = cudaMemcpy(plan->d_mu, plan->h_mu, plan->xres*plan->yres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_sigast, plan->h_sigast, plan->xres*plan->yres*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }

    if (what == SV_SYNC_ALL && plan->noutpoints>0) {
  	    err = cudaMemcpy(plan->d_outpointdata, plan->h_outpointdata, 6*plan->noutpoints*plan->nsteps*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("sakrble, %s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->d_outpoint_pos, plan->h_outpoint_pos, 3*plan->noutpoints*sizeof(int), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }


    if (what == SV_SYNC_ALL && plan->h_tsfset!=NULL) {
	    err = cudaMemcpy(plan->d_tsf_jpool_e, plan->h_tsf_jpool_e, plan->tsf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
   	    err = cudaMemcpy(plan->d_tsf_jpool_h, plan->h_tsf_jpool_h, plan->tsf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
   	    err = cudaMemcpy(plan->d_tsf_jpool_epsilon, plan->h_tsf_jpool_epsilon, plan->tsf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
   	    err = cudaMemcpy(plan->d_tsf_jpool_sigma, plan->h_tsf_jpool_sigma, plan->tsf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_tsf_jpool_mu, plan->h_tsf_jpool_mu, plan->tsf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    	err = cudaMemcpy(plan->d_tsf_jpool_sigast, plan->h_tsf_jpool_sigast, plan->tsf_jpool_size*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    	err = cudaMemcpy(plan->d_tsfset, plan->h_tsfset, 21*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    	err = cudaMemcpy(plan->d_tsf_jpvals, plan->h_tsf_jpvals, 8*sizeof(float), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));     
    }
    if (what == SV_SYNC_ALL && plan->ff_size!=0) {
	    err = cudaMemcpy(plan->d_iset, plan->h_iset, plan->iset_size*sizeof(int), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_ff_ex, plan->h_ff_ex, plan->ff_size*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) 
                printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_ff_ey, plan->h_ff_ey, plan->ff_size*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) 
                printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_ff_ez, plan->h_ff_ez, plan->ff_size*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) 
                printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_ff_hlp_ex, plan->h_ff_hlp_ex, plan->ff_size*plan->nhlps*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) 
                printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_ff_hlp_ey, plan->h_ff_hlp_ey, plan->ff_size*plan->nhlps*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) 
                printf("%s\n", cudaGetErrorString(err));
            err = cudaMemcpy(plan->d_ff_hlp_ez, plan->h_ff_hlp_ez, plan->ff_size*plan->nhlps*sizeof(float), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) 
                printf("%s\n", cudaGetErrorString(err));
    }


    if (what == SV_SYNC_ALL && plan->nsums > 0) {
        err = cudaMemcpy(plan->d_sums, plan->h_sums, plan->nsums*plan->nsteps*sizeof(float), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_epsilon, plan->h_sum_epsilon, plan->nsums*sizeof(float), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_mu, plan->h_sum_mu, plan->nsums*sizeof(float), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_sigma, plan->h_sum_sigma, plan->nsums*sizeof(float), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_sigast, plan->h_sum_sigast, plan->nsums*sizeof(float), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_mode, plan->h_sum_mode, plan->nsums*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_i0, plan->h_sum_i0, plan->nsums*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_i1, plan->h_sum_i1, plan->nsums*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_j0, plan->h_sum_j0, plan->nsums*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_sum_j1, plan->h_sum_j1, plan->nsums*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }
     if (what == SV_SYNC_ALL && plan->nmat>0) {
        err = cudaMemcpy(plan->d_mat, plan->h_mat, plan->xres*plan->yres*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_mattype, plan->h_mattype, plan->nmat*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
        err = cudaMemcpy(plan->d_mattab, plan->h_mattab, 18*plan->nmat*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));

    }
    if (what == SV_SYNC_ALL) {
        err = cudaMemcpy(plan->d_bnds, plan->h_bnds, 6*sizeof(int), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));         
    }

    err = cudaThreadSynchronize();
    if (err != cudaSuccess) 
        printf("E returned \"%s\"\n", cudaGetErrorString(err));
    //printf("gpu_copyto finished\n");
}

void gpu_copyfrom(SvGpuPlan *plan)
{
    cudaError_t err;
    int dev;
    err = cudaGetDevice(&dev);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));


    /*collect hlps for farfield, if there are any*/
      if (plan->nhlps) {
	    cudaError_t err = wrap_ffKernel_gethlps(plan);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    } 

    err = cudaMemcpy(plan->h_ex, plan->d_ex, plan->xres*plan->yres*sizeof(float), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
    err = cudaMemcpy(plan->h_ey, plan->d_ey, plan->xres*plan->yres*sizeof(float), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
    err = cudaMemcpy(plan->h_ez, plan->d_ez, plan->xres*plan->yres*sizeof(float), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
    err = cudaMemcpy(plan->h_hx, plan->d_hx, plan->xres*plan->yres*sizeof(float), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
    err = cudaMemcpy(plan->h_hy, plan->d_hy, plan->xres*plan->yres*sizeof(float), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
    err = cudaMemcpy(plan->h_hz, plan->d_hz, plan->xres*plan->yres*sizeof(float), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));

    if (plan->ff_size!=0) {
	    err = cudaMemcpy(plan->h_iset, plan->d_iset, plan->iset_size*sizeof(int), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->h_ff_ex, plan->d_ff_ex, plan->ff_size*sizeof(float), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->h_ff_ey, plan->d_ff_ey, plan->ff_size*sizeof(float), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->h_ff_ez, plan->d_ff_ez, plan->ff_size*sizeof(float), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }

    if (plan->noutpoints > 0) {
  	    err = cudaMemcpy(plan->h_outpointdata, plan->d_outpointdata, 6*plan->noutpoints*plan->nsteps*sizeof(float), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
	    err = cudaMemcpy(plan->h_outpoint_pos, plan->d_outpoint_pos, 3*plan->noutpoints*sizeof(int), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }

    if (plan->nsums > 0) {
	    err = cudaMemcpy(plan->h_sums, plan->d_sums, plan->nsums*plan->nsteps*sizeof(float), cudaMemcpyDeviceToHost);
	    if (err != cudaSuccess) 
            printf("%s\n", cudaGetErrorString(err));
    }
		
    err = cudaThreadSynchronize();
    if (err != cudaSuccess) 
        printf("E returned \"%s\"\n", cudaGetErrorString(err));

}

void gpu_ystep_e(SvGpuPlan *plan)
{
    cudaError_t err = wrap_eKernel(plan);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}

void gpu_ystep_h(SvGpuPlan *plan)
{
    cudaError_t err = wrap_hKernel(plan);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}

void gpu_liao_cpy(SvGpuPlan *plan)
{
    cudaError_t err;
    
    plan->step += 1;
    err = wrap_liao_cpybnd(plan);
    
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}

void gpu_liao_bnd(SvGpuPlan *plan)
{
    cudaError_t err = wrap_liao_applybnd(plan);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}

void gpu_srce(SvGpuPlan *plan)
{
    cudaError_t err = wrap_srceKernel(plan, plan->source_i, plan->source_j, plan->source_k, plan->source_ex, plan->source_ey, plan->source_ez);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}

void gpu_tsf_estep(SvGpuPlan *plan)
{
    cudaError_t err = wrap_tsf_e_Kernel(plan);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}

void gpu_tsf_hstep(SvGpuPlan *plan)
{
    cudaError_t err = wrap_tsf_h_Kernel(plan);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}


void gpu_mbndx_e(SvGpuPlan *plan)
{
    cudaError_t err = wrap_mbndx_e_Kernel(plan);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}

void gpu_mbndy_e(SvGpuPlan *plan)
{
    cudaError_t err = wrap_mbndy_e_Kernel(plan);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}

void gpu_mbndx_h(SvGpuPlan *plan)
{
    cudaError_t err = wrap_mbndx_h_Kernel(plan);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}

void gpu_mbndy_h(SvGpuPlan *plan)
{
    cudaError_t err = wrap_mbndy_h_Kernel(plan);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}

void gpu_tsf_jstep(SvGpuPlan *plan)
{
    cudaError_t err = wrap_tsf_jstepKernel(plan, plan->source_ex);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}

void gpu_sum(SvGpuPlan *plan)
{
    cudaError_t err = wrap_fastsumKernel(plan);//wrap_fastsumKernel(plan);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}

void gpu_opnt(SvGpuPlan *plan)
{
    cudaError_t err = wrap_outpointKernel(plan);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}

void gpu_ff(SvGpuPlan *plan)
{
    cudaError_t err = wrap_ffKernel(plan);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}
void gpu_ffa(SvGpuPlan *plan)
{
    cudaError_t err = wrap_ffKernel_hlps(plan);
    if (err != cudaSuccess) 
        printf("%s\n", cudaGetErrorString(err));
}

void *gpu_main(void *data)
{
    SvGpuAPlan *aplan = (SvGpuAPlan *)data;
    SvGpuPlan *plan = aplan->plan;
    SvGTType todo;

    while (1) {
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
	        case SV_GT_TSF_JSTEP:
	    	    //printf("J step on device %d\n", plan->device);
	    	    gpu_tsf_jstep(plan);
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
           	case SV_GT_FFA:
	    	    gpu_ffa(plan);
	    	    g_queue_pop_tail(aplan->todos);
	    	    break;		        
            case SV_GT_NOTHING:
	    	//printf("Nothing on device %d\n", plan->device);
//	    	    sleep(1);
	    	    g_queue_pop_tail(aplan->todos);
	    	break;
                default:
	    	break;	    			    
	    }
    }
}


/*allocate structures and start worker threads*/
SvGpu* sv_gpu_new(SvPool *mp, SvSet *set)
{
    int i, n, splitsize;
    cudaError_t err;
    GError *gerr=NULL;
    struct cudaDeviceProp prop;
    SvGpuPlan *plan;
    SvGpu* gpu;
    int maxthreads[MAX_GPUS];

    int ngpu = 0, ngset=0, nprescribed, sel_gpu[MAX_GPUS];

    if (set->sc.ngpu == 0) 
        return NULL;

    gpu = (SvGpu*)g_malloc(sizeof(SvGpu));

    if (set->sc.devicequery) {
	    printf("running getdevicecount\n");
	    err = cudaGetDeviceCount(&ngpu);
	    printf("finished\n");
	    if (err != cudaSuccess) 
            fprintf(stderr, "Error: %s\n", cudaGetErrorString(err));
	    if (set->sc.verbose) 
            printf("Found %d GPUs\n", ngpu);


	    for (i=0; i<ngpu; i++) {
		    err = cudaGetDeviceProperties(&prop, i);
		    if (set->sc.verbose>1) {
			    printf("The Properties of the Device with ID %d are\n", i);
			    printf("Device Name             : %s\n", prop.name);
			    printf("Device Memory Size      : %u\n", (guint)prop.totalGlobalMem);
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
        ngpu = set->sc.ngpu;   /*no test for gpus, just believe user, unsafe*/
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
    if (nprescribed==0 && set->sc.ngpu>0) { /*just pick first few devices*/
        if (set->sc.ngpu <= ngpu) 
            ngset = set->sc.ngpu;
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

    if (set->sc.verbose > 1) 
        printf("%d devices used: %d%d%d%d %d%d%d%d %d%d%d%d %d%d%d%d\n", 
		    ngset, sel_gpu[0], sel_gpu[1], sel_gpu[2], sel_gpu[3],
		    sel_gpu[4], sel_gpu[5], sel_gpu[6], sel_gpu[7],
		    sel_gpu[8], sel_gpu[9], sel_gpu[10], sel_gpu[11],
            sel_gpu[12], sel_gpu[13], sel_gpu[14], sel_gpu[15]);

    /*perform allocation*/
    gpu->ngset = ngset;
    gpu->gset = (SvGpuTSet *)g_malloc(ngset*sizeof(SvGpuTSet));
    for (n=0; n<ngset; n++) {
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

        plan->bndx0 = set->sb.bx0;
        plan->bndxn = set->sb.bxn;
        plan->bndy0 = set->sb.by0;
        plan->bndyn = set->sb.byn;

        plan->xres = set->sp.xres;
        plan->yres = set->sp.yres;
        if (ngset==1) {
            plan->xfrom = 0;
            plan->xto = set->sp.xres;
        } else {
            plan->xfrom = MAX(0, (n*set->sp.xres/ngset) - 1);
            plan->xto = MIN(set->sp.xres, ((n+1)*set->sp.xres/ngset) + 1);
        }
        plan->xres = plan->xto - plan->xfrom;

        plan->dir = 1;
        plan->dx = (float)set->sp.dx;
        plan->dy = (float)set->sp.dy;
        plan->dt = (float)set->plan.dt;
        plan->nsteps = set->sc.nsteps;
        plan->nhlps = 0;

        plan->tmmode = set->tmmode;

        if (set->sc.verbose > 1) 
            printf("Allocating %dx%d\n", plan->xres, plan->yres);
        /*
           if (ugpu>= 0 || ngpu==1) plan[n].bnd = 3;
           else if (n==0) plan[n].bnd = 1;
           else if (n==(ngpu-1)) plan[n].bnd = 2;
           else plan[n].bnd = 0;
           */
        plan->h_ex = (float *)g_malloc(plan->xres*plan->yres*sizeof(float));
        plan->h_ey = (float *)g_malloc(plan->xres*plan->yres*sizeof(float));
        plan->h_ez = (float *)g_malloc(plan->xres*plan->yres*sizeof(float));
        plan->h_hx = (float *)g_malloc(plan->xres*plan->yres*sizeof(float));
        plan->h_hy = (float *)g_malloc(plan->xres*plan->yres*sizeof(float));
        plan->h_hz = (float *)g_malloc(plan->xres*plan->yres*sizeof(float));

        if (set->plan.matmode==SV_MATMODE_FULL || set->plan.matmode==SV_MATMODE_ELECTRIC) {
            plan->h_epsilon = (float *)g_malloc(plan->xres*plan->yres*sizeof(float));
            plan->h_sigma   = (float *)g_malloc(plan->xres*plan->yres*sizeof(float));
        }
        if (set->plan.matmode==SV_MATMODE_FULL || set->plan.matmode==SV_MATMODE_MAGNETIC) {
            plan->h_mu     =  (float *)g_malloc(plan->xres*plan->yres*sizeof(float));
            plan->h_sigast =  (float *)g_malloc(plan->xres*plan->yres*sizeof(float));
        }
        if (mp->nmat) {
            plan->h_mat     = (int *)g_malloc(plan->xres*plan->yres*sizeof(int));
            plan->h_mattype = (int *)g_malloc(mp->nmat*sizeof(int));
            plan->h_mattab  = (float *)g_malloc(18*mp->nmat*sizeof(float));
            plan->nmat = mp->nmat;
            for (i=0; i<mp->nmat; i++) {
                plan->h_mattype[i] = mp->mats[i].type;
                plan->h_mattab[18*i + SV_GM_EPSILON]    = (float)mp->mats[i].epsilon;
                plan->h_mattab[18*i + SV_GM_MU]         = (float)mp->mats[i].mu;
                plan->h_mattab[18*i + SV_GM_SIGMA]      = (float)mp->mats[i].sigma;
                plan->h_mattab[18*i + SV_GM_SIGAST]     = (float)mp->mats[i].sigast;
                plan->h_mattab[18*i + SV_GM_DRUDE_OMEGA_P] = (float)mp->mats[i].drude_omega_p;
                plan->h_mattab[18*i + SV_GM_DRUDE_NU]   = (float)mp->mats[i].drude_nu;
                plan->h_mattab[18*i + SV_GM_CP3_A0]     = (float)mp->mats[i].cp3_a[0];
                plan->h_mattab[18*i + SV_GM_CP3_A1]     = (float)mp->mats[i].cp3_a[1];
                plan->h_mattab[18*i + SV_GM_CP3_A2]     = (float)mp->mats[i].cp3_a[2];
                plan->h_mattab[18*i + SV_GM_CP3_PHI0]   = (float)mp->mats[i].cp3_phi[0];
                plan->h_mattab[18*i + SV_GM_CP3_PHI1]   = (float)mp->mats[i].cp3_phi[1];
                plan->h_mattab[18*i + SV_GM_CP3_PHI2]   = (float)mp->mats[i].cp3_phi[2];
                plan->h_mattab[18*i + SV_GM_CP3_OMEGA0] = (float)mp->mats[i].cp3_omega[0];
                plan->h_mattab[18*i + SV_GM_CP3_OMEGA1] = (float)mp->mats[i].cp3_omega[1];
                plan->h_mattab[18*i + SV_GM_CP3_OMEGA2] = (float)mp->mats[i].cp3_omega[2];
                plan->h_mattab[18*i + SV_GM_CP3_GAMMA0] = (float)mp->mats[i].cp3_gamma[0];
                plan->h_mattab[18*i + SV_GM_CP3_GAMMA1] = (float)mp->mats[i].cp3_gamma[1];
                plan->h_mattab[18*i + SV_GM_CP3_GAMMA2] = (float)mp->mats[i].cp3_gamma[2];

            }
        } else {
            plan->h_mat = NULL;
            plan->h_mattype = NULL;
            plan->h_mattab = NULL;
            plan->nmat = 0;
        }

	    plan->h_isplrc = FALSE;
	    for (i=0; i<mp->nmat; i++) {
	    	if (mp->mats[i].type==SV_MAT_DRUDE || mp->mats[i].type==SV_MAT_CP3) 
	    		plan->h_isplrc = TRUE;
	    }



        plan->noutpoints = set->so.npnts;
        plan->h_outpoint_pos = (int *)g_malloc(3*plan->noutpoints*sizeof(int));
        for (i=0; i<plan->noutpoints; i++) {
            plan->h_outpoint_pos[3*i] = set->so.pnts[i].i;
            plan->h_outpoint_pos[3*i + 1] = set->so.pnts[i].j;
//            plan->h_outpoint_pos[3*i + 2] = set->so.pnts[i].k;
        }
        plan->h_outpointdata = (float *)g_malloc(6*plan->noutpoints*set->sc.nsteps*sizeof(float));
        for (i=0; i<(6*plan->noutpoints*plan->nsteps); i++) 
            plan->h_outpointdata[i] = 0;

        if (mp->src->tsf) {
            plan->h_tsfset = (float *)g_malloc(21*sizeof(float));
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
                plan->h_tsf_jpool_mu[i] = (float)mp->src->tsf->mu;
                plan->h_tsf_jpool_epsilon[i] = (float)mp->src->tsf->epsilon;
            }
            for (i=0; i<8; i++) 
                plan->h_tsf_jpvals[i] = 0;
            plan->h_tsfset[0] = (float)mp->src->tsf->i0;
            plan->h_tsfset[1] = (float)mp->src->tsf->j0;
            plan->h_tsfset[2] = 0.0f;
            plan->h_tsfset[3] = (float)mp->src->tsf->i1;
            plan->h_tsfset[4] = (float)mp->src->tsf->j1;
            plan->h_tsfset[5] = 0.0f;
            plan->h_tsfset[6] = (float)mp->src->tsf->theta;
            plan->h_tsfset[7] = (float)mp->src->tsf->phi;
            plan->h_tsfset[8] = (float)mp->src->tsf->psi;
            plan->h_tsfset[9] = (float)mp->src->tsf->corr;
            plan->h_tsfset[10] = 0.0f;
            plan->h_tsfset[11] = 0.0f;
            plan->h_tsfset[12] = 0.0f;
            plan->h_tsfset[13] = 0.0f;
            plan->h_tsfset[14] = 0.0f;
            plan->h_tsfset[15] = 0.0f;
            //16..20 are left for future for compatibility with 3d 
        } else 
            plan->h_tsfset = NULL;


        if (set->sf.nrs) {
	        plan->ff_size = 0;
            splitsize = MAX(set->sf.i1-set->sf.i0, set->sf.j1-set->sf.j0);
            
	        if (set->sf.nrs<=256) 
                plan->nhlps = (MIN(splitsize, plan->maxthreads)/set->sf.nrs);
            else 
                plan->nhlps = 0;

            if (set->sc.verbose > 1) 
                printf("We will use %d helpers for NFFF computation.\n", plan->nhlps);


	        plan->iset_size = 37 + 5*set->sf.nrs;
	        plan->h_iset = (int *)g_malloc(plan->iset_size*sizeof(int));

            plan->h_iset[0] = set->sf.i0;
            plan->h_iset[1] = set->sf.i1;
            plan->h_iset[2] = set->sf.j0;
            plan->h_iset[3] = set->sf.j1;
            plan->h_iset[4] = 0;
            plan->h_iset[5] = 0;
            /*array positions 6..35 are reserved for skipping boundary in future*/
            plan->h_iset[36] = set->sf.nrs; 
	        for (i = 0; i<set->sf.nrs; i++) {
		     //   plan->h_iset[5*i+37] = mp->farfield->ndata;
		     //   plan->ff_size += mp->farfield->ndata;
		     //   plan->h_iset[5*i+38] = (int)mp->farfield->rpoints[i].istart; //istart
		     //   plan->h_iset[5*i+39] = (int)set->sf.ri[i];
		     //   plan->h_iset[5*i+40] = (int)set->sf.rj[i];
		     //   plan->h_iset[5*i+41] = (int)set->sf.rk[i];
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
	    } else 
            plan->ff_size = 0;  

 
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
                  
            for (i=0; i<(set->so.nsums*set->sc.nsteps); i++)
                plan->h_sums[i] = 0;

            for (i=0; i<(set->so.nsums); i++) {
                plan->h_sum_epsilon[i]  = (float)set->so.sums[i].epsilon;
                plan->h_sum_mu[i]       = (float)set->so.sums[i].mu;
                plan->h_sum_sigma[i]    = (float)set->so.sums[i].sigma;
                plan->h_sum_sigast[i]   = (float)set->so.sums[i].sigast;
                plan->h_sum_mode[i]     = set->so.sums[i].component;
                plan->h_sum_i0[i] = set->so.sums[i].i0;
                plan->h_sum_i1[i] = set->so.sums[i].i1;
                plan->h_sum_j0[i] = set->so.sums[i].j0;
                plan->h_sum_j1[i] = set->so.sums[i].j1;

            }
        } else 
            plan->nsums = 0;

        plan->mb_bx0 = set->smb.bx0;
        plan->mb_by0 = set->smb.by0;
        plan->mb_bxn = set->smb.bxn;
        plan->mb_byn = set->smb.byn;
        plan->mb_bx0pos = set->smb.bx0pos;
        plan->mb_by0pos = set->smb.by0pos;
        plan->mb_bxnpos = set->smb.bxnpos;
        plan->mb_bynpos = set->smb.bynpos;


        /*start thread*/
        g_thread_create((GThreadFunc)gpu_main, (void *)(gpu->gset[n].aplan), TRUE, &gerr);
    }

    return gpu;
}

#endif

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
