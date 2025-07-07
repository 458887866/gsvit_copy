
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


/*  gpu_kernels.h : 
 *  all the algorithms that run on GPU and are compied with nvcc
 */

#ifdef UCUDA

#include <cuda_runtime_api.h>
      
/*what should be/is computed*/
typedef enum {
	SV_GT_INIT = 0,
	SV_GT_COPYTO = 1,
	SV_GT_COPYFROM = 2,
	SV_GT_YSTEPE = 3,
	SV_GT_YSTEPH = 4,
	SV_GT_LIAO_CPY = 5,
	SV_GT_LIAO_BND = 6,
	SV_GT_SRCE_POINT = 7,
	SV_GT_SRCH_POINT = 8,
        SV_GT_TSF_ESTEP      = 9,
        SV_GT_TSF_HSTEP      = 10,
        SV_GT_TSF_JSTEP    = 11,
        SV_GT_MBNDX_E    = 12,
        SV_GT_MBNDY_E    = 13,
        SV_GT_MBNDZ_E    = 14,
        SV_GT_MBNDX_H    = 15,
        SV_GT_MBNDY_H    = 16,
        SV_GT_MBNDZ_H    = 17,
        SV_GT_SUM        = 18,
        SV_GT_FF         = 19,
        SV_GT_FFA        = 20,
        SV_GT_CPML_ESTEP  = 21,
        SV_GT_CPML_HSTEP  = 22,
        SV_GT_OPNT  = 23,
        SV_GT_SF_JSTEP    = 24,
        SV_GT_PFF        = 25,
        SV_GT_TSFF_ESTEP   = 26,
        SV_GT_TSFF_HSTEP   = 27,
        SV_GT_TSFF_JSTEP   = 28,
        SV_GT_FORCE        = 29,
	SV_GT_NOTHING = 30
} SvGTType;

typedef enum {
        SV_GM_EPSILON = 0,
        SV_GM_MU = 1,
        SV_GM_SIGMA = 2,
        SV_GM_SIGAST = 3,
        SV_GM_DRUDE_OMEGA_P = 4,
        SV_GM_DRUDE_NU = 5,
        SV_GM_CP3_A0 = 6,
        SV_GM_CP3_A1 = 7,
        SV_GM_CP3_A2 = 8,
        SV_GM_CP3_PHI0 = 9,
        SV_GM_CP3_PHI1 = 10,
        SV_GM_CP3_PHI2 = 11,
        SV_GM_CP3_OMEGA0 = 12,
        SV_GM_CP3_OMEGA1 = 13,
        SV_GM_CP3_OMEGA2 = 14,
        SV_GM_CP3_GAMMA0 = 15,
        SV_GM_CP3_GAMMA1 = 16,
        SV_GM_CP3_GAMMA2 = 17,
        SV_GM_ADE_A0 = 18,
        SV_GM_ADE_A1 = 19,
        SV_GM_ADE_A2 = 20,
        SV_GM_ADE_BP0 = 21, //two members
        SV_GM_ADE_BP1 = 23,
        SV_GM_ADE_BP2 = 25, 
        SV_GM_ADE_BP3 = 27,
        SV_GM_ADE_BP4 = 29,
        SV_GM_ADE_C0 = 31,
        SV_GM_ADE_C1 = 32,
        SV_GM_ADE_C2 = 33,
        SV_GM_ADE_C3 = 34,
        SV_GM_PLRC_D_CHI = 35,
        SV_GM_PLRC_D_XI = 36,
        SV_GM_PLRC_D_DCHI = 37,
        SV_GM_PLRC_D_DXI = 38,
        SV_GM_PLRC_P_CHI = 39,
        SV_GM_PLRC_P_XI = 41,
        SV_GM_PLRC_P_DCHIR = 43,
        SV_GM_PLRC_P_DXIR = 45,
        SV_GM_PLRC_P_DCHII = 47,
        SV_GM_PLRC_P_DXII = 49,
        SV_GM_N = 51
} SvGMType;


typedef struct{
    int device;
    int xres;
    int yres;
    int zres;
    int xfrom;
    int xto;
    int step;
    int matmode;

    int bndx0;
    int bndxn;
    int bndy0;
    int bndyn;
    int bndz0;
    int bndzn;

    float dx;
    float dy;
    float dz;
    float dt;

    /*newly allocated host structures for device sync*/
    float *h_ex;
    float *h_ey;
    float *h_ez;
    float *h_hx;
    float *h_hy;
    float *h_hz;
    float *h_epsilon;
    float *h_mu;
    float *h_sigma;
    float *h_sigast;
    float *h_ff_ex;
    float *h_ff_ey;
    float *h_ff_ez;
    float *h_ff_hlp_ex;
    float *h_ff_hlp_ey;
    float *h_ff_hlp_ez;
    float *h_pff_ex;
    float *h_pff_ey;
    float *h_pff_ez;
    float *h_pff_hlp_ex;
    float *h_pff_hlp_ey;
    float *h_pff_hlp_ez;
    float *h_sums;
    float *h_sum_epsilon;
    float *h_sum_mu;
    float *h_sum_sigma;
    float *h_sum_sigast;
    int *h_sum_mode;
    int *h_sum_i0;
    int *h_sum_i1;
    int *h_sum_j0;
    int *h_sum_j1;
    int *h_sum_k0;
    int *h_sum_k1;
    int *h_iset;
    int *h_piset;
    float *h_forces;
    int *h_force_i0;
    int *h_force_i1;
    int *h_force_j0;
    int *h_force_j1;
    int *h_force_k0;
    int *h_force_k1;
    
 
    float *h_outpointdata;
    float *d_outpointdata;
    int *h_outpoint_pos;
    int *d_outpoint_pos;
    int noutpoints;
    

    int *h_bnds;
    int *h_depths;

    int *d_bnds;
    int *d_depths;

    int cpml_depth_bx0;
    int cpml_depth_bxn;
    int cpml_depth_by0;
    int cpml_depth_byn;
    int cpml_depth_bz0;
    int cpml_depth_bzn;


    float *h_cpml_be_x0;
    float *h_cpml_ce_x0;
    float *h_cpml_bh_x0;
    float *h_cpml_ch_x0;
    float *h_cpml_kappae_x0;
    float *h_cpml_kappah_x0;

    float *h_cpml_be_xn;
    float *h_cpml_ce_xn;
    float *h_cpml_bh_xn;
    float *h_cpml_ch_xn;
    float *h_cpml_kappae_xn;
    float *h_cpml_kappah_xn;

    float *h_cpml_be_y0;
    float *h_cpml_ce_y0;
    float *h_cpml_bh_y0;
    float *h_cpml_ch_y0;
    float *h_cpml_kappae_y0;
    float *h_cpml_kappah_y0;

    float *h_cpml_be_yn;
    float *h_cpml_ce_yn;
    float *h_cpml_bh_yn;
    float *h_cpml_ch_yn;
    float *h_cpml_kappae_yn;
    float *h_cpml_kappah_yn;

    float *h_cpml_be_z0;
    float *h_cpml_ce_z0;
    float *h_cpml_bh_z0;
    float *h_cpml_ch_z0;
    float *h_cpml_kappae_z0;
    float *h_cpml_kappah_z0;

    float *h_cpml_be_zn;
    float *h_cpml_ce_zn;
    float *h_cpml_bh_zn;
    float *h_cpml_ch_zn;
    float *h_cpml_kappae_zn;
    float *h_cpml_kappah_zn;




    /*pointers to device parameters*/
    float *d_ex;
    float *d_ey;
    float *d_ez;
    float *d_hx;
    float *d_hy;
    float *d_hz;
    float *d_epsilon;
    float *d_mu;
    float *d_sigma;
    float *d_sigast;
    float *d_liao_x0;
    float *d_liao_xn;
    float *d_liao_y0;
    float *d_liao_yn;
    float *d_liao_z0;
    float *d_liao_zn;
    float *d_ff_ex;  /*same integreator (in floats) as rampnt.at(m)).ex*/
    float *d_ff_ey;
    float *d_ff_ez;
    float *d_ff_hlp_ex;
    float *d_ff_hlp_ey;
    float *d_ff_hlp_ez;    
    float *d_pff_ex; 
    float *d_pff_ey;
    float *d_pff_ez;
    float *d_pff_hlp_ex;
    float *d_pff_hlp_ey;
    float *d_pff_hlp_ez;    
    float *d_sums;
    float *d_sum_epsilon;
    float *d_sum_mu;
    float *d_sum_sigma;
    float *d_sum_sigast;
    float *d_sum_accumulator;
    int *d_sum_mode;
    int *d_sum_i0;
    int *d_sum_i1;
    int *d_sum_j0;
    int *d_sum_j1;
    int *d_sum_k0;
    int *d_sum_k1;
    int nhlps;
    int nsums;
    int nsteps;
    int *d_iset;
    int *d_piset;
    int iset_size;
    int piset_size;
    int ff_size;
    int pff_size;
    float *d_forces;
    int *d_force_i0;
    int *d_force_i1;
    int *d_force_j0;
    int *d_force_j1;
    int *d_force_k0;
    int *d_force_k1;
    int nforces;

    float *h_tsf_jpool_e;
    float *h_tsf_jpool_h;
    float *h_tsf_jpool_epsilon;
    float *h_tsf_jpool_sigma;
    float *h_tsf_jpool_mu;
    float *h_tsf_jpool_sigast;
    float *h_tsf_jpvals;
    float *h_tsfset;
    float *h_tsf_fiberset;
    int tsf_jpool_size;
    float *d_tsf_jpool_e;
    float *d_tsf_jpool_h;
    float *d_tsf_jpool_epsilon;
    float *d_tsf_jpool_sigma;
    float *d_tsf_jpool_mu;
    float *d_tsf_jpool_sigast;
    float *d_tsf_jpvals;
    float *d_tsfset;
    float *d_tsf_fiberset;
    int dir;

    float *d_plrcx; //seven fields sequetially (always if any plrc is present), for debye cp and cp3 models
    float *d_plrcy;
    float *d_plrcz;
    int h_isplrc;

    float *d_pxp; //ADE storage fields, x[0], y[0], z[0], x[1], y[1], z[1] sequentially
    float *d_px;
    float *d_dpxp;
    float *d_dpx;
    float *d_exp;
    float *d_eyp;
    float *d_ezp;
    int h_isade; 

    int tsff_nxm;
    float *h_tsff_jpool_e;
    float *h_tsff_jpool_h;
    float *h_tsff_jpool_epsilon;
    float *h_tsff_jpool_sigma;
    float *h_tsff_jpool_mu;
    float *h_tsff_jpool_sigast;
    float *h_tsff_jpvals;
    float *h_tsffset;
    int tsff_jpool_size;
    float *d_tsff_jpool_e;
    float *d_tsff_jpool_h;
    float *d_tsff_jpool_epsilon;
    float *d_tsff_jpool_sigma;
    float *d_tsff_jpool_mu;
    float *d_tsff_jpool_sigast;
    float *d_tsff_jpvals;
    float *d_tsffset;
    float *tsff_data;
    float tsff_freal;
 
    float *h_sf_jpool_e;
    float *h_sf_jpool_h;
    float *h_sf_jpool_epsilon;
    float *h_sf_jpool_sigma;
    float *h_sf_jpool_mu;
    float *h_sf_jpool_sigast;
    float *h_sf_jpvals;
    float *h_sfset;
    int sf_jpool_size;
    float *d_sf_jpool_e;
    float *d_sf_jpool_h;
    float *d_sf_jpool_epsilon;
    float *d_sf_jpool_sigma;
    float *d_sf_jpool_mu;
    float *d_sf_jpool_sigast;
    float *d_sf_jpvals;
    float *d_sfset;
    float *d_sf_fiberset;
 
    float *h_cpml_p_x0;
    float *h_cpml_p_xn;
    float *h_cpml_p_y0;
    float *h_cpml_p_yn;
    float *h_cpml_p_z0;
    float *h_cpml_p_zn;


    float *d_cpml_p_x0;
    float *d_cpml_p_xn;
    float *d_cpml_p_y0;
    float *d_cpml_p_yn;
    float *d_cpml_p_z0;
    float *d_cpml_p_zn;

    float *d_cpml_be_x0;
    float *d_cpml_ce_x0;
    float *d_cpml_bh_x0;
    float *d_cpml_ch_x0;
    float *d_cpml_kappae_x0;
    float *d_cpml_kappah_x0;

    float *d_cpml_be_xn;
    float *d_cpml_ce_xn;
    float *d_cpml_bh_xn;
    float *d_cpml_ch_xn;
    float *d_cpml_kappae_xn;
    float *d_cpml_kappah_xn;

    float *d_cpml_be_y0;
    float *d_cpml_ce_y0;
    float *d_cpml_bh_y0;
    float *d_cpml_ch_y0;
    float *d_cpml_kappae_y0;
    float *d_cpml_kappah_y0;

    float *d_cpml_be_yn;
    float *d_cpml_ce_yn;
    float *d_cpml_bh_yn;
    float *d_cpml_ch_yn;
    float *d_cpml_kappae_yn;
    float *d_cpml_kappah_yn;

    float *d_cpml_be_z0;
    float *d_cpml_ce_z0;
    float *d_cpml_bh_z0;
    float *d_cpml_ch_z0;
    float *d_cpml_kappae_z0;
    float *d_cpml_kappah_z0;

    float *d_cpml_be_zn;
    float *d_cpml_ce_zn;
    float *d_cpml_bh_zn;
    float *d_cpml_ch_zn;
    float *d_cpml_kappae_zn;
    float *d_cpml_kappah_zn;

    int *d_mat;
    int *h_mat;
    int *d_mattype;
    int *h_mattype;
    float *d_mattab;
    float *h_mattab;
    int nmat;

    int source_i;
    int source_j;
    int source_k;
    float source_ex; //this is used also for tsf
    float source_ey;
    float source_ez;
    float source_hx;
    float source_hy;
    float source_hz;
    int is_ex;
    int is_ey;
    int is_ez;
    int is_hx;
    int is_hy;
    int is_hz;

    int mb_bx0;
    int mb_bxn;
    int mb_by0;
    int mb_byn;
    int mb_bz0;
    int mb_bzn;
    int mb_bx0pos;
    int mb_by0pos;
    int mb_bz0pos;
    int mb_bxnpos;
    int mb_bynpos;
    int mb_bznpos;

    int maxthreads;
} SvGpuPlan;

cudaError_t wrap_eKernel(SvGpuPlan *plan);


cudaError_t wrap_hKernel(SvGpuPlan *plan);

//cudaError_t wrap_yKernel(TGPUplan *plan);

cudaError_t wrap_liao_cpybnd(SvGpuPlan *plan);

cudaError_t wrap_liao_applybnd(SvGpuPlan *plan);

cudaError_t wrap_srceKernel(SvGpuPlan *plan, int i, int j, int k, float ex, float ey, float ez);
cudaError_t wrap_srchKernel(SvGpuPlan *plan, int i, int j, int k, float ex, float ey, float ez);

cudaError_t wrap_tsf_jstepKernel(SvGpuPlan *plan, float e);
cudaError_t wrap_sf_jstepKernel(SvGpuPlan *plan, float e);
cudaError_t wrap_tsf_e_Kernel(SvGpuPlan *plan);
cudaError_t wrap_tsf_h_Kernel(SvGpuPlan *plan);

cudaError_t wrap_tsff_jstepKernel(SvGpuPlan *plan, float e);
cudaError_t wrap_tsff_e_Kernel(SvGpuPlan *plan);
cudaError_t wrap_tsff_h_Kernel(SvGpuPlan *plan);


cudaError_t wrap_mbndx_e_Kernel(SvGpuPlan *plan);
cudaError_t wrap_mbndy_e_Kernel(SvGpuPlan *plan);
cudaError_t wrap_mbndz_e_Kernel(SvGpuPlan *plan);
cudaError_t wrap_mbndx_h_Kernel(SvGpuPlan *plan);
cudaError_t wrap_mbndy_h_Kernel(SvGpuPlan *plan);
cudaError_t wrap_mbndz_h_Kernel(SvGpuPlan *plan);

cudaError_t wrap_outpointKernel(SvGpuPlan *plan);

cudaError_t wrap_sumKernel(SvGpuPlan *plan);
cudaError_t wrap_fastsumKernel(SvGpuPlan *plan);
cudaError_t wrap_forceKernel(SvGpuPlan *plan);


cudaError_t wrap_ffKernel(SvGpuPlan *plan);
cudaError_t wrap_ffKernel_hlps(SvGpuPlan *plan);
cudaError_t wrap_ffKernel_gethlps(SvGpuPlan *plan);

cudaError_t wrap_pffKernel(SvGpuPlan *plan);
cudaError_t wrap_pffKernel_hlps(SvGpuPlan *plan);
cudaError_t wrap_pffKernel_gethlps(SvGpuPlan *plan);

cudaError_t wrap_cpml_estep(SvGpuPlan *plan);
cudaError_t wrap_cpml_hstep(SvGpuPlan *plan);


//cudaError_t wrap_fiberKernel(TGPUplan *plan);


#endif

