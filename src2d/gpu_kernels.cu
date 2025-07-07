
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


/*  gpu_kernels.cu : 
 *  all the algorithms that run on GPU and are compied with nvcc
 *  FIXME: split this into more files
 */

#if defined (_WIN32) && defined (_MSC_VER)
#include "../config_win.h"
#else
#include "../config.h"
#endif
#ifdef UCUDA

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"
//#include <cutil.h>
#include <cuda_runtime_api.h>
//#include <multithreading.h>

extern "C" {
#include "gpu_kernels.h"
}

#define CU_PI 3.14159265358979323846

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

/*values in the main field arrays*/
#define CAT(i, j) ((j)*xres + (i))

/*values stored for Liao boundary*/
#define CLX(off, i)  (off*xres + (i)) 
#define CLXP(off, i) ((off + 1)*xres + (i))
#define CLY(off, j)  (off*yres + (j)) 
#define CLYP(off, j) ((off + 1)*yres + (j))

/*various constants to help searching in parameter arrays*/
#define EX 0
#define EY 2
#define EZ 4
#define HX 6
#define HY 8
#define HZ 10

#define I0 0
#define I1 1
#define J0 2
#define J1 3
#define K0 4
#define K1 5
#define SKIPI0 6
#define SKIPI0_JMIN 7
#define SKIPI0_KMIN 8
#define SKIPI0_JMAX 9
#define SKIPI0_KMAX 10
#define SKIPIN 11
#define SKIPIN_JMIN 12
#define SKIPIN_KMIN 13
#define SKIPIN_JMAX 14
#define SKIPIN_KMAX 15

#define SKIPJ0 16
#define SKIPJ0_IMIN 17
#define SKIPJ0_KMIN 18
#define SKIPJ0_IMAX 19
#define SKIPJ0_KMAX 20
#define SKIPJN 21
#define SKIPJN_IMIN 22
#define SKIPJN_KMIN 23
#define SKIPJN_IMAX 24
#define SKIPJN_KMAX 25

#define SKIPK0 26
#define SKIPK0_IMIN 27
#define SKIPK0_JMIN 28
#define SKIPK0_IMAX 29
#define SKIPK0_JMAX 30
#define SKIPKN 31
#define SKIPKN_IMIN 32
#define SKIPKN_JMIN 33
#define SKIPKN_IMAX 34
#define SKIPKN_JMAX 35

#define BNDX0 0
#define BNDXN 1
#define BNDY0 2
#define BNDYN 3
#define BNDZ0 4
#define BNDZN 5

#define NPNTS 36
#define PNPNTS 20
#define PNDIV 200

#define GETEPS(i) (18*i + SV_GM_EPSILON)
#define GETMU(i) (18*i + SV_GM_MU)
#define GETSIGMA(i) (18*i + SV_GM_SIGMA)
#define GETSIGAST(i) (18*i + SV_GM_SIGAST)
#define GETDRUDEO(i) (18*i + SV_GM_DRUDE_OMEGA_P)
#define GETDRUDEN(i) (18*i + SV_GM_DRUDE_NU)
#define GETCPA0(i) (18*i + SV_GM_CP3_A0)
#define GETCPA1(i) (18*i + SV_GM_CP3_A1)
#define GETCPA2(i) (18*i + SV_GM_CP3_A2)
#define GETCPPHI0(i) (18*i + SV_GM_CP3_PHI0)
#define GETCPPHI1(i) (18*i + SV_GM_CP3_PHI1)
#define GETCPPHI2(i) (18*i + SV_GM_CP3_PHI2)
#define GETCPOMEGA0(i) (18*i + SV_GM_CP3_OMEGA0)
#define GETCPOMEGA1(i) (18*i + SV_GM_CP3_OMEGA1)
#define GETCPOMEGA2(i) (18*i + SV_GM_CP3_OMEGA2)
#define GETCPGAMMA0(i) (18*i + SV_GM_CP3_GAMMA0)
#define GETCPGAMMA1(i) (18*i + SV_GM_CP3_GAMMA1)
#define GETCPGAMMA2(i) (18*i + SV_GM_CP3_GAMMA2)

#define MU_0 1.256637061435917295e-6
#define EPSILON_0 8.854187817620389850e-12


#define MEM_N 0
#define LS 299792458.0
#define JCPI 0.0795774715 
/*(1/4Pi)*/

#define SV_BOUNDARY_CPML 3

__device__ float k_dcomp(int i, int j, int k, int xres, int yres, int zres, float theta, float phi,
                         int i0, int i1, int j0, int j1, int k0, int k1);

__device__ float k_gex(float field, float theta, float phi, float psi);
__device__ float k_gey(float field, float theta, float phi, float psi);
__device__ float k_gez(float field, float theta, float phi, float psi);
__device__ float k_ghx(float field, float theta, float phi, float psi);
__device__ float k_ghy(float field, float theta, float phi, float psi);
__device__ float k_ghz(float field, float theta, float phi, float psi);
__device__ float get_dval(float *line, int res, float x);


/*vacuum or any material given by voxel-by-voxel set of optical parameters*/

__global__  void
eKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_epsilon, float *d_mu, float *d_sigma, float *d_sigast, int xres, int yres, int tmmode,
         float dx, float dy, float dt, int dir)
{
    int i, j;
    float ca, cb;
    long int pos;

    i = threadIdx.x;
    j = blockIdx.x;

    if (i<=0 || j<=0) return;
    if (i>=(xres) || j>=(yres)) return;

    pos = CAT(i, j);

    ca  = (1 - d_sigma[pos]*dt/2/d_epsilon[pos])/(1 + d_sigma[pos]*dt/2/d_epsilon[pos]);
    cb  = (dt/d_epsilon[pos]/dx)/(1 + d_sigma[pos]*dt/2/d_epsilon[pos]);


    if (tmmode) {
        
        d_ez[pos] = ca*d_ez[pos] + cb*((d_hy[pos] - d_hy[CAT(i-1, j)]) -
                            (d_hx[pos] - d_hx[CAT(i, j-1)]));
    }
    else {
        d_ex[pos] = ca*d_ex[pos] + cb*(d_hz[pos] - d_hz[CAT(i, j-1)]);

        d_ey[pos] = ca*d_ey[pos] - cb*(d_hz[pos] - d_hz[CAT(i-1, j)]);
    }

}

/*any material including tabulated material (voxel-by-voxel set of tabulated material parameters)*/
__global__  void
eKernel_tab(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_epsilon, float *d_mu, float *d_sigma, float *d_sigast, 
         int *d_mat, int *d_mattype, float *d_mattab, int nmat, int matmode,
         int xres, int yres, int tmmode,
         float dx, float dy, float dt, int dir)
{
    int i, j;
    float ca, cb, sigma, epsilon;
    long int pos;
    int mattype_xm, mattype_ym, mattype;

    i = threadIdx.x;
    j = blockIdx.x;

    if (i<=0 || j<=0) return;
    if (i>=(xres) || j>=(yres)) return;

    pos = CAT(i, j);

    if (i>0 && j>0) {
	    mattype_xm = d_mattype[d_mat[CAT(i-1, j)]];
	    mattype_ym = d_mattype[d_mat[CAT(i, j-1)]];
    }

    mattype = d_mattype[d_mat[pos]];

    if (d_mat[pos]==0 && !(matmode==0 || matmode==2)) { //= mattype 0 = linear material given pixel by pixel
	    sigma = d_sigma[pos];
	    epsilon = d_epsilon[pos];
    } else if (mattype==1) { //tabulated linear material, here should be also the cp3 and drude option
	    sigma = d_mattab[GETSIGMA(d_mat[pos])];
	    epsilon = d_mattab[GETEPS(d_mat[pos])];
    } else {
	    sigma = 0;
	    epsilon = EPSILON_0;
    }
    ca  = (1 - sigma*dt/2/epsilon)/(1 + sigma*dt/2/epsilon);
    cb  = (dt/epsilon/dx)/(1 + sigma*dt/2/epsilon);

    if (tmmode) {
        
        d_ez[pos] = ca*d_ez[pos] + cb*((d_hy[pos] - d_hy[CAT(i-1, j)]) -
                            (d_hx[pos] - d_hx[CAT(i, j-1)]));
    }
    else {
        d_ex[pos] = ca*d_ex[pos] + cb*(d_hz[pos] - d_hz[CAT(i, j-1)]);

        d_ey[pos] = ca*d_ey[pos] - cb*(d_hz[pos] - d_hz[CAT(i-1, j)]);
    }
  
    /*PEC treatment*/
    if (i>0 && j>0) {
	    if ((mattype_xm!=10 && mattype==10) || (mattype_xm==10 && mattype!=10))
	    {
		    d_ey[pos] = 0;
		    d_ez[pos] = 0;
	    }

	    if ((mattype_ym!=10 && mattype==10) || (mattype_ym==10 && mattype!=10))
	    {
		    d_ex[pos] = 0;
		    d_ez[pos] = 0;
	    }
    }


}

__global__  void
eKernel_none(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         int xres, int yres, int tmmode,
         float dx, float dy, float dt, int dir)
{
    int i, j;
    float cb;
    long int pos;

    i = threadIdx.x;
    j = blockIdx.x;

    if (i<=0 || j<=0) return;
    if (i>=(xres) || j>=(yres)) return;

    pos = CAT(i, j);
    cb  = dt/EPSILON_0/dx;

    if (tmmode) {
        
        d_ez[pos] = d_ez[pos] + cb*((d_hy[pos] - d_hy[CAT(i-1, j)]) -
                            (d_hx[pos] - d_hx[CAT(i, j-1)]));
    }
    else {
        d_ex[pos] = d_ex[pos] + cb*(d_hz[pos] - d_hz[CAT(i, j-1)]);

        d_ey[pos] = d_ey[pos] - cb*(d_hz[pos] - d_hz[CAT(i-1, j)]);
    }

}

/*any linear material given voxel-by-voxel. Note that hKernel_tab (for tabulated magnetic material) is not implemented now*/
__global__  void
hKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_epsilon, float *d_mu, float *d_sigma, float *d_sigast, int xres, int yres, int tmmode,
         float dx, float dy, float dt, int dir)
{
    int i, j;
    float da, db;
    long int pos;

    i = threadIdx.x;
    j = blockIdx.x;

    if (i<0 || j<0) return;
    if (i>=(xres-1) || j>=(yres-1)) return;

    pos = CAT(i, j);

    da  = (1 - d_sigast[pos]*dt/2/d_mu[pos])/(1 + d_sigast[pos]*dt/2/d_mu[pos]);
    db  = (dt/d_mu[pos]/dx)/(1 + d_sigast[pos]*dt/2/d_mu[pos]);

    if (tmmode) {

	    d_hx[pos] = da*d_hx[pos] - db*(d_ez[CAT(i, j+1)] - d_ez[pos]);

	    d_hy[pos] = da*d_hy[pos] + db*(d_ez[CAT(i+1, j)] - d_ez[pos]);
    } else {

	    d_hz[pos] = da*d_hz[pos] + db*((d_ex[CAT(i, j+1)] - d_ex[pos]) -
			    (d_ey[CAT(i+1, j)] - d_ey[pos]));
    }


}

/*vaccum only*/
__global__  void
hKernel_none(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         int xres, int yres, int tmmode,
         float dx, float dy, float dt, int dir)
{
    int i, j;
    float db;
    long int pos;

    i = threadIdx.x;
    j = blockIdx.x;

    if (i<0 || j<0) return;
    if (i>=(xres-1) || j>=(yres-1)) return;

    pos = CAT(i, j);

    db  = dt/MU_0/dx;

    if (tmmode) {

	    d_hx[pos] = d_hx[pos] - db*(d_ez[CAT(i, j+1)] - d_ez[pos]);

	    d_hy[pos] = d_hy[pos] + db*(d_ez[CAT(i+1, j)] - d_ez[pos]);
    } else {

	    d_hz[pos] = d_hz[pos] + db*((d_ex[CAT(i, j+1)] - d_ex[pos]) -
			    (d_ey[CAT(i+1, j)] - d_ey[pos]));
    }
}

/*copy data for absorbing boundary condition*/
__global__  void
liaocpyKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_x0, float *d_xn, float *d_y0, float *d_yn, 
         int xres, int yres, int* bnds,
         float dx, float dy, float dt, int dir)
{
    int i, j;
    long int pos;

    i = threadIdx.x;
    j = blockIdx.x;

    pos = CAT(i, j);

    if (bnds[BNDX0] == 2) {
        if (i==0) {
            d_x0[CLX(EX, j)] = d_ex[pos];
            d_x0[CLX(EY, j)] = d_ey[pos];
            d_x0[CLX(EZ, j)] = d_ez[pos];
            d_x0[CLX(HX, j)] = d_hx[pos];
            d_x0[CLX(HY, j)] = d_hy[pos];
            d_x0[CLX(HZ, j)] = d_hz[pos];
        } else if (i==1) {
            d_x0[CLXP(EX, j)] = d_ex[pos];
            d_x0[CLXP(EY, j)] = d_ey[pos];
            d_x0[CLXP(EZ, j)] = d_ez[pos];
            d_x0[CLXP(HX, j)] = d_hx[pos];
            d_x0[CLXP(HY, j)] = d_hy[pos];
            d_x0[CLXP(HZ, j)] = d_hz[pos];
        }
    }

    if (bnds[BNDXN] == 2) {
        if (i==(xres-1)) {
            d_xn[CLX(EX, j)] = d_ex[pos];
            d_xn[CLX(EY, j)] = d_ey[pos];
            d_xn[CLX(EZ, j)] = d_ez[pos];
            d_xn[CLX(HX, j)] = d_hx[pos];
            d_xn[CLX(HY, j)] = d_hy[pos];
            d_xn[CLX(HZ, j)] = d_hz[pos];
        } else if (i==(xres-2)) {
            d_xn[CLXP(EX, j)] = d_ex[pos];
            d_xn[CLXP(EY, j)] = d_ey[pos];
            d_xn[CLXP(EZ, j)] = d_ez[pos];
            d_xn[CLXP(HX, j)] = d_hx[pos];
            d_xn[CLXP(HY, j)] = d_hy[pos];
            d_xn[CLXP(HZ, j)] = d_hz[pos];
        }
    }

    if (bnds[BNDY0] == 2) {
        if (j==0) {
            d_y0[CLY(EX, i)] = d_ex[pos];
            d_y0[CLY(EY, i)] = d_ey[pos];
            d_y0[CLY(EZ, i)] = d_ez[pos];
            d_y0[CLY(HX, i)] = d_hx[pos];
            d_y0[CLY(HY, i)] = d_hy[pos];
            d_y0[CLY(HZ, i)] = d_hz[pos];
        } else if (j==1) {
            d_y0[CLYP(EX, i)] = d_ex[pos];
            d_y0[CLYP(EY, i)] = d_ey[pos];
            d_y0[CLYP(EZ, i)] = d_ez[pos];
            d_y0[CLYP(HX, i)] = d_hx[pos];
            d_y0[CLYP(HY, i)] = d_hy[pos];
            d_y0[CLYP(HZ, i)] = d_hz[pos];
        }
    }

    if (bnds[BNDYN] == 2) {
        if (j==(yres-1)) {
            d_yn[CLY(EX, i)] = d_ex[pos];
            d_yn[CLY(EY, i)] = d_ey[pos];
            d_yn[CLY(EZ, i)] = d_ez[pos];
            d_yn[CLY(HX, i)] = d_hx[pos];
            d_yn[CLY(HY, i)] = d_hy[pos];
            d_yn[CLY(HZ, i)] = d_hz[pos];
        } else if (j==(yres-2)) {
            d_yn[CLYP(EX, i)] = d_ex[pos];
            d_yn[CLYP(EY, i)] = d_ey[pos];
            d_yn[CLYP(EZ, i)] = d_ez[pos];
            d_yn[CLYP(HX, i)] = d_hx[pos];
            d_yn[CLYP(HY, i)] = d_hy[pos];
            d_yn[CLYP(HZ, i)] = d_hz[pos];
        }
    }
}


/*run the absorbing boundary conditon*/
__global__  void
liaorunKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz, float *d_epsilon,
         float *d_x0, float *d_xn, float *d_y0, float *d_yn, 
         int xres, int yres, int* bnds, int matmode,
         float dx, float dy, float dt, int dir)
{
    float lssx, lssy, ind;
    float mx, my;
    int i, j;
    long int pos;

    i = threadIdx.x;
    j = blockIdx.x;

    pos = CAT(i, j);

    if (matmode==1 || matmode == 3)
        ind = sqrt(d_epsilon[pos]/EPSILON_0);
    else ind = 1;

    lssx = (dt*LS/ind-dx)/(dt*LS/ind+dx);
    lssy = (dt*LS/ind-dy)/(dt*LS/ind+dy);
    mx = (float)lssx;
    my = (float)lssy;

    if (bnds[BNDX0] == 2) {
        if (i==0) {
            d_ex[pos] = d_x0[CLXP(EX, j)] + mx*(d_ex[CAT(i+1, j)] - d_x0[CLX(EX, j)]);
            d_ey[pos] = d_x0[CLXP(EY, j)] + mx*(d_ey[CAT(i+1, j)] - d_x0[CLX(EY, j)]);
            d_ez[pos] = d_x0[CLXP(EZ, j)] + mx*(d_ez[CAT(i+1, j)] - d_x0[CLX(EZ, j)]);
        }
    }
    if (bnds[BNDXN] == 2) {
        if (i==(xres-1)) {
            d_ex[pos] = d_xn[CLXP(EX, j)] + mx*(d_ex[CAT(i-1, j)] - d_xn[CLX(EX, j)]);
            d_ey[pos] = d_xn[CLXP(EY, j)] + mx*(d_ey[CAT(i-1, j)] - d_xn[CLX(EY, j)]);
            d_ez[pos] = d_xn[CLXP(EZ, j)] + mx*(d_ez[CAT(i-1, j)] - d_xn[CLX(EZ, j)]);
            d_hx[pos] = d_xn[CLXP(HX, j)] + mx*(d_hx[CAT(i-1, j)] - d_xn[CLX(HX, j)]);
            d_hy[pos] = d_xn[CLXP(HY, j)] + mx*(d_hy[CAT(i-1, j)] - d_xn[CLX(HY, j)]);
            d_hz[pos] = d_xn[CLXP(HZ, j)] + mx*(d_hz[CAT(i-1, j)] - d_xn[CLX(HZ, j)]);
        }
    }

    if (bnds[BNDY0] == 2) {
        if (j==0) {
            d_ex[pos] = d_y0[CLYP(EX, i)] + my*(d_ex[CAT(i, j+1)] - d_y0[CLY(EX, i)]);
            d_ey[pos] = d_y0[CLYP(EY, i)] + my*(d_ey[CAT(i, j+1)] - d_y0[CLY(EY, i)]);
            d_ez[pos] = d_y0[CLYP(EZ, i)] + my*(d_ez[CAT(i, j+1)] - d_y0[CLY(EZ, i)]);
        }
    }
    if (bnds[BNDYN] == 2) {
        if (j==(yres-1)) {
            d_ex[pos] = d_yn[CLYP(EX, i)] + my*(d_ex[CAT(i, j-1)] - d_yn[CLY(EX, i)]);
            d_ey[pos] = d_yn[CLYP(EY, i)] + my*(d_ey[CAT(i, j-1)] - d_yn[CLY(EY, i)]);
            d_ez[pos] = d_yn[CLYP(EZ, i)] + my*(d_ez[CAT(i, j-1)] - d_yn[CLY(EZ, i)]);
            d_hx[pos] = d_yn[CLYP(HX, i)] + my*(d_hx[CAT(i, j-1)] - d_yn[CLY(HX, i)]);
            d_hy[pos] = d_yn[CLYP(HY, i)] + my*(d_hy[CAT(i, j-1)] - d_yn[CLY(HY, i)]);
            d_hz[pos] = d_yn[CLYP(HZ, i)] + my*(d_hz[CAT(i, j-1)] - d_yn[CLY(HZ, i)]);
        }
    }
}


/*electric field point source*/
__global__  void
srcepointKernel(float *d_ex, float *d_ey, float *d_ez, int xres, int yres, int ipos, int jpos,
         float ex, float ey, float ez, int dir)
{
    int i, j;
    long int pos;

    i = threadIdx.x;
    j = blockIdx.x;

    if (i==ipos && j==jpos) {  
       pos = CAT(i, j);
       if (ex != 0) d_ex[pos] = ex;
       if (ey != 0) d_ey[pos] = ey;
       if (ez != 0) d_ez[pos] = ez;
    }

}

/*magnetic field point source*/
__global__  void
srchpointKernel(float *d_hx, float *d_hy, float *d_hz, int xres, int yres, int ipos, int jpos,
         float hx, float hy, float hz, int dir)
{
    int i, j;
    long int pos;

    i = threadIdx.x;
    j = blockIdx.x;

    if (i==ipos && j==jpos) {  
       pos = CAT(i, j);
       if (hx != 0) d_hx[pos] = hx;
       if (hy != 0) d_hy[pos] = hy;
       if (hz != 0) d_hz[pos] = hz;
    }
}

/*set of functions for TSF source*/
/*
__inline__ __device__  float
k_dist(float x1, float y1, float z1, float x2, float y2, float z2)
{
    return sqrtf((x1 - x2)*(x1 - x2)
                + (y1 - y2)*(y1 - y2)
                + (z1 - z2)*(z1 - z2));
}

__inline__ __device__ float 
k_angvec(float x1, float y1, float z1, float x2, float y2, float z2)
{
    return ((x1*x2 + y1*y2 + z1*z2)/(sqrtf((float)(x1*x1 + y1*y1 + z1*z1)*(x2*x2 + y2*y2 + z2*z2))));
}

__device__ float 
k_angle(float x1, float y1, float z1, float x2, float y2, float z2, int side) 
{
    if (side == 0) return k_angvec(x2 - x1, y2 - y1, z2 - z1, 0, -1, 0);
    else if (side == 1) return k_angvec(x2 - x1, y2 - y1, z2 - z1, 0, 1, 0);
    else if (side == 2) return k_angvec(x2 - x1, y2 - y1, z2 - z1, -1, 0, 0);
    else if (side == 3) return k_angvec(x2 - x1, y2 - y1, z2 - z1, 1, 0, 0);
    else if (side == 4) return k_angvec(x2 - x1, y2 - y1, z2 - z1, 0, 0, -1);
    else return k_angvec(x2 - x1, y2 - y1, z2 - z1, 0, 0, 1);
}


*/

/*
__global__  void
tsfjstepKernel(float *d_jpool_e, float *d_jpool_h, float *d_jpvals,
            float *d_jpool_epsilon, float *d_jpool_mu,
            float *d_jpool_sigma, float *d_jpool_sigast, 
            float dx, float dt, float e, int n)
{
        int i;
        //copybound
        float de0, de1, dh0, dh1, den, dem, dhn, dhm;
        float dAngleMult=1;
        float dMulth=dt/(MU_0*dx)/dAngleMult;
        float dMulte=dt/(EPSILON_0*dx)/dAngleMult;         

        de0 = d_jpvals[0];
        de1 = d_jpvals[1];
        dh0 = d_jpvals[2];
        dh1 = d_jpvals[3];
        den = d_jpvals[4];
        dem = d_jpvals[5];
        dhn = d_jpvals[6];
        dhm = d_jpvals[7];

        //ystep_e
        for (i=1; i<(n); i++)
           d_jpool_e[i] += 1.0/(d_jpool_epsilon[i])*dMulte*(d_jpool_h[i-1] - d_jpool_h[i]);

        //bound
        d_jpool_e[0] = de1 + (dt*LS-dx)/(dt*LS+dx)*(d_jpool_e[1] - de0);
        d_jpool_h[0] = dh1 + (dt*LS-dx)/(dt*LS+dx)*(d_jpool_h[1] - dh0);
        d_jpool_e[n-1] = dem + (dt*LS-dx)/(dt*LS+dx)*(d_jpool_e[n-2] - den);
        d_jpool_h[n-1] = dhm + (dt*LS-dx)/(dt*LS+dx)*(d_jpool_h[n-2] - dhn);

        //copybound
        d_jpvals[0] = d_jpool_e[0];
        d_jpvals[1] = d_jpool_e[1];
        d_jpvals[2] = d_jpool_h[0];
        d_jpvals[3] = d_jpool_h[1];
        d_jpvals[4] = d_jpool_e[n-1];
        d_jpvals[5] = d_jpool_e[n-2];
        d_jpvals[6] = d_jpool_h[n-1];
        d_jpvals[7] = d_jpool_h[n-2];
      
        //applysource
        d_jpool_e[0] = e;

        //ystep_h
        for (i=0; i<(n-1); i++)
           d_jpool_h[i] += 1.0/(d_jpool_mu[i])*dMulth*(d_jpool_e[i] - d_jpool_e[i+1]);  

}
*/

/*set of functions for */ 
/*
__global__  void
mbnxeKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
          int xres, int yres, int zres, int mb_bx0, int mb_bxn, 
          int mb_bx0pos, int mb_bxnpos, int mb_by0pos, int mb_bynpos, int mb_bz0pos, int mb_bznpos, int dir)
{
    int i, j, k, pos1, pos2;
    if (dir==2) {
       i = blockIdx.x;
       j = blockIdx.y;
       k = threadIdx.x;  
    } else if (dir==0) {
       j = blockIdx.x;
       k = blockIdx.y;
       i = threadIdx.x;  
    } else {
       i = blockIdx.x;
       k = blockIdx.y;
       j = threadIdx.x;  
    }

    if ((mb_bx0 == 4 || mb_bxn == 4) && (i==mb_bx0pos && j>=mb_by0pos && j<mb_bynpos && k>=mb_bz0pos && k<mb_bznpos)) {
        pos1 = CAT(i, j, k);
        pos2 = CAT(mb_bxnpos, j, k);
        d_ey[pos2] = d_ey[pos1];
        d_ez[pos2] = d_ez[pos1];
     }
}
*/

/*
__global__  void
mbnyeKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
          int xres, int yres, int zres, int mb_by0, int mb_byn, 
          int mb_bx0pos, int mb_bxnpos, int mb_by0pos, int mb_bynpos, int mb_bz0pos, int mb_bznpos, int dir)
{
    int i, j, k, pos1, pos2;
    if (dir==2) {
       i = blockIdx.x;
       j = blockIdx.y;
       k = threadIdx.x;  
    } else if (dir==0) {
       j = blockIdx.x;
       k = blockIdx.y;
       i = threadIdx.x;  
    } else {
       i = blockIdx.x;
       k = blockIdx.y;
       j = threadIdx.x;  
    }

    if ((mb_by0 == 4 || mb_byn == 4) && (j==mb_by0pos && i>=mb_bx0pos && i<mb_bxnpos && k>=mb_bz0pos && k<mb_bznpos)) {
        pos1 = CAT(i, j, k);
        pos2 = CAT(i, mb_bynpos, k);
        d_ex[pos2] = d_ex[pos1];
        d_ez[pos2] = d_ez[pos1];
    }
}
*/

/*
__global__  void
mbnxhKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
          int xres, int yres, int zres, int mb_bx0, int mb_bxn, 
          int mb_bx0pos, int mb_bxnpos, int mb_by0pos, int mb_bynpos, int mb_bz0pos, int mb_bznpos, int dir)
{
    int i, j, k, pos1, pos2;
    if (dir==2) {
       i = blockIdx.x;
       j = blockIdx.y;
       k = threadIdx.x;  
    } else if (dir==0) {
       j = blockIdx.x;
       k = blockIdx.y;
       i = threadIdx.x;  
    } else {
       i = blockIdx.x;
       k = blockIdx.y;
       j = threadIdx.x;  
    }

    if ((mb_bx0 == 4 || mb_bxn == 4) && (i==mb_bx0pos && j>=mb_by0pos && j<mb_bynpos && k>=mb_bz0pos && k<mb_bznpos)) {
        pos1 = CAT(i-1, j, k);
        pos2 = CAT(mb_bxnpos-1, j, k);
        d_hy[pos1] = d_hy[pos2];
        d_hz[pos1] = d_hz[pos2];
     }
}
*/

/*
__global__  void
mbnyhKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
          int xres, int yres, int zres, int mb_by0, int mb_byn, 
          int mb_bx0pos, int mb_bxnpos, int mb_by0pos, int mb_bynpos, int mb_bz0pos, int mb_bznpos, int dir)
{
    int i, j, k, pos1, pos2;
    if (dir==2) {
       i = blockIdx.x;
       j = blockIdx.y;
       k = threadIdx.x;  
    } else if (dir==0) {
       j = blockIdx.x;
       k = blockIdx.y;
       i = threadIdx.x;  
    } else {
       i = blockIdx.x;
       k = blockIdx.y;
       j = threadIdx.x;  
    }

    if ((mb_by0 == 4 || mb_byn == 4) && (j==mb_by0pos && i>=mb_bx0pos && i<mb_bxnpos && k>=mb_bz0pos && k<mb_bznpos)) {
        pos1 = CAT(i, j-1, k);
        pos2 = CAT(i, mb_bynpos-1, k);
        d_hx[pos1] = d_hx[pos2];
        d_hz[pos1] = d_hz[pos2];
     }

}
*/

/*

__device__ float
k_dcomp(int i, int j, int k, int xres, int yres, int zres, float theta, float phi,
      int i0, int i1, int j0, int j1, int k0, int k1)
{
        float rx, ry, rz;
        float ax, ay, az;

        ax = sin(theta)*cos(phi);
        ay = sin(theta)*sin(phi);
        az = cos(theta);
        if (theta >= 0 && theta <= (CU_PI/2.0)) {
            if (phi >= 0 && phi <= (CU_PI/2.0)) {
                rx = i - i0;
                ry = j - j0;
                rz = k - k0;
            }
            else if (phi > (CU_PI/2.0) && phi <= CU_PI) {
                rx = i - i1;
                ry = j - j0;
                rz = k - k0;
            }
            else if (phi > CU_PI && phi <= (3.0*CU_PI/2.0)) {
                rx = i - i1;
                ry = j - j1;
                rz = k - k0;
            }
            else {
                rx = i - i0;
                ry = j - j1;
                rz = k - k0;
            }
        }
        else if (theta < CU_PI && theta > (CU_PI/2.0)) 
        {
            if (phi >= 0 && phi <= (CU_PI/2.0)) {
                rx = i - i0;
                ry = j - j0;
                rz = k - k1;
            }
            else if (phi > (CU_PI/2.0) && phi <= CU_PI) {
                rx = i - i1;
                ry = j - j0;
                rz = k - k1;
            }
            else if (phi > CU_PI && phi <= (3.0*CU_PI/2.0)) {
                rx = i - i1;
                ry = j - j1;
                rz = k - k1;
            }
            else {
                rx = i - i0;
                ry = j - j1;
                rz = k - k1;
            }
         }
    return 10 + (ax*rx + ay*ry + az*rz);
}
*/

/*
__device__ float
k_rdcomp(int i, int j, int k, int xres, int yres, int zres, float theta, float phi,
      int i0, int i1, int j0, int j1, int k0, int k1)
{
        float rx, ry, rz;
        float ax, ay, az;

        ax = sin(theta)*cos(phi);
        ay = sin(theta)*sin(phi);
        az = cos(theta);

        rx = -i+xres/2;
        ry = -j+yres/2;
        rz = -k+zres/2;

        if ((ax*rx + ay*ry + az*rz)<250)
        return 250 - (ax*rx + ay*ry + az*rz);
        else return 0;
}
*/

/*
__device__ float k_gex(float field, float theta, float phi, float psi)
{
    return field*(cos(psi)*sin(phi) - sin(psi)*cos(theta)*cos(phi));
}
__device__ float k_gey(float field, float theta, float phi, float psi)
{
    return field*(-cos(psi)*cos(phi) - sin(psi)*cos(theta)*sin(phi));
}
__device__ float k_gez(float field, float theta, float phi, float psi)
{
    return field*(sin(psi)*sin(theta));
}
__device__ float k_ghx(float field, float theta, float phi, float psi)
{
    return field*(sin(psi)*sin(phi) + cos(psi)*cos(theta)*cos(phi));
}
__device__ float k_ghy(float field, float theta, float phi, float psi)
{
    return field*(-sin(psi)*cos(phi) + cos(psi)*cos(theta)*sin(phi));
}
__device__ float k_ghz(float field, float theta, float phi, float psi)
{
    return field*(-cos(psi)*sin(theta));
}
*/

/*
__device__ float
get_dval(float *line, int res, float x)
{
    float w1, w2, w3, w4;
    int l = (int)(x);
    float a = x-(float)l;

    if (x>=1 && x<(res-1))
    {
        w1=a+1; w2=a; w3=1-a; w4=2-a;
        w1=4-8*w1+5*w1*w1-w1*w1*w1;
        w2=1-2*w2*w2+w2*w2*w2;
        w3=1-2*w3*w3+w3*w3*w3;
        w4=4-8*w4+5*w4*w4-w4*w4*w4;
        return w1*line[l-1]+w2*line[l]+w3*line[l+1]+w4*line[l+2];
    } else if ((x<1 && x>=0) || x>=(res-1) && x<(res))
    {
        return (1-a)*line[l]+a*line[l+1];
    }
    else return 0;
}
*/
/*

__global__  void
tsf_estep_aKernel(float *d_ex, float *d_ey, float *d_ez,
          float *d_hx, float *d_hy, float *d_hz,
          float *d_epsilon, float *d_mu, int matmode,
          float *d_tsfset,
          float *d_jpool_e, float *d_jpool_h, int jpool_size,
          float dx, float dt,
          int xres, int yres, int zres, int dir)
{
    int i, j, k;
    long int pos;
    float d, vh;  
    int i0 = (int)d_tsfset[0];
    int j0 = (int)d_tsfset[1];
    int k0 = (int)d_tsfset[2];
    int i1 = (int)d_tsfset[3];
    int j1 = (int)d_tsfset[4];
    int k1 = (int)d_tsfset[5];
    float theta = d_tsfset[6];
    float phi = d_tsfset[7];
    float psi = d_tsfset[8];
    float corr = d_tsfset[9];
    float epsilon;
    float gcorr = 1;
    int skip_i0 = (int)d_tsfset[10];
    int skip_in = (int)d_tsfset[11];
    int gaussian = (int)d_tsfset[16];
    float gaussian_fxpos = d_tsfset[17];
    float gaussian_fypos = d_tsfset[18];
    float gaussian_rx = d_tsfset[19];
    float gaussian_ry = d_tsfset[20];

    if (dir==2) {
       i = blockIdx.x;
       j = blockIdx.y;
       k = threadIdx.x;  
    } else if (dir==0) {
       j = blockIdx.x;
       k = blockIdx.y;
       i = threadIdx.x;  
    } else {
       i = blockIdx.x;
       k = blockIdx.y;
       j = threadIdx.x;  
    }

    pos = CAT(i, j, k);
    if (matmode == 0 || matmode == 1)
        epsilon = d_epsilon[pos];
    else epsilon = EPSILON_0;


    if ((!skip_i0) && i==(i0))
    {
        if (j>=j0 && j<=j1 && k>=k0 && k<=k1)
        {
             if (gaussian) gcorr*=gaussmult(i0-1, j, k, gaussian_fxpos, gaussian_fypos, gaussian_rx, gaussian_ry);
             d = corr*k_dcomp(i0-1, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             vh = gcorr*get_dval(d_jpool_h, jpool_size, d);

             if (j<j1) d_ey[pos] += dt/dx/epsilon*
                     k_ghz(vh, theta, phi, psi);                                         
             if (k<k1) d_ez[pos] -= dt/dx/epsilon*
                     k_ghy(vh, theta, phi, psi);                                         
        }
    }
    if ((!skip_in) && i==i1)
    {
        if (j>=j0 && j<=j1 && k>=k0 && k<=k1)
        {
             if (gaussian) gcorr*=gaussmult(i1, j, k, gaussian_fxpos, gaussian_fypos, gaussian_rx, gaussian_ry);
             d = corr*k_dcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             vh = gcorr*get_dval(d_jpool_h, jpool_size, d);

             if (j<j1) d_ey[pos] -= dt/dx/epsilon*
                     k_ghz(vh, theta, phi, psi);                                         
             if (k<k1) d_ez[pos] += dt/dx/epsilon*
                     k_ghy(vh, theta, phi, psi);                                         
        }
    }
}
*/

/*
__global__  void
tsf_estep_bKernel(float *d_ex, float *d_ey, float *d_ez,
          float *d_hx, float *d_hy, float *d_hz,
          float *d_epsilon, float *d_mu, int matmode,
          float *d_tsfset,
          float *d_jpool_e, float *d_jpool_h, int jpool_size,
          float dx, float dt,
          int xres, int yres, int zres, int dir)
{
    int i, j, k;
    long int pos;
    float d, vh, epsilon;    
    int i0 = (int)d_tsfset[0];
    int j0 = (int)d_tsfset[1];
    int k0 = (int)d_tsfset[2];
    int i1 = (int)d_tsfset[3];
    int j1 = (int)d_tsfset[4];
    int k1 = (int)d_tsfset[5];
    float theta = d_tsfset[6];
    float phi = d_tsfset[7];
    float psi = d_tsfset[8];
    float corr = d_tsfset[9];
    float gcorr = 1;
    int skip_j0 = (int)d_tsfset[12];
    int skip_jn = (int)d_tsfset[13];
    int gaussian = (int)d_tsfset[16];
    float gaussian_fxpos = d_tsfset[17];
    float gaussian_fypos = d_tsfset[18];
    float gaussian_rx = d_tsfset[19];
    float gaussian_ry = d_tsfset[20];

    if (dir==2) {
       i = blockIdx.x;
       j = blockIdx.y;
       k = threadIdx.x;  
    } else if (dir==0) {
       j = blockIdx.x;
       k = blockIdx.y;
       i = threadIdx.x;  
    } else {
       i = blockIdx.x;
       k = blockIdx.y;
       j = threadIdx.x;  
    }

    pos = CAT(i, j, k);

    if (matmode == 0 || matmode == 1)
        epsilon = d_epsilon[pos];
    else epsilon = EPSILON_0;

    if ((!skip_j0) && j==(j0))
    {
        if (i>=i0 && i<=i1 && k>=k0 && k<=k1)
        {
             if (gaussian) gcorr*=gaussmult(i, j0-1, k, gaussian_fxpos, gaussian_fypos, gaussian_rx, gaussian_ry);
             d = corr*k_dcomp(i, j0-1, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             vh = gcorr*get_dval(d_jpool_h, jpool_size, d);

             if (i<i1) d_ex[pos] -= dt/dx/epsilon*
                     k_ghz(vh, theta, phi, psi);                                         
             if (k<k1) d_ez[pos] += dt/dx/epsilon*
                     k_ghx(vh, theta, phi, psi);                                         
        }
    }
    if ((!skip_jn) && j==j1)
    {
        if (i>=i0 && i<=i1 && k>=k0 && k<=k1)
        {
             if (gaussian) gcorr*=gaussmult(i, j1, k, gaussian_fxpos, gaussian_fypos, gaussian_rx, gaussian_ry);
             d = corr*k_dcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             vh = gcorr*get_dval(d_jpool_h, jpool_size, d);

             if (i<i1) d_ex[pos] += dt/dx/epsilon*
                     k_ghz(vh, theta, phi, psi);                                         
             if (k<k1) d_ez[pos] -= dt/dx/epsilon*
                     k_ghx(vh, theta, phi, psi);                                         
        }
    }

}
*/

/*
__global__  void
tsf_hstep_aKernel(float *d_ex, float *d_ey, float *d_ez,
          float *d_hx, float *d_hy, float *d_hz,
          float *d_epsilon, float *d_mu, int matmode,
          float *d_tsfset,
          float *d_jpool_e, float *d_jpool_h, int jpool_size,
          float dx, float dt,
          int xres, int yres, int zres, int dir)
{
    int i, j, k;
    long int pos;
    float d, ve, mu;  
    int i0 = (int)d_tsfset[0];
    int j0 = (int)d_tsfset[1];
    int k0 = (int)d_tsfset[2];
    int i1 = (int)d_tsfset[3];
    int j1 = (int)d_tsfset[4];
    int k1 = (int)d_tsfset[5];
    float theta = d_tsfset[6];
    float phi = d_tsfset[7];
    float psi = d_tsfset[8];
    float corr = d_tsfset[9];
    float gcorr = 1;
    int skip_i0 = (int)d_tsfset[10];
    int skip_in = (int)d_tsfset[11];
    int gaussian = (int)d_tsfset[16];
    float gaussian_fxpos = d_tsfset[17];
    float gaussian_fypos = d_tsfset[18];
    float gaussian_rx = d_tsfset[19];
    float gaussian_ry = d_tsfset[20];




    if (dir==2) {
       i = blockIdx.x;
       j = blockIdx.y;
       k = threadIdx.x;  
    } else if (dir==0) {
       j = blockIdx.x;
       k = blockIdx.y;
       i = threadIdx.x;  
    } else {
       i = blockIdx.x;
       k = blockIdx.y;
       j = threadIdx.x;  
    }

    pos = CAT(i, j, k);
    if (matmode == 0 || matmode == 2)
        mu = d_mu[pos];
    else mu = MU_0;


    if ((!skip_i0) && i==i0)
    {
        if (j>=j0 && j<=j1 && k>=k0 && k<=k1)
        {
             if (gaussian) gcorr*=gaussmult(i, j, k, gaussian_fxpos, gaussian_fypos, gaussian_rx, gaussian_ry);
             d = corr*k_dcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
             
             ve = gcorr*get_dval(d_jpool_e, jpool_size, d);

             if (k<k1) d_hy[CAT(i-1, j, k)] -= dt/dx/mu*
                     k_gez(ve, theta, phi, psi);                                         
             if (j<j1) d_hz[CAT(i-1, j, k)] += dt/dx/mu*
                     k_gey(ve, theta, phi, psi);            
        }
    }

    if ((!skip_in) && i==i1)
    {
        if (j>=j0 && j<=j1 && k>=k0 && k<=k1)
        {
             if (gaussian) gcorr*=gaussmult(i, j, k, gaussian_fxpos, gaussian_fypos, gaussian_rx, gaussian_ry);
             d = corr*k_dcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             ve = gcorr*get_dval(d_jpool_e, jpool_size, d);

             if (k<k1) d_hy[pos] += dt/dx/mu*
                     k_gez(ve, theta, phi, psi);                                         
             if (j<j1) d_hz[pos] -= dt/dx/mu*
                     k_gey(ve, theta, phi, psi);            
        }
    }

}
*/
/*
__global__  void
tsf_hstep_bKernel(float *d_ex, float *d_ey, float *d_ez,
          float *d_hx, float *d_hy, float *d_hz, 
          float *d_epsilon, float *d_mu, int matmode,
          float *d_tsfset,
          float *d_jpool_e, float *d_jpool_h, int jpool_size,
          float dx, float dt,
          int xres, int yres, int zres, int dir)
{
    int i, j, k;
    long int pos;
    float d, ve, mu;    
    int i0 = (int)d_tsfset[0];
    int j0 = (int)d_tsfset[1];
    int k0 = (int)d_tsfset[2];
    int i1 = (int)d_tsfset[3];
    int j1 = (int)d_tsfset[4];
    int k1 = (int)d_tsfset[5];
    float theta = d_tsfset[6];
    float phi = d_tsfset[7];
    float psi = d_tsfset[8];
    float corr = d_tsfset[9];
    float gcorr = 1;
    int skip_j0 = (int)d_tsfset[12];
    int skip_jn = (int)d_tsfset[13];
    int gaussian = (int)d_tsfset[16];
    float gaussian_fxpos = d_tsfset[17];
    float gaussian_fypos = d_tsfset[18];
    float gaussian_rx = d_tsfset[19];
    float gaussian_ry = d_tsfset[20];



    if (dir==2) {
       i = blockIdx.x;
       j = blockIdx.y;
       k = threadIdx.x;  
    } else if (dir==0) {
       j = blockIdx.x;
       k = blockIdx.y;
       i = threadIdx.x;  
    } else {
       i = blockIdx.x;
       k = blockIdx.y;
       j = threadIdx.x;  
    }

    pos = CAT(i, j, k);
    if (matmode == 0 || matmode == 2)
        mu = d_mu[pos];
    else mu = MU_0;

 

    if ((!skip_j0) && j==j0)
    {
        if (i>=i0 && i<=i1 && k>=k0 && k<=k1)
        {
             if (gaussian) gcorr*=gaussmult(i, j, k, gaussian_fxpos, gaussian_fypos, gaussian_rx, gaussian_ry);
             d = corr*k_dcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             ve = gcorr*get_dval(d_jpool_e, jpool_size, d);

             if (k<k1) d_hx[CAT(i, j-1, k)] += dt/dx/mu*
                     k_gez(ve, theta, phi, psi);                                         
             if (i<i1) d_hz[CAT(i, j-1, k)] -= dt/dx/mu*
                     k_gex(ve, theta, phi, psi);            
        }
    }
    if ((!skip_jn) && j==j1)
    {
        if (i>=i0 && i<=i1 && k>=k0 && k<=k1)
        {
             if (gaussian) gcorr*=gaussmult(i, j, k, gaussian_fxpos, gaussian_fypos, gaussian_rx, gaussian_ry);
             d = corr*k_dcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             ve = gcorr*get_dval(d_jpool_e, jpool_size, d);

             if (k<k1) d_hx[pos] -= dt/dx/mu*
                     k_gez(ve, theta, phi, psi);                                         
             if (i<i1) d_hz[pos] += dt/dx/mu*
                     k_gex(ve, theta, phi, psi);            
        }
    }

}
*/


/*
__device__ void 
get_matprops(int i, int j, int k, float *epsilon, float *mu,
             int *d_mat, int *d_mattype, float *d_mattab, int matmode, int nmat,
             float *d_epsilon, float *d_mu, int xres, int yres, int zres)
{
	long int pos = CAT(i, j, k);
	int mattype = d_mattype[d_mat[pos]];  //0..normal, >1..tabulated

	if (matmode==0) {
		*epsilon = EPSILON_0;
		*mu = MU_0;
	} 
	else if (matmode == 1 || matmode==2 || matmode==3) { //= mattype 0 = linear material given pixel by pixel
		if (nmat==0 || d_mat[pos] == 0) {
			if (matmode == 2) *epsilon = EPSILON_0;
			else *epsilon = d_epsilon[pos];
			if (matmode == 1) *mu = MU_0;
			*mu = d_mu[pos];
		} else { //some tabulated material inside linear material
			*epsilon = d_mattab[GETEPS(d_mat[pos])]; 
			*mu = d_mattab[GETMU(d_mat[pos])];
		}
	} else if (nmat>0 && mattype==1) { //tabulated linear material
		*epsilon = d_mattab[GETEPS(d_mat[pos])]; 
		*mu = d_mattab[GETMU(d_mat[pos])];
	} else {
		*epsilon = EPSILON_0;
		*mu = MU_0;
	}

}
*/

/*
__global__  void
fastsumKernel_a(float *d_ex, float *d_ey, float *d_ez, int xres, int yres, int zres, 
          float *d_epsilon, float *d_mu, float *d_sum_epsilon, float *d_sum_sigma, int *d_sum_mode, 
          int *i0, int *i1, int *j0, int *j1, int *k0, int *k1,
          int nsums, float *d_acc, int dir,
          int *d_mat, int *d_mattype, float *d_mattab, int matmode, int nmat)
{
 
    int m, i, j, k, pos;
    float val, epsilon, mu;

    if (dir==2) {
       i = blockIdx.x;
       j = blockIdx.y;
       m = threadIdx.x;  

       d_acc[m*xres*yres + i*xres + j] = 0;

       if (i<i0[m] || i>=i1[m] || j<j0[m] || j>=j1[m]) return;

       
       for (k=k0[m]; k<k1[m]; k++)
       {
           pos = CAT(i, j, k);

           get_matprops(i, j, k, &epsilon, &mu, 
				    d_mat, d_mattype, d_mattab, matmode, nmat, d_epsilon, d_mu, xres, yres, zres);

            if (fabs(epsilon-d_sum_epsilon[m])<(EPSILON_0/100.0))
           {
               val = 0;
               if (d_sum_mode[m] == 0 || d_sum_mode[m] == 3 || d_sum_mode[m] == 4)
                   val += d_ex[pos]*d_ex[pos];
               if (d_sum_mode[m] == 1 || d_sum_mode[m] == 3 || d_sum_mode[m] == 4)
                   val += d_ey[pos]*d_ey[pos];
               if (d_sum_mode[m] == 2 || d_sum_mode[m] == 3 || d_sum_mode[m] == 4)
                   val += d_ez[pos]*d_ez[pos];

               if (d_sum_mode[m] == 4) d_acc[m*xres*yres + i*xres + j] += d_sum_sigma[m]*val;
               else d_acc[m*xres*yres + i*xres + j] += val;
           }
       }

    } else if (dir==0) {
       j = blockIdx.x;
       k = blockIdx.y;
       i = threadIdx.x;  
    } else {
       i = blockIdx.x;
       k = blockIdx.y;
       j = threadIdx.x;  
    }

}
*/

/*
__global__  void
fastsumKernel_b(float *d_acc, float *d_sums, int nsums, int nsteps, int step, int xres, int yres, int zres, int dir)
{
    int i, j;

    int m = threadIdx.x;
    
    if (dir==2) {
       for (i=0; i<xres; i++)
       {
          for (j=0; j<yres; j++)
          {
              d_sums[m*nsteps + step] += d_acc[m*xres*yres + i*xres + j];
          }
       }
    }
}
*/
  
/* 
__global__  void
sumKernel(float *d_ex, float *d_ey, float *d_ez, int xres, int yres, int zres, 
          float *d_epsilon, float *d_mu, float *d_sums, float *d_sum_epsilon, float *d_sum_sigma, int *d_sum_mode, 
          int *i0, int *i1, int *j0, int *j1, int *k0, int *k1,
          int nsteps, int step,
          int *d_mat, int *d_mattype, float *d_mattab, int matmode, int nmat)
{
    int i, j, k, pos;
    int m = threadIdx.x;
    float val, epsilon, mu;


    for (i=i0[m]; i<i1[m]; i++)
    {
        for (j=j0[m]; j<j1[m]; j++)
        {
            for (k=k0[m]; k<k1[m]; k++)
            {
                pos = CAT(i, j, k);
                val = 0;
    
                get_matprops(i, j, k, &epsilon, &mu, 
				    d_mat, d_mattype, d_mattab, matmode, nmat, d_epsilon, d_mu, xres, yres, zres);

                if (fabs(epsilon-d_sum_epsilon[m])<(EPSILON_0/100.0))
                {
                    if (d_sum_mode[m] == 0 || d_sum_mode[m] == 3 || d_sum_mode[m] == 4)
                        val += d_ex[pos]*d_ex[pos];
                    else if (d_sum_mode[m] == 1 || d_sum_mode[m] == 3 || d_sum_mode[m] == 4)
                        val += d_ey[pos]*d_ey[pos];
                    else if (d_sum_mode[m] == 2 || d_sum_mode[m] == 3 || d_sum_mode[m] == 4)
                        val += d_ez[pos]*d_ez[pos];

                    if (d_sum_mode[m] == 4) d_sums[m*nsteps + step] += d_sum_sigma[m]*val;
                    else d_sums[m*nsteps + step] += val;
                }
            }
        }
    }
}
*/




__global__  void
outpointKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
                              float *d_outpointdata, int *d_outpoint_pos, int nsteps, int noutpoints,
                              int step, int xres, int yres)
{
    int i;
    for (i=0; i<noutpoints; i++) 
    {
        d_outpointdata[6*(i*nsteps + step)] = d_ex[CAT(d_outpoint_pos[3*i], d_outpoint_pos[3*i+1])];
        d_outpointdata[6*(i*nsteps + step) + 1] = d_ey[CAT(d_outpoint_pos[3*i], d_outpoint_pos[3*i+1])];
        d_outpointdata[6*(i*nsteps + step) + 2] = d_ez[CAT(d_outpoint_pos[3*i], d_outpoint_pos[3*i+1])];
        d_outpointdata[6*(i*nsteps + step) + 3] = d_hx[CAT(d_outpoint_pos[3*i], d_outpoint_pos[3*i+1])];
        d_outpointdata[6*(i*nsteps + step) + 4] = d_hy[CAT(d_outpoint_pos[3*i], d_outpoint_pos[3*i+1])];
        d_outpointdata[6*(i*nsteps + step) + 5] = d_hz[CAT(d_outpoint_pos[3*i], d_outpoint_pos[3*i+1])];
    }
}



/*All the kernel wrappers called from rest of the code are here below*/
 
cudaError_t wrap_hKernel(SvGpuPlan *plan)
{
        cudaError_t err;

        int ga, b;

	b = plan->xres; //FIXME: do this dependent on problem size and graphics card computation capability
	ga = plan->yres;

	dim3 dimBlock(b);
        dim3 dimGrid(ga);

        if (plan->matmode == 0 || plan->matmode == 1)  {
 	   hKernel_none<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                                              plan->xres, plan->yres, plan->tmmode,
                                              plan->dx, plan->dy, plan->dt, plan->dir);
         }
        else 
   	   hKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                                              plan->d_epsilon, plan->d_mu, plan->d_sigma, plan->d_sigast, 
                                              plan->xres, plan->yres, plan->tmmode,
                                              plan->dx, plan->dy, plan->dt, plan->dir);
        
        err = cudaGetLastError();
        if (err) printf("H returned \"%s\"\nH calling sync\n", cudaGetErrorString(err));

        err = cudaThreadSynchronize();
        if (err) printf("H returned \"%s\"\n", cudaGetErrorString(err));
        return cudaGetLastError();
}

cudaError_t wrap_eKernel(SvGpuPlan *plan)
{
        cudaError_t err;

        //cudaEvent_t evt;
        //cudaEventCreate(&evt);
        int ga, b;

	b = plan->xres; //FIXME: do this dependent on problem size and graphics card computation capability
	ga = plan->yres;

	dim3 dimBlock(b);
        dim3 dimGrid(ga);

        /*if (plan->h_isplrc && plan->step==0) {
             printf("Initialising plrc fields");
             plrcnullKernel<<<1, 1>>>(plan->d_plrcx, plan->d_plrcy, plan->d_plrcz,
                                     plan->xres, plan->yres, plan->zres);
        } */
      
        //printf("E Launch %d %d %d %g %g %g %g\n", plan->xres, plan->yres, plan->zres, plan->dx, plan->dy, plan->dz, plan->dt);

	
	if (plan->nmat>0) {
		eKernel_tab<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz, 
                                plan->d_epsilon, plan->d_mu, plan->d_sigma, plan->d_sigast,
				plan->d_mat, plan->d_mattype, plan->d_mattab, plan->nmat, plan->matmode,
				plan->xres, plan->yres, plan->tmmode,
				plan->dx, plan->dy, plan->dt, plan->dir);
	}
	else if (plan->matmode == 0 || plan->matmode == 2)  {
		eKernel_none<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz, 
				plan->xres, plan->yres, plan->tmmode,
				plan->dx, plan->dy, plan->dt, plan->dir);
	}
	else {
		eKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
				plan->d_epsilon, plan->d_mu, plan->d_sigma, plan->d_sigast, 
				plan->xres, plan->yres, plan->tmmode,
				plan->dx, plan->dy, plan->dt, plan->dir);
	}
	err = cudaGetLastError();
//	cudaEventRecord(evt, NULL);
//	while (cudaEventQuery(evt) == cudaErrorNotReady)
//	{
//		usleep(100);
//	}
//	cudaEventDestroy(evt);

        if (err) printf("E returned \"%s\"\nE calling sync\n", cudaGetErrorString(err));
        err = cudaThreadSynchronize();
        if (err) printf("E returned \"%s\"\n", cudaGetErrorString(err));
        return cudaGetLastError();
}
cudaError_t wrap_liao_cpybnd(SvGpuPlan *plan)
{
        int ga, b;

	b = plan->xres; //FIXME: do this dependent on problem size and graphics card computation capability
	ga = plan->yres;

	dim3 dimBlock(b);
        dim3 dimGrid(ga);

        
	liaocpyKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                                              plan->d_liao_x0, plan->d_liao_xn, plan->d_liao_y0, plan->d_liao_yn, 
                                              plan->xres, plan->yres, 
                                              plan->d_bnds,
                                              plan->dx, plan->dy, plan->dt, plan->dir);

        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_liao_applybnd(SvGpuPlan *plan)
{
        int ga, b;

	b = plan->xres; //FIXME: do this dependent on problem size and graphics card computation capability
	ga = plan->yres;

	dim3 dimBlock(b);
        dim3 dimGrid(ga);
        

	liaorunKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz, plan->d_epsilon,
                                              plan->d_liao_x0, plan->d_liao_xn, plan->d_liao_y0, plan->d_liao_yn, 
                                              plan->xres, plan->yres, 
                                              plan->d_bnds, plan->matmode,
                                              plan->dx, plan->dy, plan->dt, plan->dir);

        cudaThreadSynchronize();
        return cudaGetLastError();
}


cudaError_t wrap_mbndx_e_Kernel(SvGpuPlan *plan)
{
        int ga, b;

	b = plan->xres; //FIXME: do this dependent on problem size and graphics card computation capability
	ga = plan->yres;

	dim3 dimBlock(b);
        dim3 dimGrid(ga);

        printf("Error! Not implemented yet.\n");
  /*      mbnxeKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                              plan->xres, plan->yres, plan->zres, plan->mb_bx0, plan->mb_bxn, 
                              plan->mb_bx0pos, plan->mb_bxnpos, 
                              plan->mb_by0pos, plan->mb_bynpos,
                              plan->mb_bz0pos, plan->mb_bznpos, plan->dir);
*/
        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_mbndy_e_Kernel(SvGpuPlan *plan)
{
        int ga, b;

	b = plan->xres; //FIXME: do this dependent on problem size and graphics card computation capability
	ga = plan->yres;

	dim3 dimBlock(b);
        dim3 dimGrid(ga);

        printf("Error! Not implemented yet.\n");
  /*      mbnyeKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                              plan->xres, plan->yres, plan->zres, plan->mb_by0, plan->mb_byn, 
                              plan->mb_bx0pos, plan->mb_bxnpos, 
                              plan->mb_by0pos, plan->mb_bynpos,
                              plan->mb_bz0pos, plan->mb_bznpos, plan->dir);
*/
        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_mbndz_e_Kernel(SvGpuPlan *plan)
{
        int ga, b;

	b = plan->xres; //FIXME: do this dependent on problem size and graphics card computation capability
	ga = plan->yres;

	dim3 dimBlock(b);
        dim3 dimGrid(ga);

       printf("Error! Not implemented yet.\n");
  /*      mbnzeKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                              plan->xres, plan->yres, plan->zres, plan->mb_bz0, plan->mb_bzn, 
                              plan->mb_bx0pos, plan->mb_bxnpos, 
                              plan->mb_by0pos, plan->mb_bynpos,
                              plan->mb_bz0pos, plan->mb_bznpos, plan->dir);
*/
        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_mbndx_h_Kernel(SvGpuPlan *plan)
{
        int ga, b;

	b = plan->xres; //FIXME: do this dependent on problem size and graphics card computation capability
	ga = plan->yres;

	dim3 dimBlock(b);
        dim3 dimGrid(ga);

        printf("Error! Not implemented yet.\n");
  /*      mbnxhKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                              plan->xres, plan->yres, plan->zres, plan->mb_bx0, plan->mb_bxn, 
                              plan->mb_bx0pos, plan->mb_bxnpos, 
                              plan->mb_by0pos, plan->mb_bynpos,
                              plan->mb_bz0pos, plan->mb_bznpos, plan->dir);
*/
        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_mbndy_h_Kernel(SvGpuPlan *plan)
{
        int ga, b;

	b = plan->xres; //FIXME: do this dependent on problem size and graphics card computation capability
	ga = plan->yres;

	dim3 dimBlock(b);
        dim3 dimGrid(ga);

        printf("Error! Not implemented yet.\n");
  /*      mbnyhKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                              plan->xres, plan->yres, plan->zres, plan->mb_by0, plan->mb_byn, 
                              plan->mb_bx0pos, plan->mb_bxnpos, 
                              plan->mb_by0pos, plan->mb_bynpos,
                              plan->mb_bz0pos, plan->mb_bznpos, plan->dir);
*/
        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_srceKernel(SvGpuPlan *plan, int i, int j, int k, float ex, float ey, float ez)
{
        int ga, b;

	b = plan->xres; //FIXME: do this dependent on problem size and graphics card computation capability
	ga = plan->yres;

	dim3 dimBlock(b);
        dim3 dimGrid(ga);

        srcepointKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->xres, plan->yres, i, j, ex, ey, ez, plan->dir);
        cudaThreadSynchronize();
         
        return cudaGetLastError();
}


cudaError_t wrap_fastsumKernel(SvGpuPlan *plan)
{
        int ga, b;

	b = plan->xres; //FIXME: do this dependent on problem size and graphics card computation capability
	ga = plan->yres;

	dim3 dimBlock(b);
        dim3 dimGrid(ga);

        printf("Error! Not implemented yet.\n");
  /*      fastsumKernel_a<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->xres, plan->yres, plan->zres,
                                             plan->d_epsilon, plan->d_mu, plan->d_sum_epsilon, plan->d_sum_sigma, plan->d_sum_mode,
                                             plan->d_sum_i0, plan->d_sum_i1, plan->d_sum_j0, plan->d_sum_j1, plan->d_sum_k0, plan->d_sum_k1,
                                             plan->nsums, plan->d_sum_accumulator, plan->dir,
                                             plan->d_mat, plan->d_mattype, plan->d_mattab, plan->matmode, plan->nmat);
*/
        cudaThreadSynchronize();

        dim3 dimBlockb(plan->nsums);
        //fastsumKernel_b<<<1, dimBlockb>>>(plan->d_sum_accumulator, plan->d_sums, plan->nsums, plan->nsteps, plan->step, plan->xres, plan->yres, plan->zres,
        //                          plan->dir);

        cudaThreadSynchronize();
        return cudaGetLastError();

}



cudaError_t wrap_sumKernel(SvGpuPlan *plan)
{
        dim3 dimBlock(plan->nsums);

        printf("Error! Not implemented yet.\n");
  /*      sumKernel<<<1, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->xres, plan->yres, plan->zres, 
                                   plan->d_epsilon, plan->d_mu, plan->d_sums, plan->d_sum_epsilon, plan->d_sum_sigma, plan->d_sum_mode, 
                                   plan->d_sum_i0, plan->d_sum_i1, plan->d_sum_j0, plan->d_sum_j1, plan->d_sum_k0, plan->d_sum_k1,
                                   plan->nsteps, plan->step,
                                   plan->d_mat, plan->d_mattype, plan->d_mattab, plan->matmode, plan->nmat);
*/
        cudaThreadSynchronize();
         
        return cudaGetLastError();
}


cudaError_t wrap_srchKernel(SvGpuPlan *plan, int i, int j, int k, float hx, float hy, float hz)
{
        int ga, b;

	b = plan->xres; //FIXME: do this dependent on problem size and graphics card computation capability
	ga = plan->yres;

	dim3 dimBlock(b);
        dim3 dimGrid(ga);

        srchpointKernel<<<dimGrid, dimBlock>>>(plan->d_hx, plan->d_hy, plan->d_hz, plan->xres, plan->yres, i, j, hx, hy, hz, plan->dir);
        cudaThreadSynchronize();
         
        return cudaGetLastError();
}

cudaError_t wrap_ffKernel(SvGpuPlan *plan)
{
        dim3 dimBlock(plan->maxthreads);
        dim3 dimGrid((int)ceil((float)plan->h_iset[NPNTS]/plan->maxthreads));
        printf("Error! Not implemented yet.\n");
  /*      ffKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez,
                                        plan->d_iset, 
                                        plan->d_ff_ex, plan->d_ff_ey, plan->d_ff_ez,
                                        plan->xres, plan->yres, plan->zres,
                                        plan->dx, plan->dy, plan->dz, plan->dt, 
                                        plan->step, plan->maxthreads);
        
*/
        cudaThreadSynchronize();        
        return cudaGetLastError();
}


cudaError_t wrap_ffKernel_hlps(SvGpuPlan *plan)
{
        dim3 dimBlock(plan->h_iset[NPNTS]*plan->nhlps); //this should be 512 or 1024
        printf("Error! Not implemented yet.\n");
  /*      ffKernel_a<<<1, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez,
                                        plan->d_iset, 
                                        plan->d_ff_hlp_ex, plan->d_ff_hlp_ey, plan->d_ff_hlp_ez,
                                        plan->xres, plan->yres, plan->zres,
                                        plan->dx, plan->dy, plan->dz, plan->dt, 
                                        plan->step, plan->nhlps);
*/
        cudaThreadSynchronize();    
        return cudaGetLastError();
}

cudaError_t wrap_ffKernel_gethlps(SvGpuPlan *plan)
{
        printf("Error! Not implemented yet.\n");
  /*      ffKernel_b<<<1, 1>>>(plan->d_ff_ex, plan->d_ff_ey, plan->d_ff_ez,
                             plan->d_ff_hlp_ex, plan->d_ff_hlp_ey, plan->d_ff_hlp_ez,
                             plan->d_iset, plan->nhlps);*/
        cudaThreadSynchronize();    
        return cudaGetLastError();
}



cudaError_t wrap_outpointKernel(SvGpuPlan *plan)
{
        outpointKernel<<<1, 1>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                              plan->d_outpointdata, plan->d_outpoint_pos, plan->nsteps, plan->noutpoints,
                              plan->step, plan->xres, plan->yres);              

        cudaThreadSynchronize();
        return cudaGetLastError();


}

cudaError_t wrap_tsf_jstepKernel(SvGpuPlan *plan, float e)
{
        printf("Error! Not implemented yet.\n");
  /*      tsfjstepKernel<<<1, 1>>>(plan->d_tsf_jpool_e, plan->d_tsf_jpool_h, plan->d_tsf_jpvals, 
                              plan->d_tsf_jpool_epsilon, plan->d_tsf_jpool_mu,
                              plan->d_tsf_jpool_sigma, plan->d_tsf_jpool_sigast,
                              plan->dx/plan->h_tsfset[9], plan->dt, e, plan->tsf_jpool_size);
*/
        cudaThreadSynchronize();
        return cudaGetLastError();
}


cudaError_t wrap_tsf_e_Kernel(SvGpuPlan *plan)
{
        int ga, b;

	b = plan->xres; //FIXME: do this dependent on problem size and graphics card computation capability
	ga = plan->yres;

	dim3 dimBlock(b);
        dim3 dimGrid(ga);

        printf("Error! Not implemented yet.\n");
/*
        tsf_estep_aKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, 
                                         plan->d_hx, plan->d_hy, plan->d_hz,
                                         plan->d_epsilon, plan->d_mu, plan->matmode,
                                         plan->d_tsfset,
                                         plan->d_tsf_jpool_e, plan->d_tsf_jpool_h, plan->tsf_jpool_size,
                                         plan->dx, plan->dt,
                                         plan->xres, plan->yres, plan->zres, plan->dir);
        tsf_estep_bKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, 
                                         plan->d_hx, plan->d_hy, plan->d_hz,
                                         plan->d_epsilon, plan->d_mu, plan->matmode,
                                         plan->d_tsfset,
                                         plan->d_tsf_jpool_e, plan->d_tsf_jpool_h, plan->tsf_jpool_size,
                                         plan->dx, plan->dt,
                                         plan->xres, plan->yres, plan->zres, plan->dir);
*/

        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_tsf_h_Kernel(SvGpuPlan *plan)
{
        int ga, b;

	b = plan->xres; //FIXME: do this dependent on problem size and graphics card computation capability
	ga = plan->yres;

	dim3 dimBlock(b);
        dim3 dimGrid(ga);

        printf("Error! Not implemented yet.\n");
/*
        tsf_hstep_aKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, 
                                         plan->d_hx, plan->d_hy, plan->d_hz,
                                         plan->d_epsilon, plan->d_mu, plan->matmode,
                                         plan->d_tsfset,
                                         plan->d_tsf_jpool_e, plan->d_tsf_jpool_h, plan->tsf_jpool_size,
                                         plan->dx, plan->dt,
                                         plan->xres, plan->yres, plan->zres, plan->dir);
        tsf_hstep_bKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, 
                                         plan->d_hx, plan->d_hy, plan->d_hz,
                                         plan->d_epsilon, plan->d_mu, plan->matmode,
                                         plan->d_tsfset,
                                         plan->d_tsf_jpool_e, plan->d_tsf_jpool_h, plan->tsf_jpool_size,
                                         plan->dx, plan->dt,
                                         plan->xres, plan->yres, plan->zres, plan->dir);
*/

        cudaThreadSynchronize();
        return cudaGetLastError();
}


/*
__global__ void kernel() { ... }
cudaError_t kernel_driver()
        {
        kernel<<<blocks, threads>>>();
        #ifdef NDEBUG
        return cudaSuccess;
        #else
        cudaThreadSynchronize();
        return cudaGetLastError();
        #endif
        }
*/


#endif

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
