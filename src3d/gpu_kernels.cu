
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

#define CAT(i, j, k) ((i)*zres*yres + (j)*zres + (k))

#define CLX(off, j, k)  (off*yres*zres + (j)*zres + (k)) 
#define CLXP(off, j, k) ((off + 1)*yres*zres + (j)*zres + (k))
#define CLY(off, i, k)  (off*xres*zres + (k)*xres + (i)) 
#define CLYP(off, i, k) ((off + 1)*xres*zres + (k)*xres + (i))
#define CLZ(off, i, j)  (off*xres*yres + (j)*xres + (i)) 
#define CLZP(off, i, j) ((off + 1)*xres*yres + (j)*xres + (i))

#define CX0AT(i,j,k) ((i)*zres*yres + (j)*zres + (k))
#define CX0_EYX (0)
#define CX0_EZX (depths[BNDX0]*yres*zres)
#define CX0_HYX (2*depths[BNDX0]*yres*zres)
#define CX0_HZX (3*depths[BNDX0]*yres*zres)

#define CXNAT(i,j,k) ((i)*zres*yres + (j)*zres + (k))
#define CXN_EYX (0)
#define CXN_EZX (depths[BNDXN]*yres*zres)
#define CXN_HYX (2*depths[BNDXN]*yres*zres)
#define CXN_HZX (3*depths[BNDXN]*yres*zres)

#define CY0AT(i,j,k) ((i)*zres*depths[BNDY0] + (j)*zres + (k))
#define CY0_EXY (0)
#define CY0_EZY (xres*depths[BNDY0]*zres)
#define CY0_HXY (2*xres*depths[BNDY0]*zres)
#define CY0_HZY (3*xres*depths[BNDY0]*zres)

#define CYNAT(i,j,k) ((i)*zres*depths[BNDYN] + (j)*zres + (k))
#define CYN_EXY (0)
#define CYN_EZY (xres*depths[BNDYN]*zres)
#define CYN_HXY (2*xres*depths[BNDYN]*zres)
#define CYN_HZY (3*xres*depths[BNDYN]*zres)

#define CZ0AT(i,j,k) ((i)*depths[BNDZ0]*yres + (j)*depths[BNDZ0] + (k))
#define CZ0_EXZ (0)
#define CZ0_EYZ (xres*yres*depths[BNDZ0])
#define CZ0_HXZ (2*xres*yres*depths[BNDZ0])
#define CZ0_HYZ (3*xres*yres*depths[BNDZ0])

#define CZNAT(i,j,k) ((i)*depths[BNDZN]*yres + (j)*depths[BNDZN] + (k))
#define CZN_EXZ (0)
#define CZN_EYZ (xres*yres*depths[BNDZN])
#define CZN_HXZ (2*xres*yres*depths[BNDZN])
#define CZN_HYZ (3*xres*yres*depths[BNDZN])




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

#define PIMIN 6
#define PJMIN 7
#define PIMAX 8
#define PJMAX 9

#define PSKIPK0 11
#define PSKIPK0_IMIN 12
#define PSKIPK0_JMIN 13
#define PSKIPK0_IMAX 14
#define PSKIPK0_JMAX 15
#define PSKIPKN 16
#define PSKIPKN_IMIN 17
#define PSKIPKN_JMIN 18
#define PSKIPKN_IMAX 19
#define PSKIPKN_JMAX 20

#define BNDX0 0
#define BNDXN 1
#define BNDY0 2
#define BNDYN 3
#define BNDZ0 4
#define BNDZN 5

#define NPNTS 36
#define PNPNTS 20
#define PNDIV 200
//nasleduji jednotlive body x,y,z

#define GETEPS(i) (SV_GM_N*i + SV_GM_EPSILON)
#define GETMU(i) (SV_GM_N*i + SV_GM_MU)
#define GETSIGMA(i) (SV_GM_N*i + SV_GM_SIGMA)
#define GETSIGAST(i) (SV_GM_N*i + SV_GM_SIGAST)
#define GETDRUDEO(i) (SV_GM_N*i + SV_GM_DRUDE_OMEGA_P)
#define GETDRUDEN(i) (SV_GM_N*i + SV_GM_DRUDE_NU)
#define GETCPA0(i) (SV_GM_N*i + SV_GM_CP3_A0)
#define GETCPA1(i) (SV_GM_N*i + SV_GM_CP3_A1)
#define GETCPA2(i) (SV_GM_N*i + SV_GM_CP3_A2)
#define GETCPPHI0(i) (SV_GM_N*i + SV_GM_CP3_PHI0)
#define GETCPPHI1(i) (SV_GM_N*i + SV_GM_CP3_PHI1)
#define GETCPPHI2(i) (SV_GM_N*i + SV_GM_CP3_PHI2)
#define GETCPOMEGA0(i) (SV_GM_N*i + SV_GM_CP3_OMEGA0)
#define GETCPOMEGA1(i) (SV_GM_N*i + SV_GM_CP3_OMEGA1)
#define GETCPOMEGA2(i) (SV_GM_N*i + SV_GM_CP3_OMEGA2)
#define GETCPGAMMA0(i) (SV_GM_N*i + SV_GM_CP3_GAMMA0)
#define GETCPGAMMA1(i) (SV_GM_N*i + SV_GM_CP3_GAMMA1)
#define GETCPGAMMA2(i) (SV_GM_N*i + SV_GM_CP3_GAMMA2)

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



__global__  void
peKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_epsilon, float *d_mu, float *d_sigma, float *d_sigast, int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int dir, int nthreads, int tw)
{
    int i, j, k;
    float ca, cb, hxvat, hyvat, hzvat, sigma, epsilon;
    long int pos;
    int zyres = zres*yres;

    //pos = blockIdx.x*blockIdx.y*nthreads + threadIdx.x; 
    pos = (blockIdx.x*tw + blockIdx.y)*nthreads + threadIdx.x;

    i = (int)(pos/(zyres));
    j = (pos - i*zyres)/zres;
    k = pos - i*zyres - j*zres;

    if (i<0 || j<0 || k<0) return;
    if (i>=(xres) || j>=(yres) || k>=zres) return;

    hxvat = d_hx[pos];
    hyvat = d_hy[pos];
    hzvat = d_hz[pos];
    sigma = d_sigma[pos];
    epsilon = d_epsilon[pos];

    ca  = (1 - sigma*dt/2/epsilon)/(1 + sigma*dt/2/epsilon);
    cb  = (dt/epsilon/dx)/(1 + sigma*dt/2/epsilon);

    if (j>0 && k>0)
       d_ex[pos] = ca*d_ex[pos] + cb*((hzvat - d_hz[pos-zres]) -
               (hyvat - d_hy[pos-1]));                       

    if (i>0 && k>0)
       d_ey[pos] = ca*d_ey[pos] + cb*((hxvat - d_hx[pos-1]) -  
               (hzvat - d_hz[pos-zyres]));

    if (i>0 && j>0)
       d_ez[pos] = ca*d_ez[pos] + cb*((hyvat - d_hy[pos-zyres]) -
               (hxvat - d_hx[pos-zres]));

}

__global__  void
eKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_epsilon, float *d_mu, float *d_sigma, float *d_sigast, int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k;
    float ca, cb, hxvat, hyvat, hzvat;
    long int pos;

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
 
    if (i<0 || j<0 || k<0) return;
    if (i>=(xres) || j>=(yres) || k>=zres) return;
 
    pos = CAT(i, j, k);
    hxvat = d_hx[pos];
    hyvat = d_hy[pos];
    hzvat = d_hz[pos];
    ca  = (1 - d_sigma[pos]*dt/2/d_epsilon[pos])/(1 + d_sigma[pos]*dt/2/d_epsilon[pos]);
    cb  = (dt/d_epsilon[pos]/dx)/(1 + d_sigma[pos]*dt/2/d_epsilon[pos]);

    if (j>0 && k>0) 
       d_ex[pos] = ca*d_ex[pos] + cb*((hzvat - d_hz[CAT(i, j-1, k)]) -
               (hyvat - d_hy[CAT(i, j, k-1)]));

    if (i>0 && k>0)  
       d_ey[pos] = ca*d_ey[pos] + cb*((hxvat - d_hx[CAT(i, j, k-1)]) -
               (hzvat - d_hz[CAT(i-1, j, k)]));

    if (i>0 && j>0) 
       d_ez[pos] = ca*d_ez[pos] + cb*((hyvat - d_hy[CAT(i-1, j, k)]) -
               (hxvat - d_hx[CAT(i, j-1, k)]));

}


__device__ void cmult(float a, float b, float c, float d, float *ra, float *rb)
{
    *ra = a*c - b*d;
    *rb = a*d + b*c;
}

__device__ void cdiv(float a, float b, float c, float d, float *ra, float *rb)
{
    float sub = c*c + d*d;
    *ra = (a*c + b*d)/sub;
    *rb = (b*c - a*d)/sub;
}

__device__ void
fract(float a, float gamma, float omega, float phi, float dt, float *re, float *im)
{
    float r_sub, i_sub;
    float r_sup, i_sup;

    r_sub = gamma;
    i_sub = -omega;
    r_sup = -2*a*omega*sin(phi);//was plus
    i_sup = -2*a*omega*cos(phi);

    cdiv(r_sup, i_sup, r_sub, i_sub, re, im);
}

__device__ void
dxxin0(float dt, float *re, float *im, int *d_mat, float *d_mattab, long int pos)
{
    float r_efact, i_efact;
    float r_sup, i_sup;

    r_efact = 1 - exp(-d_mattab[GETCPGAMMA0(d_mat[pos])]*dt)*cos(d_mattab[GETCPOMEGA0(d_mat[pos])]*dt);
    i_efact = -exp(-d_mattab[GETCPGAMMA0(d_mat[pos])]*dt)*sin(d_mattab[GETCPOMEGA0(d_mat[pos])]*dt);

    fract(d_mattab[GETCPA0(d_mat[pos])], d_mattab[GETCPGAMMA0(d_mat[pos])], d_mattab[GETCPOMEGA0(d_mat[pos])], d_mattab[GETCPPHI0(d_mat[pos])], dt, &r_sup, &i_sup);

    cmult(r_efact, i_efact, r_efact, i_efact, &r_efact, &i_efact);
    cmult(r_sup, i_sup, r_efact, i_efact, &r_sup, &i_sup);

    *re = r_sup;
    *im = i_sup;
}
__device__ void
dxxin1(float dt, float *re, float *im, int *d_mat, float *d_mattab, long int pos)
{
    float r_efact, i_efact;
    float r_sup, i_sup;

    r_efact = 1 - exp(-d_mattab[GETCPGAMMA1(d_mat[pos])]*dt)*cos(d_mattab[GETCPOMEGA1(d_mat[pos])]*dt);
    i_efact = -exp(-d_mattab[GETCPGAMMA1(d_mat[pos])]*dt)*sin(d_mattab[GETCPOMEGA1(d_mat[pos])]*dt);

    fract(d_mattab[GETCPA1(d_mat[pos])], d_mattab[GETCPGAMMA1(d_mat[pos])], d_mattab[GETCPOMEGA1(d_mat[pos])], d_mattab[GETCPPHI1(d_mat[pos])], dt, &r_sup, &i_sup);

    cmult(r_efact, i_efact, r_efact, i_efact, &r_efact, &i_efact);
    cmult(r_sup, i_sup, r_efact, i_efact, &r_sup, &i_sup);

    *re = r_sup;
    *im = i_sup;
}

__device__ float
xxi0(float dt, int *d_mat, float *d_mattab, long int pos)
{
    float r_efact, i_efact;
    float omega, gamma;
    float r_sup, i_sup;
    float sum_sup = 0;

    r_efact = 1 - exp(-d_mattab[GETCPGAMMA0(d_mat[pos])]*dt)*cos(d_mattab[GETCPOMEGA0(d_mat[pos])]*dt);
    i_efact = -exp(-d_mattab[GETCPGAMMA0(d_mat[pos])]*dt)*sin(d_mattab[GETCPOMEGA0(d_mat[pos])]*dt);

    fract(d_mattab[GETCPA0(d_mat[pos])], d_mattab[GETCPGAMMA0(d_mat[pos])], d_mattab[GETCPOMEGA0(d_mat[pos])], d_mattab[GETCPPHI0(d_mat[pos])], dt, &r_sup, &i_sup);
    cmult(r_sup, i_sup, r_efact, i_efact, &r_sup, &i_sup);

    sum_sup += r_sup;

    r_efact = 1 - exp(-d_mattab[GETCPGAMMA1(d_mat[pos])]*dt)*cos(d_mattab[GETCPOMEGA1(d_mat[pos])]*dt);
    i_efact = -exp(-d_mattab[GETCPGAMMA1(d_mat[pos])]*dt)*sin(d_mattab[GETCPOMEGA1(d_mat[pos])]*dt);

    fract(d_mattab[GETCPA1(d_mat[pos])], d_mattab[GETCPGAMMA1(d_mat[pos])], d_mattab[GETCPOMEGA1(d_mat[pos])], d_mattab[GETCPPHI1(d_mat[pos])], dt, &r_sup, &i_sup);
    cmult(r_sup, i_sup, r_efact, i_efact, &r_sup, &i_sup);

    sum_sup += r_sup;

    omega = d_mattab[GETDRUDEO(d_mat[pos])];
    gamma = d_mattab[GETDRUDEN(d_mat[pos])];

    return -(omega*omega/gamma/gamma)*(1-exp(-gamma*dt)) + dt*omega*omega/gamma + sum_sup;
}


__device__ void
update_rc_cpn0(float *plrc, float *iplrc, float e,
               float dt, int *d_mat, float *d_mattab, long int pos)
{
    float r_efact, i_efact;
    float r_exp, i_exp;
    float r_res, i_res;

    r_exp = exp(-d_mattab[GETCPGAMMA0(d_mat[pos])]*dt)*cos(d_mattab[GETCPOMEGA0(d_mat[pos])]*dt);
    i_exp = exp(-d_mattab[GETCPGAMMA0(d_mat[pos])]*dt)*sin(d_mattab[GETCPOMEGA0(d_mat[pos])]*dt);
    cmult(r_exp, i_exp, *plrc, *iplrc, &r_res, &i_res);

    dxxin0(dt, &r_efact, &i_efact, d_mat, d_mattab, pos);

    *plrc = r_res + r_efact*e;
    *iplrc = i_res + i_efact*e;
}
__device__ void
update_rc_cpn1(float *plrc, float *iplrc, float e,
               float dt, int *d_mat, float *d_mattab, long int pos)
{
    float r_efact, i_efact;
    float r_exp, i_exp;
    float r_res, i_res;

    r_exp = exp(-d_mattab[GETCPGAMMA1(d_mat[pos])]*dt)*cos(d_mattab[GETCPOMEGA1(d_mat[pos])]*dt);
    i_exp = exp(-d_mattab[GETCPGAMMA1(d_mat[pos])]*dt)*sin(d_mattab[GETCPOMEGA1(d_mat[pos])]*dt);
    cmult(r_exp, i_exp, *plrc, *iplrc, &r_res, &i_res);

    dxxin1(dt, &r_efact, &i_efact, d_mat, d_mattab, pos);

    *plrc = r_res + r_efact*e;
    *iplrc = i_res + i_efact*e;
}


__device__ void	
cp_update(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz, int *d_mat, int *d_mattype, float *d_mattab, int nmat, int matmode,
                  float *d_plrcx, float *d_plrcy, float *d_plrcz, int i, int j, int k, int xres, int yres, int zres, float dt, float dx)
{
    long int pos = CAT(i, j, k);
    long int shift = xres*yres*zres;
    float ca, cb, cc, xxival;

    float pex = d_ex[pos];
    float pey = d_ey[pos];
    float pez = d_ez[pos];

    float omega = d_mattab[GETDRUDEO(d_mat[pos])];
    float gamma = d_mattab[GETDRUDEN(d_mat[pos])];
    float epsilon = d_mattab[GETEPS(d_mat[pos])]/EPSILON_0;

    d_plrcx[pos] = exp(-gamma*dt)*d_plrcx[pos]
	    - (omega*omega/gamma/gamma)*(1-exp(-gamma*dt))*(1-exp(-gamma*dt))*pex;

    d_plrcy[pos] = exp(-gamma*dt)*d_plrcy[pos]
	    - (omega*omega/gamma/gamma)*(1-exp(-gamma*dt))*(1-exp(-gamma*dt))*pey;

    d_plrcz[pos] = exp(-gamma*dt)*d_plrcz[pos]
	    - (omega*omega/gamma/gamma)*(1-exp(-gamma*dt))*(1-exp(-gamma*dt))*pez;


    update_rc_cpn0(d_plrcx + shift, d_plrcx + 2*shift, pex, dt, d_mat, d_mattab, pos);
    update_rc_cpn0(d_plrcy + shift, d_plrcy + 2*shift, pey, dt, d_mat, d_mattab, pos);
    update_rc_cpn0(d_plrcz + shift, d_plrcz + 2*shift, pez, dt, d_mat, d_mattab, pos);

    update_rc_cpn0(d_plrcx + 3*shift, d_plrcx + 4*shift, pex, dt, d_mat, d_mattab, pos);
    update_rc_cpn0(d_plrcy + 3*shift, d_plrcy + 4*shift, pey, dt, d_mat, d_mattab, pos);
    update_rc_cpn0(d_plrcz + 3*shift, d_plrcz + 4*shift, pez, dt, d_mat, d_mattab, pos);


    xxival = xxi0(dt, d_mat, d_mattab, pos);
    ca = epsilon/(epsilon + xxival);
    cb = dt/(dx*EPSILON_0*(epsilon + xxival));
    cc = 1.0/(epsilon + xxival);


    d_ex[pos] = ca*pex
	    + cc*(d_plrcx[pos] + d_plrcx[pos + shift] + d_plrcx[pos + 3*shift])
	    + cb*((d_hz[pos] - d_hz[CAT(i, j-1, k)]) - (d_hy[pos] - d_hy[CAT(i,j,k-1)]));

    d_ey[pos] = ca*pey
	    + cc*(d_plrcy[pos] + d_plrcy[pos + shift] + d_plrcy[pos + 3*shift])
	    + cb*((d_hx[pos] - d_hx[CAT(i,j,k-1)]) - (d_hz[pos] - d_hz[CAT(i-1,j,k)]));

    d_ez[pos] = ca*pez
	    + cc*(d_plrcz[pos] + d_plrcz[pos + shift] + d_plrcz[pos + 3*shift])
	    + cb*((d_hy[pos] - d_hy[CAT(i-1,j,k)]) - (d_hx[pos] - d_hx[CAT(i,j-1,k)]));

}

__device__ void	
plrc_update(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz, int *d_mat, int *d_mattype, float *d_mattab, int nmat, int matmode,
                  float *d_plrcx, float *d_plrcy, float *d_plrcz, int i, int j, int k, int xres, int yres, int zres, float dt, float dx)
{
    long int pos = CAT(i, j, k);
    long int shift = xres*yres*zres;
    float ca, cb, cc, a, b, c, d, dchirdiff, dchiidiff;
    int n;

    float pex = d_ex[pos];
    float pey = d_ey[pos];
    float pez = d_ez[pos];

    float sumxi0 = d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_D_XI)]  + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_P_XI)] + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_P_XI+1)];
    float sumchi0 = d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_D_CHI)]  + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_P_CHI)] + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_P_CHI+1)];
    float sumpsix = d_plrcx[pos] + d_plrcx[pos + shift] + d_plrcx[pos + 3*shift];
    float sumpsiy = d_plrcy[pos] + d_plrcy[pos + shift] + d_plrcy[pos + 3*shift];
    float sumpsiz = d_plrcz[pos] + d_plrcz[pos + shift] + d_plrcz[pos + 3*shift];
    float epsilon = d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_EPSILON)]/EPSILON_0;

    ca = (2*EPSILON_0*(epsilon - sumxi0) - d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_SIGMA)]*dt)/(2*EPSILON_0*(epsilon - sumxi0 + sumchi0) + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_SIGMA)]*dt);
    cb = 2*EPSILON_0/(2*EPSILON_0*(epsilon - sumxi0 + sumchi0) + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_SIGMA)]*dt);
    cc = 2*dt/(2*EPSILON_0*(epsilon - sumxi0 + sumchi0) + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_SIGMA)]*dt)/dx;


    d_ex[pos] = ca*pex
	    + cb*sumpsix
	    + cc*((d_hz[pos] - d_hz[CAT(i, j-1, k)]) - (d_hy[pos] - d_hy[CAT(i,j,k-1)]));

    d_ey[pos] = ca*pey
	    + cb*sumpsiy
	    + cc*((d_hx[pos] - d_hx[CAT(i,j,k-1)]) - (d_hz[pos] - d_hz[CAT(i-1,j,k)]));

    d_ez[pos] = ca*pez
	    + cb*sumpsiz
	    + cc*((d_hy[pos] - d_hy[CAT(i-1,j,k)]) - (d_hx[pos] - d_hx[CAT(i,j-1,k)]));


    dchirdiff = (d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_D_DCHI)] - d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_D_DXI)]);

    d_plrcx[pos] = dchirdiff*d_ex[pos]
	    + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_D_DXI)]*pex + exp(-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_DRUDE_NU)]*dt)*d_plrcx[pos];

    d_plrcy[pos] = dchirdiff*d_ey[pos]
	    + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_D_DXI)]*pey + exp(-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_DRUDE_NU)]*dt)*d_plrcy[pos];

    d_plrcz[pos] = dchirdiff*d_ez[pos]
	    + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_D_DXI)]*pez + exp(-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_DRUDE_NU)]*dt)*d_plrcz[pos];


    for (n=0; n<2; n++) {
                                   

         a = exp(-dt*d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_CP3_GAMMA0 + n)])*cos(dt*d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_CP3_OMEGA0 + n)]);
         b = -exp(-dt*d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_CP3_GAMMA0 + n)])*sin(dt*d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_CP3_OMEGA0 + n)]);
         dchirdiff = (d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_P_DCHIR + n)] - d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_P_DXIR + n)]);
         dchiidiff = (d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_P_DCHII + n)] - d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_P_DXII + n)]);

         c = d_plrcx[pos + (2*n + 1)*shift];
         d = d_plrcx[pos + (2*n + 2)*shift];

         d_plrcx[pos + (2*n + 1)*shift] = dchirdiff*d_ex[pos]
                                        + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_P_DXIR + n)]*pex + a*c - b*d;
         d_plrcx[pos + (2*n + 2)*shift] = dchiidiff*d_ex[pos]
                                        + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_P_DXII + n)]*pex + a*d + b*c;


         c = d_plrcy[pos + (2*n + 1)*shift];
         d = d_plrcy[pos + (2*n + 2)*shift];


         d_plrcy[pos + (2*n + 1)*shift] = dchirdiff*d_ey[pos]
                                        + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_P_DXIR + n)]*pey + a*c - b*d;
         d_plrcy[pos + (2*n + 2)*shift] = dchiidiff*d_ey[pos]
                                        + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_P_DXII + n)]*pey + a*d + b*c;


         c = d_plrcz[pos + (2*n + 1)*shift];
         d = d_plrcz[pos + (2*n + 2)*shift];


         d_plrcz[pos + (2*n + 1)*shift] = dchirdiff*d_ez[pos]
                                        + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_P_DXIR + n)]*pez + a*c - b*d;
         d_plrcz[pos + (2*n + 2)*shift] = dchiidiff*d_ez[pos]
                                        + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_PLRC_P_DXII + n)]*pez + a*d + b*c;
    }


}
 
__device__ void	
ade_update(float *d_ex, float *d_ey, float *d_ez, float *d_exp, float *d_eyp, float *d_ezp, float *d_hx, float *d_hy, float *d_hz, int *d_mat, int *d_mattype, float *d_mattab, int nmat, int matmode,
                  float *d_px, float *d_pxp, float *d_dpx, float *d_dpxp, int i, int j, int k, int xres, int yres, int zres, float dt, float dx)
{
    int n;
    long int pos = CAT(i, j, k);
    long int shift = xres*yres*zres;
    float vdpx, vpx[2], vdpy, vpy[2], vdpz, vpz[2];

    float pex = d_ex[pos];
    float pey = d_ey[pos];
    float pez = d_ez[pos];

    float c0 = d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_C0)]/dx;
    float c1 = d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_C1)];
    float c2 = d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_C2)];
    float c3 = d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_C3)];

    d_ex[pos] = c0*((d_hz[pos] - d_hz[CAT(i, j-1, k)]) - (d_hy[pos] - d_hy[CAT(i,j,k-1)]))
               + c1*(-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_A0)]*d_dpxp[pos] + (1-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_A1)])*d_dpx[pos])
               + c1*(-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP0)]*d_pxp[pos] + (1-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP1)])*d_px[pos])
               + c1*(-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP0 + 1)]*d_pxp[pos + 3*shift] + (1-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP1 + 1)])*d_px[pos + 3*shift])
               + c2*d_exp[pos] + c3*pex;
    
    d_ey[pos] = c0*((d_hx[pos] - d_hx[CAT(i,j,k-1)]) - (d_hz[pos] - d_hz[CAT(i-1,j,k)]))
               + c1*(-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_A0)]*d_dpxp[pos + shift] + (1-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_A1)])*d_dpx[pos + shift])
               + c1*(-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP0)]*d_pxp[pos + shift] + (1-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP1)])*d_px[pos + shift])
               + c1*(-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP0 + 1)]*d_pxp[pos + 4*shift] + (1-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP1 + 1)])*d_px[pos + 4*shift])
               + c2*d_eyp[pos] + c3*pey;

    d_ez[pos] = c0*((d_hy[pos] - d_hy[CAT(i-1,j,k)]) - (d_hx[pos] - d_hx[CAT(i,j-1,k)]))
               + c1*(-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_A0)]*d_dpxp[pos + 2*shift] + (1-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_A1)])*d_dpx[pos + 2*shift])
               + c1*(-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP0)]*d_pxp[pos + 2*shift] + (1-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP1)])*d_px[pos + 2*shift])
               + c1*(-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP0 + 1)]*d_pxp[pos + 5*shift] + (1-d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP1 + 1)])*d_px[pos + 5*shift])
               + c2*d_ezp[pos] + c3*pez;
  

    /*store previous p terms*/
    vdpx = d_dpx[pos];
    vdpy = d_dpx[pos + shift];
    vdpz = d_dpx[pos + 2*shift];
    for (n=0; n<2; n++) {
	    vpx[n] = d_px[pos + 3*n*shift];
	    vpy[n] = d_px[pos + (3*n+1)*shift];
	    vpz[n] = d_px[pos + (3*n+2)*shift];
    }

    /*update all p terms*/
    d_dpx[pos] = d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_A0)]*d_dpxp[pos] + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_A1)]*d_dpx[pos]
	    + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_A2)]*(d_exp[pos] + 2*pex + d_ex[pos]);
    d_dpx[pos + shift] = d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_A0)]*d_dpxp[pos + shift] + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_A1)]*d_dpx[pos + shift]
	    + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_A2)]*(d_eyp[pos] + 2*pey + d_ey[pos]);
    d_dpx[pos + 2*shift] = d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_A0)]*d_dpxp[pos + 2*shift] + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_A1)]*d_dpx[pos + 2*shift]
	    + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_A2)]*(d_ezp[pos] + 2*pez + d_ez[pos]);


    for (n=0; n<2; n++) {

            d_px[pos + 3*n*shift]     = d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP0 + n)]*d_pxp[pos + 3*n*shift] + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP1 + n)]*d_px[pos + 3*n*shift]
	               + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP2 + n)]*d_exp[pos] + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP3 + n)]*pex + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP4 + n)]*d_ex[pos];

            d_px[pos + (3*n+1)*shift] = d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP0 + n)]*d_pxp[pos + (3*n+1)*shift] + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP1 + n)]*d_px[pos + (3*n+1)*shift]
	               + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP2 + n)]*d_eyp[pos] + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP3 + n)]*pey + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP4 + n)]*d_ey[pos];

            d_px[pos + (3*n+2)*shift] = d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP0 + n)]*d_pxp[pos + (3*n+2)*shift] + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP1 + n)]*d_px[pos + (3*n+2)*shift]
	               + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP2 + n)]*d_ezp[pos] + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP3 + n)]*pez + d_mattab[(SV_GM_N*d_mat[pos] + SV_GM_ADE_BP4 + n)]*d_ez[pos];

    }

    /*save e n-1*/
    d_exp[pos] = pex;
    d_eyp[pos] = pey;
    d_ezp[pos] = pez;

    /*save previous p terms*/
    d_dpxp[pos] = vdpx;
    d_dpxp[pos + shift] = vdpy;
    d_dpxp[pos + 2*shift] = vdpz;
    for (n=0; n<2; n++) {
	    d_pxp[pos + 3*n*shift] = vpx[n];
	    d_pxp[pos + (3*n+1)*shift] = vpy[n];
	    d_pxp[pos + (3*n+2)*shift] = vpz[n];
    }

}

__device__ float
dch0(float omegadiv, float dt)
{
    return -omegadiv*omegadiv*(1.0-exp(-LIGHT_SPEED*dt))*(1.0-exp(-LIGHT_SPEED*dt));
}

__device__ float
ch0(float omegadiv, float dt)
{
    return omegadiv*omegadiv*LIGHT_SPEED*dt - omegadiv*omegadiv*(1.0-exp(-LIGHT_SPEED*dt));
}

__device__ float
dch1(float omegadiv, float dt)
{
    return -omegadiv*omegadiv*exp(-LIGHT_SPEED*dt)*(1.0-exp(-LIGHT_SPEED*dt))*(1.0-exp(-LIGHT_SPEED*dt));
}


__device__ void	
drude_update(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz, int *d_mat, int *d_mattype, float *d_mattab, int nmat, int matmode,
                  float *d_plrcx, float *d_plrcy, float *d_plrcz, int i, int j, int k, int xres, int yres, int zres, float dt, float dx)
{
    long int pos = CAT(i, j, k);
    float ca, cb, cc;

    float pex = d_ex[pos];
    float pey = d_ey[pos];
    float pez = d_ez[pos];

    float vch0 = ch0(d_mattab[GETDRUDEO(d_mat[pos])], dt);
    float vdch1 = dch1(d_mattab[GETDRUDEO(d_mat[pos])], dt);
    float vdch0 = dch0(d_mattab[GETDRUDEO(d_mat[pos])], dt);
    float epsilon = d_mattab[GETEPS(d_mat[pos])]/EPSILON_0;

    ca = (epsilon+vdch0)/(epsilon+vch0);
    cb = dt/(epsilon+vch0)/EPSILON_0/dx;
    cc = 1.0/(epsilon+vch0);

    d_ex[pos] = ca*pex
	    + cc*(d_plrcx[pos])
	    + cb*((d_hz[pos] - d_hz[CAT(i, j-1, k)]) - (d_hy[pos] - d_hy[CAT(i,j,k-1)]));

    d_ey[pos] = ca*pey
	    + cc*(d_plrcy[pos])
	    + cb*((d_hx[pos] - d_hx[CAT(i,j,k-1)]) - (d_hz[pos] - d_hz[CAT(i-1,j,k)]));

    d_ez[pos] = ca*pez
	    + cc*(d_plrcz[pos])
	    + cb*((d_hy[pos] - d_hy[CAT(i-1,j,k)]) - (d_hx[pos] - d_hx[CAT(i,j-1,k)]));

    d_plrcx[pos] = exp(-d_mattab[GETDRUDEN(d_mat[pos])]*dt)*d_plrcx[pos]
	    + vdch1*pex;

    d_plrcy[pos] = exp(-d_mattab[GETDRUDEN(d_mat[pos])]*dt)*d_plrcy[pos]
	    + vdch1*pey;

    d_plrcz[pos] = exp(-d_mattab[GETDRUDEN(d_mat[pos])]*dt)*d_plrcz[pos]
	    + vdch1*pez;

//    d_ex[pos] = d_ey[pos] = d_ez[pos] = 0;
}

__global__  void
eKernel_tab(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_epsilon, float *d_mu, float *d_sigma, float *d_sigast,
         int *d_mat, int *d_mattype, float *d_mattab, int nmat, int matmode,
         int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int dir, float *d_sf_jpool_e, float *d_sfset, int d_sf_jpool_size,
         float *d_plrcx, float *d_plrcy, float *d_plrcz, float *d_exp, float *d_eyp, float *d_ezp, 
         float *d_px, float *d_pxp, float *d_dpx, float *d_dpxp)
{
    int i, j, k, mattype;
    int mattype_xm, mattype_ym, mattype_zm, mattype_xmym, mattype_xmzm, mattype_ymzm;
    float ca, cb, hxvat, hyvat, hzvat, sigma, epsilon, theta, phi, psi, d, ve;
    long int pos;

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
 
    if (i<0 || j<0 || k<0) return;
    if (i>=(xres) || j>=(yres) || k>=zres) return;
 
    pos = CAT(i, j, k);
    if (i>0 && j>0 && k>0) {
	    mattype_xm = d_mattype[d_mat[CAT(i-1, j, k)]];
	    mattype_ym = d_mattype[d_mat[CAT(i, j-1, k)]];
	    mattype_zm = d_mattype[d_mat[CAT(i, j, k-1)]];
	    mattype_xmym = d_mattype[d_mat[CAT(i-1, j-1, k)]];
	    mattype_xmzm = d_mattype[d_mat[CAT(i-1, j, k-1)]];
	    mattype_ymzm = d_mattype[d_mat[CAT(i, j-1, k-1)]];
    }

    mattype = d_mattype[d_mat[pos]];

    hxvat = d_hx[pos];
    hyvat = d_hy[pos];
    hzvat = d_hz[pos];
  
    if (mattype==2) {
       	drude_update(d_ex, d_ey, d_ez, d_hx, d_hy, d_hz, d_mat, d_mattype, d_mattab, nmat, matmode,
                  d_plrcx, d_plrcy, d_plrcz, i, j, k, xres, yres, zres, dt, dx); 
    }
    else if (mattype==4) { //cp model
	cp_update(d_ex, d_ey, d_ez, d_hx, d_hy, d_hz, d_mat, d_mattype, d_mattab, nmat, matmode,
                  d_plrcx, d_plrcy, d_plrcz, i, j, k, xres, yres, zres, dt, dx);
    } 
    else if (mattype==5) { //ade model
	ade_update(d_ex, d_ey, d_ez, d_exp, d_eyp, d_ezp, d_hx, d_hy, d_hz, d_mat, d_mattype, d_mattab, nmat, matmode,
                  d_px, d_pxp, d_dpx, d_dpxp, i, j, k, xres, yres, zres, dt, dx);
    } 
    else if (mattype==6) { //plrc model
	plrc_update(d_ex, d_ey, d_ez, d_hx, d_hy, d_hz, d_mat, d_mattype, d_mattab, nmat, matmode,
                  d_plrcx, d_plrcy, d_plrcz, i, j, k, xres, yres, zres, dt, dx);
    } 
    else { 
	    if (d_mat[pos]==0 && !(matmode==0 || matmode==2)) { /*= mattype 0 = linear material given pixel by pixel*/
		    sigma = d_sigma[pos];
		    epsilon = d_epsilon[pos];
	    } else if (mattype==1) { /*tabulated linear material, here should be also the cp3 and drude option*/
		    sigma = d_mattab[GETSIGMA(d_mat[pos])]; //sigma = d_mattab[GETSIGMA(mattype)];
		    epsilon = d_mattab[GETEPS(d_mat[pos])]; //epsilon = d_mattab[GETEPS(mattype)];        
	    } else {
		    sigma = 0;
		    epsilon = EPSILON_0;
	    }
	    ca  = (1 - sigma*dt/2/epsilon)/(1 + sigma*dt/2/epsilon);
	    cb  = (dt/epsilon/dx)/(1 + sigma*dt/2/epsilon);

            if (j>0 && k>0) 
   	       d_ex[pos] = ca*d_ex[pos] + cb*((hzvat - d_hz[CAT(i, j-1, k)]) -
			    (hyvat - d_hy[CAT(i, j, k-1)]));

            if (i>0 && k>0) 
	       d_ey[pos] = ca*d_ey[pos] + cb*((hxvat - d_hx[CAT(i, j, k-1)]) -
			    (hzvat - d_hz[CAT(i-1, j, k)]));

            if (i>0 && j>0) 
 	       d_ez[pos] = ca*d_ez[pos] + cb*((hyvat - d_hy[CAT(i-1, j, k)]) -
			    (hxvat - d_hx[CAT(i, j-1, k)]));

    }

    /*pec treatment*/
    if (i>0 && j>0 && k>0) {
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

	    if ((mattype_zm!=10 && mattype==10) || (mattype_zm==10 && mattype!=10))
	    {
		    d_ex[pos] = 0;
		    d_ey[pos] = 0;
	    }

	    if (mattype_xmym==10 && mattype!=10)
	    {
		    d_ez[pos] = 0;
	    }
	    if (mattype_xmzm==10 && mattype!=10)
	    {
		    d_ey[pos] = 0;
	    }
	    if (mattype_ymzm==10 && mattype!=10)
	    {
		    d_ex[pos] = 0;
	    }
    }


    /*sf source treatment*/
    if (d_sfset && (i>0 && j>0 && k>0)) { 
            theta = d_sfset[6];
            phi = d_sfset[7];
            psi = d_sfset[8];

	    d = k_dcomp(i, j, k,
			    xres, yres, zres,
			    theta, phi,    
			    0, xres, 0, yres, 0, zres);
	    ve = get_dval(d_sf_jpool_e, d_sf_jpool_size, d);

	    if ((mattype_xm!=10 && mattype==10) || (mattype_xm==10 && mattype!=10))
	    {
     
		    d_ey[pos] = -k_gey(ve, theta, phi, psi);
		    d_ez[pos] = -k_gez(ve, theta, phi, psi);
	    }

	    if ((mattype_ym!=10 && mattype==10) || (mattype_ym==10 && mattype!=10))
	    {
		    d_ex[pos] = -k_gex(ve, theta, phi, psi);
		    d_ez[pos] = -k_gez(ve, theta, phi, psi);
	    }

	    if ((mattype_zm!=10 && mattype==10) || (mattype_zm==10 && mattype!=10))
	    {
		    d_ex[pos] = -k_gex(ve, theta, phi, psi);
		    d_ey[pos] = -k_gey(ve, theta, phi, psi);
	    }

	    if (mattype_xmym==10 && mattype!=10)
	    {
		    d_ez[pos] = -k_gez(ve, theta, phi, psi);
	    }
	    if (mattype_xmzm==10 && mattype!=10)
	    {
		    d_ey[pos] = -k_gey(ve, theta, phi, psi);
	    }
	    if (mattype_ymzm==10 && mattype!=10)
	    {
		    d_ex[pos] = -k_gex(ve, theta, phi, psi);
	    }
    }



}

__global__  void
eKernel_cpml(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_epsilon, float *d_mu, float *d_sigma, float *d_sigast,
         int *bnds, int *depths,
         float *d_cpml_kappae_x0, float *d_cpml_kappae_xn, float *d_cpml_kappae_y0, float *d_cpml_kappae_yn, float *d_cpml_kappae_z0, float *d_cpml_kappae_zn,
         int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k;
    float ca, cb, hxvat, hyvat, hzvat, kappax, kappay, kappaz;

    long int pos;

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
 
    if (i<0 || j<0 || k<0) return;
    if (bnds[BNDX0] != SV_BOUNDARY_CPML && i==0) return;
    if (bnds[BNDY0] != SV_BOUNDARY_CPML && j==0) return;
    if (bnds[BNDZ0] != SV_BOUNDARY_CPML && k==0) return;

    if (i>=(xres) || j>=(yres) || k>=zres) return;

    pos = CAT(i, j, k);
   
    hxvat = d_hx[pos];
    hyvat = d_hy[pos];
    hzvat = d_hz[pos];
    
 
    ca  = (1 - d_sigma[pos]*dt/2/d_epsilon[pos])/(1 + d_sigma[pos]*dt/2/d_epsilon[pos]);
    cb  = (dt/d_epsilon[pos]/dx)/(1 + d_sigma[pos]*dt/2/d_epsilon[pos]);


    kappax = kappay = kappaz = 1;
    
    if (bnds[BNDX0] == SV_BOUNDARY_CPML && i<depths[BNDX0]) {
	    kappax = MAX(kappax, d_cpml_kappae_x0[i]);
    }
    if (bnds[BNDXN] == SV_BOUNDARY_CPML && i>=(xres - depths[BNDXN])) {
	    kappax = MAX(kappax, d_cpml_kappae_xn[i - (xres - depths[BNDXN])]);//
    }
    if (bnds[BNDY0] == SV_BOUNDARY_CPML && j<depths[BNDY0]) {
	    kappay = MAX(kappay, d_cpml_kappae_y0[j]);
    }
    if (bnds[BNDYN] == SV_BOUNDARY_CPML && j>=(yres - depths[BNDYN])) {
	    kappay = MAX(kappay, d_cpml_kappae_yn[j - (yres - depths[BNDYN])]);//
    }

    if (bnds[BNDZ0] == SV_BOUNDARY_CPML && k<depths[BNDZ0]) {
	    kappaz = MAX(kappaz, d_cpml_kappae_z0[k]);
    }
    if (bnds[BNDZN] == SV_BOUNDARY_CPML && k>=(zres - depths[BNDZN])) {
	    kappaz = MAX(kappaz, d_cpml_kappae_zn[k - (zres - depths[BNDZN])]);
    }

    if ((bnds[BNDX0] != SV_BOUNDARY_CPML && i>0 && j>0 && k>0) || (bnds[BNDX0] == SV_BOUNDARY_CPML && !(j==0 || k==0)))
       d_ex[pos] = ca*d_ex[pos] + cb*((hzvat - d_hz[CAT(i, j-1, k)])/kappay -
               (hyvat - d_hy[CAT(i, j, k-1)])/kappaz);

    if ((bnds[BNDY0] != SV_BOUNDARY_CPML && i>0 && j>0 && k>0) || (bnds[BNDY0] == SV_BOUNDARY_CPML && !(i==0 || k==0)))
       d_ey[pos] = ca*d_ey[pos] + cb*((hxvat - d_hx[CAT(i, j, k-1)])/kappaz -
               (hzvat - d_hz[CAT(i-1, j, k)])/kappax);

    if ((bnds[BNDZ0] != SV_BOUNDARY_CPML && i>0 && j>0 && k>0) || (bnds[BNDZ0] == SV_BOUNDARY_CPML && !(i==0 || j==0)))
       d_ez[pos] = ca*d_ez[pos] + cb*((hyvat - d_hy[CAT(i-1, j, k)])/kappax -
               (hxvat - d_hx[CAT(i, j-1, k)])/kappay);

}

__global__  void
eKernel_cpml_tab(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_epsilon, float *d_mu, float *d_sigma, float *d_sigast,
         int *d_mat, int *d_mattype, float *d_mattab, int nmat, int matmode,
         int *bnds, int *depths,
         float *d_cpml_kappae_x0, float *d_cpml_kappae_xn, float *d_cpml_kappae_y0, float *d_cpml_kappae_yn, float *d_cpml_kappae_z0, float *d_cpml_kappae_zn,
         int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int dir, float *d_sf_jpool_e, float *d_sfset, int d_sf_jpool_size,
         float *d_plrcx, float *d_plrcy, float *d_plrcz, float *d_exp, float *d_eyp, float *d_ezp,
         float *d_px, float *d_pxp, float *d_dpx, float *d_dpxp)
{
    int i, j, k, mattype;
    int mattype_xm, mattype_ym, mattype_zm, mattype_xmym, mattype_xmzm, mattype_ymzm;
    float ca, cb, hxvat, hyvat, hzvat, kappax, kappay, kappaz;
    float sigma, epsilon, theta, phi, psi, d, ve;

    long int pos;

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
 
    if (i<0 || j<0 || k<0) return;
    if (bnds[BNDX0] != SV_BOUNDARY_CPML && i==0) return;
    if (bnds[BNDY0] != SV_BOUNDARY_CPML && j==0) return;
    if (bnds[BNDZ0] != SV_BOUNDARY_CPML && k==0) return;

    if (i>=(xres) || j>=(yres) || k>=zres) return;

    pos = CAT(i, j, k);
    if (i>0 && j>0 && k>0) {
	    mattype_xm = d_mattype[d_mat[CAT(i-1, j, k)]];
	    mattype_ym = d_mattype[d_mat[CAT(i, j-1, k)]];
	    mattype_zm = d_mattype[d_mat[CAT(i, j, k-1)]];
	    mattype_xmym = d_mattype[d_mat[CAT(i-1, j-1, k)]];
	    mattype_xmzm = d_mattype[d_mat[CAT(i-1, j, k-1)]];
	    mattype_ymzm = d_mattype[d_mat[CAT(i, j-1, k-1)]];
    }

    mattype = d_mattype[d_mat[pos]];

    hxvat = d_hx[pos];
    hyvat = d_hy[pos];
    hzvat = d_hz[pos];
    
   if (mattype==2) {
       	drude_update(d_ex, d_ey, d_ez, d_hx, d_hy, d_hz, d_mat, d_mattype, d_mattab, nmat, matmode,
                  d_plrcx, d_plrcy, d_plrcz, i, j, k, xres, yres, zres, dt, dx); 
    }
    else if (mattype==4) { //cp model
	cp_update(d_ex, d_ey, d_ez, d_hx, d_hy, d_hz, d_mat, d_mattype, d_mattab, nmat, matmode,
                  d_plrcx, d_plrcy, d_plrcz, i, j, k, xres, yres, zres, dt, dx);
    } 
    else if (mattype==5) { //ade model
	ade_update(d_ex, d_ey, d_ez, d_exp, d_eyp, d_ezp, d_hx, d_hy, d_hz, d_mat, d_mattype, d_mattab, nmat, matmode,
                  d_px, d_pxp, d_dpx, d_dpxp, i, j, k, xres, yres, zres, dt, dx);
    } 
    else if (mattype==6) { //plrc model
	plrc_update(d_ex, d_ey, d_ez, d_hx, d_hy, d_hz, d_mat, d_mattype, d_mattab, nmat, matmode,
                  d_plrcx, d_plrcy, d_plrcz, i, j, k, xres, yres, zres, dt, dx);
    } 
    else {
       if (d_mat[pos]==0 && !(matmode==0 || matmode==2)) { /*= mattype 0 = linear material given pixel by pixel*/
            sigma = d_sigma[pos];
            epsilon = d_epsilon[pos];
       } else if (mattype==1) { /*tabulated linear material*/
            sigma = d_mattab[GETSIGMA(d_mat[pos])];
            epsilon = d_mattab[GETEPS(d_mat[pos])];        
       } else if (mattype==10) {
            sigma = 0;
            epsilon = EPSILON_0;
       }
       ca  = (1 - sigma*dt/2/epsilon)/(1 + sigma*dt/2/epsilon);
       cb  = (dt/epsilon/dx)/(1 + sigma*dt/2/epsilon);


       kappax = kappay = kappaz = 1;
    
       if (bnds[BNDX0] == SV_BOUNDARY_CPML && i<depths[BNDX0]) {
	    kappax = MAX(kappax, d_cpml_kappae_x0[i]);
       }
       if (bnds[BNDXN] == SV_BOUNDARY_CPML && i>=(xres - depths[BNDXN])) {
	    kappax = MAX(kappax, d_cpml_kappae_xn[i - (xres - depths[BNDXN])]);//
       }
       if (bnds[BNDY0] == SV_BOUNDARY_CPML && j<depths[BNDY0]) {
   	    kappay = MAX(kappay, d_cpml_kappae_y0[j]);
       }
       if (bnds[BNDYN] == SV_BOUNDARY_CPML && j>=(yres - depths[BNDYN])) {
	    kappay = MAX(kappay, d_cpml_kappae_yn[j - (yres - depths[BNDYN])]);//
       }

       if (bnds[BNDZ0] == SV_BOUNDARY_CPML && k<depths[BNDZ0]) {
	    kappaz = MAX(kappaz, d_cpml_kappae_z0[k]);
       }
       if (bnds[BNDZN] == SV_BOUNDARY_CPML && k>=(zres - depths[BNDZN])) {
	    kappaz = MAX(kappaz, d_cpml_kappae_zn[k - (zres - depths[BNDZN])]);
       }

       if ((bnds[BNDX0] != SV_BOUNDARY_CPML && i>0 && j>0 && k>0) || (bnds[BNDX0] == SV_BOUNDARY_CPML && !(j==0 || k==0)))
          d_ex[pos] = ca*d_ex[pos] + cb*((hzvat - d_hz[CAT(i, j-1, k)])/kappay -
               (hyvat - d_hy[CAT(i, j, k-1)])/kappaz);

       if ((bnds[BNDY0] != SV_BOUNDARY_CPML && i>0 && j>0 && k>0) || (bnds[BNDY0] == SV_BOUNDARY_CPML && !(i==0 || k==0)))
          d_ey[pos] = ca*d_ey[pos] + cb*((hxvat - d_hx[CAT(i, j, k-1)])/kappaz -
               (hzvat - d_hz[CAT(i-1, j, k)])/kappax);

       if ((bnds[BNDZ0] != SV_BOUNDARY_CPML && i>0 && j>0 && k>0) || (bnds[BNDZ0] == SV_BOUNDARY_CPML && !(i==0 || j==0)))
          d_ez[pos] = ca*d_ez[pos] + cb*((hyvat - d_hy[CAT(i-1, j, k)])/kappax -
               (hxvat - d_hx[CAT(i, j-1, k)])/kappay);
    }

    /*pec treatment*/
    if (i>0 && j>0 && k>0) { 
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

       if ((mattype_zm!=10 && mattype==10) || (mattype_zm==10 && mattype!=10))
       {
	    d_ex[pos] = 0;
	    d_ey[pos] = 0;
       }

       if (mattype_xmym==10 && mattype!=10)
       {
	    d_ez[pos] = 0;
       }
       if (mattype_xmzm==10 && mattype!=10)
       {
	    d_ey[pos] = 0;
       }
       if (mattype_ymzm==10 && mattype!=10)
       {
	    d_ex[pos] = 0;
       }
   }
    /*sf source treatment*/
    if (d_sfset && (i>0 && j>0 && k>0)) { 
            theta = d_sfset[6];
            phi = d_sfset[7];
            psi = d_sfset[8];

	    d = k_dcomp(i, j, k,
			    xres, yres, zres,
			    theta, phi,    
			    0, xres, 0, yres, 0, zres);
	    ve = get_dval(d_sf_jpool_e, d_sf_jpool_size, d);

	    if ((mattype_xm!=10 && mattype==10) || (mattype_xm==10 && mattype!=10))
	    {
     
		    d_ey[pos] = -k_gey(ve, theta, phi, psi);
		    d_ez[pos] = -k_gez(ve, theta, phi, psi);
	    }

	    if ((mattype_ym!=10 && mattype==10) || (mattype_ym==10 && mattype!=10))
	    {
		    d_ex[pos] = -k_gex(ve, theta, phi, psi);
		    d_ez[pos] = -k_gez(ve, theta, phi, psi);
	    }

	    if ((mattype_zm!=10 && mattype==10) || (mattype_zm==10 && mattype!=10))
	    {
		    d_ex[pos] = -k_gex(ve, theta, phi, psi);
		    d_ey[pos] = -k_gey(ve, theta, phi, psi);
	    }

	    if (mattype_xmym==10 && mattype!=10)
	    {
		    d_ez[pos] = -k_gez(ve, theta, phi, psi);
	    }
	    if (mattype_xmzm==10 && mattype!=10)
	    {
		    d_ey[pos] = -k_gey(ve, theta, phi, psi);
	    }
	    if (mattype_ymzm==10 && mattype!=10)
	    {
		    d_ex[pos] = -k_gex(ve, theta, phi, psi);
	    }
    }

 

}

__global__  void
eKernel_none(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         int xres, int yres, int zres, float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k;
    float cb, hxvat, hyvat, hzvat;
    long int pos;

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
 
    if (i<0 || j<0 || k<0) return;
    if (i>=(xres) || j>=(yres) || k>=zres) return;
 
    pos = CAT(i, j, k);
    hxvat = d_hx[pos];
    hyvat = d_hy[pos];
    hzvat = d_hz[pos];
    cb  = dt/EPSILON_0/dx;

    if (j>0 && k>0)
       d_ex[pos] = d_ex[pos] + cb*((hzvat - d_hz[CAT(i, j-1, k)]) -
               (hyvat - d_hy[CAT(i, j, k-1)]));

    if (i>0 && k>0)
       d_ey[pos] = d_ey[pos] + cb*((hxvat - d_hx[CAT(i, j, k-1)]) -
               (hzvat - d_hz[CAT(i-1, j, k)]));

    if (i>0 && j>0)
       d_ez[pos] = d_ez[pos] + cb*((hyvat - d_hy[CAT(i-1, j, k)]) -
               (hxvat - d_hx[CAT(i, j-1, k)]));

}
__global__  void
phKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_epsilon, float *d_mu, float *d_sigma, float *d_sigast, int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int dir, int nthreads, int tw)
{
    int i, j, k;
    float da, db, exvat, eyvat, ezvat, sigast, mu;
    long int pos;
    int zyres = zres*yres;

    pos = (blockIdx.x*tw + blockIdx.y)*nthreads + threadIdx.x; 

    i = (int)(pos/(zyres));
    j = (pos - i*zyres)/zres;
    k = pos - i*zyres - j*zres;

    if (i<0 || j<0 || k<0) return;
    if (i>(xres - 1) || j>(yres - 1) || k>zres - 1) return;

    exvat = d_ex[pos];
    eyvat = d_ey[pos];
    ezvat = d_ez[pos];
    sigast = d_sigast[pos];
    mu = d_mu[pos];

    da  = (1 - sigast*dt/2/mu)/(1 + sigast*dt/2/mu);
    db  = (dt/mu/dx)/(1 + sigast*dt/2/mu);

    if (j<(yres - 1) && k<(zres - 1))
       d_hx[pos] = da*d_hx[pos] + db*((d_ey[pos+1] - eyvat) -
              (d_ez[pos+zres] - ezvat));

    if (i<(xres - 1) && k<(zres - 1))
       d_hy[pos] = da*d_hy[pos] + db*((d_ez[pos+zyres] - ezvat) -
              (d_ex[pos+1] - exvat));

    if (i<(xres - 1) && j<(yres - 1))
       d_hz[pos] = da*d_hz[pos] + db*((d_ex[pos+zres] - exvat) -
              (d_ey[pos+zyres] - eyvat));

}

__global__  void
hKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_epsilon, float *d_mu, float *d_sigma, float *d_sigast, int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k;
    float da, db, exvat, eyvat, ezvat;
    long int pos;

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

    if (i<0 || j<0 || k<0) return;
    if (i>(xres - 1) || j>(yres - 1) || k>(zres-1)) return;

    pos = CAT(i, j, k);
    exvat = d_ex[pos];
    eyvat = d_ey[pos];
    ezvat = d_ez[pos];
    da  = (1 - d_sigast[pos]*dt/2/d_mu[pos])/(1 + d_sigast[pos]*dt/2/d_mu[pos]);
    db  = (dt/d_mu[pos]/dx)/(1 + d_sigast[pos]*dt/2/d_mu[pos]);

    if (j<(yres - 1) && k<(zres - 1)) 
       d_hx[pos] = da*d_hx[pos] + db*((d_ey[CAT(i, j, k+1)] - eyvat) -
              (d_ez[CAT(i, j+1, k)] - ezvat));

    if (i<(xres - 1) && k<(zres - 1))
       d_hy[pos] = da*d_hy[pos] + db*((d_ez[CAT(i+1, j, k)] - ezvat) -
              (d_ex[CAT(i, j, k+1)] - exvat));

    if (i<(xres - 1) && j<(yres - 1))
       d_hz[pos] = da*d_hz[pos] + db*((d_ex[CAT(i, j+1, k)] - exvat) -
              (d_ey[CAT(i+1, j, k)] - eyvat));

}

__global__  void
hKernel_tab(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_epsilon, float *d_mu, float *d_sigma, float *d_sigast, 
         int *d_mat, int *d_mattype, float *d_mattab, int nmat, int matmode,
         int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k, mattype;
    float da, db, exvat, eyvat, ezvat, mu, sigast;
    long int pos;

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

    if (i<0 || j<0 || k<0) return;
    if (i>(xres - 1) || j>(yres - 1) || k>(zres-1)) return;

    pos = CAT(i, j, k);

    mattype = d_mattype[d_mat[pos]];

    exvat = d_ex[pos];
    eyvat = d_ey[pos];
    ezvat = d_ez[pos];

    if (d_mat[pos]==0 && !(matmode==0 || matmode==1)) { /*= mattype 0 = linear material given pixel by pixel*/
	    sigast = d_sigast[pos];
	    mu = d_mu[pos];
    } else if (mattype==1 || mattype==2 || mattype==3 || mattype==4 || mattype==5 || mattype==6) { /*tabulated linear material, here should be also the cp3 and drude option*/
	    sigast = d_mattab[GETSIGAST(d_mat[pos])]; 
	    mu = d_mattab[GETMU(d_mat[pos])];
    } else {
	    sigast = 0;
	    mu = MU_0;
    }


    da  = (1 - sigast*dt/2/mu)/(1 + sigast*dt/2/mu);
    db  = (dt/mu/dx)/(1 + sigast*dt/2/mu);

    if (j<(yres - 1) && k<(zres - 1))  
       d_hx[pos] = da*d_hx[pos] + db*((d_ey[CAT(i, j, k+1)] - eyvat) -
              (d_ez[CAT(i, j+1, k)] - ezvat));

    if (i<(xres - 1) && k<(zres - 1))  
       d_hy[pos] = da*d_hy[pos] + db*((d_ez[CAT(i+1, j, k)] - ezvat) -
              (d_ex[CAT(i, j, k+1)] - exvat));

    if (i<(xres - 1) && j<(yres - 1))  
       d_hz[pos] = da*d_hz[pos] + db*((d_ex[CAT(i, j+1, k)] - exvat) -
              (d_ey[CAT(i+1, j, k)] - eyvat));

}
__global__  void
hKernel_cpml(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_epsilon, float *d_mu, float *d_sigma, float *d_sigast, int matmode,
         int *bnds, int *depths,
         float *d_cpml_kappah_x0, float *d_cpml_kappah_xn, float *d_cpml_kappah_y0, float *d_cpml_kappah_yn, float *d_cpml_kappah_z0, float *d_cpml_kappah_zn,
         int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k;
    float da, db, exvat, eyvat, ezvat, kappax, kappay, kappaz;
    long int pos;

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

    if (i<0 || j<0 || k<0) return;
    if (i>(xres - 1) || j>(yres - 1) || k>(zres-1)) return;

    
    kappax = kappay = kappaz = 1;
    
    if (bnds[BNDX0] == SV_BOUNDARY_CPML && i<depths[BNDX0]) {
	    kappax = MAX(kappax, d_cpml_kappah_x0[i]);
    }
    if (bnds[BNDXN] == SV_BOUNDARY_CPML && i>=(xres - depths[BNDXN])) {
	    kappax = MAX(kappax, d_cpml_kappah_xn[i - (xres - depths[BNDXN])]);//
    }
    if (bnds[BNDY0] == SV_BOUNDARY_CPML && j<depths[BNDY0]) {
	    kappay = MAX(kappay, d_cpml_kappah_y0[j]);
    }
    if (bnds[BNDYN] == SV_BOUNDARY_CPML && j>=(yres - depths[BNDYN])) {
	    kappay = MAX(kappay, d_cpml_kappah_yn[j - (yres - depths[BNDYN])]);//
    }

    if (bnds[BNDZ0] == SV_BOUNDARY_CPML && k<depths[BNDZ0]) {
	    kappaz = MAX(kappaz, d_cpml_kappah_z0[k]);
    }
    if (bnds[BNDZN] == SV_BOUNDARY_CPML && k>=(zres - depths[BNDZN])) {
	    kappaz = MAX(kappaz, d_cpml_kappah_zn[k - (zres - depths[BNDZN])]);
    }

    pos = CAT(i, j, k);
    exvat = d_ex[pos];
    eyvat = d_ey[pos];
    ezvat = d_ez[pos];

    if (matmode == 1 || matmode == 4) {
       da = 1;
       db = dt/MU_0/dx;
    } else {
       da  = (1 - d_sigast[pos]*dt/2/d_mu[pos])/(1 + d_sigast[pos]*dt/2/d_mu[pos]);
       db  = (dt/d_mu[pos]/dx)/(1 + d_sigast[pos]*dt/2/d_mu[pos]);
    }

    if (j<(yres - 1) && k<(zres - 1))  
       d_hx[pos] = da*d_hx[pos] + db*((d_ey[CAT(i, j, k+1)] - eyvat)/kappaz -
              (d_ez[CAT(i, j+1, k)] - ezvat)/kappay);

    if (i<(xres - 1) && k<(zres - 1))  
       d_hy[pos] = da*d_hy[pos] + db*((d_ez[CAT(i+1, j, k)] - ezvat)/kappax -
              (d_ex[CAT(i, j, k+1)] - exvat)/kappaz);

    if (i<(xres - 1) && j<(yres - 1))  
       d_hz[pos] = da*d_hz[pos] + db*((d_ex[CAT(i, j+1, k)] - exvat)/kappay -
              (d_ey[CAT(i+1, j, k)] - eyvat)/kappax);

}

__global__  void
hKernel_cpml_tab(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_epsilon, float *d_mu, float *d_sigma, float *d_sigast, int *d_mat, int *d_mattype, float *d_mattab, int nmat, int matmode,
         int *bnds, int *depths,
         float *d_cpml_kappah_x0, float *d_cpml_kappah_xn, float *d_cpml_kappah_y0, float *d_cpml_kappah_yn, float *d_cpml_kappah_z0, float *d_cpml_kappah_zn,
         int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k, mattype;
    float da, db, exvat, eyvat, ezvat, kappax, kappay, kappaz, mu, sigast;
    long int pos;

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

    if (i<0 || j<0 || k<0) return;
    if (i>(xres - 1) || j>(yres - 1) || k>(zres-1)) return;

    
    kappax = kappay = kappaz = 1;
    
    if (bnds[BNDX0] == SV_BOUNDARY_CPML && i<depths[BNDX0]) {
	    kappax = MAX(kappax, d_cpml_kappah_x0[i]);
    }
    if (bnds[BNDXN] == SV_BOUNDARY_CPML && i>=(xres - depths[BNDXN])) {
	    kappax = MAX(kappax, d_cpml_kappah_xn[i - (xres - depths[BNDXN])]);//
    }
    if (bnds[BNDY0] == SV_BOUNDARY_CPML && j<depths[BNDY0]) {
	    kappay = MAX(kappay, d_cpml_kappah_y0[j]);
    }
    if (bnds[BNDYN] == SV_BOUNDARY_CPML && j>=(yres - depths[BNDYN])) {
	    kappay = MAX(kappay, d_cpml_kappah_yn[j - (yres - depths[BNDYN])]);//
    }

    if (bnds[BNDZ0] == SV_BOUNDARY_CPML && k<depths[BNDZ0]) {
	    kappaz = MAX(kappaz, d_cpml_kappah_z0[k]);
    }
    if (bnds[BNDZN] == SV_BOUNDARY_CPML && k>=(zres - depths[BNDZN])) {
	    kappaz = MAX(kappaz, d_cpml_kappah_zn[k - (zres - depths[BNDZN])]);
    }

    pos = CAT(i, j, k);
    exvat = d_ex[pos];
    eyvat = d_ey[pos];
    ezvat = d_ez[pos];

    mattype = d_mattype[d_mat[pos]];

    if (d_mat[pos]==0 && !(matmode==0 || matmode==1)) { /*= mattype 0 = linear material given pixel by pixel*/
	    sigast = d_sigast[pos];
	    mu = d_mu[pos];
    } else if (mattype==1 || mattype==2 || mattype==3 || mattype==4 || mattype==5 || mattype==6) { /*tabulated linear material, here should be also the cp3 and drude option*/
            sigast = d_mattab[GETSIGAST(d_mat[pos])];
            mu = d_mattab[GETMU(d_mat[pos])];
    } else {
	    sigast = 0;
	    mu = MU_0;
    }

    da  = (1 - sigast*dt/2/mu)/(1 + sigast*dt/2/mu);
    db  = (dt/mu/dx)/(1 + sigast*dt/2/mu);

    if (j<(yres - 1) && k<(zres - 1))   
       d_hx[pos] = da*d_hx[pos] + db*((d_ey[CAT(i, j, k+1)] - eyvat)/kappaz -
              (d_ez[CAT(i, j+1, k)] - ezvat)/kappay);

    if (i<(xres - 1) && k<(zres - 1))
       d_hy[pos] = da*d_hy[pos] + db*((d_ez[CAT(i+1, j, k)] - ezvat)/kappax -
              (d_ex[CAT(i, j, k+1)] - exvat)/kappaz);

    if (i<(xres - 1) && j<(yres - 1))
       d_hz[pos] = da*d_hz[pos] + db*((d_ex[CAT(i, j+1, k)] - exvat)/kappay -
              (d_ey[CAT(i+1, j, k)] - eyvat)/kappax);

}

__global__  void
hKernel_none(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k;
    float db, exvat, eyvat, ezvat;
    long int pos;

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

    if (i<0 || j<0 || k<0) return;
    if (i>(xres - 1) || j>(yres - 1) || k>(zres-1)) return;

    pos = CAT(i, j, k);
    exvat = d_ex[pos];
    eyvat = d_ey[pos];
    ezvat = d_ez[pos];
    db  = dt/MU_0/dx;

    if (j<(yres - 1) && k<(zres - 1))
       d_hx[pos] = d_hx[pos] + db*((d_ey[CAT(i, j, k+1)] - eyvat) -
              (d_ez[CAT(i, j+1, k)] - ezvat));

    if (i<(xres - 1) && k<(zres - 1))
       d_hy[pos] = d_hy[pos] + db*((d_ez[CAT(i+1, j, k)] - ezvat) -
              (d_ex[CAT(i, j, k+1)] - exvat));

    if (i<(xres - 1) && j<(yres - 1))
       d_hz[pos] = d_hz[pos] + db*((d_ex[CAT(i, j+1, k)] - exvat) -
              (d_ey[CAT(i+1, j, k)] - eyvat));

}
/*
__global__  void
yKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_epsilon, float *d_mu, float *d_sigma, float *d_sigast, int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k;
    float da, db, exvat, eyvat, ezvat;
    float ca, cb, hxvat, hyvat, hzvat;
    long int pos;
 
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

    if (i>=0 && j>=0 && k>=0 && i<(xres-1) && j<(yres-1) && k<(zres-1)) 
    {
        exvat = d_ex[pos];
        eyvat = d_ey[pos];
        ezvat = d_ez[pos];
        da  = (1 - d_sigast[pos]*dt/2/d_mu[pos])/(1 + d_sigast[pos]*dt/2/d_mu[pos]);
        db  = (dt/d_mu[pos]/dx)/(1 + d_sigast[pos]*dt/2/d_mu[pos]);

        d_hx[pos] = hxvat = da*d_hx[pos] + db*((d_ey[CAT(i, j, k+1)] - eyvat) -
              (d_ez[CAT(i, j+1, k)] - ezvat));

        d_hy[pos] = hyvat = da*d_hy[pos] + db*((d_ez[CAT(i+1, j, k)] - ezvat) -
              (d_ex[CAT(i, j, k+1)] - exvat));

        d_hz[pos] = hzvat = da*d_hz[pos] + db*((d_ex[CAT(i, j+1, k)] - exvat) -
              (d_ey[CAT(i+1, j, k)] - eyvat));
    }
    __syncthreads();

    if (i>0 && j>0 && k>0 && i<(xres) && j<(yres) && k<(zres)) 
    {
        //hxvat = d_hx[pos];
        //hyvat = d_hy[pos];
        //hzvat = d_hz[pos];
        ca  = (1 - d_sigma[pos]*dt/2/d_epsilon[pos])/(1 + d_sigma[pos]*dt/2/d_epsilon[pos]);
        cb  = (dt/d_epsilon[pos]/dx)/(1 + d_sigma[pos]*dt/2/d_epsilon[pos]);

        d_ex[pos] = ca*d_ex[pos] + cb*((hzvat - d_hz[CAT(i, j-1, k)]) -
               (hyvat - d_hy[CAT(i, j, k-1)]));

        d_ey[pos] = ca*d_ey[pos] + cb*((hxvat - d_hx[CAT(i, j, k-1)]) -
               (hzvat - d_hz[CAT(i-1, j, k)]));

        d_ez[pos] = ca*d_ez[pos] + cb*((hyvat - d_hy[CAT(i-1, j, k)]) -
               (hxvat - d_hx[CAT(i, j-1, k)]));
    }
}
*/
__global__  void
liaocpyKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
         float *d_x0, float *d_xn, float *d_y0, float *d_yn, float *d_z0, float *d_zn,
         int xres, int yres, int zres, int* bnds,
         float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k;
    long int pos;
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

    if (bnds[BNDX0] == 2) {
        if (i==0) {
            d_x0[CLX(EX, j, k)] = d_ex[pos];
            d_x0[CLX(EY, j, k)] = d_ey[pos];
            d_x0[CLX(EZ, j, k)] = d_ez[pos];
            d_x0[CLX(HX, j, k)] = d_hx[pos];
            d_x0[CLX(HY, j, k)] = d_hy[pos];
            d_x0[CLX(HZ, j, k)] = d_hz[pos];
        } else if (i==1) {
            d_x0[CLXP(EX, j, k)] = d_ex[pos];
            d_x0[CLXP(EY, j, k)] = d_ey[pos];
            d_x0[CLXP(EZ, j, k)] = d_ez[pos];
            d_x0[CLXP(HX, j, k)] = d_hx[pos];
            d_x0[CLXP(HY, j, k)] = d_hy[pos];
            d_x0[CLXP(HZ, j, k)] = d_hz[pos];
        }
    }

    if (bnds[BNDXN] == 2) {
        if (i==(xres-1)) {
            d_xn[CLX(EX, j, k)] = d_ex[pos];
            d_xn[CLX(EY, j, k)] = d_ey[pos];
            d_xn[CLX(EZ, j, k)] = d_ez[pos];
            d_xn[CLX(HX, j, k)] = d_hx[pos];
            d_xn[CLX(HY, j, k)] = d_hy[pos];
            d_xn[CLX(HZ, j, k)] = d_hz[pos];
        } else if (i==(xres-2)) {
            d_xn[CLXP(EX, j, k)] = d_ex[pos];
            d_xn[CLXP(EY, j, k)] = d_ey[pos];
            d_xn[CLXP(EZ, j, k)] = d_ez[pos];
            d_xn[CLXP(HX, j, k)] = d_hx[pos];
            d_xn[CLXP(HY, j, k)] = d_hy[pos];
            d_xn[CLXP(HZ, j, k)] = d_hz[pos];
        }
    }

    if (bnds[BNDY0] == 2) {
        if (j==0) {
            d_y0[CLY(EX, i, k)] = d_ex[pos];
            d_y0[CLY(EY, i, k)] = d_ey[pos];
            d_y0[CLY(EZ, i, k)] = d_ez[pos];
            d_y0[CLY(HX, i, k)] = d_hx[pos];
            d_y0[CLY(HY, i, k)] = d_hy[pos];
            d_y0[CLY(HZ, i, k)] = d_hz[pos];
        } else if (j==1) {
            d_y0[CLYP(EX, i, k)] = d_ex[pos];
            d_y0[CLYP(EY, i, k)] = d_ey[pos];
            d_y0[CLYP(EZ, i, k)] = d_ez[pos];
            d_y0[CLYP(HX, i, k)] = d_hx[pos];
            d_y0[CLYP(HY, i, k)] = d_hy[pos];
            d_y0[CLYP(HZ, i, k)] = d_hz[pos];
        }
    }

    if (bnds[BNDYN] == 2) {
        if (j==(yres-1)) {
            d_yn[CLY(EX, i, k)] = d_ex[pos];
            d_yn[CLY(EY, i, k)] = d_ey[pos];
            d_yn[CLY(EZ, i, k)] = d_ez[pos];
            d_yn[CLY(HX, i, k)] = d_hx[pos];
            d_yn[CLY(HY, i, k)] = d_hy[pos];
            d_yn[CLY(HZ, i, k)] = d_hz[pos];
        } else if (j==(yres-2)) {
            d_yn[CLYP(EX, i, k)] = d_ex[pos];
            d_yn[CLYP(EY, i, k)] = d_ey[pos];
            d_yn[CLYP(EZ, i, k)] = d_ez[pos];
            d_yn[CLYP(HX, i, k)] = d_hx[pos];
            d_yn[CLYP(HY, i, k)] = d_hy[pos];
            d_yn[CLYP(HZ, i, k)] = d_hz[pos];
        }
    }

    if (bnds[BNDZ0] == 2)
    {
        if (k==0) {
            d_z0[CLZ(EX, i, j)] = d_ex[pos];
            d_z0[CLZ(EY, i, j)] = d_ey[pos];
            d_z0[CLZ(EZ, i, j)] = d_ez[pos];
            d_z0[CLZ(HX, i, j)] = d_hx[pos];
            d_z0[CLZ(HY, i, j)] = d_hy[pos];
            d_z0[CLZ(HZ, i, j)] = d_hz[pos];
        } else if (k==1) {
            d_z0[CLZP(EX, i, j)] = d_ex[pos];
            d_z0[CLZP(EY, i, j)] = d_ey[pos];
            d_z0[CLZP(EZ, i, j)] = d_ez[pos];
            d_z0[CLZP(HX, i, j)] = d_hx[pos];
            d_z0[CLZP(HY, i, j)] = d_hy[pos];
            d_z0[CLZP(HZ, i, j)] = d_hz[pos];
        } 
    }

    if (bnds[BNDZN] == 2) {
        if (k==(zres-1)) {
            d_zn[CLZ(EX, i, j)] = d_ex[pos];
            d_zn[CLZ(EY, i, j)] = d_ey[pos];
            d_zn[CLZ(EZ, i, j)] = d_ez[pos];
            d_zn[CLZ(HX, i, j)] = d_hx[pos];
            d_zn[CLZ(HY, i, j)] = d_hy[pos];
            d_zn[CLZ(HZ, i, j)] = d_hz[pos];
        } else if (k==(zres-2)) {
            d_zn[CLZP(EX, i, j)] = d_ex[pos];
            d_zn[CLZP(EY, i, j)] = d_ey[pos];
            d_zn[CLZP(EZ, i, j)] = d_ez[pos];
            d_zn[CLZP(HX, i, j)] = d_hx[pos];
            d_zn[CLZP(HY, i, j)] = d_hy[pos];
            d_zn[CLZP(HZ, i, j)] = d_hz[pos];
        }
    }



}

__global__  void
liaorunKernelx(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz, float *d_epsilon,
         float *d_x0, float *d_xn, float *d_y0, float *d_yn, float *d_z0, float *d_zn, 
         int xres, int yres, int zres, int* bnds, int matmode,
         float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k;
    float lssx, ind;
    float mx;
    long int pos;
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

    if (matmode==1 || matmode == 3)
        ind = sqrt(d_epsilon[pos]/EPSILON_0);
    else ind = 1;

    lssx = (dt*LS/ind-dx)/(dt*LS/ind+dx);
    mx = (float)lssx;

    if (bnds[BNDX0] == 2) {
        if (i==0) {
            d_ex[pos] = d_x0[CLXP(EX, j, k)] + mx*(d_ex[CAT(i+1, j, k)] - d_x0[CLX(EX, j, k)]);
            d_ey[pos] = d_x0[CLXP(EY, j, k)] + mx*(d_ey[CAT(i+1, j, k)] - d_x0[CLX(EY, j, k)]);
            d_ez[pos] = d_x0[CLXP(EZ, j, k)] + mx*(d_ez[CAT(i+1, j, k)] - d_x0[CLX(EZ, j, k)]);
        }
    }
    if (bnds[BNDXN] == 2) {
        if (i==(xres-1)) {
            d_ex[pos] = d_xn[CLXP(EX, j, k)] + mx*(d_ex[CAT(i-1, j, k)] - d_xn[CLX(EX, j, k)]);
            d_ey[pos] = d_xn[CLXP(EY, j, k)] + mx*(d_ey[CAT(i-1, j, k)] - d_xn[CLX(EY, j, k)]);
            d_ez[pos] = d_xn[CLXP(EZ, j, k)] + mx*(d_ez[CAT(i-1, j, k)] - d_xn[CLX(EZ, j, k)]);
            d_hx[pos] = d_xn[CLXP(HX, j, k)] + mx*(d_hx[CAT(i-1, j, k)] - d_xn[CLX(HX, j, k)]);
            d_hy[pos] = d_xn[CLXP(HY, j, k)] + mx*(d_hy[CAT(i-1, j, k)] - d_xn[CLX(HY, j, k)]);
            d_hz[pos] = d_xn[CLXP(HZ, j, k)] + mx*(d_hz[CAT(i-1, j, k)] - d_xn[CLX(HZ, j, k)]);
        }
    }

}
__global__  void
liaorunKernely(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz, float *d_epsilon,
         float *d_x0, float *d_xn, float *d_y0, float *d_yn, float *d_z0, float *d_zn, 
         int xres, int yres, int zres, int* bnds, int matmode,
         float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k;
    float lssy, ind;
    float my;
    long int pos;
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

    if (matmode==1 || matmode == 3)
        ind = sqrt(d_epsilon[pos]/EPSILON_0);
    else ind = 1;

    lssy = (dt*LS/ind-dy)/(dt*LS/ind+dy);
    my = (float)lssy;

    if (bnds[BNDY0] == 2) {
        if (j==0) {
            d_ex[pos] = d_y0[CLYP(EX, i, k)] + my*(d_ex[CAT(i, j+1, k)] - d_y0[CLY(EX, i, k)]);
            d_ey[pos] = d_y0[CLYP(EY, i, k)] + my*(d_ey[CAT(i, j+1, k)] - d_y0[CLY(EY, i, k)]);
            d_ez[pos] = d_y0[CLYP(EZ, i, k)] + my*(d_ez[CAT(i, j+1, k)] - d_y0[CLY(EZ, i, k)]);
        }
    }
    if (bnds[BNDYN] == 2) {
        if (j==(yres-1)) {
            d_ex[pos] = d_yn[CLYP(EX, i, k)] + my*(d_ex[CAT(i, j-1, k)] - d_yn[CLY(EX, i, k)]);
            d_ey[pos] = d_yn[CLYP(EY, i, k)] + my*(d_ey[CAT(i, j-1, k)] - d_yn[CLY(EY, i, k)]);
            d_ez[pos] = d_yn[CLYP(EZ, i, k)] + my*(d_ez[CAT(i, j-1, k)] - d_yn[CLY(EZ, i, k)]);
            d_hx[pos] = d_yn[CLYP(HX, i, k)] + my*(d_hx[CAT(i, j-1, k)] - d_yn[CLY(HX, i, k)]);
            d_hy[pos] = d_yn[CLYP(HY, i, k)] + my*(d_hy[CAT(i, j-1, k)] - d_yn[CLY(HY, i, k)]);
            d_hz[pos] = d_yn[CLYP(HZ, i, k)] + my*(d_hz[CAT(i, j-1, k)] - d_yn[CLY(HZ, i, k)]);
        }
    }


}

__global__  void
liaorunKernelz(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz, float *d_epsilon,

         float *d_x0, float *d_xn, float *d_y0, float *d_yn, float *d_z0, float *d_zn, 
         int xres, int yres, int zres, int* bnds, int matmode,
         float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k;
    float lssz, ind;
    float /*mx, my,*/ mz;
    long int pos;
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

    if (matmode==1 || matmode == 3)
        ind = sqrt(d_epsilon[pos]/EPSILON_0);
    else ind = 1;

    lssz = (dt*LS/ind-dz)/(dt*LS/ind+dz);
    /*mx = (float)lssx;*/
    /*my = (float)lssy;*/
    mz = (float)lssz;

    if (bnds[BNDZ0] == 2) {
        if (k==0) {
            d_ex[pos] = d_z0[CLZP(EX, i, j)] + mz*(d_ex[CAT(i, j, k+1)] - d_z0[CLZ(EX, i, j)]);
            d_ey[pos] = d_z0[CLZP(EY, i, j)] + mz*(d_ey[CAT(i, j, k+1)] - d_z0[CLZ(EY, i, j)]);
            d_ez[pos] = d_z0[CLZP(EZ, i, j)] + mz*(d_ez[CAT(i, j, k+1)] - d_z0[CLZ(EZ, i, j)]);
        }
    }
    if (bnds[BNDZN] == 2) {
        if (k==(zres-1)) {
            d_ex[pos] = d_zn[CLZP(EX, i, j)] + mz*(d_ex[CAT(i, j, k-1)] - d_zn[CLZ(EX, i, j)]);
            d_ey[pos] = d_zn[CLZP(EY, i, j)] + mz*(d_ey[CAT(i, j, k-1)] - d_zn[CLZ(EY, i, j)]);
            d_ez[pos] = d_zn[CLZP(EZ, i, j)] + mz*(d_ez[CAT(i, j, k-1)] - d_zn[CLZ(EZ, i, j)]);
            d_hx[pos] = d_zn[CLZP(HX, i, j)] + mz*(d_hx[CAT(i, j, k-1)] - d_zn[CLZ(HX, i, j)]);
            d_hy[pos] = d_zn[CLZP(HY, i, j)] + mz*(d_hy[CAT(i, j, k-1)] - d_zn[CLZ(HY, i, j)]);
            d_hz[pos] = d_zn[CLZP(HZ, i, j)] + mz*(d_hz[CAT(i, j, k-1)] - d_zn[CLZ(HZ, i, j)]);
        }
    }

}

__global__  void
cpmlhKernelx(float *d_ey, float *d_ez, float *d_hy, float *d_hz, float *d_mu, float *d_sigast, int matmode,
            float *d_cpml_p_x0, float *d_cpml_p_xn,
            float *d_cpml_bh_x0, float *d_cpml_ch_x0, float *d_cpml_bh_xn, float *d_cpml_ch_xn,
            int xres, int yres, int zres, int *bnds, int *depths,
            float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k, ir;
    float db, bh, ch;
    long int pos;
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

    if (i<0 || j<0 || k<0) return;
    if (i>=(xres - 1) || j>=(yres - 1) || k>=(zres-1)) return;

    if (matmode == 1 || matmode == 4) db = dt/MU_0/dx;
    else db  = (dt/d_mu[pos]/dx)/(1 + d_sigast[pos]*dt/2/d_mu[pos]);


    if (bnds[BNDX0] == SV_BOUNDARY_CPML && (i>=0 &&  i<(depths[BNDX0] - 1) && j>=1 && j<(yres-1) && k>=1 && k<(zres-1)))
    {
        bh = d_cpml_bh_x0[i];
        ch = d_cpml_ch_x0[i];

        d_cpml_p_x0[CX0AT(i,j,k) + CX0_HZX] = bh*d_cpml_p_x0[CX0AT(i,j,k) + CX0_HZX]
            + ch*(d_ey[CAT(i+1, j, k)] - d_ey[pos])/dx;

        d_cpml_p_x0[CX0AT(i,j,k) + CX0_HYX] = bh*d_cpml_p_x0[CX0AT(i,j,k) + CX0_HYX]
            + ch*(d_ez[CAT(i+1, j, k)] - d_ez[pos])/dx;

        d_hy[pos] += dx*db*d_cpml_p_x0[CX0AT(i,j,k) + CX0_HYX]; //-
        d_hz[pos] -= dx*db*d_cpml_p_x0[CX0AT(i,j,k) + CX0_HZX]; //+
    }
    if (bnds[BNDXN] == SV_BOUNDARY_CPML && (i>=(xres - depths[BNDXN] + 1) && i<(xres-1) && j>=1 && j<(yres-1) && k>=1 && k<(zres-1)))
    {
        ir = i - (xres - depths[BNDXN]);
        bh = d_cpml_bh_xn[ir];
        ch = d_cpml_ch_xn[ir];

        d_cpml_p_xn[CXNAT(ir,j,k) + CXN_HZX] = bh*d_cpml_p_xn[CXNAT(ir,j,k) + CXN_HZX]
            + ch*(d_ey[CAT(i+1, j, k)] - d_ey[pos])/dx;

        d_cpml_p_xn[CXNAT(ir,j,k) + CXN_HYX] = bh*d_cpml_p_xn[CXNAT(ir,j,k) + CXN_HYX]
            + ch*(d_ez[CAT(i+1, j, k)] - d_ez[pos])/dx;

        d_hy[pos] += dx*db*d_cpml_p_xn[CXNAT(ir,j,k) + CXN_HYX]; //-
        d_hz[pos] -= dx*db*d_cpml_p_xn[CXNAT(ir,j,k) + CXN_HZX]; //+

    }
}

__global__  void
cpmlhKernely(float *d_ex, float *d_ez, float *d_hx, float *d_hz, float *d_mu, float *d_sigast, int matmode,
            float *d_cpml_p_y0, float *d_cpml_p_yn,
            float *d_cpml_bh_y0, float *d_cpml_ch_y0, float *d_cpml_bh_yn, float *d_cpml_ch_yn,
            int xres, int yres, int zres, int *bnds, int *depths,
            float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k, jr;
    float db, bh, ch;
    long int pos;
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

    if (i<0 || j<0 || k<0) return;
    if (i>=(xres - 1) || j>=(yres - 1) || k>=(zres-1)) return;

    if (matmode == 1 || matmode == 4) db = dt/MU_0/dx;
    else
    db  = (dt/d_mu[pos]/dx)/(1 + d_sigast[pos]*dt/2/d_mu[pos]);


    if (bnds[BNDY0] == SV_BOUNDARY_CPML && (j>=0 &&  j<(depths[BNDY0] - 1) && i>=1 && i<(xres-1) && k>=1 && k<(zres-1)))
    {
        bh = d_cpml_bh_y0[j];
        ch = d_cpml_ch_y0[j];

        d_cpml_p_y0[CY0AT(i,j,k) + CY0_HXY] = bh*d_cpml_p_y0[CY0AT(i,j,k) + CY0_HXY]
            + ch*(d_ez[CAT(i, j+1, k)] - d_ez[pos])/dx;

        d_cpml_p_y0[CY0AT(i,j,k) + CY0_HZY] = bh*d_cpml_p_y0[CY0AT(i,j,k) + CY0_HZY]
            + ch*(d_ex[CAT(i, j+1, k)] - d_ex[pos])/dx;

        d_hx[pos] -= dx*db*d_cpml_p_y0[CY0AT(i,j,k) + CY0_HXY]; 
        d_hz[pos] += dx*db*d_cpml_p_y0[CY0AT(i,j,k) + CY0_HZY]; 
    }
    if (bnds[BNDYN] == SV_BOUNDARY_CPML && (j>=(yres - depths[BNDYN] + 1) && j<(yres-1) && i>=1 && i<(xres-1) && k>=1 && k<(zres-1)))
    {
        jr = j - (yres - depths[BNDYN]);
        bh = d_cpml_bh_yn[jr];
        ch = d_cpml_ch_yn[jr];

        d_cpml_p_yn[CYNAT(i,jr,k) + CYN_HXY] = bh*d_cpml_p_yn[CYNAT(i,jr,k) + CYN_HXY]
            + ch*(d_ez[CAT(i, j+1, k)] - d_ez[pos])/dx;

        d_cpml_p_yn[CYNAT(i,jr,k) + CYN_HZY] = bh*d_cpml_p_yn[CYNAT(i,jr,k) + CYN_HZY]
            + ch*(d_ex[CAT(i, j+1, k)] - d_ex[pos])/dx;

        d_hx[pos] -= dx*db*d_cpml_p_yn[CYNAT(i,jr,k) + CYN_HXY]; 
        d_hz[pos] += dx*db*d_cpml_p_yn[CYNAT(i,jr,k) + CYN_HZY]; 

    }
}
__global__  void
cpmlhKernelz(float *d_ex, float *d_ey, float *d_hx, float *d_hy, float *d_mu, float *d_sigast, int matmode,
            float *d_cpml_p_z0, float *d_cpml_p_zn,
            float *d_cpml_bh_z0, float *d_cpml_ch_z0, float *d_cpml_bh_zn, float *d_cpml_ch_zn,
            int xres, int yres, int zres, int *bnds, int *depths,
            float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k, kr;
    float db, bh, ch;
    long int pos;
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

    if (i<0 || j<0 || k<0) return;
    if (i>=(xres - 1) || j>=(yres - 1) || k>=(zres-1)) return;

    if (matmode == 1 || matmode == 4) db = dt/MU_0/dx;
    else
    db  = (dt/d_mu[pos]/dx)/(1 + d_sigast[pos]*dt/2/d_mu[pos]);


    if (bnds[BNDZ0] == SV_BOUNDARY_CPML && (k>=0 &&  k<(depths[BNDZ0] - 1) && i>=1 && i<(xres-1) && j>=1 && j<(yres-1)))
    {
        bh = d_cpml_bh_z0[k];
        ch = d_cpml_ch_z0[k];

        d_cpml_p_z0[CZ0AT(i,j,k) + CZ0_HYZ] = bh*d_cpml_p_z0[CZ0AT(i,j,k) + CZ0_HYZ]
            + ch*(d_ex[CAT(i, j, k+1)] - d_ex[pos])/dx;

        d_cpml_p_z0[CZ0AT(i,j,k) + CZ0_HXZ] = bh*d_cpml_p_z0[CZ0AT(i,j,k) + CZ0_HXZ]
            + ch*(d_ey[CAT(i, j, k+1)] - d_ey[pos])/dx;

        d_hy[pos] -= dx*db*d_cpml_p_z0[CZ0AT(i,j,k) + CZ0_HYZ]; //-
        d_hx[pos] += dx*db*d_cpml_p_z0[CZ0AT(i,j,k) + CZ0_HXZ]; //+
    }
    if (bnds[BNDZN] == SV_BOUNDARY_CPML && (k>=(zres - depths[BNDZN] + 1) && k<(zres-1) && i>=1 && i<(xres-1) && j>=1 && j<(yres-1)))
    {
        kr = k - (zres - depths[BNDZN]);
        bh = d_cpml_bh_zn[kr];
        ch = d_cpml_ch_zn[kr];

        d_cpml_p_zn[CZNAT(i,j,kr) + CZN_HYZ] = bh*d_cpml_p_zn[CZNAT(i,j,kr) + CZN_HYZ]
            + ch*(d_ex[CAT(i, j, k+1)] - d_ex[pos])/dx;

        d_cpml_p_zn[CZNAT(i,j,kr) + CZN_HXZ] = bh*d_cpml_p_zn[CZNAT(i,j,kr) + CZN_HXZ]
            + ch*(d_ey[CAT(i, j, k+1)] - d_ey[pos])/dx;

        d_hy[pos] -= dx*db*d_cpml_p_zn[CZNAT(i,j,kr) + CZN_HYZ]; //-
        d_hx[pos] += dx*db*d_cpml_p_zn[CZNAT(i,j,kr) + CZN_HXZ]; //+

    }
}

__global__  void
cpmleKernelx(float *d_ey, float *d_ez, float *d_hy, float *d_hz, float *d_epsilon, float *d_sigma,
            int *d_mat, int *d_mattype, float *d_mattab, int nmat,
            float *d_cpml_p_x0, float *d_cpml_p_xn,
            float *d_cpml_be_x0, float *d_cpml_ce_x0, float *d_cpml_be_xn, float *d_cpml_ce_xn,
            int xres, int yres, int zres, int *bnds, int *depths,
            float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k, ir;
    int mattype, mat;
    float cb, be, ce, sigma, epsilon;
    long int pos;
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

    if (nmat>0) {
	    mat = d_mat[pos];
	    mattype = d_mattype[d_mat[pos]];

	    if (mat==0) { /*= mattype 0 = linear material given pixel by pixel*/
		    sigma = d_sigma[pos];
		    epsilon = d_epsilon[pos];
	    } else if (mattype==1) { /*tabulated linear material*/
		    sigma = d_mattab[GETSIGMA(mat)];
		    epsilon = d_mattab[GETEPS(mat)];
	    } else if (mattype==10) {
		    sigma = 0;
		    epsilon = EPSILON_0;
	    } else {
		    sigma = d_sigma[pos];
		    epsilon = d_epsilon[pos];

	    }
    } else {
        sigma = d_sigma[pos];
        epsilon = d_epsilon[pos];

    }

    cb  = (dt/epsilon/dx)/(1 + sigma*dt/2/epsilon);


    if (bnds[BNDX0] == SV_BOUNDARY_CPML && (i>=1 &&  i<(depths[BNDX0]) && j>=0 && j<(yres) && k>=0 && k<(zres)))
    {
        be = d_cpml_be_x0[i];
        ce = d_cpml_ce_x0[i];

        d_cpml_p_x0[CX0AT(i,j,k) + CX0_EYX] = be*d_cpml_p_x0[CX0AT(i,j,k) + CX0_EYX]
            + ce*(d_hz[pos] - d_hz[CAT(i-1, j, k)])/dx;

        d_cpml_p_x0[CX0AT(i,j,k) + CX0_EZX] = be*d_cpml_p_x0[CX0AT(i,j,k) + CX0_EZX]
            + ce*(d_hy[pos] - d_hy[CAT(i-1, j, k)])/dx;

        if (!(j==0 || k==0))
        d_ey[pos] -= dx*cb*d_cpml_p_x0[CX0AT(i,j,k) + CX0_EYX]; 

        if (!(j==0 || k==0))
        d_ez[pos] += dx*cb*d_cpml_p_x0[CX0AT(i,j,k) + CX0_EZX]; 
    }
    if (bnds[BNDXN] == SV_BOUNDARY_CPML && (i>=(xres - depths[BNDXN] + 1) &&  i<(xres) && j>=0 && j<(yres) && k>=0 && k<(zres)))
    {
        ir = i - (xres - depths[BNDXN]);
        be = d_cpml_be_xn[ir];
        ce = d_cpml_ce_xn[ir];

        d_cpml_p_xn[CXNAT(ir,j,k) + CXN_EYX] = be*d_cpml_p_xn[CXNAT(ir,j,k) + CXN_EYX]
            + ce*(d_hz[pos] - d_hz[CAT(i-1, j, k)])/dx;

        d_cpml_p_xn[CXNAT(ir,j,k) + CXN_EZX] = be*d_cpml_p_xn[CXNAT(ir,j,k) + CXN_EZX]
            + ce*(d_hy[pos] - d_hy[CAT(i-1, j, k)])/dx;

        if (!(j==0 || k==0))
        d_ey[pos] -= dx*cb*d_cpml_p_xn[CXNAT(ir,j,k) + CXN_EYX]; 

        if (!(j==0 || k==0))
        d_ez[pos] += dx*cb*d_cpml_p_xn[CXNAT(ir,j,k) + CXN_EZX]; 
    }
}

__global__  void
cpmleKernely(float *d_ex, float *d_ez, float *d_hx, float *d_hz, float *d_epsilon, float *d_sigma,
            int *d_mat, int *d_mattype, float *d_mattab, int nmat,
            float *d_cpml_p_y0, float *d_cpml_p_yn,
            float *d_cpml_be_y0, float *d_cpml_ce_y0, float *d_cpml_be_yn, float *d_cpml_ce_yn,
            int xres, int yres, int zres, int *bnds, int *depths,
            float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k, jr;
    int mattype, mat;
    float cb, be, ce, sigma, epsilon;
    long int pos;
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
    if (nmat>0) {
	    mat = d_mat[pos];
	    mattype = d_mattype[d_mat[pos]];

	    if (mat==0) { /*= mattype 0 = linear material given pixel by pixel*/
		    sigma = d_sigma[pos];
		    epsilon = d_epsilon[pos];
	    } else if (mattype==1) { /*tabulated linear material*/
		    sigma = d_mattab[GETSIGMA(mat)];
		    epsilon = d_mattab[GETEPS(mat)];
	    } else if (mattype==10) {
		    sigma = 0;
		    epsilon = EPSILON_0;
	    } else {
		    sigma = d_sigma[pos];
		    epsilon = d_epsilon[pos];

	    }
    } else {
        sigma = d_sigma[pos];
        epsilon = d_epsilon[pos];

    }

    cb  = (dt/epsilon/dx)/(1 + sigma*dt/2/epsilon);

    //cb  = (dt/d_epsilon[pos]/dx)/(1 + d_sigma[pos]*dt/2/d_epsilon[pos]);

    if (bnds[BNDY0] == SV_BOUNDARY_CPML && (j>=1 &&  j<(depths[BNDY0]) && i>=0 && i<(xres) && k>=0 && k<(zres)))
    {
        be = d_cpml_be_y0[j];
        ce = d_cpml_ce_y0[j];

        d_cpml_p_y0[CY0AT(i,j,k) + CY0_EXY] = be*d_cpml_p_y0[CY0AT(i,j,k) + CY0_EXY]
            + ce*(d_hz[pos] - d_hz[CAT(i, j-1, k)])/dx;

        d_cpml_p_y0[CY0AT(i,j,k) + CY0_EZY] = be*d_cpml_p_y0[CY0AT(i,j,k) + CY0_EZY]
            + ce*(d_hx[pos] - d_hx[CAT(i, j-1, k)])/dx;

        if (!(j==0 || k==0))
        d_ex[pos] += dx*cb*d_cpml_p_y0[CY0AT(i,j,k) + CY0_EXY]; 

        if (!(j==0 || k==0))
        d_ez[pos] -= dx*cb*d_cpml_p_y0[CY0AT(i,j,k) + CY0_EZY]; 
    }
    if (bnds[BNDYN] == SV_BOUNDARY_CPML && (j>=(yres - depths[BNDYN] + 1) &&  j<(yres) && i>=0 && i<(xres) && k>=0 && k<(zres)))
    {
        jr = j - (yres - depths[BNDYN]);
        be = d_cpml_be_yn[jr];
        ce = d_cpml_ce_yn[jr];

        d_cpml_p_yn[CYNAT(i,jr,k) + CYN_EXY] = be*d_cpml_p_yn[CYNAT(i,jr,k) + CYN_EXY]
            + ce*(d_hz[pos] - d_hz[CAT(i, j-1, k)])/dx;

        d_cpml_p_yn[CYNAT(i,jr,k) + CYN_EZY] = be*d_cpml_p_yn[CYNAT(i,jr,k) + CYN_EZY]
            + ce*(d_hx[pos] - d_hx[CAT(i, j-1, k)])/dx;

        if (!(i==0 || k==0))
        d_ex[pos] += dx*cb*d_cpml_p_yn[CYNAT(i,jr,k) + CYN_EXY]; 

        if (!(i==0 || k==0))
        d_ez[pos] -= dx*cb*d_cpml_p_yn[CYNAT(i,jr,k) + CYN_EZY]; 
    }
}
__global__  void
cpmleKernelz(float *d_ex, float *d_ey, float *d_hx, float *d_hy, float *d_epsilon, float *d_sigma,
            int *d_mat, int *d_mattype, float *d_mattab, int nmat,
            float *d_cpml_p_z0, float *d_cpml_p_zn,
            float *d_cpml_be_z0, float *d_cpml_ce_z0, float *d_cpml_be_zn, float *d_cpml_ce_zn,
            int xres, int yres, int zres, int *bnds, int *depths,
            float dx, float dy, float dz, float dt, int dir)
{
    int i, j, k, kr;
    int mattype, mat;
    float cb, be, ce, sigma, epsilon;
    long int pos;
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
    if (nmat>0) {
	    mat = d_mat[pos];
	    mattype = d_mattype[d_mat[pos]];

	    if (mat==0) { //= mattype 0 = linear material given pixel by pixel
		    sigma = d_sigma[pos];
		    epsilon = d_epsilon[pos];
	    } else if (mattype==1) { //tabulated linear material
		    sigma = d_mattab[GETSIGMA(mat)];
		    epsilon = d_mattab[GETEPS(mat)];
	    } else if (mattype==10) {
		    sigma = 0;
		    epsilon = EPSILON_0;
	    } else {
		    sigma = d_sigma[pos];
		    epsilon = d_epsilon[pos];

	    }
    } else {
        sigma = d_sigma[pos];
        epsilon = d_epsilon[pos];

    }

    cb  = (dt/epsilon/dx)/(1 + sigma*dt/2/epsilon);

    //cb  = (dt/d_epsilon[pos]/dx)/(1 + d_sigma[pos]*dt/2/d_epsilon[pos]);

    if (bnds[BNDZ0] == SV_BOUNDARY_CPML && (k>=1 &&  k<(depths[BNDZ0]) && i>=0 && i<(xres) && j>=0 && j<(yres)))
    {
        be = d_cpml_be_z0[k];
        ce = d_cpml_ce_z0[k];

        d_cpml_p_z0[CZ0AT(i,j,k) + CZ0_EYZ] = be*d_cpml_p_z0[CZ0AT(i,j,k) + CZ0_EYZ]
            + ce*(d_hx[pos] - d_hx[CAT(i, j, k-1)])/dx;

        d_cpml_p_z0[CZ0AT(i,j,k) + CZ0_EXZ] = be*d_cpml_p_z0[CZ0AT(i,j,k) + CZ0_EXZ]
            + ce*(d_hy[pos] - d_hy[CAT(i, j, k-1)])/dx;

        if (!(i==0 || j==0))
        d_ey[pos] += dx*cb*d_cpml_p_z0[CZ0AT(i,j,k) + CZ0_EYZ]; 

        if (!(i==0 || j==0))
        d_ex[pos] -= dx*cb*d_cpml_p_z0[CZ0AT(i,j,k) + CZ0_EXZ]; 
    }
    if (bnds[BNDZN] == SV_BOUNDARY_CPML && (k>=(zres - depths[BNDZN] + 1) &&  k<(zres) && i>=0 && i<(xres) && j>=0 && j<(yres)))
    {
        kr = k - (zres - depths[BNDZN]);
        be = d_cpml_be_zn[kr];
        ce = d_cpml_ce_zn[kr];

        d_cpml_p_zn[CZNAT(i,j,kr) + CZN_EYZ] = be*d_cpml_p_zn[CZNAT(i,j,kr) + CZN_EYZ]
            + ce*(d_hx[pos] - d_hx[CAT(i, j, k-1)])/dx;

        d_cpml_p_zn[CZNAT(i,j,kr) + CZN_EXZ] = be*d_cpml_p_zn[CZNAT(i,j,kr) + CZN_EXZ]
            + ce*(d_hy[pos] - d_hy[CAT(i, j, k-1)])/dx;

        if (!(i==0 || j==0))
        d_ey[pos] += dx*cb*d_cpml_p_zn[CZNAT(i,j,kr) + CZN_EYZ]; 

        if (!(i==0 || j==0))
        d_ex[pos] -= dx*cb*d_cpml_p_zn[CZNAT(i,j,kr) + CZN_EXZ]; 
    }
}



__global__  void
srcepointKernel(float *d_ex, float *d_ey, float *d_ez, int xres, int yres, int zres, int ipos, int jpos, int kpos,
         float ex, float ey, float ez, int dir)
{
    int i, j, k;
    long int pos;
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
    
    if (i==ipos && j==jpos && k==kpos) {  
       pos = CAT(i, j, k);
       if (ex != 0) d_ex[pos] = ex;
       if (ey != 0) d_ey[pos] = ey;
       if (ez != 0) d_ez[pos] = ez;
    }

}
__global__  void
srchpointKernel(float *d_hx, float *d_hy, float *d_hz, int xres, int yres, int zres, int ipos, int jpos, int kpos,
         float hx, float hy, float hz, int dir)
{
    int i, j, k;
    long int pos;
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
    if (i==ipos && j==jpos && k==kpos) {  
       pos = CAT(i, j, k);
       if (hx != 0) d_hx[pos] = hx;
       if (hy != 0) d_hy[pos] = hy;
       if (hz != 0) d_hz[pos] = hz;
    }

}

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




// iterovat pres blok s cachovanim hodnot do __shared__ tak jak doted,
// ale inkrementovat soucasne nthreads farfield dat
// (aby kazde vlakno puvodne urcene pro konkretni bod v prostoru inkremenovalo jiny farfield bod)
// melo byt to tedy fungovat dobre az do celkoveho poctu yres*zres farfield bodu.


/*__global__  void
pffKernelx(float *d_ex, float *d_ey, float *d_ez, 
         int* iset, float *d_ff_ex, float *d_ff_ey, float *d_ff_ez,
         int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int step, int k, int nthreads)
{
    int i, j, nn, npnt, ndat, shift, istart, ii, swap;
    long int pos;
    float x, y, z;// devx, devy, devz, devxd, devyd, devzd;
    float r, theta, size, a, b, c;
    float dvadt = 0.5/dt, lsdt = LS*dt;;  
    bool skipi0, skipin;

    extern __shared__ float s[];

    j = threadIdx.x; //only one block computing all for now, more could be used only with helpers

    if (k<0 || k>=(zres) || j<0 || j>=(yres)) return;
    if (k<iset[K0] || k>=iset[K1]) return;

    skipi0 = skipin = 0;
    if (iset[SKIPI0] && iset[SKIPI0_JMIN]<=iset[J0] && iset[SKIPI0_KMIN]<=iset[K0] && iset[SKIPI0_JMAX]>=iset[J1] && iset[SKIPI0_KMAX]>=iset[K1]) skipi0 = 1;
    if (iset[SKIPIN] && iset[SKIPIN_JMIN]<=iset[J0] && iset[SKIPIN_KMIN]<=iset[K0] && iset[SKIPIN_JMAX]>=iset[J1] && iset[SKIPIN_KMAX]>=iset[K1]) skipin = 1;

    i=iset[I0];

    if (j<yres) {
	    pos = CAT(i, j, k);
	    s[6*j+0] = d_ex[pos];
	    s[6*j+1] = d_ey[pos];
	    s[6*j+2] = d_ez[pos];
	    s[6*j+3] = (d_ex[CAT(i + 1,j,k)] - d_ex[CAT(i - 1,j,k)])/2/dx;
	    s[6*j+4] = (d_ey[CAT(i + 1,j,k)] - d_ey[CAT(i - 1,j,k)])/2/dy;
	    s[6*j+5] = (d_ez[CAT(i + 1,j,k)] - d_ez[CAT(i - 1,j,k)])/2/dz;
    }

    __syncthreads();

    //find position of this point and other properties
    //iterate over x planes
    size = dx*dz;

    if (!skipi0) {
	    for (swap=0; swap<nthreads; swap++) { //swap so every pixel commits to every far field point, yres = nthreads //was yres

		    npnt = threadIdx.x - swap; //this will by default run by 0...n //was threads = yres
		    if (npnt<0) npnt += nthreads;  //wrap around //was +=yres

		    //all this works only up to nthreads far field points

		    if (npnt<0 || npnt>=iset[NPNTS]) continue;

		    ndat = iset[NPNTS + 1 + 5*npnt];
		    istart = iset[NPNTS + 2 + 5*npnt];
		    x = (float)iset[NPNTS + 3 + 5*npnt];
		    y = (float)iset[NPNTS + 4 + 5*npnt];
		    z = (float)iset[NPNTS + 5 + 5*npnt];
		    shift = 0;
		    for (ii=0; ii<npnt; ii++) 
			    shift += iset[NPNTS + 1 + 5*ii];

             
		    if (j>=iset[J0] && j<iset[J1]) {
			    if (!(iset[SKIPI0] && j>=iset[SKIPI0_JMIN] && k>=iset[SKIPI0_KMIN] && j<=iset[SKIPI0_JMAX] && k<=iset[SKIPI0_KMAX]))
			    {
				    r = k_dist(x, y, z, i, j, k);
				    theta = JCPI*size*(i-x)/r;
                                    r *= dx;
				    nn = shift - istart + (int)(step + 1 + r/lsdt);
				    a = JCPI/r*size;
				    b = -theta/(r*r);
				    c = -theta/(LS*r)*dvadt;

				    if (nn < 3 || (nn-shift)>=(ndat-3)) {
					    continue;
				    }

				    d_ff_ex[nn] += -c*s[6*j+0];
				    d_ff_ex[nn - 1] += -a*s[6*j+3] + b*s[6*j+0];
				    d_ff_ex[nn - 2] += c*s[6*j+0];

				    d_ff_ey[nn] += -c*s[6*j+1];
				    d_ff_ey[nn - 1] += -a*s[6*j+4] + b*s[6*j+1];
				    d_ff_ey[nn - 2] += c*s[6*j+1];

				    d_ff_ez[nn] += -c*s[6*j+2]; 
				    d_ff_ez[nn - 1] += -a*s[6*j+5] + b*s[6*j+2];
				    d_ff_ez[nn - 2] += c*s[6*j+2];
			    }
		    }
                   __syncthreads();
	    }
    }
    //FIXME up to this point it was tested (so it produces same results as convetional approach)
    // and optimized, but still it is slower than convetional approach (even with no helpers), at least for more than 200 farfield points


       i=iset[I1];
       pos = CAT(i, j, k);
       s[6*j+0] = d_ex[pos];
       s[6*j+1] = d_ey[pos];
       s[6*j+2] = d_ez[pos];
       s[6*j+3] = (d_ex[CAT(i + 1,j,k)] - d_ex[CAT(i - 1,j,k)])/2/dx;
       s[6*j+4] = (d_ey[CAT(i + 1,j,k)] - d_ey[CAT(i - 1,j,k)])/2/dy;
       s[6*j+5] = (d_ez[CAT(i + 1,j,k)] - d_ez[CAT(i - 1,j,k)])/2/dz;


       __syncthreads();


       if (!skipin) {
       npnt = 
       for (npnt=0; npnt<iset[NPNTS]; npnt++) {

       ndat = iset[NPNTS + 1 + 5*npnt];
       istart = iset[NPNTS + 2 + 5*npnt];
       x = (float)iset[NPNTS + 3 + 5*npnt];
       y = (float)iset[NPNTS + 4 + 5*npnt];
       z = (float)iset[NPNTS + 5 + 5*npnt];
       shift = 0;
       for (ii=0; ii<npnt; ii++) 
       shift += iset[NPNTS + 1 + 5*ii];

       if (j>=iset[J0] && j<iset[J1]) {

       if (!(iset[SKIPIN] && j>=iset[SKIPIN_JMIN] && k>=iset[SKIPIN_KMIN] && j<=iset[SKIPIN_JMAX] && k<=iset[SKIPIN_KMAX]))
       { 

       r = k_dist(x, y, z, i, j, k)*dx;
       theta = k_angle(i, j, k, x, y, z, 3);
       nn = shift - istart + (int)(step + 1 + r/(LS*dt));
       a = JCPI/r*size;
       b = -JCPI*theta/(r*r)*size;
       c = -JCPI*theta/(LS*r)*size*dvadt;
       if (nn < 3 || (nn - shift)>=(ndat-3)) {
       continue;
       }

       d_ff_ex[nn] += -c*s[6*j+0];
       d_ff_ex[nn - 1] += a*s[6*j+3] + b*s[6*j+0];
       d_ff_ex[nn - 2] += c*s[6*j+0];

       d_ff_ey[nn] += -c*s[6*j+1];
       d_ff_ey[nn - 1] += a*s[6*j+4] + b*s[6*j+1];
       d_ff_ey[nn - 2] += c*s[6*j+1];

				    d_ff_ez[nn] += -c*s[6*j+2]; 
				    d_ff_ez[nn - 1] += a*s[6*j+5] + b*s[6*j+2];
				    d_ff_ez[nn - 2] += c*s[6*j+2];
			    }
	    }
    }
}
*/
/*
__global__  void
pffKernely(float *d_ex, float *d_ey, float *d_ez, 
         int* iset, float *d_ff_ex, float *d_ff_ey, float *d_ff_ez,
         int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int step, int nthreads)
{
    int i, j, k, nn, npnt, ndat, shift, istart, ii;
    long int pos;
    float x, y, z, devx, devy, devz, devxd, devyd, devzd;
    float r, theta, size, a, b, c;
    float dvadt = 0.5/dt;  
    bool skipj0, skipjn;
    
    npnt = blockIdx.x*nthreads + threadIdx.x;
    k = npnt/xres;
    i = npnt - k*xres;
    if (k<0 || k>=(zres) || i<0 || i>=(xres)) return;

    skipj0 = skipjn = 0;
    if (iset[SKIPJ0] && iset[SKIPJ0_IMIN]<=iset[I0] && iset[SKIPJ0_KMIN]<=iset[K0] && iset[SKIPJ0_IMAX]>=iset[I1] && iset[SKIPJ0_KMAX]>=iset[K1]) skipj0 = 1;
    if (iset[SKIPJN] && iset[SKIPJN_IMIN]<=iset[I0] && iset[SKIPJN_KMIN]<=iset[K0] && iset[SKIPJN_IMAX]>=iset[I1] && iset[SKIPJN_KMAX]>=iset[K1]) skipjn = 1;

    //find position of this point and other properties
    //iterate over x planes
    size = dx*dz;

    if (!skipj0) {
	    j=iset[J0];

	    pos = CAT(i, j, k);
	    devx = d_ex[pos];
	    devy = d_ey[pos];
	    devz = d_ez[pos];
	    devxd = (d_ex[CAT(i, j + 1,k)] - d_ex[CAT(i, j - 1,k)])/2/dx;
	    devyd = (d_ey[CAT(i, j + 1,k)] - d_ey[CAT(i, j - 1,k)])/2/dy;
	    devzd = (d_ez[CAT(i, j + 1,k)] - d_ez[CAT(i, j - 1,k)])/2/dz;

	    for (npnt=0; npnt<iset[NPNTS]; npnt++) {
		    ndat = iset[NPNTS + 1 + 5*npnt];
		    istart = iset[NPNTS + 2 + 5*npnt];
		    x = (float)iset[NPNTS + 3 + 5*npnt];
		    y = (float)iset[NPNTS + 4 + 5*npnt];
		    z = (float)iset[NPNTS + 5 + 5*npnt];
		    shift = 0;
		    for (ii=0; ii<npnt; ii++) 
			    shift += iset[NPNTS + 1 + 5*ii];

	            if (!(iset[SKIPJ0] && i>=iset[SKIPJ0_IMIN] && k>=iset[SKIPJ0_KMIN] && i<=iset[SKIPJ0_IMAX] && k<=iset[SKIPJ0_KMAX]))
		    {
			    r = k_dist(x, y, z, i, j, k)*dx;
			    theta = k_angle(i, j, k, x, y, z, 0); 
			    nn = shift - istart + (int)(step + 1 + r/(LS*dt));
			    a = JCPI/r*size;
			    b = -JCPI*theta/(r*r)*size;
			    c = -JCPI*theta/(LS*r)*size*dvadt;
			    if (nn < 3 || (nn-shift)>=(ndat-3)) {
				    continue;
			    }
			    d_ff_ex[nn] += -c*devx;
			    d_ff_ex[nn - 1] += -a*devxd + b*devx;
			    d_ff_ex[nn - 2] += c*devx;

			    d_ff_ey[nn] += -c*devy;
			    d_ff_ey[nn - 1] += -a*devyd + b*devy;
			    d_ff_ey[nn - 2] += c*devy;

			    d_ff_ez[nn] += -c*devz; 
			    d_ff_ez[nn - 1] += -a*devzd + b*devz;
			    d_ff_ez[nn - 2] += c*devz;
		    }
	    }
    }
    if (!skipjn) {
	    j=iset[J1];
	    pos = CAT(i, j, k);
	    devx = d_ex[pos];
	    devy = d_ey[pos];
	    devz = d_ez[pos];
	    devxd = (d_ex[CAT(i, j + 1,k)] - d_ex[CAT(i, j - 1,k)])/2/dx;
	    devyd = (d_ey[CAT(i, j + 1,k)] - d_ey[CAT(i, j - 1,k)])/2/dy;
	    devzd = (d_ez[CAT(i, j + 1,k)] - d_ez[CAT(i, j - 1,k)])/2/dz;

	    for (npnt=0; npnt<iset[NPNTS]; npnt++) {
		    ndat = iset[NPNTS + 1 + 5*npnt];
		    istart = iset[NPNTS + 2 + 5*npnt];
		    x = (float)iset[NPNTS + 3 + 5*npnt];
		    y = (float)iset[NPNTS + 4 + 5*npnt];
		    z = (float)iset[NPNTS + 5 + 5*npnt];
		    shift = 0;
		    for (ii=0; ii<npnt; ii++) 
			    shift += iset[NPNTS + 1 + 5*ii];

                    if (!(iset[SKIPJN] && i>=iset[SKIPJN_IMIN] && k>=iset[SKIPJN_KMIN] && i<=iset[SKIPJN_IMAX] && k<=iset[SKIPJN_KMAX]))
		    { 
			    r = k_dist(x, y, z, i, j, k)*dx;
			    theta = k_angle(i, j, k, x, y, z, 1);
			    nn = shift - istart + (int)(step + 1 + r/(LS*dt));
			    a = JCPI/r*size;
			    b = -JCPI*theta/(r*r)*size;
			    c = -JCPI*theta/(LS*r)*size*dvadt;
			    if (nn < 3 || (nn - shift)>=(ndat-3)) {
				    continue;
			    }
			    d_ff_ex[nn] += -c*devx;
			    d_ff_ex[nn - 1] += a*devxd + b*devx;
			    d_ff_ex[nn - 2] += c*devx;

			    d_ff_ey[nn] += -c*devy;
			    d_ff_ey[nn - 1] += a*devyd + b*devy;
			    d_ff_ey[nn - 2] += c*devy;

			    d_ff_ez[nn] += -c*devz; 
			    d_ff_ez[nn - 1] += a*devzd + b*devz;
			    d_ff_ez[nn - 2] += c*devz;
		    }
	    }
    }
}


__global__  void
pffKernelz(float *d_ex, float *d_ey, float *d_ez, 
         int* iset, float *d_ff_ex, float *d_ff_ey, float *d_ff_ez,
         int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int step, int nthreads)
{
    int i, j, k, nn, npnt, ndat, shift, istart, ii;
    long int pos;
    float x, y, z, devx, devy, devz, devxd, devyd, devzd;
    float r, theta, size, a, b, c;
    float dvadt = 0.5/dt;  
    bool skipk0, skipkn;

    npnt = blockIdx.x*nthreads + threadIdx.x;
    j = npnt/xres;
    i = npnt - j*xres;
    if (i<0 || i>=(xres) || j<0 || j>=(yres)) return;


    skipk0 = skipkn = 0;
    if (iset[SKIPK0] && iset[SKIPK0_IMIN]<=iset[I0] && iset[SKIPK0_JMIN]<=iset[J0] && iset[SKIPK0_IMAX]>=iset[I1] && iset[SKIPK0_JMAX]>=iset[J1]) skipk0 = 1;
    if (iset[SKIPKN] && iset[SKIPKN_IMIN]<=iset[I0] && iset[SKIPKN_JMIN]<=iset[J0] && iset[SKIPKN_IMAX]>=iset[I1] && iset[SKIPKN_JMAX]>=iset[J1]) skipkn = 1;

    //find position of this point and other properties
    //iterate over x planes
    size = dx*dz;

    if (!skipk0) {
	    k=iset[K0];

	    pos = CAT(i, j, k);
	    devx = d_ex[pos];
	    devy = d_ey[pos];
	    devz = d_ez[pos];
	    devxd = (d_ex[CAT(i,j,k+1)] - d_ex[CAT(i,j,k-1)])/2/dx;
	    devyd = (d_ey[CAT(i,j,k+1)] - d_ey[CAT(i,j,k-1)])/2/dy;
	    devzd = (d_ez[CAT(i,j,k+1)] - d_ez[CAT(i,j,k-1)])/2/dz;

	    for (npnt=0; npnt<iset[NPNTS]; npnt++) {
		    ndat = iset[NPNTS + 1 + 5*npnt];
		    istart = iset[NPNTS + 2 + 5*npnt];
		    x = (float)iset[NPNTS + 3 + 5*npnt];
		    y = (float)iset[NPNTS + 4 + 5*npnt];
		    z = (float)iset[NPNTS + 5 + 5*npnt];
		    shift = 0;
		    for (ii=0; ii<npnt; ii++) 
			    shift += iset[NPNTS + 1 + 5*ii];

                    if (!(iset[SKIPK0] && i>=iset[SKIPK0_IMIN] && j>=iset[SKIPK0_JMIN] && i<=iset[SKIPK0_IMAX] && j<=iset[SKIPK0_JMAX]))
		    {
			    r = k_dist(x, y, z, i, j, k)*dx;
			    theta = k_angle(i, j, k, x, y, z, 4); 
			    nn = shift - istart + (int)(step + 1 + r/(LS*dt));
			    a = JCPI/r*size;
			    b = -JCPI*theta/(r*r)*size;
			    c = -JCPI*theta/(LS*r)*size*dvadt;
			    if (nn < 3 || (nn-shift)>=(ndat-3)) {
				    continue;
			    }
			    d_ff_ex[nn] += -c*devx;
			    d_ff_ex[nn - 1] += -a*devxd + b*devx;
			    d_ff_ex[nn - 2] += c*devx;

			    d_ff_ey[nn] += -c*devy;
			    d_ff_ey[nn - 1] += -a*devyd + b*devy;
			    d_ff_ey[nn - 2] += c*devy;

			    d_ff_ez[nn] += -c*devz; 
			    d_ff_ez[nn - 1] += -a*devzd + b*devz;
			    d_ff_ez[nn - 2] += c*devz;
		    }
	    }
    }
    if (!skipkn) {
	    k=iset[K1];
	    pos = CAT(i, j, k);
	    devx = d_ex[pos];
	    devy = d_ey[pos];
	    devz = d_ez[pos];
	    devxd = (d_ex[CAT(i,j,k+1)] - d_ex[CAT(i,j,k-1)])/2/dx;
	    devyd = (d_ey[CAT(i,j,k+1)] - d_ey[CAT(i,j,k-1)])/2/dy;
	    devzd = (d_ez[CAT(i,j,k+1)] - d_ez[CAT(i,j,k-1)])/2/dz;

	    for (npnt=0; npnt<iset[NPNTS]; npnt++) {
		    ndat = iset[NPNTS + 1 + 5*npnt];
		    istart = iset[NPNTS + 2 + 5*npnt];
		    x = (float)iset[NPNTS + 3 + 5*npnt];
		    y = (float)iset[NPNTS + 4 + 5*npnt];
		    z = (float)iset[NPNTS + 5 + 5*npnt];
		    shift = 0;
		    for (ii=0; ii<npnt; ii++) 
			    shift += iset[NPNTS + 1 + 5*ii];

                    if (!(iset[SKIPKN] && i>=iset[SKIPKN_IMIN] && j>=iset[SKIPKN_JMIN] && i<=iset[SKIPKN_IMAX] && j<=iset[SKIPKN_JMAX]))
		    { 
			    r = k_dist(x, y, z, i, j, k)*dx;
			    theta = k_angle(i, j, k, x, y, z, 5);
			    nn = shift - istart + (int)(step + 1 + r/(LS*dt));
			    a = JCPI/r*size;
			    b = -JCPI*theta/(r*r)*size;
			    c = -JCPI*theta/(LS*r)*size*dvadt;
			    if (nn < 3 || (nn - shift)>=(ndat-3)) {
				    continue;
			    }
			    d_ff_ex[nn] += -c*devx;
			    d_ff_ex[nn - 1] += a*devxd + b*devx;
			    d_ff_ex[nn - 2] += c*devx;

			    d_ff_ey[nn] += -c*devy;
			    d_ff_ey[nn - 1] += a*devyd + b*devy;
			    d_ff_ey[nn - 2] += c*devy;

			    d_ff_ez[nn] += -c*devz; 
			    d_ff_ez[nn - 1] += a*devzd + b*devz;
			    d_ff_ez[nn - 2] += c*devz;
		    }
	    }
    }
}
*/

__global__  void
ffKernel(float *d_ex, float *d_ey, float *d_ez, 
         int* iset, float *d_ff_ex, float *d_ff_ey, float *d_ff_ez,
         int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int step, int nthreads)
{
    int i, j, k,  nn, npnt, ndat, shift, istart;
    long int pos;
    float x, y, z;
    float r, theta, size, a, b, c, dev, lsdt;
    float dvadt = 0.5/dt;  
    bool skipi0, skipin, skipj0, skipjn, skipk0, skipkn;
    int zyres = zres*yres;
    
    //get far field point to be computed
    //npnt = threadIdx.x;
    npnt = blockIdx.x*nthreads + threadIdx.x;
    if (npnt>=iset[NPNTS]) return; 

    skipi0 = skipin = skipj0 = skipjn = skipk0 = skipkn = 0;
    if (iset[SKIPI0] && iset[SKIPI0_JMIN]<=iset[J0] && iset[SKIPI0_KMIN]<=iset[K0] && iset[SKIPI0_JMAX]>=iset[J1] && iset[SKIPI0_KMAX]>=iset[K1]) skipi0 = 1;
    if (iset[SKIPIN] && iset[SKIPIN_JMIN]<=iset[J0] && iset[SKIPIN_KMIN]<=iset[K0] && iset[SKIPIN_JMAX]>=iset[J1] && iset[SKIPIN_KMAX]>=iset[K1]) skipin = 1;
    if (iset[SKIPJ0] && iset[SKIPJ0_IMIN]<=iset[I0] && iset[SKIPJ0_KMIN]<=iset[K0] && iset[SKIPJ0_IMAX]>=iset[I1] && iset[SKIPJ0_KMAX]>=iset[K1]) skipj0 = 1;
    if (iset[SKIPJN] && iset[SKIPJN_IMIN]<=iset[I0] && iset[SKIPJN_KMIN]<=iset[K0] && iset[SKIPJN_IMAX]>=iset[I1] && iset[SKIPJN_KMAX]>=iset[K1]) skipjn = 1;
    if (iset[SKIPK0] && iset[SKIPK0_IMIN]<=iset[I0] && iset[SKIPK0_JMIN]<=iset[J0] && iset[SKIPK0_IMAX]>=iset[I1] && iset[SKIPK0_JMAX]>=iset[J1]) skipk0 = 1;
    if (iset[SKIPKN] && iset[SKIPKN_IMIN]<=iset[I0] && iset[SKIPKN_JMIN]<=iset[J0] && iset[SKIPKN_IMAX]>=iset[I1] && iset[SKIPKN_JMAX]>=iset[J1]) skipkn = 1;


    //find position of this point and other properties
    ndat = iset[NPNTS + 1 + 5*npnt];
    istart = iset[NPNTS + 2 + 5*npnt];
    x = (float)iset[NPNTS + 3 + 5*npnt];
    y = (float)iset[NPNTS + 4 + 5*npnt];
    z = (float)iset[NPNTS + 5 + 5*npnt];

    //determine shift in accumulator for this point
    shift = 0;
    for (i=0; i<npnt; i++) 
       shift += iset[NPNTS + 1 + 5*i];

    lsdt = LS*dt;
    //iterate over x planes
    size = dy*dz;
    
    if (!skipi0) {
	    for (j=iset[J0]; j<iset[J1]; j++)
	    {
		    for (k=iset[K0]; k<iset[K1]; k++)
		    {
			    if (!(iset[SKIPI0] && j>=iset[SKIPI0_JMIN] && k>=iset[SKIPI0_KMIN] && j<=iset[SKIPI0_JMAX] && k<=iset[SKIPI0_KMAX]))
			    {
				    i=iset[I0];
				    r = k_dist(x, y, z, i, j, k);
				    theta = JCPI*size*(i-x)/r;
                                    r *= dx;

				    nn = shift - istart + (int)(step + 1 + r/lsdt);
				    a = JCPI/r*size/2.0;
				    b = -theta/(r*r);
				    c = -theta/(LS*r)*dvadt;
				    if (nn < 3 || (nn-shift)>=(ndat-3)) {
					    continue;
				    }
				    pos = CAT(i, j, k);
				    dev = d_ex[pos];

				    d_ff_ex[nn] += -c*dev;
				    d_ff_ex[nn - 1] += -a*((d_ex[pos+zyres] - d_ex[pos-zyres])/dx) + b*dev;
				    d_ff_ex[nn - 2] += c*dev;

				    dev = d_ey[pos];
				    d_ff_ey[nn] += -c*dev;
				    d_ff_ey[nn - 1] += -a*((d_ey[pos+zyres] - d_ey[pos-zyres])/dy) + b*dev;
				    d_ff_ey[nn - 2] += c*dev;

				    dev = d_ez[pos];
				    d_ff_ez[nn] += -c*dev; 
				    d_ff_ez[nn - 1] += -a*((d_ez[pos+zyres] - d_ez[pos-zyres])/dz) + b*dev;
				    d_ff_ez[nn - 2] += c*dev;
			    }
		    }
	    }
    }

    if (!skipin) {
	    for (j=iset[J0]; j<iset[J1]; j++)
	    {
		    for (k=iset[K0]; k<iset[K1]; k++)
		    {
			    if (!(iset[SKIPIN] && j>=iset[SKIPIN_JMIN] && k>=iset[SKIPIN_KMIN] && j<=iset[SKIPIN_JMAX] && k<=iset[SKIPIN_KMAX]))
			    { 

				    i=iset[I1];
				    r = k_dist(x, y, z, i, j, k);
				    theta = JCPI*size*(x-i)/r;
                                    r *= dx;

				    nn = shift - istart + (int)(step + 1 + r/lsdt);
				    a = JCPI/r*size/2.0;
				    b = -theta/(r*r);
				    c = -theta/(LS*r)*dvadt;
				    if (nn < 3 || (nn - shift)>=(ndat-3)) {
					    continue;
				    }
				    pos = CAT(i, j, k);
				    dev = d_ex[pos];
				    d_ff_ex[nn] += -c*dev;
				    d_ff_ex[nn - 1] += a*((d_ex[pos+zyres] - d_ex[pos-zyres])/dx) + b*dev;
				    d_ff_ex[nn - 2] += c*dev;

				    dev = d_ey[pos];
				    d_ff_ey[nn] += -c*dev;
				    d_ff_ey[nn - 1] += a*((d_ey[pos+zyres] - d_ey[pos-zyres])/dy) + b*dev;
				    d_ff_ey[nn - 2] += c*dev;

				    dev = d_ez[pos];
				    d_ff_ez[nn] += -c*dev; 
				    d_ff_ez[nn - 1] += a*((d_ez[pos+zyres] - d_ez[pos-zyres])/dz) + b*dev;
				    d_ff_ez[nn - 2] += c*dev;
			    }
		    }
	    }
    }

    size = dx*dz;
    if (!skipj0) { 
	    for (i=iset[I0]; i<iset[I1]; i++)
	    {
		    for (k=iset[K0]; k<iset[K1]; k++)
		    {
			    if (!(iset[SKIPJ0] && i>=iset[SKIPJ0_IMIN] && k>=iset[SKIPJ0_KMIN] && i<=iset[SKIPJ0_IMAX] && k<=iset[SKIPJ0_KMAX]))
			    {

				    j=iset[J0];
				    r = k_dist(x, y, z, i, j, k);
				    theta = JCPI*size*(j-y)/r; 
                                    r *= dx;

				    nn = shift - istart + (int)(step + 1 + r/lsdt);
				    a = JCPI/r*size/2.0;
				    b = -theta/(r*r);
				    c = -theta/(LS*r)*dvadt;

				    if (nn < 3 || (nn-shift)>=(ndat-3)) {
					    continue;
				    }
				    pos = CAT(i, j, k);
				    dev = d_ex[pos];
				    d_ff_ex[nn] += -c*dev;
				    d_ff_ex[nn - 1] += -a*((d_ex[pos+zres] - d_ex[pos-zres])/dx) + b*dev;
				    d_ff_ex[nn - 2] += c*dev;

				    dev = d_ey[pos];
				    d_ff_ey[nn] += -c*dev;
				    d_ff_ey[nn - 1] += -a*((d_ey[pos+zres] - d_ey[pos-zres])/dy) + b*dev;
				    d_ff_ey[nn - 2] += c*dev;

				    dev = d_ez[pos];
				    d_ff_ez[nn] += -c*dev; 
				    d_ff_ez[nn - 1] += -a*((d_ez[pos+zres] - d_ez[pos-zres])/dz) + b*dev;
				    d_ff_ez[nn - 2] += c*dev;

			    }
		    }
	    }
    }
    if (!skipjn) {
	    for (i=iset[I0]; i<iset[I1]; i++)
	    {
		    for (k=iset[K0]; k<iset[K1]; k++)
		    {
			    if (!(iset[SKIPJN] && i>=iset[SKIPJN_IMIN] && k>=iset[SKIPJN_KMIN] && i<=iset[SKIPJN_IMAX] && k<=iset[SKIPJN_KMAX]))
			    {

				    j=iset[J1];
				    r = k_dist(x, y, z, i, j, k);
				    theta = JCPI*size*(y-j)/r;
                                    r *= dx;

				    nn = shift - istart + (int)(step + 1 + r/lsdt);
				    a = JCPI/r*size/2.0;
				    b = -theta/(r*r);
				    c = -theta/(LS*r)*dvadt;

				    if (nn < 3 || (nn-shift)>=(ndat-3)) {
					    continue;
				    }
				    pos = CAT(i, j, k);
				    dev = d_ex[pos];
				    d_ff_ex[nn] += -c*dev;
				    d_ff_ex[nn - 1] += a*((d_ex[pos+zres] - d_ex[pos-zres])/dx) + b*dev;
				    d_ff_ex[nn - 2] += c*dev;

				    dev = d_ey[pos];
				    d_ff_ey[nn] += -c*dev;
				    d_ff_ey[nn - 1] += a*((d_ey[pos+zres] - d_ey[pos-zres])/dy) + b*dev;
				    d_ff_ey[nn - 2] += c*dev;

				    dev = d_ez[pos];
				    d_ff_ez[nn] += -c*dev; 
				    d_ff_ez[nn - 1] += a*((d_ez[pos+zres] - d_ez[pos-zres])/dz) + b*dev;
				    d_ff_ez[nn - 2] += c*dev;

			    }

		    }     

	    }
    }
    size = dx*dy;
    if (!skipk0) {
            for (i=iset[I0]; i<iset[I1]; i++)
	    {
		    for (j=iset[J0]; j<iset[J1]; j++)
		    {
			    if (!(iset[SKIPK0] && i>=iset[SKIPK0_IMIN] && j>=iset[SKIPK0_JMIN] && i<=iset[SKIPK0_IMAX] && j<=iset[SKIPK0_JMAX]))
			    {
				    k = iset[K0];
				    r = k_dist(x, y, z, i, j, k);
				    theta = JCPI*size*(k-z)/r; 
                                    r *= dx;

				    nn = shift - istart +(int)(step + 1 + r/lsdt);
				    a = JCPI/r*size/2.0;
				    b = -theta/(r*r);
				    c = -theta/(LS*r)*dvadt;

				    if (nn < 3 || (nn-shift)>=(ndat-3)) {
					    continue;
				    }
				    pos = CAT(i, j, k);
				    dev = d_ex[pos];
				    d_ff_ex[nn] += -c*dev;
				    d_ff_ex[nn - 1] += -a*((d_ex[pos+1] - d_ex[pos-1])/dx) + b*dev;
				    d_ff_ex[nn - 2] += c*dev;

				    dev = d_ey[pos];
				    d_ff_ey[nn] += -c*dev;
				    d_ff_ey[nn - 1] += -a*((d_ey[pos+1] - d_ey[pos-1])/dy) + b*dev;
				    d_ff_ey[nn - 2] += c*dev;

				    dev = d_ez[pos];
				    d_ff_ez[nn] += -c*dev; 
				    d_ff_ez[nn - 1] += -a*((d_ez[pos+1] - d_ez[pos-1])/dz) + b*dev;
				    d_ff_ez[nn - 2] += c*dev;

			    }
		    }
	    }
    }
    if (!skipkn) {
	    for (i=iset[I0]; i<iset[I1]; i++)
	    {
		    for (j=iset[J0]; j<iset[J1]; j++)
		    {
			    if (!(iset[SKIPKN] && i>=iset[SKIPKN_IMIN] && j>=iset[SKIPKN_JMIN] && i<=iset[SKIPKN_IMAX] && j<=iset[SKIPKN_JMAX]))
			    {
				    k = iset[K1];
				    r = k_dist(x, y, z, i, j, k);
				    theta = JCPI*size*(z-k)/r;
                                    r *= dx;

				    nn = shift - istart +(int)(step + 1 + r/lsdt);
				    a = JCPI/r*size/2.0;
				    b = -theta/(r*r);
				    c = -theta/(LS*r)*dvadt;
				    if (nn < 3 || (nn-shift)>=(ndat-3)) {
					    continue;
				    }
				    pos = CAT(i, j, k);

				    dev = d_ex[pos];
				    d_ff_ex[nn] += -c*dev;
				    d_ff_ex[nn - 1] += a*((d_ex[pos+1] - d_ex[pos-1])/dx) + b*dev;
				    d_ff_ex[nn - 2] += c*dev;

				    dev = d_ey[pos];
				    d_ff_ey[nn] += -c*dev;
				    d_ff_ey[nn - 1] += a*((d_ey[pos+1] - d_ey[pos-1])/dy) + b*dev;
				    d_ff_ey[nn - 2] += c*dev;

				    dev = d_ez[pos];
				    d_ff_ez[nn] += -c*dev; 
				    d_ff_ez[nn - 1] += a*((d_ez[pos+1] - d_ez[pos-1])/dz) + b*dev;
				    d_ff_ez[nn - 2] += c*dev;

			    }
		    }
	    }
    }
}

__global__  void
ffKernel_a(float *d_ex, float *d_ey, float *d_ez, 
         int* iset, float *d_ff_ex, float *d_ff_ey, float *d_ff_ez,
         int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int step, int nhlps)
{
    int i, j, k,  nn, npnt, ndat, shift, istart, thiscase;
    long int pos;
    float x, y, z, lsdt;
    float r, theta, size, a, b, c, dev;
    float dvadt = 0.5/dt;  
    bool skipi0, skipin, skipj0, skipjn, skipk0, skipkn;
    int zyres = zres*yres; 

    skipi0 = skipin = skipj0 = skipjn = skipk0 = skipkn = 0;
    if (iset[SKIPI0] && iset[SKIPI0_JMIN]<=iset[J0] && iset[SKIPI0_KMIN]<=iset[K0] && iset[SKIPI0_JMAX]>=iset[J1] && iset[SKIPI0_KMAX]>=iset[K1]) skipi0 = 1;
    if (iset[SKIPIN] && iset[SKIPIN_JMIN]<=iset[J0] && iset[SKIPIN_KMIN]<=iset[K0] && iset[SKIPIN_JMAX]>=iset[J1] && iset[SKIPIN_KMAX]>=iset[K1]) skipin = 1;
    if (iset[SKIPJ0] && iset[SKIPJ0_IMIN]<=iset[I0] && iset[SKIPJ0_KMIN]<=iset[K0] && iset[SKIPJ0_IMAX]>=iset[I1] && iset[SKIPJ0_KMAX]>=iset[K1]) skipj0 = 1;
    if (iset[SKIPJN] && iset[SKIPJN_IMIN]<=iset[I0] && iset[SKIPJN_KMIN]<=iset[K0] && iset[SKIPJN_IMAX]>=iset[I1] && iset[SKIPJN_KMAX]>=iset[K1]) skipjn = 1;
    if (iset[SKIPK0] && iset[SKIPK0_IMIN]<=iset[I0] && iset[SKIPK0_JMIN]<=iset[J0] && iset[SKIPK0_IMAX]>=iset[I1] && iset[SKIPK0_JMAX]>=iset[J1]) skipk0 = 1;
    if (iset[SKIPKN] && iset[SKIPKN_IMIN]<=iset[I0] && iset[SKIPKN_JMIN]<=iset[J0] && iset[SKIPKN_IMAX]>=iset[I1] && iset[SKIPKN_JMAX]>=iset[J1]) skipkn = 1;


    //get far field point to be computed
    npnt = (int)threadIdx.x/(int)nhlps;    //number of computed point from npnt*nhlps threads
    thiscase = threadIdx.x - npnt*nhlps;   //what to compute 0..nhlps
    //find position of this point and other properties
    ndat = iset[NPNTS + 1 + 5*npnt];
    istart = iset[NPNTS + 2 + 5*npnt];
    x = (float)iset[NPNTS + 3 + 5*npnt];
    y = (float)iset[NPNTS + 4 + 5*npnt];
    z = (float)iset[NPNTS + 5 + 5*npnt];

    //determine shift in accumulator for this point
    shift = 0;
    for (i=0; i<npnt; i++) 
       shift += iset[NPNTS + 1 + 5*i]*nhlps;   //helpers are arranged point by point, find shift for this point
    shift += ndat*thiscase;                   //find shift for this helper

    lsdt = LS*dt;
    size = dx*dz;
    
    if (!skipi0) {
	    for (j=(iset[J0] + thiscase*(iset[J1]-iset[J0])/nhlps); j<(iset[J0] + (thiscase+1)*(iset[J1]-iset[J0])/nhlps); j++)
	    {
		    for (k=iset[K0]; k<iset[K1]; k++)
		    {
			    if (!(iset[SKIPI0] && j>=iset[SKIPI0_JMIN] && k>=iset[SKIPI0_KMIN] && j<=iset[SKIPI0_JMAX] && k<=iset[SKIPI0_KMAX]))
			    {
				    i=iset[I0];
				    r = k_dist(x, y, z, i, j, k);
				    theta = JCPI*size*(i-x)/r;
                                    r *= dx;

				    nn = shift - istart + (int)(step + 1 + r/lsdt);
				    a = JCPI/r*size/2.0;
				    b = -theta/(r*r);
				    c = -theta/(LS*r)*dvadt;
				    if (nn < 3 || (nn-shift)>=(ndat-3)) {
					    continue;
				    }
				    pos = CAT(i, j, k);
				    dev = d_ex[pos];

				    d_ff_ex[nn] += -c*dev;
				    d_ff_ex[nn - 1] += -a*((d_ex[pos+zyres] - d_ex[pos-zyres])/dx) + b*dev;
				    d_ff_ex[nn - 2] += c*dev;

				    dev = d_ey[pos];
				    d_ff_ey[nn] += -c*dev;
				    d_ff_ey[nn - 1] += -a*((d_ey[pos+zyres] - d_ey[pos-zyres])/dy) + b*dev;
				    d_ff_ey[nn - 2] += c*dev;

				    dev = d_ez[pos];
				    d_ff_ez[nn] += -c*dev; 
				    d_ff_ez[nn - 1] += -a*((d_ez[pos+zyres] - d_ez[pos-zyres])/dz) + b*dev;
				    d_ff_ez[nn - 2] += c*dev;
			    }
		    }
	    }
    }

    if (!skipin) {
	    for (j=(iset[J0] + thiscase*(iset[J1]-iset[J0])/nhlps); j<(iset[J0] + (thiscase+1)*(iset[J1]-iset[J0])/nhlps); j++)
	    {
		    for (k=iset[K0]; k<iset[K1]; k++)
		    {
			    if (!(iset[SKIPIN] && j>=iset[SKIPIN_JMIN] && k>=iset[SKIPIN_KMIN] && j<=iset[SKIPIN_JMAX] && k<=iset[SKIPIN_KMAX]))
			    { 

				    i=iset[I1];
				    r = k_dist(x, y, z, i, j, k);
				    theta = JCPI*size*(x-i)/r;
                                    r *= dx;

				    nn = shift - istart + (int)(step + 1 + r/lsdt);
				    a = JCPI/r*size/2.0;
				    b = -theta/(r*r);
				    c = -theta/(LS*r)*dvadt;
				    if (nn < 3 || (nn - shift)>=(ndat-3)) {
					    continue;
				    }
				    pos = CAT(i, j, k);
				    dev = d_ex[pos];
				    d_ff_ex[nn] += -c*dev;
				    d_ff_ex[nn - 1] += a*((d_ex[pos+zyres] - d_ex[pos-zyres])/dx) + b*dev;
				    d_ff_ex[nn - 2] += c*dev;

				    dev = d_ey[pos];
				    d_ff_ey[nn] += -c*dev;
				    d_ff_ey[nn - 1] += a*((d_ey[pos+zyres] - d_ey[pos-zyres])/dy) + b*dev;
				    d_ff_ey[nn - 2] += c*dev;

				    dev = d_ez[pos];
				    d_ff_ez[nn] += -c*dev; 
				    d_ff_ez[nn - 1] += a*((d_ez[pos+zyres] - d_ez[pos-zyres])/dz) + b*dev;
				    d_ff_ez[nn - 2] += c*dev;
			    }
		    }
	    }
    }

    size = dx*dz;
    if (!skipj0) { 
	    for (i=(iset[I0] + thiscase*(iset[I1]-iset[I0])/nhlps); i<(iset[I0] + (thiscase+1)*(iset[I1]-iset[I0])/nhlps); i++)
	    {
		    for (k=iset[K0]; k<iset[K1]; k++)
		    {
			    if (!(iset[SKIPJ0] && i>=iset[SKIPJ0_IMIN] && k>=iset[SKIPJ0_KMIN] && i<=iset[SKIPJ0_IMAX] && k<=iset[SKIPJ0_KMAX]))
			    {

				    j=iset[J0];
				    r = k_dist(x, y, z, i, j, k);
				    theta = JCPI*size*(j-y)/r; 
                                    r *= dx;

				    nn = shift - istart + (int)(step + 1 + r/lsdt);
				    a = JCPI/r*size/2.0;
				    b = -theta/(r*r);
				    c = -theta/(LS*r)*dvadt;

				    if (nn < 3 || (nn-shift)>=(ndat-3)) {
					    continue;
				    }
				    pos = CAT(i, j, k);
				    dev = d_ex[pos];
				    d_ff_ex[nn] += -c*dev;
				    d_ff_ex[nn - 1] += -a*((d_ex[pos+zres] - d_ex[pos-zres])/dx) + b*dev;
				    d_ff_ex[nn - 2] += c*dev;

				    dev = d_ey[pos];
				    d_ff_ey[nn] += -c*dev;
				    d_ff_ey[nn - 1] += -a*((d_ey[pos+zres] - d_ey[pos-zres])/dy) + b*dev;
				    d_ff_ey[nn - 2] += c*dev;

				    dev = d_ez[pos];
				    d_ff_ez[nn] += -c*dev; 
				    d_ff_ez[nn - 1] += -a*((d_ez[pos+zres] - d_ez[pos-zres])/dz) + b*dev;
				    d_ff_ez[nn - 2] += c*dev;

			    }
		    }
	    }
    }
    if (!skipjn) {
	    for (i=(iset[I0] + thiscase*(iset[I1]-iset[I0])/nhlps); i<(iset[I0] + (thiscase+1)*(iset[I1]-iset[I0])/nhlps); i++)
	    {
		    for (k=iset[K0]; k<iset[K1]; k++)
		    {
			    if (!(iset[SKIPJN] && i>=iset[SKIPJN_IMIN] && k>=iset[SKIPJN_KMIN] && i<=iset[SKIPJN_IMAX] && k<=iset[SKIPJN_KMAX]))
			    {

				    j=iset[J1];
				    r = k_dist(x, y, z, i, j, k);
				    theta = JCPI*size*(y-j)/r;
                                    r *= dx;

				    nn = shift - istart + (int)(step + 1 + r/lsdt);
				    a = JCPI/r*size/2.0;
				    b = -theta/(r*r);
				    c = -theta/(LS*r)*dvadt;

				    if (nn < 3 || (nn-shift)>=(ndat-3)) {
					    continue;
				    }
				    pos = CAT(i, j, k);
				    dev = d_ex[pos];
				    d_ff_ex[nn] += -c*dev;
				    d_ff_ex[nn - 1] += a*((d_ex[pos+zres] - d_ex[pos-zres])/dx) + b*dev;
				    d_ff_ex[nn - 2] += c*dev;

				    dev = d_ey[pos];
				    d_ff_ey[nn] += -c*dev;
				    d_ff_ey[nn - 1] += a*((d_ey[pos+zres] - d_ey[pos-zres])/dy) + b*dev;
				    d_ff_ey[nn - 2] += c*dev;

				    dev = d_ez[pos];
				    d_ff_ez[nn] += -c*dev; 
				    d_ff_ez[nn - 1] += a*((d_ez[pos+zres] - d_ez[pos-zres])/dz) + b*dev;
				    d_ff_ez[nn - 2] += c*dev;

			    }

		    }     

	    }
    }
    if (!skipk0) {
	    size = dx*dy;
	    for (i=(iset[I0] + thiscase*(iset[I1]-iset[I0])/nhlps); i<(iset[I0] + (thiscase+1)*(iset[I1]-iset[I0])/nhlps); i++)
	    {
		    for (j=iset[J0]; j<iset[J1]; j++)
		    {
			    if (!(iset[SKIPK0] && i>=iset[SKIPK0_IMIN] && j>=iset[SKIPK0_JMIN] && i<=iset[SKIPK0_IMAX] && j<=iset[SKIPK0_JMAX]))
			    {
				    k = iset[K0];
				    r = k_dist(x, y, z, i, j, k);
				    theta = JCPI*size*(k-z)/r; 
                                    r *= dx;

				    nn = shift - istart +(int)(step + 1 + r/lsdt);
				    a = JCPI/r*size/2.0;
				    b = -theta/(r*r);
				    c = -theta/(LS*r)*dvadt;

				    if (nn < 3 || (nn-shift)>=(ndat-3)) {
					    continue;
				    }
				    pos = CAT(i, j, k);
				    dev = d_ex[pos];
				    d_ff_ex[nn] += -c*dev;
				    d_ff_ex[nn - 1] += -a*((d_ex[pos+1] - d_ex[pos-1])/dx) + b*dev;
				    d_ff_ex[nn - 2] += c*dev;

				    dev = d_ey[pos];
				    d_ff_ey[nn] += -c*dev;
				    d_ff_ey[nn - 1] += -a*((d_ey[pos+1] - d_ey[pos-1])/dy) + b*dev;
				    d_ff_ey[nn - 2] += c*dev;

				    dev = d_ez[pos];
				    d_ff_ez[nn] += -c*dev; 
				    d_ff_ez[nn - 1] += -a*((d_ez[pos+1] - d_ez[pos-1])/dz) + b*dev;
				    d_ff_ez[nn - 2] += c*dev;

			    }
		    }
	    }
    }
    if (!skipkn) {
	    for (i=(iset[I0] + thiscase*(iset[I1]-iset[I0])/nhlps); i<(iset[I0] + (thiscase+1)*(iset[I1]-iset[I0])/nhlps); i++)
	    {
		    for (j=iset[J0]; j<iset[J1]; j++)
		    {
			    if (!(iset[SKIPKN] && i>=iset[SKIPKN_IMIN] && j>=iset[SKIPKN_JMIN] && i<=iset[SKIPKN_IMAX] && j<=iset[SKIPKN_JMAX]))
			    {
				    k = iset[K1];
				    r = k_dist(x, y, z, i, j, k);
				    theta = JCPI*size*(z-k)/r;
                                    r *= dx;

				    nn = shift - istart +(int)(step + 1 + r/lsdt);
				    a = JCPI/r*size/2.0;
				    b = -theta/(r*r);
				    c = -theta/(LS*r)*dvadt;
				    if (nn < 3 || (nn-shift)>=(ndat-3)) {
					    continue;
				    }
				    pos = CAT(i, j, k);

				    dev = d_ex[pos];
				    d_ff_ex[nn] += -c*dev;
				    d_ff_ex[nn - 1] += a*((d_ex[pos+1] - d_ex[pos-1])/dx) + b*dev;
				    d_ff_ex[nn - 2] += c*dev;

				    dev = d_ey[pos];
				    d_ff_ey[nn] += -c*dev;
				    d_ff_ey[nn - 1] += a*((d_ey[pos+1] - d_ey[pos-1])/dy) + b*dev;
				    d_ff_ey[nn - 2] += c*dev;

				    dev = d_ez[pos];
				    d_ff_ez[nn] += -c*dev; 
				    d_ff_ez[nn - 1] += a*((d_ez[pos+1] - d_ez[pos-1])/dz) + b*dev;
				    d_ff_ez[nn - 2] += c*dev;

				    //ipr = i;
				    //jpr = j;
				    //kpr = k;

			    }
		    }
	    }
    }

}

__global__  void
ffKernel_b(float *d_ff_ex, float *d_ff_ey, float *d_ff_ez,
         float *d_ff_hlp_ex, float *d_ff_hlp_ey, float *d_ff_hlp_ez,
         int *iset, int nhlps)
{
    int i, pnt, j, k, ndat, ashift, hshift, baseshift;

    for (i=0; i<iset[NPNTS]; i++) //iterate over points
    {
        ndat = iset[NPNTS + 1 + 5*i];
        baseshift = ashift = 0;

        for (pnt=0; pnt<i; pnt++) {
		baseshift += iset[NPNTS + 1 + 5*pnt]*nhlps;
                ashift += iset[NPNTS + 1 + 5*pnt];
        }

        for (j=0; j<ndat; j++) {
            d_ff_ex[ashift + j] = 0;
            d_ff_ey[ashift + j] = 0;
            d_ff_ez[ashift + j] = 0;
        }
        
        for (k=0; k<nhlps; k++) //iterate through helpers
        {
             hshift = baseshift + ndat*k;

             for (j=0; j<ndat; j++) //iterate through accumulator
             {
                  d_ff_ex[ashift + j] += d_ff_hlp_ex[hshift + j];
                  d_ff_ey[ashift + j] += d_ff_hlp_ey[hshift + j];
                  d_ff_ez[ashift + j] += d_ff_hlp_ez[hshift + j];                  
             }
        }
    }

}
 
__global__  void
pffKernel(float *d_ex, float *d_ey, float *d_ez, 
         int* piset, float *d_pff_ex, float *d_pff_ey, float *d_pff_ez,
         int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int step)
{
    int i, j, k,  nn, npnt, ndat, shift, istart;
    long int pos;
    float x, y, z;
    float r, theta, size, a, b, c, dev;
    float dvadt = 0.5/dt;  
    int pimin, pimax, pjmin, pjmax, pi, pj, pri, prj;

    
    //get far field point to be computed
    //npnt = threadIdx.x;
//    npnt = blockIdx.y*PNDIV + blockIdx.x;
    npnt = blockIdx.x*PNDIV + threadIdx.x; //(thread jde do PNDIV)
    if (npnt>=piset[10]) return; 

    //find position of this point and other properties
    ndat = piset[PNPNTS + 1 + 5*npnt];
    istart = piset[PNPNTS + 2 + 5*npnt];
    x = (float)piset[PNPNTS + 3 + 5*npnt];
    y = (float)piset[PNPNTS + 4 + 5*npnt];
    z = (float)piset[PNPNTS + 5 + 5*npnt];
    pimin = piset[PIMIN];
    pimax = piset[PIMAX];
    pjmin = piset[PJMIN];
    pjmax = piset[PJMAX];

    //determine shift in accumulator for this point
    shift = 0;
    for (i=0; i<npnt; i++) 
	    shift += piset[PNPNTS + 1 + 5*i];

    size = dx*dy;

    

    for (pi = pimin; pi<pimax; pi++)
    {
	    for (pj = pjmin; pj<pjmax; pj++)
	    {
		    pri = (piset[I1]-piset[I0])*pi;
		    prj = (piset[J1]-piset[J0])*pj;

		    for (i=piset[I0]; i<piset[I1]; i++)
		    {
			    for (j=piset[J0]; j<piset[J1]; j++)
			    {
				    if (!(piset[PSKIPK0] && i>=piset[PSKIPK0_IMIN] && j>=piset[PSKIPK0_JMIN] && i<=piset[PSKIPK0_IMAX] && j<=piset[PSKIPK0_JMAX]))
				    {
					    k = piset[K0];
					    r = sqrt((x-(i+pri))*(x-(i+pri)) + (y-(j+prj))*(y-(j+prj)) + (z-k)*(z-k));
					    theta = (k-z)/r;
					    r *= dx;

					    nn = shift - istart +(int)(step + 1 + r/(LS*dt));
					    a = JCPI/r*size;
					    b = -JCPI*theta/(r*r)*size;
					    c = -JCPI*theta/(LS*r)*size*dvadt;

					    if (nn < 3 || (nn-shift)>=(ndat-3)) {
						    continue;
					    }
					    pos = CAT(i, j, k);
					    dev = d_ex[pos];
					    d_pff_ex[nn] += -c*dev;
					    d_pff_ex[nn - 1] += -a*((d_ex[CAT(i,j,k+1)] - d_ex[CAT(i,j,k-1)])/2/dx) + b*dev;
					    d_pff_ex[nn - 2] += c*dev;

					    dev = d_ey[pos];
					    d_pff_ey[nn] += -c*dev;
					    d_pff_ey[nn - 1] += -a*((d_ey[CAT(i,j,k+1)] - d_ey[CAT(i,j,k-1)])/2/dy) + b*dev;
					    d_pff_ey[nn - 2] += c*dev;

					    dev = d_ez[pos];
					    d_pff_ez[nn] += -c*dev; 
					    d_pff_ez[nn - 1] += -a*((d_ez[CAT(i,j,k+1)] - d_ez[CAT(i,j,k-1)])/2/dz) + b*dev;
					    d_pff_ez[nn - 2] += c*dev;

				    }
				    if (!(piset[PSKIPKN] && i>=piset[PSKIPKN_IMIN] && j>=piset[PSKIPKN_JMIN] && i<=piset[PSKIPKN_IMAX] && j<=piset[PSKIPKN_JMAX]))
				    {
					    k = piset[K1];
					    r = sqrt((x-(i+pri))*(x-(i+pri)) + (y-(j+prj))*(y-(j+prj)) + (z-k)*(z-k));
					    theta = -(k-z)/r;
					    r *= dx;

					    nn = shift - istart +(int)(step + 1 + r/(LS*dt));
					    a = JCPI/r*size;
					    b = -JCPI*theta/(r*r)*size;
					    c = -JCPI*theta/(LS*r)*size*dvadt;
					    
                                            if (nn < 3 || (nn-shift)>=(ndat-3)) {
						    continue;
					    }
					    pos = CAT(i, j, k);
					    dev = d_ex[pos];
					    d_pff_ex[nn] += -c*dev;
					    d_pff_ex[nn - 1] += a*((d_ex[CAT(i,j,k+1)] - d_ex[CAT(i,j,k-1)])/2/dx) + b*dev;
					    d_pff_ex[nn - 2] += c*dev;

					    dev = d_ey[pos];
					    d_pff_ey[nn] += -c*dev;
					    d_pff_ey[nn - 1] += a*((d_ey[CAT(i,j,k+1)] - d_ey[CAT(i,j,k-1)])/2/dy) + b*dev;
					    d_pff_ey[nn - 2] += c*dev;

					    dev = d_ez[pos];
					    d_pff_ez[nn] += -c*dev; 
					    d_pff_ez[nn - 1] += a*((d_ez[CAT(i,j,k+1)] - d_ez[CAT(i,j,k-1)])/2/dz) + b*dev;
					    d_pff_ez[nn - 2] += c*dev;

				    }     
			    }
		    }

	    }

    }
//    d_pff_ex[shift] = npnt;
//    d_pff_ey[shift] = 12345;
//    d_pff_ez[shift] = 65432;

}

/*
__global__  void
pffsKernel(float *d_ex, float *d_ey, float *d_ez, 
		int* piset, float *d_pff_ex, float *d_pff_ey, float *d_pff_ez,
		int xres, int yres, int zres,
		float dx, float dy, float dz, float dt, int step)
{
	int i, j, k, m, nn, npnt, ndat, shift, istart, nth;
	long int pos;
	float x, y, z;
	float r, theta, size, a, b, c, dev;
	float dvadt = 0.5/dt;  
	int pimin = piset[PIMIN]; 
        int pimax = piset[PIMAX];
        int pjmin = piset[PJMIN];
        int pjmax = piset[PJMAX];
        int pi, pj, pri, prj;

        extern __shared__ float shared[];

        float* sex = (float*)shared; 
        float* sey = (float*)&sex[(pimax-pimin)*(pjmax-pjmin)*2000]; 
        float* sez = (float*)&sey[(pimax-pimin)*(pjmax-pjmin)*2000];
        int* sepos = (int*)&sez[(pimax-pimin)*(pjmax-pjmin)*2000]; 
        int* pis = (int*)&sepos[(pimax-pimin)*(pjmax-pjmin)]; 
        int* pjs = (int*)&pis[(pimax-pimin)*(pjmax-pjmin)]; 

        //FIXME, all this fails due to too small shared memory on gpu.

	//get repetion thread and far field point to be computed
	nth = threadIdx.x;
	npnt = blockIdx.y*PNDIV + blockIdx.x;
	if (npnt>=piset[10]) return; 

	//find position of this point and other properties
	ndat = piset[PNPNTS + 1 + 5*npnt];
	istart = piset[PNPNTS + 2 + 5*npnt];
	x = (float)piset[PNPNTS + 3 + 5*npnt];
	y = (float)piset[PNPNTS + 4 + 5*npnt];
	z = (float)piset[PNPNTS + 5 + 5*npnt];

	//determine shift in accumulator for this point
	shift = 0;
	for (i=0; i<npnt; i++) 
		shift += piset[PNPNTS + 1 + 5*i];

	m = 0;

	//we have (pimax-pimin)*(pjmax-pjmin) threads available
	pi = pis[nth];
	pj = pjs[nth];

	//we have (pimax-pimin)*(pjmax-pjmin) threads available
	pri = (piset[I1]-piset[I0])*pi;
	prj = (piset[J1]-piset[J0])*pj;

	r = sqrt((x-(i+pri))*(x-(i+pri)) + (y-(j+prj))*(y-(j+prj)) + (z-k)*(z-k));
	nn = shift - istart +(int)(step + 1 + r/(LS*dt));

	pis[nth] = pi;
	pjs[nth] = pj; 
	sepos[m++] = nn - 1000;
	for (i=0; i<2000; i++) {
		sex[nth*2000 + i] = 0;
		sey[nth*2000 + i] = 0;
		sez[nth*2000 + i] = 0;
	}

	size = dx*dy;

	pri = (piset[I1]-piset[I0])*pi;
	prj = (piset[J1]-piset[J0])*pj;

//FIXME: sex is shared, but 200*(pimax-pimin)*(pjmax-pjmin)

        __syncthreads();
	for (i=piset[I0]; i<piset[I1]; i++)
	{
		for (j=piset[J0]; j<piset[J1]; j++)
		{
			if (!(piset[PSKIPK0] && i>=piset[PSKIPK0_IMIN] && j>=piset[PSKIPK0_JMIN] && i<=piset[PSKIPK0_IMAX] && j<=piset[PSKIPK0_JMAX]))
			{
				k = piset[K0];
				r = sqrt((x-(i+pri))*(x-(i+pri)) + (y-(j+prj))*(y-(j+prj)) + (z-k)*(z-k));
				theta = (k-z)/r;
				r *= dx;

				nn = shift - istart +(int)(step + 1 + r/(LS*dt)) - sepos[m];
				a = JCPI/r*size;
				b = -JCPI*theta/(r*r)*size;
				c = -JCPI*theta/(LS*r)*size*dvadt;

				if (nn < 3 || (nn-shift)>=(ndat-3)) {
					continue;
				}

                                nn += nth*2000; //find appropriate accumulator

				pos = CAT(i, j, k);
				dev = d_ex[pos];
				sex[nn] += -c*dev;
				sex[nn - 1] += -a*((d_ex[CAT(i,j,k+1)] - d_ex[CAT(i,j,k-1)])/2/dx) + b*dev;
				sex[nn - 2] += c*dev;

				dev = d_ey[pos];
				sey[nn] += -c*dev;
				sey[nn - 1] += -a*((d_ey[CAT(i,j,k+1)] - d_ey[CAT(i,j,k-1)])/2/dy) + b*dev;
				sey[nn - 2] += c*dev;

				dev = d_ez[pos];
				sez[nn] += -c*dev; 
				sez[nn - 1] += -a*((d_ez[CAT(i,j,k+1)] - d_ez[CAT(i,j,k-1)])/2/dz) + b*dev;
				sez[nn - 2] += c*dev;
			}
		}
	}
        __syncthreads();

//FIXME: sex is shared, but 200*(pimax-pimin)*(pjmax-pjmin)
	if (nth == 0) {
                for (m=0; m<((pimax-pimin)*(pjmax-pjmin)); m++) 
                {
                        for (i=0; i<2000; i++) {
				d_pff_ex[sepos[m]+1000 + i] += sex[m*2000+i];  //this must be global first!!!!! or better, shared but dynamically allocated for all the threads!
                                d_pff_ey[sepos[m]+1000 + i] += sey[m*2000+i];
                                d_pff_ez[sepos[m]+1000 + i] += sez[m*2000+i];
			}
		}
	}
}
*/
__global__  void
pffKernel_a(float *d_ex, float *d_ey, float *d_ez, 
         int* piset, float *d_pff_ex, float *d_pff_ey, float *d_pff_ez,
         int xres, int yres, int zres,
         float dx, float dy, float dz, float dt, int step)
{
    int i, j, k,  nn, npnt, ndat, shift, istart, thiscase;
    long int pos;
    float x, y, z;
    float r, theta, size, a, b, c, dev;
    float dvadt = 0.5/dt;  
    int pi, pj, pri, prj;
    int nhlps = (piset[PIMAX]-piset[PIMIN])*(piset[PJMAX]-piset[PJMIN]);

    npnt = (int)threadIdx.x/(int)nhlps;
    thiscase = threadIdx.x - npnt*nhlps;   
    
    //get far field point to be computed
    //npnt = threadIdx.x;
    npnt = blockIdx.y*PNDIV + blockIdx.x;
    if (npnt>=piset[10]) return; 

    //find position of this point and other properties
    ndat = piset[PNPNTS + 1 + 5*npnt];
    istart = piset[PNPNTS + 2 + 5*npnt];
    x = (float)piset[PNPNTS + 3 + 5*npnt];
    y = (float)piset[PNPNTS + 4 + 5*npnt];
    z = (float)piset[PNPNTS + 5 + 5*npnt];
    
    //determine shift in accumulator for this point
    shift = 0;
    for (i=0; i<npnt; i++) 
	    shift += piset[PNPNTS + 1 + 5*i]*nhlps;
    shift += ndat*thiscase;

    pi = (int)(threadIdx.x/(piset[PIMAX]-piset[PIMIN]));
    pj = threadIdx.x - pi*(piset[PIMAX]-piset[PIMIN]);

    pi += piset[PIMIN];
    pj += piset[PJMIN];

    size = dx*dy;
    pri = (piset[I1]-piset[I0])*pi;
    prj = (piset[J1]-piset[J0])*pj;

    for (i=piset[I0]; i<piset[I1]; i++)
    {
	    for (j=piset[J0]; j<piset[J1]; j++)
	    {
		    if (!(piset[PSKIPK0] && i>=piset[PSKIPK0_IMIN] && j>=piset[PSKIPK0_JMIN] && i<=piset[PSKIPK0_IMAX] && j<=piset[PSKIPK0_JMAX])) 
		    {
			    k = piset[K0];
			    r = sqrt((x-(i+pri))*(x-(i+pri)) + (y-(j+prj))*(y-(j+prj)) + (z-k)*(z-k));
			    theta = (k-z)/r;
			    r *= dx;

			    nn = shift - istart +(int)(step + 1 + r/(LS*dt));
			    a = JCPI/r*size;
			    b = -JCPI*theta/(r*r)*size;
			    c = -JCPI*theta/(LS*r)*size*dvadt;

			    if (nn < 3 || (nn-shift)>=(ndat-3)) {
				    continue;
			    }
			    pos = CAT(i, j, k);
			    dev = d_ex[pos];
			    d_pff_ex[nn] += -c*dev;
			    d_pff_ex[nn - 1] += -a*((d_ex[CAT(i,j,k+1)] - d_ex[CAT(i,j,k-1)])/2/dx) + b*dev;
			    d_pff_ex[nn - 2] += c*dev;

			    dev = d_ey[pos];
			    d_pff_ey[nn] += -c*dev;
			    d_pff_ey[nn - 1] += -a*((d_ey[CAT(i,j,k+1)] - d_ey[CAT(i,j,k-1)])/2/dy) + b*dev;
			    d_pff_ey[nn - 2] += c*dev;

			    dev = d_ez[pos];
			    d_pff_ez[nn] += -c*dev; 
			    d_pff_ez[nn - 1] += -a*((d_ez[CAT(i,j,k+1)] - d_ez[CAT(i,j,k-1)])/2/dz) + b*dev;
			    d_pff_ez[nn - 2] += c*dev;
		    }
		    if (!(piset[PSKIPKN] && i>=piset[PSKIPKN_IMIN] && j>=piset[PSKIPKN_JMIN] && i<=piset[PSKIPKN_IMAX] && j<=piset[PSKIPKN_JMAX]))
		    {
			    k = piset[K1];
			    r = sqrt((x-(i+pri))*(x-(i+pri)) + (y-(j+prj))*(y-(j+prj)) + (z-k)*(z-k));
			    theta = -(k-z)/r;
			    r *= dx;

			    nn = shift - istart +(int)(step + 1 + r/(LS*dt));
			    a = JCPI/r*size;
			    b = -JCPI*theta/(r*r)*size;
			    c = -JCPI*theta/(LS*r)*size*dvadt;

			    if (nn < 3 || (nn-shift)>=(ndat-3)) {
				    continue;
			    }
			    pos = CAT(i, j, k);
			    dev = d_ex[pos];
			    d_pff_ex[nn] += -c*dev;
			    d_pff_ex[nn - 1] += a*((d_ex[CAT(i,j,k+1)] - d_ex[CAT(i,j,k-1)])/2/dx) + b*dev;
			    d_pff_ex[nn - 2] += c*dev;

			    dev = d_ey[pos];
			    d_pff_ey[nn] += -c*dev;
			    d_pff_ey[nn - 1] += a*((d_ey[CAT(i,j,k+1)] - d_ey[CAT(i,j,k-1)])/2/dy) + b*dev;
			    d_pff_ey[nn - 2] += c*dev;

			    dev = d_ez[pos];
			    d_pff_ez[nn] += -c*dev; 
			    d_pff_ez[nn - 1] += a*((d_ez[CAT(i,j,k+1)] - d_ez[CAT(i,j,k-1)])/2/dz) + b*dev;
			    d_pff_ez[nn - 2] += c*dev;
		    }    
	    }
    }
}

__global__  void
pffKernel_b(float *d_pff_ex, float *d_pff_ey, float *d_pff_ez,
         float *d_pff_hlp_ex, float *d_pff_hlp_ey, float *d_pff_hlp_ez,
         int *piset)
{
    int i, pnt, j, k, ndat, ashift, hshift, baseshift;
    int nhlps = (piset[PIMAX]-piset[PIMIN])*(piset[PJMAX]-piset[PJMIN]);

    for (i=0; i<piset[10]; i++) //iterate over points
    {
        ndat = piset[PNPNTS + 1 + 5*i];
        baseshift = ashift = 0;

        for (pnt=0; pnt<i; pnt++) {
		baseshift += piset[PNPNTS + 1 + 5*pnt]*nhlps;
                ashift += piset[PNPNTS + 1 + 5*pnt];
        }

        for (j=0; j<ndat; j++) {
            d_pff_ex[ashift + j] = 0;
            d_pff_ey[ashift + j] = 0;
            d_pff_ez[ashift + j] = 0;
        }
        
        for (k=0; k<nhlps; k++) //iterate through helpers
        {
             hshift = baseshift + ndat*k;

             for (j=0; j<ndat; j++) //iterate through accumulator
             {
                  d_pff_ex[ashift + j] += d_pff_hlp_ex[hshift + j];
                  d_pff_ey[ashift + j] += d_pff_hlp_ey[hshift + j];
                  d_pff_ez[ashift + j] += d_pff_hlp_ez[hshift + j];                  
             }
        }
    }

}

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

 
__global__  void
tsffjstepKernel(float *d_jpool_e, float *d_jpool_h, float *d_jpvals,
            float *d_jpool_epsilon, float *d_jpool_mu,
            float *d_jpool_sigma, float *d_jpool_sigast, 
            float dx, float dt, float e, int n,
            int ntsf, int pos)
{
        int i, shift = n*pos, jpshift = 8*pos;
         
        //copybound
        float de0, de1, dh0, dh1, den, dem, dhn, dhm;
        float dAngleMult=1;
        float dMulth=dt/(MU_0*dx)/dAngleMult;
        float dMulte=dt/(EPSILON_0*dx)/dAngleMult;         

        de0 = d_jpvals[0 + jpshift];
        de1 = d_jpvals[1 + jpshift];
        dh0 = d_jpvals[2 + jpshift];
        dh1 = d_jpvals[3 + jpshift];
        den = d_jpvals[4 + jpshift];
        dem = d_jpvals[5 + jpshift];
        dhn = d_jpvals[6 + jpshift];
        dhm = d_jpvals[7 + jpshift];

        //ystep_e
        for (i=1; i<(n); i++)
           d_jpool_e[i + shift] += 1.0/(d_jpool_epsilon[i + shift])*dMulte*(d_jpool_h[i-1 + shift] - d_jpool_h[i + shift]);

        //bound
        d_jpool_e[0 + shift] = de1 + (dt*LS-dx)/(dt*LS+dx)*(d_jpool_e[1 + shift] - de0);
        d_jpool_h[0 + shift] = dh1 + (dt*LS-dx)/(dt*LS+dx)*(d_jpool_h[1 + shift] - dh0);
        d_jpool_e[n-1 + shift] = dem + (dt*LS-dx)/(dt*LS+dx)*(d_jpool_e[n-2 + shift] - den);
        d_jpool_h[n-1 + shift] = dhm + (dt*LS-dx)/(dt*LS+dx)*(d_jpool_h[n-2 + shift] - dhn);


        //copybound
        d_jpvals[0 + jpshift] = d_jpool_e[0 + shift];
        d_jpvals[1 + jpshift] = d_jpool_e[1 + shift];
        d_jpvals[2 + jpshift] = d_jpool_h[0 + shift];
        d_jpvals[3 + jpshift] = d_jpool_h[1 + shift];
        d_jpvals[4 + jpshift] = d_jpool_e[n-1 + shift];
        d_jpvals[5 + jpshift] = d_jpool_e[n-2 + shift];
        d_jpvals[6 + jpshift] = d_jpool_h[n-1 + shift];
        d_jpvals[7 + jpshift] = d_jpool_h[n-2 + shift];
      
        //applysource
        d_jpool_e[0 + shift] = e;

        //ystep_h
        for (i=0; i<(n-1); i++)
           d_jpool_h[i + shift] += 1.0/(d_jpool_mu[i + shift])*dMulth*(d_jpool_e[i + shift] - d_jpool_e[i+1 + shift]);  

}

__global__  void
sfjstepKernel(float *d_jpool_e, float *d_jpool_h, float *d_jpvals,
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
__global__  void
mbnzeKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
          int xres, int yres, int zres, int mb_bz0, int mb_bzn, 
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

    if ((mb_bz0 == 4 || mb_bzn == 4) && (k==mb_bz0pos && i>=mb_bx0pos && i<mb_bxnpos && j>=mb_by0pos && j<mb_bynpos)) {
        pos1 = CAT(i, j, k);
        pos2 = CAT(i, j, mb_bznpos);
        d_ex[pos2] = d_ex[pos1];
        d_ey[pos2] = d_ey[pos1];
     }
}

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
__global__  void
mbnzhKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
          int xres, int yres, int zres, int mb_bz0, int mb_bzn, 
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

    if ((mb_bz0 == 4 || mb_bzn == 4) && (k==mb_bz0pos && i>=mb_bx0pos && i<mb_bxnpos && j>=mb_by0pos && j<mb_bynpos)) {
        pos1 = CAT(i, j, k-1);
        pos2 = CAT(i, j, mb_bznpos-1);
        d_hx[pos1] = d_hx[pos2];
        d_hy[pos1] = d_hy[pos2];
     }

}



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

__device__ float
gaussmult(int x, int y, int z, float xc, float yc, float rx, float ry)
{
    return exp(-((float)(x-xc))*((float)(x-xc))/rx/rx/2.0)*exp(-((float)(y-yc))*((float)(y-yc))/ry/ry/2.0);
}

__device__ float
fibermult(int x, int y, int z, float xc, float yc, float radius, float cutoff, float epsilon_core, float epsilon_cladding, float dx, float lambda)
{
    float v, omegazero, radval;

    v = 2*CU_PI/lambda*radius*dx*(epsilon_core - epsilon_cladding);
    omegazero = (0.65 + 1.619/v/sqrt(v) + 2.879/pow(v,6))*radius*dx;
    radval = ((x-xc)*(x-xc) + (y-yc)*(y-yc));

    if (radval < (cutoff*cutoff)) {
         return exp(-radval*dx*dx/omegazero/omegazero);
    }

    return 0;
}


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
    int fiber = (int)d_tsfset[21];
    float fiber_fxpos = d_tsfset[22];
    float fiber_fypos = d_tsfset[23];
    float fiber_radius = d_tsfset[24];
    float fiber_cutoff = d_tsfset[25];
    float fiber_core = d_tsfset[26];
    float fiber_cladding = d_tsfset[27];
    float fiber_lambda = d_tsfset[28];


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
             else if (fiber) gcorr*=fibermult(i0-1, j, k, fiber_fxpos, fiber_fypos, fiber_radius, fiber_cutoff, fiber_core, fiber_cladding, dx, fiber_lambda);

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
             else if (fiber) gcorr*=fibermult(i1, j, k, fiber_fxpos, fiber_fypos, fiber_radius, fiber_cutoff, fiber_core, fiber_cladding, dx, fiber_lambda);

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
    int fiber = (int)d_tsfset[21];
    float fiber_fxpos = d_tsfset[22];
    float fiber_fypos = d_tsfset[23];
    float fiber_radius = d_tsfset[24];
    float fiber_cutoff = d_tsfset[25];
    float fiber_core = d_tsfset[26];
    float fiber_cladding = d_tsfset[27];
    float fiber_lambda = d_tsfset[28];


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
             else if (fiber) gcorr*=fibermult(i, j0-1, k, fiber_fxpos, fiber_fypos, fiber_radius, fiber_cutoff, fiber_core, fiber_cladding, dx, fiber_lambda);

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
             else if (fiber) gcorr*=fibermult(i, j1, k, fiber_fxpos, fiber_fypos, fiber_radius, fiber_cutoff, fiber_core, fiber_cladding, dx, fiber_lambda);

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
__global__  void
tsf_estep_cKernel(float *d_ex, float *d_ey, float *d_ez,
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
    int skip_k0 = (int)d_tsfset[14];
    int skip_kn = (int)d_tsfset[15];
    int gaussian = (int)d_tsfset[16];
    float gaussian_fxpos = d_tsfset[17];
    float gaussian_fypos = d_tsfset[18];
    float gaussian_rx = d_tsfset[19];
    float gaussian_ry = d_tsfset[20];
    int fiber = (int)d_tsfset[21];
    float fiber_fxpos = d_tsfset[22];
    float fiber_fypos = d_tsfset[23];
    float fiber_radius = d_tsfset[24];
    float fiber_cutoff = d_tsfset[25];
    float fiber_core = d_tsfset[26];
    float fiber_cladding = d_tsfset[27];
    float fiber_lambda = d_tsfset[28];



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


    
    if ((!skip_k0) && k==(k0))
    {
        if (i>=i0 && i<=i1 && j>=j0 && j<=j1)
        {
             if (gaussian) gcorr*=gaussmult(i, j, k0-1, gaussian_fxpos, gaussian_fypos, gaussian_rx, gaussian_ry);
             else if (fiber) gcorr*=fibermult(i, j, k0-1, fiber_fxpos, fiber_fypos, fiber_radius, fiber_cutoff, fiber_core, fiber_cladding, dx, fiber_lambda);
    
             d = corr*k_dcomp(i, j, k0-1,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             vh = gcorr*get_dval(d_jpool_h, jpool_size, d);

             if (i<i1) d_ex[pos] += dt/dx/epsilon*
                     k_ghy(vh, theta, phi, psi);                                         
             if (j<j1) d_ey[pos] -= dt/dx/epsilon*
                     k_ghx(vh, theta, phi, psi);                                         
        }
    }
    if ((!skip_kn) && k==k1)
    {
        if (i>=i0 && i<=i1 && j>=j0 && j<=j1)
        {
             if (gaussian) gcorr*=gaussmult(i, j, k1, gaussian_fxpos, gaussian_fypos, gaussian_rx, gaussian_ry);
             else if (fiber) gcorr*=fibermult(i, j, k1, fiber_fxpos, fiber_fypos, fiber_radius, fiber_cutoff, fiber_core, fiber_cladding, dx, fiber_lambda);

             d = corr*k_dcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             vh = gcorr*get_dval(d_jpool_h, jpool_size, d);

             if (i<i1) d_ex[pos] -= dt/dx/epsilon*
                     k_ghy(vh, theta, phi, psi);                                         
             if (j<j1) d_ey[pos] += dt/dx/epsilon*
                     k_ghx(vh, theta, phi, psi);                                         
        }
    }

}

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
    int fiber = (int)d_tsfset[21];
    float fiber_fxpos = d_tsfset[22];
    float fiber_fypos = d_tsfset[23];
    float fiber_radius = d_tsfset[24];
    float fiber_cutoff = d_tsfset[25];
    float fiber_core = d_tsfset[26];
    float fiber_cladding = d_tsfset[27];
    float fiber_lambda = d_tsfset[28];





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
             else if (fiber) gcorr*=fibermult(i, j, k, fiber_fxpos, fiber_fypos, fiber_radius, fiber_cutoff, fiber_core, fiber_cladding, dx, fiber_lambda);
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
             else if (fiber) gcorr*=fibermult(i, j, k, fiber_fxpos, fiber_fypos, fiber_radius, fiber_cutoff, fiber_core, fiber_cladding, dx, fiber_lambda);
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
    int fiber = (int)d_tsfset[21];
    float fiber_fxpos = d_tsfset[22];
    float fiber_fypos = d_tsfset[23];
    float fiber_radius = d_tsfset[24];
    float fiber_cutoff = d_tsfset[25];
    float fiber_core = d_tsfset[26];
    float fiber_cladding = d_tsfset[27];
    float fiber_lambda = d_tsfset[28];




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
              else if (fiber) gcorr*=fibermult(i, j, k, fiber_fxpos, fiber_fypos, fiber_radius, fiber_cutoff, fiber_core, fiber_cladding, dx, fiber_lambda);
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
              else if (fiber) gcorr*=fibermult(i, j, k, fiber_fxpos, fiber_fypos, fiber_radius, fiber_cutoff, fiber_core, fiber_cladding, dx, fiber_lambda);
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
__global__  void
tsf_hstep_cKernel(float *d_ex, float *d_ey, float *d_ez,
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
    int skip_k0 = (int)d_tsfset[14];
    int skip_kn = (int)d_tsfset[15];
    int gaussian = (int)d_tsfset[16];
    float gaussian_fxpos = d_tsfset[17];
    float gaussian_fypos = d_tsfset[18];
    float gaussian_rx = d_tsfset[19];
    float gaussian_ry = d_tsfset[20];
    int fiber = (int)d_tsfset[21];
    float fiber_fxpos = d_tsfset[22];
    float fiber_fypos = d_tsfset[23];
    float fiber_radius = d_tsfset[24];
    float fiber_cutoff = d_tsfset[25];
    float fiber_core = d_tsfset[26];
    float fiber_cladding = d_tsfset[27];
    float fiber_lambda = d_tsfset[28];




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



    if ((!skip_k0) && k==k0)
    {
        if (i>=i0 && i<=i1 && j>=j0 && j<=j1)
        {
             if (gaussian) gcorr*=gaussmult(i, j, k, gaussian_fxpos, gaussian_fypos, gaussian_rx, gaussian_ry);
              else if (fiber) gcorr*=fibermult(i, j, k, fiber_fxpos, fiber_fypos, fiber_radius, fiber_cutoff, fiber_core, fiber_cladding, dx, fiber_lambda);
            d = corr*k_dcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             ve = gcorr*get_dval(d_jpool_e, jpool_size, d);

             if (j<j1) d_hx[CAT(i, j, k-1)] -= dt/dx/mu*
                     k_gey(ve, theta, phi, psi);                                         
             if (i<i1) d_hy[CAT(i, j, k-1)] += dt/dx/mu*
                     k_gex(ve, theta, phi, psi);            
        }
    }
    if ((!skip_kn) && k==k1)
    {
        if (i>=i0 && i<=i1 && j>=j0 && j<=j1)
        {
             if (gaussian) gcorr*=gaussmult(i, j, k, gaussian_fxpos, gaussian_fypos, gaussian_rx, gaussian_ry);
              else if (fiber) gcorr*=fibermult(i, j, k, fiber_fxpos, fiber_fypos, fiber_radius, fiber_cutoff, fiber_core, fiber_cladding, dx, fiber_lambda);
            d = corr*k_dcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             ve = gcorr*get_dval(d_jpool_e, jpool_size, d);

             if (j<j1) d_hx[pos] += dt/dx/mu*
                     k_gey(ve, theta, phi, psi);                                         
             if (i<i1) d_hy[pos] -= dt/dx/mu*
                     k_gex(ve, theta, phi, psi);            
        }
    }

}

__global__  void
tsff_estep_aKernel(float *d_ex, float *d_ey, float *d_ez,
          float *d_hx, float *d_hy, float *d_hz,
          float *d_epsilon, float *d_mu, int matmode,
          float *d_tsfset,
          float *d_jpool_e, float *d_jpool_h, int jpool_size,
          float dx, float dt,
          int xres, int yres, int zres, int dir,
          float theta, float phi, float psi, float corr, float amplitude,
          int nwave, int ntot)
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
    float epsilon;
    int skip_i0 = (int)d_tsfset[10];
    int skip_in = (int)d_tsfset[11];

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
             d = corr*k_rdcomp(i0-1, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             vh = get_dval(d_jpool_h, ntot*jpool_size, d+nwave*jpool_size)*amplitude;

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
             d = corr*k_rdcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             vh = get_dval(d_jpool_h, ntot*jpool_size, d+nwave*jpool_size)*amplitude;

             if (j<j1) d_ey[pos] -= dt/dx/epsilon*
                     k_ghz(vh, theta, phi, psi);                                         
             if (k<k1) d_ez[pos] += dt/dx/epsilon*
                     k_ghy(vh, theta, phi, psi);                                         
        }
    }
}
__global__  void
tsff_estep_bKernel(float *d_ex, float *d_ey, float *d_ez,
          float *d_hx, float *d_hy, float *d_hz,
          float *d_epsilon, float *d_mu, int matmode,
          float *d_tsfset,
          float *d_jpool_e, float *d_jpool_h, int jpool_size,
          float dx, float dt,
          int xres, int yres, int zres, int dir,
          float theta, float phi, float psi, float corr, float amplitude,
          int nwave, int ntot)
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
    int skip_j0 = (int)d_tsfset[12];
    int skip_jn = (int)d_tsfset[13];


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
             d = corr*k_rdcomp(i, j0-1, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             vh = get_dval(d_jpool_h, ntot*jpool_size, d+nwave*jpool_size)*amplitude;

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
             d = corr*k_rdcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             vh = get_dval(d_jpool_h, ntot*jpool_size, d+nwave*jpool_size)*amplitude;

             if (i<i1) d_ex[pos] += dt/dx/epsilon*
                     k_ghz(vh, theta, phi, psi);                                         
             if (k<k1) d_ez[pos] -= dt/dx/epsilon*
                     k_ghx(vh, theta, phi, psi);                                         
        }
    }

}

__global__  void
tsff_estep_cKernel(float *d_ex, float *d_ey, float *d_ez,
          float *d_hx, float *d_hy, float *d_hz,
          float *d_epsilon, float *d_mu, int matmode,
          float *d_tsfset,
          float *d_jpool_e, float *d_jpool_h, int jpool_size,
          float dx, float dt,
          int xres, int yres, int zres, int dir,
          float theta, float phi, float psi, float corr, float amplitude,
          int nwave, int ntot)
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
    int skip_k0 = (int)d_tsfset[14];
    int skip_kn = (int)d_tsfset[15];

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

   
    if ((!skip_k0) && k==(k0))
    {
        if (i>=i0 && i<=i1 && j>=j0 && j<=j1)
        {
             d = corr*k_rdcomp(i, j, k0-1,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             vh = get_dval(d_jpool_h, ntot*jpool_size, d+nwave*jpool_size)*amplitude;

             if (i<i1) d_ex[pos] += dt/dx/epsilon*
                     k_ghy(vh, theta, phi, psi);                                         
             if (j<j1) d_ey[pos] -= dt/dx/epsilon*
                     k_ghx(vh, theta, phi, psi);                                         
        }
    }

    if ((!skip_kn) && k==k1)
    {
        if (i>=i0 && i<=i1 && j>=j0 && j<=j1)
        {
             d = corr*k_rdcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             vh = get_dval(d_jpool_h, ntot*jpool_size, d+nwave*jpool_size)*amplitude;

             if (i<i1) d_ex[pos] -= dt/dx/epsilon*
                     k_ghy(vh, theta, phi, psi);                                         
             if (j<j1) d_ey[pos] += dt/dx/epsilon*
                     k_ghx(vh, theta, phi, psi);                                         
        }
    }

}

__global__  void
tsff_hstep_aKernel(float *d_ex, float *d_ey, float *d_ez,
          float *d_hx, float *d_hy, float *d_hz,
          float *d_epsilon, float *d_mu, int matmode,
          float *d_tsfset,
          float *d_jpool_e, float *d_jpool_h, int jpool_size,
          float dx, float dt,
          int xres, int yres, int zres, int dir,
          float theta, float phi, float psi, float corr, float amplitude,
          int nwave, int ntot)
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
    int skip_i0 = (int)d_tsfset[10];
    int skip_in = (int)d_tsfset[11];



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
             d = corr*k_rdcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
             
             ve = get_dval(d_jpool_e, ntot*jpool_size, d+nwave*jpool_size)*amplitude;

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
             d = corr*k_rdcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             ve = get_dval(d_jpool_e, ntot*jpool_size, d+nwave*jpool_size)*amplitude;

             if (k<k1) d_hy[pos] += dt/dx/mu*
                     k_gez(ve, theta, phi, psi);                                         
             if (j<j1) d_hz[pos] -= dt/dx/mu*
                     k_gey(ve, theta, phi, psi);            
        }
    }

}
__global__  void
tsff_hstep_bKernel(float *d_ex, float *d_ey, float *d_ez,
          float *d_hx, float *d_hy, float *d_hz, 
          float *d_epsilon, float *d_mu, int matmode,
          float *d_tsfset,
          float *d_jpool_e, float *d_jpool_h, int jpool_size,
          float dx, float dt,
          int xres, int yres, int zres, int dir,
          float theta, float phi, float psi, float corr, float amplitude,
          int nwave, int ntot)
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
    int skip_j0 = (int)d_tsfset[12];
    int skip_jn = (int)d_tsfset[13];


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
             d = corr*k_rdcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             ve = get_dval(d_jpool_e, ntot*jpool_size, d+nwave*jpool_size)*amplitude;

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
             d = corr*k_rdcomp(i, j, k,
                          xres, yres, zres,
                          theta, phi,
                          0, xres, 0, yres, 0, zres);
              
             ve = get_dval(d_jpool_e, ntot*jpool_size, d+nwave*jpool_size)*amplitude;

             if (k<k1) d_hx[pos] -= dt/dx/mu*
                     k_gez(ve, theta, phi, psi);                                         
             if (i<i1) d_hz[pos] += dt/dx/mu*
                     k_gex(ve, theta, phi, psi);            
        }
    }

}

__global__  void
tsff_hstep_cKernel(float *d_ex, float *d_ey, float *d_ez,
          float *d_hx, float *d_hy, float *d_hz,
          float *d_epsilon, float *d_mu, int matmode,
          float *d_tsffset,
          float *d_jpool_e, float *d_jpool_h, int jpool_size,
          float dx, float dt,
          int xres, int yres, int zres, int dir,
          float theta, float phi, float psi, float corr, float amplitude,
          int nwave, int ntot)
{
    int i, j, k;
    long int pos;
    float d, ve, mu;    
    int i0 = (int)d_tsffset[0];
    int j0 = (int)d_tsffset[1];
    int k0 = (int)d_tsffset[2];
    int i1 = (int)d_tsffset[3];
    int j1 = (int)d_tsffset[4];
    int k1 = (int)d_tsffset[5];
    int skip_k0 = (int)d_tsffset[14];
    int skip_kn = (int)d_tsffset[15];

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


    if ((!skip_k0) && k==k0)
    {
	    if (i>=i0 && i<=i1 && j>=j0 && j<=j1)
	    {
		    d = corr*k_rdcomp(i, j, k,
				    xres, yres, zres,
				    theta, phi,
				    0, xres, 0, yres, 0, zres);

		    ve = get_dval(d_jpool_e, ntot*jpool_size, d+nwave*jpool_size)*amplitude;


		    if (j<j1) d_hx[CAT(i, j, k-1)] -= dt/dx/mu*
			    k_gey(ve, theta, phi, psi);                                         
		    if (i<i1) d_hy[CAT(i, j, k-1)] += dt/dx/mu*
			    k_gex(ve, theta, phi, psi);            
	    }
    }
    if ((!skip_kn) && k==k1)
    {
	    if (i>=i0 && i<=i1 && j>=j0 && j<=j1)
	    {
		    d = corr*k_rdcomp(i, j, k,
				    xres, yres, zres,
				    theta, phi,
				    0, xres, 0, yres, 0, zres);

		    ve = get_dval(d_jpool_e, ntot*jpool_size, d+nwave*jpool_size)*amplitude;

		    if (j<j1) d_hx[pos] += dt/dx/mu*
			    k_gey(ve, theta, phi, psi);                                         
		    if (i<i1) d_hy[pos] -= dt/dx/mu*
			    k_gex(ve, theta, phi, psi);            
	    }
    }


}

//tsff kernel: slo by predelat to na (xres x yres x nphi), iterovat se bude pres ntheta,
//ale mozna bude efektivnejsi nechat to na (xres x yres x zres) a iterovat pres oboji jako u cpu varianty

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
	else if (matmode == 1 || matmode==2 || matmode==3) { /*= mattype 0 = linear material given pixel by pixel*/
		if (nmat==0 || d_mat[pos] == 0) {
			if (matmode == 2) *epsilon = EPSILON_0;
			else *epsilon = d_epsilon[pos];
			if (matmode == 1) *mu = MU_0;
			*mu = d_mu[pos];
		} else { //some tabulated material inside linear material
			*epsilon = d_mattab[GETEPS(d_mat[pos])]; 
			*mu = d_mattab[GETMU(d_mat[pos])];
		}
	} else if (nmat>0 && mattype==1) { /*tabulated linear material*/
		*epsilon = d_mattab[GETEPS(d_mat[pos])]; 
		*mu = d_mattab[GETMU(d_mat[pos])];
	} else {
		*epsilon = EPSILON_0;
		*mu = MU_0;
	}

}


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


__device__ void 
prod(float t[3][3], float n[3], float f[3])
{
    f[0] = t[0][0]*n[0] + t[0][1]*n[1] + t[0][2]*n[2];
    f[1] = t[1][0]*n[0] + t[1][1]*n[1] + t[1][2]*n[2];
    f[2] = t[2][0]*n[0] + t[2][1]*n[1] + t[2][2]*n[2];
}




__device__ void 
assemble_matrix(float t[3][3], float ex, float ey, float ez, 
                float hx, float hy, float hz, float epsilon, float mu)
{
    float val;

    val = (epsilon*(ex*ex+ey*ey+ez*ez) + mu*(hx*hx+hy*hy+hz*hz))/2.0;
    t[0][0] = epsilon*ex*ex + mu*hx*hx - val;
    t[1][1] = epsilon*ey*ey + mu*hy*hy - val;
    t[2][2] = epsilon*ez*ez + mu*hz*hz - val;
    t[0][1] = t[1][0] = epsilon*ex*ey + mu*hx*hy;
    t[0][2] = t[2][0] = epsilon*ex*ez + mu*hx*hz;
    t[1][2] = t[2][1] = epsilon*ey*ez + mu*hy*hz;

}


__global__  void
plrcnullKernel(float *d_plrcx, float *d_plrcy, float *d_plrcz, 
               int xres, int yres, int zres)
{
     long int k;

     for (k=0; k<(7*xres*yres*zres); k++) 
          d_plrcx[k] = d_plrcy[k] = d_plrcz[k] = 0;

          
}
__global__  void
adenullKernel(float *d_pxp, float *d_px, float *d_dpxp, float *d_dpx, 
               float *d_exp, float *d_eyp, float *d_ezp,
               int xres, int yres, int zres)
{
     long int k;

     for (k=0; k<(3*xres*yres*zres); k++)  
          d_dpxp[k] = d_dpx[k] = 0;

     for (k=0; k<(6*xres*yres*zres); k++)  
          d_pxp[k] = d_px[k] = 0;

     for (k=0; k<(xres*yres*zres); k++)  
          d_exp[k] = d_eyp[k] = d_ezp[k] = 0;

          
}

__global__  void
forceKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
          int xres, int yres, int zres, 
          float *d_epsilon, float *d_mu, float *d_forces, 
          int *i0, int *i1, int *j0, int *j1, int *k0, int *k1,
          int nsteps, int step,
          int *d_mat, int *d_mattype, float *d_mattab, int matmode, int nmat)
{
    int m = threadIdx.x;
    int i, j, k;
    float faccx, faccy, faccz;
    float t[3][3];
    float normal[3], f[3];
    float epsilon, mu, ex, ey, ez, hx, hy, hz;
    long int pos;

    faccx = faccy = faccz = 0;
    for (j=j0[m]; j<j1[m]; j++)
    {
	    for (k=k0[m]; k<k1[m]; k++)
	    {
		    pos = CAT(i0[m],j,k);
		    get_matprops(i0[m], j, k, &epsilon, &mu, 
				    d_mat, d_mattype, d_mattab, matmode, nmat, d_epsilon, d_mu, xres, yres, zres);

		    ex = d_ex[pos];
		    ey = d_ey[pos];
		    ez = d_ez[pos];
		    hx = d_hx[pos];
		    hy = d_hy[pos];
		    hz = d_hz[pos];

		    assemble_matrix(t, ex, ey, ez, hx, hy, hz, epsilon, mu);

		    normal[0] = 1; normal[1] = 0; normal[2] = 0;
		    prod(t, normal, f);

		    faccx += f[0];
		    faccy += f[1];
		    faccz += f[2];

		    pos = CAT(i1[m],j,k);
		    get_matprops(i1[m], j, k, &epsilon, &mu, 
				    d_mat, d_mattype, d_mattab, matmode, nmat, d_epsilon, d_mu, xres, yres, zres);

		    ex = d_ex[pos];
		    ey = d_ey[pos];
		    ez = d_ez[pos];
		    hx = d_hx[pos];
		    hy = d_hy[pos];
		    hz = d_hz[pos];

		    assemble_matrix(t, ex, ey, ez, hx, hy, hz, epsilon, mu);
		    normal[0] = -1; normal[1] = 0; normal[2] = 0;
		    prod(t, normal, f);
		    faccx += f[0];
		    faccy += f[1];
		    faccz += f[2];
	    }
    }
    for (i=i0[m]; i<i1[m]; i++)
    {
	    for (k=k0[m]; k<k1[m]; k++)
	    {
		    /*Local values of ex, ey, ez at x0. Normal of this facet is 1 0 0*/
		    /*T_ij*/

		    pos = CAT(i, j0[m], k);
		    get_matprops(i, j0[m], k, &epsilon, &mu, 
				    d_mat, d_mattype, d_mattab, matmode, nmat, d_epsilon, d_mu, xres, yres, zres);

		    ex = d_ex[pos];
		    ey = d_ey[pos];
		    ez = d_ez[pos];
		    hx = d_hx[pos];
		    hy = d_hy[pos];
		    hz = d_hz[pos];

		    assemble_matrix(t, ex, ey, ez, hx, hy, hz, epsilon, mu);
		    normal[0] = 0; normal[1] = 1; normal[2] = 0;
		    prod(t, normal, f);
		    faccx += f[0];
		    faccy += f[1];
		    faccz += f[2];

		    pos = CAT(i, j1[m], k);
		    get_matprops(i, j1[m], k, &epsilon, &mu, 
				    d_mat, d_mattype, d_mattab, matmode, nmat, d_epsilon, d_mu, xres, yres, zres);

		    ex = d_ex[pos];
		    ey = d_ey[pos];
		    ez = d_ez[pos];
		    hx = d_hx[pos];
		    hy = d_hy[pos];
		    hz = d_hz[pos];

		    assemble_matrix(t, ex, ey, ez, hx, hy, hz, epsilon, mu);
		    normal[0] = 0; normal[1] = -1; normal[2] = 0;
		    prod(t, normal, f);
		    faccx += f[0];
		    faccy += f[1];
		    faccz += f[2];

	    }
    }
    for (i=i0[m]; i<i1[m]; i++)
    {
	    for (j=j0[m]; j<j1[m]; j++)
	    {
		    pos = CAT(i, j, k0[m]);
		    get_matprops(i, j, k0[m], &epsilon, &mu, 
				    d_mat, d_mattype, d_mattab, matmode, nmat, d_epsilon, d_mu, xres, yres, zres);

		    ex = d_ex[pos];
		    ey = d_ey[pos];
		    ez = d_ez[pos];
		    hx = d_hx[pos];
		    hy = d_hy[pos];
		    hz = d_hz[pos];

		    assemble_matrix(t, ex, ey, ez, hx, hy, hz, epsilon, mu);
		    normal[0] = 0; normal[1] = 0; normal[2] = 1;
		    prod(t, normal, f);
		    faccx += f[0];
		    faccy += f[1];
		    faccz += f[2];

		    pos = CAT(i, j, k1[m]);
		    get_matprops(i, j, k1[m], &epsilon, &mu, 
				    d_mat, d_mattype, d_mattab, matmode, nmat, d_epsilon, d_mu, xres, yres, zres);

		    ex = d_ex[pos];
		    ey = d_ey[pos];
		    ez = d_ez[pos];
		    hx = d_hx[pos];
		    hy = d_hy[pos];
		    hz = d_hz[pos];

		    assemble_matrix(t, ex, ey, ez, hx, hy, hz, epsilon, mu);
		    normal[0] = 0; normal[1] = 0; normal[2] = -1;
		    prod(t, normal, f);
		    faccx += f[0];
		    faccy += f[1];
		    faccz += f[2];

	    }
    }

    d_forces[3*(m*nsteps + step)] = faccx;
    d_forces[3*(m*nsteps + step)+1] = faccy;
    d_forces[3*(m*nsteps + step)+2] = faccz;

}

__global__  void
outpointKernel(float *d_ex, float *d_ey, float *d_ez, float *d_hx, float *d_hy, float *d_hz,
                              float *d_outpointdata, int *d_outpoint_pos, int nsteps, int noutpoints,
                              int step, int xres, int yres, int zres)
{
    int i;
    for (i=0; i<noutpoints; i++) 
    {
        d_outpointdata[6*(i*nsteps + step)] = d_ex[CAT(d_outpoint_pos[3*i], d_outpoint_pos[3*i+1], d_outpoint_pos[3*i+2])];
        d_outpointdata[6*(i*nsteps + step) + 1] = d_ey[CAT(d_outpoint_pos[3*i], d_outpoint_pos[3*i+1], d_outpoint_pos[3*i+2])];
        d_outpointdata[6*(i*nsteps + step) + 2] = d_ez[CAT(d_outpoint_pos[3*i], d_outpoint_pos[3*i+1], d_outpoint_pos[3*i+2])];
        d_outpointdata[6*(i*nsteps + step) + 3] = d_hx[CAT(d_outpoint_pos[3*i], d_outpoint_pos[3*i+1], d_outpoint_pos[3*i+2])];
        d_outpointdata[6*(i*nsteps + step) + 4] = d_hy[CAT(d_outpoint_pos[3*i], d_outpoint_pos[3*i+1], d_outpoint_pos[3*i+2])];
        d_outpointdata[6*(i*nsteps + step) + 5] = d_hz[CAT(d_outpoint_pos[3*i], d_outpoint_pos[3*i+1], d_outpoint_pos[3*i+2])];
    }
}



/*
__global__  void
fiberKernel(float *d_ex, float *d_ey, float *d_ez,
          float *d_hx, float *d_hy, float *d_hz,
          float *d_fiberset,
          float *d_jpool_e, float *d_jpool_h, int jpool_size,
          float dx, float dt,
          int xres, int yres, int zres, int dir)
{
    int i, j, k;
    long int pos;
    int i0 = (int)d_fiberset[0];
    int j0 = (int)d_fiberset[1];
    int k0 = (int)d_fiberset[2];
    int ij = (int)d_fiberset[3];
    int jk = (int)d_fiberset[4];
    float lambda = d_fiberset[5];
    float radius = d_fiberset[6];
    float core = d_fiberset[7];
    float cladding = d_fiberset[8];
    int ijsize = 20; //FIXME
    int jksize = 20;
    float omegazero, eval, hval, radival, radval, v;
    
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

    v = 2*CU_PI/lambda*radius*dx*(core - cladding);
    omegazero = (0.65 + 1.619/v/sqrt(v) + 2.879/pow(v,6))*radius*dx;

    if (i0==-1 && j0==-1) {
         if (k==k0 && i>ij && j>jk && i<(ij+ijsize) && j<(jk+jksize)) {
                radival = sqrt((float)((i-ij-ijsize/2)*(i-ij-ijsize/2) + (j-jk-jksize/2)*(j-jk-jksize/2)));
                radval = radival*dx;

                eval = dt/dx/EPSILON_0*d_jpool_h[10];
                hval = dt/dx/MU_0*d_jpool_e[10];

                if (radival < radius) {
		    d_ex[pos] = eval*exp(-(radval*radval/omegazero/omegazero));
                    d_ey[pos] = -eval*exp(-(radval*radval/omegazero/omegazero));
                    d_ez[pos] = 0;
                    d_hx[pos] = -hval*exp(-(radval*radval/omegazero/omegazero));
                    d_hy[pos] = hval*exp(-(radval*radval/omegazero/omegazero));
                    d_hz[pos] = 0; 
                }
         }
    }

}
*/
 
cudaError_t wrap_hKernel(SvGpuPlan *plan)
{
        cudaError_t err;

        int ga, gb, b, size;
        int xsize, ysize; 

        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 
        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);
        dim3 dimpBlock(plan->maxthreads);

        size = (int)ceil(((float)plan->xres*plan->yres*plan->zres)/(float)plan->maxthreads); //total number of thread blocks necessary
        xsize = plan->xres;
        ysize = (int)ceil((float)size/((float)plan->xres));
        dim3 dimpGrid(xsize, ysize);


        if (plan->bndx0==3 || plan->bndxn==3 || plan->bndy0==3 || plan->bndyn==3 || plan->bndz0==3 || plan->bndzn==3)
        {
            if (plan->nmat>0)
               hKernel_cpml_tab<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                                                plan->d_epsilon, plan->d_mu, plan->d_sigma, plan->d_sigast, 
                                                plan->d_mat, plan->d_mattype, plan->d_mattab, plan->nmat, plan->matmode,
                                                plan->d_bnds, plan->d_depths,
                                                plan->d_cpml_kappah_x0, plan->d_cpml_kappah_xn, plan->d_cpml_kappah_y0, plan->d_cpml_kappah_yn, plan->d_cpml_kappah_z0, plan->d_cpml_kappah_zn,
                                                plan->xres, plan->yres, plan->zres,
                                                plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);


            else
               hKernel_cpml<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                                                plan->d_epsilon, plan->d_mu, plan->d_sigma, plan->d_sigast, plan->matmode,
                                                plan->d_bnds, plan->d_depths,
                                                plan->d_cpml_kappah_x0, plan->d_cpml_kappah_xn, plan->d_cpml_kappah_y0, plan->d_cpml_kappah_yn, plan->d_cpml_kappah_z0, plan->d_cpml_kappah_zn,
                                                plan->xres, plan->yres, plan->zres,
                                                plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);

        }
        else if (plan->nmat>0)
               hKernel_tab<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                                                plan->d_epsilon, plan->d_mu, plan->d_sigma, plan->d_sigast, 
                                                plan->d_mat, plan->d_mattype, plan->d_mattab, plan->nmat, plan->matmode,
                                                plan->xres, plan->yres, plan->zres,
                                                plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);

        else if (plan->matmode == 0 || plan->matmode == 1) {

               hKernel_none<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                                              plan->xres, plan->yres, plan->zres,
                                              plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);
          }
        
        else  {
   	   phKernel<<<dimpGrid, dimpBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                                              plan->d_epsilon, plan->d_mu, plan->d_sigma, plan->d_sigast, 
                                              plan->xres, plan->yres, plan->zres,
                                              plan->dx, plan->dy, plan->dz, plan->dt, plan->dir, plan->maxthreads, ysize);
    	   /*hKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                                              plan->d_epsilon, plan->d_mu, plan->d_sigma, plan->d_sigast, 
                                              plan->xres, plan->yres, plan->zres,
                                              plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);*/
        }
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

        int ga, gb, b, size, xsize, ysize;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);
        dim3 dimpBlock(plan->maxthreads);

        size = (int)ceil(((float)plan->xres*plan->yres*plan->zres)/(float)plan->maxthreads); //total number of thread blocks necessary
        xsize = plan->xres;
        ysize = (int)ceil((float)size/((float)plan->xres));
        dim3 dimpGrid(xsize, ysize);

        if (plan->h_isplrc && plan->step==0) {
             printf("Initialising plrc fields");
             plrcnullKernel<<<1, 1>>>(plan->d_plrcx, plan->d_plrcy, plan->d_plrcz,
                                     plan->xres, plan->yres, plan->zres);
        } 
        if (plan->h_isade && plan->step==0) {
             printf("Initialising ade fields");
             adenullKernel<<<1, 1>>>(plan->d_px, plan->d_pxp, plan->d_dpx, plan->d_dpxp, plan->d_exp, plan->d_eyp, plan->d_ezp,
                                     plan->xres, plan->yres, plan->zres);
        } 
       


        //printf("E Launch %d %d %d %g %g %g %g\n", plan->xres, plan->yres, plan->zres, plan->dx, plan->dy, plan->dz, plan->dt);

        if (plan->bndx0==3 || plan->bndxn==3 || plan->bndy0==3 || plan->bndyn==3 || plan->bndz0==3 || plan->bndzn==3)
	{
              if (plan->nmat>0)
		eKernel_cpml_tab<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
				plan->d_epsilon, plan->d_mu, plan->d_sigma, plan->d_sigast, 
                                plan->d_mat, plan->d_mattype, plan->d_mattab, plan->nmat, plan->matmode,
                                plan->d_bnds, plan->d_depths,
                                plan->d_cpml_kappae_x0, plan->d_cpml_kappae_xn, plan->d_cpml_kappae_y0, plan->d_cpml_kappae_yn, plan->d_cpml_kappae_z0, plan->d_cpml_kappae_zn,
				plan->xres, plan->yres, plan->zres,
				plan->dx, plan->dy, plan->dz, plan->dt, plan->dir,
                                plan->d_sf_jpool_e, plan->d_sfset, plan->sf_jpool_size,
                                plan->d_plrcx, plan->d_plrcy, plan->d_plrcz,
                                plan->d_exp, plan->d_eyp, plan->d_ezp,
                                plan->d_px, plan->d_pxp, plan->d_dpx, plan->d_dpxp);
              else
		eKernel_cpml<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
				plan->d_epsilon, plan->d_mu, plan->d_sigma, plan->d_sigast, 
                                plan->d_bnds, plan->d_depths,
                                plan->d_cpml_kappae_x0, plan->d_cpml_kappae_xn, plan->d_cpml_kappae_y0, plan->d_cpml_kappae_yn, plan->d_cpml_kappae_z0, plan->d_cpml_kappae_zn,
				plan->xres, plan->yres, plan->zres,
				plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);
        
	}
	else if (plan->nmat>0) {
		eKernel_tab<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz, 
                                plan->d_epsilon, plan->d_mu, plan->d_sigma, plan->d_sigast,
				plan->d_mat, plan->d_mattype, plan->d_mattab, plan->nmat, plan->matmode,
				plan->xres, plan->yres, plan->zres,
				plan->dx, plan->dy, plan->dz, plan->dt, plan->dir, plan->d_sf_jpool_e, plan->d_sfset, plan->sf_jpool_size,
                                plan->d_plrcx, plan->d_plrcy, plan->d_plrcz,
                                plan->d_exp, plan->d_eyp, plan->d_ezp,
                                plan->d_px, plan->d_pxp, plan->d_dpx, plan->d_dpxp);
	}
	else if (plan->matmode == 0 || plan->matmode == 2)  {
		eKernel_none<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz, 
				plan->xres, plan->yres, plan->zres,
				plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);
	}
	else {
		peKernel<<<dimpGrid, dimpBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
				plan->d_epsilon, plan->d_mu, plan->d_sigma, plan->d_sigast, 
				plan->xres, plan->yres, plan->zres,
				plan->dx, plan->dy, plan->dz, plan->dt, plan->dir, plan->maxthreads, ysize);
		/*eKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
				plan->d_epsilon, plan->d_mu, plan->d_sigma, plan->d_sigast, 
				plan->xres, plan->yres, plan->zres,
				plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);*/
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
/*
cudaError_t wrap_yKernel(TGPUplan *plan)
{
        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

	hKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                                              plan->d_epsilon, plan->d_mu, plan->d_sigma, plan->d_sigast, 
                                              plan->xres, plan->yres, plan->zres,
                                              plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);

        cudaThreadSynchronize();
   	eKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                                              plan->d_epsilon, plan->d_mu, plan->d_sigma, plan->d_sigast, 
                                              plan->xres, plan->yres, plan->zres,
                                              plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);

        cudaThreadSynchronize();
         return cudaGetLastError();
}
*/
cudaError_t wrap_liao_cpybnd(SvGpuPlan *plan)
{
        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

	liaocpyKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                                              plan->d_liao_x0, plan->d_liao_xn, plan->d_liao_y0, plan->d_liao_yn, 
                                              plan->d_liao_z0, plan->d_liao_zn,
                                              plan->xres, plan->yres, plan->zres, 
                                              plan->d_bnds,
                                              plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);

        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_liao_applybnd(SvGpuPlan *plan)
{
        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

	liaorunKernelx<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz, plan->d_epsilon,
                                              plan->d_liao_x0, plan->d_liao_xn, plan->d_liao_y0, plan->d_liao_yn, 
                                              plan->d_liao_z0, plan->d_liao_zn,
                                              plan->xres, plan->yres, plan->zres, 
                                              plan->d_bnds, plan->matmode,
                                              plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);
        cudaThreadSynchronize();
	liaorunKernely<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz, plan->d_epsilon,
                                              plan->d_liao_x0, plan->d_liao_xn, plan->d_liao_y0, plan->d_liao_yn, 
                                              plan->d_liao_z0, plan->d_liao_zn,
                                              plan->xres, plan->yres, plan->zres, 
                                              plan->d_bnds, plan->matmode,
                                              plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);
        cudaThreadSynchronize();
	liaorunKernelz<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz, plan->d_epsilon,
                                              plan->d_liao_x0, plan->d_liao_xn, plan->d_liao_y0, plan->d_liao_yn, 
                                              plan->d_liao_z0, plan->d_liao_zn,
                                              plan->xres, plan->yres, plan->zres, 
                                              plan->d_bnds, plan->matmode,
                                              plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);

        cudaThreadSynchronize();
        return cudaGetLastError();
}


cudaError_t wrap_cpml_estep(SvGpuPlan *plan)
{

        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        }

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

        cpmleKernelx<<<dimGrid, dimBlock>>>(
            plan->d_ey, plan->d_ez, plan->d_hy, plan->d_hz, plan->d_epsilon, plan->d_sigma,
            plan->d_mat, plan->d_mattype, plan->d_mattab, plan->nmat,
            plan->d_cpml_p_x0, plan->d_cpml_p_xn,
            plan->d_cpml_be_x0, plan->d_cpml_ce_x0, plan->d_cpml_be_xn, plan->d_cpml_ce_xn,
            plan->xres, plan->yres, plan->zres, plan->d_bnds, plan->d_depths,
            plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);


        cudaThreadSynchronize();

        cpmleKernely<<<dimGrid, dimBlock>>>(
            plan->d_ex, plan->d_ez, plan->d_hx, plan->d_hz, plan->d_epsilon, plan->d_sigma,
            plan->d_mat, plan->d_mattype, plan->d_mattab, plan->nmat,
            plan->d_cpml_p_y0, plan->d_cpml_p_yn,
            plan->d_cpml_be_y0, plan->d_cpml_ce_y0, plan->d_cpml_be_yn, plan->d_cpml_ce_yn,
            plan->xres, plan->yres, plan->zres, plan->d_bnds, plan->d_depths,
            plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);


        cudaThreadSynchronize();

        cpmleKernelz<<<dimGrid, dimBlock>>>(
            plan->d_ex, plan->d_ey, plan->d_hx, plan->d_hy, plan->d_epsilon, plan->d_sigma,
            plan->d_mat, plan->d_mattype, plan->d_mattab, plan->nmat,
            plan->d_cpml_p_z0, plan->d_cpml_p_zn,
            plan->d_cpml_be_z0, plan->d_cpml_ce_z0, plan->d_cpml_be_zn, plan->d_cpml_ce_zn,
            plan->xres, plan->yres, plan->zres, plan->d_bnds, plan->d_depths,
            plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);


        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_cpml_hstep(SvGpuPlan *plan)
{
        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        }

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);
        
        cpmlhKernelx<<<dimGrid, dimBlock>>>(
            plan->d_ey, plan->d_ez, plan->d_hy, plan->d_hz, plan->d_mu, plan->d_sigast, plan->matmode,
            plan->d_cpml_p_x0, plan->d_cpml_p_xn,
            plan->d_cpml_bh_x0, plan->d_cpml_ch_x0, plan->d_cpml_bh_xn, plan->d_cpml_ch_xn,
            plan->xres, plan->yres, plan->zres, plan->d_bnds, plan->d_depths, 
            plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);

        cudaThreadSynchronize();
        cpmlhKernely<<<dimGrid, dimBlock>>>(
            plan->d_ex, plan->d_ez, plan->d_hx, plan->d_hz, plan->d_mu, plan->d_sigast, plan->matmode,
            plan->d_cpml_p_y0, plan->d_cpml_p_yn,
            plan->d_cpml_bh_y0, plan->d_cpml_ch_y0, plan->d_cpml_bh_yn, plan->d_cpml_ch_yn,
            plan->xres, plan->yres, plan->zres, plan->d_bnds, plan->d_depths, 
            plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);

        cudaThreadSynchronize();

        cpmlhKernelz<<<dimGrid, dimBlock>>>(
            plan->d_ex, plan->d_ey, plan->d_hx, plan->d_hy, plan->d_mu, plan->d_sigast, plan->matmode,
            plan->d_cpml_p_z0, plan->d_cpml_p_zn,
            plan->d_cpml_bh_z0, plan->d_cpml_ch_z0, plan->d_cpml_bh_zn, plan->d_cpml_ch_zn,
            plan->xres, plan->yres, plan->zres, plan->d_bnds, plan->d_depths, 
            plan->dx, plan->dy, plan->dz, plan->dt, plan->dir);

        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_mbndx_e_Kernel(SvGpuPlan *plan)
{
        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

        mbnxeKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                              plan->xres, plan->yres, plan->zres, plan->mb_bx0, plan->mb_bxn, 
                              plan->mb_bx0pos, plan->mb_bxnpos, 
                              plan->mb_by0pos, plan->mb_bynpos,
                              plan->mb_bz0pos, plan->mb_bznpos, plan->dir);

        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_mbndy_e_Kernel(SvGpuPlan *plan)
{
        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

        mbnyeKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                              plan->xres, plan->yres, plan->zres, plan->mb_by0, plan->mb_byn, 
                              plan->mb_bx0pos, plan->mb_bxnpos, 
                              plan->mb_by0pos, plan->mb_bynpos,
                              plan->mb_bz0pos, plan->mb_bznpos, plan->dir);

        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_mbndz_e_Kernel(SvGpuPlan *plan)
{
        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

        mbnzeKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                              plan->xres, plan->yres, plan->zres, plan->mb_bz0, plan->mb_bzn, 
                              plan->mb_bx0pos, plan->mb_bxnpos, 
                              plan->mb_by0pos, plan->mb_bynpos,
                              plan->mb_bz0pos, plan->mb_bznpos, plan->dir);

        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_mbndx_h_Kernel(SvGpuPlan *plan)
{
        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

        mbnxhKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                              plan->xres, plan->yres, plan->zres, plan->mb_bx0, plan->mb_bxn, 
                              plan->mb_bx0pos, plan->mb_bxnpos, 
                              plan->mb_by0pos, plan->mb_bynpos,
                              plan->mb_bz0pos, plan->mb_bznpos, plan->dir);

        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_mbndy_h_Kernel(SvGpuPlan *plan)
{
        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

        mbnyhKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                              plan->xres, plan->yres, plan->zres, plan->mb_by0, plan->mb_byn, 
                              plan->mb_bx0pos, plan->mb_bxnpos, 
                              plan->mb_by0pos, plan->mb_bynpos,
                              plan->mb_bz0pos, plan->mb_bznpos, plan->dir);

        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_mbndz_h_Kernel(SvGpuPlan *plan)
{
        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

        mbnzhKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                              plan->xres, plan->yres, plan->zres, plan->mb_bz0, plan->mb_bzn, 
                              plan->mb_bx0pos, plan->mb_bxnpos, 
                              plan->mb_by0pos, plan->mb_bynpos,
                              plan->mb_bz0pos, plan->mb_bznpos, plan->dir);

        cudaThreadSynchronize();
        return cudaGetLastError();
}
cudaError_t wrap_srceKernel(SvGpuPlan *plan, int i, int j, int k, float ex, float ey, float ez)
{
        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

        srcepointKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->xres, plan->yres, plan->zres, i, j, k, ex, ey, ez, plan->dir);
        cudaThreadSynchronize();
         
        return cudaGetLastError();
}


cudaError_t wrap_fastsumKernel(SvGpuPlan *plan)
{
        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->nsums;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->nsums;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->nsums;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

        fastsumKernel_a<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->xres, plan->yres, plan->zres,
                                             plan->d_epsilon, plan->d_mu, plan->d_sum_epsilon, plan->d_sum_sigma, plan->d_sum_mode,
                                             plan->d_sum_i0, plan->d_sum_i1, plan->d_sum_j0, plan->d_sum_j1, plan->d_sum_k0, plan->d_sum_k1,
                                             plan->nsums, plan->d_sum_accumulator, plan->dir,
                                             plan->d_mat, plan->d_mattype, plan->d_mattab, plan->matmode, plan->nmat);

        cudaThreadSynchronize();

        dim3 dimBlockb(plan->nsums);
        fastsumKernel_b<<<1, dimBlockb>>>(plan->d_sum_accumulator, plan->d_sums, plan->nsums, plan->nsteps, plan->step, plan->xres, plan->yres, plan->zres,
                                  plan->dir);

        cudaThreadSynchronize();
        return cudaGetLastError();

}



cudaError_t wrap_sumKernel(SvGpuPlan *plan)
{
        dim3 dimBlock(plan->nsums);

        sumKernel<<<1, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->xres, plan->yres, plan->zres, 
                                   plan->d_epsilon, plan->d_mu, plan->d_sums, plan->d_sum_epsilon, plan->d_sum_sigma, plan->d_sum_mode, 
                                   plan->d_sum_i0, plan->d_sum_i1, plan->d_sum_j0, plan->d_sum_j1, plan->d_sum_k0, plan->d_sum_k1,
                                   plan->nsteps, plan->step,
                                   plan->d_mat, plan->d_mattype, plan->d_mattab, plan->matmode, plan->nmat);
        cudaThreadSynchronize();
         
        return cudaGetLastError();
}

cudaError_t wrap_forceKernel(SvGpuPlan *plan)
{
        dim3 dimBlock(plan->nforces);

        forceKernel<<<1, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz, plan->xres, plan->yres, plan->zres, 
                                   plan->d_epsilon, plan->d_mu, plan->d_forces,
                                   plan->d_force_i0, plan->d_force_i1, plan->d_force_j0, plan->d_force_j1, plan->d_force_k0, plan->d_force_k1,
                                   plan->nsteps, plan->step,
                                   plan->d_mat, plan->d_mattype, plan->d_mattab, plan->matmode, plan->nmat);
        cudaThreadSynchronize();
         
        return cudaGetLastError();
}

cudaError_t wrap_srchKernel(SvGpuPlan *plan, int i, int j, int k, float hx, float hy, float hz)
{
        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

        srchpointKernel<<<dimGrid, dimBlock>>>(plan->d_hx, plan->d_hy, plan->d_hz, plan->xres, plan->yres, plan->zres, i, j, k, hx, hy, hz, plan->dir);
        cudaThreadSynchronize();
         
        return cudaGetLastError();
}

cudaError_t wrap_ffKernel(SvGpuPlan *plan)
{
     //   int k;
     //   int res = MAX(plan->yres, plan->h_iset[36]);  //use maximum of yres OR number of far field points, assuming that this is never more than n of threads
     //   dim3 dimBlockx(res);//was yres
     //   dim3 dimGridx(1);

        //dim3 dimBlockx(plan->maxthreads);
        //dim3 dimGridx((int)ceil((float)(plan->yres*plan->zres)/plan->maxthreads));

        //dim3 dimBlocky(plan->maxthreads);
        //dim3 dimGridy((int)ceil((float)(plan->xres*plan->zres)/plan->maxthreads));

       // dim3 dimBlockz(plan->maxthreads);
        //dim3 dimGridz((int)ceil((float)(plan->xres*plan->yres)/plan->maxthreads));

/*
	for (k = 0; k<plan->zres; k++) {
		pffKernelx<<<dimGridx, dimBlockx, 6*plan->yres*sizeof(float)>>>(plan->d_ex, plan->d_ey, plan->d_ez,
				plan->d_iset, 
				plan->d_ff_ex, plan->d_ff_ey, plan->d_ff_ez,
				plan->xres, plan->yres, plan->zres,
				plan->dx, plan->dy, plan->dz, plan->dt, 
				plan->step, k, res);

		cudaThreadSynchronize();
	}
*/
/*
        pffKernely<<<dimGridy, dimBlocky>>>(plan->d_ex, plan->d_ey, plan->d_ez,
                                        plan->d_iset, 
                                        plan->d_ff_ex, plan->d_ff_ey, plan->d_ff_ez,
                                        plan->xres, plan->yres, plan->zres,
                                        plan->dx, plan->dy, plan->dz, plan->dt, 
                                        plan->step, plan->maxthreads);

        cudaThreadSynchronize();

        pffKernelz<<<dimGridz, dimBlockz>>>(plan->d_ex, plan->d_ey, plan->d_ez,
                                        plan->d_iset, 
                                        plan->d_ff_ex, plan->d_ff_ey, plan->d_ff_ez,
                                        plan->xres, plan->yres, plan->zres,
                                        plan->dx, plan->dy, plan->dz, plan->dt, 
                                        plan->step, plan->maxthreads);

        cudaThreadSynchronize();
*/




        
        dim3 dimBlock(plan->maxthreads);
        dim3 dimGrid((int)ceil((float)plan->h_iset[NPNTS]/plan->maxthreads));
        ffKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez,
                                        plan->d_iset, 
                                        plan->d_ff_ex, plan->d_ff_ey, plan->d_ff_ez,
                                        plan->xres, plan->yres, plan->zres,
                                        plan->dx, plan->dy, plan->dz, plan->dt, 
                                        plan->step, plan->maxthreads);
        

        cudaThreadSynchronize();        
        return cudaGetLastError();
}


cudaError_t wrap_ffKernel_hlps(SvGpuPlan *plan)
{
        dim3 dimBlock(plan->h_iset[NPNTS]*plan->nhlps); //this should be 512 or 1024
        ffKernel_a<<<1, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez,
                                        plan->d_iset, 
                                        plan->d_ff_hlp_ex, plan->d_ff_hlp_ey, plan->d_ff_hlp_ez,
                                        plan->xres, plan->yres, plan->zres,
                                        plan->dx, plan->dy, plan->dz, plan->dt, 
                                        plan->step, plan->nhlps);

        cudaThreadSynchronize();    
        return cudaGetLastError();
}

cudaError_t wrap_ffKernel_gethlps(SvGpuPlan *plan)
{
        ffKernel_b<<<1, 1>>>(plan->d_ff_ex, plan->d_ff_ey, plan->d_ff_ez,
                             plan->d_ff_hlp_ex, plan->d_ff_hlp_ey, plan->d_ff_hlp_ez,
                             plan->d_iset, plan->nhlps);
        cudaThreadSynchronize();    
        return cudaGetLastError();
}


cudaError_t wrap_pffKernel(SvGpuPlan *plan)
{
        //h_piset[10] je nrs


   //TODO: split this to dimBlock(PNDIV) and dimGrid(rest). 
    
/*        dim3 dimGrid(PNDIV, (int)ceil((float)plan->h_piset[10]/(float)PNDIV));

        pffKernel<<<dimGrid, 1>>>(plan->d_ex, plan->d_ey, plan->d_ez,
                                        plan->d_piset, 
                                        plan->d_pff_ex, plan->d_pff_ey, plan->d_pff_ez,
                                        plan->xres, plan->yres, plan->zres,
                                        plan->dx, plan->dy, plan->dz, plan->dt, 
                                        plan->step);*/
       /*
        dim3 dimGrid(PNDIV, (int)ceil((float)plan->h_piset[10]/(float)PNDIV));
        int memsize = 2003*(plan->h_piset[PIMAX]-plan->h_piset[PIMIN])*(plan->h_piset[PJMAX]-plan->h_piset[PJMIN])*sizeof(float);

        pffsKernel<<<dimGrid, 1, 103*9*sizeof(float)>>>
                                       (plan->d_ex, plan->d_ey, plan->d_ez,
                                        plan->d_piset, 
                                        plan->d_pff_ex, plan->d_pff_ey, plan->d_pff_ez,
                                        plan->xres, plan->yres, plan->zres,
                                        plan->dx, plan->dy, plan->dz, plan->dt, 
                                        plan->step);

*/

        dim3 dimBlock(PNDIV);
        dim3 dimGrid((int)ceil((float)plan->h_piset[10]/(float)PNDIV));

        printf("pff kernel:  %d %d\n", (int)ceil((float)plan->h_piset[10]/(float)PNDIV), PNDIV);

        pffKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez,
                                        plan->d_piset, 
                                        plan->d_pff_ex, plan->d_pff_ey, plan->d_pff_ez,
                                        plan->xres, plan->yres, plan->zres,
                                        plan->dx, plan->dy, plan->dz, plan->dt, 
                                        plan->step);


        cudaThreadSynchronize();        

        return cudaGetLastError();
}

cudaError_t wrap_pffKernel_hlps(SvGpuPlan *plan)
{
        dim3 dimGrid(PNDIV, (int)ceil((float)plan->h_piset[10]/(float)PNDIV));
        dim3 dimBlock((plan->h_piset[PIMAX]-plan->h_piset[PIMIN])*(plan->h_piset[PJMAX]-plan->h_piset[PJMIN]));


        //printf("pff wrapper: %d %d %d\n", PNDIV, (int)ceil((float)plan->h_piset[PNPNTS]/(float)PNDIV), (plan->h_piset[PIMAX]-plan->h_piset[PIMIN])*(plan->h_piset[PJMAX]-plan->h_piset[PJMIN]));
        pffKernel_a<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez,
                                        plan->d_piset, 
                                        plan->d_pff_hlp_ex, plan->d_pff_hlp_ey, plan->d_pff_hlp_ez,
                                        plan->xres, plan->yres, plan->zres,
                                        plan->dx, plan->dy, plan->dz, plan->dt, 
                                        plan->step);

        cudaThreadSynchronize();    
        return cudaGetLastError();
}

cudaError_t wrap_pffKernel_gethlps(SvGpuPlan *plan)
{
        pffKernel_b<<<1, 1>>>(plan->d_pff_ex, plan->d_pff_ey, plan->d_pff_ez,
                             plan->d_pff_hlp_ex, plan->d_pff_hlp_ey, plan->d_pff_hlp_ez,
                             plan->d_piset);
        cudaThreadSynchronize();    
        return cudaGetLastError();
}


cudaError_t wrap_outpointKernel(SvGpuPlan *plan)
{
        outpointKernel<<<1, 1>>>(plan->d_ex, plan->d_ey, plan->d_ez, plan->d_hx, plan->d_hy, plan->d_hz,
                              plan->d_outpointdata, plan->d_outpoint_pos, plan->nsteps, plan->noutpoints,
                              plan->step, plan->xres, plan->yres, plan->zres);              

        cudaThreadSynchronize();
        return cudaGetLastError();


}

cudaError_t wrap_tsf_jstepKernel(SvGpuPlan *plan, float e)
{
        tsfjstepKernel<<<1, 1>>>(plan->d_tsf_jpool_e, plan->d_tsf_jpool_h, plan->d_tsf_jpvals, 
                              plan->d_tsf_jpool_epsilon, plan->d_tsf_jpool_mu,
                              plan->d_tsf_jpool_sigma, plan->d_tsf_jpool_sigast,
                              plan->dx/plan->h_tsfset[9], plan->dt, e, plan->tsf_jpool_size);

        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_tsff_jstepKernel(SvGpuPlan *plan, float e)
{
        int i;

	for (i=0; i<plan->tsff_nxm; i++) {
		tsffjstepKernel<<<1, 1>>>(plan->d_tsff_jpool_e, plan->d_tsff_jpool_h, plan->d_tsff_jpvals, 
				plan->d_tsff_jpool_epsilon, plan->d_tsff_jpool_mu,
				plan->d_tsff_jpool_sigma, plan->d_tsff_jpool_sigast,
				plan->dx/plan->tsff_data[6*i+3], plan->dt, e, plan->tsff_jpool_size,
				plan->tsff_nxm, i);

		cudaThreadSynchronize();
	}
	return cudaGetLastError();
}

cudaError_t wrap_sf_jstepKernel(SvGpuPlan *plan, float e)
{
        sfjstepKernel<<<1, 1>>>(plan->d_sf_jpool_e, plan->d_sf_jpool_h, plan->d_sf_jpvals, 
                              plan->d_sf_jpool_epsilon, plan->d_sf_jpool_mu,
                              plan->d_sf_jpool_sigma, plan->d_sf_jpool_sigast,
                              plan->dx, plan->dt, e, plan->sf_jpool_size);

        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_tsf_e_Kernel(SvGpuPlan *plan)
{
        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);


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
        tsf_estep_cKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, 
                                         plan->d_hx, plan->d_hy, plan->d_hz,
                                         plan->d_epsilon, plan->d_mu, plan->matmode,
                                         plan->d_tsfset,
                                         plan->d_tsf_jpool_e, plan->d_tsf_jpool_h, plan->tsf_jpool_size,
                                         plan->dx, plan->dt,
                                         plan->xres, plan->yres, plan->zres, plan->dir);


        cudaThreadSynchronize();
        return cudaGetLastError();
}

cudaError_t wrap_tsf_h_Kernel(SvGpuPlan *plan)
{
        int ga, gb, b;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

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
        tsf_hstep_cKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, 
                                         plan->d_hx, plan->d_hy, plan->d_hz,
                                         plan->d_epsilon, plan->d_mu, plan->matmode,
                                         plan->d_tsfset,
                                         plan->d_tsf_jpool_e, plan->d_tsf_jpool_h, plan->tsf_jpool_size,
                                         plan->dx, plan->dt,
                                         plan->xres, plan->yres, plan->zres, plan->dir);


        cudaThreadSynchronize();
        return cudaGetLastError();
}


cudaError_t wrap_tsff_e_Kernel(SvGpuPlan *plan)
{
        int ga, gb, b, i;
        float amplitude;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

        for (i=0; i<plan->tsff_nxm; i++) {

                amplitude = (float)(1e15*plan->tsff_data[6*i + 4]*plan->tsff_data[6*i + 5]*(plan->tsff_freal/(2*CU_PI*LIGHT_SPEED))*sqrt(cos(plan->tsff_data[6*i]))*sin(plan->tsff_data[6*i]));

		tsff_estep_aKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, 
				plan->d_hx, plan->d_hy, plan->d_hz,
				plan->d_epsilon, plan->d_mu, plan->matmode,
				plan->d_tsffset,
				plan->d_tsff_jpool_e, plan->d_tsff_jpool_h, plan->tsff_jpool_size,
				plan->dx, plan->dt,
				plan->xres, plan->yres, plan->zres, plan->dir,
				plan->tsff_data[6*i], plan->tsff_data[6*i + 1],
				plan->tsff_data[6*i + 2], plan->tsff_data[6*i + 3],
				amplitude, i, plan->tsff_nxm);
		tsff_estep_bKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, 
				plan->d_hx, plan->d_hy, plan->d_hz,
				plan->d_epsilon, plan->d_mu, plan->matmode,
				plan->d_tsffset,
				plan->d_tsff_jpool_e, plan->d_tsff_jpool_h, plan->tsff_jpool_size,
				plan->dx, plan->dt,
				plan->xres, plan->yres, plan->zres, plan->dir, 
				plan->tsff_data[6*i], plan->tsff_data[6*i + 1],
				plan->tsff_data[6*i + 2], plan->tsff_data[6*i + 3],
				amplitude, i, plan->tsff_nxm);
		tsff_estep_cKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, 
				plan->d_hx, plan->d_hy, plan->d_hz,
				plan->d_epsilon, plan->d_mu, plan->matmode,
				plan->d_tsffset,
				plan->d_tsff_jpool_e, plan->d_tsff_jpool_h, plan->tsff_jpool_size,
				plan->dx, plan->dt,
				plan->xres, plan->yres, plan->zres, plan->dir,
				plan->tsff_data[6*i], plan->tsff_data[6*i + 1],
				plan->tsff_data[6*i + 2], plan->tsff_data[6*i + 3],
				amplitude, i, plan->tsff_nxm);

		cudaThreadSynchronize();
	}

	return cudaGetLastError();
}

cudaError_t wrap_tsff_h_Kernel(SvGpuPlan *plan)
{
        int ga, gb, b, i;
        float amplitude;
        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

	for (i=0; i<plan->tsff_nxm; i++) {

                amplitude = (float)(1e15*plan->tsff_data[6*i + 4]*plan->tsff_data[6*i + 5]*(plan->tsff_freal/(2*CU_PI*LIGHT_SPEED))*sqrt(cos(plan->tsff_data[6*i]))*sin(plan->tsff_data[6*i]));

		tsff_hstep_aKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, 
				plan->d_hx, plan->d_hy, plan->d_hz,
				plan->d_epsilon, plan->d_mu, plan->matmode,
				plan->d_tsffset,
				plan->d_tsff_jpool_e, plan->d_tsff_jpool_h, plan->tsff_jpool_size,
				plan->dx, plan->dt,
				plan->xres, plan->yres, plan->zres, plan->dir,
				plan->tsff_data[6*i], plan->tsff_data[6*i + 1],
				plan->tsff_data[6*i + 2], plan->tsff_data[6*i + 3],
				amplitude, i, plan->tsff_nxm);

		tsff_hstep_bKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, 
				plan->d_hx, plan->d_hy, plan->d_hz,
				plan->d_epsilon, plan->d_mu, plan->matmode,
				plan->d_tsffset,
				plan->d_tsff_jpool_e, plan->d_tsff_jpool_h, plan->tsff_jpool_size,
				plan->dx, plan->dt,
				plan->xres, plan->yres, plan->zres, plan->dir,
				plan->tsff_data[6*i], plan->tsff_data[6*i + 1],
				plan->tsff_data[6*i + 2], plan->tsff_data[6*i + 3],
				amplitude, i, plan->tsff_nxm);

		tsff_hstep_cKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, 
				plan->d_hx, plan->d_hy, plan->d_hz,
				plan->d_epsilon, plan->d_mu, plan->matmode,
				plan->d_tsffset,
				plan->d_tsff_jpool_e, plan->d_tsff_jpool_h, plan->tsff_jpool_size,
				plan->dx, plan->dt,
				plan->xres, plan->yres, plan->zres, plan->dir,
				plan->tsff_data[6*i], plan->tsff_data[6*i + 1],
				plan->tsff_data[6*i + 2], plan->tsff_data[6*i + 3],
				amplitude, i, plan->tsff_nxm);


		cudaThreadSynchronize();
	}
	return cudaGetLastError();
}

/*
cudaError_t wrap_fiberKernel(TGPUplan *plan)
{
        int ga, gb, b;

        if (plan->dir==2) {
           b = plan->zres;
           ga = plan->xres;
           gb = plan->yres;
        } else if (plan->dir==0) {
           b = plan->xres;
           ga = plan->yres;
           gb = plan->zres;
        } else {
           b = plan->yres;
           ga = plan->xres;
           gb = plan->zres;
        } 

        dim3 dimBlock(b);
        dim3 dimGrid(ga, gb);

        fiberKernel<<<dimGrid, dimBlock>>>(plan->d_ex, plan->d_ey, plan->d_ez, 
                                         plan->d_hx, plan->d_hy, plan->d_hz,
                                         plan->d_fiberset,
                                         plan->d_jpool_e, plan->d_jpool_h, plan->jpool_size,
                                         plan->dx, plan->dt,
                                         plan->xres, plan->yres, plan->zres, plan->dir);

        cudaThreadSynchronize();
        return cudaGetLastError();
}
*/
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
