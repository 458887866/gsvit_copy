
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

/*  boundary.h : 
 *  boundary conditions implementation
 *
 */

#include "boundary.h"
#include "pool.h"
#include "settings.h"
#include "constants.h"
#include <math.h>
#include <omp.h>

static void
becebhch(gint i, gint imax, gdouble sigma, gdouble kappa, gdouble a, gint m, gdouble dt,
	 gdouble *be, gdouble *ce, gdouble *bh, gdouble *ch, gdouble *vkappa, gdouble *hkappa)
{
    gdouble ma = 1;
    gdouble powdepthe = pow((imax - i)/(gdouble)(imax - 0), m);
    gdouble powdepthh = pow((imax - i + 0.5)/(gdouble)(imax - 0), m);

    gdouble sve = sigma*powdepthe;
    gdouble kve = 1 + (kappa - 1)*powdepthe;
    gdouble ave = a*pow(((i-0)/((gdouble)imax-0)),ma);
    
    gdouble svh = sigma*powdepthh;
    gdouble kvh = 1 + (kappa - 1)*powdepthh;
    gdouble avh = a*pow(((i+0.5)/((gdouble)imax-0)),ma);

    *be = exp(-(sve/kve + ave)*dt/EPSILON_0);

    if (sve == 0)
	*ce=0;
    else
	*ce = (sve/(sve + kve*ave)/kve)*((*be) - 1.0);

    *bh = exp(-(svh/kvh + avh)*dt/EPSILON_0);

    if (svh == 0)
	*ch = 0;
    else
	*ch = (svh/(svh + kvh*avh)/kvh)*((*bh) - 1.0);

    *vkappa = kve;
    *hkappa = kvh;
}

static void
becebhch0(gint i, gint imax, gdouble sigma, gdouble kappa, gdouble a, gint m, gdouble dt,
	  gdouble *be, gdouble *ce, gdouble *bh, gdouble *ch, gdouble *vkappa, gdouble *hkappa)
{
    gdouble ma = 1;
    gdouble powdepthe = pow((imax - i)/(gdouble)(imax - 0), m);
    gdouble powdepthh = pow((imax - i - 0.5)/(gdouble)(imax - 0), m);

    gdouble sve = sigma*powdepthe;
    gdouble kve = 1 + (kappa - 1)*powdepthe;
    gdouble ave = a*pow(((i-0)/((gdouble)imax-0)),ma);
    
    gdouble svh = sigma*powdepthh;
    gdouble kvh = 1 + (kappa - 1)*powdepthh;
    gdouble avh = a*pow(((i-0.5)/((gdouble)imax-0)),ma);

    *be = exp(-(sve/kve + ave)*dt/EPSILON_0);
    
    if (sve == 0)
	*ce=0;
    else
	*ce = (sve/(sve + kve*ave)/kve)*((*be) - 1.0);

    *bh = exp(-(svh/kvh + avh)*dt/EPSILON_0);
    
    if (svh == 0)
	*ch = 0;
    else
	*ch = (svh/(svh + kvh*avh)/kvh)*((*bh) - 1.0);

    *vkappa = kve;
    *hkappa = kvh;
}

SvBoundary*
sv_boundary_new(SvPool *mp, SvSet *set)
{
    gint i;
    SvBoundary* bnd = (SvBoundary*)g_malloc(sizeof(SvBoundary));

    if (set->sb.bx0 == SV_BOUNDARY_LIAO || set->sb.bxn == SV_BOUNDARY_LIAO) {
        bnd->liao.fbnx = sv_dcube_new(NPLANES, set->sp.yres, set->sp.zres, NPLANES, set->sp.yres, set->sp.zres, 1);
        if (set->sc.verbose>1) 
            printf("X Liao boundary initialized\n");
    }

    if (set->sb.by0 == SV_BOUNDARY_LIAO || set->sb.byn == SV_BOUNDARY_LIAO) {
        bnd->liao.fbny = sv_dcube_new(NPLANES, set->sp.xres, set->sp.zres, NPLANES, set->sp.yres, set->sp.zres, 1);
        if (set->sc.verbose>1) 
            printf("Y Liao boundary initialized\n");
    }

    if (set->sb.bz0 == SV_BOUNDARY_LIAO || set->sb.bzn == SV_BOUNDARY_LIAO) {
        bnd->liao.fbnz = sv_dcube_new(NPLANES, set->sp.xres, set->sp.yres, NPLANES, set->sp.yres, set->sp.zres, 1);
        if (set->sc.verbose>1) 
            printf("Z Liao boundary initialized\n");
    }

    if (set->sb.bx0 == SV_BOUNDARY_CPML) {

        bnd->cpml.peyx_x0 = sv_dcube_new(set->sb.depth_bx0, set->sp.yres, set->sp.zres, set->sb.depth_bx0, set->sp.yres, set->sp.zres, 1);
        bnd->cpml.pezx_x0 = sv_dcube_new(set->sb.depth_bx0, set->sp.yres, set->sp.zres, set->sb.depth_bx0, set->sp.yres, set->sp.zres, 1);

        bnd->cpml.phzx_x0 = sv_dcube_new(set->sb.depth_bx0, set->sp.yres, set->sp.zres, set->sb.depth_bx0, set->sp.yres, set->sp.zres, 1);
        bnd->cpml.phyx_x0 = sv_dcube_new(set->sb.depth_bx0, set->sp.yres, set->sp.zres, set->sb.depth_bx0, set->sp.yres, set->sp.zres, 1);

        bnd->cpml.be_x0 = (gdouble *) g_malloc(set->sb.depth_bx0*sizeof(gdouble));
        bnd->cpml.ce_x0 = (gdouble *) g_malloc(set->sb.depth_bx0*sizeof(gdouble));
        bnd->cpml.bh_x0 = (gdouble *) g_malloc(set->sb.depth_bx0*sizeof(gdouble));
        bnd->cpml.ch_x0 = (gdouble *) g_malloc(set->sb.depth_bx0*sizeof(gdouble)); 
        bnd->cpml.kappae_x0 = (gdouble *) g_malloc(set->sb.depth_bx0*sizeof(gdouble));
        bnd->cpml.kappah_x0 = (gdouble *) g_malloc(set->sb.depth_bx0*sizeof(gdouble));

        for (i = 0; i < set->sb.depth_bx0; i++) {
            becebhch0(i, set->sb.depth_bx0, set->sb.sigma_bx0, set->sb.kappa_bx0, set->sb.a_bx0, set->sb.m_bx0, set->plan.dt,
		              bnd->cpml.be_x0+i, bnd->cpml.ce_x0+i, bnd->cpml.bh_x0+i, bnd->cpml.ch_x0+i, bnd->cpml.kappae_x0+i, bnd->cpml.kappah_x0+i);
//	        printf("x0  %d be = %g ce = %g bh = %g ch = %g kappae = %g kappah = %g\n", i, bnd->cpml.be_x0[i], bnd->cpml.ce_x0[i], bnd->cpml.bh_x0[i], bnd->cpml.ch_x0[i], bnd->cpml.kappae_x0[i], bnd->cpml.kappah_x0[i]);
        }

        if (set->sc.verbose > 1) 
            printf("X0 CPML boundary initialized\n");
    }

    if (set->sb.bxn == SV_BOUNDARY_CPML) {

        bnd->cpml.peyx_xn = sv_dcube_new(set->sb.depth_bxn, set->sp.yres, set->sp.zres, set->sb.depth_bxn, set->sp.yres, set->sp.zres, 1);
        bnd->cpml.pezx_xn = sv_dcube_new(set->sb.depth_bxn, set->sp.yres, set->sp.zres, set->sb.depth_bxn, set->sp.yres, set->sp.zres, 1);

        bnd->cpml.phzx_xn = sv_dcube_new(set->sb.depth_bxn, set->sp.yres, set->sp.zres, set->sb.depth_bxn, set->sp.yres, set->sp.zres, 1);
        bnd->cpml.phyx_xn = sv_dcube_new(set->sb.depth_bxn, set->sp.yres, set->sp.zres, set->sb.depth_bxn, set->sp.yres, set->sp.zres, 1);

        bnd->cpml.be_xn = (gdouble *) g_malloc(set->sb.depth_bxn*sizeof(gdouble));
        bnd->cpml.ce_xn = (gdouble *) g_malloc(set->sb.depth_bxn*sizeof(gdouble));
        bnd->cpml.bh_xn = (gdouble *) g_malloc(set->sb.depth_bxn*sizeof(gdouble));
        bnd->cpml.ch_xn = (gdouble *) g_malloc(set->sb.depth_bxn*sizeof(gdouble)); 
        bnd->cpml.kappae_xn = (gdouble *) g_malloc(set->sb.depth_bxn*sizeof(gdouble));
        bnd->cpml.kappah_xn = (gdouble *) g_malloc(set->sb.depth_bxn*sizeof(gdouble));

        for (i = 0; i < set->sb.depth_bxn; i++) {
            becebhch(set->sb.depth_bxn - i - 0, set->sb.depth_bxn, set->sb.sigma_bxn, set->sb.kappa_bxn, set->sb.a_bxn, set->sb.m_bxn, set->plan.dt,
		             bnd->cpml.be_xn+i, bnd->cpml.ce_xn+i, bnd->cpml.bh_xn+i, bnd->cpml.ch_xn+i, bnd->cpml.kappae_xn+i, bnd->cpml.kappah_xn+i);

//	            printf("xn  %d be = %g ce = %g bh = %g ch = %g kappae = %g kappah = %g\n", i, bnd->cpml.be_xn[i], bnd->cpml.ce_xn[i], bnd->cpml.bh_xn[i], bnd->cpml.ch_xn[i], bnd->cpml.kappae_xn[i], bnd->cpml.kappah_xn[i]);
        }

        if (set->sc.verbose > 1) 
            printf("XN CPML boundary initialized\n");
    }

    if (set->sb.by0 == SV_BOUNDARY_CPML) {
        bnd->cpml.pexy_y0 = sv_dcube_new(set->sp.xres, set->sb.depth_by0, set->sp.zres, set->sp.xres, set->sb.depth_by0, set->sp.zres, 1);
        bnd->cpml.pezy_y0 = sv_dcube_new(set->sp.xres, set->sb.depth_by0, set->sp.zres, set->sp.xres, set->sb.depth_by0, set->sp.zres, 1);

        bnd->cpml.phxy_y0 = sv_dcube_new(set->sp.xres, set->sb.depth_by0, set->sp.zres, set->sp.xres, set->sb.depth_by0, set->sp.zres, 1);
        bnd->cpml.phzy_y0 = sv_dcube_new(set->sp.xres, set->sb.depth_by0, set->sp.zres, set->sp.xres, set->sb.depth_by0, set->sp.zres, 1);
   
        bnd->cpml.be_y0 = (gdouble *) g_malloc(set->sb.depth_by0*sizeof(gdouble));
        bnd->cpml.ce_y0 = (gdouble *) g_malloc(set->sb.depth_by0*sizeof(gdouble));
        bnd->cpml.bh_y0 = (gdouble *) g_malloc(set->sb.depth_by0*sizeof(gdouble));
        bnd->cpml.ch_y0 = (gdouble *) g_malloc(set->sb.depth_by0*sizeof(gdouble)); 
        bnd->cpml.kappae_y0 = (gdouble *) g_malloc(set->sb.depth_by0*sizeof(gdouble));
        bnd->cpml.kappah_y0 = (gdouble *) g_malloc(set->sb.depth_by0*sizeof(gdouble));

        for (i = 0; i < set->sb.depth_by0; i++) {
            becebhch0(i, set->sb.depth_by0, set->sb.sigma_by0, set->sb.kappa_by0, set->sb.a_by0, set->sb.m_by0, set->plan.dt,
		              bnd->cpml.be_y0+i, bnd->cpml.ce_y0+i, bnd->cpml.bh_y0+i, bnd->cpml.ch_y0+i, bnd->cpml.kappae_y0+i, bnd->cpml.kappah_y0+i);

  //          printf("y0 %d be = %g ce = %g bh = %g ch = %g kappae = %g kappah = %g\n", i, bnd->cpml.be_y0[i], bnd->cpml.ce_y0[i], 
  //                 bnd->cpml.bh_y0[i], bnd->cpml.ch_y0[i], bnd->cpml.kappae_y0[i], bnd->cpml.kappah_y0[i]);
        }

        if (set->sc.verbose > 1) 
            printf("Y0 CPML boundary initialized\n");
    }

    if (set->sb.byn == SV_BOUNDARY_CPML) {
        bnd->cpml.pexy_yn = sv_dcube_new(set->sp.xres, set->sb.depth_byn, set->sp.zres, set->sp.xres, set->sb.depth_byn, set->sp.zres, 1);
        bnd->cpml.pezy_yn = sv_dcube_new(set->sp.xres, set->sb.depth_byn, set->sp.zres, set->sp.xres, set->sb.depth_byn, set->sp.zres, 1);

        bnd->cpml.phxy_yn = sv_dcube_new(set->sp.xres, set->sb.depth_byn, set->sp.zres, set->sp.xres, set->sb.depth_byn, set->sp.zres, 1);
        bnd->cpml.phzy_yn = sv_dcube_new(set->sp.xres, set->sb.depth_byn, set->sp.zres, set->sp.xres, set->sb.depth_byn, set->sp.zres, 1);

        bnd->cpml.be_yn = (gdouble *) g_malloc(set->sb.depth_byn*sizeof(gdouble));
        bnd->cpml.ce_yn = (gdouble *) g_malloc(set->sb.depth_byn*sizeof(gdouble));
        bnd->cpml.bh_yn = (gdouble *) g_malloc(set->sb.depth_byn*sizeof(gdouble));
        bnd->cpml.ch_yn = (gdouble *) g_malloc(set->sb.depth_byn*sizeof(gdouble)); 
        bnd->cpml.kappae_yn = (gdouble *) g_malloc(set->sb.depth_byn*sizeof(gdouble));
        bnd->cpml.kappah_yn = (gdouble *) g_malloc(set->sb.depth_byn*sizeof(gdouble));

        for (i = 0; i < set->sb.depth_byn; i++) {
            becebhch(set->sb.depth_byn - i - 0, set->sb.depth_byn, set->sb.sigma_byn, set->sb.kappa_byn, set->sb.a_byn, set->sb.m_byn, set->plan.dt,
		             bnd->cpml.be_yn+i, bnd->cpml.ce_yn+i, bnd->cpml.bh_yn+i, bnd->cpml.ch_yn+i, bnd->cpml.kappae_yn+i, bnd->cpml.kappah_yn+i);
  //          printf("yn %d be = %g ce = %g bh = %g ch = %g kappae = %g kappah = %g\n", i, bnd->cpml.be_yn[i], bnd->cpml.ce_yn[i], 
  //          bnd->cpml.bh_yn[i], bnd->cpml.ch_yn[i], bnd->cpml.kappae_yn[i], bnd->cpml.kappah_yn[i]);
        }

        if (set->sc.verbose > 1) 
            printf("YN CPML boundary initialized\n");
    }

    if (set->sb.bz0 == SV_BOUNDARY_CPML) {
        bnd->cpml.pexz_z0 = sv_dcube_new(set->sp.xres, set->sp.yres, set->sb.depth_bz0, set->sp.xres, set->sp.yres, set->sb.depth_bz0, 1);
        bnd->cpml.peyz_z0 = sv_dcube_new(set->sp.xres, set->sp.yres, set->sb.depth_bz0, set->sp.xres, set->sp.yres, set->sb.depth_bz0, 1);

        bnd->cpml.phyz_z0 = sv_dcube_new(set->sp.xres, set->sp.yres, set->sb.depth_bz0, set->sp.xres, set->sp.yres, set->sb.depth_bz0, 1);
        bnd->cpml.phxz_z0 = sv_dcube_new(set->sp.xres, set->sp.yres, set->sb.depth_bz0, set->sp.xres, set->sp.yres, set->sb.depth_bz0, 1);

        bnd->cpml.be_z0 = (gdouble *) g_malloc(set->sb.depth_bz0*sizeof(gdouble));
        bnd->cpml.ce_z0 = (gdouble *) g_malloc(set->sb.depth_bz0*sizeof(gdouble));
        bnd->cpml.bh_z0 = (gdouble *) g_malloc(set->sb.depth_bz0*sizeof(gdouble));
        bnd->cpml.ch_z0 = (gdouble *) g_malloc(set->sb.depth_bz0*sizeof(gdouble)); 
        bnd->cpml.kappae_z0 = (gdouble *) g_malloc(set->sb.depth_bz0*sizeof(gdouble));
        bnd->cpml.kappah_z0 = (gdouble *) g_malloc(set->sb.depth_bz0*sizeof(gdouble));

        for (i = 0; i < set->sb.depth_bz0; i++) {
            becebhch0(i, set->sb.depth_bz0, set->sb.sigma_bz0, set->sb.kappa_bz0, set->sb.a_bz0, set->sb.m_bz0, set->plan.dt,
		              bnd->cpml.be_z0+i, bnd->cpml.ce_z0+i, bnd->cpml.bh_z0+i, bnd->cpml.ch_z0+i, bnd->cpml.kappae_z0+i, bnd->cpml.kappah_z0+i);
     
 //           printf("z0 %d be = %g ce = %g bh = %g ch = %g kappae = %g kappah = %g\n", i, bnd->cpml.be_z0[i], bnd->cpml.ce_z0[i], 
 //                  bnd->cpml.bh_z0[i], bnd->cpml.ch_z0[i], bnd->cpml.kappae_z0[i], bnd->cpml.kappah_z0[i]);
   }

        if (set->sc.verbose > 1) 
            printf("Z0 CPML boundary initialized\n");
    }

    if (set->sb.bzn == SV_BOUNDARY_CPML) {
        bnd->cpml.pexz_zn = sv_dcube_new(set->sp.xres, set->sp.yres, set->sb.depth_bzn, set->sp.xres, set->sp.yres, set->sb.depth_bzn, 1);
        bnd->cpml.peyz_zn = sv_dcube_new(set->sp.xres, set->sp.yres, set->sb.depth_bzn, set->sp.xres, set->sp.yres, set->sb.depth_bzn, 1);

        bnd->cpml.phyz_zn = sv_dcube_new(set->sp.xres, set->sp.yres, set->sb.depth_bzn, set->sp.xres, set->sp.yres, set->sb.depth_bzn, 1);
        bnd->cpml.phxz_zn = sv_dcube_new(set->sp.xres, set->sp.yres, set->sb.depth_bzn, set->sp.xres, set->sp.yres, set->sb.depth_bzn, 1);

        bnd->cpml.be_zn = (gdouble *) g_malloc(set->sb.depth_bzn*sizeof(gdouble));
        bnd->cpml.ce_zn = (gdouble *) g_malloc(set->sb.depth_bzn*sizeof(gdouble));
        bnd->cpml.bh_zn = (gdouble *) g_malloc(set->sb.depth_bzn*sizeof(gdouble));
        bnd->cpml.ch_zn = (gdouble *) g_malloc(set->sb.depth_bzn*sizeof(gdouble)); 
        bnd->cpml.kappae_zn = (gdouble *) g_malloc(set->sb.depth_bzn*sizeof(gdouble));
        bnd->cpml.kappah_zn = (gdouble *) g_malloc(set->sb.depth_bzn*sizeof(gdouble));

        for (i = 0; i < set->sb.depth_bzn; i++) {
            becebhch(set->sb.depth_bzn - i - 0, set->sb.depth_bzn, set->sb.sigma_bzn, set->sb.kappa_bzn, set->sb.a_bzn, set->sb.m_bzn, set->plan.dt,
		             bnd->cpml.be_zn+i, bnd->cpml.ce_zn+i, bnd->cpml.bh_zn+i, bnd->cpml.ch_zn+i, bnd->cpml.kappae_zn+i, bnd->cpml.kappah_zn+i);
    
 //         printf("zn %d be = %g ce = %g bh = %g ch = %g kappae = %g kappah = %g\n", i, bnd->cpml.be_zn[i], bnd->cpml.ce_zn[i], 
 //         bnd->cpml.bh_zn[i], bnd->cpml.ch_zn[i], bnd->cpml.kappae_zn[i], bnd->cpml.kappah_zn[i]);
        }
         
        if (set->sc.verbose > 1) 
            printf("ZN CPML boundary initialized\n");
    }
    
    return bnd;
}

void
sv_pool_boundary_hstep(SvPool *mp, SvSet *set)
{
    gint i, j, k;
    gdouble db, bh, ch;
    gint xres = set->sp.xres;
    gint yres = set->sp.yres;
    gint zres = set->sp.zres;
    gdouble dx = set->sp.dx;
    gdouble dy = set->sp.dy;
    gdouble dz = set->sp.dz;
    /*gdouble dt = set->plan.dt;*/
    gint ir, jr, kr;

    db = bh = ch = 0;
    ir = jr = kr = 0;
  
    if (set->sc.verbose > 1) {
        printf("Running boundary H step...  \n ");
        fflush(stdout);
    }

    /*material boundary*/
    if (set->smb.bx0 == SV_BOUNDARY_PERIODIC || set->smb.bxn == SV_BOUNDARY_PERIODIC) {
	/*x planes*/
        for (j = set->smb.by0pos; j < set->smb.bynpos; j++) {
            for (k = set->smb.bz0pos; k < set->smb.bznpos; k++) {
		//      mp->d->hx->data[set->smb.bx0pos-1][j][k] = mp->d->hx->data[set->smb.bxnpos-1][j][k];
                mp->d->hy->data[set->smb.bx0pos-1][j][k] = mp->d->hy->data[set->smb.bxnpos-1][j][k];
                mp->d->hz->data[set->smb.bx0pos-1][j][k] = mp->d->hz->data[set->smb.bxnpos-1][j][k];

		//      mp->d->hx->data[set->smb.bxnpos][j][k] = mp->d->hx->data[set->smb.bx0pos][j][k];
		//      mp->d->hy->data[set->smb.bxnpos][j][k] = mp->d->hy->data[set->smb.bx0pos][j][k];
		//      mp->d->hz->data[set->smb.bxnpos][j][k] = mp->d->hz->data[set->smb.bx0pos][j][k];
	        }
        }
    }
    if (set->smb.by0 == SV_BOUNDARY_PERIODIC || set->smb.byn == SV_BOUNDARY_PERIODIC) {
        /*y planes*/
        for (i = set->smb.bx0pos; i < set->smb.bxnpos; i++) {
            for (k = set->smb.bz0pos; k < set->smb.bznpos; k++) {
                mp->d->hx->data[i][set->smb.by0pos-1][k] = mp->d->hx->data[i][set->smb.bynpos-1][k];
		//      mp->d->hy->data[i][set->smb.by0pos-1][k] = mp->d->hy->data[i][set->smb.bynpos-1][k];
                mp->d->hz->data[i][set->smb.by0pos-1][k] = mp->d->hz->data[i][set->smb.bynpos-1][k];

		//      mp->d->hx->data[i][set->smb.bynpos][k] = mp->d->hx->data[i][set->smb.by0pos][k];
		//      mp->d->hy->data[i][set->smb.bynpos][k] = mp->d->hy->data[i][set->smb.by0pos][k];
		//      mp->d->hz->data[i][set->smb.bynpos][k] = mp->d->hz->data[i][set->smb.by0pos][k];
	        }
        }
    }
    if (set->smb.bz0 == SV_BOUNDARY_PERIODIC || set->smb.bzn == SV_BOUNDARY_PERIODIC) {
        /*z planes*/
        for (i = set->smb.bx0pos; i < set->smb.bxnpos; i++) {
            for (j = set->smb.by0pos; j < set->smb.bynpos; j++) {
                mp->d->hx->data[i][j][set->smb.bz0pos-1] = mp->d->hx->data[i][j][set->smb.bznpos-1];
                mp->d->hy->data[i][j][set->smb.bz0pos-1] = mp->d->hy->data[i][j][set->smb.bznpos-1];
		//      mp->d->hz->data[i][j][set->smb.bz0pos-1] = mp->d->hz->data[i][j][set->smb.bznpos-1];

		//      mp->d->hx->data[i][j][set->smb.bznpos] = mp->d->hx->data[i][j][set->smb.bz0pos];
		//      mp->d->hy->data[i][j][set->smb.bznpos] = mp->d->hy->data[i][j][set->smb.bz0pos];
		//      mp->d->hz->data[i][j][set->smb.bznpos] = mp->d->hz->data[i][j][set->smb.bz0pos];
	        }
        }
    }
    if (set->sb.bx0 == SV_BOUNDARY_CPML) {
	    for (i = 0; i < (set->sb.depth_bx0 - 1); i++) {
	    
#pragma omp parallel default(shared) private(j, k, db, bh, ch)
#pragma omp for nowait
	    
            for (j = 1; j < (yres-1); j++) {
                for (k = 1; k < (zres-1); k++) {

                    db = sv_yee_data_get_db(mp->d, set, mp->mats, mp->nmat, mp->d->dt, i, j, k);
                    bh = mp->bnd->cpml.bh_x0[i];
                    ch = mp->bnd->cpml.ch_x0[i];

                    mp->bnd->cpml.phzx_x0->data[i][j][k] = bh*mp->bnd->cpml.phzx_x0->data[i][j][k] 
                               + ch*((mp->d->ey->data[i+1][j][k] - mp->d->ey->data[i][j][k]))/dx;
                    mp->bnd->cpml.phyx_x0->data[i][j][k] = bh*mp->bnd->cpml.phyx_x0->data[i][j][k] 
                               + ch*((mp->d->ez->data[i+1][j][k] - mp->d->ez->data[i][j][k]))/dx;

                    mp->d->hy->data[i][j][k] += dx*db*mp->bnd->cpml.phyx_x0->data[i][j][k]; //-
                    mp->d->hz->data[i][j][k] -= dx*db*mp->bnd->cpml.phzx_x0->data[i][j][k]; //+
                } // k
            } // j
        } // i
    } // if


    if (set->sb.by0 == SV_BOUNDARY_CPML) {
        for (j = 0; j < (set->sb.depth_by0 - 1); j++) {
	    
#pragma omp parallel default(shared) private(i, k, db, bh, ch)
#pragma omp for nowait
	    
            for (i = 1; i < (xres-1); i++) {
                for (k = 1; k < (zres-1); k++) {

                    db = sv_yee_data_get_db(mp->d, set, mp->mats, mp->nmat, mp->d->dt, i, j, k);
                    bh = mp->bnd->cpml.bh_y0[j];
                    ch = mp->bnd->cpml.ch_y0[j];

                    mp->bnd->cpml.phxy_y0->data[i][j][k] = bh*mp->bnd->cpml.phxy_y0->data[i][j][k] 
                               + ch*((mp->d->ez->data[i][j+1][k] - mp->d->ez->data[i][j][k]))/dy;
                    mp->bnd->cpml.phzy_y0->data[i][j][k] = bh*mp->bnd->cpml.phzy_y0->data[i][j][k] 
                               + ch*((mp->d->ex->data[i][j+1][k] - mp->d->ex->data[i][j][k]))/dy;

                    mp->d->hx->data[i][j][k] -= dy*db*mp->bnd->cpml.phxy_y0->data[i][j][k];
                    mp->d->hz->data[i][j][k] += dy*db*mp->bnd->cpml.phzy_y0->data[i][j][k]; 
                } // k
            } // i
        } // j
    } // if

    if (set->sb.bz0 == SV_BOUNDARY_CPML) {
        for (k = 0; k < (set->sb.depth_bz0 - 1); k++) {
	    
#pragma omp parallel default(shared) private(i, j, db, bh, ch)
#pragma omp for nowait

            for (i = 1; i < (xres-1); i++) {
                for (j = 1; j < (yres-1); j++) {

                    db = sv_yee_data_get_db(mp->d, set, mp->mats, mp->nmat, mp->d->dt, i, j, k);
                    bh = mp->bnd->cpml.bh_z0[k];
                    ch = mp->bnd->cpml.ch_z0[k];

                    mp->bnd->cpml.phyz_z0->data[i][j][k] = bh*mp->bnd->cpml.phyz_z0->data[i][j][k] 
                               + ch*((mp->d->ex->data[i][j][k+1] - mp->d->ex->data[i][j][k]))/dz;
                    mp->bnd->cpml.phxz_z0->data[i][j][k] = bh*mp->bnd->cpml.phxz_z0->data[i][j][k] 
                               + ch*((mp->d->ey->data[i][j][k+1] - mp->d->ey->data[i][j][k]))/dz;

                    mp->d->hy->data[i][j][k] -= dz*db*mp->bnd->cpml.phyz_z0->data[i][j][k]; //-
                    mp->d->hx->data[i][j][k] += dz*db*mp->bnd->cpml.phxz_z0->data[i][j][k]; //+
                } // j
            } // i
        } // k
    } // if

    if (set->sb.bxn == SV_BOUNDARY_CPML) {
  	    for (j = 1; j < (yres-1); j++) {
		    for (k = 1; k < (zres-1); k++) {
                         mp->d->hx->data[xres - 1][j][k] = mp->d->hy->data[xres-1][j][k] = mp->d->hz->data[xres-1][j][k] = 0;
		    }
	    }
      for (i = (xres - set->sb.depth_bxn + 1); i < (xres-1); i++) {
	    
#pragma omp parallel default(shared) private(j, k, db, bh, ch, ir)
#pragma omp for nowait
	    
            for (j = 1; j < (yres-1); j++) {
                for (k = 1; k < (zres-1); k++) {
                    ir = i - (xres - set->sb.depth_bxn);

                    db = sv_yee_data_get_db(mp->d, set, mp->mats, mp->nmat, mp->d->dt, i, j, k);
                    bh = mp->bnd->cpml.bh_xn[ir];
                    ch = mp->bnd->cpml.ch_xn[ir];

                    mp->bnd->cpml.phzx_xn->data[ir][j][k] = bh*mp->bnd->cpml.phzx_xn->data[ir][j][k] 
                                + ch*((mp->d->ey->data[i+1][j][k] - mp->d->ey->data[i][j][k]))/dx;
                    mp->bnd->cpml.phyx_xn->data[ir][j][k] = bh*mp->bnd->cpml.phyx_xn->data[ir][j][k] 
                                + ch*((mp->d->ez->data[i+1][j][k] - mp->d->ez->data[i][j][k]))/dx;

                    mp->d->hy->data[i][j][k] += dx*db*mp->bnd->cpml.phyx_xn->data[ir][j][k]; //-
                    mp->d->hz->data[i][j][k] -= dx*db*mp->bnd->cpml.phzx_xn->data[ir][j][k]; //+
                } // k
            } // j
        } // i
    } // if

    if (set->sb.byn == SV_BOUNDARY_CPML) {
            for (i = 1; i < (xres-1); i++) {
                for (k = 1; k < (zres-1); k++) {
                    mp->d->hx->data[i][yres-1][k] = mp->d->hy->data[i][yres-1][k] = mp->d->hz->data[i][yres-1][k] = 0;
                }
        }
        for (j = (yres - set->sb.depth_byn + 1); j < (yres-1); j++) {
	    
#pragma omp parallel default(shared) private(i,k, db, bh, ch, jr)
#pragma omp for nowait
	    
            for (i = 1; i < (xres-1); i++) {
                for (k = 1; k < (zres-1); k++) {
                    jr = j - (yres - set->sb.depth_byn);

                    db = sv_yee_data_get_db(mp->d, set, mp->mats, mp->nmat, mp->d->dt, i, j, k);
                    bh = mp->bnd->cpml.bh_yn[jr];
                    ch = mp->bnd->cpml.ch_yn[jr];

                    mp->bnd->cpml.phxy_yn->data[i][jr][k] = bh*mp->bnd->cpml.phxy_yn->data[i][jr][k] 
                                 + ch*((mp->d->ez->data[i][j+1][k] - mp->d->ez->data[i][j][k]))/dy;
                    mp->bnd->cpml.phzy_yn->data[i][jr][k] = bh*mp->bnd->cpml.phzy_yn->data[i][jr][k] 
                                 + ch*((mp->d->ex->data[i][j+1][k] - mp->d->ex->data[i][j][k]))/dy;

                    mp->d->hx->data[i][j][k] -= dy*db*mp->bnd->cpml.phxy_yn->data[i][jr][k];
                    mp->d->hz->data[i][j][k] += dy*db*mp->bnd->cpml.phzy_yn->data[i][jr][k]; 
                } // k
            } // i
        } // j
    } // if


    if (set->sb.bzn == SV_BOUNDARY_CPML) {
        for (i = 1; i < (xres-1); i++) {
            for (j = 1; j < (yres-1); j++) {
                mp->d->hx->data[i][j][zres-1] = mp->d->hy->data[i][j][zres-1] = mp->d->hz->data[i][j][zres-1] = 0;
            }
        }
        for (k = (zres - set->sb.depth_bzn + 1); k < (zres-1); k++) {

#pragma omp parallel default(shared) private(i, j, db, bh, ch, kr)
#pragma omp for nowait
	    
            for (i = 1; i < (xres-1); i++) {
                for (j = 1; j < (yres-1); j++) {
                    kr = k - (zres - set->sb.depth_bzn);

                    db = sv_yee_data_get_db(mp->d, set, mp->mats, mp->nmat, mp->d->dt, i, j, k);
                    bh = mp->bnd->cpml.bh_zn[kr];
                    ch = mp->bnd->cpml.ch_zn[kr];

                    mp->bnd->cpml.phyz_zn->data[i][j][kr] = bh*mp->bnd->cpml.phyz_zn->data[i][j][kr] 
                             + ch*((mp->d->ex->data[i][j][k+1] - mp->d->ex->data[i][j][k]))/dz;
                    mp->bnd->cpml.phxz_zn->data[i][j][kr] = bh*mp->bnd->cpml.phxz_zn->data[i][j][kr] 
                             + ch*((mp->d->ey->data[i][j][k+1] - mp->d->ey->data[i][j][k]))/dz;

                    mp->d->hy->data[i][j][k] -= dz*db*mp->bnd->cpml.phyz_zn->data[i][j][kr]; //-
                    mp->d->hx->data[i][j][k] += dz*db*mp->bnd->cpml.phxz_zn->data[i][j][kr]; //+
                } // j
            } // i
        } // k
    } // if

    if (set->sc.verbose > 1) 
        printf("done.\n");
}

void
sv_pool_boundary_estep(SvPool *mp, SvSet *set)
{
    gint i, j, k;
    gint xres = set->sp.xres;
    gint yres = set->sp.yres;
    gint zres = set->sp.zres;
    gdouble dx = set->sp.dx;
    gdouble dy = set->sp.dy;
    gdouble dz = set->sp.dz;
    gdouble dt = set->plan.dt;
    gdouble ind;
    gdouble  be = 0, ce = 0;
    gdouble cax, cay, caz, cbx, cby, cbz;
    gint ir, jr, kr;
    gint xtbc, ytbc, ztbc;

    cax = cay = caz = cbx = cby = cbz = 0;
    ir = jr = kr = 0;
    xtbc = ytbc = ztbc = 0;

    if (set->sc.verbose > 1) {
        printf("Running boundary E step...  \n ");
        fflush(stdout);
    }

    /*pool boundary */

    if (set->sb.bx0 == SV_BOUNDARY_LIAO || set->sb.bxn == SV_BOUNDARY_LIAO) {
	    for (j = 0; j < yres; j++) {
	        for (k = 0; k < zres; k++) {
	    	    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC)
	    	        ind = sqrt(mp->d->epsilon->data[0][j][k]/EPSILON_0);
	    	    else 
                    ind = 1;

	    	    if (set->sb.bx0 == SV_BOUNDARY_LIAO) {
	    	        mp->d->ex->data[0][j][k] = mp->bnd->liao.fbnx->data[EX0P][j][k] +
	    	    	(dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx) * (mp->d->ex->data[1][j][k] - mp->bnd->liao.fbnx->data[EX00][j][k]);
	    	        mp->d->ey->data[0][j][k] = mp->bnd->liao.fbnx->data[EY0P][j][k] +
	    	    	(dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy) * (mp->d->ey->data[1][j][k] - mp->bnd->liao.fbnx->data[EY00][j][k]);
	    	        mp->d->ez->data[0][j][k] = mp->bnd->liao.fbnx->data[EZ0P][j][k] +
	    	    	(dt*LIGHT_SPEED/ind-dz)/(dt*LIGHT_SPEED/ind+dz) * (mp->d->ez->data[1][j][k] - mp->bnd->liao.fbnx->data[EZ00][j][k]);

	    	        mp->d->hx->data[0][j][k] = mp->bnd->liao.fbnx->data[HX0P][j][k] +
	    	    	(dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx) * (mp->d->hx->data[1][j][k] - mp->bnd->liao.fbnx->data[HX00][j][k]);
	    	        mp->d->hy->data[0][j][k] = mp->bnd->liao.fbnx->data[HY0P][j][k] +
	    	    	(dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy) * (mp->d->hy->data[1][j][k] - mp->bnd->liao.fbnx->data[HY00][j][k]);
	    	        mp->d->hz->data[0][j][k] = mp->bnd->liao.fbnx->data[HZ0P][j][k] +
	    	    	(dt*LIGHT_SPEED/ind-dz)/(dt*LIGHT_SPEED/ind+dz) * (mp->d->hz->data[1][j][k] - mp->bnd->liao.fbnx->data[HZ00][j][k]);
	    	    }
	    	    
	    	    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC)
	    	        ind = sqrt(mp->d->epsilon->data[xres-1][j][k]/EPSILON_0);
	    	    else ind = 1;

	    	    if (set->sb.bxn == SV_BOUNDARY_LIAO) {
	    	        mp->d->ex->data[xres-1][j][k]=mp->bnd->liao.fbnx->data[EXNP][j][k] +
	    	    	(dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx) * (mp->d->ex->data[xres-2][j][k] - mp->bnd->liao.fbnx->data[EXN0][j][k]);
	    	        mp->d->ey->data[xres-1][j][k]=mp->bnd->liao.fbnx->data[EYNP][j][k] +
	    	    	(dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy) * (mp->d->ey->data[xres-2][j][k] - mp->bnd->liao.fbnx->data[EYN0][j][k]);
	    	        mp->d->ez->data[xres-1][j][k]=mp->bnd->liao.fbnx->data[EZNP][j][k] +
	    	    	(dt*LIGHT_SPEED/ind-dz)/(dt*LIGHT_SPEED/ind+dz) * (mp->d->ez->data[xres-2][j][k] - mp->bnd->liao.fbnx->data[EZN0][j][k]);

	    	        mp->d->hx->data[xres-1][j][k]=mp->bnd->liao.fbnx->data[HXNP][j][k] +
	    	    	(dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx) * (mp->d->hx->data[xres-2][j][k] - mp->bnd->liao.fbnx->data[HXN0][j][k]);
	    	        mp->d->hy->data[xres-1][j][k]=mp->bnd->liao.fbnx->data[HYNP][j][k] +
	    	    	(dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy) * (mp->d->hy->data[xres-2][j][k] - mp->bnd->liao.fbnx->data[HYN0][j][k]);
	    	        mp->d->hz->data[xres-1][j][k]=mp->bnd->liao.fbnx->data[HZNP][j][k] +
	    	    	(dt*LIGHT_SPEED/ind-dz)/(dt*LIGHT_SPEED/ind+dz) * (mp->d->hz->data[xres-2][j][k] - mp->bnd->liao.fbnx->data[HZN0][j][k]);
	    	    }
	        } // k
	    } // j
    } // if

    if (set->sb.by0 == SV_BOUNDARY_LIAO || set->sb.byn == SV_BOUNDARY_LIAO) {
	    for (i = 0; i < xres; i++) {
	        for (k = 0; k < zres; k++) {
	    	    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC)
	    	        ind = sqrt(mp->d->epsilon->data[i][0][k]/EPSILON_0);
	    	    else 
                    ind = 1;
	    	    
	    	    if (set->sb.by0 == SV_BOUNDARY_LIAO) {
	    	        mp->d->ex->data[i][0][k]=mp->bnd->liao.fbny->data[EX0P][i][k] +
	    	    	(dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx) * (mp->d->ex->data[i][1][k] - mp->bnd->liao.fbny->data[EX00][i][k]);
	    	        mp->d->ey->data[i][0][k]=mp->bnd->liao.fbny->data[EY0P][i][k] +
	    	    	(dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy) * (mp->d->ey->data[i][1][k] - mp->bnd->liao.fbny->data[EY00][i][k]);
	    	        mp->d->ez->data[i][0][k]=mp->bnd->liao.fbny->data[EZ0P][i][k] +
	    	    	(dt*LIGHT_SPEED/ind-dz)/(dt*LIGHT_SPEED/ind+dz) * (mp->d->ez->data[i][1][k] - mp->bnd->liao.fbny->data[EZ00][i][k]);

	    	        mp->d->hx->data[i][0][k]=mp->bnd->liao.fbny->data[HX0P][i][k] +
	    	    	(dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx) * (mp->d->hx->data[i][1][k] - mp->bnd->liao.fbny->data[HX00][i][k]);
	    	        mp->d->hy->data[i][0][k]=mp->bnd->liao.fbny->data[HY0P][i][k] +
	    	    	(dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy) * (mp->d->hy->data[i][1][k] - mp->bnd->liao.fbny->data[HY00][i][k]);
	    	        mp->d->hz->data[i][0][k]=mp->bnd->liao.fbny->data[HZ0P][i][k] +
	    	    	(dt*LIGHT_SPEED/ind-dz)/(dt*LIGHT_SPEED/ind+dz) * (mp->d->hz->data[i][1][k] - mp->bnd->liao.fbny->data[HZ00][i][k]);
	    	    }
	    	    
	    	    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC)
	    	        ind = sqrt(mp->d->epsilon->data[i][yres-1][k]/EPSILON_0);
	    	    else ind = 1;
	    	    
	    	    if (set->sb.byn == SV_BOUNDARY_LIAO) {
	    	        mp->d->ex->data[i][yres-1][k]=mp->bnd->liao.fbny->data[EXNP][i][k] +
	    	    	(dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx) * (mp->d->ex->data[i][yres-2][k] - mp->bnd->liao.fbny->data[EXN0][i][k]);
	    	        mp->d->ey->data[i][yres-1][k]=mp->bnd->liao.fbny->data[EYNP][i][k] +
	    	    	(dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy) * (mp->d->ey->data[i][yres-2][k] - mp->bnd->liao.fbny->data[EYN0][i][k]);
	    	        mp->d->ez->data[i][yres-1][k]=mp->bnd->liao.fbny->data[EZNP][i][k] +
	    	    	(dt*LIGHT_SPEED/ind-dz)/(dt*LIGHT_SPEED/ind+dz) * (mp->d->ez->data[i][yres-2][k] - mp->bnd->liao.fbny->data[EZN0][i][k]);

	    	        mp->d->hx->data[i][yres-1][k]=mp->bnd->liao.fbny->data[HXNP][i][k] +
	    	    	(dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx) * (mp->d->hx->data[i][yres-2][k] - mp->bnd->liao.fbny->data[HXN0][i][k]);
	    	        mp->d->hy->data[i][yres-1][k]=mp->bnd->liao.fbny->data[HYNP][i][k] +
	    	    	(dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy) * (mp->d->hy->data[i][yres-2][k] - mp->bnd->liao.fbny->data[HYN0][i][k]);
	    	        mp->d->hz->data[i][yres-1][k]=mp->bnd->liao.fbny->data[HZNP][i][k] +
	    	    	(dt*LIGHT_SPEED/ind-dz)/(dt*LIGHT_SPEED/ind+dz) * (mp->d->hz->data[i][yres-2][k] - mp->bnd->liao.fbny->data[HZN0][i][k]);
	    	    }
	        } // k
	    } // i
    } // if
    
    if (set->sb.bz0 == SV_BOUNDARY_LIAO || set->sb.bzn == SV_BOUNDARY_LIAO) {
	    for (i = 0; i < xres; i++) {
	        for (j = 0; j < yres; j++) {
	    	    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC)
	    	        ind = sqrt(mp->d->epsilon->data[i][j][0]/EPSILON_0);
	    	    else 
                    ind = 1;
	    	    
	    	    if (set->sb.bz0 == SV_BOUNDARY_LIAO) {
	    	        mp->d->ex->data[i][j][0]=mp->bnd->liao.fbnz->data[EX0P][i][j] +
	    	    	(dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx) * (mp->d->ex->data[i][j][1] - mp->bnd->liao.fbnz->data[EX00][i][j]);
	    	        mp->d->ey->data[i][j][0]=mp->bnd->liao.fbnz->data[EY0P][i][j] +
	    	    	(dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy) * (mp->d->ey->data[i][j][1] - mp->bnd->liao.fbnz->data[EY00][i][j]);
	    	        mp->d->ez->data[i][j][0]=mp->bnd->liao.fbnz->data[EZ0P][i][j] +
	    	    	(dt*LIGHT_SPEED/ind-dz)/(dt*LIGHT_SPEED/ind+dz) * (mp->d->ez->data[i][j][1] - mp->bnd->liao.fbnz->data[EZ00][i][j]);

	    	        mp->d->hx->data[i][j][0]=mp->bnd->liao.fbnz->data[HX0P][i][j] +
	    	    	(dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx) * (mp->d->hx->data[i][j][1] - mp->bnd->liao.fbnz->data[HX00][i][j]);
	    	        mp->d->hy->data[i][j][0]=mp->bnd->liao.fbnz->data[HY0P][i][j] +
	    	    	(dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy) * (mp->d->hy->data[i][j][1] - mp->bnd->liao.fbnz->data[HY00][i][j]);
	    	        mp->d->hz->data[i][j][0]=mp->bnd->liao.fbnz->data[HZ0P][i][j] +
	    	    	(dt*LIGHT_SPEED/ind-dz)/(dt*LIGHT_SPEED/ind+dz) * (mp->d->hz->data[i][j][1] - mp->bnd->liao.fbnz->data[HZ00][i][j]);
	    	    }
	    	    
	    	    if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC)
	    	        ind = sqrt(mp->d->epsilon->data[i][j][zres-2]/EPSILON_0);
	    	    else 
                    ind = 1;
	    	    
	    	    if (set->sb.bzn == SV_BOUNDARY_LIAO) {
	    	        mp->d->ex->data[i][j][zres-1]=mp->bnd->liao.fbnz->data[EXNP][i][j] +
	    	    	(dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx) * (mp->d->ex->data[i][j][zres-2] - mp->bnd->liao.fbnz->data[EXN0][i][j]);
	    	        mp->d->ey->data[i][j][zres-1]=mp->bnd->liao.fbnz->data[EYNP][i][j] +
	    	    	(dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy) * (mp->d->ey->data[i][j][zres-2] - mp->bnd->liao.fbnz->data[EYN0][i][j]);
	    	        mp->d->ez->data[i][j][zres-1]=mp->bnd->liao.fbnz->data[EZNP][i][j] +
	    	    	(dt*LIGHT_SPEED/ind-dz)/(dt*LIGHT_SPEED/ind+dz) * (mp->d->ez->data[i][j][zres-2] - mp->bnd->liao.fbnz->data[EZN0][i][j]);

	    	        mp->d->hx->data[i][j][zres-1]=mp->bnd->liao.fbnz->data[HXNP][i][j] +
	    	    	(dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx) * (mp->d->hx->data[i][j][zres-2] - mp->bnd->liao.fbnz->data[HXN0][i][j]);
	    	        mp->d->hy->data[i][j][zres-1]=mp->bnd->liao.fbnz->data[HYNP][i][j] +
	    	    	(dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy) * (mp->d->hy->data[i][j][zres-2] - mp->bnd->liao.fbnz->data[HYN0][i][j]);
	    	        mp->d->hz->data[i][j][zres-1]=mp->bnd->liao.fbnz->data[HZNP][i][j] +
	    	    	(dt*LIGHT_SPEED/ind-dz)/(dt*LIGHT_SPEED/ind+dz) * (mp->d->hz->data[i][j][zres-2] - mp->bnd->liao.fbnz->data[HZN0][i][j]);
	    	    }
	        } // j
	    } // i
    } // if

    if (set->sb.bx0 == SV_BOUNDARY_CPML) {
        for (j = 0; j < yres; j++) {
            for (k = 0; k < zres; k++) {
                mp->d->ey->data[0][j][k] = mp->d->ez->data[0][j][k] = mp->d->ex->data[0][j][k] = 0;
            }
        }

        for (i = 1; i < set->sb.depth_bx0; i++) {
	    
#pragma omp parallel default(shared) private(j, k, cbx, cby, cbz, be, ce, xtbc, ytbc, ztbc)
#pragma omp for nowait
	    
            for (j = 0; j < yres; j++) {
                for (k = 0; k < zres; k++) {

                    sv_yee_data_get_tbc(mp->d, mp->mats, mp->nmat, set, i, j, k, &xtbc, &ytbc, &ztbc);
            
                    sv_yee_data_get_cabs(mp->d, set, mp->mats, mp->nmat, mp->d->dt, i, j, k, &cax, &cay, &caz, &cbx, &cby, &cbz);
                    be = mp->bnd->cpml.be_x0[i];
                    ce = mp->bnd->cpml.ce_x0[i];
                    
                    mp->bnd->cpml.peyx_x0->data[i][j][k] = be * mp->bnd->cpml.peyx_x0->data[i][j][k] 
                             + ce * ((mp->d->hz->data[i][j][k] - mp->d->hz->data[i-1][j][k]))/dx;
                    mp->bnd->cpml.pezx_x0->data[i][j][k] = be * mp->bnd->cpml.pezx_x0->data[i][j][k] 
                             + ce * ((mp->d->hy->data[i][j][k] - mp->d->hy->data[i-1][j][k]))/dx;

                    if (!(j == 0 || k == 0) && ytbc==1)
			            mp->d->ey->data[i][j][k] -= dx * cby * mp->bnd->cpml.peyx_x0->data[i][j][k];
                    if (!(j == 0 || k == 0) && ztbc==1)
			            mp->d->ez->data[i][j][k] += dx * cbz * mp->bnd->cpml.pezx_x0->data[i][j][k]; 
                } // k
            } // j
        } // i
    } // if


    if (set->sb.by0 == SV_BOUNDARY_CPML) {
        for (i = 0; i < xres; i++) {
            for (k = 0; k < zres; k++) {
                 mp->d->ey->data[i][0][k] = mp->d->ez->data[i][0][k] = mp->d->ex->data[i][0][k] = 0;
            }
        }

        for (j = 1; j < set->sb.depth_by0; j++) {
	    
#pragma omp parallel default(shared) private(i, k, cbx, cby, cbz, be, ce, xtbc, ytbc, ztbc)
#pragma omp for nowait
	    
            for (i = 0; i < xres; i++) {
                for (k = 0; k < zres; k++) {

                    sv_yee_data_get_tbc(mp->d, mp->mats, mp->nmat, set, i, j, k, &xtbc, &ytbc, &ztbc);

                    sv_yee_data_get_cabs(mp->d, set, mp->mats, mp->nmat, mp->d->dt, i, j, k, &cax, &cay, &caz, &cbx, &cby, &cbz);
                    be = mp->bnd->cpml.be_y0[j];
                    ce = mp->bnd->cpml.ce_y0[j];

                    mp->bnd->cpml.pexy_y0->data[i][j][k] = be*mp->bnd->cpml.pexy_y0->data[i][j][k] 
                               + ce * ((mp->d->hz->data[i][j][k] - mp->d->hz->data[i][j-1][k]))/dy;
                    mp->bnd->cpml.pezy_y0->data[i][j][k] = be*mp->bnd->cpml.pezy_y0->data[i][j][k] 
                               + ce * ((mp->d->hx->data[i][j][k] - mp->d->hx->data[i][j-1][k]))/dy;

                    if (!(i == 0 || k == 0) && xtbc==1)
			            mp->d->ex->data[i][j][k] += dy * cbx * mp->bnd->cpml.pexy_y0->data[i][j][k];
                    if (!(i == 0 || k == 0) && ztbc==1)
			            mp->d->ez->data[i][j][k] -= dy * cbz * mp->bnd->cpml.pezy_y0->data[i][j][k]; 
                } // k
            } // i
        } // j
    } // if


    if (set->sb.bz0 == SV_BOUNDARY_CPML) {

       for (i = 0; i < xres; i++) {
            for (j = 0; j < yres; j++) {
                mp->d->ey->data[i][j][0] = mp->d->ez->data[i][j][0] = mp->d->ex->data[i][j][0] = 0;
            }
       }

        for (k = 1; k < set->sb.depth_bz0; k++) {
	    
#pragma omp parallel default(shared) private(i, j, cbx, cby, cbz, be, ce, xtbc, ytbc, ztbc)
#pragma omp for nowait
	    
            for (i = 0; i < xres; i++) {
                for (j = 0; j < yres; j++) {

                    sv_yee_data_get_tbc(mp->d, mp->mats, mp->nmat, set, i, j, k, &xtbc, &ytbc, &ztbc);
 
                    sv_yee_data_get_cabs(mp->d, set, mp->mats, mp->nmat, mp->d->dt, i, j, k, &cax, &cay, &caz, &cbx, &cby, &cbz);
                    be = mp->bnd->cpml.be_z0[k];
                    ce = mp->bnd->cpml.ce_z0[k];

                    mp->bnd->cpml.peyz_z0->data[i][j][k] = be*mp->bnd->cpml.peyz_z0->data[i][j][k] 
                            + ce * ((mp->d->hx->data[i][j][k] - mp->d->hx->data[i][j][k-1]))/dz;
                    mp->bnd->cpml.pexz_z0->data[i][j][k] = be*mp->bnd->cpml.pexz_z0->data[i][j][k] 
                            + ce * ((mp->d->hy->data[i][j][k] - mp->d->hy->data[i][j][k-1]))/dz;

                    if (!(i == 0 || j == 0) && ytbc==1)
			            mp->d->ey->data[i][j][k] += dx * cby * mp->bnd->cpml.peyz_z0->data[i][j][k];
                    if (!(i == 0 || j == 0) && ztbc==1)
			            mp->d->ex->data[i][j][k] -= dx * cbx * mp->bnd->cpml.pexz_z0->data[i][j][k]; 
                } // j
            } // i
        } // k
    } // if
    if (set->sb.bxn == SV_BOUNDARY_CPML) {
        for (i = (xres - set->sb.depth_bxn + 1); i < xres; i++) {
	    
#pragma omp parallel default(shared) private(j, k, cbx, cby, cbz, be, ce, ir, xtbc, ytbc, ztbc)
#pragma omp for nowait
	    
            for (j = 0; j < yres; j++) {
                for (k = 0; k < zres; k++) {
                    ir = i - (xres - set->sb.depth_bxn);

                    sv_yee_data_get_tbc(mp->d, mp->mats, mp->nmat, set, i, j, k, &xtbc, &ytbc, &ztbc);

                    sv_yee_data_get_cabs(mp->d, set, mp->mats, mp->nmat, mp->d->dt, i, j, k, &cax, &cay, &caz, &cbx, &cby, &cbz);
                    be = mp->bnd->cpml.be_xn[ir];
                    ce = mp->bnd->cpml.ce_xn[ir];
 
                    mp->bnd->cpml.peyx_xn->data[ir][j][k] = be*mp->bnd->cpml.peyx_xn->data[ir][j][k] 
                                                              + ce * ((mp->d->hz->data[i][j][k] - mp->d->hz->data[i-1][j][k]))/dx;
                    mp->bnd->cpml.pezx_xn->data[ir][j][k] = be*mp->bnd->cpml.pezx_xn->data[ir][j][k] 
                                                              + ce * ((mp->d->hy->data[i][j][k] - mp->d->hy->data[i-1][j][k]))/dx;

                    if (!(j == 0 || k == 0) && ytbc==1)
			            mp->d->ey->data[i][j][k] -= dx * cby * mp->bnd->cpml.peyx_xn->data[ir][j][k];
                    if (!(j == 0 || k == 0) && ztbc==1)
			            mp->d->ez->data[i][j][k] += dx * cbz * mp->bnd->cpml.pezx_xn->data[ir][j][k]; 

                } // k
            } // j
        } // i
    } // if
    if (set->sb.byn == SV_BOUNDARY_CPML) {
        for (j = (yres - set->sb.depth_byn + 1); j < (yres); j++) {
	    
#pragma omp parallel default(shared) private(i, k, cbx, cby, cbz, be, ce, jr, xtbc, ytbc, ztbc)
#pragma omp for nowait
	    
            for (i = 0; i < xres; i++) {
                for (k = 0; k < zres; k++) {
                    jr = j - (yres - set->sb.depth_byn);

                    sv_yee_data_get_tbc(mp->d, mp->mats, mp->nmat, set, i, j, k, &xtbc, &ytbc, &ztbc);

                    sv_yee_data_get_cabs(mp->d, set, mp->mats, mp->nmat, mp->d->dt, i, j, k, &cax, &cay, &caz, &cbx, &cby, &cbz);
                    be = mp->bnd->cpml.be_yn[jr];
                    ce = mp->bnd->cpml.ce_yn[jr];


                    mp->bnd->cpml.pexy_yn->data[i][jr][k] = be*mp->bnd->cpml.pexy_yn->data[i][jr][k] 
                           + ce * ((mp->d->hz->data[i][j][k] - mp->d->hz->data[i][j-1][k]))/dy;
                    mp->bnd->cpml.pezy_yn->data[i][jr][k] = be*mp->bnd->cpml.pezy_yn->data[i][jr][k] 
                           + ce * ((mp->d->hx->data[i][j][k] - mp->d->hx->data[i][j-1][k]))/dy;

                    if (!(i == 0 || k == 0) && xtbc==1)
			            mp->d->ex->data[i][j][k] += dy * cbx * mp->bnd->cpml.pexy_yn->data[i][jr][k];
                    if (!(i == 0 || k == 0) && ztbc==1)
			            mp->d->ez->data[i][j][k] -= dy * cbz * mp->bnd->cpml.pezy_yn->data[i][jr][k]; 
                } // k
            } // i
        } // j
    } // if
 
    if (set->sb.bzn == SV_BOUNDARY_CPML) {
        for (k = (zres - set->sb.depth_bzn + 1); k < zres; k++) {
	    
#pragma omp parallel default(shared) private(i, j, cbx, cby, cbz, be, ce, kr, xtbc, ytbc, ztbc)
#pragma omp for nowait
	    
            for (i = 0; i < xres; i++) {
                for (j = 0; j < yres; j++) {
                    kr = k - (zres - set->sb.depth_bzn);

                    sv_yee_data_get_tbc(mp->d, mp->mats, mp->nmat, set, i, j, k, &xtbc, &ytbc, &ztbc);

                    sv_yee_data_get_cabs(mp->d, set, mp->mats, mp->nmat, mp->d->dt, i, j, k, &cax, &cay, &caz, &cbx, &cby, &cbz);
                    be = mp->bnd->cpml.be_zn[kr];
                    ce = mp->bnd->cpml.ce_zn[kr];
           
                    mp->bnd->cpml.peyz_zn->data[i][j][kr] = be*mp->bnd->cpml.peyz_zn->data[i][j][kr] 
                           + ce * ((mp->d->hx->data[i][j][k] - mp->d->hx->data[i][j][k-1]))/dz;
                    mp->bnd->cpml.pexz_zn->data[i][j][kr] = be*mp->bnd->cpml.pexz_zn->data[i][j][kr] 
                           + ce * ((mp->d->hy->data[i][j][k] - mp->d->hy->data[i][j][k-1]))/dz;

                    if (!(i == 0 || j == 0) && ytbc==1)
			            mp->d->ey->data[i][j][k] += dx * cby * mp->bnd->cpml.peyz_zn->data[i][j][kr];
                    if (!(i == 0 || j == 0) && xtbc==1)
			            mp->d->ex->data[i][j][k] -= dx * cbx * mp->bnd->cpml.pexz_zn->data[i][j][kr]; 
                } // j
            } // i
        } // k
    } // if
    
    /*material boundary*/
    if (set->smb.bx0 == SV_BOUNDARY_PERIODIC || set->smb.bxn == SV_BOUNDARY_PERIODIC) {
        for (j = set->smb.by0pos; j < set->smb.bynpos; j++) {
            for (k = set->smb.bz0pos; k < set->smb.bznpos; k++) {
		//      mp->d->ex->data[set->smb.bxnpos][j][k] = mp->d->ex->data[set->smb.bx0pos][j][k];
                mp->d->ey->data[set->smb.bxnpos][j][k] = mp->d->ey->data[set->smb.bx0pos][j][k];
                mp->d->ez->data[set->smb.bxnpos][j][k] = mp->d->ez->data[set->smb.bx0pos][j][k];

		//      mp->d->ex->data[set->smb.bxnpos+1][j][k] = mp->d->ex->data[set->smb.bx0pos+1][j][k];
		//      mp->d->ey->data[set->smb.bxnpos+1][j][k] = mp->d->ey->data[set->smb.bx0pos+1][j][k];
		//      mp->d->ez->data[set->smb.bxnpos+1][j][k] = mp->d->ez->data[set->smb.bx0pos+1][j][k];
	        }
        }
    }
    
    if (set->smb.by0 == SV_BOUNDARY_PERIODIC || set->smb.byn == SV_BOUNDARY_PERIODIC) {
        /*y planes*/
        for (i = set->smb.bx0pos; i < set->smb.bxnpos; i++) {
            for (k = set->smb.bz0pos; k < set->smb.bznpos; k++) {
                mp->d->ex->data[i][set->smb.bynpos][k] = mp->d->ex->data[i][set->smb.by0pos][k];
		//      mp->d->ey->data[i][set->smb.by0pos][k] = mp->d->ey->data[i][set->smb.bynpos][k];
                mp->d->ez->data[i][set->smb.bynpos][k] = mp->d->ez->data[i][set->smb.by0pos][k];

		//      mp->d->ex->data[i][set->smb.bynpos+1][k] = mp->d->ex->data[i][set->smb.by0pos+1][k];
		//      mp->d->ey->data[i][set->smb.bynpos+1][k] = mp->d->ey->data[i][set->smb.by0pos+1][k];
		//      mp->d->ez->data[i][set->smb.bynpos+1][k] = mp->d->ez->data[i][set->smb.by0pos+1][k];
	        }
        }
    }
    
    if (set->smb.bz0 == SV_BOUNDARY_PERIODIC || set->smb.bzn == SV_BOUNDARY_PERIODIC) {
        /*z planes*/
        for (i = set->smb.bx0pos; i < set->smb.bxnpos; i++) {
            for (j = set->smb.by0pos; j < set->smb.bynpos; j++) {
                mp->d->ex->data[i][j][set->smb.bznpos] = mp->d->ex->data[i][j][set->smb.bz0pos];
                mp->d->ey->data[i][j][set->smb.bznpos] = mp->d->ey->data[i][j][set->smb.bz0pos];
		//      mp->d->ez->data[i][j][set->smb.bznpos] = mp->d->ez->data[i][j][set->smb.bz0pos];

		//      mp->d->ex->data[i][j][set->smb.bznpos+1] = mp->d->ex->data[i][j][set->smb.bz0pos+1];
		//      mp->d->ey->data[i][j][set->smb.bznpos+1] = mp->d->ey->data[i][j][set->smb.bz0pos+1];
		//      mp->d->ez->data[i][j][set->smb.bznpos+1] = mp->d->ez->data[i][j][set->smb.bz0pos+1];
            }
        }
    }

    if (set->sc.verbose > 1) 
        printf("done.\n");
}

void
sv_pool_boundary_copy(SvPool *mp, SvSet *set)
{
    gint i, j, k;
    gint xres = set->sp.xres;
    gint yres = set->sp.yres;
    gint zres = set->sp.zres;

    if (set->sc.verbose > 1) {
        printf("Running boundary copy...   ");
        fflush(stdout);
    }

    if (set->sb.bx0 == SV_BOUNDARY_LIAO || set->sb.bxn == SV_BOUNDARY_LIAO) {
        for (j = 0; j < yres; j++) {
            for (k = 0; k < zres; k++) {
                if (set->sb.bx0 == SV_BOUNDARY_LIAO) {
                    mp->bnd->liao.fbnx->data[EX00][j][k] = mp->d->ex->data[0][j][k];
                    mp->bnd->liao.fbnx->data[EY00][j][k] = mp->d->ey->data[0][j][k];
                    mp->bnd->liao.fbnx->data[EZ00][j][k] = mp->d->ez->data[0][j][k];
                    mp->bnd->liao.fbnx->data[HX00][j][k] = mp->d->hx->data[0][j][k];
                    mp->bnd->liao.fbnx->data[HY00][j][k] = mp->d->hy->data[0][j][k];
                    mp->bnd->liao.fbnx->data[HZ00][j][k] = mp->d->hz->data[0][j][k];

                    mp->bnd->liao.fbnx->data[EX0P][j][k] = mp->d->ex->data[1][j][k];
                    mp->bnd->liao.fbnx->data[EY0P][j][k] = mp->d->ey->data[1][j][k];
                    mp->bnd->liao.fbnx->data[EZ0P][j][k] = mp->d->ez->data[1][j][k];
                    mp->bnd->liao.fbnx->data[HX0P][j][k] = mp->d->hx->data[1][j][k];
                    mp->bnd->liao.fbnx->data[HY0P][j][k] = mp->d->hy->data[1][j][k];
                    mp->bnd->liao.fbnx->data[HZ0P][j][k] = mp->d->hz->data[1][j][k];
                }
                if (set->sb.bxn == SV_BOUNDARY_LIAO) {
                    mp->bnd->liao.fbnx->data[EXN0][j][k] = mp->d->ex->data[xres-1][j][k];
                    mp->bnd->liao.fbnx->data[EYN0][j][k] = mp->d->ey->data[xres-1][j][k];
                    mp->bnd->liao.fbnx->data[EZN0][j][k] = mp->d->ez->data[xres-1][j][k];
                    mp->bnd->liao.fbnx->data[HXN0][j][k] = mp->d->hx->data[xres-1][j][k];
                    mp->bnd->liao.fbnx->data[HYN0][j][k] = mp->d->hy->data[xres-1][j][k];
                    mp->bnd->liao.fbnx->data[HZN0][j][k] = mp->d->hz->data[xres-1][j][k];

                    mp->bnd->liao.fbnx->data[EXNP][j][k] = mp->d->ex->data[xres-2][j][k];
                    mp->bnd->liao.fbnx->data[EYNP][j][k] = mp->d->ey->data[xres-2][j][k];
                    mp->bnd->liao.fbnx->data[EZNP][j][k] = mp->d->ez->data[xres-2][j][k];
                    mp->bnd->liao.fbnx->data[HXNP][j][k] = mp->d->hx->data[xres-2][j][k];
                    mp->bnd->liao.fbnx->data[HYNP][j][k] = mp->d->hy->data[xres-2][j][k];
                    mp->bnd->liao.fbnx->data[HZNP][j][k] = mp->d->hz->data[xres-2][j][k];
                }
	        } // k
        } // j
    } // if

    if (set->sb.by0 == SV_BOUNDARY_LIAO || set->sb.byn == SV_BOUNDARY_LIAO) {
        for (i = 0; i < xres; i++) {
            for (k = 0; k < zres; k++) {
                if (set->sb.by0 == SV_BOUNDARY_LIAO) {
                    mp->bnd->liao.fbny->data[EX00][i][k] = mp->d->ex->data[i][0][k];
                    mp->bnd->liao.fbny->data[EY00][i][k] = mp->d->ey->data[i][0][k];
                    mp->bnd->liao.fbny->data[EZ00][i][k] = mp->d->ez->data[i][0][k];
                    mp->bnd->liao.fbny->data[HX00][i][k] = mp->d->hx->data[i][0][k];
                    mp->bnd->liao.fbny->data[HY00][i][k] = mp->d->hy->data[i][0][k];
                    mp->bnd->liao.fbny->data[HZ00][i][k] = mp->d->hz->data[i][0][k];

                    mp->bnd->liao.fbny->data[EX0P][i][k] = mp->d->ex->data[i][1][k];
                    mp->bnd->liao.fbny->data[EY0P][i][k] = mp->d->ey->data[i][1][k];
                    mp->bnd->liao.fbny->data[EZ0P][i][k] = mp->d->ez->data[i][1][k];
                    mp->bnd->liao.fbny->data[HX0P][i][k] = mp->d->hx->data[i][1][k];
                    mp->bnd->liao.fbny->data[HY0P][i][k] = mp->d->hy->data[i][1][k];
                    mp->bnd->liao.fbny->data[HZ0P][i][k] = mp->d->hz->data[i][1][k];
                }
                if (set->sb.byn == SV_BOUNDARY_LIAO) {
                    mp->bnd->liao.fbny->data[EXN0][i][k] = mp->d->ex->data[i][yres-1][k];
                    mp->bnd->liao.fbny->data[EYN0][i][k] = mp->d->ey->data[i][yres-1][k];
                    mp->bnd->liao.fbny->data[EZN0][i][k] = mp->d->ez->data[i][yres-1][k];
                    mp->bnd->liao.fbny->data[HXN0][i][k] = mp->d->hx->data[i][yres-1][k];
                    mp->bnd->liao.fbny->data[HYN0][i][k] = mp->d->hy->data[i][yres-1][k];
                    mp->bnd->liao.fbny->data[HZN0][i][k] = mp->d->hz->data[i][yres-1][k];

                    mp->bnd->liao.fbny->data[EXNP][i][k] = mp->d->ex->data[i][yres-2][k];
                    mp->bnd->liao.fbny->data[EYNP][i][k] = mp->d->ey->data[i][yres-2][k];
                    mp->bnd->liao.fbny->data[EZNP][i][k] = mp->d->ez->data[i][yres-2][k];
                    mp->bnd->liao.fbny->data[HXNP][i][k] = mp->d->hx->data[i][yres-2][k];
                    mp->bnd->liao.fbny->data[HYNP][i][k] = mp->d->hy->data[i][yres-2][k];
                    mp->bnd->liao.fbny->data[HZNP][i][k] = mp->d->hz->data[i][yres-2][k];
		        }
            } // k
        } // i
    } // if
    
    if (set->sb.bz0 == SV_BOUNDARY_LIAO || set->sb.bzn == SV_BOUNDARY_LIAO) {
        for (i = 0; i < xres; i++) {
            for (j = 0; j < yres; j++) {
                if (set->sb.bz0 == SV_BOUNDARY_LIAO) {
                    mp->bnd->liao.fbnz->data[EX00][i][j] = mp->d->ex->data[i][j][0];
                    mp->bnd->liao.fbnz->data[EY00][i][j] = mp->d->ey->data[i][j][0];
                    mp->bnd->liao.fbnz->data[EZ00][i][j] = mp->d->ez->data[i][j][0];
                    mp->bnd->liao.fbnz->data[HX00][i][j] = mp->d->hx->data[i][j][0];
                    mp->bnd->liao.fbnz->data[HY00][i][j] = mp->d->hy->data[i][j][0];
                    mp->bnd->liao.fbnz->data[HZ00][i][j] = mp->d->hz->data[i][j][0];

                    mp->bnd->liao.fbnz->data[EX0P][i][j] = mp->d->ex->data[i][j][1];
                    mp->bnd->liao.fbnz->data[EY0P][i][j] = mp->d->ey->data[i][j][1];
                    mp->bnd->liao.fbnz->data[EZ0P][i][j] = mp->d->ez->data[i][j][1];
                    mp->bnd->liao.fbnz->data[HX0P][i][j] = mp->d->hx->data[i][j][1];
                    mp->bnd->liao.fbnz->data[HY0P][i][j] = mp->d->hy->data[i][j][1];
                    mp->bnd->liao.fbnz->data[HZ0P][i][j] = mp->d->hz->data[i][j][1];
                }
                if (set->sb.bzn == SV_BOUNDARY_LIAO) {
                    mp->bnd->liao.fbnz->data[EXN0][i][j] = mp->d->ex->data[i][j][zres-1];
                    mp->bnd->liao.fbnz->data[EYN0][i][j] = mp->d->ey->data[i][j][zres-1];
                    mp->bnd->liao.fbnz->data[EZN0][i][j] = mp->d->ez->data[i][j][zres-1];
                    mp->bnd->liao.fbnz->data[HXN0][i][j] = mp->d->hx->data[i][j][zres-1];
                    mp->bnd->liao.fbnz->data[HYN0][i][j] = mp->d->hy->data[i][j][zres-1];
                    mp->bnd->liao.fbnz->data[HZN0][i][j] = mp->d->hz->data[i][j][zres-1];

                    mp->bnd->liao.fbnz->data[EXNP][i][j] = mp->d->ex->data[i][j][zres-2];
                    mp->bnd->liao.fbnz->data[EYNP][i][j] = mp->d->ey->data[i][j][zres-2];
                    mp->bnd->liao.fbnz->data[EZNP][i][j] = mp->d->ez->data[i][j][zres-2];
                    mp->bnd->liao.fbnz->data[HXNP][i][j] = mp->d->hx->data[i][j][zres-2];
                    mp->bnd->liao.fbnz->data[HYNP][i][j] = mp->d->hy->data[i][j][zres-2];
                    mp->bnd->liao.fbnz->data[HZNP][i][j] = mp->d->hz->data[i][j][zres-2];
                }
            } // j
        } // i
    } // if

    if (set->sc.verbose > 1) 
        printf("done.\n");
}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
