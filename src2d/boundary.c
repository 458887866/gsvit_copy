
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


/*  boundary.c : 
 *  boundary conditions implementation
 *
 */

#include "boundary.h"
#include "pool.h"
#include "settings.h"
#include "constants.h"
#include <math.h>


SvBoundary* 
sv_boundary_new(SvSet *set)
{
    SvBoundary* bnd = (SvBoundary*)g_malloc(sizeof(SvBoundary));
    if (set->sb.bx0 == SV_BOUNDARY_LIAO || set->sb.bxn == SV_BOUNDARY_LIAO) {
        bnd->liao.fbnx = sv_dcube_new(NPLANES, set->sp.yres, 1, 
                                      NPLANES, set->sp.yres, 1, 1);
        if (set->sc.verbose) printf("X Liao boundary initialized\n");
    }

    if (set->sb.by0 == SV_BOUNDARY_LIAO || set->sb.byn == SV_BOUNDARY_LIAO) {
        bnd->liao.fbny = sv_dcube_new(NPLANES, set->sp.xres, 1,
                                      NPLANES, set->sp.xres, 1, 1);
        if (set->sc.verbose) printf("Y Liao boundary initialized\n");
    }

    return bnd;
}

void sv_pool_boundary_hstep(SvPool *mp, SvSet *set)
{
    //gint i, j, k;
    
    if (set->sc.verbose) {
        printf("Running boundary H step...   ");
        fflush(stdout);
    }

    /*material boundary*/
/*    if (set->smb.bx0 == SV_BOUNDARY_PERIODIC || set->smb.bxn == SV_BOUNDARY_PERIODIC) {
        for (j=set->smb.by0pos; j<set->smb.bynpos; j++)
        {
            for (k=set->smb.bz0pos; k<set->smb.bznpos; k++)
            {
                mp->hx->data[set->smb.bxnpos-1][j][k] = mp->hx->data[set->smb.bx0pos+1][j][k];
                mp->hy->data[set->smb.bxnpos-1][j][k] = mp->hy->data[set->smb.bx0pos+1][j][k];
                mp->hz->data[set->smb.bxnpos-1][j][k] = mp->hz->data[set->smb.bx0pos+1][j][k];
            }
        }
    } 
    if (set->smb.by0 == SV_BOUNDARY_PERIODIC || set->smb.byn == SV_BOUNDARY_PERIODIC) {
*/
        /*y planes*/
/*        for (i=set->smb.bx0pos; i<set->smb.bxnpos; i++)
        {
            for (k=set->smb.bz0pos; k<set->smb.bznpos; k++)
            {
                mp->hx->data[i][set->smb.bynpos-1][k] = mp->hx->data[i][set->smb.by0pos+1][k];
                mp->hy->data[i][set->smb.bynpos-1][k] = mp->hy->data[i][set->smb.by0pos+1][k];
                mp->hz->data[i][set->smb.bynpos-1][k] = mp->hz->data[i][set->smb.by0pos+1][k];
            }
        }

    }
    if (set->smb.bz0 == SV_BOUNDARY_PERIODIC || set->smb.bzn == SV_BOUNDARY_PERIODIC) {
    */
        /*z planes*/
  /*      for (i=set->smb.bx0pos; i<set->smb.bxnpos; i++)
        {
            for (j=set->smb.by0pos; j<set->smb.bynpos; j++)
            {
                mp->hx->data[i][j][set->smb.bznpos-1] = mp->hx->data[i][j][set->smb.bz0pos+1];
                mp->hy->data[i][j][set->smb.bznpos-1] = mp->hy->data[i][j][set->smb.bz0pos+1];
                mp->hz->data[i][j][set->smb.bznpos-1] = mp->hz->data[i][j][set->smb.bz0pos+1];
            }
        }
    }
*/

    if (set->sc.verbose) printf("done.\n");

}

void sv_pool_boundary_estep(SvPool *mp, SvSet *set)
{
    gint i, j;
    gint xres = set->sp.xres;
    gint yres = set->sp.yres;
    gdouble dx = set->sp.dx;
    gdouble dy = set->sp.dy;
    gdouble dt = set->plan.dt;
    gdouble ind;
    

    if (set->sc.verbose) {
        printf("Running boundary E step...   ");
        fflush(stdout);
    }

    /*pool boundary */

    if (set->sb.bx0 == SV_BOUNDARY_LIAO || set->sb.bxn == SV_BOUNDARY_LIAO) {
        for (j=0; j<yres; j++) {
                ind = sqrt(mp->epsilon->data[0+j*xres]/EPSILON_0);

                if (set->sb.bx0 == SV_BOUNDARY_LIAO) {
                    
                    mp->ex->data[0+j*xres]=mp->bnd->liao.fbnx->data[EX0P][j][0] +
                        (dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx)*
                        (mp->ex->data[1+j*xres] - mp->bnd->liao.fbnx->data[EX00][j][0]);
                    mp->ey->data[0+j*xres]=mp->bnd->liao.fbnx->data[EY0P][j][0] +
                        (dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy)*
                        (mp->ey->data[1+j*xres] - mp->bnd->liao.fbnx->data[EY00][j][0]);
                    mp->ez->data[0+j*xres]=mp->bnd->liao.fbnx->data[EZ0P][j][0] +
                        (dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx)*
                        (mp->ez->data[1+j*xres] - mp->bnd->liao.fbnx->data[EZ00][j][0]);

                }
                ind = sqrt(mp->epsilon->data[(xres-1)+xres*j]/EPSILON_0);
                if (set->sb.bxn == SV_BOUNDARY_LIAO) {

                    mp->ex->data[(xres-1)+xres*j]=mp->bnd->liao.fbnx->data[EXNP][j][0] +
                        (dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx)*
                        (mp->ex->data[(xres-2)+xres*j] - mp->bnd->liao.fbnx->data[EXN0][j][0]);
                    mp->ey->data[(xres-1)+xres*j]=mp->bnd->liao.fbnx->data[EYNP][j][0] +
                        (dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy)*
                        (mp->ey->data[(xres-2)+xres*j] - mp->bnd->liao.fbnx->data[EYN0][j][0]);
                    mp->ez->data[(xres-1)+xres*j]=mp->bnd->liao.fbnx->data[EZNP][j][0] +
                        (dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx)*
                        (mp->ez->data[(xres-2)+xres*j] - mp->bnd->liao.fbnx->data[EZN0][j][0]);

                    mp->hx->data[(xres-1)+xres*j]=mp->bnd->liao.fbnx->data[HXNP][j][0] +
                        (dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx)*
                        (mp->hx->data[(xres-2)+xres*j] - mp->bnd->liao.fbnx->data[HXN0][j][0]);
                    mp->hy->data[(xres-1)+xres*j]=mp->bnd->liao.fbnx->data[HYNP][j][0] +
                        (dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy)*
                        (mp->hy->data[(xres-2)+xres*j] - mp->bnd->liao.fbnx->data[HYN0][j][0]);
                    mp->hz->data[(xres-1)+xres*j]=mp->bnd->liao.fbnx->data[HZNP][j][0] +
                        (dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx)*
                        (mp->hz->data[(xres-2)+xres*j] - mp->bnd->liao.fbnx->data[HZN0][j][0]);

                }
        }
    }
    
    if (set->sb.by0 == SV_BOUNDARY_LIAO || set->sb.byn == SV_BOUNDARY_LIAO) {
        for (i=0; i<xres; i++) {

                ind = sqrt(mp->epsilon->data[i]/EPSILON_0);
                if (set->sb.by0 == SV_BOUNDARY_LIAO) {
                    
                    mp->ex->data[i]=mp->bnd->liao.fbny->data[EX0P][i][0] +
                        (dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx)*
                        (mp->ex->data[i+xres] - mp->bnd->liao.fbny->data[EX00][i][0]);
                    mp->ey->data[i]=mp->bnd->liao.fbny->data[EY0P][i][0] +
                        (dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy)*
                        (mp->ey->data[i+xres] - mp->bnd->liao.fbny->data[EY00][i][0]);
                    mp->ez->data[i]=mp->bnd->liao.fbny->data[EZ0P][i][0] +
                        (dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx)*
                        (mp->ez->data[i+xres] - mp->bnd->liao.fbny->data[EZ00][i][0]);

                }
                ind = sqrt(mp->epsilon->data[i+(yres-1)*xres]/EPSILON_0);
                if (set->sb.byn == SV_BOUNDARY_LIAO) {

                    mp->ex->data[i+(yres-1)*xres]=mp->bnd->liao.fbny->data[EXNP][i][0] +
                        (dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx)*
                        (mp->ex->data[i+(yres-2)*xres] - mp->bnd->liao.fbny->data[EXN0][i][0]);
                    mp->ey->data[i+(yres-1)*xres]=mp->bnd->liao.fbny->data[EYNP][i][0] +
                        (dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy)*
                        (mp->ey->data[i+(yres-2)*xres] - mp->bnd->liao.fbny->data[EYN0][i][0]);
                    mp->ez->data[i+(yres-1)*xres]=mp->bnd->liao.fbny->data[EZNP][i][0] +
                        (dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx)*
                        (mp->ez->data[i+(yres-2)*xres] - mp->bnd->liao.fbny->data[EZN0][i][0]);

                    mp->hx->data[i+(yres-1)*xres]=mp->bnd->liao.fbny->data[HXNP][i][0] +
                        (dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx)*
                        (mp->hx->data[i+(yres-2)*xres] - mp->bnd->liao.fbny->data[HXN0][i][0]);
                    mp->hy->data[i+(yres-1)*xres]=mp->bnd->liao.fbny->data[HYNP][i][0] +
                        (dt*LIGHT_SPEED/ind-dy)/(dt*LIGHT_SPEED/ind+dy)*
                        (mp->hy->data[i+(yres-2)*xres] - mp->bnd->liao.fbny->data[HYN0][i][0]);
                    mp->hz->data[i+(yres-1)*xres]=mp->bnd->liao.fbny->data[HZNP][i][0] +
                        (dt*LIGHT_SPEED/ind-dx)/(dt*LIGHT_SPEED/ind+dx)*
                        (mp->hz->data[i+(yres-2)*xres] - mp->bnd->liao.fbny->data[HZN0][i][0]);

                 }
        }
    }
    /*material boundary*/
/*    if (set->smb.bx0 == SV_BOUNDARY_PERIODIC || set->smb.bxn == SV_BOUNDARY_PERIODIC) {
        for (j=set->smb.by0pos; j<set->smb.bynpos; j++)
        {
            for (k=set->smb.bz0pos; k<set->smb.bznpos; k++)
            {
                mp->ex->data[set->smb.bx0pos+1][j][k] = mp->ex->data[set->smb.bxnpos-1][j][k];
                mp->ey->data[set->smb.bx0pos+1][j][k] = mp->ey->data[set->smb.bxnpos-1][j][k];
                mp->ez->data[set->smb.bx0pos+1][j][k] = mp->ez->data[set->smb.bxnpos-1][j][k];
            }
        }
    } 
    if (set->smb.by0 == SV_BOUNDARY_PERIODIC || set->smb.byn == SV_BOUNDARY_PERIODIC) {
*/
        /*y planes*/
/*        for (i=set->smb.bx0pos; i<set->smb.bxnpos; i++)
        {
            for (k=set->smb.bz0pos; k<set->smb.bznpos; k++)
            {
                mp->ex->data[i][set->smb.by0pos+1][k] = mp->ex->data[i][set->smb.bynpos-1][k];
                mp->ey->data[i][set->smb.by0pos+1][k] = mp->ey->data[i][set->smb.bynpos-1][k];
                mp->ez->data[i][set->smb.by0pos+1][k] = mp->ez->data[i][set->smb.bynpos-1][k];
            }
        }

    }
    if (set->smb.bz0 == SV_BOUNDARY_PERIODIC || set->smb.bzn == SV_BOUNDARY_PERIODIC) {
 */   
        /*z planes*/
 /*       for (i=set->smb.bx0pos; i<set->smb.bxnpos; i++)
        {
            for (j=set->smb.by0pos; j<set->smb.bynpos; j++)
            {
                mp->ex->data[i][j][set->smb.bz0pos+1] = mp->ex->data[i][j][set->smb.bznpos-1];
                mp->ey->data[i][j][set->smb.bz0pos+1] = mp->ey->data[i][j][set->smb.bznpos-1];
                mp->ez->data[i][j][set->smb.bz0pos+1] = mp->ez->data[i][j][set->smb.bznpos-1];
            }
        }
    }
*/
    if (set->sc.verbose) printf("done.\n");

}

void sv_pool_boundary_copy(SvPool *mp, SvSet *set)
{
    gint i, j;
    gint xres = set->sp.xres;
    gint yres = set->sp.yres;

    if (set->sc.verbose) {
        printf("Running boundary copy...   ");
        fflush(stdout);
    }


    if (set->sb.bx0 == SV_BOUNDARY_LIAO || set->sb.bxn == SV_BOUNDARY_LIAO) {
        for (j=0; j<yres; j++) {
            if (set->sb.bx0 == SV_BOUNDARY_LIAO) {
                mp->bnd->liao.fbnx->data[EX00][j][0]=mp->ex->data[0+j*xres];
                mp->bnd->liao.fbnx->data[EY00][j][0]=mp->ey->data[0+j*xres];
                mp->bnd->liao.fbnx->data[EZ00][j][0]=mp->ez->data[0+j*xres];
                mp->bnd->liao.fbnx->data[HX00][j][0]=mp->hx->data[0+j*xres];
                mp->bnd->liao.fbnx->data[HY00][j][0]=mp->hy->data[0+j*xres];
                mp->bnd->liao.fbnx->data[HZ00][j][0]=mp->hz->data[0+j*xres];

                mp->bnd->liao.fbnx->data[EX0P][j][0]=mp->ex->data[1+j*xres];
                mp->bnd->liao.fbnx->data[EY0P][j][0]=mp->ey->data[1+j*xres];
                mp->bnd->liao.fbnx->data[EZ0P][j][0]=mp->ez->data[1+j*xres];
                mp->bnd->liao.fbnx->data[HX0P][j][0]=mp->hx->data[1+j*xres];
                mp->bnd->liao.fbnx->data[HY0P][j][0]=mp->hy->data[1+j*xres];
                mp->bnd->liao.fbnx->data[HZ0P][j][0]=mp->hz->data[1+j*xres];
            }
            if (set->sb.bxn == SV_BOUNDARY_LIAO) {
                mp->bnd->liao.fbnx->data[EXN0][j][0]=mp->ex->data[xres-1+j*xres];
                mp->bnd->liao.fbnx->data[EYN0][j][0]=mp->ey->data[xres-1+j*xres];
                mp->bnd->liao.fbnx->data[EZN0][j][0]=mp->ez->data[xres-1+j*xres];
                mp->bnd->liao.fbnx->data[HXN0][j][0]=mp->hx->data[xres-1+j*xres];
                mp->bnd->liao.fbnx->data[HYN0][j][0]=mp->hy->data[xres-1+j*xres];
                mp->bnd->liao.fbnx->data[HZN0][j][0]=mp->hz->data[xres-1+j*xres];

                mp->bnd->liao.fbnx->data[EXNP][j][0]=mp->ex->data[xres-2+j*xres];
                mp->bnd->liao.fbnx->data[EYNP][j][0]=mp->ey->data[xres-2+j*xres];
                mp->bnd->liao.fbnx->data[EZNP][j][0]=mp->ez->data[xres-2+j*xres];
                mp->bnd->liao.fbnx->data[HXNP][j][0]=mp->hx->data[xres-2+j*xres];
                mp->bnd->liao.fbnx->data[HYNP][j][0]=mp->hy->data[xres-2+j*xres];
                mp->bnd->liao.fbnx->data[HZNP][j][0]=mp->hz->data[xres-2+j*xres];
            }
        }
    }

    if (set->sb.by0 == SV_BOUNDARY_LIAO || set->sb.byn == SV_BOUNDARY_LIAO) {
        for (i=0; i<xres; i++) {
            if (set->sb.by0 == SV_BOUNDARY_LIAO) {

                mp->bnd->liao.fbny->data[EX00][i][0]=mp->ex->data[i];
                mp->bnd->liao.fbny->data[EY00][i][0]=mp->ey->data[i];
                mp->bnd->liao.fbny->data[EZ00][i][0]=mp->ez->data[i];
                mp->bnd->liao.fbny->data[HX00][i][0]=mp->hx->data[i];
                mp->bnd->liao.fbny->data[HY00][i][0]=mp->hy->data[i];
                mp->bnd->liao.fbny->data[HZ00][i][0]=mp->hz->data[i];

                mp->bnd->liao.fbny->data[EX0P][i][0]=mp->ex->data[i+1*xres];
                mp->bnd->liao.fbny->data[EY0P][i][0]=mp->ey->data[i+1*xres];
                mp->bnd->liao.fbny->data[EZ0P][i][0]=mp->ez->data[i+1*xres];
                mp->bnd->liao.fbny->data[HX0P][i][0]=mp->hx->data[i+1*xres];
                mp->bnd->liao.fbny->data[HY0P][i][0]=mp->hy->data[i+1*xres];
                mp->bnd->liao.fbny->data[HZ0P][i][0]=mp->hz->data[i+1*xres];
                 
            }
            if (set->sb.byn == SV_BOUNDARY_LIAO) {
                mp->bnd->liao.fbny->data[EXN0][i][0]=mp->ex->data[i+(yres-1)*xres];
                mp->bnd->liao.fbny->data[EYN0][i][0]=mp->ey->data[i+(yres-1)*xres];
                mp->bnd->liao.fbny->data[EZN0][i][0]=mp->ez->data[i+(yres-1)*xres];
                mp->bnd->liao.fbny->data[HXN0][i][0]=mp->hx->data[i+(yres-1)*xres];
                mp->bnd->liao.fbny->data[HYN0][i][0]=mp->hy->data[i+(yres-1)*xres];
                mp->bnd->liao.fbny->data[HZN0][i][0]=mp->hz->data[i+(yres-1)*xres];

                mp->bnd->liao.fbny->data[EXNP][i][0]=mp->ex->data[i+(yres-2)*xres];
                mp->bnd->liao.fbny->data[EYNP][i][0]=mp->ey->data[i+(yres-2)*xres];
                mp->bnd->liao.fbny->data[EZNP][i][0]=mp->ez->data[i+(yres-2)*xres];
                mp->bnd->liao.fbny->data[HXNP][i][0]=mp->hx->data[i+(yres-2)*xres];
                mp->bnd->liao.fbny->data[HYNP][i][0]=mp->hy->data[i+(yres-2)*xres];
                mp->bnd->liao.fbny->data[HZNP][i][0]=mp->hz->data[i+(yres-2)*xres];
            }
        }
    }


    if (set->sc.verbose) 
        printf("done.\n");

}



/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

