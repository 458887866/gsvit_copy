
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

/*  subgrid.c : 
 *  subgrid implementation
 *
 */

#include "boundary.h"
#include "yee.h"
#include "pool.h"
#include "settings.h"
#include "constants.h"
#include "subgrid.h"
#include <math.h>
#include <omp.h>


SvSubgrid *sv_subgrid_new(SvPool *mp, SvSet *set)
{
    gint i;
    SvSg *sg;


    if (set->sg.nsg==0) return NULL;

    SvSubgrid* subgrid = (SvSubgrid*)g_malloc(sizeof(SvSubgrid));


    if (set->sc.verbose>1) printf("Allocating %d subgrid sets\n", set->sg.nsg);

    subgrid->nsubgrid = set->sg.nsg;

    subgrid->sg = (SvSg **)g_malloc(subgrid->nsubgrid*sizeof(SvSg *));

    for (i=0; i<set->sg.nsg; i++) {

        sg = (SvSg*)g_malloc(sizeof(SvSg));
        sg->ifrom = set->sg.sg[i].box_i0;
        sg->jfrom = set->sg.sg[i].box_j0;
        sg->kfrom = set->sg.sg[i].box_k0;
        sg->ito = set->sg.sg[i].box_in;
        sg->jto = set->sg.sg[i].box_jn;
        sg->kto = set->sg.sg[i].box_kn;

        sg->division = set->sg.sg[i].division;

        sg->d = sv_yee_data_new(set, set->plan.dt/set->sg.sg[i].division);
        sv_yee_data_allocate(sg->d, set, 
                             (sg->ito - sg->ifrom)*sg->division + 1,
                             (sg->jto - sg->jfrom)*sg->division + 1,
                             (sg->kto - sg->kfrom)*sg->division + 1,
                             set->sp.dx/sg->division,
                             set->sp.dy/sg->division,
                             set->sp.dz/sg->division);


        sg->phxres = (sg->ito - sg->ifrom) + 2;
        sg->phyres = (sg->jto - sg->jfrom) + 2;
        sg->phzres = (sg->kto - sg->kfrom) + 2;

        sg->nx0hy = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);
        sg->nx0hz = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);
        sg->nxnhy = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);
        sg->nxnhz = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);

        sg->ny0hx = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);
        sg->ny0hz = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);
        sg->nynhx = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);
        sg->nynhz = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);

        sg->nz0hx = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);
        sg->nz0hy = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);
        sg->nznhx = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);
        sg->nznhy = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);

        sg->ax0hy = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);
        sg->ax0hz = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);
        sg->axnhy = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);
        sg->axnhz = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);

        sg->ay0hx = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);
        sg->ay0hz = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);
        sg->aynhx = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);
        sg->aynhz = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);

        sg->az0hx = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);
        sg->az0hy = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);
        sg->aznhx = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);
        sg->aznhy = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);

        sg->px0hy = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);
        sg->px0hz = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);
        sg->pxnhy = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);
        sg->pxnhz = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);

        sg->py0hx = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);
        sg->py0hz = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);
        sg->pynhx = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);
        sg->pynhz = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);

        sg->pz0hx = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);
        sg->pz0hy = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);
        sg->pznhx = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);
        sg->pznhy = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);

        sg->ppx0hy = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);
        sg->ppx0hz = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);
        sg->ppxnhy = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);
        sg->ppxnhz = gwy_data_field_new(sg->phyres, sg->phzres, sg->phyres, sg->phzres, TRUE);

        sg->ppy0hx = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);
        sg->ppy0hz = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);
        sg->ppynhx = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);
        sg->ppynhz = gwy_data_field_new(sg->phxres, sg->phzres, sg->phxres, sg->phzres, TRUE);

        sg->ppz0hx = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);
        sg->ppz0hy = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);
        sg->ppznhx = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);
        sg->ppznhy = gwy_data_field_new(sg->phxres, sg->phyres, sg->phxres, sg->phyres, TRUE);



        if (set->sc.verbose>1) printf("Allocated subgrid %d : %dx%dx%d\n", i, sg->d->xres, sg->d->yres, sg->d->zres);

        subgrid->sg[i] = sg;
     }

    return subgrid;
}


int sv_pool_subgrid_copy(SvPool *mp, SvSet *set)
{
    gint n;
    SvSg *sg;

    gboolean debugme = FALSE;

    if (set->sc.verbose > 1 && set->sg.nsg>0) printf("subgrid copy\n");

    for (n=0; n<set->sg.nsg; n++) {

        sg = mp->sg->sg[n];

        //i planes
        gwy_data_field_copy(sg->px0hy, sg->ppx0hy, FALSE);
        gwy_data_field_copy(sg->px0hz, sg->ppx0hz, FALSE);
        gwy_data_field_copy(sg->pxnhy, sg->ppxnhy, FALSE);
        gwy_data_field_copy(sg->pxnhz, sg->ppxnhz, FALSE);

        gwy_data_field_copy(sg->ax0hy, sg->px0hy, FALSE);
        gwy_data_field_copy(sg->ax0hz, sg->px0hz, FALSE);
        gwy_data_field_copy(sg->axnhy, sg->pxnhy, FALSE);
        gwy_data_field_copy(sg->axnhz, sg->pxnhz, FALSE);

        if (debugme) printf("getting ax0hy\n");
        sv_dcube_get_datafield_area(mp->d->hy, sg->ax0hy, sg->ifrom, -1, -1, 0, sg->jfrom-1, sg->kfrom-1);
        if (debugme) printf("getting ax0hz\n");
        sv_dcube_get_datafield_area(mp->d->hz, sg->ax0hz, sg->ifrom, -1, -1, 0, sg->jfrom-1, sg->kfrom-1);
        if (debugme) printf("getting axnhy\n");
        sv_dcube_get_datafield_area(mp->d->hy, sg->axnhy, sg->ito, -1, -1, 0, sg->jfrom-1, sg->kfrom-1);
        if (debugme) printf("getting axnhz\n");
        sv_dcube_get_datafield_area(mp->d->hz, sg->axnhz, sg->ito, -1, -1, 0, sg->jfrom-1, sg->kfrom-1);
  
        //j planes
        gwy_data_field_copy(sg->py0hx, sg->ppy0hx, FALSE);
        gwy_data_field_copy(sg->py0hz, sg->ppy0hz, FALSE);
        gwy_data_field_copy(sg->pynhx, sg->ppynhx, FALSE);
        gwy_data_field_copy(sg->pynhz, sg->ppynhz, FALSE);

        gwy_data_field_copy(sg->ay0hx, sg->py0hx, FALSE);
        gwy_data_field_copy(sg->ay0hz, sg->py0hz, FALSE);
        gwy_data_field_copy(sg->aynhx, sg->pynhx, FALSE);
        gwy_data_field_copy(sg->aynhz, sg->pynhz, FALSE);

        if (debugme) printf("getting ay0hx\n");
        sv_dcube_get_datafield_area(mp->d->hx, sg->ay0hx, -1, sg->jfrom, -1, sg->ifrom-1, 0, sg->kfrom-1);
        if (debugme) printf("getting ay0hz\n");
        sv_dcube_get_datafield_area(mp->d->hz, sg->ay0hz, -1, sg->jfrom, -1, sg->ifrom-1, 0, sg->kfrom-1);
        if (debugme) printf("getting aynhx\n");
        sv_dcube_get_datafield_area(mp->d->hx, sg->aynhx, -1, sg->jto, -1, sg->ifrom-1, 0, sg->kfrom-1);
        if (debugme) printf("getting aynhz\n");
        sv_dcube_get_datafield_area(mp->d->hz, sg->aynhz, -1, sg->jto, -1, sg->ifrom-1, 0, sg->kfrom-1);
    
        //k planes
        gwy_data_field_copy(sg->pz0hx, sg->ppz0hx, FALSE);
        gwy_data_field_copy(sg->pz0hy, sg->ppz0hy, FALSE);
        gwy_data_field_copy(sg->pznhx, sg->ppznhx, FALSE);
        gwy_data_field_copy(sg->pznhy, sg->ppznhy, FALSE);

        gwy_data_field_copy(sg->az0hx, sg->pz0hx, FALSE);
        gwy_data_field_copy(sg->az0hy, sg->pz0hy, FALSE);
        gwy_data_field_copy(sg->aznhx, sg->pznhx, FALSE);
        gwy_data_field_copy(sg->aznhy, sg->pznhy, FALSE);



        if (debugme) printf("getting az0hx\n");
        sv_dcube_get_datafield_area(mp->d->hx, sg->az0hx, -1, -1, sg->kfrom, sg->ifrom-1, sg->jfrom-1, 0);
        if (debugme) printf("getting az0hy\n");
        sv_dcube_get_datafield_area(mp->d->hy, sg->az0hy, -1, -1, sg->kfrom, sg->ifrom-1, sg->jfrom-1, 0);
        if (debugme) printf("getting aznhx\n");
        sv_dcube_get_datafield_area(mp->d->hx, sg->aznhx, -1, -1, sg->kto, sg->ifrom-1, sg->jfrom-1, 0);
        if (debugme) printf("getting aznhy\n");
        sv_dcube_get_datafield_area(mp->d->hy, sg->aznhy, -1, -1, sg->kto, sg->ifrom-1, sg->jfrom-1, 0);

//        printf("getting storage sg->az0hx (1, 1) from main field %d %d %d  %g\n", sg->ifrom-1+1, sg->jfrom-1+1, sg->kfrom,
//               gwy_data_field_get_val(sg->az0hx, 1, 1));
//        printf("getting storage sg->az0hy (0, 0, 0) from main field %d %d %d  %g\n", sg->ifrom-1, sg->jfrom-1, sg->kfrom,
//               gwy_data_field_get_val(sg->az0hy, 0, 0));
      
    }
    return 0;

}



gdouble
pint(gdouble main, gdouble prev, gdouble pprev, gdouble frac)
{
    gdouble a, b;
    a = (main - pprev)/2;
    b = main + pprev - 2*prev;
    return prev + a*frac + b*frac*frac/2;
}


gdouble
pext(gdouble main, gdouble prev, gdouble pprev, gdouble frac)
{
    gdouble a, b;
    a = (main - pprev)/2;
    b = main + pprev - 2*prev;
    return prev + a*frac + b*frac*frac/2;
}

void
extrapolate_datafield(GwyDataField *next, GwyDataField *actual, GwyDataField *prev, GwyDataField *pprev, gdouble frac)
{
    gint i;
    gdouble *ndata, *adata, *pdata, *ppdata;

    ndata = gwy_data_field_get_data(next);
    adata = gwy_data_field_get_data(actual);
    pdata = gwy_data_field_get_data(prev);
    ppdata = gwy_data_field_get_data(pprev);
    
    for (i=0; i<(gwy_data_field_get_xres(next)*gwy_data_field_get_yres(next)); i++)
    {
        ndata[i] = pext(adata[i], pdata[i], ppdata[i], frac);
    }

}

void extrapolate(SvSg *sg, gdouble frac)
{


    extrapolate_datafield(sg->nx0hy, sg->ax0hy, sg->px0hy, sg->ppx0hy, frac);
    extrapolate_datafield(sg->nx0hz, sg->ax0hz, sg->px0hz, sg->ppx0hz, frac);
    extrapolate_datafield(sg->nxnhy, sg->axnhy, sg->pxnhy, sg->ppxnhy, frac);
    extrapolate_datafield(sg->nxnhz, sg->axnhz, sg->pxnhz, sg->ppxnhz, frac);

    extrapolate_datafield(sg->ny0hx, sg->ay0hx, sg->py0hx, sg->ppy0hx, frac);
    extrapolate_datafield(sg->ny0hz, sg->ay0hz, sg->py0hz, sg->ppy0hz, frac);
    extrapolate_datafield(sg->nynhx, sg->aynhx, sg->pynhx, sg->ppynhx, frac);
    extrapolate_datafield(sg->nynhz, sg->aynhz, sg->pynhz, sg->ppynhz, frac);

    extrapolate_datafield(sg->nz0hx, sg->az0hx, sg->pz0hx, sg->ppz0hx, frac);
    extrapolate_datafield(sg->nz0hy, sg->az0hy, sg->pz0hy, sg->ppz0hy, frac);
    extrapolate_datafield(sg->nznhx, sg->aznhx, sg->pznhx, sg->ppznhx, frac);
    extrapolate_datafield(sg->nznhy, sg->aznhy, sg->pznhy, sg->ppznhy, frac);

    //printf("extrapolate storage hx[1,1]: %g %g %g  %g\n",  gwy_data_field_get_val(sg->ppz0hx, 1, 1), gwy_data_field_get_val(sg->pz0hx, 1, 1), gwy_data_field_get_val(sg->az0hx, 1, 1), gwy_data_field_get_val(sg->nz0hx, 1, 1));

}

void interp(SvSg *sg, SvPool *mp, SvSet *set)
{
    gint i, j, k; //local field positions
    gdouble mdposx, mdposy, mdposz;

    gdouble factor = ceil(sg->division/2.0)/sg->division;

    gint xres = sg->d->xres;
    gint yres = sg->d->yres;
    gint zres = sg->d->zres;


    gboolean debugme = FALSE;

    //i planes
    for (j=0; j<yres; j++)
    {
       for (k=0; k<zres; k++) 
       {
            mdposy = ((gdouble)j/sg->division) + factor + 0.5;
            mdposz = ((gdouble)k/sg->division) + 1 + 0.5;

            if (debugme) printf("field iplane hy (0 %d %d) from  nx0hy %g %g\n", j, k, mdposy, mdposz);

            sg->d->hy->data[0][j][k]      = gwy_data_field_get_dval(sg->nx0hy, mdposy, mdposz, GWY_INTERPOLATION_BILINEAR);
            sg->d->hy->data[xres-1][j][k] = gwy_data_field_get_dval(sg->nxnhy, mdposy, mdposz, GWY_INTERPOLATION_BILINEAR);

            mdposy = ((gdouble)j/sg->division) + 1 + 0.5;
            mdposz = ((gdouble)k/sg->division) + factor + 0.5;

            if (debugme) printf("field iplane hz (0 %d %d) from  nx0hz %g %g\n", j, k, mdposy, mdposz);


            sg->d->hz->data[0][j][k]      = gwy_data_field_get_dval(sg->nx0hz, mdposy, mdposz, GWY_INTERPOLATION_BILINEAR);
            sg->d->hz->data[xres-1][j][k] = gwy_data_field_get_dval(sg->nxnhz, mdposy, mdposz, GWY_INTERPOLATION_BILINEAR);
        }
    }
    
    //j planes
    for (i=0; i<xres; i++)
    {
       for (k=0; k<zres; k++) 
       {
            mdposx = ((gdouble)i/sg->division) + factor + 0.5;
            mdposz = ((gdouble)k/sg->division) + 1 + 0.5;

            if (debugme) printf("field jplane hx (%d 0 %d) from  ny0hx %g %g\n", i, k, mdposx, mdposz);


            sg->d->hx->data[i][0][k]      = gwy_data_field_get_dval(sg->ny0hx, mdposx, mdposz, GWY_INTERPOLATION_BILINEAR);
            sg->d->hx->data[i][yres-1][k] = gwy_data_field_get_dval(sg->nynhx, mdposx, mdposz, GWY_INTERPOLATION_BILINEAR);

            mdposx = ((gdouble)i/sg->division) + 1 + 0.5;
            mdposz = ((gdouble)k/sg->division) + factor + 0.5;

            if (debugme) printf("field jplane hz (%d 0 %d) from  ny0hz %g %g\n", i, k, mdposx, mdposz);


            sg->d->hz->data[i][0][k]      = gwy_data_field_get_dval(sg->ny0hz, mdposx, mdposz, GWY_INTERPOLATION_BILINEAR);
            sg->d->hz->data[i][yres-1][k] = gwy_data_field_get_dval(sg->nynhz, mdposx, mdposz, GWY_INTERPOLATION_BILINEAR);
        }
    }
    
    //k planes
    for (i=0; i<xres; i++)
    {
       for (j=0; j<yres; j++)
       {
            mdposx = ((gdouble)i/sg->division) + factor + 0.5;
            mdposy = ((gdouble)j/sg->division) + 1 + 0.5;

            
            if (debugme && i<2 && j<2) printf("sg field kplane hx (%d %d 0) from  storage nz0hx at %g %g  set to %g     nz0hx[1,1] is %g  wtf %g\n", i, j, mdposx, mdposy, gwy_data_field_get_dval(sg->nz0hx, mdposx, mdposy, GWY_INTERPOLATION_BILINEAR), 
                                   gwy_data_field_get_val(sg->nz0hx, 1, 1), gwy_data_field_get_dval(sg->nz0hx, 1+0.5, 1+0.5, GWY_INTERPOLATION_BILINEAR));

            sg->d->hx->data[i][j][0]      = gwy_data_field_get_dval(sg->nz0hx, mdposx, mdposy, GWY_INTERPOLATION_BILINEAR);
            sg->d->hx->data[i][j][zres-1] = gwy_data_field_get_dval(sg->nznhx, mdposx, mdposy, GWY_INTERPOLATION_BILINEAR);

            mdposx = ((gdouble)i/sg->division) + 1 + 0.5;
            mdposy = ((gdouble)j/sg->division) + factor + 0.5;

            if (debugme) printf("field kplane hy (%d %d 0) from  nz0hy %g %g\n", i, j, mdposx, mdposy);


            sg->d->hy->data[i][j][0]      = gwy_data_field_get_dval(sg->nz0hy, mdposx, mdposy, GWY_INTERPOLATION_BILINEAR);
            sg->d->hy->data[i][j][zres-1] = gwy_data_field_get_dval(sg->nznhy, mdposx, mdposy, GWY_INTERPOLATION_BILINEAR);
       }
    }

}

void hback(SvSg *sg, SvPool *mp, SvSet *set)
{
    gint i, j, k, isg, jsg, ksg;

    for (i=sg->ifrom; i<sg->ito; i++) 
    {
       for (j=sg->jfrom; j<sg->ito; j++) 
       {
           for (k=sg->kfrom; k<sg->kto; k++) 
           {
               isg = (i-sg->ifrom)*sg->division;
               jsg = (j-sg->jfrom)*sg->division;
               ksg = (k-sg->kfrom)*sg->division;

               if (isg>=sg->d->xres) continue;
               if (jsg>=sg->d->yres) continue;
               if (ksg>=sg->d->zres) continue;


               if (isg>1 && isg<(sg->d->xres-1)) mp->d->hx->data[i][j][k] = sg->d->hx->data[isg+1][jsg][ksg];
               if (jsg>1 && jsg<(sg->d->yres-1)) mp->d->hy->data[i][j][k] = sg->d->hy->data[isg][jsg+1][ksg];
               if (ksg>1 && ksg<(sg->d->zres-1)) mp->d->hz->data[i][j][k] = sg->d->hz->data[isg][jsg][ksg+1];
            }
       }  
    }
}

void hweight(SvSg *sg, SvPool *mp, SvSet *set)
{
    gint i, j, k;

    gint xres = sg->d->xres;
    gint yres = sg->d->yres;
    gint zres = sg->d->zres;

    //i planes
    for (j=0; j<yres; j++)
    {
        for (k=0; k<zres; k++)
        {
            sg->d->hy->data[1][j][k] = 0.95*sg->d->hy->data[1][j][k] + 0.05*(sg->d->hy->data[0][j][k] + sg->d->hy->data[2][j][k])/2.0;
            sg->d->hz->data[1][j][k] = 0.95*sg->d->hz->data[1][j][k] + 0.05*(sg->d->hz->data[0][j][k] + sg->d->hz->data[2][j][k])/2.0;
            sg->d->hy->data[xres-2][j][k] = 0.95*sg->d->hy->data[xres-2][j][k] + 0.05*(sg->d->hy->data[xres-1][j][k] + sg->d->hy->data[xres-3][j][k])/2.0;
            sg->d->hz->data[xres-2][j][k] = 0.95*sg->d->hz->data[xres-2][j][k] + 0.05*(sg->d->hz->data[xres-1][j][k] + sg->d->hz->data[xres-3][j][k])/2.0;
        }
    }

    for (i=0; i<xres; i++)
    {
       for (k=0; k<zres; k++)
       {
           sg->d->hx->data[i][1][k] = 0.95*sg->d->hx->data[i][1][k] + 0.05*(sg->d->hx->data[i][0][k] + sg->d->hx->data[i][2][k])/2.0;
           sg->d->hz->data[i][1][k] = 0.95*sg->d->hz->data[i][1][k] + 0.05*(sg->d->hz->data[i][0][k] + sg->d->hz->data[i][2][k])/2.0;
           sg->d->hx->data[i][yres-2][k] = 0.95*sg->d->hx->data[i][yres-2][k] + 0.05*(sg->d->hx->data[i][yres-1][k] + sg->d->hx->data[i][yres-3][k])/2.0;
           sg->d->hz->data[i][yres-2][k] = 0.95*sg->d->hz->data[i][yres-2][k] + 0.05*(sg->d->hz->data[i][yres-1][k] + sg->d->hz->data[i][yres-3][k])/2.0;
       }
    }

    for (i=0; i<xres; i++)
    {
       for (j=0; j<yres; j++)
       {
           sg->d->hx->data[i][j][1] = 0.95*sg->d->hx->data[i][j][1] + 0.05*(sg->d->hx->data[i][j][0] + sg->d->hx->data[i][j][2])/2.0;
           sg->d->hy->data[i][j][1] = 0.95*sg->d->hy->data[i][j][1] + 0.05*(sg->d->hy->data[i][j][0] + sg->d->hy->data[i][j][2])/2.0;
           sg->d->hx->data[i][j][zres-2] = 0.95*sg->d->hx->data[i][j][zres-2] + 0.05*(sg->d->hx->data[i][j][zres-1] + sg->d->hx->data[i][j][zres-3])/2.0;
           sg->d->hy->data[i][j][zres-2] = 0.95*sg->d->hy->data[i][j][zres-2] + 0.05*(sg->d->hy->data[i][j][zres-1] + sg->d->hy->data[i][j][zres-3])/2.0;
         }
    }
}

void weight(gdouble *global, gdouble *local)
{
    gdouble val = *global;
    gdouble sval = *local;

    *global = 0.8*val + 0.2*sval;
    *local = 0.8*sval + 0.2*val;
}

void eweight(SvSg *sg, SvPool *mp, SvSet *set)
{
    gint i, j, k, in, jn, kn, sxres, syres, szres;

    sxres = sg->ito - sg->ifrom;
    syres = sg->jto - sg->jfrom;
    szres = sg->kto - sg->kfrom;

    for (i=sg->ifrom; i<=sg->ito; i++) 
    {
       for (j=sg->jfrom; j<=sg->ito; j++) 
       {
           for (k=sg->kfrom; k<=sg->kto; k++) 
           {
               in = i-sg->ifrom;
               jn = j-sg->jfrom;
               kn = k-sg->kfrom;

               //i planes
               if (in==1 && jn>0 && kn>0 && jn<(syres-1) && kn<(szres-1))
               {
                   weight(&(mp->d->ey->data[i][j][k]), &(sg->d->ey->data[in*sg->division][jn*sg->division+1][kn*sg->division]));
                   weight(&(mp->d->ez->data[i][j][k]), &(sg->d->ez->data[in*sg->division][jn*sg->division][kn*sg->division+1]));

               }
               if (in==1 && jn>=0 && kn>=0 && jn<(syres) && kn<(szres))
               {
                   weight(&(mp->d->ex->data[i][j][k]), &(sg->d->ex->data[in*sg->division][jn*sg->division+1][kn*sg->division+1]));
               }

               //j planes
               if (jn==1 && in>0 && kn>0 && in<(sxres-1) && kn<(szres-1))
               {
                   weight(&(mp->d->ex->data[i][j][k]), &(sg->d->ex->data[in*sg->division+1][jn*sg->division][kn*sg->division]));
                   weight(&(mp->d->ez->data[i][j][k]), &(sg->d->ez->data[in*sg->division][jn*sg->division][kn*sg->division+1]));

               }
               if (jn==1 && in>=0 && kn>=0 && in<(sxres) && kn<(szres))
               {
                   weight(&(mp->d->ey->data[i][j][k]), &(sg->d->ey->data[in*sg->division+1][jn*sg->division][kn*sg->division+1]));
               }


               //k planes
               if (kn==1 && in>0 && jn>0 && in<(sxres-1) && jn<(syres-1))
               {

                   weight(&(mp->d->ex->data[i][j][k]), &(sg->d->ex->data[in*sg->division+1][jn*sg->division][kn*sg->division]));
                   weight(&(mp->d->ey->data[i][j][k]), &(sg->d->ey->data[in*sg->division][jn*sg->division+1][kn*sg->division]));

               }
               if (kn==1 && in>=0 && jn>=0 && in<(sxres) && jn<(syres))
               {
                   weight(&(mp->d->ez->data[i][j][k]), &(sg->d->ez->data[in*sg->division+1][jn*sg->division+1][kn*sg->division]));
               }

            }
       }  
    }
}


int sv_pool_subgrid_step(SvPool *mp, SvSet *set)
{
    gint n, m;
    SvSg *sg;

    if (set->sc.verbose > 1 && set->sg.nsg>0) printf("subgrid step\n");

    for (n=0; n<set->sg.nsg; n++) {

        sg = mp->sg->sg[n];

        for (m=1; m<=(2*sg->division); m++) 
        {
            if (m%2 != 0) {
                sv_yee_data_ystep_e(sg->d, NULL, NULL, NULL, set, mp->mats, mp->nmat);

            }

            if (m%2 != 0 && m==sg->division) {
                eweight(sg, mp, set);
             }


            if (m%2 == 0) {
                sv_yee_data_ystep_h(sg->d, NULL, NULL, set, mp->mats, mp->nmat);

           }
            
            if (m%2 == 0) {
                extrapolate(sg, (gdouble)m/(gdouble)(2*sg->division));
                interp(sg, mp, set);
                hweight(sg, mp, set);
            }
            

        }
        hback(sg, mp, set);

    }

    return 0;
}



/*

____________________

   div 3

   e    1/6
   h    2/6
   ziskat interpolovane h z H  2/6

   e    3/6
   vymenit male e 3/6  s velkym E 1/2
   h    4/6
   ziskat interpolovane h z H  4/6

   e    5/6
   h    6/6
   ziskat interpolovane h z H  6/6

   vratit h do H

_____________________

   div 5

   e    1/10
   h    2/10
   ziskat interpolovane h z H  2/10

   e    3/10
   h    4/10
   ziskat interpolovane h z H  4/10

   e    5/10
   vymenit male e 5/10 s velkym E 1/2
   h    6/10
   ziskat interpolovane h z H  6/10

   e    7/10
   h    8/10
   ziskat interpolovane h z H  8/10

   e    9/10
   h    10/10
   ziskat interpolovane h z H  10/10


   vratit h do H

_______________________




   */





/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
