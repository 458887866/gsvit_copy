
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


 /*  pool.c :
  *  main algorithms (Yee solver, including PRLC, PEC, parts of TSF and CPML)
  */

#if defined (_WIN32) && defined (_MSC_VER)
#include "../config_win.h"
#else
#include "../config.h"
#endif
#include <glib.h>
#include <stdlib.h>
#include "yee.h"
#include "pool.h"
#include "constants.h"
#include "source.h"
#include <math.h>
#include <omp.h>
#include "subgrid.h"

SvPool*
sv_pool_new(SvSet *set)
{
    SvPool *mp = (SvPool *)g_malloc(sizeof(SvPool));

    /*alloc main fields*/
    mp->d = sv_yee_data_new(set, set->plan.dt);

    /*alloc other structures*/
    mp->src = sv_source_new();
    mp->out = sv_output_new();

    mp->nmat = 0;
    mp->mats = NULL;

    if (set->sc.verbose > 1)
        printf("Pool initialized\n");

    return mp;
}

void
sv_pool_allocate(SvSet *set, SvPool *mp)
{
    sv_yee_data_allocate(mp->d, set, set->sp.xres, set->sp.yres, set->sp.zres, set->sp.dx, set->sp.dy, set->sp.dz);

    mp->bnd = sv_boundary_new(mp, set);
    mp->farfield = sv_farfield_new(mp, set);
    mp->pfarfield = sv_pfarfield_new(mp, set);
    mp->sg = sv_subgrid_new(mp, set);

    return;
}

void
sv_pool_allocate_gpus(SvPool *mp, SvSet *set)
{
#ifdef UCUDA
    mp->gpu = sv_gpu_new(mp, set);
#endif
}

void
sv_pool_free(SvPool *mp, SvSet *set)
{
    sv_yee_data_free(mp->d, set);
}

/*Yee step: calculate H from E*/
int
sv_pool_ystep_h(SvPool *mp, SvSet *set)
{
    return sv_yee_data_ystep_h(mp->d, &(set->sb), mp->bnd, set, mp->mats, mp->nmat);
}

/*Yee step: calculate E from H*/
gint
sv_pool_ystep_e(SvPool *mp, SvSet *set)
{
    return sv_yee_data_ystep_e(mp->d, &(set->sb), mp->bnd, mp->src, set, mp->mats, mp->nmat);
}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
