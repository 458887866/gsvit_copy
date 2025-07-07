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


/*  svit.c : 
 *  main file calling all the algorithms, including main() function
 */

#include <stdio.h>
#include <glib.h>
#include <stdlib.h>
#include <string.h>
#if defined (_WIN32) && defined (_MSC_VER)
#include "../config_win.h"
#else
#include "../config.h"
#include <sys/time.h>
#endif
#include "settings.h"
#include "pool.h"
#include "plan.h"
//#include "tests.h"
#include <libgwymodule/gwymodule.h>


void puthelp()
{
    fprintf(stderr, "Gsvit2d v %s - 2D FDTD solver with GPU support\n", VERSION);
    fprintf(stderr, "Usage: ./gsvit2d parameter_file for normal computation\n\n");
    fprintf(stderr, "Tests: ./gsvit2d test 0 for checking available GPUs only\n");
    fprintf(stderr, "       ./gsvit2d test 1 for making single test on CPU and all GPUs (if there are any)\n");
    fprintf(stderr, "       ./gsvit2d test 2 for testing key algorithms performance at 100x100 voxels\n");
    fprintf(stderr, "       ./gsvit2d test 3 for testing key algorithms performance at 200x200 voxels\n");
    fprintf(stderr, "       ./gsvit2d test 12 for comparing GPU/CPU time scaling up to 200x200 voxels\n");
    fprintf(stderr, "       ./gsvit2d test 13 for comparing GPU/CPU time scaling up to 300x300 voxels\n");
    fprintf(stderr, "       ./gsvit2d test 14 for comparing GPU/CPU time scaling up to 400x400 voxels\n");
    fprintf(stderr, "       ./gsvit2d test 15 for comparing GPU/CPU time scaling up to 500x500 voxels\n");
    fprintf(stderr, "       ./gsvit2d test 20 for checking GPU time scaling up to 400x400 voxels\n");
    fprintf(stderr, "       ./gsvit2d test 3N for testing multithread speedup at N00xN00 voxels (i.e. 100-900)\n");
    fprintf(stderr, "       ./gsvit2d test e.g. 10000010 for running a specific test on CPU (see documentation)\n");
    fprintf(stderr, "Report bugs to %s\n", PACKAGE_BUGREPORT);

}


int main(int argc, char *argv[])
{
    SvSet set;
    SvPool *mp;
    gint i;
#ifndef G_OS_WIN32
    struct timeval time;
#endif
    GRand *rnd;

#if GLIB_MINOR_VERSION < 32
    g_thread_init(NULL);
#endif
#if GLIB_MINOR_VERSION < 36
    g_type_init();
#endif

    gwy_process_type_init();

    if (argc<2)
    {
        puthelp();
        return 0;
    }

#ifdef G_OS_WIN32
    rnd = g_rand_new();
    g_rand_set_seed(rnd, g_random_int() & 0x7fffffff);
#else
    gettimeofday(&time,NULL);
    rnd = g_rand_new_with_seed((time.tv_sec * 1000) + (time.tv_usec / 1000));
#endif
    set.sc.suffix = g_rand_int_range(rnd, 0, 999999);

    if (strcmp(argv[1], "test")==0) {
        if (argc<3)
        {
            printf("Error: please specify the test level (0-39)\n");
            return 1;
        }
        else {
            //run_tests(atoi(argv[2]));
            return 0;
        }
    }


    /*parse settings*/
    clear_settings(&set);
    if (parse_settings(argv[1], &set)) {
        printf("Settings error!\n"); 
        puthelp();
        return 1;
    } else {
        if (set.sc.verbose>0) printf("Parameter file loaded successfully\n");
    }
    

    /*initialize data and make plans*/
    mp = init_and_make_plan(&set);
    mp->set = &set;

#ifdef UCUDA
    if (set.plan.gpumode == SV_GPUMODE_FULL)
    {
	sv_pool_gpu_init(mp, &set);
	sv_pool_gpu_copyto(mp, &set, SV_SYNC_ALL);
    }
#endif

    /*run computation*/
    for (i=0; i<set.sc.nsteps; i++)
    {
        if (set.sc.verbose > 0) {
                printf("______________ step %d ______________\n", set.sc.step_act);
                fflush(stdout);
        }


	if (set.plan.gpumode == SV_GPUMODE_NONE)
	{	
		sv_pool_boundary_copy(mp, &set);
		sv_pool_ystep_h(mp, &set);
		sv_pool_apply_source_hstep(mp, &set);
		sv_pool_boundary_hstep(mp, &set);
		sv_pool_ystep_e(mp, &set);
		sv_pool_apply_source_estep(mp, &set);

		sv_pool_farfield(mp, &set);
		sv_pool_boundary_estep(mp, &set);
		sv_pool_getsum(mp, &set);
		if (is_output(mp, &set)) 
			sv_pool_output(mp, &set);
	} else {
#ifdef UCUDA
		sv_pool_gpu_boundary_copy(mp, &set);
		sv_pool_gpu_ystep_h(mp, &set);  
		sv_pool_gpu_apply_source_hstep(mp, &set); // uz toto vede k zahade
		sv_pool_gpu_boundary_hstep(mp, &set);
		sv_pool_gpu_ystep_e(mp, &set);  
		sv_pool_gpu_apply_source_estep(mp, &set); // uz toto vede k zahade
		sv_pool_gpu_boundary_estep(mp, &set);
                sv_pool_gpu_getsum(mp, &set);
		
		if (is_output(mp, &set) || i==(set.sc.nsteps-1)) {
		   sv_pool_gpu_copyfrom(mp, &set, SV_SYNC_EH);
		   sv_pool_output(mp, &set);
		}
#endif
	}

        set.sc.step_act++;
    }
    remove_temporary_files(&set);

    return 0;

}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

