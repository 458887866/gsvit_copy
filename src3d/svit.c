
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
#include "tests.h"
#include <libgwymodule/gwymodule.h>

void
puthelp()
{
    fprintf(stderr, "Gsvit v %s - FDTD solver with GPU support\n", VERSION);
    fprintf(stderr, "Usage: ./gsvit parameter_file for normal computation\n\n");
    fprintf(stderr, "Tests: ./gsvit test 0 for checking available GPUs only\n");
    fprintf(stderr, "       ./gsvit test 1 for making single test on CPU and all GPUs (if there are any)\n");
    fprintf(stderr, "       ./gsvit test 2 for testing key algorithms performance at 100x100x100 voxels\n");
    fprintf(stderr, "       ./gsvit test 3 for testing key algorithms performance at 200x200x200 voxels\n");
    fprintf(stderr, "       ./gsvit test 12 for comparing GPU/CPU time scaling up to 200x200x200 voxels\n");
    fprintf(stderr, "       ./gsvit test 13 for comparing GPU/CPU time scaling up to 300x300x300 voxels\n");
    fprintf(stderr, "       ./gsvit test 14 for comparing GPU/CPU time scaling up to 400x400x400 voxels\n");
    fprintf(stderr, "       ./gsvit test 15 for comparing GPU/CPU time scaling up to 500x500x500 voxels\n");
    fprintf(stderr, "       ./gsvit test 20 for checking GPU time scaling up to 400x400x400 voxels\n");
    fprintf(stderr, "       ./gsvit test 3N for testing multithread speedup at N00xN00xN00 voxels (i.e. 100-900)\n");
    fprintf(stderr, "       ./gsvit test e.g. 10000010 for running a specific test on CPU (see documentation)\n");
    fprintf(stderr, "Report bugs to %s\n", PACKAGE_BUGREPORT);
 
}

/*
static void
load_modules(void)
{
    static const gchar *const module_types[] = { "file", "layer", NULL };
    GPtrArray *module_dirs;
    const gchar *q;
    gchar *p;
    guint i;

    module_dirs = g_ptr_array_new();

    p = gwy_find_self_dir("modules");
    for (i = 0; module_types[i]; i++)
        g_ptr_array_add(module_dirs, g_build_filename(p, module_types[i], NULL));
    g_free(p);

    q = gwy_get_user_dir();
    for (i = 0; module_types[i]; i++)
        g_ptr_array_add(module_dirs, g_build_filename(q, module_types[i], NULL));

    g_ptr_array_add(module_dirs, NULL);
    gwy_module_register_modules((const gchar**)module_dirs->pdata);

    for (i = 0; module_dirs->pdata[i]; i++)
        g_free(module_dirs->pdata[i]);
    g_ptr_array_free(module_dirs, TRUE);
}*/

static gboolean
killcheck(SvSet *set)
{
    FILE *fr;
    gint val;

    if (set->sc.killer) {
        if (g_file_test(set->sc.killer, G_FILE_TEST_EXISTS)) {
            fr = fopen(set->sc.killer, "r");
            if (fr) {
                fscanf(fr, "%d", &val);
                fclose(fr);
                if (val == 1) return TRUE;
            }
        }
    }
    return FALSE;
}

int
main(int argc, char *argv[])
{
    SvSet set;
    SvPool *mp;
    gint i;

#if GLIB_MINOR_VERSION < 32
    g_thread_init(NULL);
#endif
#if GLIB_MINOR_VERSION < 36
    g_type_init();
#endif

    gwy_process_type_init();

    if (argc < 2) {
        puthelp();
        return 1;
    }

    if (strcmp(argv[1], "test") == 0) {
        if (argc < 3) {
            printf("Error: please specify the test level (0-39)\n");
            return 1;
        } else {
            run_tests(atoi(argv[2]));
            return 0;
        }
    }

    if (strcmp(argv[1], "benchmark") == 0) {
        if (argc < 3) {
            printf("Error: please specify the number of cores (0-1000)\n");
            return 1;
        } else {
            run_benchmark(atoi(argv[2]));
            return 0;
        }
    }

    if (strcmp(argv[1], "bigbenchmark") == 0) {
        if (argc < 3) {
            printf("Error: please specify the number of cores (0-1000)\n");
            return 1;
        } else {
            run_bigbenchmark(atoi(argv[2]));
            return 0;
        }
    }

    set.sc.killer = NULL;
    if (argc > 1 && strcmp(argv[1], "--killer") == 0) {
        if (argc < 3) {
            printf("Error: please specify killer file\n");
            return 1;
        } else
            set.sc.killer = g_strdup(argv[2]);
    }
    if (argc > 2 && strcmp(argv[2], "--killer") == 0) {
        if (argc < 4) {
            printf("Error: please specify killer file\n");
            return 1;
        } else
            set.sc.killer = g_strdup(argv[3]);
    }
    if (argc > 3 && strcmp(argv[3], "--killer")==0) {
        if (argc < 5) {
            printf("Error: please specify killer file\n");
            return 1;
        } else
            set.sc.killer = g_strdup(argv[4]);
    }
    //printf("killer set to %s, %d\n", set.sc.killer, argc);

    /*parse settings*/
    clear_settings(&set, FALSE);
    if (!parse_settings(argv[1], &set, TRUE)) {
        printf("Settings error!\n"); 
        puthelp();
        return 1;
    } else
        if (set.sc.verbose>0) 
            printf("Parameter file loaded successfully\n");


    /*initialize data and make plans*/
    mp = init_and_make_plan(&set);
    mp->set = &set;


    if (mp == NULL) 
        return 1;
    
#ifdef UCUDA
    if (set.plan.gpumode == SV_GPUMODE_FULL) {
	    sv_pool_gpu_init(mp, &set);
	    sv_pool_gpu_copyto(mp, &set, SV_SYNC_ALL);
    }
#endif

    /*run computation*/
    for (i = 0; i < set.sc.nsteps; i++) {
        if (set.sc.verbose > 0) {
	        printf("______________ step %d ______________\n", set.sc.step_act);
	        fflush(stdout);
        }

	    if (set.plan.gpumode == SV_GPUMODE_NONE) {	
	        sv_pool_boundary_copy(mp, &set);
	        sv_pool_ystep_h(mp, &set);               //H n+1/2
	        sv_pool_apply_source_hstep(mp, &set);    //H n+1/2 source, tsf applying Es n. before create Es n, Hs n+1/2
	        sv_pool_boundary_hstep(mp, &set);

            sv_pool_subgrid_copy(mp, &set);
            sv_pool_subgrid_step(mp, &set);

            sv_pool_ystep_e(mp, &set);

            sv_pool_boundary_estep(mp, &set);        //E n+1
	        sv_pool_apply_source_estep(mp, &set);    //E n+1 source, tsf applying Hs n+1/2

            sv_pool_getsum(mp, &set);
	        sv_pool_farfield(mp, &set);
	        sv_pool_pfarfield(mp, &set);
	        sv_pool_getpoints(mp, &set);
	        sv_pool_getforce(mp, &set);
	        if (is_output(mp, &set) || i == (set.sc.nsteps-1)) 
	    	    sv_pool_output(mp, &set);
	    } else {
#ifdef UCUDA
	        sv_pool_gpu_boundary_copy(mp, &set);
	        sv_pool_gpu_ystep_h(mp, &set);
	        sv_pool_gpu_apply_source_hstep(mp, &set);  
	        sv_pool_gpu_boundary_hstep(mp, &set);
	        sv_pool_gpu_ystep_e(mp, &set);  
	        sv_pool_gpu_apply_source_estep(mp, &set);
	        sv_pool_gpu_boundary_estep(mp, &set);
	        sv_pool_gpu_getsum(mp, &set);
	        sv_pool_gpu_apply_farfield(mp, &set);
	        sv_pool_gpu_apply_pfarfield(mp, &set);
	        sv_pool_gpu_getpoints(mp, &set);
	        sv_pool_gpu_getforce(mp, &set);

	        if (is_output(mp, &set) || i == (set.sc.nsteps-1)) {
		        sv_pool_gpu_copyfrom(mp, &set, SV_SYNC_EH);
		        sv_pool_output(mp, &set);
	        }
#endif
	    }

	    if (killcheck(&set)) {
	        if (set.sc.verbose)
	    	    printf("GSvit forced to quit by external file\n");
	        goto exit;
	    }
	    set.sc.step_act++;
    }

    printf("GSvit successfully finished.\n");

exit:
    remove_temporary_files(&set);
    sv_pool_free(mp, &set);

    return 0;
}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
