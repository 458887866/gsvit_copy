
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

/*  tests.c : 
 *  various tests and benchmarks
 */

#if defined (_WIN32) && defined (_MSC_VER)
#include "../config_win.h"
#else
#include "../config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include "tests.h"
#include "settings.h"
#include <math.h>
#include <time.h>

#ifdef UCUDA
#include <cuda_runtime_api.h>
#endif

#ifdef G_OS_WIN32
#include <windows.h>
//#define snprintf _snprintf
#else
#include <unistd.h>
#endif

typedef enum {
    TEST_SOURCE_POINT = 0,
    TEST_SOURCE_TSF = 1,
    TEST_SOURCE_SF = 2,
    TEST_SOURCE_FOCUSED = 3,
    TEST_SOURCE_LTSF = 4,
    TEST_SOURCE_LTSFF = 5 
} TestSourceType;

typedef enum {
    TEST_MATCHECK_NO = 0,
    TEST_MATCHECK_YES = 1
} TestMatcheckType;

typedef enum {
    TEST_FF_NO = 0,
    TEST_FF_YES = 1,
    TEST_FF_PFF = 2
} TestFarfieldType;

typedef enum {
    TEST_MATERIAL_NONE = 0,
    TEST_MATERIAL_ELECTRIC = 1,
    TEST_MATERIAL_MAGNETIC = 2,
    TEST_MATERIAL_LINTAB_ELECTRIC = 3,
    TEST_MATERIAL_LINTAB_MAGNETIC = 4,
    TEST_MATERIAL_LINTAB_PEC = 5,
    TEST_MATERIAL_LINTAB_DRUDE = 6,
    TEST_MATERIAL_LINTAB_CP = 7,
    TEST_MATERIAL_LINTAB_ADE = 8,
    TEST_MATERIAL_LINTAB_PLRC = 9,
} TestMaterialType;
    
static gint
get_n_cores() {
#ifdef G_OS_WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

gdouble
make_test(gint res, SvBoundaryType bnd, TestMaterialType mat, TestMatcheckType matcheck, TestFarfieldType farfield, TestSourceType source,
	  gboolean gpu, gint ugpu, gdouble *diff, gint nsteps, gint threads)
{
    gchar id[10];
    gchar *parname = "tmp_selftest.par";
    gchar *sourcename = "tmp_selftest.source";
    gchar *materialname = "tmp_selftest.material";
    gchar *testpath;
    gdouble val, damp, size, sum;
    gint pos1, pos2, i;
    gint returnstatus, skip;
    gfloat val1, val1a, val1b, val1c, val2, val2a, val2b, val2c;
    GError *error = NULL;
    GTimeVal starttv, endtv;
    gchar *outlog = NULL;
    gchar *errlog = NULL;
    gchar **argv;
    FILE *fwpar, *fwsource, *fwmaterial, *test_out, *std_out;

    argv = (gchar **)g_malloc(3*sizeof(gchar *));
    for (i = 0; i < 3; i++)
	    argv[i] = (gchar *)g_malloc(1024*sizeof(gchar));
    
#ifdef G_OS_WIN32
    sprintf(argv[0], "%s\\bin\\gsvit3d.exe", get_self_dir());
    //snprintf(argv[0], 1024*sizeof(gchar), "%s", get_self_dir());
#else
    snprintf(argv[0], 1024*sizeof(gchar), "%s", get_self_dir());
#endif
    
    sprintf(argv[1], "tmp_selftest.par");
    argv[2] = NULL;

    g_snprintf(id, sizeof(id), "%d%d%d%d%d%d", res, bnd, mat, matcheck, farfield, source);
    
#ifdef G_OS_WIN32
    testpath = g_build_path("\\", get_self_dir(), "share", "gsvit", "tests", "selftests", id, NULL);
#endif

#ifndef G_OS_WIN32
    testpath = g_build_path("/", g_path_get_dirname(g_path_get_dirname(get_self_dir())), "share", "gsvit", "tests", "selftests", id, NULL);
#endif

    if (gpu)
	    printf("GPU test: %d x %d x %d, boundary %d, material %d, matcheck %d, farfield %d, source %d (id: %s)\n", res, res, res, bnd, mat, matcheck, farfield, source, id);
    else
	    printf("CPU test: %d x %d x %d, boundary %d, material %d, matcheck %d, farfield %d, source %d (id: %s)\n", res, res, res, bnd, mat, matcheck, farfield, source, id);

    /*write basic options*/
    fwpar = fopen(parname, "w");
    if (fwpar == NULL) {
	    fprintf(stderr, "Error: cannot open file %s for writing\n", parname);
	    return -1;
    }

    if (mat == TEST_MATERIAL_LINTAB_ADE || mat == TEST_MATERIAL_LINTAB_PLRC ||  mat == TEST_MATERIAL_LINTAB_DRUDE ||  mat == TEST_MATERIAL_LINTAB_CP) 
	    fprintf(fwpar, "POOL\n%d %d %d 1e-9 1e-9 1e-9\n\n", res, res, res);
    else
	    fprintf(fwpar, "POOL\n%d %d %d 1e-6 1e-6 1e-6\n\n", res, res, res);
    fprintf(fwpar, "COMP\n%d\n\n", nsteps);

    if (mat > 0)
        if (mat == TEST_MATERIAL_LINTAB_ADE || mat == TEST_MATERIAL_LINTAB_PLRC ||  mat == TEST_MATERIAL_LINTAB_DRUDE)
	        fprintf(fwpar, "DT_MULT\n0.5\n\n");

    if (matcheck == TEST_MATCHECK_YES)
	    fprintf(fwpar, "MATMODE_CHECK\n1\n\n");
    else
	    fprintf(fwpar, "MATMODE_CHECK\n0\n\n");

    fprintf(fwpar, "THREADS\n%d\n\n", threads);

    if (mat > 0 || source == TEST_SOURCE_LTSF || source == TEST_SOURCE_LTSFF) {
	    fprintf(fwpar, "MEDIUM_VECTOR\n%s\n\n", materialname);
	    fwmaterial = fopen(materialname, "w");
	    if (mat == TEST_MATERIAL_ELECTRIC)
	        fprintf(fwmaterial, "8 %d %d %d %d %d %d 0 2 1 0 0\n", res/2+10, 30, 30, res-30, res-30, res-30);
	    else if (mat == TEST_MATERIAL_MAGNETIC)
	        fprintf(fwmaterial, "8 %d %d %d %d %d %d 0 1 2 0 0\n", res/2+10, 30, 30, res-30, res-30, res-30);
	    else if (mat == TEST_MATERIAL_LINTAB_ELECTRIC)
	        fprintf(fwmaterial, "8 %d %d %d %d %d %d 1 2 1 0 0\n", res/2+10, 30, 30, res-30, res-30, res-30);
	    else if (mat == TEST_MATERIAL_LINTAB_MAGNETIC)
	        fprintf(fwmaterial, "8 %d %d %d %d %d %d 1 1 2 0 0\n", res/2+10, 30, 30, res-30, res-30, res-30);
	    else if (mat == TEST_MATERIAL_LINTAB_PEC)
	        fprintf(fwmaterial, "8 %d %d %d %d %d %d 10\n", res/2+10, 30, 30, res-30, res-30, res-30);
	    else if (mat == TEST_MATERIAL_LINTAB_DRUDE)
	    	fprintf(fwmaterial, "8 %d %d %d %d %d %d 2 5.0 9.5 0.0987\n", res/2+10, 30, 30, res-30, res-30, res-30);
	    else if (mat == TEST_MATERIAL_LINTAB_CP)
	        fprintf(fwmaterial, "8 %d %d %d %d %d %d 4  1.1431 1.3202e16 1.0805e14  0.26698 -1.2371 3.8711e15 4.4642e14  3.0834 -1.0968 4.1684e15 2.3555e15\n", 
	    	    res/2+10, 30, 30, res-30, res-30, res-30);
	    else if (mat == TEST_MATERIAL_LINTAB_ADE)
	        fprintf(fwmaterial, "8 %d %d %d %d %d %d 5  1.1431 0 1.3202e16 1.0805e14  0.26698 -1.2371 3.8711e15 4.4642e14  3.0834 -1.0968 4.1684e15 2.3555e15\n", 
	    	    res/2+10, 30, 30, res-30, res-30, res-30);
	    else if (mat == TEST_MATERIAL_LINTAB_PLRC)
	        fprintf(fwmaterial, "8 %d %d %d %d %d %d 6  1.1431 0 1.3202e16 1.0805e14  0.26698 -1.2371 3.8711e15 4.4642e14  3.0834 -1.0968 4.1684e15 2.3555e15\n", 
	    	    res/2+10, 30, 30, res-30, res-30, res-30);
            else if (source == TEST_SOURCE_LTSF || source == TEST_SOURCE_LTSFF) 
	            fprintf(fwmaterial, "8 %d %d %d %d %d %d 1 4 1 0 0\n", 0, 0, res/2-15, res, res, res/2+15);
		
	    fclose(fwmaterial);
    }

    skip = 100;
    if (farfield == TEST_FF_NO) {
        if (mat == TEST_MATERIAL_LINTAB_CP || mat == TEST_MATERIAL_LINTAB_PLRC || mat == TEST_MATERIAL_LINTAB_ADE ||  mat == TEST_MATERIAL_LINTAB_DRUDE)
	        fprintf(fwpar, "OUT_POINT\nEx %d %d %d %d %s\n\n", skip, res/2, res/2+30, res/2+30, id);
        else
	        fprintf(fwpar, "OUT_POINT\nEx %d %d %d %d %s\n\n", skip, res/2+30, res/2+30, res/2+30, id);
    } else if (farfield == TEST_FF_YES)
	    fprintf(fwpar, "NFFF\n%d %d %d %d %d %d 1 10000000\n\nNFFF_RAMAHI_POINT\n%d %d %d %s\n\n", 16, 16, 16, res-16, res-16, res-16, res-12, res-12, res-12, id);
    else {
	    fprintf(fwpar, "PERIODIC_NFFF\n%d %d %d %d %d %d -1 -1 2 2\n\nPERIODIC_NFFF_RAMAHI_POINT\n%d %d %d %s\n\n", 26, 26, 26, res-26, res-26, res-26, res-12, res-12, res-12, id);
	    fprintf(fwpar, "MBOUNDARY_X0\nperiodic %d\n\nMBOUNDARY_Y0\nperiodic %d\n\nMBOUNDARY_XN\nperiodic %d\n\nMBOUNDARY_YN\nperiodic %d\n\n", 26, 26, res-26, res-26); 
    }

    if (gpu)
        fprintf(fwpar, "GPU\n1\n\nUGPU\n%d\n\n", ugpu);

    /*for debug only*/
    if (0) {
	    fprintf(fwpar, "OUT_IMAGE\nEx 20 %d -1 -1 exx logscale\n\n", res/2);
	    fprintf(fwpar, "OUT_IMAGE\nEx 20 -1 %d -1 exy logscale\n\n", res/2);
	    fprintf(fwpar, "OUT_IMAGE\nEx 20 -1 -1 %d exz logscale\n\n", res/2);

	    fprintf(fwpar, "OUT_FILE\n%s.gwy\n\n", id);

	    fprintf(fwpar, "VERBOSE\n3\n\n");
    }

    /*sources*/
    if (source == TEST_SOURCE_POINT) {
	    fwsource = fopen(sourcename, "w");
	    size = 10;

	    fprintf(fwsource, "%d\n", nsteps);
	    for (i = 0; i < nsteps; i++) {
	        val = 5 * sin(0.25*(double)i);
	        //if (i<100) val *= 0.01*((double)(i));
	        damp = ((double)i*1e-6 - 4*size) / sqrt((gdouble)2.0) / size;
	        val *= exp(-damp*damp);
	        fprintf(fwsource, "%d %g %g %g %g %g %g\n", i, val, 0.0, 0.0, 0.0, 0.0, 0.0);
	    }
	    fclose(fwsource);
	    fprintf(fwpar, "SOURCE_POINT\n%d %d %d 0 %s\n\n", res/2, res/2, res/2, sourcename);
    } else if (source == TEST_SOURCE_TSF) {
	    if (farfield == TEST_FF_PFF)
	        fprintf(fwpar, "SOURCE_TSF\n20 20 20 %d %d %d 0 0 1.57079632679 2 20e-6 10 5\n\n", res-20, res-20, res-20);
	    else
	        fprintf(fwpar, "SOURCE_TSF\n20 20 20 %d %d %d 1.57079632679 1.57079632679 0 2 20e-6 10 5\n\n", res-20, res-20, res-20);
    } else if (source == TEST_SOURCE_SF) {
	    fprintf(fwpar, "SOURCE_SF\n1.57079632679 1.57079632679 0 2 20e-6 10 5\n\n");
    } else if (source == TEST_SOURCE_FOCUSED) {
	    fprintf(fwpar, "SOURCE_TSFF\n18 18 18 %d %d %d 0.4 1000 0 5 5 2 20e-6 100 5\n\n", res-18, res-18, res-18);
    } else if (source == TEST_SOURCE_LTSF) {
	    fprintf(fwpar, "SOURCE_LTSF\n18 18 18 %d %d %d 0.707 0.3 0 2 %d 4 1 0 0 %d 1 1 0 0 2 20e-6 100 5\n\n", res-18, res-18, res-18, res/2-15, res/2+15);
    } else if (source == TEST_SOURCE_LTSFF) {
	    fprintf(fwpar, "SOURCE_LTSFF\n18 18 18 %d %d %d 0.4 1000 0 5 5 2 %d 4 1 0 0 %d 1 1 0 0 2 20e-6 100 5\n\n", res-18, res-18, res-18, res/2-15, res/2+15);
    }


    /*boundary conditions*/
    if (bnd == SV_BOUNDARY_LIAO)
	    fprintf(fwpar, "BOUNDARY_ALL\nliao\n\n");
    else if (bnd == SV_BOUNDARY_CPML)
	    fprintf(fwpar, "BOUNDARY_ALL\ncpml 10 3 -1 0.03 4\n\n");

    fclose(fwpar);

    printf("Running test (%s %s)...\n", argv[0], argv[1]);

    g_get_current_time(&starttv);
    if (g_spawn_sync(NULL, argv, NULL, G_SPAWN_SEARCH_PATH, NULL, NULL, &outlog, &errlog, &returnstatus, &error)); // ??? (RS)
    else {
	    printf("Error spawning process. %s\n", error->message);
	    return -1;
    }
    g_get_current_time(&endtv);

    if (g_strstr_len(outlog, 2000, "failure") || g_strstr_len(outlog, 2000, "ERROR") || g_strstr_len(outlog, 2000, "failed")) {
	    printf("Error while running test (most typically CUDA launch failure).\n");
	    return -1;
    }

    if (g_strstr_len(errlog, 2000, "failure") || g_strstr_len(errlog, 2000, "ERROR") || g_strstr_len(errlog, 2000, "failed")) {
	    printf("Error while running test (most typically lack of memory)\n");
	    return -1;
    }

    if (diff != NULL) {
	    test_out = fopen(id, "r");
	    if (test_out == NULL) {
	        printf("Error: cannot open test result\n");
	        return -1;
	    }

	    std_out = fopen(testpath, "r");
	    if (std_out == NULL) {
	        printf("Error: canot open standard result at (%s)\n", testpath);
	        return -1;
	    }

	    *diff = 0;
	    sum = 0;
	    for (i = 0; i < (nsteps); i++) {
            if (!farfield) {
	    	    fscanf(test_out, "%d %f", &pos1, &val1);
	    	    fscanf(std_out, "%d %f", &pos2, &val2);
            } else {
                fscanf(test_out, "%d %f", &pos1, &val1a);
                fscanf(test_out, "%f", &val1b);
                fscanf(test_out, "%f", &val1c);
                fscanf(std_out, "%d %f", &pos2, &val2a);
                fscanf(std_out, "%f", &val2b);
                fscanf(std_out, "%f", &val2c);
                val1 = (gfloat)(sqrt(val1a*val1a + val1b*val1b + val1c*val1c));
                val2 = (gfloat)(sqrt(val2a*val2a + val2b*val2b + val2c*val2c));
            }

	        if (pos1 != pos2) {
	    	    printf("Error: Standard result at %s seems to belong to a different task.\n", testpath);
	    	    return -1;
	        }
	        
	        *diff += (val2-val1)*(val2-val1);
	        sum += val1*val1;
	        //printf("%d %g %g %g\n", pos1, val1, val2, val2 - val1);
	    } // i
	    
	    fclose(test_out);
	    fclose(std_out);
    } // if DIFF!= NULL

    if (diff != NULL) {
	    printf("Test time %g s, difference from expected result %g (%g %%)\n",
	        (gdouble)(endtv.tv_sec - starttv.tv_sec)+(gdouble)(endtv.tv_usec - starttv.tv_usec)/G_USEC_PER_SEC, *diff, (*diff)/sum*100.0);
    }

    for (i = 0; i < 3; i++)
	    g_free(argv[i]);
    g_free(argv); 

    return (gdouble)(endtv.tv_sec - starttv.tv_sec) + (gdouble)(endtv.tv_usec - starttv.tv_usec) / G_USEC_PER_SEC;
}

void
run_tests(gint level)
{
    gdouble cputime, gputime, mttime;
    gdouble result;
    gdouble cpu_time[11];
    gdouble gpu_time[11];
    gchar request[10], resrequest[3], *r;
    gint i, j, k, l, m, n, ncores;
    gint res, size;
    
#ifdef UCUDA
    cudaError_t err;
    struct cudaDeviceProp prop;
    gint ngpu = 0;
    gboolean rungpu = 1;
#endif

#ifndef UCUDA
    gboolean rungpu = 0;
    gint ngpu = 0;
#endif

    printf("Running GSvit tests at level %d\n", level);
    ncores = get_n_cores();
    printf("System has %d cores\n", ncores);
    
#ifdef G_OS_WIN32
    printf("Gsvit is installed at %s\\bin\n", get_self_dir());
#else
    printf("Gsvit is installed at %s\n", g_path_get_dirname(get_self_dir()));
#endif
  
#ifdef UCUDA
    printf("Program is compiled with GPU support\nSearching for available GPUs...\n");

    err = cudaGetDeviceCount(&ngpu);
    if (err){
	    printf("Error: %s\n", cudaGetErrorString(err));
	    ngpu = 0;
    }
    printf("Found %d GPUs\n", ngpu);

    for (i = 0; i < ngpu; i++) {
        err = cudaGetDeviceProperties(&prop, i);
        printf("The Properties of the Device with ID %d are\n", i);
        printf("Device Name             : %s\n", prop.name);
        printf("Device Memory Size      : %u\n", (guint)prop.totalGlobalMem);
        printf("Block Shared memory size: %d\n", (gint)prop.sharedMemPerBlock);
        printf("Max grid size           : %dx%dx%d\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
        printf("Max threads dim         : %dx%dx%d\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
    }
#endif

    if (!rungpu)
	    printf("Program is compiled without GPU support\n");

    if (level == 1) {
	    printf("%d simple tests scheduled in this test level\n", ngpu+1);
	    printf("Test 1:\n");
	    cputime = make_test(100, SV_BOUNDARY_NONE, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, 1);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
		        gputime = make_test(100, SV_BOUNDARY_NONE, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 201, 1);
		        printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
    } // if level == 1

    if (level == 2) {
        printf("%d tests scheduled in this test level\n", (rungpu+2)*16);


	    printf("___________________________________________________\nTest 1: free space, no BC, nothing special\n");
	    cputime = make_test(100, SV_BOUNDARY_NONE, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, 1);
	    mttime = make_test(100, SV_BOUNDARY_NONE, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_NONE, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 201, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }



	    printf("___________________________________________________\nTest 2: electric material, matmode checking on, Liao BC\n");
	    cputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, 1);
        mttime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, -1);
        printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 201, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }

	
	    printf("___________________________________________________\nTest 3: electric tabulated material, matmode checking on, Liao BC\n");
	    cputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, 1);
        mttime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, -1);
        printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 201, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }

	    printf("___________________________________________________\nTest 4: magnetic material, matmode checking on, Liao BC\n");
	    cputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_MAGNETIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, 1);
        mttime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_MAGNETIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, -1);
        printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_MAGNETIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 201, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }

	    
	    printf("___________________________________________________\nTest 5: magnetic tabulated material, matmode checking on, Liao BC\n");
	    cputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_MAGNETIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, 1);
        mttime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_MAGNETIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, -1);
        printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_MAGNETIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 201, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }

	
	    printf("___________________________________________________\nTest 6: CPML bc, far field calculation\n");
	    cputime = make_test(100, SV_BOUNDARY_CPML, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_YES, TEST_SOURCE_POINT, 0, 0, &result, 201, 1);
	    mttime = make_test(100, SV_BOUNDARY_CPML, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_YES, TEST_SOURCE_POINT, 0, 0, &result, 201, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_CPML, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_YES, TEST_SOURCE_POINT, 1, i, &result, 201, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	    
	    printf("___________________________________________________\nTest 7: TSF source, CPML bc, PEC, far field calculation\n");
	    cputime = make_test(100, SV_BOUNDARY_CPML, TEST_MATERIAL_LINTAB_PEC, TEST_MATCHECK_NO, TEST_FF_YES, TEST_SOURCE_TSF, 0, 0, &result, 201, 1);
	    mttime = make_test(100, SV_BOUNDARY_CPML, TEST_MATERIAL_LINTAB_PEC, TEST_MATCHECK_NO, TEST_FF_YES, TEST_SOURCE_TSF, 0, 0, &result, 201, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_CPML, TEST_MATERIAL_LINTAB_PEC, TEST_MATCHECK_NO, TEST_FF_YES, TEST_SOURCE_TSF, 1, i, &result, 201, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	    
	    printf("___________________________________________________\nTest 8: TSF source, CPML bc, far field calculation, electric material, matmode checking on\n");
	    cputime = make_test(100, SV_BOUNDARY_CPML, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_YES, TEST_SOURCE_TSF, 0, 0, &result, 201, 1);
	    mttime = make_test(100, SV_BOUNDARY_CPML, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_YES, TEST_SOURCE_TSF, 0, 0, &result, 201, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	    	for (i=0; i<ngpu; i++) {
	    		gputime = make_test(100, SV_BOUNDARY_CPML, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_YES, TEST_SOURCE_TSF, 1, i, &result, 201, 1);
	    		printf("GPU %d speedup: %g\n", i, cputime/gputime);
	    	}
	    }
	
	    printf("___________________________________________________\nTest 9: SF source, PEC, Liao BC\n");
	    cputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_PEC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_SF, 0, 0, &result, 201, 1);
	    mttime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_PEC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_SF, 0, 0, &result, 201, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);        
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_PEC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_SF, 1, i, &result, 201, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }

	    printf("___________________________________________________\nTest 10: Debye material\n");
	    cputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_DRUDE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, 1);
        mttime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_DRUDE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, -1);
        printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_DRUDE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 201, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	
	    printf("___________________________________________________\nTest 11: CP material\n");
	    cputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_CP, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, 1);
        mttime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_CP, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, -1);
        printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_CP, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 201, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	    
	    printf("___________________________________________________\nTest 12: ADE material\n");
	    cputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_ADE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, 1);
            mttime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_ADE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, -1);
            printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_ADE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 201, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	    
	    printf("___________________________________________________\nTest 13: PLRC material\n");
	    cputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_PLRC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, 1);
        mttime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_PLRC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 201, -1);
        printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_PLRC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 201, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }

	    printf("___________________________________________________\nTest 14: TSFF source, Liao bc, PEC\n");
	    cputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_FOCUSED, 0, 0, &result, 801, 1);
        mttime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_FOCUSED, 0, 0, &result, 801, -1);
        printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i= 0 ; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_FOCUSED, 1, i, &result, 801, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }


	    printf("___________________________________________________\nTest 15: TSF source, periodic BC, periodic far field calculation, PEC\n");
	    cputime = make_test(100, SV_BOUNDARY_NONE, TEST_MATERIAL_LINTAB_PEC, TEST_MATCHECK_NO, TEST_FF_PFF, TEST_SOURCE_TSF, 0, 0, &result, 201, 1);
	    mttime = make_test(100, SV_BOUNDARY_NONE, TEST_MATERIAL_LINTAB_PEC, TEST_MATCHECK_NO, TEST_FF_PFF, TEST_SOURCE_TSF, 0, 0, &result, 201, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_NONE, TEST_MATERIAL_LINTAB_PEC, TEST_MATCHECK_NO, TEST_FF_PFF, TEST_SOURCE_TSF, 1, i, &result, 201, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }


	    printf("___________________________________________________\nTest 16: LTSF source, Liao bc, PEC\n");
	    cputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_LTSF, 0, 0, &result, 801, 1);
            mttime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_LTSF, 0, 0, &result, 801, -1);
            printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i= 0 ; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_LTSF, 1, i, &result, 801, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
/*
	    printf("___________________________________________________\nTest 17: LTSFFF source, Liao bc, PEC\n");
	    cputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_LTSFF, 0, 0, &result, 801, 1);
            mttime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_LTSFF, 0, 0, &result, 801, -1);
            printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i= 0 ; i < ngpu; i++) {
	    	    gputime = make_test(100, SV_BOUNDARY_LIAO, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_LTSFF, 1, i, &result, 801, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
*/


    } // if level == 2

    if (level == 3) {
        printf("%d tests scheduled in this test level\n", (rungpu+2)*15);

        
	    printf("___________________________________________________\nTest 1: free space, no BC, nothing special\n");
	    cputime = make_test(200, SV_BOUNDARY_NONE, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, 1);
	    mttime = make_test(200, SV_BOUNDARY_NONE, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(200, SV_BOUNDARY_NONE, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 401, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	    
	    printf("___________________________________________________\nTest 2: electric material, matmode checking on, Liao BC\n");
	    cputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, 1);
	    mttime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 401, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	    
	    printf("___________________________________________________\nTest 3: electric tabulated material, matmode checking on, Liao BC\n");
	    cputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, 1);
	    mttime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 401, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }

	    printf("___________________________________________________\nTest 4: magnetic material, matmode checking on, Liao BC\n");
	    cputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_MAGNETIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, 1);
	    mttime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_MAGNETIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_MAGNETIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 401, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	    
	    printf("___________________________________________________\nTest 5: magnetic tabulated material, matmode checking on, Liao BC\n");
	    cputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_MAGNETIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, 1);
	    mttime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_MAGNETIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_MAGNETIC, TEST_MATCHECK_YES, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 401, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	    

	    printf("___________________________________________________\nTest 6: CPML bc, far field calculation\n");
	    cputime = make_test(200, SV_BOUNDARY_CPML, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_YES, TEST_SOURCE_POINT, 0, 0, &result, 401, 1);
	    mttime = make_test(200, SV_BOUNDARY_CPML, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_YES, TEST_SOURCE_POINT, 0, 0, &result, 401, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(200, SV_BOUNDARY_CPML, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_YES, TEST_SOURCE_POINT, 1, i, &result, 401, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	    
	    printf("___________________________________________________\nTest 7: TSF source, CPML bc, PEC, far field calculation\n");
	    cputime = make_test(200, SV_BOUNDARY_CPML, TEST_MATERIAL_LINTAB_PEC, TEST_MATCHECK_NO, TEST_FF_YES, TEST_SOURCE_TSF, 0, 0, &result, 401, 1);
	    mttime = make_test(200, SV_BOUNDARY_CPML, TEST_MATERIAL_LINTAB_PEC, TEST_MATCHECK_NO, TEST_FF_YES, TEST_SOURCE_TSF, 0, 0, &result, 401, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(200, SV_BOUNDARY_CPML, TEST_MATERIAL_LINTAB_PEC, TEST_MATCHECK_NO, TEST_FF_YES, TEST_SOURCE_TSF, 1, i, &result, 401, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	    
	    printf("___________________________________________________\nTest 8: TSF source, CPML bc, far field calculation, electric material, matmode checking on\n");
	    cputime = make_test(200, SV_BOUNDARY_CPML, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_YES, TEST_SOURCE_TSF, 0, 0, &result, 401, 1);
	    mttime = make_test(200, SV_BOUNDARY_CPML, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_YES, TEST_SOURCE_TSF, 0, 0, &result, 401, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(200, SV_BOUNDARY_CPML, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_YES, TEST_FF_YES, TEST_SOURCE_TSF, 1, i, &result, 401, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	    
	    printf("___________________________________________________\nTest 9: SF source, PEC, Liao BC\n");
	    cputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_PEC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_SF, 0, 0, &result, 401, 1);
	    mttime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_PEC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_SF, 0, 0, &result, 401, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_PEC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_SF, 1, i, &result, 401, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	    
	    printf("___________________________________________________\nTest 10: Debye material\n");
	    cputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_DRUDE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, 1);
	    mttime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_DRUDE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_DRUDE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 401, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	   
 
	    printf("___________________________________________________\nTest 11: CP material\n");
	    cputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_CP, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, 1);
	    mttime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_CP, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_CP, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 401, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	    
	    printf("___________________________________________________\nTest 12: ADE material\n");
	    cputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_ADE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, 1);
	    mttime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_ADE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_ADE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 401, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	    
	    printf("___________________________________________________\nTest 13: PLRC material\n");
	    cputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_PLRC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, 1);
	    mttime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_PLRC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, &result, 401, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_LINTAB_PLRC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, i, &result, 401, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }

	    printf("___________________________________________________\nTest 14: TSFF source, Liao bc, PEC\n");
	    cputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_FOCUSED, 0, 0, &result, 801, 1);
	    mttime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_FOCUSED, 0, 0, &result, 801, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(200, SV_BOUNDARY_LIAO, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_FOCUSED, 1, i, &result, 801, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
	    
	    printf("___________________________________________________\nTest 15: TSF source, periodic BC, periodic far field calculation, PEC\n");
	    cputime = make_test(200, SV_BOUNDARY_NONE, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_PFF, TEST_SOURCE_TSF, 0, 0, &result, 401, 1);
	    mttime = make_test(200, SV_BOUNDARY_NONE, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_PFF, TEST_SOURCE_TSF, 0, 0, &result, 401, -1);
	    printf("Multithreading speedup: %g\n", cputime/mttime);
	    if (rungpu) {
	        for (i = 0; i < ngpu; i++) {
	    	    gputime = make_test(200, SV_BOUNDARY_NONE, TEST_MATERIAL_NONE, TEST_MATCHECK_NO, TEST_FF_PFF, TEST_SOURCE_TSF, 1, i, &result, 401, 1);
	    	    printf("GPU %d speedup: %g\n", i, cputime/gputime);
	        }
	    }
    } // if level == 3

    /*
    if (level==3) {
        printf("%d tests scheduled in this test level\n", (rungpu+1)*(TEST_MATCHECK_YES+1)*(SV_BOUNDARY_CPML+1)*(SV_MAT_CP+1)*(TEST_FF_YES+1));
	
	for (i=0; i<=rungpu; i++) {   //iterate over gpu/nongpu
	    for (j=0; j<=TEST_MATCHECK_YES; j++) {  //iterate over checking material option
	        for (k=0; k<=SV_BOUNDARY_CPML; k++) { //iterate over different boundary conditions
		    for (l=0; l<=SV_MAT_CP; l++) { //iterate over different materials in the pool
		        for (m=0; m<=TEST_FF_YES; m++) { //iterate over farfield point instead of local one
			    for (n=0; n<=TEST_SOURCE_TSF; n++) { //iterate over different sources

			        cputime = make_test(100, k, l, j, m, n, 0, 0, &result);
				if (rungpu) {
				    gputime = make_test(100, k, l, j, m, n, 1, 0, &result);
				    printf("GPU %d speedup: %g\n", i, cputime/gputime);
				}
			    } // n
			} // m
		    } // l
		} // k               
	    } // j           
	} // i

    } // if 
    */
    
    if (level == 12) {
        printf("11 tests scheduled in this test level\n");
        for (i = 0; i < 11; i++) {
	        printf("Test %d:\n", i);
	        cpu_time[i] = make_test(100 + 10*i, SV_BOUNDARY_NONE, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, NULL, 201, 1);
	        if (rungpu) {
		        gpu_time[i] = make_test(100 + 10*i, SV_BOUNDARY_NONE, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, 0, NULL, 201, 1);
		        printf("GPU speedup: %g\n", cpu_time[i]/gpu_time[i]);
	        }
        }
        printf("---------------------------------------------------------\n");
        printf("Test summary:");
        if (rungpu)
	        printf("#cube size     cpu time     gpu time    speedup\n");
        else
	        printf("#cube size     cpu time\n");
        for (i = 0; i < 11; i++) {
	        if (rungpu)
		        printf("%d %g %g %g\n", 100 + 10*i, cpu_time[i], gpu_time[i], cpu_time[i]/gpu_time[i]);
	        else
		        printf("%d %g\n", 100 + 10*i, cpu_time[i]);
	    }
    } // if level == 12

    if (level == 13) {
        printf("11 tests scheduled in this test level\n");
        for (i = 0; i < 11; i++) {
	        printf("Test %d:\n", i);
	        cpu_time[i] = make_test(100 + 20*i, SV_BOUNDARY_NONE, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, NULL, 201, 1);
	        if (rungpu) {
		        gpu_time[i] = make_test(100 + 20*i, SV_BOUNDARY_NONE, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, 0, NULL, 201, 1);
		        printf("GPU speedup: %g\n", cpu_time[i]/gpu_time[i]);
	        }
        }
        printf("---------------------------------------------------------\n");
        printf("Test summary:");
        if (rungpu)
	        printf("#cube size     cpu time     gpu time    speedup\n");
        else
	        printf("#cube size     cpu time\n");
        for (i = 0; i < 11; i++) {
	        if (rungpu)
		        printf("%d %g %g %g\n", 100 + 20*i, cpu_time[i], gpu_time[i], cpu_time[i]/gpu_time[i]);
	        else
		        printf("%d %g\n", 100 + 20*i, cpu_time[i]);
	    }
    } // if level == 13

    if (level == 14) {
        printf("11 tests scheduled in this test level\n");
        for (i = 0; i < 11; i++) {
	        printf("Test %d:\n", i);
	        cpu_time[i] = make_test(100 + 30*i, SV_BOUNDARY_NONE, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, NULL, 201, 1);
	        if (rungpu) {
		        gpu_time[i] = make_test(100 + 30*i, SV_BOUNDARY_NONE, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, 0, NULL, 201, 1);
		        printf("GPU speedup: %g\n", cpu_time[i]/gpu_time[i]);
	        }
        }
        printf("---------------------------------------------------------\n");
        printf("Test summary:");
        if (rungpu)
	        printf("#cube size     cpu time     gpu time    speedup\n");
        else
	        printf("#cube size     cpu time\n");
        for (i = 0; i < 11; i++) {
	        if (rungpu)
		        printf("%d %g %g %g\n", 100 + 30*i, cpu_time[i], gpu_time[i], cpu_time[i]/gpu_time[i]);
	        else
		        printf("%d %g\n", 100 + 30*i, cpu_time[i]);
	    }
    } // if level == 14

    if (level == 15) {
        printf("11 tests scheduled in this test level\n");
        for (i = 0; i < 11; i++) {
	        printf("Test %d:\n", i);
	        cpu_time[i] = make_test(100 + 40*i, SV_BOUNDARY_NONE, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, NULL, 201, 1);
	        if (rungpu) {
		        gpu_time[i] = make_test(100 + 40*i, SV_BOUNDARY_NONE, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, 0, NULL, 201, 1);
		        printf("GPU speedup: %g\n", cpu_time[i]/gpu_time[i]);
	        }
        }
        printf("---------------------------------------------------------\n");
        printf("Test summary:");
        if (rungpu)
	        printf("#cube size     cpu time     gpu time    speedup\n");
        else
	        printf("#cube size     cpu time\n");
        for (i = 0; i < 11; i++) {
	        if (rungpu)
		        printf("%d %g %g %g\n", 100 + 40*i, cpu_time[i], gpu_time[i], cpu_time[i]/gpu_time[i]);
	        else
		        printf("%d %g\n", 100 + 40*i, cpu_time[i]);
	    }
    } // if level == 15

    if (level == 20) {
        printf("11 tests scheduled in this test level\n");
        for (i = 0; i < 11; i++) {
	        printf("Test %d:\n", i);
	        if (rungpu)
		        gpu_time[i] = make_test(100 + 30*i, SV_BOUNDARY_NONE, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 1, 0, NULL, 201, 1);
        }
        printf("---------------------------------------------------------\n");
        printf("Test summary:");
        if (rungpu)
	        printf("#cube size     gpu time\n");
        for (i = 0; i < 11; i++) {
	        if (rungpu)
		        printf("%d %g\n", 100 + 30*i, gpu_time[i]);
	    }
    } // if level == 20

    if (level > 30 && level < 40) { 
        printf("%d tests scheduled in this test level\n", ncores);
        size = (level-30)*100;
        for (i = 1; i <= ncores; i++) {
	        printf("Test %d:\n", i);
	        cpu_time[i] = make_test(size, SV_BOUNDARY_NONE, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, NULL, 201, i);
	        printf("Result for %d threads: time %g,  speedup %g\n", i, cpu_time[i], cpu_time[1]/cpu_time[i]);
        }
        printf("---------------------------------------------------------\n");
        printf("Test summary for cube size %dx%dx%d\n", size, size, size);
        printf("#cores     cpu time   speedup\n");
        for (i = 1; i <= ncores; i++) 
	        printf("%d %g %g\n", i, cpu_time[i], cpu_time[1]/cpu_time[i]);
    } // if 30 < level < 40

    if (level > 10000) {
	    g_snprintf(request, sizeof(request), "%d", level);
	    resrequest[0] = request[0];
	    resrequest[1] = request[1];
	    resrequest[2] = request[2];
	    res = atoi(resrequest);            
	    r = g_strdup_printf("%c", request[3]); j = atoi(r);
	    r = g_strdup_printf("%c", request[4]); k = atoi(r);
	    r = g_strdup_printf("%c", request[5]); l = atoi(r);
	    r = g_strdup_printf("%c", request[6]); m = atoi(r);
	    r = g_strdup_printf("%c", request[7]); n = atoi(r);
	    printf("1 test scheduled in this test level: resolution %d boundary %d, material %d, matcheck %d, farfield %d, source %d, gpu %d, ugpu %d\n", res, j, k, l, m, n, 0, 0);
	    if (res > 101)
	        cputime = make_test(res, (SvBoundaryType)j, (TestMaterialType)k, (TestMatcheckType)l, (TestFarfieldType)m, (TestSourceType)n, 0, 0, &result, 301, 1);
	    else
	        cputime = make_test(res, (SvBoundaryType)j, (TestMaterialType)k, (TestMatcheckType)l, (TestFarfieldType)m, (TestSourceType)n, 0, 0, &result, 201, 1);
    } // if level > 10000
}

void
run_benchmark(gint level)
{
    gdouble cputime;
    gint ncores;
    gint size = 900;

    printf("Running GSvit benchmark for %d cores\n", level);
    ncores = get_n_cores();
    printf("System has %d cores\n", ncores);
    printf("Gsvit is installed at %s\n", g_path_get_dirname(get_self_dir()));
    printf("Running test...\n");
  
    cputime = make_test(size, SV_BOUNDARY_NONE, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, NULL, 201, level);
    printf("Test result for cube size %dx%dx%d\n", size, size, size);
    printf("#cores   cpu time\n");
    printf("%d %g\n", level, cputime);
}

void
run_bigbenchmark(gint level)
{
    gdouble cputime;
    gint ncores;
    gint size = 4000;

    printf("Running GSvit very large scale benchmark for %d cores\n", level);
    ncores = get_n_cores();
    printf("System has %d cores\n", ncores);
    printf("Gsvit is installed at %s\n", g_path_get_dirname(get_self_dir()));
    printf("Running test...\n");
  
    cputime = make_test(size, SV_BOUNDARY_NONE, TEST_MATERIAL_ELECTRIC, TEST_MATCHECK_NO, TEST_FF_NO, TEST_SOURCE_POINT, 0, 0, NULL, 201, level);
    printf("Test result for cube size %dx%dx%d\n", size, size, size);
    printf("#cores   cpu time\n");
    printf("%d %g\n", level, cputime);
}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
