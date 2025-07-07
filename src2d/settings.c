
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


/*  settings.c : 
 *  all the main structures that are used for driving computation, 
 *  generally data loaded from different input files
 */


#include <stdlib.h>
#include <stdio.h>
#include <glib.h>
#include <math.h>
#include "settings.h"
#include "constants.h"
#include <string.h>
#ifndef G_OS_WIN32
#include <unistd.h>
#endif


typedef enum {
   SRC_DIRECT = 0,
   SRC_SINE = 1,
   SRC_PULSE = 2
} SrcType;

typedef enum {
   SRCF_FULL = 0,
   SRCF_TSF = 1,
   SRCF_LTSF = 2,
   SRCF_ARB = 3
} SrcFormat;

gchar*
get_self_dir()
{
#ifndef G_OS_WIN32
    gchar buf[1024];
    memset(buf, 0, sizeof(buf));
    readlink("/proc/self/exe", buf, sizeof(buf)-1);
    return g_strdup(buf);
#endif
    
#ifdef G_OS_WIN32 
#if (GLIB_CHECK_VERSION(2, 16, 0))
    return g_win32_get_package_installation_directory_of_module(NULL);
#else
    return g_win32_get_package_installation_directory(NULL, NULL);
#endif
#endif

    return NULL;
}

gint
get_int(FILE *F, gint *val, gchar *key)
{
    gchar value[100], *err;

    fscanf(F, "%s", value);;
    *val = strtol(value, &err, 10);

    if (*err != '\0')
    {
        fprintf(stderr, "Error parsing integer value (key %s)\n", key);
        return 1;
    }
    return 0;
}

gint
get_double(FILE *F, gdouble *val, gchar *key)
{
    gchar value[100], *err;

    fscanf(F, "%s", value);;
    *val = strtod(value, &err);

    if (*err != '\0')
    {
        fprintf(stderr, "Error parsing double value (key %s)\n", key);
        return 1;
    }
    return 0;
}

gint
get_component(FILE *fr, gint *val, gchar *key)
{
    gchar component[20];
    fscanf(fr, "%s", component);

    if (strcmp(component, "Ex") == 0 || strcmp(component, "ex") == 0)
        *val = SV_COMP_EX;
    else if (strcmp(component, "Ey") == 0 || strcmp(component, "ey") == 0)
        *val = SV_COMP_EY;
    else if (strcmp(component, "Ez") == 0 || strcmp(component, "ez") == 0)
        *val = SV_COMP_EZ;
    else if (strcmp(component, "Hx") == 0 || strcmp(component, "hx") == 0)
        *val = SV_COMP_HX;
    else if (strcmp(component, "Hy") == 0 || strcmp(component, "hy") == 0)
        *val = SV_COMP_HY;
    else if (strcmp(component, "Hz") == 0 || strcmp(component, "hz") == 0)
        *val = SV_COMP_HZ;
    else if (strcmp(component, "All") == 0 || strcmp(component, "all") == 0)
        *val = SV_COMP_ALL;
    else if (strcmp(component, "Cur") == 0 || strcmp(component, "cur") == 0)
        *val = SV_COMP_CUR;
    else if (strcmp(component, "Epsilon") == 0 || strcmp(component, "epsilon") == 0)
        *val = SV_COMP_EPSILON;
    else if (strcmp(component, "Sigma") == 0 || strcmp(component, "sigma") == 0)
        *val = SV_COMP_SIGMA;
    else if (strcmp(component, "Mu") == 0 || strcmp(component, "mu") == 0)
        *val = SV_COMP_MU;
    else if (strcmp(component, "Sigast") == 0 || strcmp(component, "sigast") == 0)
        *val = SV_COMP_SIGAST;
    else
    {
        fprintf(stderr, "Error parsing field component (key %s)\n", key);
        return 1;
    }
    return 0;
}

gint
get_scomponent(FILE *fr, gint *val, gchar *key)
{
    gchar component[5];
    fscanf(fr, "%s", component);

    if (strcmp(component, "Ex") == 0 || strcmp(component, "ex") == 0)
        *val = SV_SUM_EX;
    else if (strcmp(component, "Ey") == 0 || strcmp(component, "ey") == 0)
        *val = SV_SUM_EY;
    else if (strcmp(component, "Ez") == 0 || strcmp(component, "ez") == 0)
        *val = SV_SUM_EZ;
    else if (strcmp(component, "All") == 0 || strcmp(component, "all") == 0)
        *val = SV_SUM_ALL;
    else if (strcmp(component, "Abs") == 0 || strcmp(component, "abs") == 0)
        *val = SV_SUM_ABS;
    else
    {
        fprintf(stderr, "Error parsing field component (key %s)\n", key);
        return 1;
    }
    return 0;
}

gint
write_source(SrcFormat format, SrcType type, gint nsteps, gdouble wavelength, gdouble size, gdouble dt, gint component, gchar *filename, gdouble amplitude, gdouble theta, gdouble phi)
{
    gint i, dampval = 100;
    gdouble val, wavefreq, damp, ax, ay, az;
    FILE *fw = fopen(filename, "w");

    if (!fw) {
        fprintf(stderr, "Error: cannot write source file: %s\n", filename);
        return 1;
    }

    if (format==SRCF_LTSF) {
        dt/=3;
        nsteps *= 3;
        size *= 3;
        dampval *= 3;
    }

    wavefreq = LIGHT_SPEED*2*G_PI/wavelength;
    if (type==SRC_PULSE) size *= dt;

    fprintf(fw, "%d\n", nsteps);
    for (i=0; i<nsteps; i++)
    {
        val = amplitude*sin(wavefreq*dt*(gdouble)i);
        if (i<dampval) val *= (1.0/(gdouble)dampval)*((gdouble)(i));

        if (type==SRC_PULSE) {
           damp = ((gdouble)i*dt - 4.0*size)/sqrt((gdouble)2.0)/size;
           val *= exp(-damp*damp);
        }
        if (format==SRCF_FULL) {
            if (component==SV_COMP_EX) fprintf(fw, "%d %g %g %g %g %g %g\n",
                     i, val, 0.0, 0.0, 0.0, 0.0, 0.0);
            else if (component==SV_COMP_EY) fprintf(fw, "%d %g %g %g %g %g %g\n",
                     i, 0.0, val, 0.0, 0.0, 0.0, 0.0);
            else if (component==SV_COMP_EZ) fprintf(fw, "%d %g %g %g %g %g %g\n",
                     i, 0.0, 0.0, val, 0.0, 0.0, 0.0);
            else if (component==SV_COMP_HX) fprintf(fw, "%d %g %g %g %g %g %g\n",
                     i, 0.0, 0.0, 0.0, val, 0.0, 0.0);
            else if (component==SV_COMP_HY) fprintf(fw, "%d %g %g %g %g %g %g\n",
                     i, 0.0, 0.0, 0.0, 0.0, val, 0.0);
            else if (component==SV_COMP_HZ) fprintf(fw, "%d %g %g %g %g %g %g\n",
                     i, 0.0, 0.0, 0.0, 0.0, 0.0, val);

        }
        if (format==SRCF_ARB) {
            ax = sin(theta)*cos(phi);
            ay = sin(theta)*sin(phi);
            az = cos(theta);

            fprintf(fw, "%d %g %g %g %g %g %g\n",
                     i, ax*val, ay*val, az*val, 0.0, 0.0, 0.0);
        } else {
            fprintf(fw, "%d %g\n", i, val);

        }

    }

    fclose(fw);

   // if (format==SRCF_TSF) fprintf(fw, "%d\n");

    return 0;
}



int clear_settings(SvSet *set)
{
    gint i;

    /*no sources by default*/
    set->ss.npnts = 0;
    set->ss.tsf.i0 = set->ss.tsf.j0 = 0;
    set->ss.tsf.i1 = set->ss.tsf.j1 = 0;
    set->ss.tsf.epsilon = 1;
    set->ss.tsf.mu = 1;

    /*no output by default*/
    set->so.npnts = 0;
    set->so.nsums = 0;
    set->so.nplns = 0;
    set->so.nimgs = 0;
    set->sc.verbose = 0;
    set->sc.step_act = 0;
    /*no material data by default*/
    set->sm.in_material = NULL;
    set->sm.in_vector = NULL;
    /*none boundary condition by default (including material boundary)*/
    set->sb.bx0 = set->sb.bxn = SV_BOUNDARY_NONE;
    set->sb.by0 = set->sb.byn = SV_BOUNDARY_NONE;
    set->smb.bx0 = set->smb.bxn = SV_BOUNDARY_NONE;
    set->smb.by0 = set->smb.byn = SV_BOUNDARY_NONE;
    set->smb.bx0pos = set->smb.by0pos = 0;
    set->smb.bxnpos = set->smb.bynpos = -1;
    /*te mode by default*/
    set->tmmode = 1;
    /*none nfff data by default*/
    set->sf.savefile = NULL;
    set->sf.i0 = set->sf.i1 = set->sf.j0 = set->sf.j1 = 0;
    set->sf.nrs = 0;
    set->sc.nthreads = -1;


    /*use full material properties algorithm by default*/
    set->sm.matmode_check = 0;
    /*do not use GPU by default*/
    set->sc.ngpu = 0;
    set->so.outfile = g_strdup_printf("out.gwy");
    for (i=0; i<MAX_GPUS; i++) set->sc.ugpu[i] = 0;

    return 0;
}

void remove_temporary_files(SvSet *set)
{
    gint i;

    for (i=0; i<set->ss.npnts; i++)
    {
	    if (set->ss.pnts[i].mode != 0) {
		    printf("Removing file %s\n", set->ss.pnts[i].filename);

		    //preliminary hack
		    remove(set->ss.pnts[i].filename);
	    }
    }
    if ((set->ss.tsf.i0 + set->ss.tsf.j0
         + set->ss.tsf.i1 + set->ss.tsf.j1)!=0)
    {
          if (set->ss.tsf.mode != 0) {
                    printf("Removing file %s\n", set->ss.tsf.filename);

                    //preliminary hack
                    remove(set->ss.tsf.filename);
            }

    }


}

int parse_settings(gchar *filename, SvSet *set)
{
    gchar key[100];
    gchar value[20];
    gchar buffer[100];
    gchar tmpsourcefile[256];
    gint val;
    gint dimensionality = 3; //assuming 3d version by default
    gdouble amplitude, wavelength, theta, phi, size, dt;
    gint srctype;

    FILE *fr = fopen(filename, "r");
    int n_source_pnts_size = 0;
    int n_output_pnts_size = 0;
    int n_output_imgs_size = 0;
    int n_output_plns_size = 0;
    int n_output_sums_size = 0;

    if (!fr) 
        return 1;
 
    /*inits*/
    clear_settings(set);

    /*parse file*/
    fscanf(fr, "%s", key);

    while ((strstr(key, "#") != NULL)) {
	    fgets(buffer, 100, fr);
        fscanf(fr, "%s", key);
    }

    if (strcmp(key, "DIMENSIONALITY") == 0) {
	    if (get_int(fr, &dimensionality, key)) {
            goto exit;
        }
	    if (dimensionality != 2) {
            fprintf(stderr, "Error: Wrong parameter file dimensionality, it should be 2 for 2D version of Gvit"); 
            goto exit;
        }
    } else {
        fprintf(stderr, "Error: Missing file dimensionality, it should be 2 for 2D version of Gvit"); 
        goto exit;
    }

    while (fscanf(fr, "%s", key) != EOF)
    {
//        printf("key: %s\n", key);

        if (strcmp(key, "VERBOSE") == 0) {
            if (get_int(fr, &(set->sc.verbose), key)) 
                goto exit;
        } else if (strcmp(key, "POOL") == 0) {
            if (get_int(fr, &(set->sp.xres), key)) 
                goto exit;
            if (get_int(fr, &(set->sp.yres), key)) 
                goto exit;
            if (get_double(fr, &(set->sp.dx), key)) 
                goto exit;
            if (get_double(fr, &(set->sp.dy), key)) 
                goto exit;
        } else if (strcmp(key, "COMP") == 0) {
            if (get_int(fr, &(set->sc.nsteps), key)) 
                goto exit;
        } else if (strcmp(key, "GPU") == 0) {
            if (get_int(fr, &(set->sc.ngpu), key)) 
                goto exit;
        } else if (strcmp(key, "TMMODE") == 0) {
            if (get_int(fr, &(set->tmmode), key)) 
                goto exit;
            if (set->tmmode) 
                printf("TM mode simulation.\n");
            else 
                printf("TE mode simulation.\n");
        } else if (strcmp(key, "UGPU") == 0) {
            if (get_int(fr, &val, key)) 
                goto exit;
            set->sc.ugpu[val] = 1;
        } else if (strcmp(key, "THREADS") == 0) {
            if (get_int(fr, &(set->sc.nthreads), key)) {
                goto exit;
            }
        } else if (strcmp(key, "MATMODE_CHECK") == 0) {
            if (get_int(fr, &(set->sm.matmode_check), key)) 
                goto exit;
        } else if (strcmp(key, "OUT_FILE") == 0) {
            fscanf(fr, "%s", buffer);
            set->so.outfile = g_strdup(buffer);
        } else if (strcmp(key, "MEDIUM_LINEAR") == 0) { 
            fscanf(fr, "%s", buffer);
            set->sm.in_material = g_strdup(buffer);
            if (!g_file_test(set->sm.in_material, G_FILE_TEST_EXISTS)) {
                fprintf(stderr, "Error: File %s does not exist.\n", set->sm.in_material);
                goto exit;
            }
        } else if (strcmp(key, "MEDIUM_VECTOR") == 0) {
            fscanf(fr, "%s", buffer);
            set->sm.in_vector = g_strdup(buffer);
            if (!g_file_test(set->sm.in_vector, G_FILE_TEST_EXISTS)) {
                fprintf(stderr, "Error: File %s does not exist.\n", set->sm.in_vector);
                goto exit;
            }
        } else if (strcmp(key, "SOURCE_POINT") == 0) {
            /*alloc or realloc the necessary structure*/
            if (set->ss.npnts==0) {
                n_source_pnts_size = 10;
                set->ss.pnts = (SvSrcPoint *)g_malloc(n_source_pnts_size*sizeof(SvSrcPoint));
            } else if (set->ss.npnts >= n_source_pnts_size) { /*realloc*/
                n_source_pnts_size += 10;
                set->ss.pnts = (SvSrcPoint *)g_realloc(set->ss.pnts, n_source_pnts_size*sizeof(SvSrcPoint));
            }

            /*load parameters*/
            if (get_int(fr, &(set->ss.pnts[set->ss.npnts].i), key)) {
                goto exit;
            }
            if (get_int(fr, &(set->ss.pnts[set->ss.npnts].j), key)) {
                goto exit;
            }
            if (get_int(fr, &(srctype), key)) {
                goto exit;
            }
            set->ss.pnts[set->ss.npnts].mode = srctype;
            if (srctype==0) {  //load from file
                fscanf(fr, "%100[^\n\r]", buffer);
                set->ss.pnts[set->ss.npnts].filename=g_strstrip(g_strdup(buffer));
                if (set->sc.verbose>1) 
                    printf("Source will be loaded from file %s\n", set->ss.pnts[set->ss.npnts].filename);
                if (!g_file_test(set->ss.pnts[set->ss.npnts].filename, G_FILE_TEST_EXISTS)) {
                    fprintf(stderr, "Error: File %s does not exist.\n", set->ss.pnts[set->ss.npnts].filename);
                    goto exit;
               }
           } else if (srctype==1 || srctype==2) { //generate simple sources
                if (set->sc.nsteps == 0) {
                    fprintf(stderr, "Error: COMP parameter (number of steps) should precede generated source command.\n");
                    goto exit;
                }
                if (get_double(fr, &wavelength, key)) 
                    goto exit;
                if (srctype==2 && get_double(fr, &size, key)) 
                    goto exit;
                if (get_double(fr, &amplitude, key)) 
                    goto exit;

                if (get_double(fr, &theta, key)) 
                    goto exit;
                if (get_double(fr, &phi, key)) 
                    goto exit;

                set->ss.pnts[set->ss.npnts].source_amplitude = amplitude;
                set->ss.pnts[set->ss.npnts].theta = theta;
                set->ss.pnts[set->ss.npnts].phi = phi;

                set->ss.pnts[set->ss.npnts].source_wl = wavelength;
                set->ss.pnts[set->ss.npnts].length = size;

                printf("Generating %g %g\n", wavelength, amplitude);

                dt = 1.0/LIGHT_SPEED/sqrt(1.0/set->sp.dx/set->sp.dx+1.0/set->sp.dy/set->sp.dy);
                g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "tmp_source_%06d", set->sc.suffix);
                if (write_source(SRCF_ARB, (SrcType)srctype, set->sc.nsteps, wavelength, size, dt, 0, tmpsourcefile, amplitude, theta, phi)) 
                    goto exit;

                set->ss.pnts[set->ss.npnts].filename = g_strdup(tmpsourcefile);
            } else {
                fprintf(stderr, "Error: unknown source type\n");
                goto exit;
            }

            set->ss.npnts++;
        } else if (strcmp(key, "SOURCE_TSF") == 0) {
            if (get_int(fr, &(set->ss.tsf.i0), key)) 
                goto exit;
            if (get_int(fr, &(set->ss.tsf.j0), key)) 
                goto exit;
            if (get_int(fr, &(set->ss.tsf.i1), key)) 
                goto exit;
            if (get_int(fr, &(set->ss.tsf.j1), key)) 
                goto exit;
            if (get_double(fr, &(set->ss.tsf.phi), key)) 
                goto exit;
            if (get_double(fr, &(set->ss.tsf.psi), key)) 
                goto exit;
            set->ss.tsf.theta = G_PI/2.0;
            if (get_int(fr, &(srctype), key)) 
                goto exit;

            set->ss.tsf.mode = srctype;
            if (srctype==0) {  //load from file
                fscanf(fr, "%100[^\n\r]", buffer);
                set->ss.tsf.filename=g_strstrip(g_strdup(buffer));
                if (set->sc.verbose>1) 
                    printf("Source will be loaded from file %s\n", set->ss.tsf.filename);
            } else if (srctype==1 || srctype==2) { //generate simple sources
                if (set->sc.nsteps == 0) {
                    fprintf(stderr, "Error: COMP parameter (number of steps) should precede generated source command.\n");
                    goto exit;
                }
                if (get_double(fr, &wavelength, key))
                    goto exit;
                if (srctype==2 && get_double(fr, &size, key)) 
                    goto exit;
                if (get_double(fr, &amplitude, key)) 
                    goto exit;

                set->ss.tsf.source_wl = wavelength;
                set->ss.tsf.length = size;
                set->ss.tsf.source_amplitude = amplitude;

                dt = 1.0/LIGHT_SPEED/sqrt(1.0/set->sp.dx/set->sp.dx+1.0/set->sp.dy/set->sp.dy);
                g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "tmp_source_%06d", set->sc.suffix);
                if (write_source(SRCF_TSF, (SrcType)srctype, set->sc.nsteps, wavelength, size, dt, 0, tmpsourcefile, amplitude, 0, 0)) 
                    goto exit;

                set->ss.tsf.filename=g_strdup(tmpsourcefile);
            } else {
                fprintf(stderr, "Error: unknown source type\n");
                goto exit;
            }
        } else if (strcmp(key, "TSFSOURCE_MATERIAL") == 0) {
            if (get_double(fr, &(set->ss.tsf.epsilon), key)) 
                goto exit;
            if (get_double(fr, &(set->ss.tsf.mu), key)) 
                goto exit;
            if (get_double(fr, &(set->ss.tsf.sigma), key)) 
                goto exit;
            if (get_double(fr, &(set->ss.tsf.sigast), key)) 
                goto exit;
        } else if (strcmp(key, "OUT_POINT") == 0) {
            /*alloc or realloc the necessary structure*/
            if (set->so.npnts==0) {
                n_output_pnts_size = 10;
                set->so.pnts = (SvOutputPar *)g_malloc(n_output_pnts_size*sizeof(SvOutputPar));
            } else if (set->so.npnts >= n_output_pnts_size) { /*realloc*/
                n_output_pnts_size += 10;
                set->so.pnts = (SvOutputPar *)g_realloc(set->so.pnts, n_output_pnts_size*sizeof(SvOutputPar));
            }

            /*load parameters*/
            if (get_component(fr, &(set->so.pnts[set->so.npnts].component), key)) 
                goto exit;
            if (get_int(fr, &(set->so.pnts[set->so.npnts].skip), key)) 
                goto exit;
            if (get_int(fr, &(set->so.pnts[set->so.npnts].i), key)) 
                goto exit;
            if (get_int(fr, &(set->so.pnts[set->so.npnts].j), key)) 
                goto exit;
            fscanf(fr, "%s", (set->so.pnts[set->so.npnts].filebase));
            set->so.npnts++;
        } else if (strcmp(key, "OUT_IMAGE") == 0) {
            /*alloc or realloc the necessary structure*/
            if (set->so.nimgs==0) {
                n_output_imgs_size = 10;
                set->so.imgs = (SvOutputPar *)g_malloc(n_output_imgs_size*sizeof(SvOutputPar));
            } else if (set->so.nimgs >= n_output_imgs_size) { /*realloc*/
                n_output_imgs_size += 10;
                set->so.imgs = (SvOutputPar *)g_realloc(set->so.imgs, n_output_imgs_size*sizeof(SvOutputPar));
            }

            /*load parameters*/
            if (get_component(fr, &(set->so.imgs[set->so.nimgs].component), key)) 
                goto exit;
            if (get_int(fr, &(set->so.imgs[set->so.nimgs].skip), key)) 
                goto exit;
           
            fscanf(fr, "%100[^\n\r]", buffer);
            sprintf(set->so.imgs[set->so.nimgs].filebase, "%s", g_strstrip(g_strdup(buffer)));

            set->so.nimgs++;
        } else if (strcmp(key, "OUT_PLANE") == 0) {
            /*alloc or realloc the necessary structure*/
            if (set->so.nplns==0) {
                n_output_plns_size = 10;
                set->so.plns = (SvOutputPar *)g_malloc(n_output_plns_size*sizeof(SvOutputPar));
            } else if (set->so.nplns >= n_output_plns_size) { /*realloc*/
                n_output_plns_size += 10;
                set->so.plns = (SvOutputPar *)g_realloc(set->so.plns, n_output_plns_size*sizeof(SvOutputPar));
            }

            /*load parameters*/
            if (get_component(fr, &(set->so.plns[set->so.nplns].component), key)) 
                goto exit;
            if (get_int(fr, &(set->so.plns[set->so.nplns].skip), key)) 
                goto exit;
            if (get_int(fr, &(set->so.plns[set->so.nplns].start), key)) 
                goto exit;
            if (get_int(fr, &(set->so.plns[set->so.nplns].stop), key)) 
                goto exit;
            if (get_int(fr, &(set->so.plns[set->so.nplns].logscale), key)) 
                goto exit;
            
            fscanf(fr, "%100[^\n\r]", buffer);
            sprintf(set->so.plns[set->so.nplns].filebase, "%s", g_strstrip(g_strdup(buffer)));

            set->so.nplns++;
        } else if (strcmp(key, "OUT_SUM") == 0) {
            /*alloc or realloc the necessary structure*/
            if (set->so.nsums==0) {
                n_output_sums_size = 10;
                set->so.sums = (SvOutputSum *)g_malloc(n_output_sums_size*sizeof(SvOutputSum));
            } else if (set->so.nsums >= n_output_sums_size) { /*realloc*/
                n_output_sums_size += 10;
                set->so.sums = (SvOutputSum *)g_realloc(set->so.sums, n_output_sums_size*sizeof(SvOutputSum));
            }

            /*load parameters*/
            if (get_scomponent(fr, &(set->so.sums[set->so.nsums].component), key)) 
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].skip), key)) 
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].i0), key)) 
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].j0), key)) 
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].i1), key)) 
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].j1), key)) 
                goto exit;
            if (get_double(fr, &(set->so.sums[set->so.nsums].epsilon), key)) 
                goto exit;
            if (get_double(fr, &(set->so.sums[set->so.nsums].mu), key)) 
                goto exit;
            if (get_double(fr, &(set->so.sums[set->so.nsums].sigma), key)) 
                goto exit;
            if (get_double(fr, &(set->so.sums[set->so.nsums].sigast), key)) 
                goto exit;
            fscanf(fr, "%s", (set->so.sums[set->so.nsums].filename));
            set->so.sums[set->so.nsums].epsilon *= EPSILON_0;
            set->so.sums[set->so.nsums].mu *= MU_0;
            set->so.nsums++;
        } else if (strcmp(key, "MBOUNDARY_X0") == 0) {
            fscanf(fr, "%s", value);
            if (get_int(fr, &(set->smb.bx0pos), key)) 
                goto exit;
            if (strcmp(value, "periodic") == 0) {
                set->smb.bx0 = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        }
        else if (strcmp(key, "MBOUNDARY_XN") == 0) {
            fscanf(fr, "%s", value);
            if (get_int(fr, &(set->smb.bxnpos), key)) 
                goto exit;
            if (strcmp(value, "periodic") == 0) {
                set->smb.bxn = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "MBOUNDARY_Y0") == 0) {
            fscanf(fr, "%s", value);
            if (get_int(fr, &(set->smb.by0pos), key)) 
                goto exit;
            if (strcmp(value, "periodic") == 0) {
                set->smb.by0 = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        }
        else if (strcmp(key, "MBOUNDARY_YN") == 0) {
            fscanf(fr, "%s", value);
            if (get_int(fr, &(set->smb.bynpos), key)) 
                goto exit;
            if (strcmp(value, "periodic") == 0) {
                set->smb.byn = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "MBOUNDARY_ALL") == 0) {
            fscanf(fr, "%s", value);
            if (get_int(fr, &(set->smb.bx0pos), key)) 
                goto exit;
            if (get_int(fr, &(set->smb.bxnpos), key)) 
                goto exit;
            if (get_int(fr, &(set->smb.by0pos), key)) 
                goto exit;
            if (get_int(fr, &(set->smb.bynpos), key)) 
                goto exit;
            if (strcmp(value, "periodic") == 0) {
                set->smb.bx0 = set->smb.bxn = set->smb.by0 = set->smb.byn = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "BOUNDARY_ALL") == 0) {
            fscanf(fr, "%s", value);
            if (strcmp(value, "none") == 0) {
                set->sb.bx0 = set->sb.bxn = set->sb.by0 = set->sb.byn = SV_BOUNDARY_NONE;
            } else if (strcmp(value, "pec") == 0) {
                set->sb.bx0 = set->sb.bxn = set->sb.by0 = set->sb.byn = SV_BOUNDARY_PEC;
            } else if (strcmp(value, "liao") == 0) {
                set->sb.bx0 = set->sb.bxn = set->sb.by0 = set->sb.byn = SV_BOUNDARY_LIAO;
            } else {
                fprintf(stderr, "Error: unsupported boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "NFFF_SAVE") == 0) {
            if (get_int(fr, &(set->sf.i0), key)) 
                goto exit;
            if (get_int(fr, &(set->sf.j0), key)) 
                goto exit;
            if (get_int(fr, &(set->sf.i1), key)) 
                goto exit;
            if (get_int(fr, &(set->sf.j1), key)) 
                goto exit;
            fscanf(fr, "%s", buffer);
            set->sf.savefile = g_strdup(buffer);
            set->sf.save = 1;
            set->sf.ramahi = 0;
        } else if (strcmp(key, "NFFF_RSAVE") == 0) {
            if (get_int(fr, &(set->sf.i0), key)) 
                goto exit;
            if (get_int(fr, &(set->sf.j0), key)) 
                goto exit;
            if (get_int(fr, &(set->sf.i1), key)) 
                goto exit;
            if (get_int(fr, &(set->sf.j1), key)) 
                goto exit;
            fscanf(fr, "%s", buffer);
            set->sf.savefile = g_strdup(buffer);
            set->sf.save = 1;
            set->sf.ramahi = 1;
        }
         else if (strcmp(key, "NFFF_LOAD") == 0) {
            fscanf(fr, "%s", buffer);
            set->sf.savefile = g_strdup(buffer);
            set->sf.save = 0;
        }
      }

    if (set->smb.bxnpos == -1) 
        set->smb.bxnpos = set->sp.xres;
    if (set->smb.bynpos == -1) 
        set->smb.bynpos = set->sp.yres;
  
    return 0;

exit:
    fclose(fr);
    return 1;
}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
