
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
#include <libprocess/gwyprocess.h>
#include "settings.h"
#include "constants.h"
#include "plan.h"
#include "global.h"
#include <math.h>
#include <string.h>
#ifndef G_OS_WIN32
#include <sys/time.h>
#include <unistd.h>
#endif

typedef enum {
    SRC_DIRECT = 0,
    SRC_SINE = 1,
    SRC_PULSE = 2,
    SRC_BROADBAND = 3
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
    readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    return g_strdup(buf);
#else
#if (GLIB_CHECK_VERSION(2, 16, 0))
    return g_win32_get_package_installation_directory_of_module(NULL);
#else
    return g_win32_get_package_installation_directory(NULL, NULL);
#endif
#endif
    /*    LPWSTR buf_win;
        gchar* buf_utf8;
        buf_win = (LPWSTR)malloc (1024*sizeof(WCHAR));
        if (0 != GetModuleFileNameW(NULL, buf_win, 1024)) {
        //PathRemoveFileSpecW(buf_win);
        if ((buf_utf8 = g_utf16_to_utf8((gunichar2*)buf_win, -1, NULL, NULL, NULL)) != NULL)
        return g_strdup(buf_utf8);*/

    return NULL;
}

gint
get_self_dir_ex(gchar* buf, gint size)
{
#ifndef G_OS_WIN32
    readlink("/proc/self/exe", buf, size - 1);
    sprintf(buf, "%s", g_path_get_dirname(buf));
#else
#if (GLIB_CHECK_VERSION(2, 16, 0))
    sprintf(buf, "%s", g_win32_get_package_installation_directory_of_module(NULL));
#else
    sprintf(buf, "%s", g_win32_get_package_installation_directory(NULL, NULL));
#endif
#endif

    return 1;
}

/* read settings from file functions */

gint
get_int(FILE *F, gint *val, gchar *key)
{
    gchar value[100], *err;

    fscanf(F, "%s", value);
    *val = strtol(value, &err, 10);

    if (*err != '\0') {
        *val = (gint)strtod(value, &err);   /* this is to load integers in exponential format i.e. 2e7 */
        if (*err != '\0') {
            fprintf(stderr, "Error parsing integer value (key %s)\n", key);
            return 1;
        }
    }

    return 0;
}

gint
get_double(FILE *F, gdouble *val, gchar *key)
{
    gchar value[100], *err;

    fscanf(F, "%s", value);
    *val = strtod(value, &err);

    if (*err != '\0') {
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
    else if (strcmp(component, "AllFFT") == 0 || strcmp(component, "allfft") == 0)
        *val = SV_COMP_ALLFFT;
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
    else if (strcmp(component, "Material") == 0 || strcmp(component, "material") == 0)
        *val = SV_COMP_MAT;
    else {
        fprintf(stderr, "Error parsing field component (key %s)\n", key);
        return 1;
    }

    return 0;
}

gint
get_cube_component(FILE *fr, gint *val, gchar *key)
{
    gchar component[20];
    fscanf(fr, "%s", component);

    if (strcmp(component, "Ex") == 0 || strcmp(component, "ex") == 0)
        *val = SV_OVOLUME_EX;
    else if (strcmp(component, "Ey") == 0 || strcmp(component, "ey") == 0)
        *val = SV_OVOLUME_EY;
    else if (strcmp(component, "Ez") == 0 || strcmp(component, "ez") == 0)
        *val = SV_OVOLUME_EZ;
    else if (strcmp(component, "Hx") == 0 || strcmp(component, "hx") == 0)
        *val = SV_OVOLUME_HX;
    else if (strcmp(component, "Hy") == 0 || strcmp(component, "hy") == 0)
        *val = SV_OVOLUME_HY;
    else if (strcmp(component, "Hz") == 0 || strcmp(component, "hz") == 0)
        *val = SV_OVOLUME_HZ;
    else if (strcmp(component, "All") == 0 || strcmp(component, "all") == 0)
        *val = SV_OVOLUME_ALL;
    else if (strcmp(component, "Epsilon") == 0 || strcmp(component, "epsilon") == 0)
        *val = SV_OVOLUME_EPSILON;
    else if (strcmp(component, "Sigma") == 0 || strcmp(component, "sigma") == 0)
        *val = SV_OVOLUME_SIGMA;
    else if (strcmp(component, "Mu") == 0 || strcmp(component, "mu") == 0)
        *val = SV_OVOLUME_MU;
    else if (strcmp(component, "Sigast") == 0 || strcmp(component, "sigast") == 0)
        *val = SV_OVOLUME_SIGAST;
    else if (strcmp(component, "Material") == 0 || strcmp(component, "material") == 0)
        *val = SV_OVOLUME_MAT;
    else if (strcmp(component, "Mattype") == 0 || strcmp(component, "mattype") == 0)
        *val = SV_OVOLUME_MATTYPE;
    else if (strcmp(component, "Abs") == 0 || strcmp(component, "abs") == 0)
        *val = SV_OVOLUME_ABS;
    else if (strcmp(component, "Sumall") == 0 || strcmp(component, "sumall") == 0)
        *val = SV_OVOLUME_SUMALL;
    else if (strcmp(component, "Maxall") == 0 || strcmp(component, "maxall") == 0)
        *val = SV_OVOLUME_MAXALL;
    else {
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
    else if (strcmp(component, "Max") == 0 || strcmp(component, "max") == 0)
        *val = SV_SUM_MAX;
    else {
        fprintf(stderr, "Error parsing field component (key %s)\n", key);
        return 1;
    }

    return 0;
}

/* read settings from string functions */

gint
sget_int(const gchar *str, gint *val, gchar *key)
{
    gchar value[100], *err;

    sscanf(str, "%s", value);
    *val = strtol(value, &err, 10);

    if (*err != '\0') {
        *val = (gint)strtod(value, &err);   /* this is to load integers in exponential format i.e. 2e7 */
        if (*err != '\0') {
            fprintf(stderr, "Error parsing integer value (key %s)\n", key);
            return 1;
        }
    }

    return 0;
}

gint
sget_double(const gchar *str, gdouble *val, gchar *key)
{
    gchar value[100], *err;

    sscanf(str, "%s", value);
    *val = strtod(value, &err);

    if (*err != '\0') {
        fprintf(stderr, "Error parsing double value (key %s)\n", key);
        return 1;
    }

    return 0;
}

gint
sget_component(const gchar *str, gint *val, gchar *key)
{
    gchar component[20];
    sscanf(str, "%s", component);

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
    else if (strcmp(component, "AllFFT") == 0 || strcmp(component, "allfft") == 0)
        *val = SV_COMP_ALLFFT;
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
    else if (strcmp(component, "Material") == 0 || strcmp(component, "material") == 0)
        *val = SV_COMP_MAT;
    else {
        fprintf(stderr, "Error parsing field component (key %s)\n", key);
        return 1;
    }

    return 0;
}

gint
sget_cube_component(const gchar *str, gint *val, gchar *key)
{
    gchar component[20];
    sscanf(str, "%s", component);

    if (strcmp(component, "Ex") == 0 || strcmp(component, "ex") == 0)
        *val = SV_OVOLUME_EX;
    else if (strcmp(component, "Ey") == 0 || strcmp(component, "ey") == 0)
        *val = SV_OVOLUME_EY;
    else if (strcmp(component, "Ez") == 0 || strcmp(component, "ez") == 0)
        *val = SV_OVOLUME_EZ;
    else if (strcmp(component, "Hx") == 0 || strcmp(component, "hx") == 0)
        *val = SV_OVOLUME_HX;
    else if (strcmp(component, "Hy") == 0 || strcmp(component, "hy") == 0)
        *val = SV_OVOLUME_HY;
    else if (strcmp(component, "Hz") == 0 || strcmp(component, "hz") == 0)
        *val = SV_OVOLUME_HZ;
    else if (strcmp(component, "All") == 0 || strcmp(component, "all") == 0)
        *val = SV_OVOLUME_ALL;
    else if (strcmp(component, "Epsilon") == 0 || strcmp(component, "epsilon") == 0)
        *val = SV_OVOLUME_EPSILON;
    else if (strcmp(component, "Sigma") == 0 || strcmp(component, "sigma") == 0)
        *val = SV_OVOLUME_SIGMA;
    else if (strcmp(component, "Mu") == 0 || strcmp(component, "mu") == 0)
        *val = SV_OVOLUME_MU;
    else if (strcmp(component, "Sigast") == 0 || strcmp(component, "sigast") == 0)
        *val = SV_OVOLUME_SIGAST;
    else if (strcmp(component, "Material") == 0 || strcmp(component, "material") == 0)
        *val = SV_OVOLUME_MAT;
    else if (strcmp(component, "Mattype") == 0 || strcmp(component, "mattype") == 0)
        *val = SV_OVOLUME_MATTYPE;
    else if (strcmp(component, "Abs") == 0 || strcmp(component, "abs") == 0)
        *val = SV_OVOLUME_ABS;
    else if (strcmp(component, "Sumall") == 0 || strcmp(component, "sumall") == 0)
        *val = SV_OVOLUME_SUMALL;
    else if (strcmp(component, "Maxall") == 0 || strcmp(component, "maxall") == 0)
        *val = SV_OVOLUME_MAXALL;
    else {
        fprintf(stderr, "Error parsing field component (key %s)\n", key);
        return 1;
    }

    return 0;
}

gint
sget_scomponent(const gchar *str, gint *val, gchar *key)
{
    gchar component[5];
    sscanf(str, "%s", component);

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
    else if (strcmp(component, "Max") == 0 || strcmp(component, "max") == 0)
        *val = SV_SUM_MAX;
    else {
        fprintf(stderr, "Error parsing field component (key %s)\n", key);
        return 1;
    }

    return 0;
}


void generate_broadband_signal(gdouble *vals, gint nsteps, gdouble wavelength, gdouble wlspan, gdouble dt, gdouble amplitude)
{
    gint i, istart, istop, dampval = 100;
    gdouble phase;
    GwyDataLine *in_real = gwy_data_line_new(nsteps, nsteps, 1);
    GwyDataLine *in_im = gwy_data_line_new(nsteps, nsteps, 1);
    GwyDataLine *out_real = gwy_data_line_new(nsteps, nsteps, 1);
    GwyDataLine *out_im = gwy_data_line_new(nsteps, nsteps, 1);
    GRand *rnd;


    istop = (gint)(LIGHT_SPEED*nsteps*dt / (wavelength - wlspan / 2));
    istart = (gint)(LIGHT_SPEED*nsteps*dt / (wavelength + wlspan / 2));

    rnd = g_rand_new();
    g_rand_set_seed(rnd, 1);
    for (i = istart; i < istop; i++) {
        phase = g_rand_double_range(rnd, 0, G_PI);
        gwy_data_line_set_val(in_real, i, cos(phase));
        gwy_data_line_set_val(in_im, i, sin(phase));
    }

    gwy_data_line_fft(in_real, in_im, out_real, out_im,
                      GWY_WINDOWING_NONE,
                      GWY_TRANSFORM_DIRECTION_BACKWARD,
                      GWY_INTERPOLATION_BILINEAR,
                      TRUE, 2);

    for (i = 0; i < nsteps; i++) {
        vals[i] = gwy_data_line_get_val(out_real, i);
        if (i < dampval) vals[i] *= (1.0 / (gdouble)dampval) * ((gdouble)(i));
    }

    gwy_object_unref(in_real);
    gwy_object_unref(in_im);
    gwy_object_unref(out_real);
    gwy_object_unref(out_im);

    g_rand_free(rnd);
}

gint
write_source(SrcFormat format, SrcType type, gint nsteps, gdouble wavelength, gdouble wlspan, gint pulse_width, gdouble dt,
             gint component, gchar *filename, gdouble amplitude, gdouble theta, gdouble phi)
{
    gint i, j, dampval = 100, shiftval = 1000;
    gdouble *vals;
    gdouble wavefreq, damp, ax, ay, az;
    gdouble dpulse_width;

    FILE *fw = fopen(filename, "w");

    if (!fw) {
        fprintf(stderr, "Error: cannot write source file: %s\n", filename);
        return 1;
    }

    if (format == SRCF_LTSF) {
        dt /= 3;
        nsteps *= 3;
        pulse_width *= 3;
        dampval *= 3;
    }

    wavefreq = LIGHT_SPEED * 2 * G_PI / wavelength;
    if (type == SRC_PULSE) 
        dpulse_width = pulse_width * dt;


    vals = (gdouble *)g_malloc(nsteps * sizeof(gdouble));

    if (type == SRC_BROADBAND)
        generate_broadband_signal(vals, nsteps, wavelength, wlspan, dt, amplitude);
    else {
        for (i = 0; i < nsteps; i++) {
            if (type == SRC_BROADBAND) {

                /*
                   if (i<(period/2))
                   wl = wavelength - wlspan/2 + (double)i*wlspan/(double)(period/2);
                   else if (i<period)
                   wl = wavelength - wlspan/2 + (double)(period - i)*wlspan/(double)(period/2);
                   else amplitude = 0;

                   /            printf("wl: %d %g\n", i, wl);

                   wavefreq = LIGHT_SPEED*2*G_PI/wl;
                   val = amplitude * sin(wavefreq*dt*(gdouble)i);
                 */

                vals[i] = 0;
                shiftval = 100000;
                for (j = 0; j < 50; j++) {
                    wavefreq = LIGHT_SPEED * 2 * G_PI / (wavelength - wlspan / 2 + (double)j*wlspan / 50.0);
                    vals[i] += 0.1 * amplitude * sin(wavefreq*dt*(gdouble)(i + shiftval));
                }

            } else {
                vals[i] = amplitude * sin(wavefreq*dt*(gdouble)i); //was only 5, 5e7 for particle
            }

            if (i < dampval)
                vals[i] *= (1.0 / (gdouble)dampval) * ((gdouble)(i));

            if (type == SRC_PULSE) {
                damp = ((gdouble)i*dt - 4.0*dpulse_width) / sqrt((gdouble)2.0) / dpulse_width;
                vals[i] *= exp(-damp*damp);
            }
        }
    }


    fprintf(fw, "%d\n", nsteps);
    for (i = 0; i < nsteps; i++) {
        if (format == SRCF_FULL) {
            if (component == SV_COMP_EX)
                fprintf(fw, "%d %g %g %g %g %g %g\n", i, vals[i], 0.0, 0.0, 0.0, 0.0, 0.0);
            else if (component == SV_COMP_EY)
                fprintf(fw, "%d %g %g %g %g %g %g\n", i, 0.0, vals[i], 0.0, 0.0, 0.0, 0.0);
            else if (component == SV_COMP_EZ)
                fprintf(fw, "%d %g %g %g %g %g %g\n", i, 0.0, 0.0, vals[i], 0.0, 0.0, 0.0);
            else if (component == SV_COMP_HX)
                fprintf(fw, "%d %g %g %g %g %g %g\n", i, 0.0, 0.0, 0.0, vals[i], 0.0, 0.0);
            else if (component == SV_COMP_HY)
                fprintf(fw, "%d %g %g %g %g %g %g\n", i, 0.0, 0.0, 0.0, 0.0, vals[i], 0.0);
            else if (component == SV_COMP_HZ)
                fprintf(fw, "%d %g %g %g %g %g %g\n", i, 0.0, 0.0, 0.0, 0.0, 0.0, vals[i]);
        }
        if (format == SRCF_ARB) {
            ax = sin(theta)*cos(phi);
            ay = sin(theta)*sin(phi);
            az = cos(theta);

            fprintf(fw, "%d %g %g %g %g %g %g\n", i, ax*vals[i], ay*vals[i], az*vals[i], 0.0, 0.0, 0.0);
        } else
            fprintf(fw, "%d %g\n", i, vals[i]);
    }
    fclose(fw);

    g_free(vals);

    // if (format==SRCF_TSF) 
    //    fprintf(fw, "%d\n");

    return 0;
}

void
delete_point_sources(SvSet *set)
{
    int i;

    if ((set->ss.npnts_allocated) > 0 && (NULL != set->ss.pnts)) {
        for (i = 0; i < set->ss.npnts; i++) {
            g_free(set->ss.pnts[i].source_filename);
            set->ss.pnts[i].source_filename = NULL;
            set->ss.npnts = 0;
            set->ss.pnts = NULL;
            set->ss.npnts_allocated = 0;
        }

        g_free(set->ss.pnts);
        set->ss.pnts = NULL;
    }
}

void
clear_settings(SvSet *set, gboolean remove)
{
    gint i;

#ifndef G_OS_WIN32
    struct timeval time;
#endif

    /* delete settings */
    if (remove)
        delete_point_sources(set);

    /* inits */

    /* source points - no sources by default*/
    set->ss.npnts = 0;
    set->ss.pnts = NULL;
    set->ss.npnts_allocated = 0;
    set->ss.nlocals = 0;

    /* SF */
    set->ss.sf.layered_epsilon = SOURCE_LAYERED_EPS;
    set->ss.sf.layered_mu = SOURCE_LAYERED_MU;
    set->ss.sf.layered_sigma = SOURCE_LAYERED_SIGMA;
    set->ss.sf.layered_sigast = SOURCE_LAYERED_SIGAST;

    if (remove)
        g_free(set->ss.sf.source_filename);
    set->ss.sf.source_filename = NULL;
    set->ss.sf.source_mode = SOURCE_MODE;
    set->ss.sf.source_amplitude = SOURCE_AMPLITUDE;
    set->ss.sf.source_wl = SOURCE_WAVELENGTH;
    set->ss.sf.source_wlspan = SOURCE_WAVELENGTHSPAN;
    set->ss.sf.source_pulsewidth = SOURCE_PULSE_WIDTH;

    set->ss.sf.ia_theta = SOURCE_IA_THETADEG;
    set->ss.sf.ia_phi = SOURCE_IA_PHIDEG;
    set->ss.sf.ia_psi = SOURCE_IA_PSIDEG;

    set->ss.sf.valid = FALSE;

    /* TSF */
    set->ss.tsf.box_i0 = set->ss.tsf.box_j0 = set->ss.tsf.box_k0 = SOURCE_BOX_VERTEX0;
    set->ss.tsf.box_in = set->ss.tsf.box_jn = set->ss.tsf.box_kn = SOURCE_BOX_VERTEXN;

    set->ss.tsf.layered_epsilon = SOURCE_LAYERED_EPS;
    set->ss.tsf.layered_mu = SOURCE_LAYERED_MU;
    set->ss.tsf.layered_sigma = SOURCE_LAYERED_SIGMA;
    set->ss.tsf.layered_sigast = SOURCE_LAYERED_SIGAST;

    set->ss.tsf.gaussian = SOURCE_GAUSSIAN_MULT_ENABLE;
    set->ss.tsf.gaussian_fxpos = SOURCE_XXX_MULT_CENTER;
    set->ss.tsf.gaussian_fypos = SOURCE_XXX_MULT_CENTER;
    set->ss.tsf.gaussian_rx = SOURCE_XXX_MULT_RADIUS;
    set->ss.tsf.gaussian_ry = SOURCE_XXX_MULT_RADIUS;

    set->ss.tsf.radial = SOURCE_RADIAL_MULT_ENABLE;
    set->ss.tsf.radial_fxpos = SOURCE_XXX_MULT_CENTER;
    set->ss.tsf.radial_fypos = SOURCE_XXX_MULT_CENTER;
    set->ss.tsf.radial_rx = SOURCE_XXX_MULT_RADIUS;
    set->ss.tsf.radial_ry = SOURCE_XXX_MULT_RADIUS;

    set->ss.tsf.fiber = SOURCE_FIBER_MULT_ENABLE;
    set->ss.tsf.fiber_fxpos = SOURCE_XXX_MULT_CENTER;
    set->ss.tsf.fiber_fypos = SOURCE_XXX_MULT_CENTER;
    set->ss.tsf.fiber_radius = SOURCE_XXX_MULT_RADIUS;
    set->ss.tsf.fiber_cutoff = SOURCE_FIBER_MULT_CUTOFF;
    set->ss.tsf.fiber_epsilon_core = SOURCE_FIBER_MULT_EPSILON_CORE;
    set->ss.tsf.fiber_epsilon_cladding = SOURCE_FIBER_MULT_EPSILON_CLADDING;

    if (remove)
        g_free(set->ss.tsf.source_filename);
    set->ss.tsf.source_filename = NULL;
    set->ss.tsf.source_mode = SOURCE_MODE;
    set->ss.tsf.source_wl = SOURCE_WAVELENGTH;
    set->ss.tsf.source_wlspan = SOURCE_WAVELENGTHSPAN;
    set->ss.tsf.source_amplitude = SOURCE_AMPLITUDE;
    set->ss.tsf.source_pulsewidth = SOURCE_PULSE_WIDTH;

    set->ss.tsf.box_boundary_skipi0 = set->ss.tsf.box_boundary_skipj0 = set->ss.tsf.box_boundary_skipk0 = SOURCE_BOUNDARY_SKIP;
    set->ss.tsf.box_boundary_skipin = set->ss.tsf.box_boundary_skipjn = set->ss.tsf.box_boundary_skipkn = SOURCE_BOUNDARY_SKIP;
    set->ss.tsf.box_boundary_skipdepth = SOURCE_BOUNDARY_SKIP_DEPTH;

    set->ss.tsf.valid = SOURCE_TSF_VALID;
    /* TSF */

    /* LTSF */
    set->ss.ltsf.box_i0 = set->ss.ltsf.box_j0 = set->ss.ltsf.box_k0 = SOURCE_BOX_VERTEX0;
    set->ss.ltsf.box_in = set->ss.ltsf.box_jn = set->ss.ltsf.box_kn = SOURCE_BOX_VERTEXN;

    if (remove) {
        for(i=0; i < set->ss.ltsf.layered_count; i++)
            g_free(set->ss.ltsf.layered_material[i]);
    }
    memset(set->ss.ltsf.layered_zpos, 0, sizeof(set->ss.ltsf.layered_zpos));
    memset(set->ss.ltsf.layered_epsilon, 0, sizeof(set->ss.ltsf.layered_epsilon));
    memset(set->ss.ltsf.layered_mu, 0, sizeof(set->ss.ltsf.layered_mu));
    memset(set->ss.ltsf.layered_sigma, 0, sizeof(set->ss.ltsf.layered_sigma));
    memset(set->ss.ltsf.layered_sigast, 0, sizeof(set->ss.ltsf.layered_sigast));
    memset(set->ss.ltsf.layered_material, 0, sizeof(set->ss.ltsf.layered_material));
    set->ss.ltsf.layered_count = 1;
    set->ss.ltsf.layered_zpos[0] = SOURCE_LAYERED_ZPOS;
    set->ss.ltsf.layered_epsilon[0] = SOURCE_LAYERED_EPS;
    set->ss.ltsf.layered_mu[0] = SOURCE_LAYERED_MU;
    set->ss.ltsf.layered_sigma[0] = SOURCE_LAYERED_SIGMA;
    set->ss.ltsf.layered_sigast[0] = SOURCE_LAYERED_SIGAST;
    set->ss.ltsf.layered_material[0] = g_strdup(SOURCE_LAYERED_MATERIAL);

    set->ss.ltsf.gaussian = SOURCE_GAUSSIAN_MULT_ENABLE;
    set->ss.ltsf.gaussian_fxpos = SOURCE_XXX_MULT_CENTER;
    set->ss.ltsf.gaussian_fypos = SOURCE_XXX_MULT_CENTER;
    set->ss.ltsf.gaussian_rx = SOURCE_XXX_MULT_RADIUS;
    set->ss.ltsf.gaussian_ry = SOURCE_XXX_MULT_RADIUS;

    set->ss.ltsf.radial = SOURCE_RADIAL_MULT_ENABLE;
    set->ss.ltsf.radial_fxpos = SOURCE_XXX_MULT_CENTER;
    set->ss.ltsf.radial_fypos = SOURCE_XXX_MULT_CENTER;
    set->ss.ltsf.radial_rx = SOURCE_XXX_MULT_RADIUS;
    set->ss.ltsf.radial_ry = SOURCE_XXX_MULT_RADIUS;

    set->ss.ltsf.fiber = SOURCE_FIBER_MULT_ENABLE;
    set->ss.ltsf.fiber_fxpos = SOURCE_XXX_MULT_CENTER;
    set->ss.ltsf.fiber_fypos = SOURCE_XXX_MULT_CENTER;
    set->ss.ltsf.fiber_radius = SOURCE_XXX_MULT_RADIUS;
    set->ss.ltsf.fiber_cutoff = SOURCE_FIBER_MULT_CUTOFF;
    set->ss.ltsf.fiber_epsilon_core = SOURCE_FIBER_MULT_EPSILON_CORE;
    set->ss.ltsf.fiber_epsilon_cladding = SOURCE_FIBER_MULT_EPSILON_CLADDING;

    if (remove)
        g_free(set->ss.ltsf.source_filename);
    set->ss.ltsf.source_filename = NULL;

    set->ss.ltsf.box_boundary_skipi0 = set->ss.ltsf.box_boundary_skipj0 = set->ss.ltsf.box_boundary_skipk0 = SOURCE_BOUNDARY_SKIP;
    set->ss.ltsf.box_boundary_skipin = set->ss.ltsf.box_boundary_skipjn = set->ss.ltsf.box_boundary_skipkn = SOURCE_BOUNDARY_SKIP;
    set->ss.ltsf.box_boundary_skipdepth = SOURCE_BOUNDARY_SKIP_DEPTH;

    set->ss.ltsf.valid = SOURCE_LTSF_VALID;
    /* LTSF */

    /* TSFF */
    set->ss.tsff.box_i0 = set->ss.tsff.box_j0 = set->ss.tsff.box_k0 = SOURCE_BOX_VERTEX0;
    set->ss.tsff.box_in = set->ss.tsff.box_jn = set->ss.tsff.box_kn = SOURCE_BOX_VERTEXN;
    set->ss.tsff.xshift = SOURCE_TSFF_SHIFT;
    set->ss.tsff.yshift = SOURCE_TSFF_SHIFT;
    set->ss.tsff.zshift = SOURCE_TSFF_SHIFT;
    set->ss.tsff.source_mode = SOURCE_MODE;
    set->ss.tsff.source_wl = SOURCE_WAVELENGTH;
    set->ss.tsff.source_wlspan = SOURCE_WAVELENGTHSPAN;

    if (remove)
        g_free(set->ss.tsff.source_filename);
    set->ss.tsff.source_filename = NULL;

    set->ss.tsff.box_boundary_skipi0 = set->ss.tsff.box_boundary_skipj0 = set->ss.tsff.box_boundary_skipk0 = SOURCE_BOUNDARY_SKIP;
    set->ss.tsff.box_boundary_skipin = set->ss.tsff.box_boundary_skipjn = set->ss.tsff.box_boundary_skipkn = SOURCE_BOUNDARY_SKIP;
    set->ss.tsff.box_boundary_skipdepth = SOURCE_BOUNDARY_SKIP_DEPTH;

    set->ss.tsff.valid = SOURCE_TSFF_VALID;
    /* TSFF */

    /* LTSFF */
    set->ss.ltsff.box_i0 = set->ss.ltsff.box_j0 = set->ss.ltsff.box_k0 = SOURCE_BOX_VERTEX0;
    set->ss.ltsff.box_in = set->ss.ltsff.box_jn = set->ss.ltsff.box_kn = SOURCE_BOX_VERTEXN;

    if (remove) {
        for (i = 0; i < set->ss.ltsff.layered_count; i++)
            g_free(set->ss.ltsff.layered_material[i]);
    }

    memset(set->ss.ltsff.layered_zpos, 0, sizeof(set->ss.ltsff.layered_zpos));
    memset(set->ss.ltsff.layered_epsilon, 0, sizeof(set->ss.ltsff.layered_epsilon));
    memset(set->ss.ltsff.layered_mu, 0, sizeof(set->ss.ltsff.layered_mu));
    memset(set->ss.ltsff.layered_sigma, 0, sizeof(set->ss.ltsff.layered_sigma));
    memset(set->ss.ltsff.layered_sigast, 0, sizeof(set->ss.ltsff.layered_sigast));
    memset(set->ss.ltsff.layered_material, 0, sizeof(set->ss.ltsff.layered_material));
    set->ss.ltsff.layered_count = 1;
    set->ss.ltsff.layered_zpos[0] = SOURCE_LAYERED_ZPOS;
    set->ss.ltsff.layered_epsilon[0] = SOURCE_LAYERED_EPS;
    set->ss.ltsff.layered_mu[0] = SOURCE_LAYERED_MU;
    set->ss.ltsff.layered_sigma[0] = SOURCE_LAYERED_SIGMA;
    set->ss.ltsff.layered_sigast[0] = SOURCE_LAYERED_SIGAST;
    set->ss.ltsff.layered_material[0] = g_strdup(SOURCE_LAYERED_MATERIAL);

    if (remove)
        g_free(set->ss.ltsff.source_filename);
    set->ss.ltsff.source_filename = NULL;
    set->ss.ltsff.fast = 0;

    set->ss.ltsff.box_boundary_skipi0 = set->ss.ltsff.box_boundary_skipj0 = set->ss.ltsff.box_boundary_skipk0 = SOURCE_BOUNDARY_SKIP;
    set->ss.ltsff.box_boundary_skipin = set->ss.ltsff.box_boundary_skipjn = set->ss.ltsff.box_boundary_skipkn = SOURCE_BOUNDARY_SKIP;
    set->ss.ltsff.box_boundary_skipdepth = SOURCE_BOUNDARY_SKIP_DEPTH;

    set->ss.ltsff.valid = SOURCE_LTSFF_VALID;
    /* LTSFF */


    set->ss.lambda_min = set->ss.lambda_max = set->ss.lambda_center = -1;

    set->ss.ext.filebase_ex = NULL;
    set->ss.ext.filebase_ey = NULL;
    set->ss.ext.filebase_ez = NULL;
    set->ss.ext.filebase_hx = NULL;
    set->ss.ext.filebase_hy = NULL;
    set->ss.ext.filebase_hz = NULL;

    /*no output by default*/
    set->so.npnts = 0;
    set->so.npnts_allocated = 0;
    set->so.pnts = NULL;

    set->so.nplns = 0;
    set->so.nplns_allocated = 0;
    set->so.plns = NULL;

    set->so.nsums = 0;
    set->so.nsums_allocated = 0;
    set->so.sums = NULL;

    set->so.nimgs = 0;
    set->so.nimgs_allocated = 0;
    set->so.imgs = NULL;

    set->so.nsubgrid_imgs = 0;
    set->so.nsubgrid_imgs_allocated = 0;
    set->so.subgrid_imgs = NULL;

    set->so.ncubs = 0;
    set->so.ncubs_allocated = 0;
    set->so.cubs = NULL;

    set->so.nforces = 0;
    set->so.nforces_allocated = 0;
    set->so.forces = NULL;

    set->so.savespectrum = 0;

    set->sc.verbose = 0;
    set->sc.step_act = 0;

    /*no material data by default*/
    set->sm.in_voxel_filename = NULL;
    set->sm.in_vector_filename = NULL;
    set->sm.in_vtk_data_filename = NULL;

    /*use full material properties algorithm by default*/
    set->sm.matmode_check = 0;

    /*do not try to smooth materials by default*/
    set->sm.smooth = 0;

    /*no grow modifier by default*/
    set->sm.ngrowths = 0;
    for (i = 0; i < MATERIAL_COUNT; i++) {
        set->sm.grow_i0[i] = set->sm.grow_j0[i] = set->sm.grow_k0[i] = set->sm.grow_in[i] = set->sm.grow_jn[i] = set->sm.grow_kn[i] = 0;
        set->sm.grow_addindex[i] = set->sm.grow_attachindex[i] = set->sm.grow_subsampling[i] = set->sm.grow_nsteps[i] = set->sm.grow_seed[i] = 0;
        set->sm.grow_mobility[i] = set->sm.grow_probability[i] = 0.0f;
        set->sm.grow_skipi0[i] = set->sm.grow_skipin[i] = set->sm.grow_skipj0[i] = set->sm.grow_skipjn[i] = set->sm.grow_skipk0[i] = set->sm.grow_skipkn[i] = MATERIAL_GROW_SKIP;
    }

    /*no roughness modifier by default*/
    set->sm.nroughens = 0;   
    for (i = 0; i < MATERIAL_COUNT; i++) {
        set->sm.rough_radius_peak[i] = set->sm.rough_radius_span[i] = set->sm.rough_iterations[i] = 0;
        set->sm.rough_probability[i] = 0.0f;
        set->sm.rough_matindex[i] = set->sm.rough_voidindex[i] = set->sm.rough_seed[i] = 0;
    }

    /*no spectral modifier by default*/
    set->sm.nspectrals = 0;
    for (i = 0; i < MATERIAL_COUNT; i++) {
        set->sm.spectral_sigma[i] = set->sm.spectral_t[i] = set->sm.spectral_seed[i] = set->sm.spectral_matindex[i] = 0;
    }

    /*no expression modifier by default*/
    set->sm.nexprs = 0;
    for (i = 0; i < MATERIAL_COUNT; i++) {
        set->sm.expr_i0[i] = set->sm.expr_j0[i] = set->sm.expr_k0[i] = set->sm.expr_in[i] = set->sm.expr_jn[i] = set->sm.expr_kn[i] = 0;
        set->sm.expr_matindex[i] = set->sm.expr_voidindex[i] = set->sm.expr_maxdist[i] = set->sm.expr_distmode[i] = 0;
        memset(set->sm.expr_expr, 0, sizeof(gchar)*MATERIAL_EXPR_EXPR_CHARS);
    }

    /*no mathematic expression modifier by default*/
    set->sm.nexprs = 0;

    /*no material crop by default*/
    set->sm.crop = FALSE;    


    /*use local files for materials if there are any*/
    set->sm.localfiles = TRUE;

    /*none boundary condition by default (including material boundary)*/
    set->sb.bx0 = set->sb.bxn = SV_BOUNDARY_NONE;
    set->sb.by0 = set->sb.byn = SV_BOUNDARY_NONE;
    set->sb.bz0 = set->sb.bzn = SV_BOUNDARY_NONE;
    set->smb.bx0 = set->smb.bxn = SV_BOUNDARY_NONE;
    set->smb.by0 = set->smb.byn = SV_BOUNDARY_NONE;
    set->smb.bz0 = set->smb.bzn = SV_BOUNDARY_NONE;
    set->smb.bx0pos = set->smb.by0pos = set->smb.bz0pos = 0;
    set->smb.bxnpos = set->smb.bynpos = set->smb.bznpos = -1;

    set->sb.depth_bx0 = 10;
    set->sb.m_bx0 = 3;
    set->sb.sigma_bx0 = -1;
    set->sb.a_bx0 = 0.03;
    set->sb.kappa_bx0 = 4;

    set->sb.depth_bxn = 10;
    set->sb.m_bxn = 3;
    set->sb.sigma_bxn = -1;
    set->sb.a_bxn = 0.03;
    set->sb.kappa_bxn = 4;

    set->sb.depth_by0 = 10;
    set->sb.m_by0 = 3;
    set->sb.sigma_by0 = -1;
    set->sb.a_by0 = 0.03;
    set->sb.kappa_by0 = 4;

    set->sb.depth_byn = 10;
    set->sb.m_byn = 3;
    set->sb.sigma_byn = -1;
    set->sb.a_byn = 0.03;
    set->sb.kappa_byn = 4;

    set->sb.depth_bz0 = 10;
    set->sb.m_bz0 = 3;
    set->sb.sigma_bz0 = -1;
    set->sb.a_bz0 = 0.03;
    set->sb.kappa_bz0 = 4;

    set->sb.depth_bzn = 10;
    set->sb.m_bzn = 3;
    set->sb.sigma_bzn = -1;
    set->sb.a_bzn = 0.03;
    set->sb.kappa_bzn = 4;

    /*no farfield points by default*/
    set->sf.nrs = 0;
    set->sf.ri = set->sf.rj = set->sf.rk = NULL;
    set->sf.box_i0 = set->sf.box_j0 = set->sf.box_k0 = NFFF_BOX_VERTEX0;
    set->sf.box_in = set->sf.box_jn = set->sf.box_kn = NFFF_BOX_VERTEXN;
    set->sf.box_boundary_skipi0 = set->sf.box_boundary_skipj0 = set->sf.box_boundary_skipk0 = NFFF_BOX_BOUNDARY_SKIP;
    set->sf.box_boundary_skipin = set->sf.box_boundary_skipjn = set->sf.box_boundary_skipkn = NFFF_BOX_BOUNDARY_SKIP;
    set->sf.skipi0_jmin = NFFF_SKIPI0;
    set->sf.skipi0_kmin = NFFF_SKIPI0;
    set->sf.skipi0_jmax = NFFF_SKIPI0;
    set->sf.skipi0_kmax = NFFF_SKIPI0;
    set->sf.skipin_jmin = NFFF_SKIPIN;
    set->sf.skipin_kmin = NFFF_SKIPIN;
    set->sf.skipin_jmax = NFFF_SKIPIN;
    set->sf.skipin_kmax = NFFF_SKIPIN;
    set->sf.skipj0_imin = NFFF_SKIPJ0;
    set->sf.skipj0_kmin = NFFF_SKIPJ0;
    set->sf.skipj0_imax = NFFF_SKIPJ0;
    set->sf.skipj0_kmax = NFFF_SKIPJ0;
    set->sf.skipjn_imin = NFFF_SKIPJN;
    set->sf.skipjn_kmin = NFFF_SKIPJN;
    set->sf.skipjn_imax = NFFF_SKIPJN;
    set->sf.skipjn_kmax = NFFF_SKIPJN;
    set->sf.skipk0_imin = NFFF_SKIPK0;
    set->sf.skipk0_jmin = NFFF_SKIPK0;
    set->sf.skipk0_imax = NFFF_SKIPK0;
    set->sf.skipk0_jmax = NFFF_SKIPK0;
    set->sf.skipkn_imin = NFFF_SKIPKN;
    set->sf.skipkn_jmin = NFFF_SKIPKN;
    set->sf.skipkn_imax = NFFF_SKIPKN;
    set->sf.skipkn_jmax = NFFF_SKIPKN;

    set->sf.nsets = 0;
    set->sf.nareas = 0;
    set->sf.nsquares = 0;

    //set->sf.valid = FALSE;

    set->spf.nrs = 0;
    set->spf.postprocess = 0;
    set->spf.ppstart = 0;
    set->spf.ri = set->spf.rj = set->spf.rk = NULL;
    set->spf.box_i0 = set->spf.box_in = set->spf.box_j0 = PNFFF_BOX_VERTEX0;
    set->spf.box_jn = set->spf.box_k0 = set->spf.box_kn = PNFFF_BOX_VERTEXN;
    set->spf.box_boundary_skipk0 = set->spf.box_boundary_skipkn = PNFFF_BOX_BOUNDARY_SKIP;
    set->spf.skipk0_imin = set->spf.skipk0_jmin = set->spf.skipk0_imax = set->spf.skipk0_jmax = PNFFF_SKIPK0;
    set->spf.skipkn_imin = set->spf.skipkn_jmin = set->spf.skipkn_imax = set->spf.skipkn_jmax = PNFFF_SKIPKN;
    set->spf.pimin = PNFFF_INTEGRATION_XMIN;
    set->spf.pimax = PNFFF_INTEGRATION_XMAX;
    set->spf.pjmin = PNFFF_INTEGRATION_YMIN;    
    set->spf.pjmax = PNFFF_INTEGRATION_YMAX;

    set->spf.nsets = 0;
    set->spf.nareas = 0;

    //set->spf.valid = FALSE;

    set->sc.nsteps = BASIC_STEPS;
    //set->sc.suffix = 0;
    set->sc.dtmult = 1;
    /*do not use GPU by default*/
    set->sc.usegpu = BASIC_USEGPU;
    /*maximum threads by default*/
    set->sc.nthreads = BASIC_NTHREADS;
    set->so.outfile = g_strdup_printf(BASIC_OUTFILE);
    for (i = 0; i < MAX_GPUS; i++)
        set->sc.ugpu[i] = BASIC_UGPUINDEX;
    /*run devicequery for gpus by default*/
    set->sc.devicequery = 1;

    set->sg.nsg = 0;

    GRand *rnd;
#ifdef G_OS_WIN32
    rnd = g_rand_new();
    g_rand_set_seed(rnd, g_random_int() & 0x7fffffff);
#else
    gettimeofday(&time, NULL);
    rnd = g_rand_new_with_seed((time.tv_sec * 1000) + (time.tv_usec / 1000));
#endif

    set->sc.suffix = g_rand_int_range(rnd, 0, 999999);
} /* clear_settings */

void
clear_settings_mat(SvSetMat *set, gboolean remove)
{
    if (remove) {
        if (set->spheres)
            g_array_free(set->spheres, TRUE);
        if (set->cones)
            g_array_free(set->cones, TRUE);
        if (set->rcones)
            g_array_free(set->rcones, TRUE);
        if (set->cylinders)
            g_array_free(set->cylinders, TRUE);
        if (set->voxels)
            g_array_free(set->voxels, TRUE);
        if (set->tetrahedrons)
            g_array_free(set->tetrahedrons, TRUE);
        if (set->gwydds)
            g_array_free(set->gwydds, TRUE);
    }

    memset(set, 0, sizeof(SvSetMat));
}


void
source_point_alloc(SvSet *set)
{
    gint i;

    /*alloc or realloc the necessary structure*/
    if (set->ss.npnts_allocated == 0) {
        set->ss.npnts_allocated = 10;
        set->ss.pnts = (SvSrcPoint *)g_malloc(set->ss.npnts_allocated * sizeof(SvSrcPoint));
    } else if (set->ss.npnts >= set->ss.npnts_allocated) { /*realloc*/
        set->ss.npnts_allocated += 10;
        set->ss.pnts = (SvSrcPoint *)g_realloc(set->ss.pnts, set->ss.npnts_allocated * sizeof(SvSrcPoint));
    }

    for (i = set->ss.npnts; i < set->ss.npnts_allocated; i++) {
        set->ss.pnts[i].point_origin_position_i = 50;
        set->ss.pnts[i].point_origin_position_j = 50;
        set->ss.pnts[i].point_origin_position_k = 50;
        set->ss.pnts[i].point_origin_theta = 0.0f;
        set->ss.pnts[i].point_origin_phi = 0.0f;

        set->ss.pnts[i].source_mode = SOURCE_MODE;
        set->ss.pnts[i].source_filename = NULL;
        set->ss.pnts[i].source_wl = SOURCE_WAVELENGTH;
        set->ss.pnts[i].source_wlspan = SOURCE_WAVELENGTHSPAN;
        set->ss.pnts[i].source_amplitude = SOURCE_AMPLITUDE;
        set->ss.pnts[i].source_pulsewidth = SOURCE_PULSE_WIDTH;
    }
}


void
output_alloc(SvSet *set, SvType type)
{
    /*alloc or realloc the necessary structure*/

    if (SV_OUTTYPE_POINT == type) {
        if (set->so.npnts_allocated == 0) {
            set->so.npnts_allocated = 10;
            set->so.pnts = (SvOutputPar *)g_malloc(set->so.npnts_allocated * sizeof(SvOutputPar));
        } else if (set->so.npnts >= set->so.npnts_allocated) { /*realloc*/
            set->so.npnts_allocated += 10;
            set->so.pnts = (SvOutputPar *)g_realloc(set->so.pnts, set->so.npnts_allocated * sizeof(SvOutputPar));
        }
    }

    if (SV_OUTTYPE_IMAGE == type) {
        if (set->so.nimgs_allocated == 0) {
            set->so.nimgs_allocated = 10;
            set->so.imgs = (SvOutputPar *)g_malloc(set->so.nimgs_allocated * sizeof(SvOutputPar));
        } else if (set->so.nimgs >= set->so.nimgs_allocated) { /*realloc*/
            set->so.nimgs_allocated += 10;
            set->so.imgs = (SvOutputPar *)g_realloc(set->so.imgs, set->so.nimgs_allocated * sizeof(SvOutputPar));
        }
    }

    if (SV_OUTTYPE_SUBGRID_IMAGE == type) {
        if (set->so.nsubgrid_imgs_allocated == 0) {
            set->so.nsubgrid_imgs_allocated = 10;
            set->so.subgrid_imgs = (SvOutputPar *)g_malloc(set->so.nsubgrid_imgs_allocated * sizeof(SvOutputPar));
        } else if (set->so.nsubgrid_imgs >= set->so.nsubgrid_imgs_allocated) { /*realloc*/
            set->so.nsubgrid_imgs_allocated += 10;
            set->so.subgrid_imgs = (SvOutputPar *)g_realloc(set->so.subgrid_imgs, set->so.nsubgrid_imgs_allocated * sizeof(SvOutputPar));
        }
    }

    if (SV_OUTTYPE_PLANE == type) {
        if (set->so.nplns_allocated == 0) {
            set->so.nplns_allocated = 10;
            set->so.plns = (SvOutputPar *)g_malloc(set->so.nplns_allocated * sizeof(SvOutputPar));
        } else if (set->so.nplns >= set->so.nplns_allocated) { /*realloc*/
            set->so.nplns_allocated += 10;
            set->so.plns = (SvOutputPar *)g_realloc(set->so.plns, set->so.nplns_allocated * sizeof(SvOutputPar));
        }
    }

    if (SV_OUTTYPE_VOLUME == type) {
        if (set->so.ncubs_allocated == 0) {
            set->so.ncubs_allocated = 10;
            set->so.cubs = (SvOutputPar *)g_malloc(set->so.ncubs_allocated * sizeof(SvOutputPar));
        } else if (set->so.ncubs >= set->so.ncubs_allocated) { /*realloc*/
            set->so.ncubs_allocated += 10;
            set->so.cubs = (SvOutputPar *)g_realloc(set->so.cubs, set->so.ncubs_allocated * sizeof(SvOutputPar));
        }
    }

    if (SV_OUTTYPE_SUM == type || SV_OUTTYPE_SUMTAB == type) {
        if (set->so.nsums_allocated == 0) {
            set->so.nsums_allocated = 10;
            set->so.sums = (SvOutputSum *)g_malloc(set->so.nsums_allocated * sizeof(SvOutputSum));
        } else if (set->so.nsums >= set->so.nsums_allocated) { /*realloc*/
            set->so.nsums_allocated += 10;
            set->so.sums = (SvOutputSum *)g_realloc(set->so.sums, set->so.nsums_allocated * sizeof(SvOutputSum));
        }
    }

    if (SV_OUTTYPE_FORCE == type) {
        if (set->so.nforces_allocated == 0) {
            set->so.nforces_allocated = 10;
            set->so.forces = (SvOutputForce *)g_malloc(set->so.nforces_allocated * sizeof(SvOutputForce));
        } else if (set->so.nforces >= set->so.nforces_allocated) { /*realloc*/
            set->so.nforces_allocated += 10;
            set->so.forces = (SvOutputForce *)g_realloc(set->so.forces, set->so.nforces_allocated * sizeof(SvOutputForce));
        }
    }
}

void
init_settings(SvSet *set, SvType type)
{
    /*inits*/

    gint i;

    if (SV_SRCTYPE_POINT == type) {
        //g_assert(set->ss.pnts == NULL);        

        source_point_alloc(set);

        set->ss.npnts++;
        //data->is_psrc[data->set.ss.npnts-1] = TRUE;     /* set visibility */
    } else if (SV_SRCTYPE_SF == type) {
        /* 1. Source properties */
        set->ss.sf.source_mode = SOURCE_MODE;
        g_free(set->ss.sf.source_filename);
        set->ss.sf.source_filename = NULL;
        set->ss.sf.source_wl = SOURCE_WAVELENGTH;
        set->ss.sf.source_wlspan = SOURCE_WAVELENGTHSPAN;
        set->ss.sf.source_amplitude = SOURCE_AMPLITUDE;
        set->ss.sf.source_pulsewidth = SOURCE_PULSE_WIDTH;

        /* 2. Incident angle */
        set->ss.sf.ia_theta = SOURCE_IA_THETADEG;
        set->ss.sf.ia_phi = SOURCE_IA_PHIDEG;
        set->ss.sf.ia_psi = SOURCE_IA_PSIDEG;

        set->ss.sf.layered_epsilon = SOURCE_LAYERED_EPS;
        set->ss.sf.layered_mu = SOURCE_LAYERED_MU;
        set->ss.sf.layered_sigma = SOURCE_LAYERED_SIGMA;
        set->ss.sf.layered_sigast = SOURCE_LAYERED_SIGAST;

        set->ss.sf.valid = TRUE;        
    } else if (SV_SRCTYPE_TSF == type) {
        /* 0a. Box properties */
        set->ss.tsf.box_i0 = set->ss.tsf.box_j0 = set->ss.tsf.box_k0 = SOURCE_BOX_VERTEX0;
        set->ss.tsf.box_in = set->ss.tsf.box_jn = set->ss.tsf.box_kn = SOURCE_BOX_VERTEXN;
        set->ss.tsf.box_boundary_skipi0 = set->ss.tsf.box_boundary_skipj0 = set->ss.tsf.box_boundary_skipk0 = SOURCE_BOUNDARY_SKIP;
        set->ss.tsf.box_boundary_skipin = set->ss.tsf.box_boundary_skipjn = set->ss.tsf.box_boundary_skipkn = SOURCE_BOUNDARY_SKIP;
        set->ss.tsf.box_boundary_skipdepth = SOURCE_BOUNDARY_SKIP_DEPTH;

        /* 1. Source properties */
        set->ss.tsf.source_mode = SOURCE_MODE;
        g_free(set->ss.tsf.source_filename);
        set->ss.tsf.source_filename = NULL;
        set->ss.tsf.source_wl = SOURCE_WAVELENGTH;
        set->ss.tsf.source_wlspan = SOURCE_WAVELENGTHSPAN;
        set->ss.tsf.source_amplitude = SOURCE_AMPLITUDE;
        set->ss.tsf.source_pulsewidth = SOURCE_PULSE_WIDTH;

        /* 2. Incident angle */
        set->ss.tsf.ia_theta = SOURCE_IA_THETADEG;
        set->ss.tsf.ia_phi = SOURCE_IA_PHIDEG;
        set->ss.tsf.ia_psi = SOURCE_IA_PSIDEG;
        
        set->ss.tsf.layered_epsilon = SOURCE_LAYERED_EPS;
        set->ss.tsf.layered_mu = SOURCE_LAYERED_MU;
        set->ss.tsf.layered_sigma = SOURCE_LAYERED_SIGMA;
        set->ss.tsf.layered_sigast = SOURCE_LAYERED_SIGAST;

        set->ss.tsf.gaussian = SOURCE_GAUSSIAN_MULT_ENABLE;
        set->ss.tsf.gaussian_fxpos = SOURCE_XXX_MULT_CENTER;
        set->ss.tsf.gaussian_fypos = SOURCE_XXX_MULT_CENTER;
        set->ss.tsf.gaussian_rx = SOURCE_XXX_MULT_RADIUS;
        set->ss.tsf.gaussian_ry = SOURCE_XXX_MULT_RADIUS;

        set->ss.tsf.radial = SOURCE_RADIAL_MULT_ENABLE;
        set->ss.tsf.radial_fxpos = SOURCE_XXX_MULT_CENTER;
        set->ss.tsf.radial_fypos = SOURCE_XXX_MULT_CENTER;
        set->ss.tsf.radial_rx = SOURCE_XXX_MULT_RADIUS;
        set->ss.tsf.radial_ry = SOURCE_XXX_MULT_RADIUS;

        set->ss.tsf.fiber = SOURCE_FIBER_MULT_ENABLE;
        set->ss.tsf.fiber_radius = SOURCE_XXX_MULT_RADIUS;
        set->ss.tsf.fiber_fxpos = SOURCE_XXX_MULT_CENTER;
        set->ss.tsf.fiber_fypos = SOURCE_XXX_MULT_CENTER;
        set->ss.tsf.fiber_cutoff = SOURCE_FIBER_MULT_CUTOFF;
        set->ss.tsf.fiber_epsilon_core = SOURCE_FIBER_MULT_EPSILON_CORE;
        set->ss.tsf.fiber_epsilon_cladding = SOURCE_FIBER_MULT_EPSILON_CLADDING;

        set->ss.tsf.valid = TRUE;
    } else if (SV_SRCTYPE_LTSF == type) {
        /* 0a. Box properties */
        set->ss.ltsf.box_i0 = set->ss.ltsf.box_j0 = set->ss.ltsf.box_k0 = SOURCE_BOX_VERTEX0;
        set->ss.ltsf.box_in = set->ss.ltsf.box_jn = set->ss.ltsf.box_kn = SOURCE_BOX_VERTEXN;
        set->ss.ltsf.box_boundary_skipi0 = set->ss.ltsf.box_boundary_skipj0 = set->ss.ltsf.box_boundary_skipk0 = SOURCE_BOUNDARY_SKIP;
        set->ss.ltsf.box_boundary_skipin = set->ss.ltsf.box_boundary_skipjn = set->ss.ltsf.box_boundary_skipkn = SOURCE_BOUNDARY_SKIP;
        set->ss.ltsf.box_boundary_skipdepth = SOURCE_BOUNDARY_SKIP_DEPTH;

        /* 1. Source properties */
        set->ss.ltsf.source_mode = SOURCE_MODE;
        g_free(set->ss.ltsf.source_filename);
        set->ss.ltsf.source_filename = NULL;
        set->ss.ltsf.source_wl = SOURCE_WAVELENGTH;
        set->ss.ltsf.source_wlspan = SOURCE_WAVELENGTHSPAN;
        set->ss.ltsf.source_amplitude = SOURCE_AMPLITUDE;
        set->ss.ltsf.source_pulsewidth = SOURCE_PULSE_WIDTH;

        /* 2. Incident angle */
        set->ss.ltsf.ia_theta = SOURCE_IA_THETADEG;
        set->ss.ltsf.ia_phi = SOURCE_IA_PHIDEG;
        set->ss.ltsf.ia_psi = SOURCE_IA_PSIDEG;

        /* 4. Layered source properties */
        for (i = 0; i < set->ss.ltsf.layered_count; i++)
            g_free(set->ss.ltsf.layered_material[i]);

        set->ss.ltsf.layered_zpos[0] = SOURCE_LAYERED_ZPOS;
        set->ss.ltsf.layered_epsilon[0] = SOURCE_LAYERED_EPS;
        set->ss.ltsf.layered_mu[0] = SOURCE_LAYERED_MU;
        set->ss.ltsf.layered_sigma[0] = SOURCE_LAYERED_SIGMA;
        set->ss.ltsf.layered_sigast[0] = SOURCE_LAYERED_SIGAST;
        set->ss.ltsf.layered_material[0] = g_strdup(SOURCE_LAYERED_MATERIAL);

        set->ss.ltsf.layered_zpos[1] = SOURCE_LAYERED_ZPOS;
        set->ss.ltsf.layered_epsilon[1] = SOURCE_LAYERED_EPS;
        set->ss.ltsf.layered_mu[1] = SOURCE_LAYERED_MU;
        set->ss.ltsf.layered_sigma[1] = SOURCE_LAYERED_SIGMA;
        set->ss.ltsf.layered_sigast[1] = SOURCE_LAYERED_SIGAST;
        set->ss.ltsf.layered_material[1] = g_strdup(SOURCE_LAYERED_MATERIAL);

        set->ss.ltsf.layered_count = 2;

        set->ss.ltsf.gaussian = SOURCE_GAUSSIAN_MULT_ENABLE;
        set->ss.ltsf.gaussian_fxpos = SOURCE_XXX_MULT_CENTER;
        set->ss.ltsf.gaussian_fypos = SOURCE_XXX_MULT_CENTER;
        set->ss.ltsf.gaussian_rx = SOURCE_XXX_MULT_RADIUS;
        set->ss.ltsf.gaussian_ry = SOURCE_XXX_MULT_RADIUS;

        set->ss.ltsf.radial = SOURCE_RADIAL_MULT_ENABLE;
        set->ss.ltsf.radial_fxpos = SOURCE_XXX_MULT_CENTER;
        set->ss.ltsf.radial_fypos = SOURCE_XXX_MULT_CENTER;
        set->ss.ltsf.radial_rx = SOURCE_XXX_MULT_RADIUS;
        set->ss.ltsf.radial_ry = SOURCE_XXX_MULT_RADIUS;

        set->ss.ltsf.fiber = SOURCE_FIBER_MULT_ENABLE;
        set->ss.ltsf.fiber_fxpos = SOURCE_XXX_MULT_CENTER;
        set->ss.ltsf.fiber_fypos = SOURCE_XXX_MULT_CENTER;
        set->ss.ltsf.fiber_radius = SOURCE_XXX_MULT_RADIUS;
        set->ss.ltsf.fiber_cutoff = SOURCE_FIBER_MULT_CUTOFF;
        set->ss.ltsf.fiber_epsilon_core = SOURCE_FIBER_MULT_EPSILON_CORE;
        set->ss.ltsf.fiber_epsilon_cladding = SOURCE_FIBER_MULT_EPSILON_CLADDING;

        set->ss.ltsf.valid = TRUE;
    } else if (SV_SRCTYPE_TSFF == type) {
        /* 0a. Box properties */
        set->ss.tsff.box_i0 = set->ss.tsff.box_j0 = set->ss.tsff.box_k0 = SOURCE_BOX_VERTEX0;
        set->ss.tsff.box_in = set->ss.tsff.box_jn = set->ss.tsff.box_kn = SOURCE_BOX_VERTEXN;
        set->ss.tsff.box_boundary_skipi0 = set->ss.tsff.box_boundary_skipj0 = set->ss.tsff.box_boundary_skipk0 = SOURCE_BOUNDARY_SKIP;
        set->ss.tsff.box_boundary_skipin = set->ss.tsff.box_boundary_skipjn = set->ss.tsff.box_boundary_skipkn = SOURCE_BOUNDARY_SKIP;
        set->ss.tsff.box_boundary_skipdepth = SOURCE_BOUNDARY_SKIP_DEPTH;

        /* 1. Source properties */
        set->ss.tsff.source_mode = 1;
        g_free(set->ss.tsff.source_filename);
        set->ss.tsff.source_mode = SOURCE_MODE;
        set->ss.tsff.source_wl = SOURCE_WAVELENGTH;
        set->ss.tsff.source_wlspan = SOURCE_WAVELENGTHSPAN;
        set->ss.tsff.source_amplitude = SOURCE_AMPLITUDE;
        set->ss.tsff.source_pulsewidth = SOURCE_PULSE_WIDTH;

        /* 3. Focused source properties */
        set->ss.tsff.focused_thetamax = SOURCE_FS_THETAMAXDEG;  /* aperture max angle [rad] */
        set->ss.tsff.focused_fdist = SOURCE_FS_FDIST;           /* focal distance */
        set->ss.tsff.focused_pol = SOURCE_FS_POLARISATIONDEG;   /* polarization angle [rad]*/
        set->ss.tsff.focused_nip = SOURCE_FS_NIP;               /* N integration points */
        set->ss.tsff.focused_mip = SOURCE_FS_MIP;               /* M integration points */

        set->ss.tsff.xshift = SOURCE_TSFF_SHIFT;
        set->ss.tsff.yshift = SOURCE_TSFF_SHIFT;
        set->ss.tsff.zshift = SOURCE_TSFF_SHIFT;

        set->ss.tsff.valid = TRUE;
    } else if (SV_SRCTYPE_LTSFF == type) {
        /* 0a. Box properties */
        set->ss.ltsff.box_i0 = set->ss.ltsff.box_j0 = set->ss.ltsff.box_k0 = SOURCE_BOX_VERTEX0;
        set->ss.ltsff.box_in = set->ss.ltsff.box_jn = set->ss.ltsff.box_kn = SOURCE_BOX_VERTEXN;
        set->ss.ltsff.box_boundary_skipi0 = set->ss.ltsff.box_boundary_skipj0 = set->ss.ltsff.box_boundary_skipk0 = SOURCE_BOUNDARY_SKIP;
        set->ss.ltsff.box_boundary_skipin = set->ss.ltsff.box_boundary_skipjn = set->ss.ltsff.box_boundary_skipkn = SOURCE_BOUNDARY_SKIP;
        set->ss.ltsff.box_boundary_skipdepth = SOURCE_BOUNDARY_SKIP_DEPTH;

        /* 1. Source properties */
        g_free(set->ss.ltsff.source_filename);
        set->ss.ltsff.source_mode = SOURCE_MODE;
        set->ss.ltsff.source_wl = SOURCE_WAVELENGTH;
        set->ss.ltsff.source_wlspan = SOURCE_WAVELENGTHSPAN;
        set->ss.ltsff.source_amplitude = SOURCE_AMPLITUDE;
        set->ss.ltsff.source_pulsewidth = SOURCE_PULSE_WIDTH;

        /* 3. Focused source properties */
        set->ss.ltsff.focused_thetamax = SOURCE_FS_THETAMAXDEG; /* aperture max angle [rad] */
        set->ss.ltsff.focused_fdist = SOURCE_FS_FDIST;          /* focal distance */
        set->ss.ltsff.focused_pol = SOURCE_FS_POLARISATIONDEG;  /* polarization angle [rad]*/
        set->ss.ltsff.focused_nip = SOURCE_FS_NIP;              /* N integration points */
        set->ss.ltsff.focused_mip = SOURCE_FS_MIP;              /* M integration points */

        /* 4. Layered source properties */
        for (i = 0; i < set->ss.ltsff.layered_count; i++)
            g_free(set->ss.ltsff.layered_material[i]);

        set->ss.ltsff.layered_zpos[0] = SOURCE_LAYERED_ZPOS;
        set->ss.ltsff.layered_epsilon[0] = SOURCE_LAYERED_EPS;
        set->ss.ltsff.layered_mu[0] = SOURCE_LAYERED_MU;
        set->ss.ltsff.layered_sigma[0] = SOURCE_LAYERED_SIGMA;
        set->ss.ltsff.layered_sigast[0] = SOURCE_LAYERED_SIGAST;
        set->ss.ltsff.layered_material[0] = g_strdup(SOURCE_LAYERED_MATERIAL);

        set->ss.ltsff.layered_zpos[1] = SOURCE_LAYERED_ZPOS;
        set->ss.ltsff.layered_epsilon[1] = SOURCE_LAYERED_EPS;
        set->ss.ltsff.layered_mu[1] = SOURCE_LAYERED_MU;
        set->ss.ltsff.layered_sigma[1] = SOURCE_LAYERED_SIGMA;
        set->ss.ltsff.layered_sigast[1] = SOURCE_LAYERED_SIGAST;
        set->ss.ltsff.layered_material[1] = g_strdup(SOURCE_LAYERED_MATERIAL);

        set->ss.ltsff.layered_count = 2;

        set->ss.ltsff.fast = SOURCE_LTSFF_FAST;

        set->ss.ltsff.valid = TRUE;
    } else if (SV_OUTTYPE_POINT == type) {
        set->so.npnts++;
        output_alloc(set, SV_OUTTYPE_POINT);

        set->so.pnts[set->so.npnts - 1].i = OUTPUT_POINT_POS;
        set->so.pnts[set->so.npnts - 1].j = OUTPUT_POINT_POS;
        set->so.pnts[set->so.npnts - 1].k = OUTPUT_POINT_POS;
        set->so.pnts[set->so.npnts - 1].step = OUTPUT_STEP;
        set->so.pnts[set->so.npnts - 1].component = OUTPUT_COMPONENT;
        memset(set->so.pnts[set->so.npnts - 1].filebase, 0, sizeof(set->so.pnts[set->so.npnts - 1].filebase));
        g_snprintf(set->so.pnts[set->so.npnts - 1].filebase, sizeof(set->so.pnts[set->so.npnts - 1].filebase), GENERAL_FILENAME);
        set->so.pnts[set->so.npnts - 1].start = OUTPUT_START;
        set->so.pnts[set->so.npnts - 1].stop = OUTPUT_STOP;
        set->so.pnts[set->so.npnts - 1].format = OUTPUT_FORMAT;
    } else if (SV_OUTTYPE_IMAGE == type) {
        set->so.nimgs++;
        output_alloc(set, SV_OUTTYPE_IMAGE);

        set->so.imgs[set->so.nimgs - 1].i = OUTPUT_IMAGE_POSI;
        set->so.imgs[set->so.nimgs - 1].j = OUTPUT_IMAGE_POSJ;
        set->so.imgs[set->so.nimgs - 1].k = OUTPUT_IMAGE_POSK;
        set->so.imgs[set->so.nimgs - 1].step = OUTPUT_STEP;
        set->so.imgs[set->so.nimgs - 1].component = OUTPUT_COMPONENT;
        memset(set->so.imgs[set->so.nimgs - 1].filebase, 0, sizeof(set->so.imgs[set->so.nimgs - 1].filebase));
        g_snprintf(set->so.imgs[set->so.nimgs - 1].filebase, sizeof(set->so.imgs[set->so.nimgs - 1].filebase), GENERAL_FILENAME);
        set->so.imgs[set->so.nimgs - 1].start = OUTPUT_START;
        set->so.imgs[set->so.nimgs - 1].stop = OUTPUT_STOP;
        set->so.imgs[set->so.nimgs - 1].format = OUTPUT_FORMAT;
    } else if (SV_OUTTYPE_SUBGRID_IMAGE == type) {
        set->so.nsubgrid_imgs++;
        output_alloc(set, SV_OUTTYPE_SUBGRID_IMAGE);

        set->so.subgrid_imgs[set->so.nsubgrid_imgs - 1].i = OUTPUT_IMAGE_POSI;
        set->so.subgrid_imgs[set->so.nsubgrid_imgs - 1].j = OUTPUT_IMAGE_POSJ;
        set->so.subgrid_imgs[set->so.nsubgrid_imgs - 1].k = OUTPUT_IMAGE_POSK;
        set->so.subgrid_imgs[set->so.nsubgrid_imgs - 1].step = OUTPUT_STEP;
        set->so.subgrid_imgs[set->so.nsubgrid_imgs - 1].component = OUTPUT_COMPONENT;
        memset(set->so.subgrid_imgs[set->so.nsubgrid_imgs - 1].filebase, 0, sizeof(set->so.subgrid_imgs[set->so.nsubgrid_imgs - 1].filebase));
        g_snprintf(set->so.subgrid_imgs[set->so.nsubgrid_imgs - 1].filebase, sizeof(set->so.subgrid_imgs[set->so.nsubgrid_imgs - 1].filebase), GENERAL_FILENAME);
        set->so.subgrid_imgs[set->so.nsubgrid_imgs - 1].start = OUTPUT_START;
        set->so.subgrid_imgs[set->so.nsubgrid_imgs - 1].stop = OUTPUT_STOP;
        set->so.subgrid_imgs[set->so.nsubgrid_imgs - 1].format = OUTPUT_FORMAT;
    } else if (SV_OUTTYPE_PLANE == type) {
        set->so.nplns++;
        output_alloc(set, SV_OUTTYPE_PLANE);

        set->so.plns[set->so.nplns - 1].i = OUTPUT_PLANE_POSI;
        set->so.plns[set->so.nplns - 1].j = OUTPUT_PLANE_POSJ;
        set->so.plns[set->so.nplns - 1].k = OUTPUT_PLANE_POSK;
        set->so.plns[set->so.nplns - 1].step = OUTPUT_STEP;
        set->so.plns[set->so.nplns - 1].component = OUTPUT_COMPONENT;
        memset(set->so.plns[set->so.nplns - 1].filebase, 0, sizeof(set->so.plns[set->so.nplns - 1].filebase));
        g_snprintf(set->so.plns[set->so.nplns - 1].filebase, sizeof(set->so.plns[set->so.nplns - 1].filebase), GENERAL_FILENAME);
        set->so.plns[set->so.nplns - 1].start = OUTPUT_START;
        set->so.plns[set->so.nplns - 1].stop = OUTPUT_STOP;
        set->so.plns[set->so.nplns - 1].format = OUTPUT_FORMAT;
    } else if (SV_OUTTYPE_VOLUME == type) {
        set->so.ncubs++;
        output_alloc(set, SV_OUTTYPE_VOLUME);

        set->so.cubs[set->so.ncubs - 1].i = OUTPUT_VOLUME_POSI;
        set->so.cubs[set->so.ncubs - 1].j = OUTPUT_VOLUME_POSJ;
        set->so.cubs[set->so.ncubs - 1].k = OUTPUT_VOLUME_POSK;
        set->so.cubs[set->so.ncubs - 1].step = OUTPUT_STEP;
        set->so.cubs[set->so.ncubs - 1].component = OUTPUT_COMPONENT;
        memset(set->so.cubs[set->so.ncubs - 1].filebase, 0, sizeof(set->so.cubs[set->so.ncubs - 1].filebase));
        g_snprintf(set->so.cubs[set->so.ncubs - 1].filebase, sizeof(set->so.cubs[set->so.ncubs - 1].filebase), GENERAL_FILENAME);
        set->so.cubs[set->so.ncubs - 1].start = OUTPUT_START;
        set->so.cubs[set->so.ncubs - 1].stop = OUTPUT_STOP;
        set->so.cubs[set->so.ncubs - 1].format = OUTPUT_FORMAT;
    } else if (SV_OUTTYPE_SUM == type || SV_OUTTYPE_SUMTAB == type) {
        set->so.nsums++;
        output_alloc(set, SV_OUTTYPE_SUM);

        set->so.sums[set->so.nsums - 1].box_i0 = OUTPUT_SUM_VERTEX0;
        set->so.sums[set->so.nsums - 1].box_j0 = OUTPUT_SUM_VERTEX0;
        set->so.sums[set->so.nsums - 1].box_k0 = OUTPUT_SUM_VERTEX0;
        set->so.sums[set->so.nsums - 1].box_in = OUTPUT_SUM_VERTEXN;
        set->so.sums[set->so.nsums - 1].box_jn = OUTPUT_SUM_VERTEXN;
        set->so.sums[set->so.nsums - 1].box_kn = OUTPUT_SUM_VERTEXN;

        set->so.sums[set->so.nsums - 1].component = OUTPUT_COMPONENT;
        set->so.sums[set->so.nsums - 1].step = OUTPUT_STEP;

        set->so.sums[set->so.nsums - 1].layered_epsilon = OUTPUT_SUM_EPS;
        set->so.sums[set->so.nsums - 1].layered_mu = OUTPUT_SUM_MU;
        set->so.sums[set->so.nsums - 1].layered_sigma = OUTPUT_SUM_SIGMA;
        set->so.sums[set->so.nsums - 1].layered_sigast = OUTPUT_SUM_SIGAST;

        set->so.sums[set->so.nsums - 1].stringbased = 0;
        memset(set->so.sums[set->so.nsums - 1].string, 0, sizeof(set->so.sums[set->so.nsums - 1].string));

        memset(set->so.sums[set->so.nsums - 1].filename, 0, sizeof(set->so.sums[set->so.nsums - 1].filename));
        g_snprintf(set->so.sums[set->so.nsums - 1].filename, sizeof(set->so.sums[set->so.nsums - 1].filename), GENERAL_FILENAME);
    } else if (SV_OUTTYPE_FORCE == type) {
        set->so.nforces++;
        output_alloc(set, SV_OUTTYPE_FORCE);

        set->so.forces[set->so.nforces - 1].box_i0 = OUTPUT_FORCE_VERTEX0;
        set->so.forces[set->so.nforces - 1].box_j0 = OUTPUT_FORCE_VERTEX0;
        set->so.forces[set->so.nforces - 1].box_k0 = OUTPUT_FORCE_VERTEX0;
        set->so.forces[set->so.nforces - 1].box_in = OUTPUT_FORCE_VERTEXN;
        set->so.forces[set->so.nforces - 1].box_jn = OUTPUT_FORCE_VERTEXN;
        set->so.forces[set->so.nforces - 1].box_kn = OUTPUT_FORCE_VERTEXN;

        set->so.forces[set->so.nforces - 1].step = OUTPUT_STEP;

        memset(set->so.forces[set->so.nforces - 1].filename, 0, sizeof(set->so.forces[set->so.nforces - 1].filename));
        g_snprintf(set->so.forces[set->so.nforces - 1].filename, sizeof(set->so.forces[set->so.nforces - 1].filename), GENERAL_FILENAME);
    } else if (SV_TYPE_MATERIAL_GROW == type) {
        set->sm.ngrowths++;
        set->sm.grow_i0[set->sm.ngrowths - 1] = MATERIAL_GROW_I0;
        set->sm.grow_j0[set->sm.ngrowths - 1] = MATERIAL_GROW_J0;
        set->sm.grow_k0[set->sm.ngrowths - 1] = MATERIAL_GROW_K0;
        set->sm.grow_in[set->sm.ngrowths - 1] = set->sp.xres - MATERIAL_GROW_INBORDER;
        set->sm.grow_jn[set->sm.ngrowths - 1] = set->sp.yres - MATERIAL_GROW_JNBORDER;
        set->sm.grow_kn[set->sm.ngrowths - 1] = set->sp.zres - MATERIAL_GROW_KNBORDER;

        set->sm.grow_skipi0[set->sm.ngrowths - 1] = MATERIAL_GROW_SKIP;
        set->sm.grow_skipj0[set->sm.ngrowths - 1] = MATERIAL_GROW_SKIP;
        set->sm.grow_skipk0[set->sm.ngrowths - 1] = MATERIAL_GROW_SKIP;
        set->sm.grow_skipin[set->sm.ngrowths - 1] = MATERIAL_GROW_SKIP;
        set->sm.grow_skipjn[set->sm.ngrowths - 1] = MATERIAL_GROW_SKIP;
        set->sm.grow_skipkn[set->sm.ngrowths - 1] = MATERIAL_GROW_SKIP;

        set->sm.grow_addindex[set->sm.ngrowths - 1] = MATERIAL_GROW_ADDINDEX;
        set->sm.grow_attachindex[set->sm.ngrowths - 1] = MATERIAL_GROW_ATTACHINDEX;
        set->sm.grow_subsampling[set->sm.ngrowths - 1] = MATERIAL_GROW_SUBSAMPLING;
        set->sm.grow_seed[set->sm.ngrowths - 1] = MATERIAL_GROW_SEED;
        set->sm.grow_nsteps[set->sm.ngrowths - 1] = MATERIAL_GROW_NSTEPS;
        set->sm.grow_mobility[set->sm.ngrowths - 1] = MATERIAL_GROW_MOBILITY;
        set->sm.grow_probability[set->sm.ngrowths - 1] = MATERIAL_GROW_PROBABILITY;

        set->sm.localfiles = TRUE;
    } else if (SV_TYPE_MATERIAL_ROUGHNESS == type) {
        set->sm.nroughens++;
        set->sm.rough_radius_peak[set->sm.nroughens - 1] = MATERIAL_ROUGHNESS_RADIUSPEAK;
        set->sm.rough_radius_span[set->sm.nroughens - 1] = MATERIAL_ROUGHNESS_RADIUSSPAN;
        set->sm.rough_iterations[set->sm.nroughens - 1] = MATERIAL_ROUGHNESS_ITERATIONS;
        set->sm.rough_probability[set->sm.nroughens - 1] = MATERIAL_ROUGHNESS_PROBABILITY;
        set->sm.rough_matindex[set->sm.nroughens - 1] = MATERIAL_ROUGHNESS_MATINDEX;
        set->sm.rough_voidindex[set->sm.nroughens - 1] = MATERIAL_ROUGHNESS_VOIDINDEX;
        set->sm.rough_seed[set->sm.nroughens - 1] = MATERIAL_ROUGHNESS_SEED;
    } else if (SV_TYPE_MATERIAL_SPECTRAL == type) {
        set->sm.nspectrals++;
        set->sm.spectral_sigma[set->sm.nspectrals - 1] = MATERIAL_SPECTRAL_SIGMA;
        set->sm.spectral_t[set->sm.nspectrals - 1] = MATERIAL_SPECTRAL_T;
        set->sm.spectral_matindex[set->sm.nspectrals - 1] = MATERIAL_SPECTRAL_MATINDEX;
        set->sm.spectral_seed[set->sm.nspectrals - 1] = MATERIAL_SPECTRAL_SEED;
    } else if (SV_TYPE_MATERIAL_EXPRESSION == type) {
        set->sm.nexprs++;
        set->sm.expr_i0[set->sm.nexprs - 1] = MATERIAL_EXPR_I0;
        set->sm.expr_j0[set->sm.nexprs - 1] = MATERIAL_EXPR_J0;
        set->sm.expr_k0[set->sm.nexprs - 1] = MATERIAL_EXPR_K0;
        set->sm.expr_in[set->sm.nexprs - 1] = set->sp.xres - MATERIAL_EXPR_INBORDER;
        set->sm.expr_jn[set->sm.nexprs - 1] = set->sp.yres - MATERIAL_EXPR_JNBORDER;
        set->sm.expr_kn[set->sm.nexprs - 1] = set->sp.zres - MATERIAL_EXPR_KNBORDER;
        set->sm.expr_matindex[set->sm.nexprs - 1] = MATERIAL_EXPR_MATINDEX;
        set->sm.expr_voidindex[set->sm.nexprs - 1] = MATERIAL_EXPR_VOIDINDEX;
        set->sm.expr_maxdist[set->sm.nexprs - 1] = MATERIAL_EXPR_MAXDIST;
        set->sm.expr_distmode[set->sm.nexprs - 1] = MATERIAL_EXPR_DISTMODE;
        memset(set->sm.expr_expr[set->sm.nexprs - 1], 0, sizeof(gchar)*MATERIAL_EXPR_EXPR_CHARS);
        g_snprintf(set->sm.expr_expr[set->sm.nexprs - 1], sizeof(gchar)*MATERIAL_EXPR_EXPR_CHARS, MATERIAL_EXPR_EXPR);
    } else if (SV_TYPE_NFFF_POINT == type) {
        set->sf.box_i0 = set->sf.box_j0 = set->sf.box_k0 = NFFF_BOX_VERTEX0;
        set->sf.box_in = set->sf.box_jn = set->sf.box_kn = NFFF_BOX_VERTEXN;
        set->sf.box_boundary_skipi0 = NFFF_BOX_BOUNDARY_SKIP;
        set->sf.skipi0_jmin = set->sf.skipi0_kmin = set->sf.skipi0_jmax = set->sf.skipi0_kmax = NFFF_SKIPI0;
        set->sf.box_boundary_skipin = NFFF_BOX_BOUNDARY_SKIP;
        set->sf.skipin_jmin = set->sf.skipin_kmin = set->sf.skipin_jmax = set->sf.skipin_kmax = NFFF_SKIPIN;
        set->sf.box_boundary_skipj0 = NFFF_BOX_BOUNDARY_SKIP;
        set->sf.skipj0_imin = set->sf.skipj0_kmin = set->sf.skipj0_imax = set->sf.skipj0_kmax = NFFF_SKIPJ0;
        set->sf.box_boundary_skipjn = NFFF_BOX_BOUNDARY_SKIP;
        set->sf.skipjn_imin = set->sf.skipjn_kmin = set->sf.skipjn_imax = set->sf.skipjn_kmax = NFFF_SKIPJN;
        set->sf.box_boundary_skipk0 = NFFF_BOX_BOUNDARY_SKIP;
        set->sf.skipk0_imin = set->sf.skipk0_jmin = set->sf.skipk0_imax = set->sf.skipk0_jmax = NFFF_SKIPK0;
        set->sf.box_boundary_skipkn = NFFF_BOX_BOUNDARY_SKIP;
        set->sf.skipkn_imin = set->sf.skipkn_jmin = set->sf.skipkn_imax = set->sf.skipkn_jmax = NFFF_SKIPKN;

        set->sf.nrs++;
        if (set->sf.nrs == 1) {
            set->sf.source_filename = (gchar **)g_malloc(set->sf.nrs * sizeof(gchar *));
            set->sf.ri = (gint *)g_malloc(set->sf.nrs * sizeof(gint));
            set->sf.rj = (gint *)g_malloc(set->sf.nrs * sizeof(gint));
            set->sf.rk = (gint *)g_malloc(set->sf.nrs * sizeof(gint));
            set->sf.individual = (gint *)g_malloc(set->sf.nrs * sizeof(gint));
        } else {
            set->sf.ri = (gint *)g_realloc(set->sf.ri, set->sf.nrs * sizeof(gint));
            set->sf.rj = (gint *)g_realloc(set->sf.rj, set->sf.nrs * sizeof(gint));
            set->sf.rk = (gint *)g_realloc(set->sf.rk, set->sf.nrs * sizeof(gint));
            set->sf.individual = (gint *)g_realloc(set->sf.individual, set->sf.nrs * sizeof(gint));
            set->sf.source_filename = (gchar **)g_realloc(set->sf.source_filename, set->sf.nrs * sizeof(gchar *));
        }

        set->sf.source_filename[set->sf.nrs - 1] = g_strdup(GENERAL_FILENAME);
        set->sf.ri[set->sf.nrs - 1] = NFFFP_POS_I;
        set->sf.rj[set->sf.nrs - 1] = NFFFP_POS_J;
        set->sf.rk[set->sf.nrs - 1] = NFFFP_POS_K;
        set->sf.individual[set->sf.nrs - 1] = NFFFP_INDIVIDUAL;

        //set->sf.valid = TRUE;
    } else if (SV_TYPE_NFFF_AREA == type) {
        set->sf.nareas++;
        set->sf.area_thetares[set->sf.nareas - 1] = NFFFA_THETARES;
        set->sf.area_phires[set->sf.nareas - 1] = NFFFA_PHIRES;
        set->sf.area_radius[set->sf.nareas - 1] = NFFFA_RADIUS;
        set->sf.area_thetafrom[set->sf.nareas - 1] = NFFFA_THETAFROMDEG  * G_PI / 180;
        set->sf.area_phifrom[set->sf.nareas - 1] = NFFFA_PHIFROMDEG * G_PI / 180;
        set->sf.area_thetato[set->sf.nareas - 1] = NFFFA_THETATODEG  * G_PI / 180;
        set->sf.area_phito[set->sf.nareas - 1] = NFFFA_PHITODEG * G_PI / 180;
        set->sf.area_savefile[set->sf.nareas - 1] = NFFFA_SAVEFILE;

        //set->sf.valid = TRUE;
    } else if (SV_TYPE_PNFFF_POINT == type) {
        set->spf.box_i0 = set->spf.box_j0 = set->spf.box_k0 = PNFFF_BOX_VERTEX0;
        set->spf.box_in = set->spf.box_jn = set->spf.box_kn = PNFFF_BOX_VERTEXN;
        set->spf.box_boundary_skipk0 = PNFFF_BOX_BOUNDARY_SKIP;
        set->spf.skipk0_imin = set->spf.skipk0_jmin = set->spf.skipk0_imax = set->spf.skipk0_jmax = PNFFF_SKIPK0;
        set->spf.box_boundary_skipkn = PNFFF_BOX_BOUNDARY_SKIP;
        set->spf.skipkn_imin = set->spf.skipkn_jmin = set->spf.skipkn_imax = set->spf.skipkn_jmax = PNFFF_SKIPKN;

        set->spf.nrs++;
        if (set->spf.nrs == 1) {
            set->spf.source_filename = (gchar **)g_malloc(set->spf.nrs * sizeof(gchar *));
            set->spf.ri = (gint *)g_malloc(set->spf.nrs * sizeof(gint));
            set->spf.rj = (gint *)g_malloc(set->spf.nrs * sizeof(gint));
            set->spf.rk = (gint *)g_malloc(set->spf.nrs * sizeof(gint));
            set->spf.individual = (gint *)g_malloc(set->spf.nrs * sizeof(gint));
        } else {
            set->spf.ri = (gint *)g_realloc(set->spf.ri, set->spf.nrs * sizeof(gint));
            set->spf.rj = (gint *)g_realloc(set->spf.rj, set->spf.nrs * sizeof(gint));
            set->spf.rk = (gint *)g_realloc(set->spf.rk, set->spf.nrs * sizeof(gint));
            set->spf.individual = (gint *)g_realloc(set->spf.individual, set->spf.nrs * sizeof(gint));
            set->spf.source_filename = (gchar **)g_realloc(set->spf.source_filename, set->spf.nrs * sizeof(gchar *));
        }

        set->spf.source_filename[set->spf.nrs - 1] = g_strdup(GENERAL_FILENAME);
        set->spf.ri[set->spf.nrs - 1] = PNFFFP_POS_I;
        set->spf.rj[set->spf.nrs - 1] = PNFFFP_POS_J;
        set->spf.rk[set->spf.nrs - 1] = PNFFFP_POS_K;
        set->spf.individual[set->spf.nrs - 1] = PNFFFP_INDIVIDUAL;

        set->spf.pimin = PNFFF_INTEGRATION_XMIN;
        set->spf.pimax = PNFFF_INTEGRATION_XMAX;
        set->spf.pjmin = PNFFF_INTEGRATION_YMIN;
        set->spf.pjmax = PNFFF_INTEGRATION_YMAX;

        //set->spf.valid = TRUE;
    } else if (SV_TYPE_PNFFF_AREA == type) {
        set->spf.nareas++;
        set->spf.area_thetares[set->spf.nareas - 1] = PNFFFA_THETARES;
        set->spf.area_phires[set->spf.nareas - 1] = PNFFFA_PHIRES;
        set->spf.area_radius[set->spf.nareas - 1] = PNFFFA_RADIUS;
        set->spf.area_thetafrom[set->spf.nareas - 1] = PNFFFA_THETAFROMDEG  * G_PI / 180;
        set->spf.area_phifrom[set->spf.nareas - 1] = PNFFFA_PHIFROMDEG * G_PI / 180;
        set->spf.area_thetato[set->spf.nareas - 1] = PNFFFA_THETATODEG * G_PI / 180;
        set->spf.area_phito[set->spf.nareas - 1] = PNFFFA_PHITODEG * G_PI / 180;
        set->spf.area_savefile[set->spf.nareas - 1] = PNFFFA_SAVEFILE;

        //set->spf.valid = TRUE;
    }
} /* init_settings */

void
init_add_settings_mat(SvSetMat *set, SvTypeMat type)
{    
    /*inits*/

    SvMatProp mat;
    init_settings_mat_prop(&mat);

    if (SV_TYPE_MAT_SPHERE == type) {
        SvSphere sphere;
        sphere.pnt1[0] = MAT_SPHERE_CENTER_I;
        sphere.pnt1[1] = MAT_SPHERE_CENTER_J;
        sphere.pnt1[2] = MAT_SPHERE_CENTER_K;
        sphere.radius = MAT_SPHERE_RADIUS;
        sphere.n = set->ngtot++;
        memcpy(&sphere.mat, &mat, sizeof(SvMatProp));
        g_array_append_val(set->spheres, sphere);
    } else if (SV_TYPE_MAT_VOXEL == type) {
        SvVoxel voxel;
        voxel.pnt1[0] = MAT_VOXEL_I0;
        voxel.pnt1[1] = MAT_VOXEL_J0;
        voxel.pnt1[2] = MAT_VOXEL_K0;
        voxel.pnt2[0] = MAT_VOXEL_IN;
        voxel.pnt2[1] = MAT_VOXEL_JN;
        voxel.pnt2[2] = MAT_VOXEL_KN;
        voxel.n = set->ngtot++;
        memcpy(&voxel.mat, &mat, sizeof(SvMatProp));
        g_array_append_val(set->voxels, voxel);
    } else if (SV_TYPE_MAT_CYLINDER == type) {
        SvCylinder cyl;
        cyl.pnt1[0] = MAT_CYLINDER_I0;
        cyl.pnt1[1] = MAT_CYLINDER_J0;
        cyl.pnt1[2] = MAT_CYLINDER_K0;
        cyl.pnt2[0] = MAT_CYLINDER_IN;
        cyl.pnt2[1] = MAT_CYLINDER_JN;
        cyl.pnt2[2] = MAT_CYLINDER_KN;
        cyl.radius = MAT_CYLINDER_RADIUS;
        cyl.n = set->ngtot++;
        memcpy(&cyl.mat, &mat, sizeof(SvMatProp));
        g_array_append_val(set->cylinders, cyl);
    } else if (SV_TYPE_MAT_CONE == type) {
        SvCone cone;
        cone.pnt1[0] = MAT_CONE_I0;
        cone.pnt1[1] = MAT_CONE_J0;
        cone.pnt1[2] = MAT_CONE_K0;
        cone.pnt2[0] = MAT_CONE_IN;
        cone.pnt2[1] = MAT_CONE_JN;
        cone.pnt2[2] = MAT_CONE_KN;
        cone.radius = MAT_CONE_RADIUS;
        cone.n = set->ngtot++;
        memcpy(&cone.mat, &mat, sizeof(SvMatProp));
        g_array_append_val(set->cones, cone);
    } else if (SV_TYPE_MAT_RCONE == type) {
        SvRCone rcone;
        rcone.pnt1[0] = MAT_RCONE_I0;
        rcone.pnt1[1] = MAT_RCONE_J0;
        rcone.pnt1[2] = MAT_RCONE_K0;
        rcone.pnt2[0] = MAT_RCONE_IN;
        rcone.pnt2[1] = MAT_RCONE_JN;
        rcone.pnt2[2] = MAT_RCONE_KN;
        rcone.radius1 = MAT_RCONE_RADIUS1;
        rcone.radius2 = MAT_RCONE_RADIUS2;
        rcone.n = set->ngtot++;
        memcpy(&rcone.mat, &mat, sizeof(SvMatProp));
        g_array_append_val(set->rcones, rcone);
    } else if (SV_TYPE_MAT_GWYDD == type) {
    } else if (SV_TYPE_MAT_MESH == type) {
    } else if (SV_TYPE_MAT_TETRAHEDRON == type) {
        SvTetrahedron tthn;
        tthn.pnt1[0] = MAT_TTHN_I1;
        tthn.pnt1[1] = MAT_TTHN_J1;
        tthn.pnt1[2] = MAT_TTHN_K1;
        tthn.pnt2[0] = MAT_TTHN_I2;
        tthn.pnt2[1] = MAT_TTHN_J2;
        tthn.pnt2[2] = MAT_TTHN_K2;
        tthn.pnt3[0] = MAT_TTHN_I3;
        tthn.pnt3[1] = MAT_TTHN_J3;
        tthn.pnt3[2] = MAT_TTHN_K3;
        tthn.pnt4[0] = MAT_TTHN_I4;
        tthn.pnt4[1] = MAT_TTHN_J4;
        tthn.pnt4[2] = MAT_TTHN_K4;
        tthn.setpart = 0;
        tthn.n = set->ngtot++;
        memcpy(&tthn.mat, &mat, sizeof(SvMatProp));
        g_array_append_val(set->tetrahedrons, tthn);
    }
} /* init_settings_mat */

void init_settings_mat_prop(SvMatProp *mat)
{
    /*inits*/

    //SvMatProp mat;
    memset(mat, 0, sizeof(SvMatProp));
    mat->epsilon = 1;
    mat->sigma = 0;
    mat->mu = 1;
    mat->sigast = 0;
    mat->drude_omega_p = 0;
    mat->drude_nu = 0;
    mat->type = 0;
    //mat->material = g_strdup(MAT_PROP_MATERIAL_NONE);
} /* init_settings_mat */

gboolean
is_source_valid(SvSet *set, SvType type)
{
    if (SV_SRCTYPE_SF == type)
        return set->ss.sf.valid;
    if (SV_SRCTYPE_TSF == type)
        return set->ss.tsf.valid;
    else if (SV_SRCTYPE_TSFF == type)
        return set->ss.tsff.valid;
    else if (SV_SRCTYPE_LTSF == type)
        return set->ss.ltsf.valid;
    else if (SV_SRCTYPE_LTSFF == type)
        return set->ss.ltsff.valid;
    else if (SV_SRCTYPE_EXT == type)
        return set->ss.ext.valid;

    return FALSE;

    /*    if (SV_SRCTYPE_SF == type)
            return !set->ss.sf.valid;
        if(SV_SRCTYPE_TSF == type)
            return ((set->ss.tsf.i0 + set->ss.tsf.j0 + set->ss.tsf.k0 + set->ss.tsf.i1 + set->ss.tsf.j1 + set->ss.tsf.k1)==0);
        else if (SV_SRCTYPE_TSFF == type)
            return ((set->ss.tsff.i0 + set->ss.tsff.j0 + set->ss.tsff.k0 + set->ss.tsff.i1 + set->ss.tsff.j1 + set->ss.tsff.k1)==0);
        else if (SV_SRCTYPE_LTSF == type)
            return ((set->ss.ltsf.i0 + set->ss.ltsf.j0 + set->ss.ltsf.k0 + set->ss.ltsf.i1 + set->ss.ltsf.j1 + set->ss.ltsf.k1)==0);
        else if (SV_SRCTYPE_LTSFF == type)
            return ((set->ss.ltsff.i0 + set->ss.ltsff.j0 + set->ss.ltsff.k0 + set->ss.ltsff.i1 + set->ss.ltsff.j1 + set->ss.ltsff.k1)==0);

        return TRUE;*/
}

void
remove_temporary_files(SvSet *set)
{
    gint i;

    for (i = 0; i < set->ss.npnts; i++) {
        if (set->ss.pnts[i].source_mode != 0) {
            printf("Removing file %s\n", set->ss.pnts[i].source_filename);
            remove(set->ss.pnts[i].source_filename);
        }
    }
    if (set->ss.sf.valid) {
        if (set->ss.sf.source_mode != 0) {
            printf("Removing file %s\n", set->ss.sf.source_filename);
            remove(set->ss.sf.source_filename);
        }
    }
    if (set->ss.tsf.valid) {
        if (set->ss.tsf.source_mode != 0) {
            printf("Removing file %s\n", set->ss.tsf.source_filename);
            remove(set->ss.tsf.source_filename);
        }
    }
    if (set->ss.tsff.valid) {
        if (set->ss.tsff.source_mode != 0) {
            printf("Removing file %s\n", set->ss.tsff.source_filename);
            remove(set->ss.tsff.source_filename);
        }
    }
    if (set->ss.ltsf.valid) {
        if (set->ss.ltsf.source_mode != 0) {
            printf("Removing file %s\n", set->ss.ltsf.source_filename);
            remove(set->ss.ltsf.source_filename);
        }
    }
    if (set->ss.ltsff.valid) {
        if (set->ss.ltsff.source_mode != 0) {
            printf("Removing file %s\n", set->ss.ltsff.source_filename);
            remove(set->ss.ltsff.source_filename);
        }
    }
}

int
read_settings(gchar *filename, gchar *buffer)
{
    FILE *fh;
    gint size;
    gint result;

    fh = fopen(filename, "r");
    if (fh == NULL) { 
        fprintf(stderr, "Error: Cannot open %s.\n", filename);
        return 0;
    }

    // obtain file size:
    fseek(fh, 0, SEEK_END);
    size = ftell(fh);
    rewind(fh);

    // allocate memory to contain the whole file:
    buffer = (gchar*)g_malloc(sizeof(gchar)*size);
    if (buffer == NULL) { 
        fprintf(stderr, "Error: Cannot allocate buffer to read %s.\n", filename);
    }

    // copy the file into the buffer:
    result = fread(buffer, 1, size, fh);
    if (result != size) { 
        fprintf(stderr, "Error: Cannot read %s.\n", filename);
    }

    return 1;
}

#if(0)
gint
parse_settings_from_string(const gchar *str, SvSet *set, gboolean called_from_gsvit)
{
    if (str == NULL)
        return 0;

    gchar key[100] = {0};
    gchar value[20] = {0};
    gchar buffer[200] = {0};
    gchar tmpsourcefile[512] = {0};
    gint i, j, n, val, thetares, phires, radius;
    gdouble thetafrom, phifrom, thetato, phito, theta, phi;
    gint srctype, savefile, dimensionality;
    gdouble wavelength = SOURCE_WAVELENGTH, dt, amplitude = SOURCE_AMPLITUDE, wlspan = SOURCE_WAVELENGTHSPAN;
    gint pulse_width = SOURCE_PULSE_WIDTH;
    gint ijres, jkres, orientation;
    gdouble ijfrom, jkfrom, ijto, jkto, distance;

    gint n_source_locals_size = 0;
    gint n_abs = 0;

    /*parse file*/
    while (sscanf(str, "%s", key) != EOF) {
        /*if (strstr(key, "#") != NULL) {
            fgets(buffer, 100, fr);
            continue;
        }*/
        //        printf("key: %s\n", key);

        if (strcmp(key, "DIMENSIONALITY") == 0) {
            if (sget_int(str, &dimensionality, key)) {
                goto exit;
            }
            if (dimensionality != 3) {
                fprintf(stderr, "Error: Wrong parameter file dimensionality, it should be 3 for 3D version of Gvit");
                goto exit;
            }
        }
        if (strcmp(key, "VERBOSE") == 0) {
            if (sget_int(str, &(set->sc.verbose), key))
                goto exit;
        } else if (strcmp(key, "POOL") == 0) {
            if (sget_int(str, &(set->sp.xres), key))
                goto exit;
            if (sget_int(str, &(set->sp.yres), key))
                goto exit;
            if (sget_int(str, &(set->sp.zres), key))
                goto exit;
            if (sget_double(str, &(set->sp.dx), key))
                goto exit;
            if (sget_double(str, &(set->sp.dy), key))
                goto exit;
            if (sget_double(str, &(set->sp.dz), key))
                goto exit;
        } else if (strcmp(key, "COMP") == 0) {
            if (sget_int(str, &(set->sc.nsteps), key))
                goto exit;
        } else if (strcmp(key, "DT_MULT") == 0) {
            if (sget_double(str, &(set->sc.dtmult), key))
                goto exit;
        } else if (strcmp(key, "GPU") == 0) {
            if (sget_int(str, &(set->sc.usegpu), key))
                goto exit;
        } else if (strcmp(key, "UGPU") == 0) {
            if (sget_int(str, &val, key))
                goto exit;
            set->sc.ugpu[val] = 1;
        } else if (strcmp(key, "GPU_QUERY") == 0) {
            if (sget_int(str, &(set->sc.devicequery), key))
                goto exit;
        } else if (strcmp(key, "THREADS") == 0) {
            if (sget_int(str, &(set->sc.nthreads), key))
                goto exit;
        } else if (strcmp(key, "MATMODE_CHECK") == 0) {
            if (sget_int(str, &(set->sm.matmode_check), key))
                goto exit;
        } else if (strcmp(key, "MEDIUM_SMOOTH") == 0) {
            if (sget_int(str, &(set->sm.smooth), key))
                goto exit;
        } else if (strcmp(key, "MEDIUM_CROP") == 0) {
            if (sget_int(str, &(set->sm.crop_i0), key))
                goto exit;
            if (sget_int(str, &(set->sm.crop_j0), key))
                goto exit;
            if (sget_int(str, &(set->sm.crop_k0), key))
                goto exit;
            if (sget_int(str, &(set->sm.crop_in), key))
                goto exit;
            if (sget_int(str, &(set->sm.crop_jn), key))
                goto exit;
            if (sget_int(str, &(set->sm.crop_kn), key))
                goto exit;
            if (sget_double(str, &(set->sm.crop_epsilon), key))
                goto exit;
            if (sget_double(str, &(set->sm.crop_mu), key))
                goto exit;
            if (sget_double(str, &(set->sm.crop_sigma), key))
                goto exit;
            if (sget_double(str, &(set->sm.crop_sigast), key))
                goto exit;
            set->sm.crop = TRUE;
        } else if (strcmp(key, "MEDIUM_ROUGHEN") == 0) {
            if (sget_int(str, &(set->sm.rough_radius_peak[set->sm.nroughens]), key))
                goto exit;
            if (sget_int(str, &(set->sm.rough_radius_span[set->sm.nroughens]), key))
                goto exit;
            if (sget_int(str, &(set->sm.rough_iterations[set->sm.nroughens]), key))
                goto exit;
            if (sget_double(str, &(set->sm.rough_probability[set->sm.nroughens]), key))
                goto exit;
            if (sget_int(str, &(set->sm.rough_matindex[set->sm.nroughens]), key))
                goto exit;
            if (sget_int(str, &(set->sm.rough_voidindex[set->sm.nroughens]), key))
                goto exit;
            if (sget_int(str, &(set->sm.rough_seed[set->sm.nroughens]), key))
                goto exit;
            set->sm.nroughens++;
        } else if (strcmp(key, "MEDIUM_GROW") == 0) {
            if (sget_int(str, &(set->sm.grow_i0[set->sm.ngrowths]), key))
                goto exit;
            if (sget_int(str, &(set->sm.grow_j0[set->sm.ngrowths]), key))
                goto exit;
            if (sget_int(str, &(set->sm.grow_k0[set->sm.ngrowths]), key))
                goto exit;
            if (sget_int(str, &(set->sm.grow_in[set->sm.ngrowths]), key))
                goto exit;
            if (sget_int(str, &(set->sm.grow_jn[set->sm.ngrowths]), key))
                goto exit;
            if (sget_int(str, &(set->sm.grow_kn[set->sm.ngrowths]), key))
                goto exit;
            if (sget_int(str, &(set->sm.grow_addindex[set->sm.ngrowths]), key))
                goto exit;
            if (sget_int(str, &(set->sm.grow_attachindex[set->sm.ngrowths]), key))
                goto exit;
            if (sget_int(str, &(set->sm.grow_subsampling[set->sm.ngrowths]), key))
                goto exit;
            if (sget_int(str, &(set->sm.grow_nsteps[set->sm.ngrowths]), key))
                goto exit;
            if (sget_double(str, &(set->sm.grow_mobility[set->sm.ngrowths]), key))
                goto exit;
            if (sget_double(str, &(set->sm.grow_probability[set->sm.ngrowths]), key))
                goto exit;
            if (sget_int(str, &(set->sm.grow_seed[set->sm.ngrowths]), key))
                goto exit;
            set->sm.ngrowths++;
        } else if (strcmp(key, "MEDIUM_GROW_SKIP_FACE") == 0) {
            if (set->sm.ngrowths == 0) {
                fprintf(stderr, "Error: MEDIUM_GROW_SKIP_FACE should follow MEDIUM_GROW direcrive\n");
                goto exit;
            }

            sscanf(str, "%s", value);
            if (strcmp(value, "i0") == 0)
                set->sm.grow_skipi0[set->sm.ngrowths - 1] = 1;
            else if (strcmp(value, "in") == 0)
                set->sm.grow_skipin[set->sm.ngrowths - 1] = 1;
            else if (strcmp(value, "j0") == 0)
                set->sm.grow_skipj0[set->sm.ngrowths - 1] = 1;
            else if (strcmp(value, "jn") == 0)
                set->sm.grow_skipjn[set->sm.ngrowths - 1] = 1;
            else if (strcmp(value, "k0") == 0)
                set->sm.grow_skipk0[set->sm.ngrowths - 1] = 1;
            else if (strcmp(value, "kn") == 0)
                set->sm.grow_skipkn[set->sm.ngrowths - 1] = 1;
        } else if (strcmp(key, "MEDIUM_EXPRESSION") == 0) {
            if (sget_int(str, &(set->sm.expr_i0[set->sm.nexprs]), key))
                goto exit;
            if (sget_int(str, &(set->sm.expr_j0[set->sm.nexprs]), key))
                goto exit;
            if (sget_int(str, &(set->sm.expr_k0[set->sm.nexprs]), key))
                goto exit;
            if (sget_int(str, &(set->sm.expr_in[set->sm.nexprs]), key))
                goto exit;
            if (sget_int(str, &(set->sm.expr_jn[set->sm.nexprs]), key))
                goto exit;
            if (sget_int(str, &(set->sm.expr_kn[set->sm.nexprs]), key))
                goto exit;
            if (sget_int(str, &(set->sm.expr_matindex[set->sm.nexprs]), key))
                goto exit;
            if (sget_int(str, &(set->sm.expr_voidindex[set->sm.nexprs]), key))
                goto exit;
            if (sget_int(str, &(set->sm.expr_maxdist[set->sm.nexprs]), key))
                goto exit;
            if (sget_int(str, &(set->sm.expr_distmode[set->sm.nexprs]), key))
                goto exit;
             sscanf(str, "%200[^\n\r]", buffer);
            set->sm.expr_expr[set->sm.nexprs] = g_strstrip(g_strdup(buffer));
            set->sm.nexprs++;
        } else if (strcmp(key, "MEDIUM_SPECTRAL") == 0) {
            if (sget_double(str, &(set->sm.spectral_sigma[set->sm.nspectrals]), key))
                goto exit;
            if (sget_double(str, &(set->sm.spectral_t[set->sm.nspectrals]), key))
                goto exit;
            if (sget_int(str, &(set->sm.spectral_matindex[set->sm.nspectrals]), key))
                goto exit;
            if (sget_int(str, &(set->sm.spectral_seed[set->sm.nspectrals]), key))
                goto exit;
            set->sm.nspectrals++;
        } else if (strcmp(key, "OUT_FILE") == 0) {
                      //fgets(buffer, 100, str);
            sscanf(str, "%s", buffer);
            sscanf(str, "%100[^\n\r]", buffer);
            set->so.outfile = g_strstrip(g_strdup(buffer));
        } else if (strcmp(key, "MEDIUM_LINEAR") == 0) {
            //fgets(buffer, 100, str);
            sscanf(str, "%s", buffer);
            sscanf(str, "%100[^\n\r]", buffer);
            set->sm.in_voxel_filename = g_strstrip(g_strdup(buffer));
            if (!g_file_test(set->sm.in_voxel_filename, G_FILE_TEST_EXISTS)) {
                fprintf(stderr, "Error: File %s does not exist.\n", set->sm.in_voxel_filename);
                goto exit;
            }
        } else if (strcmp(key, "MEDIUM_VECTOR") == 0) {
            //fgets(buffer, 100, str);
            sscanf(str, "%s", buffer);
            sscanf(str, "%100[^\n\r]", buffer);
            set->sm.in_vector_filename = g_strstrip(g_strdup(buffer));
            if (!g_file_test(set->sm.in_vector_filename, G_FILE_TEST_EXISTS)) {
                fprintf(stderr, "Error: File %s does not exist.\n", set->sm.in_vector_filename);
                goto exit;
            }
        } else if (strcmp(key, "MEDIUM_VTK_DATA") == 0) {
            //fgets(buffer, 100, str);
            sscanf(str, "%s", buffer);
            sscanf(str, "%100[^\n\r]", buffer);
            set->sm.in_vtk_data_filename = g_strstrip(g_strdup(buffer));
            if (!g_file_test(set->sm.in_vtk_data_filename, G_FILE_TEST_EXISTS)) {
                fprintf(stderr, "Error: File %s does not exist.\n", set->sm.in_vtk_data_filename);
                goto exit;
            }
        } else if (strcmp(key, "SUBGRID") == 0) {
            if (sget_int(str, &(set->sg.sg[set->sg.nsg].box_i0), key))
                goto exit;
            if (sget_int(str, &(set->sg.sg[set->sg.nsg].box_j0), key))
                goto exit;
            if (sget_int(str, &(set->sg.sg[set->sg.nsg].box_k0), key))
                goto exit;
            if (sget_int(str, &(set->sg.sg[set->sg.nsg].box_in), key))
                goto exit;
            if (sget_int(str, &(set->sg.sg[set->sg.nsg].box_jn), key))
                goto exit;
            if (sget_int(str, &(set->sg.sg[set->sg.nsg].box_kn), key))
                goto exit;
            if (sget_int(str, &(set->sg.sg[set->sg.nsg].division), key))
                goto exit;
            set->sg.nsg++;
        } else if (strcmp(key, "SOURCE_POINT") == 0) {
            source_point_alloc(set);

            /*load parameters*/
            if (sget_int(str, &(set->ss.pnts[set->ss.npnts].point_origin_position_i), key))
                goto exit;
            if (sget_int(str, &(set->ss.pnts[set->ss.npnts].point_origin_position_j), key))
                goto exit;
            if (sget_int(str, &(set->ss.pnts[set->ss.npnts].point_origin_position_k), key))
                goto exit;
            if (sget_int(str, &(srctype), key))
                goto exit;

            //set->ss.pnts[set->ss.npnts].point_origin_theta = theta;
            //set->ss.pnts[set->ss.npnts].point_origin_phi = phi;

            set->ss.pnts[set->ss.npnts].source_mode = srctype;
            if (srctype == 0) {  //load strom file
                sscanf(str, "%100[^\n\r]", buffer);
                g_free(set->ss.pnts[set->ss.npnts].source_filename);
                set->ss.pnts[set->ss.npnts].source_filename = g_strstrip(g_strdup(buffer));
                if (set->sc.verbose > 1)
                    printf("Source will be loaded from file %s\n", set->ss.pnts[set->ss.npnts].source_filename);
                if ((TRUE == called_from_gsvit) && !g_file_test(set->ss.pnts[set->ss.npnts].source_filename, G_FILE_TEST_EXISTS)) {
                    fprintf(stderr, "Error: File %s does not exist.\n", set->ss.pnts[set->ss.npnts].source_filename);
                    goto exit;
                }
            } else if (srctype == 1 || srctype == 2 || srctype == 3) { //generate simple sources
                if (set->sc.nsteps == 0) {
                    fprintf(stderr, "Error: COMP parameter (number of steps) should precede generated source command.\n");
                    goto exit;
                }
                wlspan = 0;
                if (sget_double(str, &wavelength, key))
                    goto exit;
                if (srctype == 2 && sget_int(str, &pulse_width, key))
                    goto exit;
                if (srctype == 3 && sget_double(str, &wlspan, key))
                    goto exit;
                if (sget_double(str, &amplitude, key))
                    goto exit;

                if (sget_double(str, &theta, key))
                    goto exit;
                if (sget_double(str, &phi, key))
                    goto exit;

                set->ss.pnts[set->ss.npnts].source_amplitude = amplitude;
                set->ss.pnts[set->ss.npnts].source_wl = wavelength;
                set->ss.pnts[set->ss.npnts].source_wlspan = wlspan;
                set->ss.pnts[set->ss.npnts].source_pulsewidth = pulse_width;

                set->ss.pnts[set->ss.npnts].point_origin_theta = theta;
                set->ss.pnts[set->ss.npnts].point_origin_phi = phi;

                dt = 1.0 / LIGHT_SPEED / sqrt(1.0 / set->sp.dx / set->sp.dx + 1.0 / set->sp.dy / set->sp.dy + 1.0 / set->sp.dz / set->sp.dz) * set->sc.dtmult;


                if (called_from_gsvit) {
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "tmp_source_%06d_%06d", set->ss.npnts, set->sc.suffix);
#ifdef G_OS_WIN32
                    /*gchar* filename_tmp = NULL;
                    filename_tmp = g_build_filename(g_get_home_dir(), ".xsvit", tmpsourcefile, NULL);
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "%s", filename_tmp);
                    g_free(filename_tmp);*/
#endif
                    if (write_source(SRCF_ARB, (SrcType)srctype, set->sc.nsteps, wavelength, wlspan, pulse_width, dt, 0, tmpsourcefile, amplitude, theta, phi))
                        goto exit;

                    g_free(set->ss.pnts[set->ss.npnts].source_filename);
                    set->ss.pnts[set->ss.npnts].source_filename = g_strdup(tmpsourcefile);
                }
            } else {
                fprintf(stderr, "Error: unknown source type\n");
                goto exit;
            }

            set->ss.npnts++;
        } else if (strcmp(key, "SOURCE_EXT") == 0) {
            if (sget_int(str, &(set->ss.ext.i), key))
                goto exit;
            if (sget_int(str, &(set->ss.ext.j), key))
                goto exit;
            if (sget_int(str, &(set->ss.ext.k), key))
                goto exit;
            if (sget_int(str, &(set->ss.ext.ijstart), key))
                goto exit;
            if (sget_int(str, &(set->ss.ext.jkstart), key))
                goto exit;
            if (sget_int(str, &(set->ss.ext.extxres), key))
                goto exit;
            if (sget_int(str, &(set->ss.ext.extyres), key))
                goto exit;
            if (sget_int(str, &(set->ss.ext.iextfrom), key))
                goto exit;
            if (sget_int(str, &(set->ss.ext.jextfrom), key))
                goto exit;
            if (sget_int(str, &(set->ss.ext.iextto), key))
                goto exit;
            if (sget_int(str, &(set->ss.ext.jextto), key))
                goto exit;
            if (sget_int(str, &(set->ss.ext.shift), key))
                goto exit;

            sscanf(str, "%s", buffer); /*no spaces allowed here*/
            set->ss.ext.filebase_ex = g_strdup(buffer);
            sscanf(str, "%s", buffer); /*no spaces allowed here*/
            set->ss.ext.filebase_ey = g_strdup(buffer);
            sscanf(str, "%s", buffer); /*no spaces allowed here*/
            set->ss.ext.filebase_ez = g_strdup(buffer);

            sscanf(str, "%s", buffer); /*no spaces allowed here*/
            set->ss.ext.filebase_hx = g_strdup(buffer);
            sscanf(str, "%s", buffer); /*no spaces allowed here*/
            set->ss.ext.filebase_hy = g_strdup(buffer);
            sscanf(str, "%s", buffer); /*no spaces allowed here*/
            set->ss.ext.filebase_hz = g_strdup(buffer);
            set->ss.ext.valid = TRUE;
        } else if (strcmp(key, "SOURCE_LOCAL") == 0) {
            /*alloc or realloc the necessary structure*/
            if (set->ss.nlocals == 0) {
                n_source_locals_size = 10;
                set->ss.locals = (SvSrcLocal *)g_malloc(n_source_locals_size * sizeof(SvSrcLocal));
            } else if (set->ss.nlocals >= n_source_locals_size) { /*realloc*/
                n_source_locals_size += 10;
                set->ss.locals = (SvSrcLocal *)g_realloc(set->ss.locals, n_source_locals_size * sizeof(SvSrcLocal));
            }

            /*load parameters*/
            if (sget_int(str, &(set->ss.locals[set->ss.nlocals].box_i0), key))
                goto exit;
            if (sget_int(str, &(set->ss.locals[set->ss.nlocals].box_j0), key))
                goto exit;
            if (sget_int(str, &(set->ss.locals[set->ss.nlocals].box_k0), key))
                goto exit;
            if (sget_int(str, &(set->ss.locals[set->ss.nlocals].box_in), key))
                goto exit;
            if (sget_int(str, &(set->ss.locals[set->ss.nlocals].box_jn), key))
                goto exit;
            if (sget_int(str, &(set->ss.locals[set->ss.nlocals].box_kn), key))
                goto exit;
            if (sget_double(str, &(set->ss.locals[set->ss.nlocals].layered_epsilon), key))
                goto exit;
            if (sget_double(str, &(set->ss.locals[set->ss.nlocals].density), key))
                goto exit;
            if (sget_double(str, &(set->ss.locals[set->ss.nlocals].strength), key))
                goto exit;
            if (sget_double(str, &(set->ss.locals[set->ss.nlocals].lambda_peak), key))
                goto exit;
            if (sget_double(str, &(set->ss.locals[set->ss.nlocals].lambda_width), key))
                goto exit;
            if (sget_int(str, &(set->ss.locals[set->ss.nlocals].source_mode), key))
                goto exit;
            if (sget_int(str, &(set->ss.locals[set->ss.nlocals].startfrom), key))
                goto exit;

            //set->ss.locals[set->ss.nlocals].layered_epsilon *= EPSILON_0;
            set->ss.nlocals++;
        } else if (strcmp(key, "SOURCE_TSF") == 0) {
            if (sget_int(str, &(set->ss.tsf.box_i0), key))
                goto exit;
            if (sget_int(str, &(set->ss.tsf.box_j0), key))
                goto exit;
            if (sget_int(str, &(set->ss.tsf.box_k0), key))
                goto exit;
            if (sget_int(str, &(set->ss.tsf.box_in), key))
                goto exit;
            if (sget_int(str, &(set->ss.tsf.box_jn), key))
                goto exit;
            if (sget_int(str, &(set->ss.tsf.box_kn), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.ia_theta), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.ia_phi), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.ia_psi), key))
                goto exit;
            if (sget_int(str, &(srctype), key))
                goto exit;
            set->ss.tsf.source_mode = srctype;
            if (srctype == 0) {  //load from file
                sscanf(str, "%100[^\n\r]", buffer);
                g_free(set->ss.tsf.source_filename);
                set->ss.tsf.source_filename = g_strstrip(g_strdup(buffer));
                if (set->sc.verbose > 1)
                    printf("Source will be loaded from file %s\n", set->ss.tsf.source_filename);
            } else if (srctype == 1 || srctype == 2 || srctype == 3) { //generate simple sources
                if (set->sc.nsteps == 0) {
                    fprintf(stderr, "Error: COMP parameter (number of steps) should precede generated source command.\n");
                    goto exit;
                }
                wlspan = 0;
                if (sget_double(str, &wavelength, key))
                    goto exit;
                if (srctype == 2 && sget_int(str, &pulse_width, key))
                    goto exit;
                if (srctype == 3 && sget_double(str, &wlspan, key))
                    goto exit;
                if (sget_double(str, &amplitude, key))
                    goto exit;

                set->ss.tsf.source_wl = wavelength;
                set->ss.tsf.source_wlspan = wlspan;
                set->ss.tsf.source_pulsewidth = pulse_width;
                set->ss.tsf.source_amplitude = amplitude;

                dt = 1.0 / LIGHT_SPEED / sqrt(1.0 / set->sp.dx / set->sp.dx + 1.0 / set->sp.dy / set->sp.dy + 1.0 / set->sp.dz / set->sp.dz) * set->sc.dtmult;

                if (called_from_gsvit) {
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "tmp_source_tsf_%06d", set->sc.suffix);
                    if (write_source(SRCF_TSF, (SrcType)srctype, set->sc.nsteps, wavelength, wlspan, pulse_width, dt, 0, tmpsourcefile, amplitude, 0, 0))
                        goto exit;

                    g_free(set->ss.tsf.source_filename);
                    set->ss.tsf.source_filename = g_strdup(tmpsourcefile);
                }
            } else {
                fprintf(stderr, "Error: unknown source type\n");
                goto exit;
            }
            set->ss.tsf.valid = TRUE;
        } else if (strcmp(key, "TSF_SKIP") == 0) {
            sscanf(str, "%s", value);
            if (strcmp(value, "i0") == 0)
                set->ss.tsf.box_boundary_skipi0 = 1;
            else if (strcmp(value, "in") == 0)
                set->ss.tsf.box_boundary_skipin = 1;
            else if (strcmp(value, "j0") == 0)
                set->ss.tsf.box_boundary_skipj0 = 1;
            else if (strcmp(value, "jn") == 0)
                set->ss.tsf.box_boundary_skipjn = 1;
            else if (strcmp(value, "k0") == 0)
                set->ss.tsf.box_boundary_skipk0 = 1;
            else if (strcmp(value, "kn") == 0)
                set->ss.tsf.box_boundary_skipkn = 1;
            else if (strcmp(value, "depth") == 0) {
                if (sget_int(str, &(set->ss.tsf.box_boundary_skipdepth), key))
                    goto exit;
            }
        } else if (strcmp(key, "TSF_GAUSSIAN_Z") == 0) {
            set->ss.tsf.gaussian = 1;
            if (sget_double(str, &(set->ss.tsf.gaussian_fxpos), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.gaussian_fypos), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.gaussian_rx), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.gaussian_ry), key))
                goto exit;
        } else if (strcmp(key, "TSF_RADIAL_Z") == 0) {
            set->ss.tsf.radial = 1;
            if (sget_double(str, &(set->ss.tsf.radial_fxpos), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.radial_fypos), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.radial_rx), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.radial_ry), key))
                goto exit;
        } else if (strcmp(key, "TSF_FIBER_Z") == 0) {
            set->ss.tsf.fiber = 1;
            if (sget_double(str, &(set->ss.tsf.fiber_fxpos), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.fiber_fypos), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.fiber_radius), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.fiber_cutoff), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.fiber_epsilon_core), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.fiber_epsilon_cladding), key))
                goto exit;
        } else if (strcmp(key, "SOURCE_LTSF") == 0) {
            if (sget_int(str, &(set->ss.ltsf.box_i0), key))
                goto exit;
            if (sget_int(str, &(set->ss.ltsf.box_j0), key))
                goto exit;
            if (sget_int(str, &(set->ss.ltsf.box_k0), key))
                goto exit;
            if (sget_int(str, &(set->ss.ltsf.box_in), key))
                goto exit;
            if (sget_int(str, &(set->ss.ltsf.box_jn), key))
                goto exit;
            if (sget_int(str, &(set->ss.ltsf.box_kn), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsf.ia_theta), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsf.ia_phi), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsf.ia_psi), key))
                goto exit;
            if (sget_int(str, &(set->ss.ltsf.layered_count), key))
                goto exit;
            for (i = 0; i < set->ss.ltsf.layered_count; i++) {
                if (sget_int(str, &(set->ss.ltsf.layered_zpos[i]), key))
                    goto exit;
                if (sget_double(str, &(set->ss.ltsf.layered_epsilon[i]), key))
                    goto exit;
                if (sget_double(str, &(set->ss.ltsf.layered_mu[i]), key))
                    goto exit;
                if (sget_double(str, &(set->ss.ltsf.layered_sigma[i]), key))
                    goto exit;
                if (sget_double(str, &(set->ss.ltsf.layered_sigast[i]), key))
                    goto exit;
            }

            if (sget_int(str, &(srctype), key))
                goto exit;
            set->ss.ltsf.source_mode = srctype;

            if (srctype == 0) {  //load from file
                sscanf(str, "%100[^\n\r]", buffer);
                set->ss.ltsf.source_filename = g_strstrip(g_strdup(buffer));

                if (set->sc.verbose > 1)
                    printf("Source will be loaded from file %s\n", set->ss.ltsf.source_filename);
            } else if (srctype == 1 || srctype == 2 || srctype == 3) { //generate simple sources
                if (set->sc.nsteps == 0) {
                    fprintf(stderr, "Error: COMP parameter (number of steps) should precede generated source command.\n");
                    goto exit;
                }
                wlspan = 0;
                if (sget_double(str, &wavelength, key))
                    goto exit;
                if (srctype == 2 && sget_int(str, &pulse_width, key))
                    goto exit;
                if (srctype == 3 && sget_double(str, &wlspan, key))
                    goto exit;
                if (sget_double(str, &amplitude, key))
                    goto exit;

                set->ss.ltsf.source_wl = wavelength;
                set->ss.ltsf.source_wlspan = wlspan;
                set->ss.ltsf.source_pulsewidth = pulse_width;
                set->ss.ltsf.source_amplitude = amplitude;

                dt = 1.0 / LIGHT_SPEED / sqrt(1.0 / set->sp.dx / set->sp.dx + 1.0 / set->sp.dy / set->sp.dy + 1.0 / set->sp.dz / set->sp.dz) * set->sc.dtmult;

                if (called_from_gsvit) {
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "tmp_source_ltsf_%06d", set->sc.suffix);
                    if (write_source(SRCF_LTSF, (SrcType)srctype, set->sc.nsteps, wavelength, wlspan, pulse_width, dt, 0, tmpsourcefile, amplitude, 0, 0))
                        goto exit;

                    g_free(set->ss.ltsf.source_filename);
                    set->ss.ltsf.source_filename = g_strdup(tmpsourcefile);
                }
            } else {
                fprintf(stderr, "Error: unknown source type\n");
                goto exit;
            }
            set->ss.ltsf.valid = TRUE;
        } else if (strcmp(key, "LTSF_GAUSSIAN") == 0) {
            set->ss.ltsf.gaussian = 1;
            if (sget_double(str, &(set->ss.ltsf.gaussian_fxpos), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsf.gaussian_fypos), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsf.gaussian_rx), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsf.gaussian_ry), key))
                goto exit;
        } else if (strcmp(key, "LTSF_RADIAL") == 0) {
            set->ss.ltsf.radial = 1;
            if (sget_double(str, &(set->ss.ltsf.radial_fxpos), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsf.radial_fypos), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsf.radial_rx), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsf.radial_ry), key))
                goto exit;
        } else if (strcmp(key, "LTSF_FIBER") == 0) {
            set->ss.ltsf.fiber = 1;
            if (sget_double(str, &(set->ss.ltsf.fiber_fxpos), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsf.fiber_fypos), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsf.fiber_radius), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsf.fiber_cutoff), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsf.fiber_epsilon_core), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsf.fiber_epsilon_cladding), key))
                goto exit;
        } else if (strcmp(key, "LTSF_SKIP") == 0) {
            sscanf(str, "%s", value);
            if (strcmp(value, "i0") == 0)
                set->ss.ltsf.box_boundary_skipi0 = 1;
            else if (strcmp(value, "in") == 0)
                set->ss.ltsf.box_boundary_skipin = 1;
            else if (strcmp(value, "j0") == 0)
                set->ss.ltsf.box_boundary_skipj0 = 1;
            else if (strcmp(value, "jn") == 0)
                set->ss.ltsf.box_boundary_skipjn = 1;
            else if (strcmp(value, "k0") == 0)
                set->ss.ltsf.box_boundary_skipk0 = 1;
            else if (strcmp(value, "kn") == 0)
                set->ss.ltsf.box_boundary_skipkn = 1;
            else if (strcmp(value, "depth") == 0) {
                if (sget_int(str, &(set->ss.ltsf.box_boundary_skipdepth), key))
                    goto exit;
            }
        } else if (strcmp(key, "SOURCE_LTSFF") == 0) {
            if (sget_int(str, &(set->ss.ltsff.box_i0), key))
                goto exit;
            if (sget_int(str, &(set->ss.ltsff.box_j0), key))
                goto exit;
            if (sget_int(str, &(set->ss.ltsff.box_k0), key))
                goto exit;
            if (sget_int(str, &(set->ss.ltsff.box_in), key))
                goto exit;
            if (sget_int(str, &(set->ss.ltsff.box_jn), key))
                goto exit;
            if (sget_int(str, &(set->ss.ltsff.box_kn), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsff.focused_thetamax), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsff.focused_fdist), key))
                goto exit;
            if (sget_double(str, &(set->ss.ltsff.focused_pol), key))
                goto exit;
            if (sget_int(str, &(set->ss.ltsff.focused_nip), key))
                goto exit;
            if (sget_int(str, &(set->ss.ltsff.focused_mip), key))
                goto exit;
            if (sget_int(str, &(set->ss.ltsff.layered_count), key))
                goto exit;
            for (i = 0; i < set->ss.ltsff.layered_count; i++) {
                if (sget_int(str, &(set->ss.ltsff.layered_zpos[i]), key))
                    goto exit;
                if (sget_double(str, &(set->ss.ltsff.layered_epsilon[i]), key))
                    goto exit;
                if (sget_double(str, &(set->ss.ltsff.layered_mu[i]), key))
                    goto exit;
                if (sget_double(str, &(set->ss.ltsff.layered_sigma[i]), key))
                    goto exit;
                if (sget_double(str, &(set->ss.ltsff.layered_sigast[i]), key))
                    goto exit;
            }

            if (sget_int(str, &(srctype), key))
                goto exit;
            set->ss.ltsff.source_mode = srctype;

            if (srctype == 0) {  //load from file
                sscanf(str, "%100[^\n\r]", buffer);
                set->ss.ltsff.source_filename = g_strstrip(g_strdup(buffer));
                if (set->sc.verbose > 1)
                    printf("Source will be loaded from file %s\n", set->ss.ltsff.source_filename);
            } else if (srctype == 1 || srctype == 2 || srctype == 3) { //generate simple sources
                if (set->sc.nsteps == 0) {
                    fprintf(stderr, "Error: COMP parameter (number of steps) should precede generated source command.\n");
                    goto exit;
                }
                wlspan = 0;
                if (sget_double(str, &wavelength, key))
                    goto exit;
                if (srctype == 2 && sget_int(str, &pulse_width, key))
                    goto exit;
                if (srctype == 3 && sget_double(str, &wlspan, key))
                    goto exit;
                if (sget_double(str, &amplitude, key))
                    goto exit;

                set->ss.ltsff.source_wl = wavelength;
                set->ss.ltsff.source_wlspan = wlspan;
                set->ss.ltsff.source_pulsewidth = pulse_width;
                set->ss.ltsff.source_amplitude = amplitude;

                dt = 1.0 / LIGHT_SPEED / sqrt(1.0 / set->sp.dx / set->sp.dx + 1.0 / set->sp.dy / set->sp.dy + 1.0 / set->sp.dz / set->sp.dz) * set->sc.dtmult;

                if (called_from_gsvit) {
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "tmp_source_ltsff_%06d", set->sc.suffix);
                    if (write_source(SRCF_LTSF, (SrcType)srctype, set->sc.nsteps, wavelength, wlspan, pulse_width, dt, 0, tmpsourcefile, amplitude, 0, 0))
                        goto exit;

                    g_free(set->ss.ltsff.source_filename);
                    set->ss.ltsff.source_filename = g_strdup(tmpsourcefile);
                }
            } else {
                fprintf(stderr, "Error: unknown source type\n");
                goto exit;
            }
            set->ss.ltsff.valid = TRUE;
        } else if (strcmp(key, "LTSFF_FAST") == 0) {
            if (sget_int(str, &(set->ss.ltsff.fast), key))
                goto exit;
        } else if (strcmp(key, "LTSFF_SKIP") == 0) {
            sscanf(str, "%s", value);
            if (strcmp(value, "i0") == 0)
                set->ss.ltsff.box_boundary_skipi0 = 1;
            else if (strcmp(value, "in") == 0)
                set->ss.ltsff.box_boundary_skipin = 1;
            else if (strcmp(value, "j0") == 0)
                set->ss.ltsff.box_boundary_skipj0 = 1;
            else if (strcmp(value, "jn") == 0)
                set->ss.ltsff.box_boundary_skipjn = 1;
            else if (strcmp(value, "k0") == 0)
                set->ss.ltsff.box_boundary_skipk0 = 1;
            else if (strcmp(value, "kn") == 0)
                set->ss.ltsff.box_boundary_skipkn = 1;
            else if (strcmp(value, "depth") == 0) {
                if (sget_int(str, &(set->ss.ltsff.box_boundary_skipdepth), key))
                    goto exit;
            }
        } else if (strcmp(key, "SOURCE_TSFF") == 0) {
            if (sget_int(str, &(set->ss.tsff.box_i0), key))
                goto exit;
            if (sget_int(str, &(set->ss.tsff.box_j0), key))
                goto exit;
            if (sget_int(str, &(set->ss.tsff.box_k0), key))
                goto exit;
            if (sget_int(str, &(set->ss.tsff.box_in), key))
                goto exit;
            if (sget_int(str, &(set->ss.tsff.box_jn), key))
                goto exit;
            if (sget_int(str, &(set->ss.tsff.box_kn), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsff.focused_thetamax), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsff.focused_fdist), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsff.focused_pol), key))
                goto exit;
            if (sget_int(str, &(set->ss.tsff.focused_nip), key))
                goto exit;
            if (sget_int(str, &(set->ss.tsff.focused_mip), key))
                goto exit;
            if (sget_int(str, &(srctype), key))
                goto exit;
            set->ss.tsff.source_mode = srctype;

            if (srctype == 0) {  //load from file
                sscanf(str, "%100[^\n\r]", buffer);
                set->ss.tsff.source_filename = g_strstrip(g_strdup(buffer));

                if (set->sc.verbose > 1)
                    printf("Source will be loaded from file %s\n", set->ss.tsff.source_filename);
            } else if (srctype == 1 || srctype == 2 || srctype == 3) { //generate simple sources
                if (set->sc.nsteps == 0) {
                    fprintf(stderr, "Error: COMP parameter (number of steps) should precede generated source command.\n");
                    goto exit;
                }
                wlspan = 0;
                if (sget_double(str, &wavelength, key))
                    goto exit;
                if (srctype == 2 && sget_int(str, &pulse_width, key))
                    goto exit;
                if (srctype == 3 && sget_double(str, &wlspan, key))
                    goto exit;
                if (sget_double(str, &amplitude, key))
                    goto exit;

                set->ss.tsff.source_wl = wavelength;
                set->ss.tsff.source_wlspan = wlspan;
                set->ss.tsff.source_pulsewidth = pulse_width;
                set->ss.tsff.source_amplitude = amplitude;

                dt = 1.0 / LIGHT_SPEED / sqrt(1.0 / set->sp.dx / set->sp.dx + 1.0 / set->sp.dy / set->sp.dy + 1.0 / set->sp.dz / set->sp.dz) * set->sc.dtmult;

                if (called_from_gsvit) {
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "tmp_source_tsff_%06d", set->sc.suffix);
                    if (write_source(SRCF_TSF, (SrcType)srctype, set->sc.nsteps, wavelength, wlspan, pulse_width, dt, 0, tmpsourcefile, amplitude, 0, 0))
                        goto exit;

                    g_free(set->ss.tsff.source_filename);
                    set->ss.tsff.source_filename = g_strdup(tmpsourcefile);
                }
            } else {
                fprintf(stderr, "Error: unknown source type\n");
                goto exit;
            }
            set->ss.tsff.valid = TRUE;
        } else if (strcmp(key, "TSFF_SHIFT") == 0) {
            if (sget_double(str, &(set->ss.tsff.xshift), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsff.yshift), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsff.zshift), key))
                goto exit;
        } else if (strcmp(key, "TSFF_SKIP") == 0) {
            sscanf(str, "%s", value);
            if (strcmp(value, "i0") == 0)
                set->ss.tsff.box_boundary_skipi0 = 1;
            else if (strcmp(value, "in") == 0)
                set->ss.tsff.box_boundary_skipin = 1;
            else if (strcmp(value, "j0") == 0)
                set->ss.tsff.box_boundary_skipj0 = 1;
            else if (strcmp(value, "jn") == 0)
                set->ss.tsff.box_boundary_skipjn = 1;
            else if (strcmp(value, "k0") == 0)
                set->ss.tsff.box_boundary_skipk0 = 1;
            else if (strcmp(value, "kn") == 0)
                set->ss.tsff.box_boundary_skipkn = 1;
            else if (strcmp(value, "depth") == 0) {
                if (sget_int(str, &(set->ss.tsff.box_boundary_skipdepth), key))
                    goto exit;
            }
        } else if (strcmp(key, "SOURCE_SF") == 0) {
            if (sget_double(str, &(set->ss.sf.ia_theta), key))
                goto exit;
            if (sget_double(str, &(set->ss.sf.ia_phi), key))
                goto exit;
            if (sget_double(str, &(set->ss.sf.ia_psi), key))
                goto exit;
            if (sget_int(str, &(srctype), key))
                goto exit;
            set->ss.sf.source_mode = srctype;

            if (srctype == 0) {  //load from file
                sscanf(str, "%100[^\n\r]", buffer);
                set->ss.sf.source_filename = g_strstrip(g_strdup(buffer));

                if (set->sc.verbose > 1)
                    printf("Source will be loaded from file %s\n", set->ss.sf.source_filename);
            } else if (srctype == 1 || srctype == 2 || srctype == 3) { //generate simple sources
                if (set->sc.nsteps == 0) {
                    fprintf(stderr, "Error: COMP parameter (number of steps) should precede generated source command.\n");
                    goto exit;
                }
                wlspan = 0;
                if (sget_double(str, &wavelength, key))
                    goto exit;
                if (srctype == 2 && sget_int(str, &pulse_width, key))
                    goto exit;
                if (srctype == 3 && sget_double(str, &wlspan, key))
                    goto exit;
                if (sget_double(str, &amplitude, key))
                    goto exit;

                set->ss.sf.source_wl = wavelength;
                set->ss.sf.source_wlspan = wlspan;
                set->ss.sf.source_pulsewidth = pulse_width;
                set->ss.sf.source_amplitude = amplitude;

                dt = 1.0 / LIGHT_SPEED / sqrt(1.0 / set->sp.dx / set->sp.dx + 1.0 / set->sp.dy / set->sp.dy + 1.0 / set->sp.dz / set->sp.dz) * set->sc.dtmult;

                if (called_from_gsvit) {
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "tmp_source_sf_%06d", set->sc.suffix);

#ifdef G_OS_WIN32
                    /*gchar* filename_tmp = NULL;
                    filename_tmp = g_build_filename(g_get_home_dir(), ".xsvit", tmpsourcefile, NULL);
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "%s", filename_tmp);
                    g_free(filename_tmp);*/
#endif

                    if (write_source(SRCF_TSF, (SrcType)srctype, set->sc.nsteps, wavelength, wlspan, pulse_width, dt, 0, tmpsourcefile, amplitude, 0, 0))
                        goto exit;

                    g_free(set->ss.sf.source_filename);
                    set->ss.sf.source_filename = g_strdup(tmpsourcefile);
                }

                //if (write_source(SRCF_TSF, (SrcType)srctype, set->sc.nsteps, wavelength, wlspan, size, dt, 0, "tmpsource", amplitude, 0, 0)) 
                //    goto exit;                  
                //set->ss.sf.source_filename = g_strdup("tmpsource");                
            } else {
                fprintf(stderr, "Error: unknown source type\n");
                goto exit;
            }
            set->ss.sf.valid = TRUE;
        } else if (strcmp(key, "TSFSOURCE_MATERIAL") == 0) {
            if (sget_double(str, &(set->ss.tsf.layered_epsilon), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.layered_mu), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.layered_sigma), key))
                goto exit;
            if (sget_double(str, &(set->ss.tsf.layered_sigast), key))
                goto exit;
        } else if (strcmp(key, "SOURCE_WAVELENGTH") == 0) {
            if (sget_double(str, &(set->ss.lambda_min), key))
                goto exit;
            if (sget_double(str, &(set->ss.lambda_center), key))
                goto exit;
            if (sget_double(str, &(set->ss.lambda_max), key))
                goto exit;
        } else if (strcmp(key, "OUT_POINT") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_POINT);
            //    	    if (set->so.npnts == 0) {
            //    		    n_output_pnts_size = 10;
            //    		    set->so.pnts = (SvOutputPar *)g_malloc(n_output_pnts_size*sizeof(SvOutputPar));
            //    	    } else if (set->so.npnts >= n_output_pnts_size) { /*realloc*/
            //    		    n_output_pnts_size += 10;

            //    		    set->so.pnts = (SvOutputPar *)g_realloc(set->so.pnts, n_output_pnts_size*sizeof(SvOutputPar));
            //    	    }


            /*load parameters*/
            if (sget_component(str, &(set->so.pnts[set->so.npnts].component), key))
                goto exit;
            if (sget_int(str, &(set->so.pnts[set->so.npnts].step), key))
                goto exit;
            if (sget_int(str, &(set->so.pnts[set->so.npnts].i), key))
                goto exit;
            if (sget_int(str, &(set->so.pnts[set->so.npnts].j), key))
                goto exit;
            if (sget_int(str, &(set->so.pnts[set->so.npnts].k), key))
                goto exit;

            if (set->so.pnts[set->so.npnts].component == SV_COMP_ALLFFT) {
                if (sget_int(str, &(set->so.pnts[set->so.npnts].start), key))
                    goto exit;
                if (sget_int(str, &(set->so.pnts[set->so.npnts].stop), key))
                    goto exit;

                if (set->so.pnts[set->so.npnts].start == -1)
                    set->so.pnts[set->so.npnts].start = 0;
                if (set->so.pnts[set->so.npnts].stop == -1)
                    set->so.pnts[set->so.npnts].stop = set->sc.nsteps;

                if (set->so.pnts[set->so.npnts].start > (set->sc.nsteps - 4)) {
                    set->so.pnts[set->so.npnts].start = set->sc.nsteps - 4;
                    fprintf(stderr, "Warning: requesting FFT output outside of the data range, cropped to %d\n", set->so.pnts[set->so.npnts].start);
                }

                if (set->so.pnts[set->so.npnts].stop > (set->sc.nsteps)) {
                    set->so.pnts[set->so.npnts].stop = set->sc.nsteps;
                    fprintf(stderr, "Warning: requesting FFT output outside of the data range, cropped to %d\n", set->so.pnts[set->so.npnts].stop);
                }

                set->so.savespectrum = 1;
            }

            sscanf(str, "%100[^\n\r]", buffer);
            sprintf(set->so.pnts[set->so.npnts].filebase, "%s", g_strstrip(g_strdup(buffer)));
            set->so.npnts++;
        } else if (strcmp(key, "OUT_IMAGE") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_IMAGE);
            /*load parameters*/
            if (sget_component(str, &(set->so.imgs[set->so.nimgs].component), key))
                goto exit;
            if (sget_int(str, &(set->so.imgs[set->so.nimgs].step), key))
                goto exit;
            if (sget_int(str, &(set->so.imgs[set->so.nimgs].i), key))
                goto exit;
            if (sget_int(str, &(set->so.imgs[set->so.nimgs].j), key))
                goto exit;
            if (sget_int(str, &(set->so.imgs[set->so.nimgs].k), key))
                goto exit;

            sscanf(str, "%100[^\n\r]", buffer);
            sprintf(set->so.imgs[set->so.nimgs].filebase, "%s", g_strstrip(g_strdup(buffer)));

            set->so.nimgs++;
        } else if (strcmp(key, "OUT_SUBGRID_IMAGE") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_SUBGRID_IMAGE);
            /*load parameters*/
            if (sget_component(str, &(set->so.subgrid_imgs[set->so.nsubgrid_imgs].component), key))
                goto exit;
            if (sget_int(str, &(set->so.subgrid_imgs[set->so.nsubgrid_imgs].step), key))
                goto exit;
            if (sget_int(str, &(set->so.subgrid_imgs[set->so.nsubgrid_imgs].i), key))
                goto exit;
            if (sget_int(str, &(set->so.subgrid_imgs[set->so.nsubgrid_imgs].j), key))
                goto exit;
            if (sget_int(str, &(set->so.subgrid_imgs[set->so.nsubgrid_imgs].k), key))
                goto exit;
            if (sget_int(str, &(set->so.subgrid_imgs[set->so.nsubgrid_imgs].format), key)) //this is the numebr of subgrid
                goto exit;

            sscanf(str, "%100[^\n\r]", buffer);
            sprintf(set->so.subgrid_imgs[set->so.nsubgrid_imgs].filebase, "%s", g_strstrip(g_strdup(buffer)));

            set->so.nsubgrid_imgs++;
        } else if (strcmp(key, "OUT_PLANE") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_PLANE);
            //    	    if (set->so.nplns == 0) {
            //    		    n_output_plns_size = 10;
            //    		    set->so.plns = (SvOutputPar *)g_malloc(n_output_plns_size*sizeof(SvOutputPar));
            //    	    } else if (set->so.nplns >= n_output_plns_size) { /*realloc*/
            //    		    n_output_plns_size += 10;
            //    		    set->so.plns = (SvOutputPar *)g_realloc(set->so.plns, n_output_plns_size*sizeof(SvOutputPar));
            //    	    }

            /*load parameters*/
            if (sget_component(str, &(set->so.plns[set->so.nplns].component), key))
                goto exit;
            if (sget_int(str, &(set->so.plns[set->so.nplns].step), key))
                goto exit;
            if (sget_int(str, &(set->so.plns[set->so.nplns].start), key))
                goto exit;
            if (sget_int(str, &(set->so.plns[set->so.nplns].stop), key))
                goto exit;
            if (sget_int(str, &(set->so.plns[set->so.nplns].format), key))
                goto exit; //here it means binary/ascii
            if (sget_int(str, &(set->so.plns[set->so.nplns].i), key))
                goto exit;
            if (sget_int(str, &(set->so.plns[set->so.nplns].j), key))
                goto exit;
            if (sget_int(str, &(set->so.plns[set->so.nplns].k), key))
                goto exit;

            sscanf(str, "%100[^\n\r]", buffer);
            sprintf(set->so.plns[set->so.nplns].filebase, "%s", g_strstrip(g_strdup(buffer)));

            set->so.nplns++;
        } else if (strcmp(key, "OUT_VOLUME") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_VOLUME);
            //    	    if (set->so.ncubs == 0) {
            //    		    n_output_cubs_size = 10;
            //    		    set->so.cubs = (SvOutputPar *)g_malloc(n_output_cubs_size*sizeof(SvOutputPar));
            //    	    } else if (set->so.ncubs >= n_output_cubs_size) { /*realloc*/
            //    		    n_output_cubs_size += 10;
            //    		    set->so.cubs = (SvOutputPar *)g_realloc(set->so.cubs, n_output_cubs_size*sizeof(SvOutputPar));
            //    	    }

            /*load parameters*/
            if (sget_cube_component(str, &(set->so.cubs[set->so.ncubs].component), key))
                goto exit;
            if (sget_int(str, &(set->so.cubs[set->so.ncubs].step), key))
                goto exit;
            if (sget_int(str, &(set->so.cubs[set->so.ncubs].start), key))
                goto exit;
            if (sget_int(str, &(set->so.cubs[set->so.ncubs].stop), key))
                goto exit;
            if (sget_int(str, &(set->so.cubs[set->so.ncubs].format), key))
                goto exit; //here it means binary/ascii

            if (set->so.cubs[set->so.ncubs].start == -1)
                set->so.cubs[set->so.ncubs].start = 0;
            if (set->so.cubs[set->so.ncubs].stop == -1)
                set->so.cubs[set->so.ncubs].stop = set->sc.nsteps;

            sscanf(str, "%100[^\n\r]", buffer);
            sprintf(set->so.cubs[set->so.ncubs].filebase, "%s", g_strstrip(g_strdup(buffer)));

            if (set->so.cubs[set->so.ncubs].component == SV_OVOLUME_ABS || set->so.cubs[set->so.ncubs].component == SV_OVOLUME_SUMALL || set->so.cubs[set->so.ncubs].component == SV_OVOLUME_MAXALL)
                set->so.cubs[set->so.ncubs].k = n_abs++; //position of data 

            set->so.ncubs++;
        } else if (strcmp(key, "OUT_SUM") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_SUM);
            //    	    if (set->so.nsums == 0) {
            //    		    n_output_sums_size = 10;
            //    		    set->so.sums = (SvOutputSum *)g_malloc(n_output_sums_size*sizeof(SvOutputSum));
            //    	    } else if (set->so.nsums >= n_output_sums_size) { /*realloc*/
            //    		    n_output_sums_size += 10;
            //    		    set->so.sums = (SvOutputSum *)g_realloc(set->so.sums, n_output_sums_size*sizeof(SvOutputSum));
            //    	    }

            /*load parameters*/
            if (sget_scomponent(str, &(set->so.sums[set->so.nsums].component), key))
                goto exit;
            if (sget_int(str, &(set->so.sums[set->so.nsums].step), key))
                goto exit;
            if (sget_int(str, &(set->so.sums[set->so.nsums].box_i0), key))
                goto exit;
            if (sget_int(str, &(set->so.sums[set->so.nsums].box_j0), key))
                goto exit;
            if (sget_int(str, &(set->so.sums[set->so.nsums].box_k0), key))
                goto exit;
            if (sget_int(str, &(set->so.sums[set->so.nsums].box_in), key))
                goto exit;
            if (sget_int(str, &(set->so.sums[set->so.nsums].box_jn), key))
                goto exit;
            if (sget_int(str, &(set->so.sums[set->so.nsums].box_kn), key))
                goto exit;
            if (sget_double(str, &(set->so.sums[set->so.nsums].layered_epsilon), key))
                goto exit;
            if (sget_double(str, &(set->so.sums[set->so.nsums].layered_mu), key))
                goto exit;
            if (sget_double(str, &(set->so.sums[set->so.nsums].layered_sigma), key))
                goto exit;
            if (sget_double(str, &(set->so.sums[set->so.nsums].layered_sigast), key))
                goto exit;
            sscanf(str, "%100[^\n\r]", buffer);
            sprintf(set->so.sums[set->so.nsums].filename, "%s", g_strstrip(g_strdup(buffer)));

            /*if (set->so.sums[set->so.nsums].layered_epsilon != -1)
            set->so.sums[set->so.nsums].layered_epsilon *= EPSILON_0;
            set->so.sums[set->so.nsums].layered_mu *= MU_0;*/

            set->so.sums[set->so.nsums].stringbased = 0;
            set->so.nsums++;
        } else if (strcmp(key, "OUT_SUMTAB") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_SUMTAB);
            //    	    if (set->so.nsums == 0) {
            //    		    n_output_sums_size = 10;
            //    		    set->so.sums = (SvOutputSum *)g_malloc(n_output_sums_size*sizeof(SvOutputSum));
            //    	    } else if (set->so.nsums >= n_output_sums_size) { /*realloc*/
            //    		    n_output_sums_size += 10;
            //    		    set->so.sums = (SvOutputSum *)g_realloc(set->so.sums, n_output_sums_size*sizeof(SvOutputSum));
            //    	    }

            /*load parameters*/
            if (sget_scomponent(str, &(set->so.sums[set->so.nsums].component), key))
                goto exit;
            if (sget_int(str, &(set->so.sums[set->so.nsums].step), key))
                goto exit;
            if (sget_int(str, &(set->so.sums[set->so.nsums].box_i0), key))
                goto exit;
            if (sget_int(str, &(set->so.sums[set->so.nsums].box_j0), key))
                goto exit;
            if (sget_int(str, &(set->so.sums[set->so.nsums].box_k0), key))
                goto exit;
            if (sget_int(str, &(set->so.sums[set->so.nsums].box_in), key))
                goto exit;
            if (sget_int(str, &(set->so.sums[set->so.nsums].box_jn), key))
                goto exit;
            if (sget_int(str, &(set->so.sums[set->so.nsums].box_kn), key))
                goto exit;
            sscanf(str, "%s", buffer);
            sprintf(set->so.sums[set->so.nsums].string, "%s", buffer);

            sscanf(str, "%100[^\n\r]", buffer);
            sprintf(set->so.sums[set->so.nsums].filename, "%s", g_strstrip(g_strdup(buffer)));

            set->so.sums[set->so.nsums].stringbased = 1;
            set->so.nsums++;
        } else if (strcmp(key, "OUT_FORCE") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_FORCE);
            //if (set->so.nforces == 0) {
            //    n_output_forces_size = 10;
            //    set->so.forces = (SvOutputForce *)g_malloc(n_output_forces_size*sizeof(SvOutputForce));
            //} else if (set->so.nforces >= n_output_forces_size) { /*realloc*/
            //    n_output_forces_size += 10;
            //    set->so.forces = (SvOutputForce *)g_realloc(set->so.forces, n_output_forces_size*sizeof(SvOutputForce));
            //}

            /*load parameters*/
            if (sget_int(str, &(set->so.forces[set->so.nforces].step), key))
                goto exit;
            if (sget_int(str, &(set->so.forces[set->so.nforces].box_i0), key))
                goto exit;
            if (sget_int(str, &(set->so.forces[set->so.nforces].box_j0), key))
                goto exit;
            if (sget_int(str, &(set->so.forces[set->so.nforces].box_k0), key))
                goto exit;
            if (sget_int(str, &(set->so.forces[set->so.nforces].box_in), key))
                goto exit;
            if (sget_int(str, &(set->so.forces[set->so.nforces].box_jn), key))
                goto exit;
            if (sget_int(str, &(set->so.forces[set->so.nforces].box_kn), key))
                goto exit;
            sscanf(str, "%100[^\n\r]", buffer);
            sprintf(set->so.forces[set->so.nforces].filename, "%s", g_strstrip(g_strdup(buffer)));
            set->so.nforces++;
        } else if (strcmp(key, "MBOUNDARY_X0") == 0) {
            sscanf(str, "%s", value);
            if (sget_int(str, &(set->smb.bx0pos), key))
                goto exit;
            if (strcmp(value, "periodic") == 0) {
                set->smb.bx0 = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "MBOUNDARY_XN") == 0) {
            sscanf(str, "%s", value);
            if (sget_int(str, &(set->smb.bxnpos), key)) {
                goto exit;
            }
            if (strcmp(value, "periodic") == 0) {
                set->smb.bxn = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "MBOUNDARY_Y0") == 0) {
            sscanf(str, "%s", value);
            if (sget_int(str, &(set->smb.by0pos), key)) {
                goto exit;
            }
            if (strcmp(value, "periodic") == 0) {
                set->smb.by0 = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "MBOUNDARY_YN") == 0) {
            sscanf(str, "%s", value);
            if (sget_int(str, &(set->smb.bynpos), key))
                goto exit;
            if (strcmp(value, "periodic") == 0) {
                set->smb.byn = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "MBOUNDARY_Z0") == 0) {
            sscanf(str, "%s", value);
            if (sget_int(str, &(set->smb.bz0pos), key))
                goto exit;
            if (strcmp(value, "periodic") == 0) {
                set->smb.bz0 = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "MBOUNDARY_ZN") == 0) {
            sscanf(str, "%s", value);
            if (sget_int(str, &(set->smb.bznpos), key))
                goto exit;
            if (strcmp(value, "periodic") == 0) {
                set->smb.bzn = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "MBOUNDARY_ALL") == 0) {
            sscanf(str, "%s", value);
            if (sget_int(str, &(set->smb.bx0pos), key))
                goto exit;
            if (sget_int(str, &(set->smb.bxnpos), key))
                goto exit;
            if (sget_int(str, &(set->smb.by0pos), key))
                goto exit;
            if (sget_int(str, &(set->smb.bynpos), key))
                goto exit;
            if (sget_int(str, &(set->smb.bz0pos), key))
                goto exit;
            if (sget_int(str, &(set->smb.bznpos), key))
                goto exit;
            if (strcmp(value, "periodic") == 0) {
                set->smb.bx0 = set->smb.bxn = set->smb.by0 = set->smb.byn = set->smb.bz0 = set->smb.bzn = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "BOUNDARY_ALL") == 0) {
            sscanf(str, "%s", value);
            if (strcmp(value, "none") == 0) {
                set->sb.bx0 = set->sb.bxn = set->sb.by0 = set->sb.byn = set->sb.bz0 = set->sb.bzn = SV_BOUNDARY_NONE;
            } else if (strcmp(value, "pec") == 0) {
                set->sb.bx0 = set->sb.bxn = set->sb.by0 = set->sb.byn = set->sb.bz0 = set->sb.bzn = SV_BOUNDARY_PEC;
            } else if (strcmp(value, "liao") == 0) {
                set->sb.bx0 = set->sb.bxn = set->sb.by0 = set->sb.byn = set->sb.bz0 = set->sb.bzn = SV_BOUNDARY_LIAO;
            } else if (strcmp(value, "cpml") == 0) {
                set->sb.bx0 = set->sb.bxn = set->sb.by0 = set->sb.byn = set->sb.bz0 = set->sb.bzn = SV_BOUNDARY_CPML;

                if (sget_int(str, &(set->sb.depth_bx0), key))
                    goto exit;
                set->sb.depth_bxn = set->sb.depth_by0 = set->sb.depth_byn = set->sb.depth_bz0 = set->sb.depth_bzn = set->sb.depth_bx0;

                if (sget_int(str, &(set->sb.m_bx0), key))
                    goto exit;
                set->sb.m_bxn = set->sb.m_by0 = set->sb.m_byn = set->sb.m_bz0 = set->sb.m_bzn = set->sb.m_bx0;

                if (sget_double(str, &(set->sb.sigma_bx0), key))
                    goto exit;
                if ((TRUE == called_from_gsvit) && (set->sb.sigma_bx0 == -1)) {
                    set->sb.sigma_bx0 = 0.8 * (gdouble)set->sb.m_bx0 / (377 * set->sp.dx);
                    if (set->sc.verbose > 1)
                        printf("CPML sigma set to optimum value of %g\n", set->sb.sigma_bx0);
                }
                set->sb.sigma_bxn = set->sb.sigma_by0 = set->sb.sigma_byn = set->sb.sigma_bz0 = set->sb.sigma_bzn = set->sb.sigma_bx0;

                if (sget_double(str, &(set->sb.a_bx0), key))
                    goto exit;
                set->sb.a_bxn = set->sb.a_by0 = set->sb.a_byn = set->sb.a_bz0 = set->sb.a_bzn = set->sb.a_bx0;

                if (sget_double(str, &(set->sb.kappa_bx0), key))
                    goto exit;
                set->sb.kappa_bxn = set->sb.kappa_by0 = set->sb.kappa_byn = set->sb.kappa_bz0 = set->sb.kappa_bzn = set->sb.kappa_bx0;
            } else {
                fprintf(stderr, "Error: unsupported boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "BOUNDARY_X0") == 0) {
            sscanf(str, "%s", value);
            if (strcmp(value, "none") == 0)
                set->sb.bx0 = SV_BOUNDARY_NONE;
            else if (strcmp(value, "pec") == 0)
                set->sb.bx0 = SV_BOUNDARY_PEC;
            else if (strcmp(value, "liao") == 0)
                set->sb.bx0 = SV_BOUNDARY_LIAO;
            else if (strcmp(value, "cpml") == 0) {
                set->sb.bx0 = SV_BOUNDARY_CPML;

                if (sget_int(str, &(set->sb.depth_bx0), key))
                    goto exit;
                if (sget_int(str, &(set->sb.m_bx0), key))
                    goto exit;
                if (sget_double(str, &(set->sb.sigma_bx0), key))
                    goto exit;
                if ((TRUE == called_from_gsvit) && (set->sb.sigma_bx0 == -1)) {
                    set->sb.sigma_bx0 = 0.8 * (gdouble)set->sb.m_bx0 / (377 * set->sp.dx);
                    if (set->sc.verbose > 1)
                        printf("CPML sigma set to optimum value of %g\n", set->sb.sigma_bx0);
                }
                if (sget_double(str, &(set->sb.a_bx0), key))
                    goto exit;
                if (sget_double(str, &(set->sb.kappa_bx0), key))
                    goto exit;
            } else {
                fprintf(stderr, "Error: unsupported boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "BOUNDARY_XN") == 0) {
            sscanf(str, "%s", value);
            if (strcmp(value, "none") == 0)
                set->sb.bxn = SV_BOUNDARY_NONE;
            else if (strcmp(value, "pec") == 0)
                set->sb.bxn = SV_BOUNDARY_PEC;
            else if (strcmp(value, "liao") == 0)
                set->sb.bxn = SV_BOUNDARY_LIAO;
            else if (strcmp(value, "cpml") == 0) {
                set->sb.bxn = SV_BOUNDARY_CPML;

                if (sget_int(str, &(set->sb.depth_bxn), key))
                    goto exit;
                if (sget_int(str, &(set->sb.m_bxn), key))
                    goto exit;
                if (sget_double(str, &(set->sb.sigma_bxn), key))
                    goto exit;
                if ((TRUE == called_from_gsvit) && (set->sb.sigma_bxn == -1)) {
                    set->sb.sigma_bxn = 0.8 * (gdouble)set->sb.m_bxn / (377 * set->sp.dx);
                    if (set->sc.verbose > 1)
                        printf("CPML sigma set to optimum value of %g\n", set->sb.sigma_bxn);
                }
                if (sget_double(str, &(set->sb.a_bxn), key))
                    goto exit;
                if (sget_double(str, &(set->sb.kappa_bxn), key))
                    goto exit;
            } else {
                fprintf(stderr, "Error: unsupported boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "BOUNDARY_Y0") == 0) {
            sscanf(str, "%s", value);
            if (strcmp(value, "none") == 0)
                set->sb.by0 = SV_BOUNDARY_NONE;
            else if (strcmp(value, "pec") == 0)
                set->sb.by0 = SV_BOUNDARY_PEC;
            else if (strcmp(value, "liao") == 0)
                set->sb.by0 = SV_BOUNDARY_LIAO;
            else if (strcmp(value, "cpml") == 0) {
                set->sb.by0 = SV_BOUNDARY_CPML;

                if (sget_int(str, &(set->sb.depth_by0), key))
                    goto exit;
                if (sget_int(str, &(set->sb.m_by0), key))
                    goto exit;
                if (sget_double(str, &(set->sb.sigma_by0), key))
                    goto exit;
                if ((TRUE == called_from_gsvit) && (set->sb.sigma_by0 == -1)) {
                    set->sb.sigma_by0 = 0.8 * (gdouble)set->sb.m_by0 / (377 * set->sp.dy);
                    if (set->sc.verbose > 1)
                        printf("CPML sigma set to optimum value of %g\n", set->sb.sigma_by0);
                }
                if (sget_double(str, &(set->sb.a_by0), key))
                    goto exit;
                if (sget_double(str, &(set->sb.kappa_by0), key))
                    goto exit;
            } else {
                fprintf(stderr, "Error: unsupported boundary condition (%s)\n", value);
                goto exit;
            }

        } else if (strcmp(key, "BOUNDARY_YN") == 0) {
            sscanf(str, "%s", value);
            if (strcmp(value, "none") == 0)
                set->sb.byn = SV_BOUNDARY_NONE;
            else if (strcmp(value, "pec") == 0)
                set->sb.byn = SV_BOUNDARY_PEC;
            else if (strcmp(value, "liao") == 0)
                set->sb.byn = SV_BOUNDARY_LIAO;
            else if (strcmp(value, "cpml") == 0) {
                set->sb.byn = SV_BOUNDARY_CPML;

                if (sget_int(str, &(set->sb.depth_byn), key))
                    goto exit;
                if (sget_int(str, &(set->sb.m_byn), key))
                    goto exit;
                if (sget_double(str, &(set->sb.sigma_byn), key))
                    goto exit;
                if ((TRUE == called_from_gsvit) && (set->sb.sigma_byn == -1)) {
                    set->sb.sigma_byn = 0.8 * (gdouble)set->sb.m_byn / (377 * set->sp.dy);
                    if (set->sc.verbose > 1)
                        printf("CPML sigma set to optimum value of %g\n", set->sb.sigma_byn);
                }
                if (sget_double(str, &(set->sb.a_byn), key))
                    goto exit;
                if (sget_double(str, &(set->sb.kappa_byn), key))
                    goto exit;
            } else {
                fprintf(stderr, "Error: unsupported boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "BOUNDARY_Z0") == 0) {
            sscanf(str, "%s", value);
            if (strcmp(value, "none") == 0)
                set->sb.bz0 = SV_BOUNDARY_NONE;
            else if (strcmp(value, "pec") == 0)
                set->sb.bz0 = SV_BOUNDARY_PEC;
            else if (strcmp(value, "liao") == 0)
                set->sb.bz0 = SV_BOUNDARY_LIAO;
            else if (strcmp(value, "cpml") == 0) {
                set->sb.bz0 = SV_BOUNDARY_CPML;

                if (sget_int(str, &(set->sb.depth_bz0), key))
                    goto exit;
                if (sget_int(str, &(set->sb.m_bz0), key))
                    goto exit;
                if (sget_double(str, &(set->sb.sigma_bz0), key))
                    goto exit;
                if ((TRUE == called_from_gsvit) && (set->sb.sigma_bz0 == -1)) {
                    set->sb.sigma_bz0 = 0.8 * (gdouble)set->sb.m_bz0 / (377 * set->sp.dz);
                    if (set->sc.verbose > 1)
                        printf("CPML sigma set to optimum value of %g\n", set->sb.sigma_bz0);
                }
                if (sget_double(str, &(set->sb.a_bz0), key))
                    goto exit;
                if (sget_double(str, &(set->sb.kappa_bz0), key))
                    goto exit;
            } else {
                fprintf(stderr, "Error: unsupported boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "BOUNDARY_ZN") == 0) {
            sscanf(str, "%s", value);
            if (strcmp(value, "none") == 0)
                set->sb.bzn = SV_BOUNDARY_NONE;
            else if (strcmp(value, "pec") == 0)
                set->sb.bzn = SV_BOUNDARY_PEC;
            else if (strcmp(value, "liao") == 0)
                set->sb.bzn = SV_BOUNDARY_LIAO;
            else if (strcmp(value, "cpml") == 0) {
                set->sb.bzn = SV_BOUNDARY_CPML;

                if (sget_int(str, &(set->sb.depth_bzn), key))
                    goto exit;
                if (sget_int(str, &(set->sb.m_bzn), key))
                    goto exit;
                if (sget_double(str, &(set->sb.sigma_bzn), key))
                    goto exit;
                if ((TRUE == called_from_gsvit) && (set->sb.sigma_bzn == -1)) {
                    set->sb.sigma_bzn = 0.8 * (gdouble)set->sb.m_bzn / (377 * set->sp.dz);
                    if (set->sc.verbose > 1)
                        printf("CPML sigma set to optimum value of %g\n", set->sb.sigma_bzn);
                }
                if (sget_double(str, &(set->sb.a_bzn), key))
                    goto exit;
                if (sget_double(str, &(set->sb.kappa_bzn), key))
                    goto exit;
            } else {
                fprintf(stderr, "Error: unsupported boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "NFFF") == 0) {
            if (sget_int(str, &(set->sf.box_i0), key))
                goto exit;
            if (sget_int(str, &(set->sf.box_j0), key))
                goto exit;
            if (sget_int(str, &(set->sf.box_k0), key))
                goto exit;
            if (sget_int(str, &(set->sf.box_in), key))
                goto exit;
            if (sget_int(str, &(set->sf.box_jn), key))
                goto exit;
            if (sget_int(str, &(set->sf.box_kn), key))
                goto exit;
            if (sget_int(str, &(set->sf.sumfrom), key))
                goto exit;
            if (sget_int(str, &(set->sf.sumto), key))
                goto exit;

            //set->sf.valid = TRUE;
        } else if (strcmp(key, "NFFF_SKIP") == 0) {
            sscanf(str, "%s", value);
            if (strcmp(value, "i0") == 0) {
                set->sf.box_boundary_skipi0 = 1;
                if (sget_int(str, &(set->sf.skipi0_jmin), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipi0_kmin), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipi0_jmax), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipi0_kmax), key))
                    goto exit;
            } else if (strcmp(value, "in") == 0) {
                set->sf.box_boundary_skipin = 1;
                if (sget_int(str, &(set->sf.skipin_jmin), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipin_kmin), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipin_jmax), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipin_kmax), key))
                    goto exit;
            } else if (strcmp(value, "j0") == 0) {
                set->sf.box_boundary_skipj0 = 1;
                if (sget_int(str, &(set->sf.skipj0_imin), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipj0_kmin), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipj0_imax), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipj0_kmax), key))
                    goto exit;
            } else if (strcmp(value, "jn") == 0) {
                set->sf.box_boundary_skipjn = 1;
                if (sget_int(str, &(set->sf.skipjn_imin), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipjn_kmin), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipjn_imax), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipjn_kmax), key))
                    goto exit;
            } else if (strcmp(value, "k0") == 0) {
                set->sf.box_boundary_skipk0 = 1;
                if (sget_int(str, &(set->sf.skipk0_imin), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipk0_jmin), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipk0_imax), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipk0_jmax), key))
                    goto exit;
            } else if (strcmp(value, "kn") == 0) {
                set->sf.box_boundary_skipkn = 1;
                if (sget_int(str, &(set->sf.skipkn_imin), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipkn_jmin), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipkn_imax), key))
                    goto exit;
                if (sget_int(str, &(set->sf.skipkn_jmax), key))
                    goto exit;
            }
        } else if (strcmp(key, "NFFF_RAMAHI_POINT") == 0) {
            if (set->sf.nrs == 0) {
                set->sf.nrs = 1;
                set->sf.ri = (gint *)g_malloc(set->sf.nrs * sizeof(gint));
                set->sf.rj = (gint *)g_malloc(set->sf.nrs * sizeof(gint));
                set->sf.rk = (gint *)g_malloc(set->sf.nrs * sizeof(gint));
                set->sf.individual = (gint *)g_malloc(set->sf.nrs * sizeof(gint));
                set->sf.source_filename = (gchar **)g_malloc(set->sf.nrs * sizeof(gchar *));
            } else { /*realloc*/
                set->sf.nrs += 1;
                set->sf.ri = (gint *)g_realloc(set->sf.ri, set->sf.nrs * sizeof(gint));
                set->sf.rj = (gint *)g_realloc(set->sf.rj, set->sf.nrs * sizeof(gint));
                set->sf.rk = (gint *)g_realloc(set->sf.rk, set->sf.nrs * sizeof(gint));
                set->sf.individual = (gint *)g_realloc(set->sf.individual, set->sf.nrs * sizeof(gint));
                set->sf.source_filename = (gchar **)g_realloc(set->sf.source_filename, set->sf.nrs * sizeof(gchar *));
            }
            if (sget_int(str, &(set->sf.ri[set->sf.nrs - 1]), key))
                goto exit;
            if (sget_int(str, &(set->sf.rj[set->sf.nrs - 1]), key))
                goto exit;
            if (sget_int(str, &(set->sf.rk[set->sf.nrs - 1]), key))
                goto exit;
            sscanf(str, "%100[^\n\r]", buffer);
            set->sf.source_filename[set->sf.nrs - 1] = g_strstrip(g_strdup(buffer));
            set->sf.individual[set->sf.nrs - 1] = 0;
        } else if (strcmp(key, "PERIODIC_NFFF") == 0) {
            if (sget_int(str, &(set->spf.box_i0), key))
                goto exit;
            if (sget_int(str, &(set->spf.box_j0), key))
                goto exit;
            if (sget_int(str, &(set->spf.box_k0), key))
                goto exit;
            if (sget_int(str, &(set->spf.box_in), key))
                goto exit;
            if (sget_int(str, &(set->spf.box_jn), key))
                goto exit;
            if (sget_int(str, &(set->spf.box_kn), key))
                goto exit;
            if (sget_int(str, &(set->spf.pimin), key))
                goto exit;
            if (sget_int(str, &(set->spf.pjmin), key))
                goto exit;
            if (sget_int(str, &(set->spf.pimax), key))
                goto exit;
            if (sget_int(str, &(set->spf.pjmax), key))
                goto exit;

            //set->spf.valid = TRUE;
        } else if (strcmp(key, "PERIODIC_NFFF_POSTPROCESS") == 0) {
            if (sget_int(str, &(set->spf.postprocess), key))
                goto exit;
            if (sget_int(str, &(set->spf.ppstart), key))
                goto exit;
        } else if (strcmp(key, "PERIODIC_NFFF_SKIP") == 0) {
            sscanf(str, "%s", value);
            if (strcmp(value, "k0") == 0) {
                set->spf.box_boundary_skipk0 = 1;
                if (sget_int(str, &(set->spf.skipk0_imin), key))
                    goto exit;
                if (sget_int(str, &(set->spf.skipk0_jmin), key))
                    goto exit;
                if (sget_int(str, &(set->spf.skipk0_imax), key))
                    goto exit;
                if (sget_int(str, &(set->spf.skipk0_jmax), key))
                    goto exit;
            } else if (strcmp(value, "kn") == 0) {
                set->spf.box_boundary_skipkn = 1;
                if (sget_int(str, &(set->spf.skipkn_imin), key))
                    goto exit;
                if (sget_int(str, &(set->spf.skipkn_jmin), key))
                    goto exit;
                if (sget_int(str, &(set->spf.skipkn_imax), key))
                    goto exit;
                if (sget_int(str, &(set->spf.skipkn_jmax), key))
                    goto exit;
            }
        } else if (strcmp(key, "PERIODIC_NFFF_RAMAHI_POINT") == 0) {
            if (set->spf.nrs == 0) {
                set->spf.nrs = 1;
                set->spf.ri = (gint *)g_malloc(set->spf.nrs * sizeof(gint));
                set->spf.rj = (gint *)g_malloc(set->spf.nrs * sizeof(gint));
                set->spf.rk = (gint *)g_malloc(set->spf.nrs * sizeof(gint));
                set->spf.individual = (gint *)g_malloc(set->spf.nrs * sizeof(gint));
                set->spf.source_filename = (gchar **)g_malloc(set->spf.nrs * sizeof(gchar *));
            } else { /*realloc*/
                set->spf.nrs += 1;
                set->spf.ri = (gint *)g_realloc(set->spf.ri, set->spf.nrs * sizeof(gint));
                set->spf.rj = (gint *)g_realloc(set->spf.rj, set->spf.nrs * sizeof(gint));
                set->spf.rk = (gint *)g_realloc(set->spf.rk, set->spf.nrs * sizeof(gint));
                set->spf.individual = (gint *)g_realloc(set->spf.individual, set->spf.nrs * sizeof(gint));
                set->spf.source_filename = (gchar **)g_realloc(set->spf.source_filename, set->spf.nrs * sizeof(gchar *));
            }
            if (sget_int(str, &(set->spf.ri[set->spf.nrs - 1]), key))
                goto exit;
            if (sget_int(str, &(set->spf.rj[set->spf.nrs - 1]), key))
                goto exit;
            if (sget_int(str, &(set->spf.rk[set->spf.nrs - 1]), key))
                goto exit;
            sscanf(str, "%100[^\n\r]", buffer);
            set->spf.source_filename[set->spf.nrs - 1] = g_strstrip(g_strdup(buffer));
            set->spf.individual[set->spf.nrs - 1] = 0; /*not part of any set*/
        } else if (strcmp(key, "PERIODIC_NFFF_SPHERICAL_AREA") == 0) {

            if (sget_int(str, &thetares, key))
                goto exit;
            if (sget_int(str, &phires, key))
                goto exit;
            if (sget_int(str, &radius, key))
                goto exit;
            if (sget_double(str, &thetafrom, key))
                goto exit;
            if (sget_double(str, &phifrom, key))
                goto exit;
            if (sget_double(str, &thetato, key))
                goto exit;
            if (sget_double(str, &phito, key))
                goto exit;
            if (sget_int(str, &savefile, key))
                goto exit;

            if (set->spf.nareas < 20) {  //for xsvit only
                set->spf.area_thetares[set->spf.nareas] = thetares;
                set->spf.area_phires[set->spf.nareas] = phires;
                set->spf.area_radius[set->spf.nareas] = radius;
                set->spf.area_thetafrom[set->spf.nareas] = thetafrom;
                set->spf.area_phifrom[set->spf.nareas] = phifrom;
                set->spf.area_thetato[set->spf.nareas] = thetato;
                set->spf.area_phito[set->spf.nareas] = phito;
                set->spf.area_savefile[set->spf.nareas] = savefile;
                set->spf.nareas++;
            }

            if (TRUE == called_from_gsvit) {
                if (set->spf.nrs == 0) {
                    set->spf.nrs = thetares*phires;
                    set->spf.ri = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->spf.rj = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->spf.rk = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->spf.individual = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->spf.source_filename = (gchar **)g_malloc(thetares*phires * sizeof(gchar *));
                } else {
                    set->spf.nrs += thetares*phires;
                    set->spf.ri = (gint *)g_realloc(set->spf.ri, set->spf.nrs * sizeof(gint));
                    set->spf.rj = (gint *)g_realloc(set->spf.rj, set->spf.nrs * sizeof(gint));
                    set->spf.rk = (gint *)g_realloc(set->spf.rk, set->spf.nrs * sizeof(gint));
                    set->spf.individual = (gint *)g_realloc(set->spf.individual, set->spf.nrs * sizeof(gint));
                    set->spf.source_filename = (gchar **)g_realloc(set->spf.source_filename, set->spf.nrs * sizeof(gchar *));
                }

                n = set->spf.nrs - thetares*phires;

                if (set->spf.nsets == 0) {
                    set->spf.nsets = 1;    /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                    set->spf.setxres = (gint *)g_malloc((set->spf.nsets + 1) * sizeof(gint));
                    set->spf.setyres = (gint *)g_malloc((set->spf.nsets + 1) * sizeof(gint));
                } else {
                    set->spf.nsets += 1;
                    set->spf.setxres = (gint *)g_realloc(set->spf.setxres, (set->spf.nsets + 1) * sizeof(gint));
                    set->spf.setyres = (gint *)g_realloc(set->spf.setyres, (set->spf.nsets + 1) * sizeof(gint));
                }

                set->spf.setxres[set->spf.nsets] = thetares;
                set->spf.setyres[set->spf.nsets] = phires;

                for (i = 0; i < thetares; i++) {
                    for (j = 0; j < phires; j++) {
                        if (thetares > 1)
                            theta = (thetafrom + (gdouble)(i)*(thetato - thetafrom) / (thetares - 1));
                        else
                            theta = thetafrom;
                        if (phires > 1)
                            phi = (phifrom + (gdouble)(j)*(phito - phifrom) / (phires - 1));
                        else
                            phi = phifrom;

                        set->spf.ri[n] = (gint)(radius*sin(theta)*cos(phi)) + set->sp.xres / 2;
                        set->spf.rj[n] = (gint)(radius*sin(theta)*sin(phi)) + set->sp.yres / 2;
                        set->spf.rk[n] = (gint)(radius*cos(theta)) + set->sp.zres / 2;

                        set->spf.individual[n] = set->spf.nsets; /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                        if (!savefile)
                            set->spf.source_filename[n] = NULL;
                        else
                            set->spf.source_filename[n] = g_strdup_printf("xnff_set%03d_%03d_%03d", set->spf.nsets, i, j);

                        n++;
                    } // j
                } // i
            }
        } else if (strcmp(key, "NFFF_SPHERICAL_AREA") == 0) {
            if (sget_int(str, &thetares, key))
                goto exit;
            if (sget_int(str, &phires, key))
                goto exit;
            if (sget_int(str, &radius, key))
                goto exit;
            if (sget_double(str, &thetafrom, key))
                goto exit;
            if (sget_double(str, &phifrom, key))
                goto exit;
            if (sget_double(str, &thetato, key))
                goto exit;
            if (sget_double(str, &phito, key))
                goto exit;
            if (sget_int(str, &savefile, key))
                goto exit;

            if (set->sf.nareas < 20) {  //for xsvit only
                set->sf.area_thetares[set->sf.nareas] = thetares;
                set->sf.area_phires[set->sf.nareas] = phires;
                set->sf.area_radius[set->sf.nareas] = radius;
                set->sf.area_thetafrom[set->sf.nareas] = thetafrom;
                set->sf.area_phifrom[set->sf.nareas] = phifrom;
                set->sf.area_thetato[set->sf.nareas] = thetato;
                set->sf.area_phito[set->sf.nareas] = phito;
                set->sf.area_savefile[set->sf.nareas] = savefile;
                set->sf.nareas++;
            }

            if (TRUE == called_from_gsvit) {
                if (set->sf.nrs == 0) {
                    set->sf.nrs = thetares*phires;
                    set->sf.ri = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->sf.rj = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->sf.rk = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->sf.individual = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->sf.source_filename = (gchar **)g_malloc(thetares*phires * sizeof(gchar *));
                } else {
                    set->sf.nrs += thetares*phires;
                    set->sf.ri = (gint *)g_realloc(set->sf.ri, set->sf.nrs * sizeof(gint));
                    set->sf.rj = (gint *)g_realloc(set->sf.rj, set->sf.nrs * sizeof(gint));
                    set->sf.rk = (gint *)g_realloc(set->sf.rk, set->sf.nrs * sizeof(gint));
                    set->sf.individual = (gint *)g_realloc(set->sf.individual, set->sf.nrs * sizeof(gint));
                    set->sf.source_filename = (gchar **)g_realloc(set->sf.source_filename, set->sf.nrs * sizeof(gchar *));
                }
                n = set->sf.nrs - thetares*phires;

                if (set->sf.nsets == 0) {
                    set->sf.nsets = 1;    /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                    set->sf.setxres = (gint *)g_malloc((set->sf.nsets + 1) * sizeof(gint));
                    set->sf.setyres = (gint *)g_malloc((set->sf.nsets + 1) * sizeof(gint));
                } else {
                    set->sf.nsets += 1;
                    set->sf.setxres = (gint *)g_realloc(set->sf.setxres, (set->sf.nsets + 1) * sizeof(gint));
                    set->sf.setyres = (gint *)g_realloc(set->sf.setyres, (set->sf.nsets + 1) * sizeof(gint));
                }
                set->sf.setxres[set->sf.nsets] = thetares;
                set->sf.setyres[set->sf.nsets] = phires;

                for (i = 0; i < thetares; i++) {
                    for (j = 0; j < phires; j++) {
                        if (thetares > 1)
                            theta = (thetafrom + (gdouble)(i)*(thetato - thetafrom) / (thetares - 1));
                        else
                            theta = thetafrom;
                        if (phires > 1)
                            phi = (phifrom + (gdouble)(j)*(phito - phifrom) / (phires - 1));
                        else
                            phi = phifrom;

                        set->sf.ri[n] = (gint)(radius*sin(theta)*cos(phi)) + set->sp.xres / 2;
                        set->sf.rj[n] = (gint)(radius*sin(theta)*sin(phi)) + set->sp.yres / 2;
                        set->sf.rk[n] = (gint)(radius*cos(theta)) + set->sp.zres / 2;

                        set->sf.individual[n] = set->sf.nsets; /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                        if (!savefile)
                            set->sf.source_filename[n] = NULL;
                        else
                            set->sf.source_filename[n] = g_strdup_printf("xnff_set%03d_%03d_%03d", set->sf.nsets, i, j);
                        n++;
                    } // j
                } // i
            }
        } else if (strcmp(key, "NFFF_PLANAR_AREA") == 0) {
            if (sget_int(str, &ijres, key))
                goto exit;
            if (sget_int(str, &jkres, key))
                goto exit;
            if (sget_double(str, &ijfrom, key))
                goto exit;
            if (sget_double(str, &jkfrom, key))
                goto exit;
            if (sget_double(str, &ijto, key))
                goto exit;
            if (sget_double(str, &jkto, key))
                goto exit;
            if (sget_int(str, &orientation, key))
                goto exit;
            if (sget_double(str, &distance, key))
                goto exit;
            if (sget_int(str, &savefile, key))
                goto exit;

            // printf("or %d\n", orientation);
            if (set->sf.nsquares < 20) {  //for xsvit only
                set->sf.square_ijres[set->sf.nsquares] = thetares;
                set->sf.square_jkres[set->sf.nsquares] = phires;
                set->sf.square_distance[set->sf.nsquares] = distance;
                set->sf.square_ijfrom[set->sf.nsquares] = thetafrom;
                set->sf.square_jkfrom[set->sf.nsquares] = phifrom;
                set->sf.square_ijto[set->sf.nsquares] = thetato;
                set->sf.square_jkto[set->sf.nsquares] = phito;
                set->sf.square_orientation[set->sf.nsquares] = orientation;
                set->sf.square_savefile[set->sf.nsquares] = savefile;
                set->sf.nsquares++;
            }

            if (TRUE == called_from_gsvit) {
                if (set->sf.nrs == 0) {
                    set->sf.nrs = ijres*jkres;
                    set->sf.ri = (gint *)g_malloc(ijres*jkres * sizeof(gint));
                    set->sf.rj = (gint *)g_malloc(ijres*jkres * sizeof(gint));
                    set->sf.rk = (gint *)g_malloc(ijres*jkres * sizeof(gint));
                    set->sf.individual = (gint *)g_malloc(ijres*jkres * sizeof(gint));
                    set->sf.source_filename = (gchar **)g_malloc(ijres*jkres * sizeof(gchar *));
                } else {
                    set->sf.nrs += ijres*jkres;
                    set->sf.ri = (gint *)g_realloc(set->sf.ri, set->sf.nrs * sizeof(gint));
                    set->sf.rj = (gint *)g_realloc(set->sf.rj, set->sf.nrs * sizeof(gint));
                    set->sf.rk = (gint *)g_realloc(set->sf.rk, set->sf.nrs * sizeof(gint));
                    set->sf.individual = (gint *)g_realloc(set->sf.individual, set->sf.nrs * sizeof(gint));
                    set->sf.source_filename = (gchar **)g_realloc(set->sf.source_filename, set->sf.nrs * sizeof(gchar *));
                }
                n = set->sf.nrs - ijres*jkres;

                if (set->sf.nsets == 0) {
                    set->sf.nsets = 1;    /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                    set->sf.setxres = (gint *)g_malloc((set->sf.nsets + 1) * sizeof(gint));
                    set->sf.setyres = (gint *)g_malloc((set->sf.nsets + 1) * sizeof(gint));
                } else {
                    set->sf.nsets += 1;
                    set->sf.setxres = (gint *)g_realloc(set->sf.setxres, (set->sf.nsets + 1) * sizeof(gint));
                    set->sf.setyres = (gint *)g_realloc(set->sf.setyres, (set->sf.nsets + 1) * sizeof(gint));
                }

                set->sf.setxres[set->sf.nsets] = ijres;
                set->sf.setyres[set->sf.nsets] = jkres;

                for (i = 0; i < ijres; i++) {
                    for (j = 0; j < jkres; j++) {
                        if (orientation == 0) { //x orientation
                            set->sf.ri[n] = (gint)(distance / set->sp.dx) + set->sp.xres / 2;
                            set->sf.rj[n] = (gint)((ijfrom + (ijto - ijfrom) * ((gdouble)i) / ((gdouble)ijres - 1)) / set->sp.dy) + set->sp.yres / 2;
                            set->sf.rk[n] = (gint)((jkfrom + (jkto - jkfrom) * ((gdouble)j) / ((gdouble)jkres - 1)) / set->sp.dz) + set->sp.zres / 2;

                            //                      printf("OUT_POINT\nAll 1 %d %d %d point_%.3d_%.3d\n\n", set->sf.ri[n], set->sf.rj[n], set->sf.rk[n], i, j); 
                            //                      printf("NFFF_RAMAHI_POINT\n%d %d %d fpoint_%.3d_%.3d\n\n", set->sf.ri[n], set->sf.rj[n], set->sf.rk[n], i, j); 

                            //printf("%d %d %d\n", set->sf.ri[n], set->sf.rj[n], set->sf.rk[n]);

                            set->sf.individual[n] = set->sf.nsets; /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                            if (!savefile)
                                set->sf.source_filename[n] = NULL;
                            else
                                set->sf.source_filename[n] = g_strdup_printf("xnff_set%03d_%03d_%03d", set->sf.nsets, i, j);
                            n++;
                        } else if (orientation == 1) { //y orientation
                            set->sf.ri[n] = (gint)((ijfrom + (ijto - ijfrom) * ((gdouble)i) / ((gdouble)ijres - 1)) / set->sp.dx) + set->sp.xres / 2;
                            set->sf.rj[n] = (gint)(distance / set->sp.dy) + set->sp.yres / 2;
                            set->sf.rk[n] = (gint)((jkfrom + (jkto - jkfrom) * ((gdouble)j) / ((gdouble)jkres - 1)) / set->sp.dz) + set->sp.zres / 2;

                            set->sf.individual[n] = set->sf.nsets; /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                            if (!savefile)
                                set->sf.source_filename[n] = NULL;
                            else
                                set->sf.source_filename[n] = g_strdup_printf("xnff_set%03d_%03d_%03d", set->sf.nsets, i, j);
                            n++;
                        } else if (orientation == 2) { //z orientation
                            set->sf.ri[n] = (gint)((ijfrom + (ijto - ijfrom) * ((gdouble)i) / ((gdouble)ijres - 1)) / set->sp.dx) + set->sp.xres / 2;
                            set->sf.rj[n] = (gint)((jkfrom + (jkto - jkfrom) * ((gdouble)j) / ((gdouble)jkres - 1)) / set->sp.dy) + set->sp.yres / 2;
                            set->sf.rk[n] = (gint)(distance / set->sp.dz) + set->sp.zres / 2;

                            set->sf.individual[n] = set->sf.nsets; /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                            if (!savefile)
                                set->sf.source_filename[n] = NULL;
                            else
                                set->sf.source_filename[n] = g_strdup_printf("xnff_set%03d_%03d_%03d", set->sf.nsets, i, j);
                            n++;
                        } else {
                            fprintf(stderr, "Error: unsupported NFFF area orientation (%d)\n", orientation);
                            goto exit;
                        }
                    } // j
                } // i
            }
        } else {
            fprintf(stderr, "Error: unsupported or unknown key (%s)\n", key);
            goto exit;
        }
    }

    if (set->smb.bxnpos == -1)
        set->smb.bxnpos = set->sp.xres;
    if (set->smb.bynpos == -1)
        set->smb.bynpos = set->sp.yres;
    if (set->smb.bznpos == -1)
        set->smb.bznpos = set->sp.zres;

    return 1;

exit:
    return 0;
}


gint
parse_settings_from_file(gchar *filename, SvSet *set, gboolean called_from_gsvit)
{
    gint result = 0;
    gchar *str = NULL;

    read_settings(filename, str);
    result = parse_settings_from_string(str, set, called_from_gsvit);

    g_free(str);

    return result;
}
#endif


#if(1)
int
parse_settings(gchar *filename, SvSet *set, gboolean called_from_gsvit)
{
    gchar key[100] = {0};
    gchar value[20] = {0};
    gchar buffer[200] = {0};
    gchar tmpsourcefile[512] = {0};
    gint i, j, n, val, thetares, phires, radius;
    gdouble thetafrom, phifrom, thetato, phito, theta, phi;
    gint srctype, savefile, dimensionality;
    gdouble wavelength = SOURCE_WAVELENGTH, dt, amplitude = SOURCE_AMPLITUDE, wlspan = SOURCE_WAVELENGTHSPAN;
    gint pulse_width = SOURCE_PULSE_WIDTH;
    gint ijres, jkres, orientation;
    gdouble ijfrom, jkfrom, ijto, jkto, distance;

    FILE *fr = fopen(filename, "r");
    gint n_source_locals_size = 0;
    gint n_abs = 0;

    if (!fr)
        return 1;

    //clear_settings(set);

    /*parse file*/
    while (fscanf(fr, "%s", key) != EOF) {
        if (strstr(key, "#") != NULL) {
            fgets(buffer, 100, fr);
            continue;
        }
        //        printf("key: %s\n", key);

        if (strcmp(key, "DIMENSIONALITY") == 0) {
            if (get_int(fr, &dimensionality, key)) {
                goto exit;
            }
            if (dimensionality != 3) {
                fprintf(stderr, "Error: Wrong parameter file dimensionality, it should be 3 for 3D version of Gvit\n");
                goto exit;
            }
        }
        if (strcmp(key, "VERBOSE") == 0) {
            if (get_int(fr, &(set->sc.verbose), key))
                goto exit;
        } else if (strcmp(key, "POOL") == 0) {
            if (get_int(fr, &(set->sp.xres), key))
                goto exit;
            if (get_int(fr, &(set->sp.yres), key))
                goto exit;
            if (get_int(fr, &(set->sp.zres), key))
                goto exit;
            if (get_double(fr, &(set->sp.dx), key))
                goto exit;
            if (get_double(fr, &(set->sp.dy), key))
                goto exit;
            if (get_double(fr, &(set->sp.dz), key))
                goto exit;
        } else if (strcmp(key, "COMP") == 0) {
            if (get_int(fr, &(set->sc.nsteps), key))
                goto exit;
        } else if (strcmp(key, "DT_MULT") == 0) {
            if (get_double(fr, &(set->sc.dtmult), key))
                goto exit;
        } else if (strcmp(key, "GPU") == 0) {
            if (get_int(fr, &(set->sc.usegpu), key))
                goto exit;
        } else if (strcmp(key, "UGPU") == 0) {
            if (get_int(fr, &val, key))
                goto exit;
            set->sc.ugpu[val] = 1;
        } else if (strcmp(key, "GPU_QUERY") == 0) {
            if (get_int(fr, &(set->sc.devicequery), key))
                goto exit;
        } else if (strcmp(key, "THREADS") == 0) {
            if (get_int(fr, &(set->sc.nthreads), key))
                goto exit;
        } else if (strcmp(key, "MATMODE_CHECK") == 0) {
            if (get_int(fr, &(set->sm.matmode_check), key))
                goto exit;
        } else if (strcmp(key, "MEDIUM_SMOOTH") == 0) {
            if (get_int(fr, &(set->sm.smooth), key))
                goto exit;
        } else if (strcmp(key, "MEDIUM_CROP") == 0) {
            if (get_int(fr, &(set->sm.crop_i0), key))
                goto exit;
            if (get_int(fr, &(set->sm.crop_j0), key))
                goto exit;
            if (get_int(fr, &(set->sm.crop_k0), key))
                goto exit;
            if (get_int(fr, &(set->sm.crop_in), key))
                goto exit;
            if (get_int(fr, &(set->sm.crop_jn), key))
                goto exit;
            if (get_int(fr, &(set->sm.crop_kn), key))
                goto exit;
            if (get_double(fr, &(set->sm.crop_epsilon), key))
                goto exit;
            if (get_double(fr, &(set->sm.crop_mu), key))
                goto exit;
            if (get_double(fr, &(set->sm.crop_sigma), key))
                goto exit;
            if (get_double(fr, &(set->sm.crop_sigast), key))
                goto exit;
            set->sm.crop = TRUE;        
        } else if (strcmp(key, "MEDIUM_GROW") == 0) {
            if (get_int(fr, &(set->sm.grow_i0[set->sm.ngrowths]), key))
                goto exit;
            if (get_int(fr, &(set->sm.grow_j0[set->sm.ngrowths]), key))
                goto exit;
            if (get_int(fr, &(set->sm.grow_k0[set->sm.ngrowths]), key))
                goto exit;
            if (get_int(fr, &(set->sm.grow_in[set->sm.ngrowths]), key))
                goto exit;
            if (get_int(fr, &(set->sm.grow_jn[set->sm.ngrowths]), key))
                goto exit;
            if (get_int(fr, &(set->sm.grow_kn[set->sm.ngrowths]), key))
                goto exit;
            if (get_int(fr, &(set->sm.grow_addindex[set->sm.ngrowths]), key))
                goto exit;
            if (get_int(fr, &(set->sm.grow_attachindex[set->sm.ngrowths]), key))
                goto exit;
            if (get_int(fr, &(set->sm.grow_subsampling[set->sm.ngrowths]), key))
                goto exit;
            if (get_int(fr, &(set->sm.grow_nsteps[set->sm.ngrowths]), key))
                goto exit;
            if (get_double(fr, &(set->sm.grow_mobility[set->sm.ngrowths]), key))
                goto exit;
            if (get_double(fr, &(set->sm.grow_probability[set->sm.ngrowths]), key))
                goto exit;
            if (get_int(fr, &(set->sm.grow_seed[set->sm.ngrowths]), key))
                goto exit;
            set->sm.ngrowths++;
        } else if (strcmp(key, "MEDIUM_GROW_SKIP_FACE") == 0) {
            if (set->sm.ngrowths == 0) {
                fprintf(stderr, "Error: MEDIUM_GROW_SKIP_FACE should follow MEDIUM_GROW direcrive\n");
                goto exit;
            }

            fscanf(fr, "%s", value);
            if (strcmp(value, "i0") == 0)
                set->sm.grow_skipi0[set->sm.ngrowths - 1] = 1;
            else if (strcmp(value, "in") == 0)
                set->sm.grow_skipin[set->sm.ngrowths - 1] = 1;
            else if (strcmp(value, "j0") == 0)
                set->sm.grow_skipj0[set->sm.ngrowths - 1] = 1;
            else if (strcmp(value, "jn") == 0)
                set->sm.grow_skipjn[set->sm.ngrowths - 1] = 1;
            else if (strcmp(value, "k0") == 0)
                set->sm.grow_skipk0[set->sm.ngrowths - 1] = 1;
            else if (strcmp(value, "kn") == 0)
                set->sm.grow_skipkn[set->sm.ngrowths - 1] = 1;
        } else if (strcmp(key, "MEDIUM_ROUGHEN") == 0) {
            if (get_int(fr, &(set->sm.rough_radius_peak[set->sm.nroughens]), key))
                goto exit;
            if (get_int(fr, &(set->sm.rough_radius_span[set->sm.nroughens]), key))
                goto exit;
            if (get_int(fr, &(set->sm.rough_iterations[set->sm.nroughens]), key))
                goto exit;
            if (get_double(fr, &(set->sm.rough_probability[set->sm.nroughens]), key))
                goto exit;
            if (get_int(fr, &(set->sm.rough_matindex[set->sm.nroughens]), key))
                goto exit;
            if (get_int(fr, &(set->sm.rough_voidindex[set->sm.nroughens]), key))
                goto exit;
            if (get_int(fr, &(set->sm.rough_seed[set->sm.nroughens]), key))
                goto exit;
            set->sm.nroughens++;
        } else if (strcmp(key, "MEDIUM_SPECTRAL") == 0) {
            if (get_double(fr, &(set->sm.spectral_sigma[set->sm.nspectrals]), key))
                goto exit;
            if (get_double(fr, &(set->sm.spectral_t[set->sm.nspectrals]), key))
                goto exit;
            if (get_int(fr, &(set->sm.spectral_matindex[set->sm.nspectrals]), key))
                goto exit;
            if (get_int(fr, &(set->sm.spectral_seed[set->sm.nspectrals]), key))
                goto exit;
            set->sm.nspectrals++;
        } else if (strcmp(key, "MEDIUM_EXPRESSION") == 0) {
            if (get_int(fr, &(set->sm.expr_i0[set->sm.nexprs]), key))
                goto exit;
            if (get_int(fr, &(set->sm.expr_j0[set->sm.nexprs]), key))
                goto exit;
            if (get_int(fr, &(set->sm.expr_k0[set->sm.nexprs]), key))
                goto exit;
            if (get_int(fr, &(set->sm.expr_in[set->sm.nexprs]), key))
                goto exit;
            if (get_int(fr, &(set->sm.expr_jn[set->sm.nexprs]), key))
                goto exit;
            if (get_int(fr, &(set->sm.expr_kn[set->sm.nexprs]), key))
                goto exit;
            if (get_int(fr, &(set->sm.expr_matindex[set->sm.nexprs]), key))
                goto exit;
            if (get_int(fr, &(set->sm.expr_voidindex[set->sm.nexprs]), key))
                goto exit;
            if (get_int(fr, &(set->sm.expr_maxdist[set->sm.nexprs]), key))
                goto exit;
            if (get_int(fr, &(set->sm.expr_distmode[set->sm.nexprs]), key))
                goto exit;
             fscanf(fr, "%200[^\n\r]", buffer);
            //set->sm.expr_expr[set->sm.nexprs] = g_strstrip(g_strdup(buffer));
            g_snprintf(set->sm.expr_expr[set->sm.nexprs], sizeof(gchar)*MATERIAL_EXPR_EXPR_CHARS, g_strstrip(g_strdup(buffer)));
            set->sm.nexprs++;        
         } else if (strcmp(key, "OUT_FILE") == 0) {
            fgets(buffer, 100, fr);
            fscanf(fr, "%100[^\n\r]", buffer);
            set->so.outfile = g_strstrip(g_strdup(buffer));
        } else if (strcmp(key, "MEDIUM_LINEAR") == 0) {
            fgets(buffer, 100, fr);
            fscanf(fr, "%100[^\n\r]", buffer);
            set->sm.in_voxel_filename = g_strstrip(g_strdup(buffer));
            if (!g_file_test(set->sm.in_voxel_filename, G_FILE_TEST_EXISTS)) {
                fprintf(stderr, "Error: File %s does not exist.\n", set->sm.in_voxel_filename);
                goto exit;
            }
        } else if (strcmp(key, "MEDIUM_VECTOR") == 0) {
            fgets(buffer, 100, fr);
            fscanf(fr, "%100[^\n\r]", buffer);
            set->sm.in_vector_filename = g_strstrip(g_strdup(buffer));
            if (!g_file_test(set->sm.in_vector_filename, G_FILE_TEST_EXISTS)) {
                fprintf(stderr, "Error: File %s does not exist.\n", set->sm.in_vector_filename);
                goto exit;
            }
        } else if (strcmp(key, "MEDIUM_VTK_DATA") == 0) {
            fgets(buffer, 100, fr);
            fscanf(fr, "%100[^\n\r]", buffer);
            set->sm.in_vtk_data_filename = g_strstrip(g_strdup(buffer));
            if (!g_file_test(set->sm.in_vtk_data_filename, G_FILE_TEST_EXISTS)) {
                fprintf(stderr, "Error: File %s does not exist.\n", set->sm.in_vtk_data_filename);
                goto exit;
            }
        } 
        else if (strcmp(key, "SUBGRID") == 0) {
            if (get_int(fr, &(set->sg.sg[set->sg.nsg].box_i0), key))
                goto exit;
            if (get_int(fr, &(set->sg.sg[set->sg.nsg].box_j0), key))
                goto exit;
            if (get_int(fr, &(set->sg.sg[set->sg.nsg].box_k0), key))
                goto exit;
            if (get_int(fr, &(set->sg.sg[set->sg.nsg].box_in), key))
                goto exit;
            if (get_int(fr, &(set->sg.sg[set->sg.nsg].box_jn), key))
                goto exit;
            if (get_int(fr, &(set->sg.sg[set->sg.nsg].box_kn), key))
                goto exit;
            if (get_int(fr, &(set->sg.sg[set->sg.nsg].division), key))
                goto exit;
            set->sg.nsg++;
        }
         else if (strcmp(key, "SOURCE_POINT") == 0) {
            source_point_alloc(set);

            /*load parameters*/
            if (get_int(fr, &(set->ss.pnts[set->ss.npnts].point_origin_position_i), key))
                goto exit;
            if (get_int(fr, &(set->ss.pnts[set->ss.npnts].point_origin_position_j), key))
                goto exit;
            if (get_int(fr, &(set->ss.pnts[set->ss.npnts].point_origin_position_k), key))
                goto exit;
            if (get_int(fr, &(srctype), key))
                goto exit;

            //set->ss.pnts[set->ss.npnts].point_origin_theta = theta;
            //set->ss.pnts[set->ss.npnts].point_origin_phi = phi;

            set->ss.pnts[set->ss.npnts].source_mode = srctype;
            if (srctype == 0) {  //load from file
                fscanf(fr, "%100[^\n\r]", buffer);
                g_free(set->ss.pnts[set->ss.npnts].source_filename);
                set->ss.pnts[set->ss.npnts].source_filename = g_strstrip(g_strdup(buffer));
                if (set->sc.verbose > 1)
                    printf("Source will be loaded from file %s\n", set->ss.pnts[set->ss.npnts].source_filename);
                if ((TRUE == called_from_gsvit) && !g_file_test(set->ss.pnts[set->ss.npnts].source_filename, G_FILE_TEST_EXISTS)) {
                    fprintf(stderr, "Error: File %s does not exist.\n", set->ss.pnts[set->ss.npnts].source_filename);
                    goto exit;
                }
            } else if (srctype == 1 || srctype == 2 || srctype == 3) { //generate simple sources
                if (set->sc.nsteps == 0) {
                    fprintf(stderr, "Error: COMP parameter (number of steps) should precede generated source command.\n");
                    goto exit;
                }
                wlspan = 0;
                if (get_double(fr, &wavelength, key))
                    goto exit;
                if (srctype == 2 && get_int(fr, &pulse_width, key))
                    goto exit;
                if (srctype == 3 && get_double(fr, &wlspan, key))
                    goto exit;
                if (get_double(fr, &amplitude, key))
                    goto exit;

                if (get_double(fr, &theta, key))
                    goto exit;
                if (get_double(fr, &phi, key))
                    goto exit;

                set->ss.pnts[set->ss.npnts].source_amplitude = amplitude;
                set->ss.pnts[set->ss.npnts].source_wl = wavelength;
                set->ss.pnts[set->ss.npnts].source_wlspan = wlspan;
                set->ss.pnts[set->ss.npnts].source_pulsewidth = pulse_width;

                set->ss.pnts[set->ss.npnts].point_origin_theta = theta;
                set->ss.pnts[set->ss.npnts].point_origin_phi = phi;

                dt = 1.0 / LIGHT_SPEED / sqrt(1.0 / set->sp.dx / set->sp.dx + 1.0 / set->sp.dy / set->sp.dy + 1.0 / set->sp.dz / set->sp.dz) * set->sc.dtmult;


                if (called_from_gsvit) {
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "tmp_source_%06d_%06d", set->ss.npnts, set->sc.suffix);
#ifdef G_OS_WIN32
                    /*gchar* filename_tmp = NULL;
                    filename_tmp = g_build_filename(g_get_home_dir(), ".xsvit", tmpsourcefile, NULL);
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "%s", filename_tmp);
                    g_free(filename_tmp);*/
#endif
                    if (write_source(SRCF_ARB, (SrcType)srctype, set->sc.nsteps, wavelength, wlspan, pulse_width, dt, 0, tmpsourcefile, amplitude, theta, phi))
                        goto exit;

                    g_free(set->ss.pnts[set->ss.npnts].source_filename);
                    set->ss.pnts[set->ss.npnts].source_filename = g_strdup(tmpsourcefile);
                }
            } else {
                fprintf(stderr, "Error: unknown source type\n");
                goto exit;
            }

            set->ss.npnts++;
        } else if (strcmp(key, "SOURCE_EXT") == 0) {
            if (get_int(fr, &(set->ss.ext.i), key))
                goto exit;
            if (get_int(fr, &(set->ss.ext.j), key))
                goto exit;
            if (get_int(fr, &(set->ss.ext.k), key))
                goto exit;
            if (get_int(fr, &(set->ss.ext.ijstart), key))
                goto exit;
            if (get_int(fr, &(set->ss.ext.jkstart), key))
                goto exit;
            if (get_int(fr, &(set->ss.ext.extxres), key))
                goto exit;
            if (get_int(fr, &(set->ss.ext.extyres), key))
                goto exit;
            if (get_int(fr, &(set->ss.ext.iextfrom), key))
                goto exit;
            if (get_int(fr, &(set->ss.ext.jextfrom), key))
                goto exit;
            if (get_int(fr, &(set->ss.ext.iextto), key))
                goto exit;
            if (get_int(fr, &(set->ss.ext.jextto), key))
                goto exit;
            if (get_int(fr, &(set->ss.ext.shift), key))
                goto exit;

            fscanf(fr, "%s", buffer); /*no spaces allowed here*/
            set->ss.ext.filebase_ex = g_strdup(buffer);
            fscanf(fr, "%s", buffer); /*no spaces allowed here*/
            set->ss.ext.filebase_ey = g_strdup(buffer);
            fscanf(fr, "%s", buffer); /*no spaces allowed here*/
            set->ss.ext.filebase_ez = g_strdup(buffer);

            fscanf(fr, "%s", buffer); /*no spaces allowed here*/
            set->ss.ext.filebase_hx = g_strdup(buffer);
            fscanf(fr, "%s", buffer); /*no spaces allowed here*/
            set->ss.ext.filebase_hy = g_strdup(buffer);
            fscanf(fr, "%s", buffer); /*no spaces allowed here*/
            set->ss.ext.filebase_hz = g_strdup(buffer);
            set->ss.ext.valid = TRUE;
        } else if (strcmp(key, "SOURCE_LOCAL") == 0) {
            /*alloc or realloc the necessary structure*/
            if (set->ss.nlocals == 0) {
                n_source_locals_size = 10;
                set->ss.locals = (SvSrcLocal *)g_malloc(n_source_locals_size * sizeof(SvSrcLocal));
            } else if (set->ss.nlocals >= n_source_locals_size) { /*realloc*/
                n_source_locals_size += 10;
                set->ss.locals = (SvSrcLocal *)g_realloc(set->ss.locals, n_source_locals_size * sizeof(SvSrcLocal));
            }

            /*load parameters*/
            if (get_int(fr, &(set->ss.locals[set->ss.nlocals].box_i0), key))
                goto exit;
            if (get_int(fr, &(set->ss.locals[set->ss.nlocals].box_j0), key))
                goto exit;
            if (get_int(fr, &(set->ss.locals[set->ss.nlocals].box_k0), key))
                goto exit;
            if (get_int(fr, &(set->ss.locals[set->ss.nlocals].box_in), key))
                goto exit;
            if (get_int(fr, &(set->ss.locals[set->ss.nlocals].box_jn), key))
                goto exit;
            if (get_int(fr, &(set->ss.locals[set->ss.nlocals].box_kn), key))
                goto exit;
            if (get_double(fr, &(set->ss.locals[set->ss.nlocals].layered_epsilon), key))
                goto exit;
            if (get_double(fr, &(set->ss.locals[set->ss.nlocals].density), key))
                goto exit;
            if (get_double(fr, &(set->ss.locals[set->ss.nlocals].strength), key))
                goto exit;
            if (get_double(fr, &(set->ss.locals[set->ss.nlocals].lambda_peak), key))
                goto exit;
            if (get_double(fr, &(set->ss.locals[set->ss.nlocals].lambda_width), key))
                goto exit;
            if (get_int(fr, &(set->ss.locals[set->ss.nlocals].source_mode), key))
                goto exit;
            if (get_int(fr, &(set->ss.locals[set->ss.nlocals].startfrom), key))
                goto exit;

            //set->ss.locals[set->ss.nlocals].layered_epsilon *= EPSILON_0;
            set->ss.nlocals++;
        } else if (strcmp(key, "SOURCE_TSF") == 0) {
            if (get_int(fr, &(set->ss.tsf.box_i0), key))
                goto exit;
            if (get_int(fr, &(set->ss.tsf.box_j0), key))
                goto exit;
            if (get_int(fr, &(set->ss.tsf.box_k0), key))
                goto exit;
            if (get_int(fr, &(set->ss.tsf.box_in), key))
                goto exit;
            if (get_int(fr, &(set->ss.tsf.box_jn), key))
                goto exit;
            if (get_int(fr, &(set->ss.tsf.box_kn), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.ia_theta), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.ia_phi), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.ia_psi), key))
                goto exit;
            if (get_int(fr, &(srctype), key))
                goto exit;
            set->ss.tsf.source_mode = srctype;
            if (srctype == 0) {  //load from file
                fscanf(fr, "%100[^\n\r]", buffer);
                g_free(set->ss.tsf.source_filename);
                set->ss.tsf.source_filename = g_strstrip(g_strdup(buffer));
                if (set->sc.verbose > 1)
                    printf("Source will be loaded from file %s\n", set->ss.tsf.source_filename);
            } else if (srctype == 1 || srctype == 2 || srctype == 3) { //generate simple sources
                if (set->sc.nsteps == 0) {
                    fprintf(stderr, "Error: COMP parameter (number of steps) should precede generated source command.\n");
                    goto exit;
                }
                wlspan = 0;
                if (get_double(fr, &wavelength, key))
                    goto exit;
                if (srctype == 2 && get_int(fr, &pulse_width, key))
                    goto exit;
                if (srctype == 3 && get_double(fr, &wlspan, key))
                    goto exit;
                if (get_double(fr, &amplitude, key))
                    goto exit;

                set->ss.tsf.source_wl = wavelength;
                set->ss.tsf.source_wlspan = wlspan;
                set->ss.tsf.source_pulsewidth = pulse_width;
                set->ss.tsf.source_amplitude = amplitude;

                dt = 1.0 / LIGHT_SPEED / sqrt(1.0 / set->sp.dx / set->sp.dx + 1.0 / set->sp.dy / set->sp.dy + 1.0 / set->sp.dz / set->sp.dz) * set->sc.dtmult;

                if (called_from_gsvit) {
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "tmp_source_tsf_%06d", set->sc.suffix);
                    if (write_source(SRCF_TSF, (SrcType)srctype, set->sc.nsteps, wavelength, wlspan, pulse_width, dt, 0, tmpsourcefile, amplitude, 0, 0))
                        goto exit;

                    g_free(set->ss.tsf.source_filename);
                    set->ss.tsf.source_filename = g_strdup(tmpsourcefile);
                }
            } else {
                fprintf(stderr, "Error: unknown source type\n");
                goto exit;
            }
            set->ss.tsf.valid = TRUE;
        } else if (strcmp(key, "TSF_SKIP") == 0) {
            fscanf(fr, "%s", value);
            if (strcmp(value, "i0") == 0)
                set->ss.tsf.box_boundary_skipi0 = 1;
            else if (strcmp(value, "in") == 0)
                set->ss.tsf.box_boundary_skipin = 1;
            else if (strcmp(value, "j0") == 0)
                set->ss.tsf.box_boundary_skipj0 = 1;
            else if (strcmp(value, "jn") == 0)
                set->ss.tsf.box_boundary_skipjn = 1;
            else if (strcmp(value, "k0") == 0)
                set->ss.tsf.box_boundary_skipk0 = 1;
            else if (strcmp(value, "kn") == 0)
                set->ss.tsf.box_boundary_skipkn = 1;
            else if (strcmp(value, "depth") == 0) {
                if (get_int(fr, &(set->ss.tsf.box_boundary_skipdepth), key))
                    goto exit;
            }
        } else if (strcmp(key, "TSF_GAUSSIAN_Z") == 0) {
            set->ss.tsf.gaussian = 1;
            if (get_double(fr, &(set->ss.tsf.gaussian_fxpos), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.gaussian_fypos), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.gaussian_rx), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.gaussian_ry), key))
                goto exit;
        } else if (strcmp(key, "TSF_RADIAL_Z") == 0) {
            set->ss.tsf.radial = 1;
            if (get_double(fr, &(set->ss.tsf.radial_fxpos), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.radial_fypos), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.radial_rx), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.radial_ry), key))
                goto exit;
        } else if (strcmp(key, "TSF_FIBER_Z") == 0) {
            set->ss.tsf.fiber = 1;
            if (get_double(fr, &(set->ss.tsf.fiber_fxpos), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.fiber_fypos), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.fiber_radius), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.fiber_cutoff), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.fiber_epsilon_core), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.fiber_epsilon_cladding), key))
                goto exit;
        } else if (strcmp(key, "SOURCE_LTSF") == 0) {
            if (get_int(fr, &(set->ss.ltsf.box_i0), key))
                goto exit;
            if (get_int(fr, &(set->ss.ltsf.box_j0), key))
                goto exit;
            if (get_int(fr, &(set->ss.ltsf.box_k0), key))
                goto exit;
            if (get_int(fr, &(set->ss.ltsf.box_in), key))
                goto exit;
            if (get_int(fr, &(set->ss.ltsf.box_jn), key))
                goto exit;
            if (get_int(fr, &(set->ss.ltsf.box_kn), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsf.ia_theta), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsf.ia_phi), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsf.ia_psi), key))
                goto exit;
            if (get_int(fr, &(set->ss.ltsf.layered_count), key))
                goto exit;
            for (i = 0; i < set->ss.ltsf.layered_count; i++) {
                gint mat99 = 0;
                long int pos = 0;
                if (get_int(fr, &(set->ss.ltsf.layered_zpos[i]), key))
                    goto exit;

                pos = ftell(fr);
                if (get_int(fr, &mat99, key))
                    goto exit;
                if (mat99 == 99) {
                    fscanf(fr, "%s", buffer);
                    g_free(set->ss.ltsf.layered_material[i]);
                    set->ss.ltsf.layered_material[i] = g_strdup(buffer);
                } else {
                    fseek(fr, pos, SEEK_SET);

                    if (get_double(fr, &(set->ss.ltsf.layered_epsilon[i]), key))
                        goto exit;
                    if (get_double(fr, &(set->ss.ltsf.layered_mu[i]), key))
                        goto exit;
                    if (get_double(fr, &(set->ss.ltsf.layered_sigma[i]), key))
                        goto exit;
                    if (get_double(fr, &(set->ss.ltsf.layered_sigast[i]), key))
                        goto exit;
                }                
            }

            if (get_int(fr, &(srctype), key))
                goto exit;
            set->ss.ltsf.source_mode = srctype;

            if (srctype == 0) {  //load from file
                fscanf(fr, "%100[^\n\r]", buffer);
                set->ss.ltsf.source_filename = g_strstrip(g_strdup(buffer));

                if (set->sc.verbose > 1)
                    printf("Source will be loaded from file %s\n", set->ss.ltsf.source_filename);
            } else if (srctype == 1 || srctype == 2 || srctype == 3) { //generate simple sources
                if (set->sc.nsteps == 0) {
                    fprintf(stderr, "Error: COMP parameter (number of steps) should precede generated source command.\n");
                    goto exit;
                }
                wlspan = 0;
                if (get_double(fr, &wavelength, key))
                    goto exit;
                if (srctype == 2 && get_int(fr, &pulse_width, key))
                    goto exit;
                if (srctype == 3 && get_double(fr, &wlspan, key))
                    goto exit;
                if (get_double(fr, &amplitude, key))
                    goto exit;

                set->ss.ltsf.source_wl = wavelength;
                set->ss.ltsf.source_wlspan = wlspan;
                set->ss.ltsf.source_pulsewidth = pulse_width;
                set->ss.ltsf.source_amplitude = amplitude;

                dt = 1.0 / LIGHT_SPEED / sqrt(1.0 / set->sp.dx / set->sp.dx + 1.0 / set->sp.dy / set->sp.dy + 1.0 / set->sp.dz / set->sp.dz) * set->sc.dtmult;

                if (called_from_gsvit) {
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "tmp_source_ltsf_%06d", set->sc.suffix);
                    if (write_source(SRCF_LTSF, (SrcType)srctype, set->sc.nsteps, wavelength, wlspan, pulse_width, dt, 0, tmpsourcefile, amplitude, 0, 0))
                        goto exit;

                    g_free(set->ss.ltsf.source_filename);
                    set->ss.ltsf.source_filename = g_strdup(tmpsourcefile);
                }
            } else {
                fprintf(stderr, "Error: unknown source type\n");
                goto exit;
            }
            set->ss.ltsf.valid = TRUE;
        } else if (strcmp(key, "LTSF_GAUSSIAN") == 0) {
            set->ss.ltsf.gaussian = 1;
            if (get_double(fr, &(set->ss.ltsf.gaussian_fxpos), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsf.gaussian_fypos), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsf.gaussian_rx), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsf.gaussian_ry), key))
                goto exit;
        } else if (strcmp(key, "LTSF_RADIAL") == 0) {
            set->ss.ltsf.radial = 1;
            if (get_double(fr, &(set->ss.ltsf.radial_fxpos), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsf.radial_fypos), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsf.radial_rx), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsf.radial_ry), key))
                goto exit;
        } else if (strcmp(key, "LTSF_FIBER") == 0) {
            set->ss.ltsf.fiber = 1;
            if (get_double(fr, &(set->ss.ltsf.fiber_fxpos), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsf.fiber_fypos), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsf.fiber_radius), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsf.fiber_cutoff), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsf.fiber_epsilon_core), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsf.fiber_epsilon_cladding), key))
                goto exit;
        } else if (strcmp(key, "LTSF_SKIP") == 0) {
            fscanf(fr, "%s", value);
            if (strcmp(value, "i0") == 0)
                set->ss.ltsf.box_boundary_skipi0 = 1;
            else if (strcmp(value, "in") == 0)
                set->ss.ltsf.box_boundary_skipin = 1;
            else if (strcmp(value, "j0") == 0)
                set->ss.ltsf.box_boundary_skipj0 = 1;
            else if (strcmp(value, "jn") == 0)
                set->ss.ltsf.box_boundary_skipjn = 1;
            else if (strcmp(value, "k0") == 0)
                set->ss.ltsf.box_boundary_skipk0 = 1;
            else if (strcmp(value, "kn") == 0)
                set->ss.ltsf.box_boundary_skipkn = 1;
            else if (strcmp(value, "depth") == 0) {
                if (get_int(fr, &(set->ss.ltsf.box_boundary_skipdepth), key))
                    goto exit;
            }
        } else if (strcmp(key, "SOURCE_LTSFF") == 0) {
            if (get_int(fr, &(set->ss.ltsff.box_i0), key))
                goto exit;
            if (get_int(fr, &(set->ss.ltsff.box_j0), key))
                goto exit;
            if (get_int(fr, &(set->ss.ltsff.box_k0), key))
                goto exit;
            if (get_int(fr, &(set->ss.ltsff.box_in), key))
                goto exit;
            if (get_int(fr, &(set->ss.ltsff.box_jn), key))
                goto exit;
            if (get_int(fr, &(set->ss.ltsff.box_kn), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsff.focused_thetamax), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsff.focused_fdist), key))
                goto exit;
            if (get_double(fr, &(set->ss.ltsff.focused_pol), key))
                goto exit;
            if (get_int(fr, &(set->ss.ltsff.focused_nip), key))
                goto exit;
            if (get_int(fr, &(set->ss.ltsff.focused_mip), key))
                goto exit;
            if (get_int(fr, &(set->ss.ltsff.layered_count), key))
                goto exit;
            for (i = 0; i < set->ss.ltsff.layered_count; i++) {
                gint mat99 = 0;
                long int pos = 0;
                if (get_int(fr, &(set->ss.ltsff.layered_zpos[i]), key))
                    goto exit;

                pos = ftell(fr);
                if (get_int(fr, &mat99, key))
                    goto exit;
                if (mat99 == 99) {
                    fscanf(fr, "%s", buffer);
                    g_free(set->ss.ltsff.layered_material[i]);
                    set->ss.ltsff.layered_material[i] = g_strdup(buffer);
                }
                else {
                    fseek(fr, pos, SEEK_SET);

                    if (get_double(fr, &(set->ss.ltsff.layered_epsilon[i]), key))
                        goto exit;
                    if (get_double(fr, &(set->ss.ltsff.layered_mu[i]), key))
                        goto exit;
                    if (get_double(fr, &(set->ss.ltsff.layered_sigma[i]), key))
                        goto exit;
                    if (get_double(fr, &(set->ss.ltsff.layered_sigast[i]), key))
                        goto exit;
                }
            }

            if (get_int(fr, &(srctype), key))
                goto exit;
            set->ss.ltsff.source_mode = srctype;

            if (srctype == 0) {  //load from file
                fscanf(fr, "%100[^\n\r]", buffer);
                set->ss.ltsff.source_filename = g_strstrip(g_strdup(buffer));
                if (set->sc.verbose > 1)
                    printf("Source will be loaded from file %s\n", set->ss.ltsff.source_filename);
            } else if (srctype == 1 || srctype == 2 || srctype == 3) { //generate simple sources
                if (set->sc.nsteps == 0) {
                    fprintf(stderr, "Error: COMP parameter (number of steps) should precede generated source command.\n");
                    goto exit;
                }
                wlspan = 0;
                if (get_double(fr, &wavelength, key))
                    goto exit;
                if (srctype == 2 && get_int(fr, &pulse_width, key))
                    goto exit;
                if (srctype == 3 && get_double(fr, &wlspan, key))
                    goto exit;
                if (get_double(fr, &amplitude, key))
                    goto exit;

                set->ss.ltsff.source_wl = wavelength;
                set->ss.ltsff.source_wlspan = wlspan;
                set->ss.ltsff.source_pulsewidth = pulse_width;
                set->ss.ltsff.source_amplitude = amplitude;

                dt = 1.0 / LIGHT_SPEED / sqrt(1.0 / set->sp.dx / set->sp.dx + 1.0 / set->sp.dy / set->sp.dy + 1.0 / set->sp.dz / set->sp.dz) * set->sc.dtmult;

                if (called_from_gsvit) {
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "tmp_source_ltsff_%06d", set->sc.suffix);
                    if (write_source(SRCF_LTSF, (SrcType)srctype, set->sc.nsteps, wavelength, wlspan, pulse_width, dt, 0, tmpsourcefile, amplitude, 0, 0))
                        goto exit;

                    g_free(set->ss.ltsff.source_filename);
                    set->ss.ltsff.source_filename = g_strdup(tmpsourcefile);
                }
            } else {
                fprintf(stderr, "Error: unknown source type\n");
                goto exit;
            }
            set->ss.ltsff.valid = TRUE;
        } else if (strcmp(key, "LTSFF_FAST") == 0) {
            if (get_int(fr, &(set->ss.ltsff.fast), key))
                goto exit;
        } else if (strcmp(key, "LTSFF_SKIP") == 0) {
            fscanf(fr, "%s", value);
            if (strcmp(value, "i0") == 0)
                set->ss.ltsff.box_boundary_skipi0 = 1;
            else if (strcmp(value, "in") == 0)
                set->ss.ltsff.box_boundary_skipin = 1;
            else if (strcmp(value, "j0") == 0)
                set->ss.ltsff.box_boundary_skipj0 = 1;
            else if (strcmp(value, "jn") == 0)
                set->ss.ltsff.box_boundary_skipjn = 1;
            else if (strcmp(value, "k0") == 0)
                set->ss.ltsff.box_boundary_skipk0 = 1;
            else if (strcmp(value, "kn") == 0)
                set->ss.ltsff.box_boundary_skipkn = 1;
            else if (strcmp(value, "depth") == 0) {
                if (get_int(fr, &(set->ss.ltsff.box_boundary_skipdepth), key))
                    goto exit;
            }
        } else if (strcmp(key, "SOURCE_TSFF") == 0) {
            if (get_int(fr, &(set->ss.tsff.box_i0), key))
                goto exit;
            if (get_int(fr, &(set->ss.tsff.box_j0), key))
                goto exit;
            if (get_int(fr, &(set->ss.tsff.box_k0), key))
                goto exit;
            if (get_int(fr, &(set->ss.tsff.box_in), key))
                goto exit;
            if (get_int(fr, &(set->ss.tsff.box_jn), key))
                goto exit;
            if (get_int(fr, &(set->ss.tsff.box_kn), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsff.focused_thetamax), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsff.focused_fdist), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsff.focused_pol), key))
                goto exit;
            if (get_int(fr, &(set->ss.tsff.focused_nip), key))
                goto exit;
            if (get_int(fr, &(set->ss.tsff.focused_mip), key))
                goto exit;
            if (get_int(fr, &(srctype), key))
                goto exit;
            set->ss.tsff.source_mode = srctype;

            if (srctype == 0) {  //load from file
                fscanf(fr, "%100[^\n\r]", buffer);
                set->ss.tsff.source_filename = g_strstrip(g_strdup(buffer));

                if (set->sc.verbose > 1)
                    printf("Source will be loaded from file %s\n", set->ss.tsff.source_filename);
            } else if (srctype == 1 || srctype == 2 || srctype == 3) { //generate simple sources
                if (set->sc.nsteps == 0) {
                    fprintf(stderr, "Error: COMP parameter (number of steps) should precede generated source command.\n");
                    goto exit;
                }
                wlspan = 0;
                if (get_double(fr, &wavelength, key))
                    goto exit;
                if (srctype == 2 && get_int(fr, &pulse_width, key))
                    goto exit;
                if (srctype == 3 && get_double(fr, &wlspan, key))
                    goto exit;
                if (get_double(fr, &amplitude, key))
                    goto exit;

                set->ss.tsff.source_wl = wavelength;
                set->ss.tsff.source_wlspan = wlspan;
                set->ss.tsff.source_pulsewidth = pulse_width;
                set->ss.tsff.source_amplitude = amplitude;

                dt = 1.0 / LIGHT_SPEED / sqrt(1.0 / set->sp.dx / set->sp.dx + 1.0 / set->sp.dy / set->sp.dy + 1.0 / set->sp.dz / set->sp.dz) * set->sc.dtmult;

                if (called_from_gsvit) {
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "tmp_source_tsff_%06d", set->sc.suffix);
                    if (write_source(SRCF_TSF, (SrcType)srctype, set->sc.nsteps, wavelength, wlspan, pulse_width, dt, 0, tmpsourcefile, amplitude, 0, 0))
                        goto exit;

                    g_free(set->ss.tsff.source_filename);
                    set->ss.tsff.source_filename = g_strdup(tmpsourcefile);
                }
            } else {
                fprintf(stderr, "Error: unknown source type\n");
                goto exit;
            }
            set->ss.tsff.valid = TRUE;
        } else if (strcmp(key, "TSFF_SHIFT") == 0) {
            if (get_double(fr, &(set->ss.tsff.xshift), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsff.yshift), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsff.zshift), key))
                goto exit;
        } else if (strcmp(key, "TSFF_SKIP") == 0) {
            fscanf(fr, "%s", value);
            if (strcmp(value, "i0") == 0)
                set->ss.tsff.box_boundary_skipi0 = 1;
            else if (strcmp(value, "in") == 0)
                set->ss.tsff.box_boundary_skipin = 1;
            else if (strcmp(value, "j0") == 0)
                set->ss.tsff.box_boundary_skipj0 = 1;
            else if (strcmp(value, "jn") == 0)
                set->ss.tsff.box_boundary_skipjn = 1;
            else if (strcmp(value, "k0") == 0)
                set->ss.tsff.box_boundary_skipk0 = 1;
            else if (strcmp(value, "kn") == 0)
                set->ss.tsff.box_boundary_skipkn = 1;
            else if (strcmp(value, "depth") == 0) {
                if (get_int(fr, &(set->ss.tsff.box_boundary_skipdepth), key))
                    goto exit;
            }
        } else if (strcmp(key, "SOURCE_SF") == 0) {
            if (get_double(fr, &(set->ss.sf.ia_theta), key))
                goto exit;
            if (get_double(fr, &(set->ss.sf.ia_phi), key))
                goto exit;
            if (get_double(fr, &(set->ss.sf.ia_psi), key))
                goto exit;
            if (get_int(fr, &(srctype), key))
                goto exit;
            set->ss.sf.source_mode = srctype;

            if (srctype == 0) {  //load from file
                fscanf(fr, "%100[^\n\r]", buffer);
                set->ss.sf.source_filename = g_strstrip(g_strdup(buffer));

                if (set->sc.verbose > 1)
                    printf("Source will be loaded from file %s\n", set->ss.sf.source_filename);
            } else if (srctype == 1 || srctype == 2 || srctype == 3) { //generate simple sources
                if (set->sc.nsteps == 0) {
                    fprintf(stderr, "Error: COMP parameter (number of steps) should precede generated source command.\n");
                    goto exit;
                }
                wlspan = 0;
                if (get_double(fr, &wavelength, key))
                    goto exit;
                if (srctype == 2 && get_int(fr, &pulse_width, key))
                    goto exit;
                if (srctype == 3 && get_double(fr, &wlspan, key))
                    goto exit;
                if (get_double(fr, &amplitude, key))
                    goto exit;

                set->ss.sf.source_wl = wavelength;
                set->ss.sf.source_wlspan = wlspan;
                set->ss.sf.source_pulsewidth = pulse_width;
                set->ss.sf.source_amplitude = amplitude;

                dt = 1.0 / LIGHT_SPEED / sqrt(1.0 / set->sp.dx / set->sp.dx + 1.0 / set->sp.dy / set->sp.dy + 1.0 / set->sp.dz / set->sp.dz) * set->sc.dtmult;

                if (called_from_gsvit) {
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "tmp_source_sf_%06d", set->sc.suffix);

#ifdef G_OS_WIN32
                    /*gchar* filename_tmp = NULL;
                    filename_tmp = g_build_filename(g_get_home_dir(), ".xsvit", tmpsourcefile, NULL);
                    g_snprintf(tmpsourcefile, sizeof(tmpsourcefile), "%s", filename_tmp);
                    g_free(filename_tmp);*/
#endif

                    if (write_source(SRCF_TSF, (SrcType)srctype, set->sc.nsteps, wavelength, wlspan, pulse_width, dt, 0, tmpsourcefile, amplitude, 0, 0))
                        goto exit;

                    g_free(set->ss.sf.source_filename);
                    set->ss.sf.source_filename = g_strdup(tmpsourcefile);
                }

                //if (write_source(SRCF_TSF, (SrcType)srctype, set->sc.nsteps, wavelength, wlspan, size, dt, 0, "tmpsource", amplitude, 0, 0)) 
                //    goto exit;                  
                //set->ss.sf.source_filename = g_strdup("tmpsource");                
            } else {
                fprintf(stderr, "Error: unknown source type\n");
                goto exit;
            }
            set->ss.sf.valid = TRUE;
        } else if (strcmp(key, "TSFSOURCE_MATERIAL") == 0) {
            if (get_double(fr, &(set->ss.tsf.layered_epsilon), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.layered_mu), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.layered_sigma), key))
                goto exit;
            if (get_double(fr, &(set->ss.tsf.layered_sigast), key))
                goto exit;
        } else if (strcmp(key, "SOURCE_WAVELENGTH") == 0) {
            if (get_double(fr, &(set->ss.lambda_min), key))
                goto exit;
            if (get_double(fr, &(set->ss.lambda_center), key))
                goto exit;
            if (get_double(fr, &(set->ss.lambda_max), key))
                goto exit;
        } else if (strcmp(key, "OUT_POINT") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_POINT);
            //    	    if (set->so.npnts == 0) {
            //    		    n_output_pnts_size = 10;
            //    		    set->so.pnts = (SvOutputPar *)g_malloc(n_output_pnts_size*sizeof(SvOutputPar));
            //    	    } else if (set->so.npnts >= n_output_pnts_size) { /*realloc*/
            //    		    n_output_pnts_size += 10;
            
            //    		    set->so.pnts = (SvOutputPar *)g_realloc(set->so.pnts, n_output_pnts_size*sizeof(SvOutputPar));
            //    	    }


                        /*load parameters*/
            if (get_component(fr, &(set->so.pnts[set->so.npnts].component), key))
                goto exit;
            if (get_int(fr, &(set->so.pnts[set->so.npnts].step), key))
                goto exit;
            if (get_int(fr, &(set->so.pnts[set->so.npnts].i), key))
                goto exit;
            if (get_int(fr, &(set->so.pnts[set->so.npnts].j), key))
                goto exit;
            if (get_int(fr, &(set->so.pnts[set->so.npnts].k), key))
                goto exit;

            if (set->so.pnts[set->so.npnts].component == SV_COMP_ALLFFT) {
                if (get_int(fr, &(set->so.pnts[set->so.npnts].start), key))
                    goto exit;
                if (get_int(fr, &(set->so.pnts[set->so.npnts].stop), key))
                    goto exit;

                if (set->so.pnts[set->so.npnts].start == -1)
                    set->so.pnts[set->so.npnts].start = 0;
                if (set->so.pnts[set->so.npnts].stop == -1)
                    set->so.pnts[set->so.npnts].stop = set->sc.nsteps;

                if (set->so.pnts[set->so.npnts].start > (set->sc.nsteps - 4)) {
                    set->so.pnts[set->so.npnts].start = set->sc.nsteps - 4;
                    fprintf(stderr, "Warning: requesting FFT output outside of the data range, cropped to %d\n", set->so.pnts[set->so.npnts].start);
                }

                if (set->so.pnts[set->so.npnts].stop > (set->sc.nsteps)) {
                    set->so.pnts[set->so.npnts].stop = set->sc.nsteps;
                    fprintf(stderr, "Warning: requesting FFT output outside of the data range, cropped to %d\n", set->so.pnts[set->so.npnts].stop);
                }
                set->so.savespectrum = 1;
            }

            fscanf(fr, "%100[^\n\r]", buffer);
            sprintf(set->so.pnts[set->so.npnts].filebase, "%s", g_strstrip(g_strdup(buffer)));
            set->so.npnts++;
        } else if (strcmp(key, "OUT_IMAGE") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_IMAGE);
            /*load parameters*/
            if (get_component(fr, &(set->so.imgs[set->so.nimgs].component), key))
                goto exit;
            if (get_int(fr, &(set->so.imgs[set->so.nimgs].step), key))
                goto exit;
            if (get_int(fr, &(set->so.imgs[set->so.nimgs].i), key))
                goto exit;
            if (get_int(fr, &(set->so.imgs[set->so.nimgs].j), key))
                goto exit;
            if (get_int(fr, &(set->so.imgs[set->so.nimgs].k), key))
                goto exit;

            fscanf(fr, "%100[^\n\r]", buffer);
            sprintf(set->so.imgs[set->so.nimgs].filebase, "%s", g_strstrip(g_strdup(buffer)));

            set->so.nimgs++;
        }   else if (strcmp(key, "OUT_SUBGRID_IMAGE") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_SUBGRID_IMAGE);
            /*load parameters*/
            if (get_component(fr, &(set->so.subgrid_imgs[set->so.nsubgrid_imgs].component), key))
                goto exit;
            if (get_int(fr, &(set->so.subgrid_imgs[set->so.nsubgrid_imgs].step), key))
                goto exit;
            if (get_int(fr, &(set->so.subgrid_imgs[set->so.nsubgrid_imgs].i), key))
                goto exit;
            if (get_int(fr, &(set->so.subgrid_imgs[set->so.nsubgrid_imgs].j), key))
                goto exit;
            if (get_int(fr, &(set->so.subgrid_imgs[set->so.nsubgrid_imgs].k), key))
                goto exit;
            if (get_int(fr, &(set->so.subgrid_imgs[set->so.nsubgrid_imgs].format), key)) //this is the numebr of subgrid
                goto exit;
 
            fscanf(fr, "%100[^\n\r]", buffer);
            sprintf(set->so.subgrid_imgs[set->so.nsubgrid_imgs].filebase, "%s", g_strstrip(g_strdup(buffer)));

            set->so.nsubgrid_imgs++;
        } else if (strcmp(key, "OUT_PLANE") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_PLANE);
            //    	    if (set->so.nplns == 0) {
            //    		    n_output_plns_size = 10;
            //    		    set->so.plns = (SvOutputPar *)g_malloc(n_output_plns_size*sizeof(SvOutputPar));
            //    	    } else if (set->so.nplns >= n_output_plns_size) { /*realloc*/
            //    		    n_output_plns_size += 10;
            //    		    set->so.plns = (SvOutputPar *)g_realloc(set->so.plns, n_output_plns_size*sizeof(SvOutputPar));
            //    	    }

                        /*load parameters*/
            if (get_component(fr, &(set->so.plns[set->so.nplns].component), key))
                goto exit;
            if (get_int(fr, &(set->so.plns[set->so.nplns].step), key))
                goto exit;
            if (get_int(fr, &(set->so.plns[set->so.nplns].start), key))
                goto exit;
            if (get_int(fr, &(set->so.plns[set->so.nplns].stop), key))
                goto exit;
            if (get_int(fr, &(set->so.plns[set->so.nplns].format), key))
                goto exit; //here it means binary/ascii
            if (get_int(fr, &(set->so.plns[set->so.nplns].i), key))
                goto exit;
            if (get_int(fr, &(set->so.plns[set->so.nplns].j), key))
                goto exit;
            if (get_int(fr, &(set->so.plns[set->so.nplns].k), key))
                goto exit;

            fscanf(fr, "%100[^\n\r]", buffer);
            sprintf(set->so.plns[set->so.nplns].filebase, "%s", g_strstrip(g_strdup(buffer)));

            set->so.nplns++;
        } else if (strcmp(key, "OUT_VOLUME") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_VOLUME);
            //    	    if (set->so.ncubs == 0) {
            //    		    n_output_cubs_size = 10;
            //    		    set->so.cubs = (SvOutputPar *)g_malloc(n_output_cubs_size*sizeof(SvOutputPar));
            //    	    } else if (set->so.ncubs >= n_output_cubs_size) { /*realloc*/
            //    		    n_output_cubs_size += 10;
            //    		    set->so.cubs = (SvOutputPar *)g_realloc(set->so.cubs, n_output_cubs_size*sizeof(SvOutputPar));
            //    	    }

                        /*load parameters*/
            if (get_cube_component(fr, &(set->so.cubs[set->so.ncubs].component), key))
                goto exit;
            if (get_int(fr, &(set->so.cubs[set->so.ncubs].step), key))
                goto exit;
            if (get_int(fr, &(set->so.cubs[set->so.ncubs].start), key))
                goto exit;
            if (get_int(fr, &(set->so.cubs[set->so.ncubs].stop), key))
                goto exit;
            if (get_int(fr, &(set->so.cubs[set->so.ncubs].format), key))
                goto exit; //here it means binary/ascii

            if (set->so.cubs[set->so.ncubs].start == -1)
                set->so.cubs[set->so.ncubs].start = 0;
            if (set->so.cubs[set->so.ncubs].stop == -1)
                set->so.cubs[set->so.ncubs].stop = set->sc.nsteps;

            fscanf(fr, "%100[^\n\r]", buffer);
            sprintf(set->so.cubs[set->so.ncubs].filebase, "%s", g_strstrip(g_strdup(buffer)));

            if (set->so.cubs[set->so.ncubs].component == SV_OVOLUME_ABS || set->so.cubs[set->so.ncubs].component == SV_OVOLUME_SUMALL || set->so.cubs[set->so.ncubs].component == SV_OVOLUME_MAXALL)
                set->so.cubs[set->so.ncubs].k = n_abs++; //position of data 

            set->so.ncubs++;
        } else if (strcmp(key, "OUT_SUM") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_SUM);
            //    	    if (set->so.nsums == 0) {
            //    		    n_output_sums_size = 10;
            //    		    set->so.sums = (SvOutputSum *)g_malloc(n_output_sums_size*sizeof(SvOutputSum));
            //    	    } else if (set->so.nsums >= n_output_sums_size) { /*realloc*/
            //    		    n_output_sums_size += 10;
            //    		    set->so.sums = (SvOutputSum *)g_realloc(set->so.sums, n_output_sums_size*sizeof(SvOutputSum));
            //    	    }

                        /*load parameters*/
            if (get_scomponent(fr, &(set->so.sums[set->so.nsums].component), key))
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].step), key))
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].box_i0), key))
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].box_j0), key))
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].box_k0), key))
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].box_in), key))
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].box_jn), key))
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].box_kn), key))
                goto exit;
            if (get_double(fr, &(set->so.sums[set->so.nsums].layered_epsilon), key))
                goto exit;
            if (get_double(fr, &(set->so.sums[set->so.nsums].layered_mu), key))
                goto exit;
            if (get_double(fr, &(set->so.sums[set->so.nsums].layered_sigma), key))
                goto exit;
            if (get_double(fr, &(set->so.sums[set->so.nsums].layered_sigast), key))
                goto exit;
            fscanf(fr, "%100[^\n\r]", buffer);
            sprintf(set->so.sums[set->so.nsums].filename, "%s", g_strstrip(g_strdup(buffer)));

            /*if (set->so.sums[set->so.nsums].layered_epsilon != -1)
                set->so.sums[set->so.nsums].layered_epsilon *= EPSILON_0;
            set->so.sums[set->so.nsums].layered_mu *= MU_0;*/

            set->so.sums[set->so.nsums].stringbased = 0;
            set->so.nsums++;
        } else if (strcmp(key, "OUT_SUMTAB") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_SUMTAB);
            //    	    if (set->so.nsums == 0) {
            //    		    n_output_sums_size = 10;
            //    		    set->so.sums = (SvOutputSum *)g_malloc(n_output_sums_size*sizeof(SvOutputSum));
            //    	    } else if (set->so.nsums >= n_output_sums_size) { /*realloc*/
            //    		    n_output_sums_size += 10;
            //    		    set->so.sums = (SvOutputSum *)g_realloc(set->so.sums, n_output_sums_size*sizeof(SvOutputSum));
            //    	    }

                        /*load parameters*/
            if (get_scomponent(fr, &(set->so.sums[set->so.nsums].component), key))
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].step), key))
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].box_i0), key))
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].box_j0), key))
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].box_k0), key))
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].box_in), key))
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].box_jn), key))
                goto exit;
            if (get_int(fr, &(set->so.sums[set->so.nsums].box_kn), key))
                goto exit;
            fscanf(fr, "%s", buffer);
            sprintf(set->so.sums[set->so.nsums].string, "%s", buffer);

            fscanf(fr, "%100[^\n\r]", buffer);
            sprintf(set->so.sums[set->so.nsums].filename, "%s", g_strstrip(g_strdup(buffer)));

            set->so.sums[set->so.nsums].stringbased = 1;
            set->so.nsums++;
        } else if (strcmp(key, "OUT_FORCE") == 0) {
            /*alloc or realloc the necessary structure*/
            output_alloc(set, SV_OUTTYPE_FORCE);
            //if (set->so.nforces == 0) {
            //    n_output_forces_size = 10;
            //    set->so.forces = (SvOutputForce *)g_malloc(n_output_forces_size*sizeof(SvOutputForce));
            //} else if (set->so.nforces >= n_output_forces_size) { /*realloc*/
            //    n_output_forces_size += 10;
            //    set->so.forces = (SvOutputForce *)g_realloc(set->so.forces, n_output_forces_size*sizeof(SvOutputForce));
            //}

            /*load parameters*/
            if (get_int(fr, &(set->so.forces[set->so.nforces].step), key))
                goto exit;
            if (get_int(fr, &(set->so.forces[set->so.nforces].box_i0), key))
                goto exit;
            if (get_int(fr, &(set->so.forces[set->so.nforces].box_j0), key))
                goto exit;
            if (get_int(fr, &(set->so.forces[set->so.nforces].box_k0), key))
                goto exit;
            if (get_int(fr, &(set->so.forces[set->so.nforces].box_in), key))
                goto exit;
            if (get_int(fr, &(set->so.forces[set->so.nforces].box_jn), key))
                goto exit;
            if (get_int(fr, &(set->so.forces[set->so.nforces].box_kn), key))
                goto exit;
            fscanf(fr, "%100[^\n\r]", buffer);
            sprintf(set->so.forces[set->so.nforces].filename, "%s", g_strstrip(g_strdup(buffer)));
            set->so.nforces++;
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
        } else if (strcmp(key, "MBOUNDARY_XN") == 0) {
            fscanf(fr, "%s", value);
            if (get_int(fr, &(set->smb.bxnpos), key)) {
                goto exit;
            }
            if (strcmp(value, "periodic") == 0) {
                set->smb.bxn = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "MBOUNDARY_Y0") == 0) {
            fscanf(fr, "%s", value);
            if (get_int(fr, &(set->smb.by0pos), key)) {
                goto exit;
            }
            if (strcmp(value, "periodic") == 0) {
                set->smb.by0 = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "MBOUNDARY_YN") == 0) {
            fscanf(fr, "%s", value);
            if (get_int(fr, &(set->smb.bynpos), key))
                goto exit;
            if (strcmp(value, "periodic") == 0) {
                set->smb.byn = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "MBOUNDARY_Z0") == 0) {
            fscanf(fr, "%s", value);
            if (get_int(fr, &(set->smb.bz0pos), key))
                goto exit;
            if (strcmp(value, "periodic") == 0) {
                set->smb.bz0 = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "MBOUNDARY_ZN") == 0) {
            fscanf(fr, "%s", value);
            if (get_int(fr, &(set->smb.bznpos), key))
                goto exit;
            if (strcmp(value, "periodic") == 0) {
                set->smb.bzn = SV_BOUNDARY_PERIODIC;
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
            if (get_int(fr, &(set->smb.bz0pos), key))
                goto exit;
            if (get_int(fr, &(set->smb.bznpos), key))
                goto exit;
            if (strcmp(value, "periodic") == 0) {
                set->smb.bx0 = set->smb.bxn = set->smb.by0 = set->smb.byn = set->smb.bz0 = set->smb.bzn = SV_BOUNDARY_PERIODIC;
            } else {
                fprintf(stderr, "Error: unsupported medium boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "BOUNDARY_ALL") == 0) {
            fscanf(fr, "%s", value);
            if (strcmp(value, "none") == 0) {
                set->sb.bx0 = set->sb.bxn = set->sb.by0 = set->sb.byn = set->sb.bz0 = set->sb.bzn = SV_BOUNDARY_NONE;
            } else if (strcmp(value, "pec") == 0) {
                set->sb.bx0 = set->sb.bxn = set->sb.by0 = set->sb.byn = set->sb.bz0 = set->sb.bzn = SV_BOUNDARY_PEC;
            } else if (strcmp(value, "liao") == 0) {
                set->sb.bx0 = set->sb.bxn = set->sb.by0 = set->sb.byn = set->sb.bz0 = set->sb.bzn = SV_BOUNDARY_LIAO;
            } else if (strcmp(value, "cpml") == 0) {
                set->sb.bx0 = set->sb.bxn = set->sb.by0 = set->sb.byn = set->sb.bz0 = set->sb.bzn = SV_BOUNDARY_CPML;

                if (get_int(fr, &(set->sb.depth_bx0), key))
                    goto exit;
                set->sb.depth_bxn = set->sb.depth_by0 = set->sb.depth_byn = set->sb.depth_bz0 = set->sb.depth_bzn = set->sb.depth_bx0;

                if (get_int(fr, &(set->sb.m_bx0), key))
                    goto exit;
                set->sb.m_bxn = set->sb.m_by0 = set->sb.m_byn = set->sb.m_bz0 = set->sb.m_bzn = set->sb.m_bx0;

                if (get_double(fr, &(set->sb.sigma_bx0), key))
                    goto exit;
                if ((TRUE == called_from_gsvit) && (set->sb.sigma_bx0 == -1)) {
                    set->sb.sigma_bx0 = 0.8 * (gdouble)set->sb.m_bx0 / (377 * set->sp.dx);
                    if (set->sc.verbose > 1)
                        printf("CPML sigma set to optimum value of %g\n", set->sb.sigma_bx0);
                }
                set->sb.sigma_bxn = set->sb.sigma_by0 = set->sb.sigma_byn = set->sb.sigma_bz0 = set->sb.sigma_bzn = set->sb.sigma_bx0;

                if (get_double(fr, &(set->sb.a_bx0), key))
                    goto exit;
                set->sb.a_bxn = set->sb.a_by0 = set->sb.a_byn = set->sb.a_bz0 = set->sb.a_bzn = set->sb.a_bx0;

                if (get_double(fr, &(set->sb.kappa_bx0), key))
                    goto exit;
                set->sb.kappa_bxn = set->sb.kappa_by0 = set->sb.kappa_byn = set->sb.kappa_bz0 = set->sb.kappa_bzn = set->sb.kappa_bx0;
            } else {
                fprintf(stderr, "Error: unsupported boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "BOUNDARY_X0") == 0) {
            fscanf(fr, "%s", value);
            if (strcmp(value, "none") == 0)
                set->sb.bx0 = SV_BOUNDARY_NONE;
            else if (strcmp(value, "pec") == 0)
                set->sb.bx0 = SV_BOUNDARY_PEC;
            else if (strcmp(value, "liao") == 0)
                set->sb.bx0 = SV_BOUNDARY_LIAO;
            else if (strcmp(value, "cpml") == 0) {
                set->sb.bx0 = SV_BOUNDARY_CPML;

                if (get_int(fr, &(set->sb.depth_bx0), key))
                    goto exit;
                if (get_int(fr, &(set->sb.m_bx0), key))
                    goto exit;
                if (get_double(fr, &(set->sb.sigma_bx0), key))
                    goto exit;
                if ((TRUE == called_from_gsvit) && (set->sb.sigma_bx0 == -1)) {
                    set->sb.sigma_bx0 = 0.8 * (gdouble)set->sb.m_bx0 / (377 * set->sp.dx);
                    if (set->sc.verbose > 1)
                        printf("CPML sigma set to optimum value of %g\n", set->sb.sigma_bx0);
                }
                if (get_double(fr, &(set->sb.a_bx0), key))
                    goto exit;
                if (get_double(fr, &(set->sb.kappa_bx0), key))
                    goto exit;
            } else {
                fprintf(stderr, "Error: unsupported boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "BOUNDARY_XN") == 0) {
            fscanf(fr, "%s", value);
            if (strcmp(value, "none") == 0)
                set->sb.bxn = SV_BOUNDARY_NONE;
            else if (strcmp(value, "pec") == 0)
                set->sb.bxn = SV_BOUNDARY_PEC;
            else if (strcmp(value, "liao") == 0)
                set->sb.bxn = SV_BOUNDARY_LIAO;
            else if (strcmp(value, "cpml") == 0) {
                set->sb.bxn = SV_BOUNDARY_CPML;

                if (get_int(fr, &(set->sb.depth_bxn), key))
                    goto exit;
                if (get_int(fr, &(set->sb.m_bxn), key))
                    goto exit;
                if (get_double(fr, &(set->sb.sigma_bxn), key))
                    goto exit;
                if ((TRUE == called_from_gsvit) && (set->sb.sigma_bxn == -1)) {
                    set->sb.sigma_bxn = 0.8 * (gdouble)set->sb.m_bxn / (377 * set->sp.dx);
                    if (set->sc.verbose > 1)
                        printf("CPML sigma set to optimum value of %g\n", set->sb.sigma_bxn);
                }
                if (get_double(fr, &(set->sb.a_bxn), key))
                    goto exit;
                if (get_double(fr, &(set->sb.kappa_bxn), key))
                    goto exit;
            } else {
                fprintf(stderr, "Error: unsupported boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "BOUNDARY_Y0") == 0) {
            fscanf(fr, "%s", value);
            if (strcmp(value, "none") == 0)
                set->sb.by0 = SV_BOUNDARY_NONE;
            else if (strcmp(value, "pec") == 0)
                set->sb.by0 = SV_BOUNDARY_PEC;
            else if (strcmp(value, "liao") == 0)
                set->sb.by0 = SV_BOUNDARY_LIAO;
            else if (strcmp(value, "cpml") == 0) {
                set->sb.by0 = SV_BOUNDARY_CPML;

                if (get_int(fr, &(set->sb.depth_by0), key))
                    goto exit;
                if (get_int(fr, &(set->sb.m_by0), key))
                    goto exit;
                if (get_double(fr, &(set->sb.sigma_by0), key))
                    goto exit;
                if ((TRUE == called_from_gsvit) && (set->sb.sigma_by0 == -1)) {
                    set->sb.sigma_by0 = 0.8 * (gdouble)set->sb.m_by0 / (377 * set->sp.dy);
                    if (set->sc.verbose > 1)
                        printf("CPML sigma set to optimum value of %g\n", set->sb.sigma_by0);
                }
                if (get_double(fr, &(set->sb.a_by0), key))
                    goto exit;
                if (get_double(fr, &(set->sb.kappa_by0), key))
                    goto exit;
            } else {
                fprintf(stderr, "Error: unsupported boundary condition (%s)\n", value);
                goto exit;
            }

        } else if (strcmp(key, "BOUNDARY_YN") == 0) {
            fscanf(fr, "%s", value);
            if (strcmp(value, "none") == 0)
                set->sb.byn = SV_BOUNDARY_NONE;
            else if (strcmp(value, "pec") == 0)
                set->sb.byn = SV_BOUNDARY_PEC;
            else if (strcmp(value, "liao") == 0)
                set->sb.byn = SV_BOUNDARY_LIAO;
            else if (strcmp(value, "cpml") == 0) {
                set->sb.byn = SV_BOUNDARY_CPML;

                if (get_int(fr, &(set->sb.depth_byn), key))
                    goto exit;
                if (get_int(fr, &(set->sb.m_byn), key))
                    goto exit;
                if (get_double(fr, &(set->sb.sigma_byn), key))
                    goto exit;
                if ((TRUE == called_from_gsvit) && (set->sb.sigma_byn == -1)) {
                    set->sb.sigma_byn = 0.8 * (gdouble)set->sb.m_byn / (377 * set->sp.dy);
                    if (set->sc.verbose > 1)
                        printf("CPML sigma set to optimum value of %g\n", set->sb.sigma_byn);
                }
                if (get_double(fr, &(set->sb.a_byn), key))
                    goto exit;
                if (get_double(fr, &(set->sb.kappa_byn), key))
                    goto exit;
            } else {
                fprintf(stderr, "Error: unsupported boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "BOUNDARY_Z0") == 0) {
            fscanf(fr, "%s", value);
            if (strcmp(value, "none") == 0)
                set->sb.bz0 = SV_BOUNDARY_NONE;
            else if (strcmp(value, "pec") == 0)
                set->sb.bz0 = SV_BOUNDARY_PEC;
            else if (strcmp(value, "liao") == 0)
                set->sb.bz0 = SV_BOUNDARY_LIAO;
            else if (strcmp(value, "cpml") == 0) {
                set->sb.bz0 = SV_BOUNDARY_CPML;

                if (get_int(fr, &(set->sb.depth_bz0), key))
                    goto exit;
                if (get_int(fr, &(set->sb.m_bz0), key))
                    goto exit;
                if (get_double(fr, &(set->sb.sigma_bz0), key))
                    goto exit;
                if ((TRUE == called_from_gsvit) && (set->sb.sigma_bz0 == -1)) {
                    set->sb.sigma_bz0 = 0.8 * (gdouble)set->sb.m_bz0 / (377 * set->sp.dz);
                    if (set->sc.verbose > 1)
                        printf("CPML sigma set to optimum value of %g\n", set->sb.sigma_bz0);
                }
                if (get_double(fr, &(set->sb.a_bz0), key))
                    goto exit;
                if (get_double(fr, &(set->sb.kappa_bz0), key))
                    goto exit;
            } else {
                fprintf(stderr, "Error: unsupported boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "BOUNDARY_ZN") == 0) {
            fscanf(fr, "%s", value);
            if (strcmp(value, "none") == 0)
                set->sb.bzn = SV_BOUNDARY_NONE;
            else if (strcmp(value, "pec") == 0)
                set->sb.bzn = SV_BOUNDARY_PEC;
            else if (strcmp(value, "liao") == 0)
                set->sb.bzn = SV_BOUNDARY_LIAO;
            else if (strcmp(value, "cpml") == 0) {
                set->sb.bzn = SV_BOUNDARY_CPML;

                if (get_int(fr, &(set->sb.depth_bzn), key))
                    goto exit;
                if (get_int(fr, &(set->sb.m_bzn), key))
                    goto exit;
                if (get_double(fr, &(set->sb.sigma_bzn), key))
                    goto exit;
                if ((TRUE == called_from_gsvit) && (set->sb.sigma_bzn == -1)) {
                    set->sb.sigma_bzn = 0.8 * (gdouble)set->sb.m_bzn / (377 * set->sp.dz);
                    if (set->sc.verbose > 1)
                        printf("CPML sigma set to optimum value of %g\n", set->sb.sigma_bzn);
                }
                if (get_double(fr, &(set->sb.a_bzn), key))
                    goto exit;
                if (get_double(fr, &(set->sb.kappa_bzn), key))
                    goto exit;
            } else {
                fprintf(stderr, "Error: unsupported boundary condition (%s)\n", value);
                goto exit;
            }
        } else if (strcmp(key, "NFFF") == 0) {
            if (get_int(fr, &(set->sf.box_i0), key))
                goto exit;
            if (get_int(fr, &(set->sf.box_j0), key))
                goto exit;
            if (get_int(fr, &(set->sf.box_k0), key))
                goto exit;
            if (get_int(fr, &(set->sf.box_in), key))
                goto exit;
            if (get_int(fr, &(set->sf.box_jn), key))
                goto exit;
            if (get_int(fr, &(set->sf.box_kn), key))
                goto exit;
            if (get_int(fr, &(set->sf.sumfrom), key))
                goto exit;
            if (get_int(fr, &(set->sf.sumto), key))
                goto exit;

            //set->sf.valid = TRUE;
        } else if (strcmp(key, "NFFF_SKIP") == 0) {
            fscanf(fr, "%s", value);
            if (strcmp(value, "i0") == 0) {
                set->sf.box_boundary_skipi0 = 1;
                if (get_int(fr, &(set->sf.skipi0_jmin), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipi0_kmin), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipi0_jmax), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipi0_kmax), key))
                    goto exit;
            } else if (strcmp(value, "in") == 0) {
                set->sf.box_boundary_skipin = 1;
                if (get_int(fr, &(set->sf.skipin_jmin), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipin_kmin), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipin_jmax), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipin_kmax), key))
                    goto exit;
            } else if (strcmp(value, "j0") == 0) {
                set->sf.box_boundary_skipj0 = 1;
                if (get_int(fr, &(set->sf.skipj0_imin), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipj0_kmin), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipj0_imax), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipj0_kmax), key))
                    goto exit;
            } else if (strcmp(value, "jn") == 0) {
                set->sf.box_boundary_skipjn = 1;
                if (get_int(fr, &(set->sf.skipjn_imin), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipjn_kmin), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipjn_imax), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipjn_kmax), key))
                    goto exit;
            } else if (strcmp(value, "k0") == 0) {
                set->sf.box_boundary_skipk0 = 1;
                if (get_int(fr, &(set->sf.skipk0_imin), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipk0_jmin), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipk0_imax), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipk0_jmax), key))
                    goto exit;
            } else if (strcmp(value, "kn") == 0) {
                set->sf.box_boundary_skipkn = 1;
                if (get_int(fr, &(set->sf.skipkn_imin), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipkn_jmin), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipkn_imax), key))
                    goto exit;
                if (get_int(fr, &(set->sf.skipkn_jmax), key))
                    goto exit;
            }
        } else if (strcmp(key, "NFFF_RAMAHI_POINT") == 0) {
            if (set->sf.nrs == 0) {
                set->sf.nrs = 1;
                set->sf.ri = (gint *)g_malloc(set->sf.nrs * sizeof(gint));
                set->sf.rj = (gint *)g_malloc(set->sf.nrs * sizeof(gint));
                set->sf.rk = (gint *)g_malloc(set->sf.nrs * sizeof(gint));
                set->sf.individual = (gint *)g_malloc(set->sf.nrs * sizeof(gint));
                set->sf.source_filename = (gchar **)g_malloc(set->sf.nrs * sizeof(gchar *));
            } else { /*realloc*/
                set->sf.nrs += 1;
                set->sf.ri = (gint *)g_realloc(set->sf.ri, set->sf.nrs * sizeof(gint));
                set->sf.rj = (gint *)g_realloc(set->sf.rj, set->sf.nrs * sizeof(gint));
                set->sf.rk = (gint *)g_realloc(set->sf.rk, set->sf.nrs * sizeof(gint));
                set->sf.individual = (gint *)g_realloc(set->sf.individual, set->sf.nrs * sizeof(gint));
                set->sf.source_filename = (gchar **)g_realloc(set->sf.source_filename, set->sf.nrs * sizeof(gchar *));
            }
            if (get_int(fr, &(set->sf.ri[set->sf.nrs - 1]), key))
                goto exit;
            if (get_int(fr, &(set->sf.rj[set->sf.nrs - 1]), key))
                goto exit;
            if (get_int(fr, &(set->sf.rk[set->sf.nrs - 1]), key))
                goto exit;
            fscanf(fr, "%100[^\n\r]", buffer);
            set->sf.source_filename[set->sf.nrs - 1] = g_strstrip(g_strdup(buffer));
            set->sf.individual[set->sf.nrs - 1] = 0;
        } else if (strcmp(key, "PERIODIC_NFFF") == 0) {
            if (get_int(fr, &(set->spf.box_i0), key))
                goto exit;
            if (get_int(fr, &(set->spf.box_j0), key))
                goto exit;
            if (get_int(fr, &(set->spf.box_k0), key))
                goto exit;
            if (get_int(fr, &(set->spf.box_in), key))
                goto exit;
            if (get_int(fr, &(set->spf.box_jn), key))
                goto exit;
            if (get_int(fr, &(set->spf.box_kn), key))
                goto exit;
            if (get_int(fr, &(set->spf.pimin), key))
                goto exit;
            if (get_int(fr, &(set->spf.pjmin), key))
                goto exit;
            if (get_int(fr, &(set->spf.pimax), key))
                goto exit;
            if (get_int(fr, &(set->spf.pjmax), key))
                goto exit;

            //set->spf.valid = TRUE;
        } else if (strcmp(key, "PERIODIC_NFFF_POSTPROCESS") == 0) {
            if (get_int(fr, &(set->spf.postprocess), key))
                goto exit;
            if (get_int(fr, &(set->spf.ppstart), key))
                goto exit;
        } else if (strcmp(key, "PERIODIC_NFFF_SKIP") == 0) {
            fscanf(fr, "%s", value);
            if (strcmp(value, "k0") == 0) {
                set->spf.box_boundary_skipk0 = 1;
                if (get_int(fr, &(set->spf.skipk0_imin), key))
                    goto exit;
                if (get_int(fr, &(set->spf.skipk0_jmin), key))
                    goto exit;
                if (get_int(fr, &(set->spf.skipk0_imax), key))
                    goto exit;
                if (get_int(fr, &(set->spf.skipk0_jmax), key))
                    goto exit;
            } else if (strcmp(value, "kn") == 0) {
                set->spf.box_boundary_skipkn = 1;
                if (get_int(fr, &(set->spf.skipkn_imin), key))
                    goto exit;
                if (get_int(fr, &(set->spf.skipkn_jmin), key))
                    goto exit;
                if (get_int(fr, &(set->spf.skipkn_imax), key))
                    goto exit;
                if (get_int(fr, &(set->spf.skipkn_jmax), key))
                    goto exit;
            }
        } else if (strcmp(key, "PERIODIC_NFFF_RAMAHI_POINT") == 0) {
            if (set->spf.nrs == 0) {
                set->spf.nrs = 1;
                set->spf.ri = (gint *)g_malloc(set->spf.nrs * sizeof(gint));
                set->spf.rj = (gint *)g_malloc(set->spf.nrs * sizeof(gint));
                set->spf.rk = (gint *)g_malloc(set->spf.nrs * sizeof(gint));
                set->spf.individual = (gint *)g_malloc(set->spf.nrs * sizeof(gint));
                set->spf.source_filename = (gchar **)g_malloc(set->spf.nrs * sizeof(gchar *));
            } else { /*realloc*/
                set->spf.nrs += 1;
                set->spf.ri = (gint *)g_realloc(set->spf.ri, set->spf.nrs * sizeof(gint));
                set->spf.rj = (gint *)g_realloc(set->spf.rj, set->spf.nrs * sizeof(gint));
                set->spf.rk = (gint *)g_realloc(set->spf.rk, set->spf.nrs * sizeof(gint));
                set->spf.individual = (gint *)g_realloc(set->spf.individual, set->spf.nrs * sizeof(gint));
                set->spf.source_filename = (gchar **)g_realloc(set->spf.source_filename, set->spf.nrs * sizeof(gchar *));
            }
            if (get_int(fr, &(set->spf.ri[set->spf.nrs - 1]), key))
                goto exit;
            if (get_int(fr, &(set->spf.rj[set->spf.nrs - 1]), key))
                goto exit;
            if (get_int(fr, &(set->spf.rk[set->spf.nrs - 1]), key))
                goto exit;
            fscanf(fr, "%100[^\n\r]", buffer);
            set->spf.source_filename[set->spf.nrs - 1] = g_strstrip(g_strdup(buffer));
            set->spf.individual[set->spf.nrs - 1] = 0; /*not part of any set*/
        } else if (strcmp(key, "PERIODIC_NFFF_SPHERICAL_AREA") == 0) {

            if (get_int(fr, &thetares, key))
                goto exit;
            if (get_int(fr, &phires, key))
                goto exit;
            if (get_int(fr, &radius, key))
                goto exit;
            if (get_double(fr, &thetafrom, key))
                goto exit;
            if (get_double(fr, &phifrom, key))
                goto exit;
            if (get_double(fr, &thetato, key))
                goto exit;
            if (get_double(fr, &phito, key))
                goto exit;
            if (get_int(fr, &savefile, key))
                goto exit;

            if (set->spf.nareas < 20) {  //for xsvit only
                set->spf.area_thetares[set->spf.nareas] = thetares;
                set->spf.area_phires[set->spf.nareas] = phires;
                set->spf.area_radius[set->spf.nareas] = radius;
                set->spf.area_thetafrom[set->spf.nareas] = thetafrom;
                set->spf.area_phifrom[set->spf.nareas] = phifrom;
                set->spf.area_thetato[set->spf.nareas] = thetato;
                set->spf.area_phito[set->spf.nareas] = phito;
                set->spf.area_savefile[set->spf.nareas] = savefile;
                set->spf.nareas++;
            }

            if (TRUE == called_from_gsvit) {
                if (set->spf.nrs == 0) {
                    set->spf.nrs = thetares*phires;
                    set->spf.ri = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->spf.rj = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->spf.rk = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->spf.individual = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->spf.source_filename = (gchar **)g_malloc(thetares*phires * sizeof(gchar *));
                } else {
                    set->spf.nrs += thetares*phires;
                    set->spf.ri = (gint *)g_realloc(set->spf.ri, set->spf.nrs * sizeof(gint));
                    set->spf.rj = (gint *)g_realloc(set->spf.rj, set->spf.nrs * sizeof(gint));
                    set->spf.rk = (gint *)g_realloc(set->spf.rk, set->spf.nrs * sizeof(gint));
                    set->spf.individual = (gint *)g_realloc(set->spf.individual, set->spf.nrs * sizeof(gint));
                    set->spf.source_filename = (gchar **)g_realloc(set->spf.source_filename, set->spf.nrs * sizeof(gchar *));
                }

                n = set->spf.nrs - thetares*phires;
                printf("n %d  %d %d\n", n, thetares, phires);

                if (set->spf.nsets == 0) {
                    set->spf.nsets = 1;    /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                    set->spf.setxres = (gint *)g_malloc((set->spf.nsets + 1) * sizeof(gint));
                    set->spf.setyres = (gint *)g_malloc((set->spf.nsets + 1) * sizeof(gint));
                } else {
                    set->spf.nsets += 1;
                    set->spf.setxres = (gint *)g_realloc(set->spf.setxres, (set->spf.nsets + 1) * sizeof(gint));
                    set->spf.setyres = (gint *)g_realloc(set->spf.setyres, (set->spf.nsets + 1) * sizeof(gint));
                }

                set->spf.setxres[set->spf.nsets] = thetares;
                set->spf.setyres[set->spf.nsets] = phires;

                for (i = 0; i < thetares; i++) {
                    for (j = 0; j < phires; j++) {
                        if (thetares > 1)
                            theta = (thetafrom + (gdouble)(i)*(thetato - thetafrom) / (thetares - 1));
                        else
                            theta = thetafrom;
                        if (phires > 1)
                            phi = (phifrom + (gdouble)(j)*(phito - phifrom) / (phires - 1));
                        else
                            phi = phifrom;

                        set->spf.ri[n] = (gint)(radius*sin(theta)*cos(phi)) + set->sp.xres / 2;
                        set->spf.rj[n] = (gint)(radius*sin(theta)*sin(phi)) + set->sp.yres / 2;
                        set->spf.rk[n] = (gint)(radius*cos(theta)) + set->sp.zres / 2;

                        set->spf.individual[n] = set->spf.nsets; /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                        if (!savefile)
                            set->spf.source_filename[n] = NULL;
                        else
                            set->spf.source_filename[n] = g_strdup_printf("xnff_set%03d_%03d_%03d", set->spf.nsets, i, j);

                        n++;
                    } // j
                } // i
            }
        } else if (strcmp(key, "NFFF_SPHERICAL_AREA") == 0) {
            if (get_int(fr, &thetares, key))
                goto exit;
            if (get_int(fr, &phires, key))
                goto exit;
            if (get_int(fr, &radius, key))
                goto exit;
            if (get_double(fr, &thetafrom, key))
                goto exit;
            if (get_double(fr, &phifrom, key))
                goto exit;
            if (get_double(fr, &thetato, key))
                goto exit;
            if (get_double(fr, &phito, key))
                goto exit;
            if (get_int(fr, &savefile, key))
                goto exit;

            if (set->sf.nareas < 20) {  //for xsvit only
                set->sf.area_thetares[set->sf.nareas] = thetares;
                set->sf.area_phires[set->sf.nareas] = phires;
                set->sf.area_radius[set->sf.nareas] = radius;
                set->sf.area_thetafrom[set->sf.nareas] = thetafrom;
                set->sf.area_phifrom[set->sf.nareas] = phifrom;
                set->sf.area_thetato[set->sf.nareas] = thetato;
                set->sf.area_phito[set->sf.nareas] = phito;
                set->sf.area_savefile[set->sf.nareas] = savefile;
                set->sf.nareas++;
            }

            if (TRUE == called_from_gsvit) {
                if (set->sf.nrs == 0) {
                    set->sf.nrs = thetares*phires;
                    set->sf.ri = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->sf.rj = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->sf.rk = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->sf.individual = (gint *)g_malloc(thetares*phires * sizeof(gint));
                    set->sf.source_filename = (gchar **)g_malloc(thetares*phires * sizeof(gchar *));
                } else {
                    set->sf.nrs += thetares*phires;
                    set->sf.ri = (gint *)g_realloc(set->sf.ri, set->sf.nrs * sizeof(gint));
                    set->sf.rj = (gint *)g_realloc(set->sf.rj, set->sf.nrs * sizeof(gint));
                    set->sf.rk = (gint *)g_realloc(set->sf.rk, set->sf.nrs * sizeof(gint));
                    set->sf.individual = (gint *)g_realloc(set->sf.individual, set->sf.nrs * sizeof(gint));
                    set->sf.source_filename = (gchar **)g_realloc(set->sf.source_filename, set->sf.nrs * sizeof(gchar *));
                }
                n = set->sf.nrs - thetares*phires;

                if (set->sf.nsets == 0) {
                    set->sf.nsets = 1;    /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                    set->sf.setxres = (gint *)g_malloc((set->sf.nsets + 1) * sizeof(gint));
                    set->sf.setyres = (gint *)g_malloc((set->sf.nsets + 1) * sizeof(gint));
                } else {
                    set->sf.nsets += 1;
                    set->sf.setxres = (gint *)g_realloc(set->sf.setxres, (set->sf.nsets + 1) * sizeof(gint));
                    set->sf.setyres = (gint *)g_realloc(set->sf.setyres, (set->sf.nsets + 1) * sizeof(gint));
                }
                set->sf.setxres[set->sf.nsets] = thetares;
                set->sf.setyres[set->sf.nsets] = phires;

                for (i = 0; i < thetares; i++) {
                    for (j = 0; j < phires; j++) {
                        if (thetares > 1)
                            theta = (thetafrom + (gdouble)(i)*(thetato - thetafrom) / (thetares - 1));
                        else
                            theta = thetafrom;
                        if (phires > 1)
                            phi = (phifrom + (gdouble)(j)*(phito - phifrom) / (phires - 1));
                        else
                            phi = phifrom;

                        set->sf.ri[n] = (gint)(radius*sin(theta)*cos(phi)) + set->sp.xres / 2;
                        set->sf.rj[n] = (gint)(radius*sin(theta)*sin(phi)) + set->sp.yres / 2;
                        set->sf.rk[n] = (gint)(radius*cos(theta)) + set->sp.zres / 2;

                        set->sf.individual[n] = set->sf.nsets; /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                        if (!savefile)
                            set->sf.source_filename[n] = NULL;
                        else
                            set->sf.source_filename[n] = g_strdup_printf("xnff_set%03d_%03d_%03d", set->sf.nsets, i, j);
                        n++;
                    } // j
                } // i
            }
        } else if (strcmp(key, "NFFF_PLANAR_AREA") == 0) {
            if (get_int(fr, &ijres, key))
                goto exit;
            if (get_int(fr, &jkres, key))
                goto exit;
            if (get_double(fr, &ijfrom, key))
                goto exit;
            if (get_double(fr, &jkfrom, key))
                goto exit;
            if (get_double(fr, &ijto, key))
                goto exit;
            if (get_double(fr, &jkto, key))
                goto exit;
            if (get_int(fr, &orientation, key))
                goto exit;
            if (get_double(fr, &distance, key))
                goto exit;
            if (get_int(fr, &savefile, key))
                goto exit;

            // printf("or %d\n", orientation);
            if (set->sf.nsquares < 20) {  //for xsvit only
                set->sf.square_ijres[set->sf.nsquares] = thetares;
                set->sf.square_jkres[set->sf.nsquares] = phires;
                set->sf.square_distance[set->sf.nsquares] = distance;
                set->sf.square_ijfrom[set->sf.nsquares] = thetafrom;
                set->sf.square_jkfrom[set->sf.nsquares] = phifrom;
                set->sf.square_ijto[set->sf.nsquares] = thetato;
                set->sf.square_jkto[set->sf.nsquares] = phito;
                set->sf.square_orientation[set->sf.nsquares] = orientation;
                set->sf.square_savefile[set->sf.nsquares] = savefile;
                set->sf.nsquares++;
            }

            if (TRUE == called_from_gsvit) {
                if (set->sf.nrs == 0) {
                    set->sf.nrs = ijres*jkres;
                    set->sf.ri = (gint *)g_malloc(ijres*jkres * sizeof(gint));
                    set->sf.rj = (gint *)g_malloc(ijres*jkres * sizeof(gint));
                    set->sf.rk = (gint *)g_malloc(ijres*jkres * sizeof(gint));
                    set->sf.individual = (gint *)g_malloc(ijres*jkres * sizeof(gint));
                    set->sf.source_filename = (gchar **)g_malloc(ijres*jkres * sizeof(gchar *));
                } else {
                    set->sf.nrs += ijres*jkres;
                    set->sf.ri = (gint *)g_realloc(set->sf.ri, set->sf.nrs * sizeof(gint));
                    set->sf.rj = (gint *)g_realloc(set->sf.rj, set->sf.nrs * sizeof(gint));
                    set->sf.rk = (gint *)g_realloc(set->sf.rk, set->sf.nrs * sizeof(gint));
                    set->sf.individual = (gint *)g_realloc(set->sf.individual, set->sf.nrs * sizeof(gint));
                    set->sf.source_filename = (gchar **)g_realloc(set->sf.source_filename, set->sf.nrs * sizeof(gchar *));
                }
                n = set->sf.nrs - ijres*jkres;

                if (set->sf.nsets == 0) {
                    set->sf.nsets = 1;    /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                    set->sf.setxres = (gint *)g_malloc((set->sf.nsets + 1) * sizeof(gint));
                    set->sf.setyres = (gint *)g_malloc((set->sf.nsets + 1) * sizeof(gint));
                } else {
                    set->sf.nsets += 1;
                    set->sf.setxres = (gint *)g_realloc(set->sf.setxres, (set->sf.nsets + 1) * sizeof(gint));
                    set->sf.setyres = (gint *)g_realloc(set->sf.setyres, (set->sf.nsets + 1) * sizeof(gint));
                }

                set->sf.setxres[set->sf.nsets] = ijres;
                set->sf.setyres[set->sf.nsets] = jkres;

                for (i = 0; i < ijres; i++) {
                    for (j = 0; j < jkres; j++) {
                        if (orientation == 0) { //x orientation
                            set->sf.ri[n] = (gint)(distance / set->sp.dx) + set->sp.xres / 2;
                            set->sf.rj[n] = (gint)((ijfrom + (ijto - ijfrom) * ((gdouble)i) / ((gdouble)ijres - 1)) / set->sp.dy) + set->sp.yres / 2;
                            set->sf.rk[n] = (gint)((jkfrom + (jkto - jkfrom) * ((gdouble)j) / ((gdouble)jkres - 1)) / set->sp.dz) + set->sp.zres / 2;

                            //                      printf("OUT_POINT\nAll 1 %d %d %d point_%.3d_%.3d\n\n", set->sf.ri[n], set->sf.rj[n], set->sf.rk[n], i, j); 
                            //                      printf("NFFF_RAMAHI_POINT\n%d %d %d fpoint_%.3d_%.3d\n\n", set->sf.ri[n], set->sf.rj[n], set->sf.rk[n], i, j); 

                                                    //printf("%d %d %d\n", set->sf.ri[n], set->sf.rj[n], set->sf.rk[n]);

                            set->sf.individual[n] = set->sf.nsets; /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                            if (!savefile)
                                set->sf.source_filename[n] = NULL;
                            else
                                set->sf.source_filename[n] = g_strdup_printf("xnff_set%03d_%03d_%03d", set->sf.nsets, i, j);
                            n++;
                        } else if (orientation == 1) { //y orientation
                            set->sf.ri[n] = (gint)((ijfrom + (ijto - ijfrom) * ((gdouble)i) / ((gdouble)ijres - 1)) / set->sp.dx) + set->sp.xres / 2;
                            set->sf.rj[n] = (gint)(distance / set->sp.dy) + set->sp.yres / 2;
                            set->sf.rk[n] = (gint)((jkfrom + (jkto - jkfrom) * ((gdouble)j) / ((gdouble)jkres - 1)) / set->sp.dz) + set->sp.zres / 2;

                            set->sf.individual[n] = set->sf.nsets; /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                            if (!savefile)
                                set->sf.source_filename[n] = NULL;
                            else
                                set->sf.source_filename[n] = g_strdup_printf("xnff_set%03d_%03d_%03d", set->sf.nsets, i, j);
                            n++;
                        } else if (orientation == 2) { //z orientation
                            set->sf.ri[n] = (gint)((ijfrom + (ijto - ijfrom) * ((gdouble)i) / ((gdouble)ijres - 1)) / set->sp.dx) + set->sp.xres / 2;
                            set->sf.rj[n] = (gint)((jkfrom + (jkto - jkfrom) * ((gdouble)j) / ((gdouble)jkres - 1)) / set->sp.dy) + set->sp.yres / 2;
                            set->sf.rk[n] = (gint)(distance / set->sp.dz) + set->sp.zres / 2;

                            set->sf.individual[n] = set->sf.nsets; /*we are counting from 1 - as the individual data are 0 and sets start from 1*/
                            if (!savefile)
                                set->sf.source_filename[n] = NULL;
                            else
                                set->sf.source_filename[n] = g_strdup_printf("xnff_set%03d_%03d_%03d", set->sf.nsets, i, j);
                            n++;
                        } else {
                            fprintf(stderr, "Error: unsupported NFFF area orientation (%d)\n", orientation);
                            goto exit;
                        }
                    } // j
                } // i
            }
        } else {
            fprintf(stderr, "Error: unsupported or unknown key (%s)\n", key);
            goto exit;
        }
    }

    if (set->smb.bxnpos == -1)
        set->smb.bxnpos = set->sp.xres;
    if (set->smb.bynpos == -1)
        set->smb.bynpos = set->sp.yres;
    if (set->smb.bznpos == -1)
        set->smb.bznpos = set->sp.zres;

    fclose(fr);

    return 1;

exit:
    fclose(fr);
    return 0;
} /* parse_settings */
#endif

int
parse_settings_mat(gchar *filename, SvSet *set_par, SvSetMat *set_mat, gboolean called_from_gsvit)
{
    gint type, testout, xres, yres, zres, j, k;
    gdouble xreal, yreal, zreal, val;
    GError* glib_error = 0;
    FILE *fr = NULL;
    gint ret;
    SvSphere sx;
    SvVoxel vx;
    SvCylinder cx;
    SvCone cnx;
    SvRCone rcnx;
    SvTetrahedron tx;
    SvGwydd gx;
    //gint ngtot = 0;
    SvMatProp mat;
    GLfloat xshift, yshift, zshift, xmult, ymult, zmult;
    gchar buff[256], matfile[256], filebase[256];
    gint attribute_pos, attribute_val;
    gint i, gwydd_channel;
    gdouble mask;
    SvDCube *buf;
    gdouble min, max;



    set_mat->nmats = 0;
    set_mat->materials = (gchar **)g_malloc(100 * sizeof(gchar *)); //FIXME: limited size array
    set_mat->overpos = (gint *)g_malloc(100 * sizeof(gint));
    for (i = 0; i < 100; i++) {
        set_mat->overpos[i] = -1;
        set_mat->materials[i] = NULL;
    }

    if ((testout = test_script_file(filename)) == SV_MATMODE_DATABASE) {
        if (override_matfile_from_database(filename, 600e-9, 600e-9, 600e-9, 12345, &(set_mat->nmats), set_mat->materials, set_mat->overpos, TRUE)) { //FIXME, here should be set->sm.localfile
            fprintf(stderr, "Error parsing material from database.\n");
            //g_snprintf(buff, sizeof(buff), "Error parsing material from database.");
            //gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
        } else {
            fprintf(stderr, "Material loaded from database (%d entries listed).\n", set_mat->nmats);
            //g_snprintf(buff, sizeof(buff), "Material loaded from database (%d entries listed).", set->nmats);
            //gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);

            /*for (i=0; i<xgc->nmats; i++)
            printf("Material %d: %s to object %d\n", i, xgc->materials[i], xgc->overpos[i]); */
        }

        //FIXME! use real wavelengths
        g_snprintf(matfile, sizeof(matfile), "tmp_matfile_%06d", 12345);
        filename = g_strdup(matfile);
    }
    if (testout == SV_MATMODE_ERROR) {
        printf("Error in material file\n");

        fprintf(stderr, "Error parsing material file.\n");
        //g_snprintf(buff, sizeof(buff), "Error parsing material file.");
        //gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);

        return 0;
    }


    fr = fopen(filename, "r");
    if (fr == NULL) {
        fprintf(stderr, "Error: cannot open material vector file for loading.\n");
        g_snprintf(buff, sizeof(buff), "Error parsing material file, see console for output\n");
        //gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
    } else {
        /* alloc data*/
        alloc_set_mat(set_mat);

        while ((ret = fscanf(fr, "%d", &type)) != EOF) { //TODO use proper locale handling functions here
            if (ret != 1) {
                fprintf(stderr, "Error parsing material file\n");
                goto exit;
            }

            switch (type) {
            case SV_ENTITY_SPHERE:
                if (scan_point(fr, sx.pnt1, sx.pnt1 + 1, sx.pnt1 + 2)) {
                    fprintf(stderr, "Error: cannot scan sphere center from material file\n");
                    goto exit;
                }
                if (scan_radius(fr, &sx.radius)) {
                    fprintf(stderr, "Error: cannot scan sphere radius from material file\n");
                    goto exit;
                }
                if (scan_material(fr, &sx.mat, 1)) {
                    fprintf(stderr, "Error: cannot scan sphere material from material file\n");
                    goto exit;
                }
                //if (sx.mat.type!=0) sv_pool_add_material(mp, &sx.mat);
                sx.n = set_mat->ngtot++;
                g_array_append_val(set_mat->spheres, sx);
                break;

            case SV_ENTITY_CYLINDER:
                if (scan_point(fr, cx.pnt1, cx.pnt1 + 1, cx.pnt1 + 2)) {
                    fprintf(stderr, "Error: cannot scan cylinder center1 from material file\n");
                    goto exit;
                }
                if (scan_point(fr, cx.pnt2, cx.pnt2 + 1, cx.pnt2 + 2)) {
                    fprintf(stderr, "Error: cannot scan cylinder center2 from material file\n");
                    goto exit;
                }
                if (scan_radius(fr, &cx.radius)) {
                    fprintf(stderr, "Error: cannot scan cylinder radius from material file\n");
                    goto exit;
                }
                if (scan_material(fr, &cx.mat, 1)) {
                    fprintf(stderr, "Error: cannot scan cylinder material from material file\n");
                    goto exit;
                }
                //if (cx.mat.type!=0) sv_pool_add_material(mp, &cx.mat);
                cx.n = set_mat->ngtot++;
                g_array_append_val(set_mat->cylinders, cx);
                break;

            case SV_ENTITY_CONE:
                if (scan_point(fr, cnx.pnt1, cnx.pnt1 + 1, cnx.pnt1 + 2)) {
                    fprintf(stderr, "Error: cannot scan cone center1 from material file\n");
                    goto exit;
                }
                if (scan_point(fr, cnx.pnt2, cnx.pnt2 + 1, cnx.pnt2 + 2)) {
                    fprintf(stderr, "Error: cannot scan cone center2 from material file\n");
                    goto exit;
                }
                if (scan_radius(fr, &cnx.radius)) {
                    fprintf(stderr, "Error: cannot scan cone radius from material file\n");
                    goto exit;
                }
                if (scan_material(fr, &cnx.mat, 1)) {
                    fprintf(stderr, "Error: cannot scan cone material from material file\n");
                    goto exit;
                }
                //if (cnx.mat.type!=0) sv_pool_add_material(mp, &cnx.mat);
                cnx.n = set_mat->ngtot++;
                g_array_append_val(set_mat->cones, cnx);
                break;

            case SV_ENTITY_RCONE:
                if (scan_point(fr, rcnx.pnt1, rcnx.pnt1 + 1, rcnx.pnt1 + 2)) {
                    fprintf(stderr, "Error: cannot scan cut cone center1 from material file\n");
                    goto exit;
                }
                if (scan_point(fr, rcnx.pnt2, rcnx.pnt2 + 1, rcnx.pnt2 + 2)) {
                    fprintf(stderr, "Error: cannot scan cut cone center2 from material file\n");
                    goto exit;
                }
                if (scan_radius(fr, &rcnx.radius1)) {
                    fprintf(stderr, "Error: cannot scan cut cone radius1 from material file\n");
                    goto exit;
                }
                if (scan_radius(fr, &rcnx.radius2)) {
                    fprintf(stderr, "Error: cannot scan cut cone radius2 from material file\n");
                    goto exit;
                }
                if (scan_material(fr, &rcnx.mat, 1)) {
                    fprintf(stderr, "Error: cannot scan cut cone material from material file\n");
                    goto exit;
                }
                //if (rcnx.mat.type!=0) sv_pool_add_material(mp, &rcnx.mat);
                rcnx.n = set_mat->ngtot++;
                g_array_append_val(set_mat->rcones, rcnx);
                break;

            case SV_ENTITY_VOXEL:
                if (scan_point(fr, vx.pnt1, vx.pnt1 + 1, vx.pnt1 + 2)) {
                    fprintf(stderr, "Error: cannot scan voxel start vertex from material file\n");
                    goto exit;
                }
                if (scan_point(fr, vx.pnt2, vx.pnt2 + 1, vx.pnt2 + 2)) {
                    fprintf(stderr, "Error: cannot scan voxel end vertex from material file\n");
                    goto exit;
                }
                if (scan_material(fr, &vx.mat, 1)) {
                    fprintf(stderr, "Error: cannot scan voxel material from material file\n");
                    goto exit;
                }
                //if (vx.mat.type!=0) sv_pool_add_material(mp, &vx.mat);
                vx.n = set_mat->ngtot++;
                g_array_append_val(set_mat->voxels, vx);
                break;

            case SV_ENTITY_TETRAHEDRON:
                if (scan_point(fr, tx.pnt1, tx.pnt1 + 1, tx.pnt1 + 2)) {
                    fprintf(stderr, "Error: cannot scan tetrahedron point1 from material file\n");
                    goto exit;
                }
                if (scan_point(fr, tx.pnt2, tx.pnt2 + 1, tx.pnt2 + 2)) {
                    fprintf(stderr, "Error: cannot scan tetrahedron point2 from material file\n");
                    goto exit;
                }
                if (scan_point(fr, tx.pnt3, tx.pnt3 + 1, tx.pnt3 + 2)) {
                    fprintf(stderr, "Error: cannot scan tetrahedron point3 from material file\n");
                    goto exit;
                }
                if (scan_point(fr, tx.pnt4, tx.pnt4 + 1, tx.pnt4 + 2)) {
                    fprintf(stderr, "Error: cannot scan tetrahedron point4 from material file\n");
                    goto exit;
                }
                if (scan_material(fr, &tx.mat, 1)) {
                    fprintf(stderr, "Error: cannot scan tetrahedron material from material file\n");
                    goto exit;
                }
                //if (tx.mat.type!=0) sv_pool_add_material(mp, &tx.mat);
                tx.setpart = 0;
                tx.n = set_mat->ngtot++;
                g_array_append_val(set_mat->tetrahedrons, tx);
                break;

            case SV_ENTITY_TETGEN:
                fscanf(fr, "%s", filebase);             //TODO this should be locale dependent
                if (scan_int(fr, &attribute_pos)) {
                    fprintf(stderr, "Error: cannot scan tetrahedral mesh attribute_number from material file\n");
                    goto exit;
                }
                if (scan_int(fr, &attribute_val)) {
                    fprintf(stderr, "Error: cannot scan tetrahedral mesh material_index from material file\n");
                    goto exit;
                }
                if (scan_point(fr, &xshift, &yshift, &zshift)) {
                    fprintf(stderr, "Error: cannot scan tetrahedral mesh shift from material file\n");
                    goto exit;
                }
                if (scan_point(fr, &xmult, &ymult, &zmult)) {
                    fprintf(stderr, "Error: cannot scan tetrahedral mesh mult from material file\n");
                    goto exit;
                }
                if (scan_material(fr, &mat, 1)) {
                    fprintf(stderr, "Error: cannot scan tetrahedral mesh material from material file\n");
                    goto exit;
                }
                //if (mat.type!=0) sv_pool_add_material(mp, &mat);

                set_mat->tetgen_attribute_pos[set_mat->ntetgens] = attribute_pos;
                set_mat->tetgen_attribute_val[set_mat->ntetgens] = attribute_val;
                set_mat->tetgen_start[set_mat->ntetgens] = set_mat->ngtot;
                set_mat->tetgen_xshift[set_mat->ntetgens] = xshift;
                set_mat->tetgen_yshift[set_mat->ntetgens] = yshift;
                set_mat->tetgen_zshift[set_mat->ntetgens] = zshift;
                set_mat->tetgen_xmult[set_mat->ntetgens] = xmult;
                set_mat->tetgen_ymult[set_mat->ntetgens] = ymult;
                set_mat->tetgen_zmult[set_mat->ntetgens] = zmult;
                g_snprintf(set_mat->tetgen_filebase[set_mat->ntetgens], 100, "%s", filebase);
                if (load_tets(filebase, set_mat->tetrahedrons, attribute_pos, attribute_val, xshift, yshift, zshift, xmult, ymult, zmult, 3, &mat, &set_mat->ngtot) != 0) {
                    fprintf(stderr, "Error parsing tetrahedral mesh from %s\n", filebase);
                    //g_snprintf(buff, sizeof(buff), "Error parsing tetrahedral mesh from %s", filebase);
                    //gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
                    goto exit;
                }
                else {
                    fprintf(stderr, "Loaded tetrahedral mesh from %s\n", filebase);
                    //g_snprintf(buff, sizeof(buff), "Loaded tetrahedral mesh from %s", filebase);
                    //gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
                }
                set_mat->tetgen_end[set_mat->ntetgens] = set_mat->ngtot;
                set_mat->ntetgens++;
                break;

            case SV_ENTITY_GWYDDION:
                fscanf(fr, "%s", filebase);             //TODO this should be locale dependent
                g_snprintf(gx.filebase, 256, "%s", filebase);
                if (scan_int(fr, &gwydd_channel)) {
                    fprintf(stderr, "Error: cannot scan Gwyddion field channel from material file\n");
                    goto exit;
                }
                gx.channel = gwydd_channel;
                if (scan_int(fr, &gx.mask)) {
                    fprintf(stderr, "Error: cannot scan Gwyddion field mask from material file\n");
                    goto exit;
                }
                if (scan_int(fr, &gx.i)) {
                    fprintf(stderr, "Error: cannot scan Gwyddion field i from material file\n");
                    goto exit;
                }
                if (scan_int(fr, &gx.j)) {
                    fprintf(stderr, "Error: cannot scan Gwyddion field j from material file\n");
                    goto exit;
                }
                if (scan_int(fr, &gx.k)) {
                    fprintf(stderr, "Error: cannot scan Gwyddion field k from material file\n");
                    goto exit;
                }
                if (scan_point(fr, &gx.xoffset, &gx.yoffset, &gx.zoffset)) {
                    fprintf(stderr, "Error: cannot scan Gwyddion field offset from material file\n");
                    goto exit;
                }
                if (scan_int(fr, &gx.depth)) {
                    fprintf(stderr, "Error: cannot scan Gwyddion field depth from material file\n");
                    goto exit;
                }
                if (scan_material(fr, &gx.mat, 1)) {
                    fprintf(stderr, "Error: cannot scan Gwyddion field material from material file\n");
                    goto exit;
                }
                //if (mat.type!=0) sv_pool_add_material(mp, &gx.mat);

                load_gwydd(filebase, &gx, 0, gwydd_channel);

                if (gx.dfield == NULL || (!GWY_IS_DATA_FIELD(gx.dfield))) {
                    fprintf(stderr, "Error: no valid Gwyddion datafield\n");
                }
                else {
                    //                              printf("loaded %d %d\n", gwy_data_field_get_xres(gx.dfield),  gwy_data_field_get_yres(gx.dfield));
                    gx.n = set_mat->ngtot++;
                    g_array_append_val(set_mat->gwydds, gx);

                    if (set_mat->gwyddata)
                        sv_dcube_free(set_mat->gwyddata);
                    set_mat->gwyddata = sv_dcube_new(set_par->sp.xres, set_par->sp.yres, set_par->sp.zres,
                        set_par->sp.xres, set_par->sp.yres, set_par->sp.zres, 1);
                    buf = sv_dcube_new(set_par->sp.xres, set_par->sp.yres, set_par->sp.zres,
                        set_par->sp.xres, set_par->sp.yres, set_par->sp.zres, 1);

                    xres = set_par->sp.xres;
                    xreal = set_par->sp.dx*xres;

                    yres = set_par->sp.yres;
                    yreal = set_par->sp.dy*yres;

                    zres = set_par->sp.zres;
                    zreal = set_par->sp.dz*zres;

                    min = gwy_data_field_get_min(gx.dfield);
                    max = gwy_data_field_get_max(gx.dfield);

                    if (gx.i != -1) {
                        set_mat->gwydd_ymin[set_mat->gwydds->len - 1] = (-gx.yoffset) / set_par->sp.dy;
                        set_mat->gwydd_zmin[set_mat->gwydds->len - 1] = (-gx.zoffset) / set_par->sp.dz;
                        set_mat->gwydd_ymax[set_mat->gwydds->len - 1] = ((gdouble)yres*set_par->sp.dy - gx.yoffset) / set_par->sp.dy;
                        set_mat->gwydd_zmax[set_mat->gwydds->len - 1] = ((gdouble)zres*set_par->sp.dz - gx.zoffset) / set_par->sp.dz;
                        set_mat->gwydd_xmin[set_mat->gwydds->len - 1] = gx.i + (gdouble)xres / xreal*min;
                        set_mat->gwydd_xmax[set_mat->gwydds->len - 1] = gx.i + (gdouble)xres / xreal*max;
                    }
                    if (gx.j != -1) {
                        set_mat->gwydd_xmin[set_mat->gwydds->len - 1] = (-gx.xoffset) / set_par->sp.dx;
                        set_mat->gwydd_zmin[set_mat->gwydds->len - 1] = (-gx.zoffset) / set_par->sp.dz;
                        set_mat->gwydd_xmax[set_mat->gwydds->len - 1] = ((gdouble)xres*set_par->sp.dx - gx.xoffset) / set_par->sp.dx;
                        set_mat->gwydd_zmax[set_mat->gwydds->len - 1] = ((gdouble)zres*set_par->sp.dz - gx.zoffset) / set_par->sp.dz;
                        set_mat->gwydd_ymin[set_mat->gwydds->len - 1] = gx.j + (gdouble)yres / yreal*min;
                        set_mat->gwydd_ymax[set_mat->gwydds->len - 1] = gx.j + (gdouble)yres / yreal*max;
                    }
                    if (gx.k != -1) {
                        set_mat->gwydd_xmin[set_mat->gwydds->len - 1] = (-gx.xoffset) / set_par->sp.dx;
                        set_mat->gwydd_ymin[set_mat->gwydds->len - 1] = (-gx.yoffset) / set_par->sp.dy;
                        set_mat->gwydd_xmax[set_mat->gwydds->len - 1] = ((gdouble)xres*set_par->sp.dx - gx.xoffset) / set_par->sp.dx;
                        set_mat->gwydd_ymax[set_mat->gwydds->len - 1] = ((gdouble)yres*set_par->sp.dy - gx.yoffset) / set_par->sp.dy;
                        set_mat->gwydd_zmin[set_mat->gwydds->len - 1] = gx.k + (gdouble)zres / zreal*min;
                        set_mat->gwydd_zmax[set_mat->gwydds->len - 1] = gx.k + (gdouble)zres / zreal*max;
                    }
                    set_mat->gwydd_nvx[set_mat->gwydds->len - 1] = 0;

                    for (i = 0; i < xres; i++) { //x
                        for (j = 0; j < yres; j++) { //y
                            for (k = 0; k < zres; k++) { //z
                                mask = 0;

                                if (gx.i != -1 && ggwy_data_field_inside(gx.dfield, (gint)(j - gx.yoffset / set_par->sp.dy), (gint)(k - gx.zoffset / set_par->sp.dz))) {
                                    val = gx.i + (gdouble)xres / xreal*(gwy_data_field_get_dval_real(gx.dfield, (gdouble)j*set_par->sp.dy - gx.yoffset, (gdouble)k*set_par->sp.dz - gx.zoffset, GWY_INTERPOLATION_BILINEAR) - gx.xoffset);
                                    if (gx.mask == 1 && gx.mfield)
                                        mask = gwy_data_field_get_dval_real(gx.mfield, (gdouble)j*set_par->sp.dy - gx.xoffset, (gdouble)k*set_par->sp.dz - gx.yoffset, GWY_INTERPOLATION_ROUND);
                                    else if (gx.mask == -1 && gx.mfield)
                                        mask = fabs(gwy_data_field_get_dval_real(gx.mfield, (gdouble)j*set_par->sp.dy - gx.xoffset, (gdouble)k*set_par->sp.dz - gx.yoffset, GWY_INTERPOLATION_ROUND) - 1);

                                    if (gx.depth > 0 && mask == 0 && i > val && i < (val + gx.depth)) {
                                        buf->data[i][j][k] = (gdouble)i / (gdouble)xres;
                                        if (set_mat->gwydd_xmin[set_mat->gwydds->len - 1] > i)
                                            set_mat->gwydd_xmin[set_mat->gwydds->len - 1] = i;
                                        if (set_mat->gwydd_xmax[set_mat->gwydds->len - 1] < i)
                                            set_mat->gwydd_xmax[set_mat->gwydds->len - 1] = i;
                                        set_mat->gwydd_nvx[set_mat->gwydds->len - 1]++;
                                    }
                                }

                                if (gx.j != -1 && ggwy_data_field_inside(gx.dfield, (gint)(i - gx.xoffset / set_par->sp.dx), (gint)(k - gx.zoffset / set_par->sp.dz))) {
                                    val = gx.j + (gdouble)yres / yreal*(gwy_data_field_get_dval_real(gx.dfield, (gdouble)i*set_par->sp.dx - gx.xoffset, (gdouble)k*set_par->sp.dz - gx.zoffset, GWY_INTERPOLATION_BILINEAR) - gx.yoffset);
                                    if (gx.mask == 1 && gx.mfield)
                                        mask = gwy_data_field_get_dval_real(gx.mfield, (gdouble)i*set_par->sp.dx - gx.xoffset, (gdouble)k*set_par->sp.dz - gx.yoffset, GWY_INTERPOLATION_ROUND);
                                    else if (gx.mask == -1 && gx.mfield)
                                        mask = fabs(gwy_data_field_get_dval_real(gx.mfield, (gdouble)i*set_par->sp.dx - gx.xoffset, (gdouble)k*set_par->sp.dz - gx.yoffset, GWY_INTERPOLATION_ROUND) - 1);


                                    if (gx.depth > 0 && mask == 0 && j > val && j < (val + gx.depth)) {
                                        buf->data[i][j][k] = (gdouble)j / (gdouble)yres;
                                        if (set_mat->gwydd_ymin[set_mat->gwydds->len - 1] > k)
                                            set_mat->gwydd_ymin[set_mat->gwydds->len - 1] = j;
                                        if (set_mat->gwydd_ymax[set_mat->gwydds->len - 1] < k)
                                            set_mat->gwydd_ymax[set_mat->gwydds->len - 1] = j;
                                        set_mat->gwydd_nvx[set_mat->gwydds->len - 1]++;
                                    }
                                }

                                if (gx.k != -1 && ggwy_data_field_inside(gx.dfield, (gint)(i - gx.xoffset / set_par->sp.dx), (gint)(j - gx.yoffset / set_par->sp.dy))) {
                                    val = gx.k + (gdouble)zres / zreal*(gwy_data_field_get_dval_real(gx.dfield, (gdouble)i*set_par->sp.dx - gx.xoffset, (gdouble)j*set_par->sp.dy - gx.yoffset, GWY_INTERPOLATION_BILINEAR) - gx.zoffset);
                                    if (gx.mask == 1 && gx.mfield)
                                        mask = gwy_data_field_get_dval_real(gx.mfield, (gdouble)i*set_par->sp.dx - gx.xoffset, (gdouble)j*set_par->sp.dy - gx.yoffset, GWY_INTERPOLATION_ROUND);
                                    else if (gx.mask == -1 && gx.mfield)
                                        mask = fabs(gwy_data_field_get_dval_real(gx.mfield, (gdouble)i*set_par->sp.dx - gx.xoffset, (gdouble)j*set_par->sp.dy - gx.yoffset, GWY_INTERPOLATION_ROUND) - 1);

                                    //                                                printf("%d %g %g\n", k, val, (val+gx.depth));
                                    if (gx.depth > 0 && k > val && k < (val + gx.depth)) {
                                        buf->data[i][j][k] = (gdouble)k / (gdouble)zres;
                                        if (set_mat->gwydd_zmin[set_mat->gwydds->len - 1] > k) set_mat->gwydd_zmin[set_mat->gwydds->len - 1] = k;
                                        if (set_mat->gwydd_zmax[set_mat->gwydds->len - 1] < k) set_mat->gwydd_zmax[set_mat->gwydds->len - 1] = k;
                                        set_mat->gwydd_nvx[set_mat->gwydds->len - 1]++;
                                    }
                                }
                            }
                        }
                    }
                    /*simplify*/

                    for (i = 1; i < xres - 1; i++) { //x
                        for (j = 1; j < yres - 1; j++) { //y
                            for (k = 1; k < zres - 1; k++) { //z
                                if (((gx.i != -1) && (((j / 2) % 10) == 0 || ((k / 2) % 10) == 0))) {
                                    if (buf->data[i][j][k] > 0 &&
                                        !(buf->data[i - 1][j][k] > 0 && buf->data[i + 1][j][k] > 0
                                            && buf->data[i][j - 1][k] > 0 && buf->data[i][j + 1][k] > 0
                                            && buf->data[i][j][k - 1] > 0 && buf->data[i][j][k + 1] > 0))
                                        set_mat->gwyddata->data[i][j][k] = buf->data[i][j][k];
                                }
                                else if (((gx.j != -1) && (((i / 2) % 10) == 0 || ((k / 2) % 10) == 0))) {
                                    if (buf->data[i][j][k] > 0 &&
                                        !(buf->data[i - 1][j][k] > 0 && buf->data[i + 1][j][k] > 0
                                            && buf->data[i][j - 1][k] > 0 && buf->data[i][j + 1][k] > 0
                                            && buf->data[i][j][k - 1] > 0 && buf->data[i][j][k + 1] > 0))
                                        set_mat->gwyddata->data[i][j][k] = buf->data[i][j][k];
                                }
                                else if (((gx.k != -1) && (((i / 2) % 10) == 0 || ((j / 2) % 10) == 0))) {
                                    if (buf->data[i][j][k] > 0 &&
                                        !(buf->data[i - 1][j][k] > 0 && buf->data[i + 1][j][k] > 0
                                            && buf->data[i][j - 1][k] > 0 && buf->data[i][j + 1][k] > 0
                                            && buf->data[i][j][k - 1] > 0 && buf->data[i][j][k + 1] > 0))
                                        set_mat->gwyddata->data[i][j][k] = buf->data[i][j][k];
                                }
                            }
                        }
                    }
                    sv_dcube_free(buf);
                }
                break;


            default:
                fprintf(stderr, "Error: unknown object type in material file (%d)\n", type);                
                goto exit;
            }
        }

        fclose(fr);
        fr = NULL;
    }

    return 1;

exit:
    fclose(fr);
    fr = NULL;

    return 0;
} /* parse_settings_mat */

/*write all the settings to a file*/
int
write_settings(gchar *filename, SvSet *set)
{
    FILE *fw;

    fw = fopen(filename, "w");

    //    fprintf(fw, "POOL\n%d %d %d %g %g %g\n\n", );

    fclose(fw);

    return 0;
}


gboolean
check_cube_inside(GString *message, gchar *description, gchar *box, gint x0, gint y0, gint z0, gint x1, gint y1, gint z1, gint i0, gint j0, gint k0, gint i1, gint j1, gint k1)
{
    gboolean isok = TRUE;


    if (i0 < x0 || i0 >= x1) {
        g_string_append_printf(message, "<b>Error:</b> %s x0 boundary is outside of %s.\n", description, box);
        isok = FALSE;
    }
    if (i0 == x0 || i0 == (x1 - 1)) {
        g_string_append_printf(message, "<b>Error:</b> %s x0 boundary is at boundary of %s.\n", description, box);
        isok = FALSE;
    }
    if (i1 < x0 || i1 >= x1) {
        g_string_append_printf(message, "<b>Error:</b> %s x1 boundary is outside of %s.\n", description, box);
        isok = FALSE;
    }
    if (i1 == x0 || i1 == (x1 - 1)) {
        g_string_append_printf(message, "<b>Error:</b> %s x1 boundary is at boundary of %s.\n", description, box);
        isok = FALSE;
    }
    if (i0 >= i1) {
        g_string_append_printf(message, "<b>Error:</b> %s size is zero or negative in x direction\n", description);
        isok = FALSE;
    }

    if (j0 < y0 || j0 >= y1) {
        g_string_append_printf(message, "<b>Error:</b> %s y0 boundary is outside of %s.\n", description, box);
        isok = FALSE;
    }
    if (j0 == y0 || j0 == (y1 - 1)) {
        g_string_append_printf(message, "<b>Error:</b> %s y0 boundary is at boundary of %s.\n", description, box);
        isok = FALSE;
    }
    if (j1 < y0 || j1 >= y1) {
        g_string_append_printf(message, "<b>Error:</b> %s y1 boundary is outside of %s.\n", description, box);
        isok = FALSE;
    }
    if (j1 == y0 || j1 == (y1 - 1)) {
        g_string_append_printf(message, "<b>Error:</b> %s y1 boundary is at boundary of %s.\n", description, box);
        isok = FALSE;
    }
    if (j0 >= j1) {
        g_string_append_printf(message, "<b>Error:</b> %s size is zero or negative in x direction\n", description);
        isok = FALSE;
    }

    if (k0 < z0 || k0 >= z1) {
        g_string_append_printf(message, "<b>Error:</b> %s z0 boundary is outside of %s.\n", description, box);
        isok = FALSE;
    }
    if (k0 == z0 || k0 == (z1 - 1)) {
        g_string_append_printf(message, "<b>Error:</b> %s z0 boundary is at boundary of %s.\n", description, box);
        isok = FALSE;
    }
    if (k1 < z0 || k1 >= z1) {
        g_string_append_printf(message, "<b>Error:</b> %s z1 boundary is outside of %s.\n", description, box);
        isok = FALSE;
    }
    if (k1 == z0 || k1 == (z1 - 1)) {
        g_string_append_printf(message, "<b>Error:</b> %s z1 boundary is at boundary of %s.\n", description, box);
        isok = FALSE;
    }
    if (k0 >= k1) {
        g_string_append_printf(message, "<b>Error:</b> %s size is zero or negative in x direction\n", description);
        isok = FALSE;
    }

    return isok;
}


gboolean
check_point_inside(GString *message, gchar *description, gchar *box, gint x0, gint y0, gint z0, gint x1, gint y1, gint z1, gint i, gint j, gint k)
{
    gboolean isok = TRUE;

    if (i < x0 || i >= x1) {
        g_string_append_printf(message, "<b>Error:</b> %s location is outside of %s.\n", description, box);
        isok = FALSE;
    }
    if (i == x0 || i == (x1 - 1)) {
        g_string_append_printf(message, "<b>Warning:</b> %s location is at boundary of %s.\n", description, box);
        isok = FALSE;
    }

    if (j < y0 || j >= y1) {
        g_string_append_printf(message, "<b>Error:</b> %s location is outside of %s.\n", description, box);
        isok = FALSE;
    }
    if (j == y0 || j == (y1 - 1)) {
        g_string_append_printf(message, "<b>Warning:</b> %s location is at boundary of %s.\n", description, box);
        isok = FALSE;
    }

    if (k < z0 || k >= z1) {
        g_string_append_printf(message, "<b>Error:</b> %s location is outside of %s.\n", description, box);
        isok = FALSE;
    }
    if (k == z0 || i == (z1 - 1)) {
        g_string_append_printf(message, "<b>Warning:</b> %s location is at boundary of %s.\n", description, box);
        isok = FALSE;
    }

    return isok;
}


/* check all the settings consistency */
GString *
sv_check(SvSet *set, GString *message, gboolean *ok)
{
    gint i;
    gboolean isok = TRUE;

    for (i = 0; i < set->ss.npnts; i++) {

        isok &= check_point_inside(message, "Point source:", "the computational domain",
                                   0, 0, 0, set->sp.xres, set->sp.yres, set->sp.zres, set->ss.pnts[i].point_origin_position_i, set->ss.pnts[i].point_origin_position_j, set->ss.pnts[i].point_origin_position_k);
    }

    if (set->ss.tsf.valid) {

        isok &= check_cube_inside(message, "TSF source:", "the computational domain", 0, 0, 0,
                                  set->sp.xres, set->sp.yres, set->sp.zres,
                                  set->ss.tsf.box_i0, set->ss.tsf.box_j0, set->ss.tsf.box_k0, set->ss.tsf.box_in, set->ss.tsf.box_jn, set->ss.tsf.box_kn);
    }

    if (set->ss.ltsf.valid) {
        isok &= check_cube_inside(message, "LTSF source:", "the computational domain", 0, 0, 0,
                                  set->sp.xres, set->sp.yres, set->sp.zres,
                                  set->ss.ltsf.box_i0, set->ss.ltsf.box_j0, set->ss.ltsf.box_k0, set->ss.ltsf.box_in, set->ss.ltsf.box_jn, set->ss.ltsf.box_kn);
    }

    if (set->ss.tsff.valid) {
        isok &= check_cube_inside(message, "TSFF source:", "the computational domain", 0, 0, 0,
                                  set->sp.xres, set->sp.yres, set->sp.zres,
                                  set->ss.tsff.box_i0, set->ss.tsff.box_j0, set->ss.tsff.box_k0, set->ss.tsff.box_in, set->ss.tsff.box_jn, set->ss.tsff.box_kn);
    }

    if (set->ss.ltsff.valid) {
        isok &= check_cube_inside(message, "LTSFF source:", "the computational domain", 0, 0, 0,
                                  set->sp.xres, set->sp.yres, set->sp.zres,
                                  set->ss.ltsff.box_i0, set->ss.ltsff.box_j0, set->ss.ltsff.box_k0, set->ss.ltsff.box_in, set->ss.ltsff.box_jn, set->ss.ltsff.box_kn);
    }

    for (i = 0; i < set->so.npnts; i++) {
        isok &= check_point_inside(message, "Point output:", "the computational domain",
                                   0, 0, 0, set->sp.xres, set->sp.yres, set->sp.zres, set->so.pnts[i].i, set->so.pnts[i].j, set->so.pnts[i].k);
    }


    if (isok)
        message = g_string_append(message, "OK.");
    else
        message = g_string_append(message, "Failed.");


    *ok = isok;

    return message;
}

void alloc_set_par_mat(SvSet *set_par, SvSetMat *set_mat, gboolean init_data)
{
    if (init_data) {
        set_par->sp.xres = set_par->sp.yres = set_par->sp.zres = COMP_DOMAIN_SIZE;
        set_par->sp.dx = set_par->sp.dy = set_par->sp.dz = COMP_DOMAIN_SPACING;
        set_par->sc.nsteps = BASIC_STEPS;
        set_par->sc.verbose = BASIC_VERBOSE;
        set_par->sc.nthreads = BASIC_NTHREADS;
        set_par->sc.usegpu = BASIC_USEGPU;
        memset(set_par->sc.ugpu, 0, sizeof(set_par->sc.ugpu));
    }

    alloc_set_mat(set_mat);
}

void alloc_set_mat(SvSetMat *set_mat)
{
    if (set_mat->spheres)
        g_array_free(set_mat->spheres, TRUE);
    if (set_mat->cones)
        g_array_free(set_mat->cones, TRUE);
    if (set_mat->rcones)
        g_array_free(set_mat->rcones, TRUE);
    if (set_mat->cylinders)
        g_array_free(set_mat->cylinders, TRUE);
    if (set_mat->voxels)
        g_array_free(set_mat->voxels, TRUE);
    if (set_mat->tetrahedrons)
        g_array_free(set_mat->tetrahedrons, TRUE);
    if (set_mat->gwydds)
        g_array_free(set_mat->gwydds, TRUE);

    set_mat->spheres = g_array_new(FALSE, FALSE, sizeof(SvSphere));
    set_mat->cones = g_array_new(FALSE, FALSE, sizeof(SvCone));
    set_mat->rcones = g_array_new(FALSE, FALSE, sizeof(SvRCone));
    set_mat->cylinders = g_array_new(FALSE, FALSE, sizeof(SvCylinder));
    set_mat->voxels = g_array_new(FALSE, FALSE, sizeof(SvVoxel));
    set_mat->tetrahedrons = g_array_new(FALSE, FALSE, sizeof(SvTetrahedron));
    set_mat->gwydds = g_array_new(FALSE, FALSE, sizeof(SvGwydd));
}


/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
