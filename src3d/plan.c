
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


 /*  plan.c :
  *  creating plan for computation (e.g. which materials to use), loading
  *  data files and allocating data
  */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <glib.h>
#include <libprocess/gwyprocess.h>
#include <libgwyddion/gwyexpr.h>
#include "global.h"
#include "plan.h"
#include "pool.h"
#include "settings.h"
#include "constants.h"
#include "output.h"
#include <math.h>
#include <glib.h>
#include <omp.h>
#include "modifiers.h"

#if defined (_WIN32) && defined (_MSC_VER)
#include "../config_win.h"
#else
#include "../config.h"
#include <sys/time.h>
#endif

typedef struct {
    gint i;
    gint j;
    gint k;
} SmCoords;

typedef struct {
    SmCoords coords;
    SvMatProp material[3][3][3]; // pole ukazatelu
    gint cells_written;
    gboolean was_mixed;
} Submesh;

typedef struct {
    SvPool *mp;
    SvSet *set;
} MpSetData;

/*
 gint type; gdouble epsilon; mu; sigma; sigast; drude_omega_p; drude_nu; cp3_a[3]; cp3_phi[3]; cp3_omega[3]; cp3_gamma[3]; ade_a0; ade_a1; ade_a2;
 ade_bp0[2]; ade_bp1[2]; ade_bp2[2]; ade_bp3[2]; ade_bp4[2]; ade_c0; ade_c1; ade_c2; ade_c3; plrc_d_chi; plrc_d_xi; plrc_d_dchi; plrc_d_dxi;
 plrc_p_chi[2]; plrc_p_xi[2]; plrc_p_dchir[2]; plrc_p_dxir[2]; plrc_p_dchii[2]; plrc_p_dxii[2]; gint pos;
*/
const SvMatProp mat_vacuum_linear = {SV_MAT_LINEAR, 1, 1, 0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, 0, 0, 0,
                     {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, 0, 0, 0, 0, 0, 0, 0, 0,
                     {0, 0}, {0, 0}, {0, 0}, {0,0}, {0, 0}, {0, 0}, -1};
const SvMatProp mat_vacuum_lintab = {SV_MAT_LINTAB, 1, 1, 0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, 0, 0, 0,
                     {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, 0, 0, 0, 0, 0, 0, 0, 0,
                     {0, 0}, {0, 0}, {0, 0}, {0,0}, {0, 0}, {0, 0}, 0};

void
gwy_app_data_browser_add(GwyContainer *data);

static gint
get_n_cores()
{
#ifdef G_OS_WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

gboolean
scan_point(FILE *fr, GLfloat *px, GLfloat *py, GLfloat *pz)
{
    GLfloat x, y, z;
    if (fscanf(fr, "%f", &x) != 1)
        return 1;  //TODO this should be locale dependent
    if (fscanf(fr, "%f", &y) != 1)
        return 1;
    if (fscanf(fr, "%f", &z) != 1)
        return 1;
    *px = x;
    *py = y;
    *pz = z;
    // printf("scan point finished\n");
    return 0;
}

gboolean
scan_radius(FILE *fr, GLfloat *pradius)
{
    GLfloat radius;
    if (fscanf(fr, "%f", &radius) != 1)
        return 1; //TODO this should be locale dependent
    *pradius = radius;
    // printf("scan radius finished\n");
    return 0;
}

gboolean
scan_int(FILE *fr, gint *pint)
{
    gint val;
    if (fscanf(fr, "%d", &val) != 1)
        return 1; //TODO this should be locale dependent
    *pint = val;
    // printf("scan integer finished\n");
    return 0;
}


/*check whether there is electric or magnetic material in material file*/
SvMatMode
test_material_file(gchar *filename)
{
    gboolean is_e = 0;
    gboolean is_m = 0;
    int i, j, k, xres, yres, zres;
    gfloat value;
    FILE *fr = fopen(filename, "rb");

    if (fr == NULL)
        return SV_MATMODE_NONE;

    if (fread(&xres, sizeof(gint), 1, fr) != 1) {
        fclose(fr);
        return SV_MATMODE_NONE;
    }
    if (fread(&yres, sizeof(gint), 1, fr) != 1) {
        fclose(fr);
        return SV_MATMODE_NONE;
    }
    if (fread(&zres, sizeof(gint), 1, fr) != 1) {
        fclose(fr);
        return SV_MATMODE_NONE;
    }
    for (i = 0; i < xres; i++) {
        for (j = 0; j < yres; j++) {
            for (k = 0; k < zres; k++) {
                if (fread(&value, sizeof(gfloat), 1, fr) != 1) {
                    fclose(fr);
                    return SV_MATMODE_NONE;
                }
                if (value != 1) {
                    is_e = 1;
                    break;
                }
            }
            if (is_e) break;
        }
        if (is_e) break;
    }
    if (!is_e) {
        for (i = 0; i < xres; i++) {
            for (j = 0; j < yres; j++) {
                for (k = 0; k < zres; k++) {
                    if (fread(&value, sizeof(gfloat), 1, fr) != 1) {
                        fclose(fr);
                        return SV_MATMODE_NONE;
                    }
                    if (value != 0) {
                        is_e = 1;
                        break;
                    }
                }
                if (is_e)
                    break;
            }
            if (is_e)
                break;
        }
    }

    for (i = 0; i < xres; i++) {
        for (j = 0; j < yres; j++) {
            for (k = 0; k < zres; k++) {
                if (fread(&value, sizeof(gfloat), 1, fr) != 1) {
                    fclose(fr);
                    return SV_MATMODE_NONE;
                }
                if (value != 1) {
                    is_m = 1;
                    break;
                }
            }
            if (is_m)
                break;
        }
        if (is_m)
            break;
    }
    if (!is_m) {
        for (i = 0; i < xres; i++) {
            for (j = 0; j < yres; j++) {
                for (k = 0; k < zres; k++) {
                    if (fread(&value, sizeof(gfloat), 1, fr) != 1) {
                        fclose(fr);
                        return SV_MATMODE_NONE;
                    }
                    if (value != 0) {
                        is_m = 1;
                        break;
                    }
                }
                if (is_m)
                    break;
            }
            if (is_m)
                break;
        }
    }
    fclose(fr);

    if ((is_e + is_m) == 0)
        return SV_MATMODE_NONE;
    if (is_e == 0 && is_m != 0)
        return SV_MATMODE_MAGNETIC;
    if (is_e != 0 && is_m == 0)
        return SV_MATMODE_ELECTRIC;

    return SV_MATMODE_FULL;
}

gboolean
test_material(FILE *fr, gboolean *is_e, gboolean *is_m, gboolean *is_tab)
{
    float epsilon, mu, sigma, sigast, omegap, nu, a, omega, phi, gamma;
    int mattype;
    char string[256];


    if (fscanf(fr, "%d", &mattype) != 1)
        return 1;

    if (mattype == SV_MAT_LINEAR || mattype == SV_MAT_LINTAB) {
        if (fscanf(fr, "%f", &epsilon) != 1)
            return 1;
        if (fscanf(fr, "%f", &mu) != 1)
            return 1;
        if (fscanf(fr, "%f", &sigma) != 1)
            return 1;
        if (fscanf(fr, "%f", &sigast) != 1)
            return 1;
        if (mattype == SV_MAT_LINEAR) { //for lintab don't alloc material fields
            if (epsilon != 1.0 || sigma > 0) *is_e = 1;
            if (mu != 1.0 || sigast > 0) *is_m = 1;
        }
    } else if (mattype == SV_MAT_DRUDE) {
        if (fscanf(fr, "%f", &epsilon) != 1)
            return 1;
        if (fscanf(fr, "%f", &omegap) != 1)
            return 1;
        if (fscanf(fr, "%f", &nu) != 1)
            return 1;
        *is_e = 1;
        /*do nothing, if only drude is here we dont need material fields*/
    } else if (mattype == SV_MAT_CP) {
        if (fscanf(fr, "%f", &epsilon) != 1)
            return 1;
        if (fscanf(fr, "%f", &omegap) != 1)
            return 1;
        if (fscanf(fr, "%f", &nu) != 1)
            return 1;
        if (fscanf(fr, "%f", &a) != 1)
            return 1;
        if (fscanf(fr, "%f", &phi) != 1)
            return 1;
        if (fscanf(fr, "%f", &omega) != 1)
            return 1;
        if (fscanf(fr, "%f", &gamma) != 1)
            return 1;
        if (fscanf(fr, "%f", &a) != 1)
            return 1;
        if (fscanf(fr, "%f", &phi) != 1)
            return 1;
        if (fscanf(fr, "%f", &omega) != 1)
            return 1;
        if (fscanf(fr, "%f", &gamma) != 1)
            return 1;
        *is_e = 1;
    } else if (mattype == SV_MAT_ADE || mattype == SV_MAT_PLRC) {
        if (fscanf(fr, "%f", &epsilon) != 1)
            return 1;
        if (fscanf(fr, "%f", &sigma) != 1)
            return 1;
        if (fscanf(fr, "%f", &omegap) != 1)
            return 1;
        if (fscanf(fr, "%f", &nu) != 1)
            return 1;
        if (fscanf(fr, "%f", &a) != 1)
            return 1;
        if (fscanf(fr, "%f", &phi) != 1)
            return 1;
        if (fscanf(fr, "%f", &omega) != 1)
            return 1;
        if (fscanf(fr, "%f", &gamma) != 1)
            return 1;
        if (fscanf(fr, "%f", &a) != 1)
            return 1;
        if (fscanf(fr, "%f", &phi) != 1)
            return 1;
        if (fscanf(fr, "%f", &omega) != 1)
            return 1;
        if (fscanf(fr, "%f", &gamma) != 1)
            return 1;
        *is_e = 1;
    } else if (mattype == SV_MAT_PEC) {
        /*do nothing*/
    } else if (mattype == SV_MAT_CP3) {
        if (fscanf(fr, "%f", &epsilon) != 1)
            return 1;
        if (fscanf(fr, "%f", &sigma) != 1)
            return 1;
        if (fscanf(fr, "%f", &a) != 1)
            return 1;
        if (fscanf(fr, "%f", &phi) != 1)
            return 1;
        if (fscanf(fr, "%f", &omega) != 1)
            return 1;
        if (fscanf(fr, "%f", &gamma) != 1)
            return 1;
        if (fscanf(fr, "%f", &a) != 1)
            return 1;
        if (fscanf(fr, "%f", &phi) != 1)
            return 1;
        if (fscanf(fr, "%f", &omega) != 1)
            return 1;
        if (fscanf(fr, "%f", &gamma) != 1)
            return 1;
        if (fscanf(fr, "%f", &a) != 1)
            return 1;
        if (fscanf(fr, "%f", &phi) != 1)
            return 1;
        if (fscanf(fr, "%f", &omega) != 1)
            return 1;
        if (fscanf(fr, "%f", &gamma) != 1)
            return 1;
        *is_e = 1;
    } else if (mattype == SV_MAT_DATABASE) {
        fscanf(fr, "%s", string);
        *is_tab = 1;
    } else {
        fprintf(stderr, "Error: test: Unknown material type? (%d)\n", mattype);
        return 1;
    }

    return 0;
}

gboolean
subst_mat(gchar *string, gdouble *epsilon, gdouble *mu, gdouble *sigma, gdouble *sigast, gdouble wfrom, gdouble wcenter, gdouble wto, gboolean localfiles)
{
    gdouble wav;
    char stype[100];
    char *nl;
    float swfrom, swto, par1, par2, par3, par4, par5, par6, par7, wf, nf, kf;
    FILE *fs = NULL;
    gchar *path = NULL;
    gdouble n, k;
    char *description = (char *)g_malloc(256 * sizeof(char));
    gint i, ndata, dir = 0;
    gdouble mult = 1e-6, wprev = 0, wnext = 0, nprev = 0, nnext = 0, kprev = 0, knext = 0;
    gdouble *wfs, *nfs, *kfs;

    if (!string)
        return 1;
    printf("mat:%s:\n", string);

    if (localfiles) {
        fs = fopen(string, "r"); 
    }

    if (!fs) {
        path = get_spectra_path(string);
        fs = fopen(path, "r");
    }

    if (fs == NULL) {
        printf("no database file for this material (%s)!\n", string);
        if (path) fprintf(stderr, "Error: cannot open data from optical database: %s\n", path);
        return 1;
    } else {
        printf("file %s successfully opened\n", path);
        fscanf(fs, "%s", stype);
        fscanf(fs, "%f", &swfrom);
        fscanf(fs, "%f", &swto);
        description = fgets(description, 256, fs);
        description = fgets(description, 256, fs);
        if (!description)
            return 1;
        //getline(&description, &len, fs);
        //getline(&description, &len, fs);
        nl = strrchr(description, '\r');
        if (nl)
            *nl = '\0';
        nl = strrchr(description, '\n');
        if (nl)
            *nl = '\0';

        swfrom *= (float)1.0e-6; //spectra ranges are in microns
        swto *= (float)1.0e-6;

        printf("Spectrum type %s, range %g %g micrometers (%s) used\n", stype, swfrom, swto, description);
        if (swfrom > wfrom || swto < wto)
            fprintf(stderr, "Warning: your source wavelength (%g to %g, center at %g) does not match database spectrum range (%g, %g)\n", wfrom, wcenter, wto, swfrom, swto);

        printf("trying to match spectrum type '%s'\n", stype);
        if (strcmp(stype, "sellmeier") == 0) {
            fscanf(fs, "%g", &par1);
            fscanf(fs, "%g", &par2);
            fscanf(fs, "%g", &par3);
            fscanf(fs, "%g", &par4);
            fscanf(fs, "%g", &par5);
            fscanf(fs, "%g", &par6);
            fscanf(fs, "%g", &par7);

            wav = wcenter*1e6; //wavelength in um
            n = sqrt((double)(par1 + par2*wav*wav / (wav*wav - par3) + par4*wav*wav / (wav*wav - par5) + par6*wav*wav / (wav*wav - par7)));

            *epsilon = n*n;
            *sigma = 0;
            *mu = 1;
            *sigast = 0;
        } else if (strcmp(stype, "cauchy") == 0) {
            fscanf(fs, "%g", &par1);
            fscanf(fs, "%g", &par2);
            fscanf(fs, "%g", &par3);

            wav = wcenter*1e6; //wavelength in um
            n = par1 + par2 / wav / wav + par3 / wav / wav / wav / wav;

            *epsilon = n*n;
            *sigma = 0;
            *mu = 1;
            *sigast = 0;
        } else if (strcmp(stype, "kasarova") == 0) {
            fscanf(fs, "%g", &par1);
            fscanf(fs, "%g", &par2);
            fscanf(fs, "%g", &par3);
            fscanf(fs, "%g", &par4);
            fscanf(fs, "%g", &par5);
            fscanf(fs, "%g", &par6);

            wav = wcenter*1e6; //wavelength in um
            n = sqrt((double)(par1 + par2*wav*wav + par3 / wav / wav + par4 / wav / wav / wav / wav + par5 / wav / wav / wav / wav / wav / wav + par6 / wav / wav / wav / wav / wav / wav / wav / wav));

            *epsilon = n*n;
            *sigma = 0;
            *mu = 1;
            *sigast = 0;
        } else if (strcmp(stype, "nk") == 0) {
            fscanf(fs, "%d", &ndata);

            wfs = (gdouble *)g_malloc(ndata * sizeof(gdouble));
            nfs = (gdouble *)g_malloc(ndata * sizeof(gdouble));
            kfs = (gdouble *)g_malloc(ndata * sizeof(gdouble));

            for (i = 0; i < (ndata); i++) {
                fscanf(fs, "%g", &wf);
                fscanf(fs, "%g", &nf);
                fscanf(fs, "%g", &kf);

                if (i == 0) {
                    if (fabs(wf*1e-6 - swfrom) < 1e-10) {
                        dir = 1;
                        mult = 1e-6;
                    } else if (fabs(wf*1e-6 - swto) < 1e-10) {
                        dir = 0;
                        mult = 1e-6;
                    } else if (fabs(wf*1e-9 - swfrom) < 1e-10) {
                        dir = 1;
                        mult = 1e-9;
                    } else if (fabs(wf*1e-9 - swto) < 1e-10) {
                        dir = 0;
                        mult = 1e-9;
                    } else {
                        fprintf(stderr, "Database entry error: wavelength range (%g %g) does not match the spectrum value %g\n", swfrom, swto, wf);
                        return 1;
                    }
                }

                wfs[i] = wf*mult;
                nfs[i] = nf;
                kfs[i] = kf;
            }

            if (dir == 1 && wcenter < wfs[0]) {
                wprev = wfs[0];
                nprev = nfs[0];
                kprev = kfs[0];
                wnext = wfs[1];
                nnext = nfs[1];
                knext = kfs[1];
                fprintf(stderr, "Warning: Spectral data extrapolated at lower side of spectrum\n");
            } else if (dir == 0 && wcenter < wfs[ndata - 1]) {
                wprev = wfs[ndata - 1];
                nprev = nfs[ndata - 1];
                kprev = kfs[ndata - 1];
                wnext = wfs[ndata - 2];
                nnext = nfs[ndata - 2];
                knext = kfs[ndata - 2];
                fprintf(stderr, "Warning: Spectral data extrapolated at lower side of spectrum\n");
            } else {
                for (i = 0; i < (ndata - 1); i++) {
                    if (dir == 1 && wfs[i] <= wcenter && wfs[i + 1] >= wcenter) { //data in table going upwards
                        wprev = wfs[i];
                        nprev = nfs[i];
                        kprev = kfs[i];
                        wnext = wfs[i + 1];
                        nnext = nfs[i + 1];
                        knext = kfs[i + 1];
                        break;
                    } else if (dir == 0 && wfs[i] >= wcenter && wfs[i + 1] <= wcenter) { //data going backwards
                        wprev = wfs[i + 1];
                        nprev = nfs[i + 1];
                        kprev = kfs[i + 1];
                        wnext = wfs[i];
                        nnext = nfs[i];
                        knext = kfs[i];
                        break;
                    }
                }
            }
            n = nprev + (nnext - nprev)*(wcenter - wprev) / (wnext - wprev);
            k = kprev + (knext - kprev)*(wcenter - wprev) / (wnext - wprev);

            *epsilon = n*n - k*k;
            *sigma = 4 * G_PI*n*k*LIGHT_SPEED / wcenter*EPSILON_0;
            *mu = 1;
            *sigast = 0;

            g_free(wfs);
            g_free(nfs);
            g_free(kfs);
        } else {
            fprintf(stderr, "Error: unsupported spectrum type (%s)\n", stype);
            return 1;
        }
    }

    fclose(fs);

    return 0;
}

/*check whether there is electric or magnetic material in material script file*/
SvMatMode
test_script_file(gchar *filename)
{
    gint type, ibuf, ret;
    GLfloat buf;
    gchar filebase[256];
    FILE *fr;
    gboolean is_e = 0;
    gboolean is_m = 0;
    gboolean is_tab = 0;


    fr = fopen(filename, "r");
    if (fr == NULL) {
        fprintf(stderr, "Error: cannot open file %s for reading\n", filename);
        return SV_MATMODE_ERROR;
    }

    while ((ret = fscanf(fr, "%d", &type)) != EOF) { //TODO use proper locale handling functions here
        if (!ret) {
            printf("Error: something strange in material input file\n");
            goto exit;
        }

        switch (type) {
            case SV_ENTITY_SPHERE:
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (scan_radius(fr, &buf))
                    goto exit;
                if (test_material(fr, &is_e, &is_m, &is_tab))
                    goto exit;
                break;

            case SV_ENTITY_CYLINDER:
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (scan_radius(fr, &buf))
                    goto exit;
                if (test_material(fr, &is_e, &is_m, &is_tab))
                    goto exit;
                break;

            case SV_ENTITY_CONE:
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (scan_radius(fr, &buf))
                    goto exit;
                if (test_material(fr, &is_e, &is_m, &is_tab))
                    goto exit;
                break;

            case SV_ENTITY_RCONE:
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (scan_radius(fr, &buf))
                    goto exit;
                if (scan_radius(fr, &buf))
                    goto exit;
                if (test_material(fr, &is_e, &is_m, &is_tab))
                    goto exit;
                break;

            case SV_ENTITY_VOXEL:
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (test_material(fr, &is_e, &is_m, &is_tab))
                    goto exit;
                break;

            case SV_ENTITY_TETRAHEDRON:
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (test_material(fr, &is_e, &is_m, &is_tab))
                    goto exit;
                break;

            case SV_ENTITY_TETGEN:
                fscanf(fr, "%s", filebase);             //TODO this should be locale dependent
                if (scan_radius(fr, &buf))
                    goto exit;
                if (scan_radius(fr, &buf))
                    goto exit;
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (test_material(fr, &is_e, &is_m, &is_tab))
                    goto exit;
                break;

            case SV_ENTITY_GWYDDION:
                fscanf(fr, "%s", filebase);             //TODO this should be locale dependent
                if (scan_int(fr, &ibuf))
                    goto exit;
                if (scan_int(fr, &ibuf))
                    goto exit;
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (scan_point(fr, &buf, &buf, &buf))
                    goto exit;
                if (scan_int(fr, &ibuf))
                    goto exit;
                if (test_material(fr, &is_e, &is_m, &is_tab))
                    goto exit;
                break;

            default:
                break;
        }
    }

    fclose(fr);

    if (is_tab)
        return SV_MATMODE_DATABASE;

    if ((is_e + is_m) == 0)
        return SV_MATMODE_NONE;
    if (is_e == 0 && is_m != 0)
        return SV_MATMODE_MAGNETIC;
    if (is_e != 0 && is_m == 0)
        return SV_MATMODE_ELECTRIC;

    return SV_MATMODE_NONE;

exit:
    fclose(fr);
    return SV_MATMODE_ERROR;
}

/*replace the database material with appropriate (e.g. interpolated) values*/
gboolean
adjust_material(FILE *fr, FILE *fw, gdouble wfrom, gdouble wcenter, gdouble wto, gint *nmaterials, gchar **materials, gint pos, gint *overpos, gboolean localfiles)
{
    float epsilon, mu, sigma, sigast, omegap, nu, a, omega, phi, gamma;
    int mattype;
    gdouble wav;
    char stype[100];
    char string[256], *nl;
    float swfrom, swto, par1, par2, par3, par4, par5, par6, par7, wf, nf, kf;
    FILE *fs = NULL;
    gchar *path = NULL;
    gdouble n, k, eps, sig;
    char *description = (char *)g_malloc(256 * sizeof(char));
    gint i, ndata, dir = 0;
    gdouble mult = 1e-6, wprev = 0, wnext = 0, nprev = 0, nnext = 0, kprev = 0, knext = 0;
    gdouble *wfs, *nfs, *kfs;

    if (fscanf(fr, "%d", &mattype) != 1)
        return 1;

    if (mattype == SV_MAT_LINEAR || mattype == SV_MAT_LINTAB) {
        fprintf(fw, "%d ", mattype);
        if (fscanf(fr, "%f", &epsilon) != 1)
            return 1;
        fprintf(fw, "%g ", epsilon);
        if (fscanf(fr, "%f", &mu) != 1)
            return 1;
        fprintf(fw, "%g ", mu);
        if (fscanf(fr, "%f", &sigma) != 1)
            return 1;
        fprintf(fw, "%g ", sigma);
        if (fscanf(fr, "%f", &sigast) != 1)
            return 1;
        fprintf(fw, "%g\n", sigast);
    } else if (mattype == SV_MAT_DRUDE) {
        fprintf(fw, "%d ", SV_MAT_DRUDE);
        if (fscanf(fr, "%f", &epsilon) != 1)
            return 1;
        fprintf(fw, "%g ", epsilon);
        if (fscanf(fr, "%f", &omegap) != 1)
            return 1;
        fprintf(fw, "%g ", omegap);
        if (fscanf(fr, "%f", &nu) != 1)
            return 1;
        fprintf(fw, "%g\n", nu);
    } else if (mattype == SV_MAT_CP) {
        fprintf(fw, "%d ", SV_MAT_CP);
        if (fscanf(fr, "%f", &epsilon) != 1)
            return 1;
        fprintf(fw, "%g ", epsilon);
        if (fscanf(fr, "%f", &omegap) != 1)
            return 1;
        fprintf(fw, "%g ", omegap);
        if (fscanf(fr, "%f", &nu) != 1)
            return 1;
        fprintf(fw, "%g ", nu);
        if (fscanf(fr, "%f", &a) != 1)
            return 1;
        fprintf(fw, "%g ", a);
        if (fscanf(fr, "%f", &phi) != 1)
            return 1;
        fprintf(fw, "%g ", phi);
        if (fscanf(fr, "%f", &omega) != 1)
            return 1;
        fprintf(fw, "%g ", omega);
        if (fscanf(fr, "%f", &gamma) != 1)
            return 1;
        fprintf(fw, "%g ", gamma);
        if (fscanf(fr, "%f", &a) != 1)
            return 1;
        fprintf(fw, "%g ", a);
        if (fscanf(fr, "%f", &phi) != 1)
            return 1;
        fprintf(fw, "%g ", phi);
        if (fscanf(fr, "%f", &omega) != 1)
            return 1;
        fprintf(fw, "%g ", omega);
        if (fscanf(fr, "%f", &gamma) != 1)
            return 1;
        fprintf(fw, "%g\n", gamma);
    } else if (mattype == SV_MAT_ADE || mattype == SV_MAT_PLRC) {
        fprintf(fw, "%d ", SV_MAT_ADE);
        if (fscanf(fr, "%f", &epsilon) != 1)
            return 1;
        fprintf(fw, "%g ", epsilon);
        if (fscanf(fr, "%g", &sigma) != 1)
            return 1;
        fprintf(fw, "%g ", sigma);
        if (fscanf(fr, "%f", &omegap) != 1)
            return 1;
        fprintf(fw, "%g ", omegap);
        if (fscanf(fr, "%f", &nu) != 1)
            return 1;
        fprintf(fw, "%g ", nu);
        if (fscanf(fr, "%f", &a) != 1)
            return 1;
        fprintf(fw, "%g ", a);
        if (fscanf(fr, "%f", &phi) != 1)
            return 1;
        fprintf(fw, "%g ", phi);
        if (fscanf(fr, "%f", &omega) != 1)
            return 1;
        fprintf(fw, "%g ", omega);
        if (fscanf(fr, "%f", &gamma) != 1)
            return 1;
        fprintf(fw, "%g ", gamma);
        if (fscanf(fr, "%f", &a) != 1)
            return 1;
        fprintf(fw, "%g ", a);
        if (fscanf(fr, "%f", &phi) != 1)
            return 1;
        fprintf(fw, "%g ", phi);
        if (fscanf(fr, "%f", &omega) != 1)
            return 1;
        fprintf(fw, "%g ", omega);
        if (fscanf(fr, "%f", &gamma) != 1)
            return 1;
        fprintf(fw, "%g\n", gamma);
    } else if (mattype == SV_MAT_PEC) {
        /*do nothing*/
        fprintf(fw, "%d\n", SV_MAT_PEC);
    } else if (mattype == SV_MAT_CP3) {
        fprintf(fw, "%d ", SV_MAT_CP3);
        if (fscanf(fr, "%g", &epsilon) != 1)
            return 1;
        fprintf(fw, "%g ", epsilon);
        if (fscanf(fr, "%g", &sigma) != 1)
            return 1;
        fprintf(fw, "%g ", sigma);
        if (fscanf(fr, "%f", &a) != 1)
            return 1;
        fprintf(fw, "%g ", a);
        if (fscanf(fr, "%f", &phi) != 1)
            return 1;
        fprintf(fw, "%g ", phi);
        if (fscanf(fr, "%f", &omega) != 1)
            return 1;
        fprintf(fw, "%g ", omega);
        if (fscanf(fr, "%f", &gamma) != 1)
            return 1;
        fprintf(fw, "%g ", gamma);
        if (fscanf(fr, "%f", &a) != 1)
            return 1;
        fprintf(fw, "%g ", a);
        if (fscanf(fr, "%f", &phi) != 1)
            return 1;
        fprintf(fw, "%g ", phi);
        if (fscanf(fr, "%f", &omega) != 1)
            return 1;
        fprintf(fw, "%g ", omega);
        if (fscanf(fr, "%f", &gamma) != 1)
            return 1;
        fprintf(fw, "%g ", gamma);
        if (fscanf(fr, "%f", &a) != 1)
            return 1;
        fprintf(fw, "%g ", a);
        if (fscanf(fr, "%f", &phi) != 1)
            return 1;
        fprintf(fw, "%g ", phi);
        if (fscanf(fr, "%f", &omega) != 1)
            return 1;
        fprintf(fw, "%g ", omega);
        if (fscanf(fr, "%f", &gamma) != 1)
            return 1;
        fprintf(fw, "%g\n", gamma);
    } else if (mattype == SV_MAT_DATABASE) {
        if (fscanf(fr, "%s", string) != 1)
            return 1;
        printf("database entry to search: %s\n", string);

        if (materials != NULL && overpos != NULL) { //for use with XSvit only
            if ((*nmaterials) < 99) {
                materials[*nmaterials] = g_strdup(string);
                overpos[*nmaterials] = pos;
                (*nmaterials) += 1;
            }
        }


       if (localfiles) 
           fs = fopen(string, "r"); 

       if (!fs) {
           path = get_spectra_path(string);
           fs = fopen(path, "r");
        }
        if (fs == NULL) {
            printf("no database file for this material (%s)!\n", string);
            if (path) 
                fprintf(stderr, "Error: cannot open data from optical database: %s\n", path);
            return 1;
        } else {
            printf("file successfully opened\n");
            fscanf(fs, "%s", stype);
            fscanf(fs, "%f", &swfrom);
            fscanf(fs, "%f", &swto);
            description = fgets(description, 256, fs);
            description = fgets(description, 256, fs);
            if (!description)
                return 1;
            //getline(&description, &len, fs);
            //getline(&description, &len, fs);
            nl = strrchr(description, '\r');
            if (nl)
                *nl = '\0';
            nl = strrchr(description, '\n');
            if (nl)
                *nl = '\0';

            swfrom *= (float)1.0e-6; //spectra ranges are in microns
            swto *= (float)1.0e-6;

            printf("Spectrum type %s, range %g %g micrometers (%s) used\n", stype, swfrom, swto, description);
            if (swfrom > wfrom || swto < wto)
                fprintf(stderr, "Warning: your source wavelength (%g to %g, center at %g) does not match database spectrum range (%g, %g)\n", wfrom, wcenter, wto, swfrom, swto);

            printf("trying to match spectrum type '%s'\n", stype);
            if (strcmp(stype, "drude") == 0) {
                fprintf(fw, "%d ", SV_MAT_DRUDE);
                fscanf(fs, "%f", &epsilon); fprintf(fw, "%g ", epsilon);
                fscanf(fs, "%f", &omegap); fprintf(fw, "%g ", omegap);
                fscanf(fs, "%f", &nu); fprintf(fw, "%g\n", nu);
            } else if (strcmp(stype, "cp") == 0) {
                fprintf(fw, "%d ", SV_MAT_CP);
                fscanf(fs, "%f", &epsilon); fprintf(fw, "%g ", epsilon);
                fscanf(fs, "%f", &omegap); fprintf(fw, "%g ", omegap);
                fscanf(fs, "%f", &nu); fprintf(fw, "%g ", nu);
                fscanf(fs, "%f", &a); fprintf(fw, "%g ", a);
                fscanf(fs, "%f", &phi); fprintf(fw, "%g ", phi);
                fscanf(fs, "%f", &omega); fprintf(fw, "%g ", omega);
                fscanf(fs, "%f", &gamma); fprintf(fw, "%g ", gamma);
                fscanf(fs, "%f", &a); fprintf(fw, "%g ", a);
                fscanf(fs, "%f", &phi); fprintf(fw, "%g ", phi);
                fscanf(fs, "%f", &omega); fprintf(fw, "%g ", omega);
                fscanf(fs, "%f", &gamma); fprintf(fw, "%g\n", gamma);
            } else if (strcmp(stype, "ade") == 0) {
                fprintf(fw, "%d ", SV_MAT_ADE);
                fscanf(fs, "%f", &epsilon); fprintf(fw, "%g ", epsilon);
                fscanf(fs, "%g", &sigma); fprintf(fw, "%g ", sigma);
                fscanf(fs, "%f", &omegap); fprintf(fw, "%g ", omegap);
                fscanf(fs, "%f", &nu); fprintf(fw, "%g ", nu);
                fscanf(fs, "%f", &a); fprintf(fw, "%g ", a);
                fscanf(fs, "%f", &phi); fprintf(fw, "%g ", phi);
                fscanf(fs, "%f", &omega); fprintf(fw, "%g ", omega);
                fscanf(fs, "%f", &gamma); fprintf(fw, "%g ", gamma);
                fscanf(fs, "%f", &a); fprintf(fw, "%g ", a);
                fscanf(fs, "%f", &phi); fprintf(fw, "%g ", phi);
                fscanf(fs, "%f", &omega); fprintf(fw, "%g ", omega);
                fscanf(fs, "%f", &gamma); fprintf(fw, "%g\n", gamma);
            } else if (strcmp(stype, "plrc") == 0) {
                fprintf(fw, "%d ", SV_MAT_PLRC);
                fscanf(fs, "%f", &epsilon); fprintf(fw, "%g ", epsilon);
                fscanf(fs, "%g", &sigma); fprintf(fw, "%g ", sigma);
                fscanf(fs, "%f", &omegap); fprintf(fw, "%g ", omegap);
                fscanf(fs, "%f", &nu); fprintf(fw, "%g ", nu);
                fscanf(fs, "%f", &a); fprintf(fw, "%g ", a);
                fscanf(fs, "%f", &phi); fprintf(fw, "%g ", phi);
                fscanf(fs, "%f", &omega); fprintf(fw, "%g ", omega);
                fscanf(fs, "%f", &gamma); fprintf(fw, "%g ", gamma);
                fscanf(fs, "%f", &a); fprintf(fw, "%g ", a);
                fscanf(fs, "%f", &phi); fprintf(fw, "%g ", phi);
                fscanf(fs, "%f", &omega); fprintf(fw, "%g ", omega);
                fscanf(fs, "%f", &gamma); fprintf(fw, "%g\n", gamma);
            } else if (strcmp(stype, "cp3") == 0) {
                fprintf(fw, "%d ", SV_MAT_CP3);
                fscanf(fs, "%g", &epsilon); fprintf(fw, "%g ", epsilon);
                fscanf(fs, "%g", &sigma); fprintf(fw, "%g ", sigma);
                fscanf(fs, "%f", &a); fprintf(fw, "%g ", a);
                fscanf(fs, "%f", &phi); fprintf(fw, "%g ", phi);
                fscanf(fs, "%f", &omega); fprintf(fw, "%g ", omega);
                fscanf(fs, "%f", &gamma); fprintf(fw, "%g ", gamma);
                fscanf(fs, "%f", &a); fprintf(fw, "%g ", a);
                fscanf(fs, "%f", &phi); fprintf(fw, "%g ", phi);
                fscanf(fs, "%f", &omega); fprintf(fw, "%g ", omega);
                fscanf(fs, "%f", &gamma); fprintf(fw, "%g ", gamma);
                fscanf(fs, "%f", &a); fprintf(fw, "%g ", a);
                fscanf(fs, "%f", &phi); fprintf(fw, "%g ", phi);
                fscanf(fs, "%f", &omega); fprintf(fw, "%g ", omega);
                fscanf(fs, "%f", &gamma); fprintf(fw, "%g\n", gamma);
            } else if (strcmp(stype, "sellmeier") == 0) {
                fprintf(fw, "%d ", SV_MAT_LINTAB);
                fscanf(fs, "%g", &par1);
                fscanf(fs, "%g", &par2);
                fscanf(fs, "%g", &par3);
                fscanf(fs, "%g", &par4);
                fscanf(fs, "%g", &par5);
                fscanf(fs, "%g", &par6);
                fscanf(fs, "%g", &par7);

                wav = wcenter*1e6; //wavelength in um
                n = sqrt((double)(par1 + par2*wav*wav / (wav*wav - par3) + par4*wav*wav / (wav*wav - par5) + par6*wav*wav / (wav*wav - par7)));

                eps = n*n;

                fprintf(fw, "%g 1 0 0\n", eps);
            } else if (strcmp(stype, "cauchy") == 0) {
                fprintf(fw, "%d ", SV_MAT_LINTAB);
                fscanf(fs, "%g", &par1);
                fscanf(fs, "%g", &par2);
                fscanf(fs, "%g", &par3);

                wav = wcenter*1e6; //wavelength in um
                n = par1 + par2 / wav / wav + par3 / wav / wav / wav / wav;
                eps = n*n;

                fprintf(fw, "%g 1 0 0\n", eps);
            } else if (strcmp(stype, "kasarova") == 0) {
                fprintf(fw, "%d ", SV_MAT_LINTAB);
                fscanf(fs, "%g", &par1);
                fscanf(fs, "%g", &par2);
                fscanf(fs, "%g", &par3);
                fscanf(fs, "%g", &par4);
                fscanf(fs, "%g", &par5);
                fscanf(fs, "%g", &par6);

                wav = wcenter*1e6; //wavelength in um
                n = sqrt((double)(par1 + par2*wav*wav + par3 / wav / wav + par4 / wav / wav / wav / wav + par5 / wav / wav / wav / wav / wav / wav + par6 / wav / wav / wav / wav / wav / wav / wav / wav));
                eps = n*n;

                fprintf(fw, "%g 1 0 0\n", eps);
            } else if (strcmp(stype, "nk") == 0) {
                fprintf(fw, "%d ", SV_MAT_LINTAB);
                fscanf(fs, "%d", &ndata);

                wfs = (gdouble *)g_malloc(ndata * sizeof(gdouble));
                nfs = (gdouble *)g_malloc(ndata * sizeof(gdouble));
                kfs = (gdouble *)g_malloc(ndata * sizeof(gdouble));

                for (i = 0; i < (ndata); i++) {
                    fscanf(fs, "%g", &wf);
                    fscanf(fs, "%g", &nf);
                    fscanf(fs, "%g", &kf);

                    if (i == 0) {
                        if (fabs(wf*1e-6 - swfrom) < 1e-10) {
                            dir = 1;
                            mult = 1e-6;
                        } else if (fabs(wf*1e-6 - swto) < 1e-10) {
                            dir = 0;
                            mult = 1e-6;
                        } else if (fabs(wf*1e-9 - swfrom) < 1e-10) {
                            dir = 1;
                            mult = 1e-9;
                        } else if (fabs(wf*1e-9 - swto) < 1e-10) {
                            dir = 0;
                            mult = 1e-9;
                        } else {
                            fprintf(stderr, "Database entry error: wavelength range (%g %g) does not match the spectrum value %g\n", swfrom, swto, wf);
                            return 1;
                        }
                    }

                    wfs[i] = wf*mult;
                    nfs[i] = nf;
                    kfs[i] = kf;
                    //                    printf("%g %g %g %d\n", wfs[i], nfs[i], kfs[i], dir);
                }
                if (dir == 1 && wcenter < wfs[0]) {
                    wprev = wfs[0];
                    nprev = nfs[0];
                    kprev = kfs[0];
                    wnext = wfs[1];
                    nnext = nfs[1];
                    knext = kfs[1];
                    fprintf(stderr, "Warning: Spectral data extrapolated at lower side of spectrum\n");
                } else if (dir == 0 && wcenter < wfs[ndata - 1]) {
                    wprev = wfs[ndata - 1];
                    nprev = nfs[ndata - 1];
                    kprev = kfs[ndata - 1];
                    wnext = wfs[ndata - 2];
                    nnext = nfs[ndata - 2];
                    knext = kfs[ndata - 2];
                    fprintf(stderr, "Warning: Spectral data extrapolated at lower side of spectrum\n");
                } else {
                    for (i = 0; i < (ndata - 1); i++) {
                        if (dir == 1 && wfs[i] <= wcenter && wfs[i + 1] >= wcenter) { //data in table going upwards
                            wprev = wfs[i];
                            nprev = nfs[i];
                            kprev = kfs[i];
                            wnext = wfs[i + 1];
                            nnext = nfs[i + 1];
                            knext = kfs[i + 1];
                            break;
                        } else if (dir == 0 && wfs[i] >= wcenter && wfs[i + 1] <= wcenter) { //data going backwards
                            wprev = wfs[i + 1];
                            nprev = nfs[i + 1];
                            kprev = kfs[i + 1];
                            wnext = wfs[i];
                            nnext = nfs[i];
                            knext = kfs[i];
                            break;
                        }

                    }
                }
                n = nprev + (nnext - nprev)*(wcenter - wprev) / (wnext - wprev);
                k = kprev + (knext - kprev)*(wcenter - wprev) / (wnext - wprev);

                eps = n*n - k*k;
                sig = 4 * G_PI*n*k*LIGHT_SPEED / wcenter*EPSILON_0; //no 4 Pi before formerly

                printf("nk: %g %g %g %g for wavelength %g\n", n, k, eps, sig, wcenter);
                fprintf(fw, "%g 1 %g 0\n", eps, sig);

                g_free(wfs);
                g_free(nfs);
                g_free(kfs);
            } else {
                fprintf(stderr, "Error: wrong spectrum type (%s)\n", stype);
                return 1;
            }
        }

        fclose(fs);
    } else {
        fprintf(stderr, "Error: adj: Unknown material type?\n");
        return 1;
    }

    return 0;
}

/*override the material file if we work with database, using wavelength span wf-wt and central wavelength wc*/
gboolean
override_matfile_from_database(gchar *filename, gdouble wf, gdouble wc, gdouble wt, gint suffix, gint *nmaterials, gchar **materials, gint *overpos, gboolean localfiles)
{
    gint type, ret;
    GLfloat buf1, buf2, buf3;
    FILE *fr, *fw;
    gchar filebase[256];
    gchar matfile[256];
    gint pos = 0;

    fr = fopen(filename, "r");
    g_snprintf(matfile, sizeof(matfile), "tmp_matfile_%06d", suffix);
    fw = fopen(matfile, "w");
    if (fr == NULL) {
        fprintf(stderr, "Error: cannot open file %s for reading\n", filename);
        return 1;
    }
    if (fw == NULL) {
        fprintf(stderr, "Error: cannot open file %s for writing\n", matfile);
        return 1;
    }

    while ((ret = fscanf(fr, "%d", &type)) != EOF) { //TODO use proper locale handling functions here
        if (!ret) {
            printf("Error: something strange in material input file\n");
            fclose(fr);
            fclose(fw);
            return 1;
        }

        fprintf(fw, "%d ", type);
        switch (type) {
            case SV_ENTITY_SPHERE:
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (scan_radius(fr, &buf1))
                    goto exit;
                fprintf(fw, "%g ", buf1);
                if (adjust_material(fr, fw, wf, wc, wt, nmaterials, materials, pos, overpos, localfiles))
                    goto exit;
                break;

            case SV_ENTITY_CYLINDER:
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (scan_radius(fr, &buf1))
                    goto exit;
                fprintf(fw, "%g ", buf1);
                if (adjust_material(fr, fw, wf, wc, wt, nmaterials, materials, pos, overpos, localfiles))
                    goto exit;
                break;

            case SV_ENTITY_CONE:
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (scan_radius(fr, &buf1))
                    goto exit;
                fprintf(fw, "%g ", buf1);
                if (adjust_material(fr, fw, wf, wc, wt, nmaterials, materials, pos, overpos, localfiles))
                    goto exit;
                break;

            case SV_ENTITY_RCONE:
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (scan_radius(fr, &buf1))
                    goto exit;
                fprintf(fw, "%g ", buf1);
                if (scan_radius(fr, &buf1))
                    goto exit;
                fprintf(fw, "%g ", buf1);
                if (adjust_material(fr, fw, wf, wc, wt, nmaterials, materials, pos, overpos, localfiles))
                    goto exit;
                break;

            case SV_ENTITY_VOXEL:
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (adjust_material(fr, fw, wf, wc, wt, nmaterials, materials, pos, overpos, localfiles))
                    goto exit;
                break;

            case SV_ENTITY_TETRAHEDRON:
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (adjust_material(fr, fw, wf, wc, wt, nmaterials, materials, pos, overpos, localfiles))
                    goto exit;
                break;

            case SV_ENTITY_TETGEN:
                if (fscanf(fr, "%s", filebase) != 1)
                    goto exit;
                fprintf(fw, "%s ", filebase);
                if (scan_radius(fr, &buf1))
                    goto exit;
                fprintf(fw, "%d ", (gint)buf1);
                if (scan_radius(fr, &buf1))
                    goto exit;
                fprintf(fw, "%d ", (gint)buf1);
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (adjust_material(fr, fw, wf, wc, wt, nmaterials, materials, pos, overpos, localfiles))
                    goto exit;
                break;

            case SV_ENTITY_GWYDDION:
                if (fscanf(fr, "%s", filebase) != 1)
                    goto exit;
                fprintf(fw, "%s ", filebase);
                if (scan_radius(fr, &buf1))
                    goto exit;
                fprintf(fw, "%d ", (gint)buf1);
                if (scan_radius(fr, &buf1))
                    goto exit;
                fprintf(fw, "%d ", (gint)buf1);
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%d %d %d ", (gint)buf1, (gint)buf2, (gint)buf3);
                if (scan_point(fr, &buf1, &buf2, &buf3))
                    goto exit;
                fprintf(fw, "%g %g %g ", buf1, buf2, buf3);
                if (scan_radius(fr, &buf1))
                    goto exit;
                fprintf(fw, "%d ", (gint)buf1);
                if (adjust_material(fr, fw, wf, wc, wt, nmaterials, materials, pos, overpos, localfiles))
                    goto exit;
                break;
            default:
                break;
        }
        pos++;
    }

    fclose(fr);
    fclose(fw);

    return 0;

exit:
    fclose(fr);
    fclose(fw);
    return 1;
}

gboolean
estimate_lambda(gdouble *data, gint ndata, gdouble dt, gdouble *lambda_min, gdouble *lambda_max, gint savespectrum)
{
    GwyDataLine *dline, *spectrum;
    GwyDataLine *iin, *rout, *iout;
    gdouble *linedata, *rdata, *idata;
    gdouble maximum = -G_MAXDOUBLE;
    gdouble wl;
    gint i, res;
    FILE *fw;

    if (ndata<10) return FALSE;

    dline = gwy_data_line_new(ndata, ndata*dt, FALSE);
    linedata = gwy_data_line_get_data(dline);

    for (i = 0; i < ndata; i++) {
        linedata[i] = data[i];
    }
    gwy_data_line_data_changed(dline);

    spectrum = gwy_data_line_new_alike(dline, TRUE);

    /*************** code pasted from Gwyddion because of a bug there ***************/
    res = dline->res;
    iin = gwy_data_line_new_alike(dline, TRUE);
    rout = gwy_data_line_new_alike(dline, FALSE);
    iout = gwy_data_line_new_alike(dline, FALSE);
    gwy_data_line_resample(spectrum, res / 2, GWY_INTERPOLATION_NONE);

    gwy_data_line_fft(dline, iin, rout, iout,
                      GWY_WINDOWING_NONE,
                      GWY_TRANSFORM_DIRECTION_FORWARD,
                      GWY_INTERPOLATION_BILINEAR,
                      TRUE, 2);

    linedata = spectrum->data;
    rdata = rout->data;
    idata = iout->data;

    /* Calculate modulus */
    for (i = 0; i < res / 2; i++) {
        linedata[i] = (rdata[i] * rdata[i] + idata[i] * idata[i]);
    }
    gwy_data_line_data_changed(spectrum);

    g_object_unref(rout);
    g_object_unref(iin);
    g_object_unref(iout);

    /************** end of code from Gwyddion *****************/

    ndata = gwy_data_line_get_res(spectrum);
    maximum = gwy_data_line_get_max(spectrum);
    linedata = gwy_data_line_get_data(spectrum);

    if (maximum == 0) {
        g_object_unref(dline);
        g_object_unref(spectrum);
        return FALSE;
    }

    *lambda_min = 1e10;
    *lambda_max = -1e10;
    for (i = 1; i < ndata; i++) {
       if (linedata[i] >= (maximum / 2)) {
           wl = LIGHT_SPEED / ((gdouble)i / ((gdouble)ndata*dt*2.0));
           *lambda_min = MIN(*lambda_min, wl);
           *lambda_max = MAX(*lambda_max, wl);
       }
    }
 
    if (savespectrum)
    {
       fw = fopen("xspectrum.txt", "w");
       for (i = 1; i < ndata; i++) {
           fprintf(fw, "%g %g\n", LIGHT_SPEED / ((gdouble)i / ((gdouble)ndata*dt*2.0)), sqrt(linedata[i]));
       }
       fclose(fw);
    }

    g_object_unref(dline);
    g_object_unref(spectrum);
    return TRUE;
}

void
get_gl_coeffs(gdouble *xs, gdouble *ws, gint n)
{
    gdouble z1, z, pp, p1, p2, p3;
    gint i, j;

    for (i = 0; i < ((n + 1) / 2); i++) {
        z = cos(G_PI*(i + 0.75) / (n + 0.5));
        do {
            p1 = 1.0;
            p2 = 0.0;
            for (j = 0; j < n; j++) {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0*j + 1.0)*z*p2 - (gdouble)j*p3) / ((gdouble)j + 1);
            }
            pp = n*(z*p1 - p2) / (z*z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (fabs(z - z1) > 1e-15);
        xs[i] = -z;
        xs[n - 1 - i] = z;
        ws[i] = 2.0 / ((1.0 - z*z)*pp*pp);
        ws[n - 1 - i] = ws[i];
    }
}

gboolean
load_sources(SvPool *mp, SvSet *set)
{
    FILE *fr;
    gint i, k, ndata, pos, i1, im;
    gdouble ex, ey, ez, hx, hy, hz, lmin = G_MAXDOUBLE, lmax = -G_MAXDOUBLE;
    gdouble *sum, *xsm, *wsm, *xsn, *wsn;
    gdouble dist;
    gdouble ax, ay, az, rx, ry, rz;

    if (set->ss.npnts)
        mp->src->sp = (SvSourcePoint *)g_malloc(set->ss.npnts * sizeof(SvSourcePoint));

    for (i = 0; i < set->ss.npnts; i++) {
        fr = fopen(set->ss.pnts[i].source_filename, "r");
        if (fr == NULL) {
            fprintf(stderr, "Error: cannot open file %s for reading\n", set->ss.pnts[i].source_filename);
            continue;
        }

        if (fscanf(fr, "%d", &ndata) != 1) {
            fprintf(stderr, "Error: cannot read number of values from file %s\n", set->ss.pnts[i].source_filename);
            fclose(fr);
            return 1;
        }
        if (ndata <= 0) {
            fprintf(stderr, "Error: zero source steps in %s\n", set->ss.pnts[i].source_filename);
            fclose(fr);
            return 1;
        }
        if (ndata < set->sc.nsteps)
            fprintf(stderr, "Warning: less data (%d) in \"%s\" than scheduled computation steps (%d)\n", ndata, set->ss.pnts[i].source_filename, set->sc.nsteps);
        mp->src->sp[i].i = set->ss.pnts[i].point_origin_position_i;
        mp->src->sp[i].j = set->ss.pnts[i].point_origin_position_j;
        mp->src->sp[i].k = set->ss.pnts[i].point_origin_position_k;
        mp->src->sp[i].sdata.ndata = ndata;
        mp->src->sp[i].sdata.layered_zpos = (gint *)g_malloc(ndata * sizeof(gint));
        mp->src->sp[i].sdata.ex = (gdouble *)g_malloc(ndata * sizeof(gdouble));
        mp->src->sp[i].sdata.ey = (gdouble *)g_malloc(ndata * sizeof(gdouble));
        mp->src->sp[i].sdata.ez = (gdouble *)g_malloc(ndata * sizeof(gdouble));
        mp->src->sp[i].sdata.hx = (gdouble *)g_malloc(ndata * sizeof(gdouble));
        mp->src->sp[i].sdata.hy = (gdouble *)g_malloc(ndata * sizeof(gdouble));
        mp->src->sp[i].sdata.hz = (gdouble *)g_malloc(ndata * sizeof(gdouble));
        sum = (gdouble *)g_malloc(ndata * sizeof(gdouble));

        for (k = 0; k < ndata; k++) {
            if (get_int(fr, &pos, "PointSource")) {
                fclose(fr);
                return 1;
            }
            if (get_double(fr, &ex, "PointSource")) {
                fclose(fr);
                return 1;
            }
            if (get_double(fr, &ey, "PointSource")) {
                fclose(fr);
                return 1;
            }
            if (get_double(fr, &ez, "PointSource")) {
                fclose(fr);
                return 1;
            }
            if (get_double(fr, &hx, "PointSource")) {
                fclose(fr);
                return 1;
            }
            if (get_double(fr, &hy, "PointSource")) {
                fclose(fr);
                return 1;
            }
            if (get_double(fr, &hz, "PointSource")) {
                fclose(fr);
                return 1;
            }

            mp->src->sp[i].sdata.layered_zpos[k] = pos;
            mp->src->sp[i].sdata.ex[k] = ex;
            mp->src->sp[i].sdata.ey[k] = ey;
            mp->src->sp[i].sdata.ez[k] = ez;
            mp->src->sp[i].sdata.hx[k] = hx;
            mp->src->sp[i].sdata.hy[k] = hy;
            mp->src->sp[i].sdata.hz[k] = hz;
            sum[k] = (ex + ey + ez);
        }

        fclose(fr);
        if (estimate_lambda(sum, ndata, set->plan.dt, &lmin, &lmax, set->so.savespectrum)) {
            if (set->sc.verbose > 1)
                printf("Estimated wavelength range: (%g, %g) m  after file %s\n", lmin, lmax, set->ss.pnts[i].source_filename);
        } else if (set->sc.verbose > 1)
            printf("Automated detection of wavelength failed on file %s\n", set->ss.pnts[i].source_filename);
    }

    if (set->ss.tsf.valid) {
        mp->src->tsf = (SvSourceTSF *)g_malloc(sizeof(SvSourceTSF));
        mp->src->tsf->box_i0 = set->ss.tsf.box_i0;
        mp->src->tsf->box_j0 = set->ss.tsf.box_j0;
        mp->src->tsf->box_k0 = set->ss.tsf.box_k0;
        mp->src->tsf->box_in = set->ss.tsf.box_in;
        mp->src->tsf->box_jn = set->ss.tsf.box_jn;
        mp->src->tsf->box_kn = set->ss.tsf.box_kn;
        mp->src->tsf->box_boundary_skipi0 = set->ss.tsf.box_boundary_skipi0;
        mp->src->tsf->box_boundary_skipj0 = set->ss.tsf.box_boundary_skipj0;
        mp->src->tsf->box_boundary_skipk0 = set->ss.tsf.box_boundary_skipk0;
        mp->src->tsf->box_boundary_skipin = set->ss.tsf.box_boundary_skipin;
        mp->src->tsf->box_boundary_skipjn = set->ss.tsf.box_boundary_skipjn;
        mp->src->tsf->box_boundary_skipkn = set->ss.tsf.box_boundary_skipkn;
        mp->src->tsf->ia_theta = set->ss.tsf.ia_theta;
        mp->src->tsf->ia_phi = set->ss.tsf.ia_phi;
        mp->src->tsf->ia_psi = set->ss.tsf.ia_psi;
        mp->src->tsf->layered_epsilon = set->ss.tsf.layered_epsilon;
        mp->src->tsf->layered_mu = set->ss.tsf.layered_mu;
        mp->src->tsf->layered_sigma = set->ss.tsf.layered_sigma;
        mp->src->tsf->layered_sigast = set->ss.tsf.layered_sigast;

        fr = fopen(set->ss.tsf.source_filename, "r");
        if (fr == NULL) {
            fprintf(stderr, "Cannot load source file %s\n", set->ss.tsf.source_filename);
            return FALSE;
        }
        if (fscanf(fr, "%d", &ndata) != 1) {
            fprintf(stderr, "Error: cannot read number of values from file %s\n", set->ss.tsf.source_filename);
            fclose(fr);
            return 1;
        }
        if (ndata <= 0) {
            fprintf(stderr, "Error: zero source steps in %s\n", set->ss.tsf.source_filename);
            fclose(fr);
            return 1;
        }
        if (ndata < set->sc.nsteps)
            fprintf(stderr, "Warning: less data (%d) in \"%s\" than scheduled computation steps (%d) in TSF source\n", ndata, set->ss.pnts[i].source_filename, set->sc.nsteps);
        mp->src->tsf->ndata = ndata;
        mp->src->tsf->layered_zpos = (gint *)g_malloc(ndata * sizeof(gint));
        mp->src->tsf->e = (gdouble *)g_malloc(ndata * sizeof(gdouble));

        for (k = 0; k < ndata; k++) {
            if (get_int(fr, &pos, "TSF Source")) {
                fclose(fr);
                return 1;
            }
            if (get_double(fr, &ex, "TSF Source")) {
                fclose(fr);
                return 1;
            }

            mp->src->tsf->layered_zpos[k] = pos;
            mp->src->tsf->e[k] = ex;
        }
        if (estimate_lambda(mp->src->tsf->e, ndata, set->plan.dt, &lmin, &lmax, set->so.savespectrum)) {
            if (set->sc.verbose > 1)
                printf("Estimated wavelength range: (%g, %g) m  after file %s\n", lmin, lmax, set->ss.tsf.source_filename);
        } else if (set->sc.verbose > 1)
            printf("Automated detection of wavelength failed on file %s\n", set->ss.tsf.source_filename);

        mp->src->tsf->corr = 1.0 / sqrt(pow(sin(mp->src->tsf->ia_theta), 4)
                                        *(pow(cos(mp->src->tsf->ia_phi), 4) + pow(sin(mp->src->tsf->ia_phi), 4))
                                        + pow(cos(mp->src->tsf->ia_theta), 4));
        if (set->sc.verbose > 1)
            printf("Angular correction factor for TSF: %g\n", mp->src->tsf->corr);


        mp->src->tsf->jp = sv_1dpool_new((gint)(sqrt((gdouble)(set->sp.xres*set->sp.xres + set->sp.yres*set->sp.yres + set->sp.zres*set->sp.zres)) + 250),
                                         set->sp.dx / mp->src->tsf->corr, set->plan.dt);

        sv_1dpool_set_source_pos(mp->src->tsf->jp, 0);
        sv_1dpool_set_material(mp->src->tsf->jp, mp->src->tsf->layered_epsilon, mp->src->tsf->layered_mu);

        fclose(fr);
    }

    if (set->ss.ltsf.valid) {
        mp->src->ltsf = (SvSourceLTSF *)g_malloc(sizeof(SvSourceLTSF));
        mp->src->ltsf->box_i0 = set->ss.ltsf.box_i0;
        mp->src->ltsf->box_j0 = set->ss.ltsf.box_j0;
        mp->src->ltsf->box_k0 = set->ss.ltsf.box_k0;
        mp->src->ltsf->box_in = set->ss.ltsf.box_in;
        mp->src->ltsf->box_jn = set->ss.ltsf.box_jn;
        mp->src->ltsf->box_kn = set->ss.ltsf.box_kn;
        mp->src->ltsf->box_boundary_skipi0 = set->ss.ltsf.box_boundary_skipi0;
        mp->src->ltsf->box_boundary_skipj0 = set->ss.ltsf.box_boundary_skipj0;
        mp->src->ltsf->box_boundary_skipk0 = set->ss.ltsf.box_boundary_skipk0;
        mp->src->ltsf->box_boundary_skipin = set->ss.ltsf.box_boundary_skipin;
        mp->src->ltsf->box_boundary_skipjn = set->ss.ltsf.box_boundary_skipjn;
        mp->src->ltsf->box_boundary_skipkn = set->ss.ltsf.box_boundary_skipkn;
        mp->src->ltsf->ia_theta = set->ss.ltsf.ia_theta;
        mp->src->ltsf->ia_phi = set->ss.ltsf.ia_phi;
        mp->src->ltsf->ia_psi = set->ss.ltsf.ia_psi;
        mp->src->ltsf->layered_count = set->ss.ltsf.layered_count;
        for (k = 0; k < mp->src->ltsf->layered_count; k++) {
            mp->src->ltsf->lpos[k] = set->ss.ltsf.layered_zpos[k];
            mp->src->ltsf->layered_epsilon[k] = set->ss.ltsf.layered_epsilon[k];
            mp->src->ltsf->layered_mu[k] = set->ss.ltsf.layered_mu[k];
            mp->src->ltsf->layered_sigma[k] = set->ss.ltsf.layered_sigma[k];
            mp->src->ltsf->layered_sigast[k] = set->ss.ltsf.layered_sigast[k];
        }

        fr = fopen(set->ss.ltsf.source_filename, "r");
        if (fr == NULL) {
            fprintf(stderr, "Cannot load source file %s\n", set->ss.ltsf.source_filename);
            return FALSE;
        }
        if (fscanf(fr, "%d", &ndata) != 1) {
            fprintf(stderr, "Error: cannot read number of values from file %s\n", set->ss.ltsf.source_filename);
            fclose(fr);
            return 1;
        }

        if (ndata <= 0) {
            fprintf(stderr, "Error: zero source steps in %s\n", set->ss.ltsf.source_filename);
            fclose(fr);
            return 1;
        }
        if (ndata < (set->sc.nsteps * 3))
            fprintf(stderr, "Warning: less data (%d) in \"%s\" than scheduled computation steps (%d) in LTSF source (need to be three times more)\n", ndata, set->ss.pnts[i].source_filename, set->sc.nsteps);
        mp->src->ltsf->ndata = ndata;
        mp->src->ltsf->layered_zpos = (gint *)g_malloc(ndata * sizeof(gint));
        mp->src->ltsf->e = (gdouble *)g_malloc(ndata * sizeof(gdouble));

        for (k = 0; k < ndata; k++) {
            if (get_int(fr, &pos, "LTSF Source")) {
                fclose(fr);
                return 1;
            }
            if (get_double(fr, &ex, "LTSF Source")) {
                fclose(fr);
                return 1;
            }

            mp->src->ltsf->layered_zpos[k] = pos;
            mp->src->ltsf->e[k] = ex;
        }
        if (estimate_lambda(mp->src->ltsf->e, ndata, set->plan.dt, &lmin, &lmax, set->so.savespectrum)) {
            lmin /= 3.0;
            lmax /= 3.0;

            if (set->sc.verbose > 1)
                printf("Estimated wavelength range: (%g, %g) m  after file %s\n", lmin, lmax, set->ss.ltsf.source_filename);
        } else if (set->sc.verbose > 1)
            printf("Automated detection of wavelength failed on file %s\n", set->ss.ltsf.source_filename);

        mp->src->ltsf->corr = 1.0 / sqrt(pow(sin(mp->src->ltsf->ia_theta), 4)
                                         *(pow(cos(mp->src->ltsf->ia_phi), 4) + pow(sin(mp->src->ltsf->ia_phi), 4))
                                         + pow(cos(mp->src->ltsf->ia_theta), 4));
        //printf("LTSF correction factor %g\n", mp->src->ltsf->corr);
        if (set->sc.verbose > 1)
            printf("Angular correction factor for LTSF: %g\n", mp->src->ltsf->corr);
        mp->src->ltsf->jp = sv_zpool_new(set->sp.zres, set->sp.dx / mp->src->ltsf->corr, set->plan.dt, set->sc.nsteps);
        mp->src->ltsf->timeshift = 12;

        sv_zpool_set_source_pos(mp->src->ltsf->jp, 0);
        sv_zpool_set_correction(mp->src->ltsf->jp, mp->src->ltsf->corr);
        sv_zpool_set_angle(mp->src->ltsf->jp, mp->src->ltsf->ia_theta, mp->src->ltsf->ia_psi);
        sv_zpool_set_material(mp->src->ltsf->jp, mp->src->ltsf->layered_count, mp->src->ltsf->lpos, mp->src->ltsf->layered_epsilon, mp->src->ltsf->layered_mu, mp->src->ltsf->layered_sigma, mp->src->ltsf->layered_sigast);

        fclose(fr);
    }

    if (set->ss.ltsff.valid) {
        mp->src->ltsff = (SvSourceLTSFF *)g_malloc(sizeof(SvSourceLTSFF));
        mp->src->ltsff->box_i0 = set->ss.ltsff.box_i0;
        mp->src->ltsff->box_j0 = set->ss.ltsff.box_j0;
        mp->src->ltsff->box_k0 = set->ss.ltsff.box_k0;
        mp->src->ltsff->box_in = set->ss.ltsff.box_in;
        mp->src->ltsff->box_jn = set->ss.ltsff.box_jn;
        mp->src->ltsff->box_kn = set->ss.ltsff.box_kn;
        mp->src->ltsff->box_boundary_skipi0 = set->ss.ltsff.box_boundary_skipi0;
        mp->src->ltsff->box_boundary_skipj0 = set->ss.ltsff.box_boundary_skipj0;
        mp->src->ltsff->box_boundary_skipk0 = set->ss.ltsff.box_boundary_skipk0;
        mp->src->ltsff->box_boundary_skipin = set->ss.ltsff.box_boundary_skipin;
        mp->src->ltsff->box_boundary_skipjn = set->ss.ltsff.box_boundary_skipjn;
        mp->src->ltsff->box_boundary_skipkn = set->ss.ltsff.box_boundary_skipkn;
        mp->src->ltsff->focused_thetamax = set->ss.ltsff.focused_thetamax;
        mp->src->ltsff->focused_fdist = set->ss.ltsff.focused_fdist;
        mp->src->ltsff->focused_nip = set->ss.ltsff.focused_nip;
        mp->src->ltsff->focused_mip = set->ss.ltsff.focused_mip;
        mp->src->ltsff->layered_count = set->ss.ltsff.layered_count;
        for (k = 0; k < mp->src->ltsff->layered_count; k++) {
            mp->src->ltsff->lpos[k] = set->ss.ltsff.layered_zpos[k];
            mp->src->ltsff->layered_epsilon[k] = set->ss.ltsff.layered_epsilon[k];
            mp->src->ltsff->layered_mu[k] = set->ss.ltsff.layered_mu[k];
            mp->src->ltsff->layered_sigma[k] = set->ss.ltsff.layered_sigma[k];
            mp->src->ltsff->layered_sigast[k] = set->ss.ltsff.layered_sigast[k];
        }

        fr = fopen(set->ss.ltsff.source_filename, "r");
        if (fr == NULL) {
            fprintf(stderr, "Cannot load source file %s\n", set->ss.ltsff.source_filename);
            return FALSE;
        }
        if (fscanf(fr, "%d", &ndata) != 1) {
            fprintf(stderr, "Error: cannot read number of values from file %s\n", set->ss.ltsff.source_filename);
            fclose(fr);
            return 1;
        }

        if (ndata <= 0) {
            fprintf(stderr, "Error: zero source steps in %s\n", set->ss.ltsff.source_filename);
            fclose(fr);
            return 1;
        }
        if (ndata < set->sc.nsteps)
            fprintf(stderr, "Warning: less data (%d) in \"%s\" than scheduled computation steps (%d) in TSFF source\n", ndata, set->ss.pnts[i].source_filename, set->sc.nsteps);
        mp->src->ltsff->ndata = ndata;
        mp->src->ltsff->layered_zpos = (gint *)g_malloc(ndata * sizeof(gint));
        mp->src->ltsff->e = (gdouble *)g_malloc(ndata * sizeof(gdouble));

        for (k = 0; k < ndata; k++) {
            if (get_int(fr, &pos, "LTSFF Source")) {
                fclose(fr);
                return 1;
            }
            if (get_double(fr, &ex, "LTSFF Source")) {
                fclose(fr);
                return 1;
            }

            mp->src->ltsff->layered_zpos[k] = pos;
            mp->src->ltsff->e[k] = ex;//(ex-pex)/set->plan.dt;
        }
        if (estimate_lambda(mp->src->ltsff->e, ndata, set->plan.dt, &lmin, &lmax, set->so.savespectrum)) {
            lmin /= 3.0;
            lmax /= 3.0;

            if (set->sc.verbose > 1)
                printf("Estimated wavelength range: (%g, %g) m  after file %s\n", lmin, lmax, set->ss.ltsff.source_filename);
        } else if (set->sc.verbose > 1) printf("Automated detection of wavelength failed on file %s\n", set->ss.ltsff.source_filename);

        mp->src->ltsff->jp = (SvZPool ***)g_malloc(mp->src->ltsff->focused_nip * sizeof(SvZPool **));
        mp->src->ltsff->ia_theta = (gdouble **)g_malloc(mp->src->ltsff->focused_nip * sizeof(gdouble *));
        mp->src->ltsff->ia_phi = (gdouble **)g_malloc(mp->src->ltsff->focused_nip * sizeof(gdouble *));
        mp->src->ltsff->ia_psi = (gdouble **)g_malloc(mp->src->ltsff->focused_nip * sizeof(gdouble *));
        mp->src->ltsff->corr = (gdouble **)g_malloc(mp->src->ltsff->focused_nip * sizeof(gdouble *));
        mp->src->ltsff->an = (gdouble **)g_malloc(mp->src->ltsff->focused_nip * sizeof(gdouble *));
        mp->src->ltsff->bm = (gdouble **)g_malloc(mp->src->ltsff->focused_nip * sizeof(gdouble *));
        mp->src->ltsff->timeshift = (gdouble **)g_malloc(mp->src->ltsff->focused_nip * sizeof(gdouble *));

        xsn = (gdouble *)g_malloc(mp->src->ltsff->focused_nip * sizeof(gdouble));
        wsn = (gdouble *)g_malloc(mp->src->ltsff->focused_nip * sizeof(gdouble));
        xsm = (gdouble *)g_malloc(mp->src->ltsff->focused_mip * sizeof(gdouble));
        wsm = (gdouble *)g_malloc(mp->src->ltsff->focused_mip * sizeof(gdouble));

        get_gl_coeffs(xsm, wsm, mp->src->ltsff->focused_mip);
        get_gl_coeffs(xsn, wsn, mp->src->ltsff->focused_nip);

        for (i1 = 0; i1 < mp->src->ltsff->focused_nip; i1++) {
            mp->src->ltsff->jp[i1] = (SvZPool **)g_malloc(mp->src->ltsff->focused_mip * sizeof(SvZPool *));
            mp->src->ltsff->ia_theta[i1] = (gdouble *)g_malloc(mp->src->ltsff->focused_mip * sizeof(gdouble));
            mp->src->ltsff->ia_phi[i1] = (gdouble *)g_malloc(mp->src->ltsff->focused_mip * sizeof(gdouble));
            mp->src->ltsff->ia_psi[i1] = (gdouble *)g_malloc(mp->src->ltsff->focused_mip * sizeof(gdouble));
            mp->src->ltsff->corr[i1] = (gdouble *)g_malloc(mp->src->ltsff->focused_mip * sizeof(gdouble));
            mp->src->ltsff->an[i1] = (gdouble *)g_malloc(mp->src->ltsff->focused_mip * sizeof(gdouble));
            mp->src->ltsff->bm[i1] = (gdouble *)g_malloc(mp->src->ltsff->focused_mip * sizeof(gdouble));
            mp->src->ltsff->timeshift[i1] = (gdouble *)g_malloc(mp->src->ltsff->focused_mip * sizeof(gdouble));

            for (im = 0; im < mp->src->ltsff->focused_mip; im++) {
                /*// extended midpoint
                  mp->src->ltsff->theta[in][im] = (in + 0.5)*mp->src->ltsff->thetamax/mp->src->ltsff->nint;
                  mp->src->ltsff->phi[in][im] = mp->src->ltsff->psi[in][im] = (im + 0.5)*2*G_PI/mp->src->ltsff->mint;
                  mp->src->ltsff->an[in][im] = mp->src->ltsff->thetamax/mp->src->ltsff->nint;
                  mp->src->ltsff->bm[in][im] = 2.0*G_PI/mp->src->ltsff->mint;
                */

                /* //gauss legendre
                   mp->src->ltsff->theta[in][im] = (xsn[in]+1)*mp->src->ltsff->thetamax/2.0;
                   mp->src->ltsff->phi[in][im] =  mp->src->ltsff->psi[in][im] = (xsm[im]+1)*G_PI;
                   mp->src->ltsff->an[in][im] = wsn[in]*mp->src->ltsff->thetamax/2;
                   mp->src->ltsff->bm[in][im] = G_PI*wsm[im];*/

                   //mixed - GL in theta but midpoint in phi (which is periodic)
                mp->src->ltsff->ia_theta[i1][im] = (xsn[i1] + 1)*mp->src->ltsff->focused_thetamax / 2.0;
                mp->src->ltsff->ia_phi[i1][im] = (im + 0.5) * 2 * G_PI / mp->src->ltsff->focused_mip;
                mp->src->ltsff->ia_psi[i1][im] = mp->src->ltsff->ia_phi[i1][im] + set->ss.ltsff.focused_pol;//no add

                mp->src->ltsff->an[i1][im] = wsn[i1] * mp->src->ltsff->focused_thetamax / 2;
                mp->src->ltsff->bm[i1][im] = 2.0*G_PI / mp->src->ltsff->focused_mip;

                mp->src->ltsff->corr[i1][im] = 1.0 / sqrt(pow(sin(mp->src->ltsff->ia_theta[i1][im]), 4)
                                                          *(pow(cos(mp->src->ltsff->ia_phi[i1][im]), 4) + pow(sin(mp->src->ltsff->ia_phi[i1][im]), 4))
                                                          + pow(cos(mp->src->ltsff->ia_theta[i1][im]), 4));

                mp->src->ltsff->jp[i1][im] = sv_zpool_new(set->sp.zres, set->sp.dx / mp->src->ltsff->corr[i1][im], set->plan.dt, set->sc.nsteps);

                ax = sin(mp->src->ltsff->ia_theta[i1][im])*cos(mp->src->ltsff->ia_phi[i1][im]);
                ay = sin(mp->src->ltsff->ia_theta[i1][im])*sin(mp->src->ltsff->ia_phi[i1][im]);
                az = cos(mp->src->ltsff->ia_theta[i1][im]);

                if (mp->src->ltsff->ia_phi[i1][im] >= 0 && mp->src->ltsff->ia_phi[i1][im] <= (G_PI / 2.0)) {
                    rx = set->sp.xres / 2 - mp->src->ltsff->box_i0;
                    ry = set->sp.yres / 2 - mp->src->ltsff->box_j0;
                    rz = set->sp.zres / 2 - mp->src->ltsff->box_k0;
                } else if (mp->src->ltsff->ia_phi[i1][im] > (G_PI / 2.0) && mp->src->ltsff->ia_phi[i1][im] <= G_PI) {
                    rx = -(set->sp.xres / 2 - mp->src->ltsff->box_i0);
                    ry = -(set->sp.yres / 2 - mp->src->ltsff->box_jn);
                    rz = set->sp.zres / 2 - mp->src->ltsff->box_k0;
                } else if (mp->src->ltsff->ia_phi[i1][im] > G_PI && mp->src->ltsff->ia_phi[i1][im] <= (3.0*G_PI / 2.0)) {
                    rx = set->sp.xres / 2 - mp->src->ltsff->box_in;
                    ry = set->sp.yres / 2 - mp->src->ltsff->box_jn;
                    rz = set->sp.zres / 2 - mp->src->ltsff->box_k0;
                } else {
                    rx = set->sp.xres / 2 - mp->src->ltsff->box_i0;
                    ry = set->sp.yres / 2 - mp->src->ltsff->box_jn;
                    rz = set->sp.zres / 2 - mp->src->ltsff->box_k0;
                }
                dist = sqrt((set->sp.xres / 2 - mp->src->ltsff->box_i0)*(set->sp.xres / 2 - mp->src->ltsff->box_i0)
                            + (set->sp.yres / 2 - mp->src->ltsff->box_j0)*(set->sp.yres / 2 - mp->src->ltsff->box_j0)
                            + (set->sp.zres / 2 - mp->src->ltsff->box_k0)*(set->sp.zres / 2 - mp->src->ltsff->box_k0));

                //printf("theta %g phi %g, a (%g %g %g)  r (%g %g %g) angle %g  corr %g %g\n", 
                //       mp->src->ltsff->theta[in][im], mp->src->ltsff->phi[in][im], ax, ay, az, rx, ry, rz,
                //       acos((ax*rx + ay*ry + az*rz)/(sqrt(rx*rx + ry*ry + rz*rz)*sqrt(ax*ax+ay*ay+az*az))), 
                //       (dist - (ax*rx + ay*ry + az*rz))/(set->plan.dt*LIGHT_SPEED/set->sp.dx), mp->src->ltsff->timeshift[in][im]);

                mp->src->ltsff->timeshift[i1][im] = 3 * (dist - (ax*rx + ay*ry + az*rz)) / (set->plan.dt*LIGHT_SPEED / set->sp.dx);

                sv_zpool_set_source_pos(mp->src->ltsff->jp[i1][im], 0);
                sv_zpool_set_correction(mp->src->ltsff->jp[i1][im], mp->src->ltsff->corr[i1][im]);
                sv_zpool_set_angle(mp->src->ltsff->jp[i1][im], mp->src->ltsff->ia_theta[i1][im], mp->src->ltsff->ia_psi[i1][im]);
                sv_zpool_set_material(mp->src->ltsff->jp[i1][im], mp->src->ltsff->layered_count, mp->src->ltsff->lpos, mp->src->ltsff->layered_epsilon, mp->src->ltsff->layered_mu, mp->src->ltsff->layered_sigma, mp->src->ltsff->layered_sigast);
            }
        }

        fclose(fr);

    }

    if (set->ss.tsff.valid) {
        mp->src->tsff = (SvSourceTSFF *)g_malloc(sizeof(SvSourceTSFF));
        mp->src->tsff->box_i0 = set->ss.tsff.box_i0;
        mp->src->tsff->box_j0 = set->ss.tsff.box_j0;
        mp->src->tsff->box_k0 = set->ss.tsff.box_k0;
        mp->src->tsff->box_in = set->ss.tsff.box_in;
        mp->src->tsff->box_jn = set->ss.tsff.box_jn;
        mp->src->tsff->box_kn = set->ss.tsff.box_kn;
        mp->src->tsff->box_boundary_skipi0 = set->ss.tsff.box_boundary_skipi0;
        mp->src->tsff->box_boundary_skipj0 = set->ss.tsff.box_boundary_skipj0;
        mp->src->tsff->box_boundary_skipk0 = set->ss.tsff.box_boundary_skipk0;
        mp->src->tsff->box_boundary_skipin = set->ss.tsff.box_boundary_skipin;
        mp->src->tsff->box_boundary_skipjn = set->ss.tsff.box_boundary_skipjn;
        mp->src->tsff->box_boundary_skipkn = set->ss.tsff.box_boundary_skipkn;
        mp->src->tsff->focused_thetamax = set->ss.tsff.focused_thetamax;
        mp->src->tsff->focused_fdist = set->ss.tsff.focused_fdist;
        mp->src->tsff->focused_nip = set->ss.tsff.focused_nip;
        mp->src->tsff->focused_mip = set->ss.tsff.focused_mip;
        mp->src->tsff->layered_epsilon = 1;//set->ss.tsff.layered_epsilon;
        mp->src->tsff->layered_mu = 1;//set->ss.tsff.layered_mu;
        mp->src->tsff->layered_sigma = 0;//set->ss.tsff.layered_sigma;
        mp->src->tsff->layered_sigast = 0;//set->ss.tsff.layered_sigast;

        fr = fopen(set->ss.tsff.source_filename, "r");
        if (fr == NULL) {
            fprintf(stderr, "Cannot load source file %s\n", set->ss.tsff.source_filename);
            return FALSE;
        }

        if (fscanf(fr, "%d", &ndata) != 1) {
            fprintf(stderr, "Error: cannot read number of values from file %s\n", set->ss.tsff.source_filename);
            fclose(fr);
            return 1;
        }
        if (ndata <= 0) {
            fprintf(stderr, "Error: zero source steps in %s\n", set->ss.tsff.source_filename);
            fclose(fr);
            return 1;
        }
        if (ndata < set->sc.nsteps)
            fprintf(stderr, "Warning: less data (%d) in \"%s\" than scheduled computation steps (%d) in TSFF source\n", ndata, set->ss.pnts[i].source_filename, set->sc.nsteps);
        mp->src->tsff->ndata = ndata;
        mp->src->tsff->layered_zpos = (gint *)g_malloc(ndata * sizeof(gint));
        mp->src->tsff->e = (gdouble *)g_malloc(ndata * sizeof(gdouble));

        for (k = 0; k < ndata; k++) {
            if (get_int(fr, &pos, "TSFF Source")) {
                fclose(fr);
                return 1;
            }
            if (get_double(fr, &ex, "TSFF Source")) {
                fclose(fr);
                return 1;
            }

            mp->src->tsff->layered_zpos[k] = pos;
            mp->src->tsff->e[k] = ex;//(ex-pex)/set->plan.dt;
        }
        if (estimate_lambda(mp->src->tsff->e, ndata, set->plan.dt, &lmin, &lmax, set->so.savespectrum)) {
            if (set->sc.verbose > 1)
                printf("Estimated wavelength range: (%g, %g) m  after file %s\n", lmin, lmax, set->ss.tsff.source_filename);
        } else if (set->sc.verbose > 1)
            printf("Automated detection of wavelength failed on file %s\n", set->ss.tsff.source_filename);

        mp->src->tsff->jp = (Sv1DPool ***)g_malloc(mp->src->tsff->focused_nip * sizeof(Sv1DPool **));
        mp->src->tsff->ia_theta = (gdouble **)g_malloc(mp->src->tsff->focused_nip * sizeof(gdouble *));
        mp->src->tsff->ia_phi = (gdouble **)g_malloc(mp->src->tsff->focused_nip * sizeof(gdouble *));
        mp->src->tsff->ia_psi = (gdouble **)g_malloc(mp->src->tsff->focused_nip * sizeof(gdouble *));
        mp->src->tsff->corr = (gdouble **)g_malloc(mp->src->tsff->focused_nip * sizeof(gdouble *));
        mp->src->tsff->an = (gdouble **)g_malloc(mp->src->tsff->focused_nip * sizeof(gdouble *));
        mp->src->tsff->bm = (gdouble **)g_malloc(mp->src->tsff->focused_nip * sizeof(gdouble *));

        xsn = (gdouble *)g_malloc(mp->src->tsff->focused_nip * sizeof(gdouble));
        wsn = (gdouble *)g_malloc(mp->src->tsff->focused_nip * sizeof(gdouble));
        xsm = (gdouble *)g_malloc(mp->src->tsff->focused_mip * sizeof(gdouble));
        wsm = (gdouble *)g_malloc(mp->src->tsff->focused_mip * sizeof(gdouble));

        get_gl_coeffs(xsm, wsm, mp->src->tsff->focused_mip);
        get_gl_coeffs(xsn, wsn, mp->src->tsff->focused_nip);

        for (i1 = 0; i1 < mp->src->tsff->focused_nip; i1++) {
            mp->src->tsff->jp[i1] = (Sv1DPool **)g_malloc(mp->src->tsff->focused_mip * sizeof(Sv1DPool *));
            mp->src->tsff->ia_theta[i1] = (gdouble *)g_malloc(mp->src->tsff->focused_mip * sizeof(gdouble));
            mp->src->tsff->ia_phi[i1] = (gdouble *)g_malloc(mp->src->tsff->focused_mip * sizeof(gdouble));
            mp->src->tsff->ia_psi[i1] = (gdouble *)g_malloc(mp->src->tsff->focused_mip * sizeof(gdouble));
            mp->src->tsff->corr[i1] = (gdouble *)g_malloc(mp->src->tsff->focused_mip * sizeof(gdouble));
            mp->src->tsff->an[i1] = (gdouble *)g_malloc(mp->src->tsff->focused_mip * sizeof(gdouble));
            mp->src->tsff->bm[i1] = (gdouble *)g_malloc(mp->src->tsff->focused_mip * sizeof(gdouble));

            for (im = 0; im < mp->src->tsff->focused_mip; im++) {
                /*// extended midpoint
                mp->src->tsff->theta[in][im] = (in + 0.5)*mp->src->tsff->thetamax/mp->src->tsff->nint;
                mp->src->tsff->phi[in][im] = mp->src->tsff->psi[in][im] = (im + 0.5)*2*G_PI/mp->src->tsff->mint;
                mp->src->tsff->an[in][im] = mp->src->tsff->thetamax/mp->src->tsff->nint;
                mp->src->tsff->bm[in][im] = 2.0*G_PI/mp->src->tsff->mint;
                */

                /* //gauss legendre
                mp->src->tsff->theta[in][im] = (xsn[in]+1)*mp->src->tsff->thetamax/2.0;
                mp->src->tsff->phi[in][im] =  mp->src->tsff->psi[in][im] = (xsm[im]+1)*G_PI;
                mp->src->tsff->an[in][im] = wsn[in]*mp->src->tsff->thetamax/2;
                mp->src->tsff->bm[in][im] = G_PI*wsm[im];*/

                //mixed - GL in theta but midpoint in phi (which is periodic)
                mp->src->tsff->ia_theta[i1][im] = (xsn[i1] + 1)*mp->src->tsff->focused_thetamax / 2.0;
                mp->src->tsff->ia_phi[i1][im] = (im + 0.5) * 2 * G_PI / mp->src->tsff->focused_mip;
                mp->src->tsff->ia_psi[i1][im] = mp->src->tsff->ia_phi[i1][im] + set->ss.tsff.focused_pol;//no add
                mp->src->tsff->an[i1][im] = wsn[i1] * mp->src->tsff->focused_thetamax / 2;
                mp->src->tsff->bm[i1][im] = 2.0*G_PI / mp->src->tsff->focused_mip;

                mp->src->tsff->corr[i1][im] = 1.0 / sqrt(pow(sin(mp->src->tsff->ia_theta[i1][im]), 4)
                                                         *(pow(cos(mp->src->tsff->ia_phi[i1][im]), 4) + pow(sin(mp->src->tsff->ia_phi[i1][im]), 4))
                                                         + pow(cos(mp->src->tsff->ia_theta[i1][im]), 4));

                mp->src->tsff->jp[i1][im] = sv_1dpool_new((gint)(sqrt((gdouble)(set->sp.xres*set->sp.xres + set->sp.yres*set->sp.yres + set->sp.zres*set->sp.zres)) + 250),
                                                          set->sp.dx / mp->src->tsff->corr[i1][im], set->plan.dt);

                sv_1dpool_set_source_pos(mp->src->tsff->jp[i1][im], 0);

                if (set->sc.verbose > 1) printf("set material: %g %g\n", mp->src->tsff->layered_epsilon, mp->src->tsff->layered_mu);
                sv_1dpool_set_material(mp->src->tsff->jp[i1][im], mp->src->tsff->layered_epsilon, mp->src->tsff->layered_mu);

                if (set->sc.verbose > 1) printf("angles: %g %g %g %d\n", mp->src->tsff->ia_theta[i1][im], mp->src->tsff->ia_phi[i1][im], mp->src->tsff->ia_psi[i1][im], (gint)(sqrt((gdouble)(set->sp.xres*set->sp.xres + set->sp.yres*set->sp.yres + set->sp.zres*set->sp.zres))) + 250);
            }
        }

        fclose(fr);
    }
    //if (set->ss.sf.source_filename) {
    if (set->ss.sf.valid) {
        mp->src->sf = (SvSourceSF *)g_malloc(sizeof(SvSourceSF));
        mp->src->sf->ia_theta = set->ss.sf.ia_theta;
        mp->src->sf->ia_phi = set->ss.sf.ia_phi;
        mp->src->sf->ia_psi = set->ss.sf.ia_psi;
        mp->src->sf->layered_epsilon = set->ss.sf.layered_epsilon;
        mp->src->sf->layered_mu = set->ss.sf.layered_mu;
        mp->src->sf->layered_sigma = set->ss.sf.layered_sigma;
        mp->src->sf->layered_sigast = set->ss.sf.layered_sigast;

        fr = fopen(set->ss.sf.source_filename, "r");
        if (fr == NULL) {
            fprintf(stderr, "Cannot load source file %s\n", set->ss.sf.source_filename);
            return FALSE;
        }
        if (fscanf(fr, "%d", &ndata) != 1) {
            fprintf(stderr, "Error: cannot read number of values from file %s\n", set->ss.sf.source_filename);
            fclose(fr);
            return 1;
        }

        if (ndata <= 0) {
            fprintf(stderr, "Error: zero source steps in %s\n", set->ss.sf.source_filename);
            fclose(fr);
            return 1;
        }
        if (ndata < set->sc.nsteps)
            fprintf(stderr, "Warning: less data (%d) in \"%s\" than scheduled computation steps (%d) in SF source\n", ndata, set->ss.pnts[i].source_filename, set->sc.nsteps);
        mp->src->sf->ndata = ndata;
        mp->src->sf->layered_zpos = (gint *)g_malloc(ndata * sizeof(gint));
        mp->src->sf->e = (gdouble *)g_malloc(ndata * sizeof(gdouble));

        for (k = 0; k < ndata; k++) {
            if (get_int(fr, &pos, "TSF Source")) {
                fclose(fr);
                return 1;
            }
            if (get_double(fr, &ex, "TSF Source")) {
                fclose(fr);
                return 1;
            }

            mp->src->sf->layered_zpos[k] = pos;
            mp->src->sf->e[k] = ex;
        }
        if (estimate_lambda(mp->src->sf->e, ndata, set->plan.dt, &lmin, &lmax, set->so.savespectrum)) {
            if (set->sc.verbose > 1)
                printf("Estimated wavelength range: (%g, %g) m  after file %s\n", lmin, lmax, set->ss.sf.source_filename);
        } else if (set->sc.verbose > 1)
            printf("Automated detection of wavelength failed on file %s\n", set->ss.sf.source_filename);

        mp->src->sf->jp = sv_1dpool_new((gint)sqrt((gdouble)(set->sp.xres*set->sp.xres + set->sp.yres*set->sp.yres + set->sp.zres*set->sp.zres) + 10), set->sp.dx, set->plan.dt);

        sv_1dpool_set_source_pos(mp->src->sf->jp, 0);
        sv_1dpool_set_material(mp->src->sf->jp, mp->src->sf->layered_epsilon, mp->src->sf->layered_mu);

        fclose(fr);
    }

    if (set->ss.lambda_min == -1) {
        if (lmin > 0 && lmin < 1e6 && lmax > 0 && lmax < 1e6 && lmin <= lmax) {
            set->ss.lambda_min = lmin;
            set->ss.lambda_max = lmax;
            set->ss.lambda_center = (lmin + lmax) / 2.0;
            if (set->sc.verbose > 1)
                printf("Total estimated source wavelength: range (%g, %g) m, center at %g m\n", set->ss.lambda_min, set->ss.lambda_max, set->ss.lambda_center);
        } else {
            if (set->sc.verbose > 1)
                printf("Warning: no automatically set wavelength range. Predefined materials optical properties may not work.\n");
        }
    } else {
        if (set->sc.verbose > 1)
            printf("Automatically detected wavelength overrriden by user: range (%g, %g) m, center at %g m\n", set->ss.lambda_min, set->ss.lambda_max, set->ss.lambda_center);
    }

    return 0;
}

static GLfloat
dist(GLfloat *a, GLfloat *b)
{
    return (GLfloat)sqrt((a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1]) + (a[2] - b[2])*(a[2] - b[2]));
}

static GLfloat
dot(GLfloat *a, GLfloat *b, GLfloat *c)
{
    return (b[0] - a[0])*(c[0] - b[0]) + (b[1] - a[1])*(c[1] - b[1]) + (b[2] - a[2])*(c[2] - b[2]);
}

static GLfloat
mag(GLfloat aa, GLfloat ab, GLfloat ac)
{
    return (GLfloat)sqrt(aa*aa + ab*ab + ac*ac);
}

static GLfloat
ltop(GLfloat *a, GLfloat *b, GLfloat *c)
{
    GLfloat a1 = a[0] - b[0];
    GLfloat a2 = a[1] - b[1];
    GLfloat a3 = a[2] - b[2];
    GLfloat b1 = a[0] - c[0];
    GLfloat b2 = a[1] - c[1];
    GLfloat b3 = a[2] - c[2];

    return mag(a2*b3 - a3*b2, a3*b1 - a1*b3, a1*b2 - a2*b1) / mag(c[0] - b[0], c[1] - b[1], c[2] - b[2]);
}

static GLfloat
linedist(GLfloat *a, GLfloat *b, GLfloat *c)
{
    // gdouble dist, dot1, dot2;

    if (dot(b, c, a) > 0 || dot(c, b, a) > 0)
        return 1e6;
    else
        return (GLfloat)fabs(ltop(a, b, c));

    /*
    dist = ltop(a,b,c);
    dot1 = dot(b,c,a);
    if (dot1 > 0) return 1e6;
    dot2 = dot(c,b,a);
    if (dot2 > 0) return 1e6;
    return fabs(dist);
    */
}

static GLfloat
dett(GLfloat a11, GLfloat a12, GLfloat a13,
     GLfloat a21, GLfloat a22, GLfloat a23,
     GLfloat a31, GLfloat a32, GLfloat a33)
{
    return a11*a22*a33 - a11*a23*a32 + a12*a23*a31 - a12*a21*a33 + a13*a21*a32 - a13*a22*a31;
}

static GLfloat
determinant(GLfloat *a, GLfloat *b, GLfloat *c, GLfloat *d)
{
    GLfloat da, db, dc, dd;

    da = dett(b[1], b[2], b[3], c[1], c[2], c[3], d[1], d[2], d[3]);
    db = dett(b[0], b[2], b[3], c[0], c[2], c[3], d[0], d[2], d[3]);
    dc = dett(b[0], b[1], b[3], c[0], c[1], c[3], d[0], d[1], d[3]);
    dd = dett(b[0], b[1], b[2], c[0], c[1], c[2], d[0], d[1], d[2]);

    return a[0] * da - a[1] * db + a[2] * dc - a[3] * dd;
}

// is_in_...      x,y,z - actual coordinates of a point

static gboolean
is_in_sphere(GLfloat x, GLfloat y, GLfloat z, GLfloat pnt1[3], GLfloat radius)
{
    if (((x - pnt1[0])*(x - pnt1[0]) + (y - pnt1[1])*(y - pnt1[1]) + (z - pnt1[2])*(z - pnt1[2])) <= radius*radius)
        return 1;
    else
        return 0;
}

static gboolean
is_in_cylinder(GLfloat x, GLfloat y, GLfloat z, GLfloat pnt1[3], GLfloat pnt2[3], GLfloat radius)
{
    GLfloat pntx[3];
    pntx[0] = x;
    pntx[1] = y;
    pntx[2] = z;

    if (linedist(pntx, pnt1, pnt2) < radius)
        return 1;
    else
        return 0;
}

static gboolean
is_in_cone(GLfloat x, GLfloat y, GLfloat z, GLfloat pnt1[3], GLfloat pnt2[3], GLfloat radius)
{
    GLfloat pntx[3];
    GLfloat length, pos;
    pntx[0] = x;
    pntx[1] = y;
    pntx[2] = z;
    length = dist(pnt1, pnt2);
    pos = dist(pntx, pnt1);

    if (linedist(pntx, pnt1, pnt2) < (pos*radius / length))
        return 1;

    return 0;
}

static gboolean
is_in_rcone(GLfloat x, GLfloat y, GLfloat z, GLfloat pnt1[3], GLfloat pnt2[3], GLfloat radius1, GLfloat radius2)
{
    GLfloat pntx[3];
    GLfloat length, pos;
    GLfloat pnt1to2_nrm[3], pnt1tox[3], pntx_proj[3];
    GLfloat proj_coeff;
    pntx[0] = x;
    pntx[1] = y;
    pntx[2] = z;
    length = dist(pnt1, pnt2);

    // normovany smerovy vektor usecky 1-2
    pnt1to2_nrm[0] = (pnt2[0] - pnt1[0]) / length;
    pnt1to2_nrm[1] = (pnt2[1] - pnt1[1]) / length;
    pnt1to2_nrm[2] = (pnt2[2] - pnt1[2]) / length;

    pnt1tox[0] = x - pnt1[0];
    pnt1tox[1] = y - pnt1[1];
    pnt1tox[2] = z - pnt1[2];

    proj_coeff = pnt1tox[0] * pnt1to2_nrm[0] + pnt1tox[1] * pnt1to2_nrm[1] + pnt1tox[2] * pnt1to2_nrm[2];

    // projekce bodu pntx na usecku pnt1-pnt2 = pnt1 + projekce vektoru 1tox na 1to2
    pntx_proj[0] = pnt1[0] + proj_coeff * pnt1to2_nrm[0];
    pntx_proj[1] = pnt1[1] + proj_coeff * pnt1to2_nrm[1];
    pntx_proj[2] = pnt1[2] + proj_coeff * pnt1to2_nrm[2];

    pos = proj_coeff;

    if (pos >= 0 && pos <= length && dist(pntx, pntx_proj) <= (radius1 + pos*(radius2 - radius1) / length))
        return 1;

    return 0;
}

static gboolean
is_in_tetrahedron(GLfloat x, GLfloat y, GLfloat z, GLfloat *p1, GLfloat *p2, GLfloat *p3, GLfloat *p4)
{
    GLfloat d0[4], d1[4], d2[4], d3[4], dn[4];
    GLfloat rdet[5];

    //printf("Th: %g %g %g,   %g %g %g,   %g %g %g,   %g %g %g\n", p1[0], p1[1], p1[2], p1[3], p2[0], p2[1], p2[2], p2[3], p3[0], p3[1], p3[2], p3[3], p4[0], p4[1], p4[2], p4[3]);

    /*det 0*/
    d0[0] = p1[0]; d0[1] = p1[1]; d0[2] = p1[2]; d0[3] = 1;
    d1[0] = p2[0]; d1[1] = p2[1]; d1[2] = p2[2]; d1[3] = 1;
    d2[0] = p3[0]; d2[1] = p3[1]; d2[2] = p3[2]; d2[3] = 1;
    d3[0] = p4[0]; d3[1] = p4[1]; d3[2] = p4[2]; d3[3] = 1;
    dn[0] = x; dn[1] = y; dn[2] = z; dn[3] = 1;

    rdet[0] = determinant(d0, d1, d2, d3);
    rdet[1] = determinant(dn, d1, d2, d3);
    rdet[2] = determinant(d0, dn, d2, d3);
    rdet[3] = determinant(d0, d1, dn, d3);
    rdet[4] = determinant(d0, d1, d2, dn);

    if (rdet[0] == 0)
        fprintf(stderr, "Warning: degenerated tetrahedron\n");
    if ((rdet[0] >= 0 && rdet[1] >= 0 && rdet[2] >= 0 && rdet[3] >= 0 && rdet[4] >= 0) ||
        (rdet[0] <= 0 && rdet[1] <= 0 && rdet[2] <= 0 && rdet[3] <= 0 && rdet[4] <= 0))
        return 1;
    else
        return 0;
}

static void
sv_pool_set_data(SvPool *mp, SvMatProp mat, SvSet *set, gint i, gint j, gint k)
{
    if (mat.type == SV_MAT_LINEAR) {
        if (set->plan.matmode == SV_MATMODE_ELECTRIC || set->plan.matmode == SV_MATMODE_FULL) {
            mp->d->epsilon->data[i][j][k] = (gfloat)(mat.epsilon * EPSILON_0);
            mp->d->sigma->data[i][j][k] = (gfloat)mat.sigma;
        }
        if (set->plan.matmode == SV_MATMODE_MAGNETIC || set->plan.matmode == SV_MATMODE_FULL) {
            mp->d->mu->data[i][j][k] = (gfloat)(mat.mu * MU_0);
            mp->d->sigast->data[i][j][k] = (gfloat)(mat.sigast);
        }
        if (mp->d->mat)
            mp->d->mat->data[i][j][k] = 0; // valgrind: mp->d->mat neinicializovano; pool.c : sv_pool_new
    } else
        mp->d->mat->data[i][j][k] = mat.pos;
}
static void
sv_yee_set_data(SvYeeData *d, SvMatProp mat, SvSet *set, gint i, gint j, gint k)
{
    if (i<0 || j<0 || k<0 || i>=d->xres || j>=d->yres || k>=d->zres) return;

    if (mat.type == SV_MAT_LINEAR) {
        if (set->plan.matmode == SV_MATMODE_ELECTRIC || set->plan.matmode == SV_MATMODE_FULL) {
            d->epsilon->data[i][j][k] = (gfloat)(mat.epsilon * EPSILON_0);
            d->sigma->data[i][j][k] = (gfloat)mat.sigma;
        }
        if (set->plan.matmode == SV_MATMODE_MAGNETIC || set->plan.matmode == SV_MATMODE_FULL) {
            d->mu->data[i][j][k] = (gfloat)(mat.mu * MU_0);
            d->sigast->data[i][j][k] = (gfloat)(mat.sigast);
        }
        if (d->mat)
            d->mat->data[i][j][k] = 0;
    } else
        d->mat->data[i][j][k] = mat.pos;
}

/*
get all the relevant material information at single point:
- matindex: -1 if there are no tabulated materials at all or we are not at tabulated material point,
             otherwise it is the position in tabulated material list, enough unique identifier if there are tabulated materials
         (in this case, 0 stands for vacuum)
- epsilon, sigma, mu, sigast: filled only if matindex == -1: this means that these four values are the only material parameters available
- note that epsilon, mu are relative values

- for simple pickoff of epsilon, sigma, mu or sigast, use sv_pool_get_epsilon, sv_pool_get_sigma, etc.

*/

static void
sv_pool_get_matinfo(SvPool *mp, SvSet *set, gint i, gint j, gint k, gint *matindex, gdouble *epsilon, gdouble *mu, gdouble *sigma, gdouble *sigast)
{

    if (mp->d->mat) { // there are tabulated materials (--> vacuum at position 0)
        *matindex = mp->d->mat->data[i][j][k];
    } else { // no tabulated materials -> SV_MAT_LINEAR
        *matindex = -1;

        if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_ELECTRIC) { // no tabulated material, allocated local properties
            *epsilon = mp->d->epsilon->data[i][j][k] / EPSILON_0; // mp->d->epsilon stores absolute values
            *sigma = mp->d->sigma->data[i][j][k];
        } else { // no material info at all, this is vacuum
            *epsilon = 1;
            *sigma = 0;
        }

        if (set->plan.matmode == SV_MATMODE_FULL || set->plan.matmode == SV_MATMODE_MAGNETIC) { // no tabulated material, allocated local properties
            *mu = mp->d->mu->data[i][j][k] / MU_0; // mp->d->mu stores absolute values
            *sigast = mp->d->sigast->data[i][j][k];
        } else { // no material info at all, this is vacuum
            *mu = 1;
            *sigast = 0;
        }
    }
}

static gboolean
almost_eq(gdouble a, gdouble b)
{
    const gdouble TOL = 1e-3;

    if (fabs(a) >= 1.0 && fabs(b) >= 1.0) // only use relative comparison for "sufficiently large" values
        return (fabs(a / b - 1.0) < TOL);
    else
        return (fabs(a - b) < TOL);
}

static gint
sv_pool_get_matindex(SvPool *mp, SvMatProp mat) // try to find the first material with given properties; return -1 if not found
{
    gint i, matindex;
    gboolean found;
    SvMatProp tmpmat;

    matindex = -1;
    found = FALSE;

    for (i = 0; (i < mp->nmat) && (found == FALSE); i++) {
        tmpmat = mp->mats[i];
        if (tmpmat.type == mat.type) {
            if (mat.type == SV_MAT_PEC ||
                (almost_eq(tmpmat.epsilon, mat.epsilon) && almost_eq(tmpmat.mu, mat.mu) && almost_eq(tmpmat.sigma, mat.sigma) && almost_eq(tmpmat.sigast, mat.sigast) &&
                 almost_eq(tmpmat.drude_omega_p, mat.drude_omega_p) && almost_eq(tmpmat.drude_nu, mat.drude_nu) &&
                 almost_eq(tmpmat.cp3_a[0], mat.cp3_a[0]) && almost_eq(tmpmat.cp3_a[1], mat.cp3_a[1]) && almost_eq(tmpmat.cp3_phi[2], mat.cp3_a[2]) &&
                 almost_eq(tmpmat.cp3_phi[0], mat.cp3_phi[0]) && almost_eq(tmpmat.cp3_phi[1], mat.cp3_phi[1]) && almost_eq(tmpmat.cp3_phi[2], mat.cp3_phi[2]) &&
                 almost_eq(tmpmat.cp3_omega[0], mat.cp3_omega[0]) && almost_eq(tmpmat.cp3_omega[1], mat.cp3_omega[1]) && almost_eq(tmpmat.cp3_omega[2], mat.cp3_omega[2]) &&
                 almost_eq(tmpmat.cp3_gamma[0], mat.cp3_gamma[0]) && almost_eq(tmpmat.cp3_gamma[1], mat.cp3_gamma[1]) && almost_eq(tmpmat.cp3_gamma[2], mat.cp3_gamma[2]))) {

                // this is too tolerant, as each quantity has a different order of magnitude (esp. note that perm***ies are absolute, i.e. very small)
                /*
                  fabs(tmpmat.epsilon - mat.epsilon) + fabs(tmpmat.mu - mat.mu) + fabs(tmpmat.sigma - mat.sigma) + fabs(tmpmat.sigast - mat.sigast ) +
                  fabs(tmpmat.drude_omega_p - mat.drude_omega_p) + fabs(tmpmat.drude_nu - mat.drude_nu) +
                  fabs(tmpmat.cp3_a[0] - mat.cp3_a[0]) + fabs(tmpmat.cp3_a[1] - mat.cp3_a[1]) + fabs(tmpmat.cp3_a[2] - mat.cp3_a[2]) +
                  fabs(tmpmat.cp3_phi[0] - mat.cp3_phi[0]) + fabs(tmpmat.cp3_phi[1] - mat.cp3_phi[1]) + fabs(tmpmat.cp3_phi[2] - mat.cp3_phi[2]) +
                  fabs(tmpmat.cp3_omega[0] - mat.cp3_omega[0]) + fabs(tmpmat.cp3_omega[1] - mat.cp3_omega[1]) + fabs(tmpmat.cp3_omega[2] - mat.cp3_omega[2]) +
                  fabs(tmpmat.cp3_gamma[0] - mat.cp3_gamma[0]) + fabs(tmpmat.cp3_gamma[1] - mat.cp3_gamma[1]) + fabs(tmpmat.cp3_gamma[2] - mat.cp3_gamma[2])
                  < TOL) {
                */
                matindex = i;
                found = TRUE;
            } // if
        } // if
    } // for

    return matindex; // OK, vraci spravne hodnoty
}

#if(0)
/* given material description (see sv_pool_matinfo) it checks if the same material is at given point (i,j,k).*/
static gboolean
sv_pool_is_there_material(SvPool *mp, SvSet *set, gint i, gint j, gint k, gint matindex, gdouble epsilon, gdouble mu, gdouble sigma, gdouble sigast)
{
    int local_matindex;
    double local_epsilon, local_mu, local_sigma, local_sigast;
    double diff = 1e-10;

    sv_pool_get_matinfo(mp, set, i, j, k, &local_matindex, &local_epsilon, &local_mu, &local_sigma, &local_sigast);

    if (matindex >= 0) {
        if (matindex == local_matindex)
            return TRUE;
    } else {
        if ((fabs(epsilon - local_epsilon) < diff) && (fabs(mu - local_mu) < diff) && (fabs(sigma - local_sigma) < diff) && (fabs(sigast - local_sigast) < diff))
            return TRUE;
    }
    return FALSE;

}
#endif

void
sv_pool_add_material(SvPool *mp, SvMatProp *mat, gboolean multme)
{
    gint i;

    // original version, remove if sv_pool_get_matindex works fine
    /*
      for (i=1; i < mp->nmat; i++) {
      if (((mat->type==SV_MAT_LINEAR || mat->type==SV_MAT_LINTAB) && (mp->mats[i].type == mat->type)
      && (fabs(mat->epsilon - mp->mats[i].epsilon) < diff))
      ||
      ((mat->type==SV_MAT_DRUDE || mat->type==SV_MAT_CP3 || mat->type==SV_MAT_CP || mat->type==SV_MAT_ADE || mat->type==SV_MAT_PLRC) && (mp->mats[i].type == mat->type)
      && (fabs(mat->epsilon - mp->mats[i].epsilon) < diff))
      || (mat->type==SV_MAT_PEC && mp->mats[i].type==SV_MAT_PEC)
      ) {
      printf("We already have this material as No %d\n", i);
      mat->pos = i;
      neww=FALSE;
      }
      }
    */

    /*zero index material it is used to determine whether use tabulated material at given point or no*/
    // RS - set it to be like vacuum
    if (mp->nmat == 0) {
        mp->mats = (SvMatProp *)g_realloc(mp->mats, sizeof(SvMatProp));
        mp->mats[0] = mat_vacuum_lintab;
        /*
            mp->mats[0].epsilon = 1;
            mp->mats[0].mu = 1;
            mp->mats[0].sigma = 0;
            mp->mats[0].sigast = 0;
            mp->mats[0].type = SV_MAT_LINTAB; // puvodne 0
        mp->mats[0].pos = 0;
        */
        mp->nmat++;
    }

    mp->mats = (SvMatProp *)g_realloc(mp->mats, (mp->nmat + 1) * sizeof(SvMatProp));
    mp->mats[mp->nmat].epsilon = mat->epsilon;
    mp->mats[mp->nmat].mu = mat->mu;
    mp->mats[mp->nmat].sigma = mat->sigma;
    mp->mats[mp->nmat].sigast = mat->sigast;
    mp->mats[mp->nmat].drude_omega_p = mat->drude_omega_p;
    mp->mats[mp->nmat].drude_nu = mat->drude_nu;
    mp->mats[mp->nmat].type = mat->type;

    if (multme) { // workaround for cropping after everything is loaded
        mp->mats[mp->nmat].epsilon *= EPSILON_0;
        mp->mats[mp->nmat].mu *= MU_0;
    }

    for (i = 0; i < 3; i++) {
        mp->mats[mp->nmat].cp3_a[i] = mat->cp3_a[i];
        mp->mats[mp->nmat].cp3_phi[i] = mat->cp3_phi[i];
        mp->mats[mp->nmat].cp3_omega[i] = mat->cp3_omega[i];
        mp->mats[mp->nmat].cp3_gamma[i] = mat->cp3_gamma[i];
    }

    mp->mats[mp->nmat].ade_a0 = mat->ade_a0;
    mp->mats[mp->nmat].ade_a1 = mat->ade_a1;
    mp->mats[mp->nmat].ade_a2 = mat->ade_a2;

    mp->mats[mp->nmat].ade_bp0[0] = mat->ade_bp0[0];
    mp->mats[mp->nmat].ade_bp0[1] = mat->ade_bp0[1];
    mp->mats[mp->nmat].ade_bp1[0] = mat->ade_bp1[0];
    mp->mats[mp->nmat].ade_bp1[1] = mat->ade_bp1[1];
    mp->mats[mp->nmat].ade_bp2[0] = mat->ade_bp2[0];
    mp->mats[mp->nmat].ade_bp2[1] = mat->ade_bp2[1];
    mp->mats[mp->nmat].ade_bp3[0] = mat->ade_bp3[0];
    mp->mats[mp->nmat].ade_bp3[1] = mat->ade_bp3[1];
    mp->mats[mp->nmat].ade_bp4[0] = mat->ade_bp4[0];
    mp->mats[mp->nmat].ade_bp4[1] = mat->ade_bp4[1];

    mp->mats[mp->nmat].ade_c0 = mat->ade_c0;
    mp->mats[mp->nmat].ade_c1 = mat->ade_c1;
    mp->mats[mp->nmat].ade_c2 = mat->ade_c2;
    mp->mats[mp->nmat].ade_c3 = mat->ade_c3;

    mp->mats[mp->nmat].plrc_d_chi = mat->plrc_d_chi;
    mp->mats[mp->nmat].plrc_d_xi = mat->plrc_d_xi;
    mp->mats[mp->nmat].plrc_d_dchi = mat->plrc_d_dchi;
    mp->mats[mp->nmat].plrc_d_dxi = mat->plrc_d_dxi;

    mp->mats[mp->nmat].plrc_p_chi[0] = mat->plrc_p_chi[0];
    mp->mats[mp->nmat].plrc_p_xi[0] = mat->plrc_p_xi[0];
    mp->mats[mp->nmat].plrc_p_dchir[0] = mat->plrc_p_dchir[0];
    mp->mats[mp->nmat].plrc_p_dxir[0] = mat->plrc_p_dxir[0];
    mp->mats[mp->nmat].plrc_p_dchii[0] = mat->plrc_p_dchii[0];
    mp->mats[mp->nmat].plrc_p_dxii[0] = mat->plrc_p_dxii[0];

    mp->mats[mp->nmat].plrc_p_chi[1] = mat->plrc_p_chi[1];
    mp->mats[mp->nmat].plrc_p_xi[1] = mat->plrc_p_xi[1];
    mp->mats[mp->nmat].plrc_p_dchir[1] = mat->plrc_p_dchir[1];
    mp->mats[mp->nmat].plrc_p_dxir[1] = mat->plrc_p_dxir[1];
    mp->mats[mp->nmat].plrc_p_dchii[1] = mat->plrc_p_dchii[1];
    mp->mats[mp->nmat].plrc_p_dxii[1] = mat->plrc_p_dxii[1];

    mp->mats[mp->nmat].pos = mp->nmat;

    switch (mat->type) {
        case SV_MAT_LINTAB:
            printf("material %d of type %d added with epsilon %f, sigma %f, mu %f, sigast %f\n", mp->nmat, mp->mats[mp->nmat].type, mat->epsilon, mat->sigma, mat->mu, mat->sigast);
            break;

        case SV_MAT_DRUDE:
            printf("material %d of type %d added with epsilon %f, omega_p %f, nu %f\n", mp->nmat, mp->mats[mp->nmat].type, mat->epsilon, mat->drude_omega_p, mat->drude_nu);
            break;

        case SV_MAT_CP3:
            printf("material %d of type %d added with epsilon %f, sigma %f,\n                                a0 %f, phi0 %f, omega0 %f, gamma0 %f,\n                                a1 %f, phi1 %f, omega1 %f, gamma1 %f,\n                                a2 %f, phi2 %f, omega2 %f, gamma2 %f\n", mp->nmat, mp->mats[mp->nmat].type, mat->epsilon, mat->sigma, mat->cp3_a[0], mat->cp3_phi[0], mat->cp3_omega[0], mat->cp3_gamma[0], mat->cp3_a[1], mat->cp3_phi[1], mat->cp3_omega[1], mat->cp3_gamma[1], mat->cp3_a[2], mat->cp3_phi[2], mat->cp3_omega[2], mat->cp3_gamma[2]);
            break;

        case SV_MAT_CP:
            printf("material %d of type %d added with epsilon %f, omega_p %f, gamma %f,\n                                a0 %f, phi0 %f, omega0 %f, gamma0 %f,\n                                a1 %f, phi1 %f, omega1 %f, gamma1 %f\n", mp->nmat, mp->mats[mp->nmat].type, mat->epsilon, mat->drude_omega_p, mat->drude_nu, mat->cp3_a[0], mat->cp3_phi[0], mat->cp3_omega[0], mat->cp3_gamma[0], mat->cp3_a[1], mat->cp3_phi[1], mat->cp3_omega[1], mat->cp3_gamma[1]); // !!!!!!!!!!!!!!!
            break;

        case SV_MAT_ADE:
            printf("material %d of type %d added with epsilon %f, sigma %f, omega_p %f, gamma%f,\n                                a0 %f, phi0 %f, omega0 %f, gamma0 %f,\n                                a1 %f, phi1 %f, omega1 %f, gamma1 %f\n", mp->nmat, mp->mats[mp->nmat].type, mat->epsilon, mat->sigma, mat->drude_omega_p, mat->drude_nu, mat->cp3_a[0], mat->cp3_phi[0], mat->cp3_omega[0], mat->cp3_gamma[0], mat->cp3_a[1], mat->cp3_phi[1], mat->cp3_omega[1], mat->cp3_gamma[1]); // !!!!!!!!
            /*
            calculate_ade_a(&mat, mat.drude_omega_p, mat.drude_nu, data->set->plan.dt);
            calculate_ade_bp(0, &mat, mat.cp3_a[0], mat.cp3_phi[0], mat.cp3_omega[0], mat.cp3_gamma[0], data->set->plan.dt);
            calculate_ade_bp(1, &mat, mat.cp3_a[1], mat.cp3_phi[1], mat.cp3_omega[1], mat.cp3_gamma[1], data->set->plan.dt);
            calculate_ade_c(&mat, mat.epsilon, mat.sigma, data->set->plan.dt);
            */
            break;

        case SV_MAT_PLRC: // NOT WORKING?
            printf("material %d of type %d added with epsilon %f, sigma %f, omega_p %f, gamma %f,\n                                a0 %f, phi0 %f, omega0 %f, gamma0 %f,\n                                a1 %f, phi1 %f, omega1 %f, gamma1 %f\n", mp->nmat, mp->mats[mp->nmat].type, mat->epsilon, mat->sigma, mat->drude_omega_p, mat->drude_nu, mat->cp3_a[0], mat->cp3_phi[0], mat->cp3_omega[0], mat->cp3_gamma[0], mat->cp3_a[1], mat->cp3_phi[1], mat->cp3_omega[1], mat->cp3_gamma[1]); // !!!!!!!!
            /*
            calculate_plrc_d(&mat, mat.drude_omega_p, mat.drude_nu, data->set->plan.dt);
            calculate_plrc_p(0, &mat, data->set->plan.dt);
            calculate_plrc_p(1, &mat, data->set->plan.dt);
            */
            break;

        case SV_MAT_PEC:
            printf("material %d of type %d added\n", mp->nmat, mp->mats[mp->nmat].type);
            break;

        default: // DATABASE -- TODO;
            break;
    }


    //printf("material %d of type %d added with eps %f, sigma %f, mu %f, sigast %f\n", mp->nmat, mp->mats[mp->nmat].type, mat->epsilon, mat->sigma, mat->mu, mat->sigast);
    // TODO rework according to material type

    mat->pos = mp->nmat;

    mp->nmat++;
}

static gint
submeshcmp(Submesh *sm1, Submesh *sm2)
{
    if (sm1->coords.i < sm2->coords.i)
        return -1;
    else if (sm1->coords.i == sm2->coords.i) {
        if (sm1->coords.j < sm2->coords.j)
            return -1;
        else if (sm1->coords.j == sm2->coords.j) {
            if (sm1->coords.k < sm2->coords.k)
                return -1;
            else if (sm1->coords.k == sm2->coords.k)
                return 0;
            else
                return +1; // sm1->i == sm2->i, sm1->j == sm2->j, sm1->k > sm2->k
        } else
            return +1; // sm1->i == sm2->i, sm1->j > sm2->j
    } else
        return +1; // sm1->i > sm2->i
}

static Submesh*
init_submesh(gint i, gint j, gint k, GTree *submeshes, SvPool *mp, SvSet *set)
{
    Submesh *sm;
    SvMatProp mat;
    SmCoords sm_coords = {i, j, k};
    gint l, m, n;

    sm = g_tree_lookup(submeshes, &sm_coords); // try to find the submesh among the mixed cells (and fill-in with respective data)
    if (sm == NULL) { // if the cell is not mixed, initialize working submesh using current cell values
        sm = g_malloc(sizeof(Submesh));
        sm->coords.i = i;
        sm->coords.j = j;
        sm->coords.k = k;
        sm->cells_written = 0;
        sm->was_mixed = FALSE;

        mat = mat_vacuum_linear; // init all material values with something harmless;
        sv_pool_get_matinfo(mp, set, i, j, k, &mat.pos, &mat.epsilon, &mat.mu, &mat.sigma, &mat.sigast); // determine background material
        if (mat.pos != -1) // tabulated material found
            mat = mp->mats[mat.pos];

        for (l = 0; l < 3; l++) {
            for (m = 0; m < 3; m++) {
                for (n = 0; n < 3; n++)
                    sm->material[l][m][n] = mat;
            } // m
        } // l
    } else { // cell mixed, note it
        sm->was_mixed = TRUE;
    }

    return sm;
}

static void
add_mat_to_submesh_sphere(gint i, gint j, gint k, GLfloat pnt1[3], GLfloat radius, Submesh *sm, SvMatProp mat)
{
    const GLfloat a = 1.0f / 3; // size of a little box
    gint l, m, n;
    sm->cells_written = 0;
    for (l = 0; l < 3; l++) {
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                if (is_in_sphere(i + (0.5f + l)*a, j + (0.5f + m)*a, k + (0.5f + n)*a, pnt1, radius)) {
                    sm->cells_written++;
                    sm->material[l][m][n] = mat;
                }
            } // n
        } // m
    } // l
}

static void
add_mat_to_submesh_cylinder(gint i, gint j, gint k, GLfloat pnt1[3], GLfloat pnt2[3], GLfloat radius, Submesh *sm, SvMatProp mat)
{
    const GLfloat a = 1.0f / 3; // size of a little box
    gint l, m, n;
    sm->cells_written = 0;
    for (l = 0; l < 3; l++) {
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                if (is_in_cylinder(i + (0.5f + l)*a, j + (0.5f + m)*a, k + (0.5f + n)*a, pnt1, pnt2, radius)) {
                    sm->cells_written++;
                    sm->material[l][m][n] = mat;
                }
            } // n
        } // m
    } // l
}

static void
add_mat_to_submesh_cone(gint i, gint j, gint k, GLfloat pnt1[3], GLfloat pnt2[3], GLfloat radius, Submesh *sm, SvMatProp mat)
{
    const GLfloat a = 1.0f / 3; // size of a little box
    gint l, m, n;
    sm->cells_written = 0;
    for (l = 0; l < 3; l++) {
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                if (is_in_cone(i + (0.5f + l)*a, j + (0.5f + m)*a, k + (0.5f + n)*a, pnt1, pnt2, radius)) {
                    sm->cells_written++;
                    sm->material[l][m][n] = mat;
                }
            } // n
        } // m
    } // l
}

static void
add_mat_to_submesh_rcone(gint i, gint j, gint k, GLfloat pnt1[3], GLfloat pnt2[3], GLfloat radius1, GLfloat radius2, Submesh *sm, SvMatProp mat)
{
    const GLfloat a = 1.0f / 3; // size of a little box
    gint l, m, n;
    sm->cells_written = 0;
    for (l = 0; l < 3; l++) {
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                if (is_in_rcone(i + (0.5f + l)*a, j + (0.5f + m)*a, k + (0.5f + n)*a, pnt1, pnt2, radius1, radius2)) {
                    sm->cells_written++;
                    sm->material[l][m][n] = mat;
                }
            } // n
        } // m
    } // l
}

static void
add_mat_to_submesh_tetrahedron(gint i, gint j, gint k, GLfloat pnt1[3], GLfloat pnt2[3], GLfloat pnt3[3], GLfloat pnt4[3], Submesh *sm, SvMatProp mat)
{
    const GLfloat a = 1.0f / 3; // size of a little box
    gint l, m, n;
    sm->cells_written = 0;
    for (l = 0; l < 3; l++) {
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                if (is_in_tetrahedron(i + (0.5f + l)*a, j + (0.5f + m)*a, k + (0.5f + n)*a, pnt1, pnt2, pnt3, pnt4)) {
                    sm->cells_written++;
                    sm->material[l][m][n] = mat;
                }
            } // n
        } // m
    } // l
}

static void
check_and_note_mixed(GTree *submeshes, Submesh *sm, SvPool *mp, SvSet *set)
{
    SvMatProp smmat;
    gint l, m, n;
    gint smmat_cnt = 0;

    if (sm->cells_written == 0) { // no material added in previous step
        if (sm->was_mixed == FALSE) // not in the tree -> just free it
            g_free(sm);
    }

    else if (sm->cells_written == 27) { // whole cell overwritten -> set data, remove from tree if needed
        sv_pool_set_data(mp, sm->material[0][0][0], set, sm->coords.i, sm->coords.j, sm->coords.k); // same material in all cells
        if (sm->was_mixed)
            g_tree_remove(submeshes, &sm->coords);
        g_free(sm);
    }

    else { // some material was added, not covering the whole cell

        // check whether one material covers the whole cell
        smmat = sm->material[0][0][0];
        for (l = 0; l < 3; l++) {
            for (m = 0; m < 3; m++) {
                for (n = 0; n < 3; n++) {
                    // equality test might be better for SV_MAT_LINEAR; at present only epsilon is compared
                    if ((smmat.pos >= 0 && sm->material[l][m][n].pos == smmat.pos) || (smmat.pos == -1 && sm->material[l][m][n].epsilon == smmat.epsilon))
                        smmat_cnt++;
                }
            }
        }

        if (smmat_cnt == 27) { // whole cell filled with the same material -> set data, remove from tree if needed
            sv_pool_set_data(mp, sm->material[0][0][0], set, sm->coords.i, sm->coords.j, sm->coords.k); // same material in all cells
            if (sm->was_mixed) // found in the tree --> remove it from the tree
                g_tree_remove(submeshes, &sm->coords);
            g_free(sm);
        } else { // mix of materials in the cell
            if (sm->was_mixed == FALSE) // not yet in the tree -> add it to the tree;
                g_tree_insert(submeshes, &sm->coords, sm);
        }
    }
}

static void
calculate_ade_a(SvMatProp *mat, gdouble omegap, gdouble nu, gdouble dt)
{
    mat->ade_a0 = (nu*dt - 2) / (nu*dt + 2);
    mat->ade_a1 = 4 / (nu*dt + 2);
    mat->ade_a2 = dt*dt*EPSILON_0*omegap*omegap / (2 * (nu*dt + 2));
}

static void
calculate_ade_bp(gint i, SvMatProp *mat, gdouble a, gdouble phi, gdouble omega, gdouble gamma, gdouble dt)
{
    mat->ade_bp0[i] = -(omega*omega*dt*dt + (2 - gamma*dt)*(2 - gamma*dt)) / ((omega*omega*dt*dt + (2 + gamma*dt)*(2 + gamma*dt)));
    mat->ade_bp1[i] = 2 * (4 - (gamma*gamma + omega*omega)*dt*dt) / ((omega*omega*dt*dt + (2 + gamma*dt)*(2 + gamma*dt)));
    mat->ade_bp2[i] = 2 * EPSILON_0*a*omega*dt*(omega*dt*cos(phi) + (2 - gamma*dt)*sin(phi)) / ((omega*omega*dt*dt + (2 + gamma*dt)*(2 + gamma*dt)));
    mat->ade_bp3[i] = 4 * EPSILON_0*a*omega*dt*(omega*dt*cos(phi) - gamma*dt*sin(phi)) / ((omega*omega*dt*dt + (2 + gamma*dt)*(2 + gamma*dt)));
    mat->ade_bp4[i] = 2 * EPSILON_0*a*omega*dt*(omega*dt*cos(phi) - (2 + gamma*dt)*sin(phi)) / ((omega*omega*dt*dt + (2 + gamma*dt)*(2 + gamma*dt)));
}

static void
calculate_ade_c(SvMatProp *mat, gdouble epsilon, gdouble sigma, gdouble dt)
{
    mat->ade_c0 = dt / (sigma*dt / 2 + mat->ade_a2 + mat->ade_bp4[0] + mat->ade_bp4[1] + EPSILON_0*epsilon);
    mat->ade_c1 = 1 / (sigma*dt / 2 + mat->ade_a2 + mat->ade_bp4[0] + mat->ade_bp4[1] + EPSILON_0*epsilon);
    mat->ade_c2 = -(mat->ade_a2 + mat->ade_bp2[0] + mat->ade_bp2[1]) / (sigma*dt / 2 + mat->ade_a2 + mat->ade_bp4[0] + mat->ade_bp4[1] + EPSILON_0*epsilon);
    mat->ade_c3 = -(sigma*dt / 2 + 2 * mat->ade_a2 + mat->ade_bp3[0] + mat->ade_bp3[1] - EPSILON_0*epsilon)
        / (sigma*dt / 2 + mat->ade_a2 + mat->ade_bp4[0] + mat->ade_bp4[1] + EPSILON_0*epsilon);
}

static void
calculate_plrc_d(SvMatProp *mat, gdouble omegap, gdouble nu, gdouble dt)
{
    mat->plrc_d_chi = (omegap*omegap) / (nu*nu)*(exp(-nu*dt) + nu*dt - 1);
    mat->plrc_d_xi = mat->plrc_d_chi*(1 / (1 - exp(nu*dt)) + 1 / (nu*dt)) - (omegap*omegap) / (nu*nu)*(1 - nu*dt / 2 * (1 / tanh(nu*dt / 2)));
    mat->plrc_d_dchi = -(omegap*omegap) / (nu*nu)*(1 - exp(-nu*dt))*(1 - exp(-nu*dt));
    mat->plrc_d_dxi = mat->plrc_d_dchi*(1 / (1 - exp(nu*dt)) + 1 / (nu*dt));
}

static void
cdiv(gdouble a, gdouble b, gdouble c, gdouble d, gdouble *ra, gdouble *rb)
// (a+bi) / (c+di) = ((ac+bd) + (bc-ad)i) / (c^2+d^2)
{
    gdouble sub = c*c + d*d;
    *ra = (a*c + b*d) / sub;
    *rb = (b*c - a*d) / sub;
}

static void
cmult(gdouble a, gdouble b, gdouble c, gdouble d, gdouble *ra, gdouble *rb)
// (a+bi) * (c+di) = (ac-bd) + (bc+ad)i
{
    *ra = a*c - b*d;
    *rb = a*d + b*c;
}

static void
calculate_plrc_p(gint i, SvMatProp *mat, gdouble dt)
{
    //get_p_values(mat.cp3_a[0], mat.cp3_phi[0], mat.cp3_omega[0], mat.cp3_gamma[0], data->set->plan.dt,
//		     &(mat.plrc_p_chi[0]), &(mat.plrc_p_xi[0]), &(mat.plrc_p_dchir[0]), &(mat.plrc_p_dchii[0]), &(mat.plrc_p_dxir[0]), &(mat.plrc_p_dxii[0]));


    gdouble ra, ia, rb, ib, ic;

    ra = -2 * mat->cp3_a[i] * mat->cp3_omega[i] * sin(mat->cp3_phi[i]);
    ia = 2 * mat->cp3_a[i] * mat->cp3_omega[i] * cos(mat->cp3_phi[i]);

    rb = 1 - exp(-dt * mat->cp3_gamma[i]) * cos(dt * mat->cp3_omega[i]);
    ib = exp(-dt * mat->cp3_gamma[i]) * sin(dt * mat->cp3_omega[i]);

    cmult(ra, ia, rb, ib, &ra, &ia);

    /*eq. 12a*/
    cdiv(ra, ia, mat->cp3_gamma[i], mat->cp3_omega[i], &mat->plrc_p_chi[i], &ic);

    /*eq. 12a->14a*/
    cmult(mat->plrc_p_chi[i], ic, rb, ib, &mat->plrc_p_dchir[i], &mat->plrc_p_dchii[i]);

    /*eq. 12b*/
    cdiv(mat->plrc_p_chi[i], ic, 1 - exp(dt * mat->cp3_gamma[i]) * cos(dt * mat->cp3_omega[i]), -exp(dt * mat->cp3_gamma[i]) * sin(dt * mat->cp3_omega[i]), &ra, &ia);
    cdiv(mat->plrc_p_chi[i], ic, dt * mat->cp3_gamma[i], dt*mat->cp3_omega[i], &rb, &ib);
    mat->plrc_p_xi[i] = ra + rb;

    /*eq. 14a->14b*/
    cdiv(mat->plrc_p_dchir[i], mat->plrc_p_dchii[i],
         1 - exp(dt * mat->cp3_gamma[i]) * cos(dt * mat->cp3_omega[i]), -exp(dt * mat->cp3_gamma[i]) * sin(dt * mat->cp3_omega[i]), &ra, &ia);
    cdiv(mat->plrc_p_dchir[i], mat->plrc_p_dchii[i], dt * mat->cp3_gamma[i], dt * mat->cp3_omega[i], &rb, &ib);
    mat->plrc_p_dxir[i] = ra + rb;
    mat->plrc_p_dxii[i] = ia + ib;
}

#if(0)
static void
get_p_values(gdouble a, gdouble phi, gdouble omega, gdouble gamma, gdouble dt, gdouble *chi0, gdouble *xi0, gdouble *dchir, gdouble *dchii, gdouble *dxir, gdouble *dxii)
{
    gdouble ra, rb, rc, rd, ia, ib, ic, id;

    ra = -2 * a*omega*sin(phi);
    ia = 2 * a*omega*cos(phi);
    rb = 1 - exp(-dt*gamma)*cos(dt*omega);
    ib = exp(-dt*gamma)*sin(dt*omega);
    cmult(ra, ia, rb, ib, &ra, &ia);

    /*eq. 12a*/
    cdiv(ra, ia, gamma, omega, &rc, &ic);
    *chi0 = rc;

    /*eq. 12a->14a*/
    cmult(rc, ic, rb, ib, &rd, &id);
    *dchir = rd;
    *dchii = id;

    /*eq. 12b*/
    cdiv(rc, ic, 1 - exp(dt*gamma)*cos(dt*omega), -exp(dt*gamma)*sin(dt*omega), &rd, &id);
    *xi0 = rd;
    cdiv(rc, ic, dt*gamma, dt*omega, &rd, &id);
    *xi0 = (*xi0) + rd;

    rc = *dchir;
    ic = *dchii;

    /*eq. 14a->14b*/
    cdiv(rc, ic, 1 - exp(dt*gamma)*cos(dt*omega), -exp(dt*gamma)*sin(dt*omega), &rd, &id);
    *dxir = rd;
    *dxii = id;
    cdiv(rc, ic, dt*gamma, dt*omega, &rd, &id);
    *dxir = (*dxir) + rd;
    *dxii = (*dxii) + id;
}
#endif

static void
register_and_set_material(SvPool *mp, SvSet *set, SvMatProp mat, gint i, gint j, gint k)
{
    gint pos;
    pos = sv_pool_get_matindex(mp, mat);
    if (pos == -1)
        sv_pool_add_material(mp, &mat, FALSE);
    else
        mat.pos = pos;
    sv_pool_set_data(mp, mat, set, i, j, k);
}

static gboolean
resolve_mixed(SmCoords *coords, Submesh *sm, MpSetData *data)
{
    gint i, l, m, n;
    SvMatProp mat, smmat;
    gint cnt_drude, cnt_cp3, cnt_cp, cnt_ade, cnt_plrc, cnt_pec, cnt_max, mat_max;
    gdouble sum_epsilon, sum_sigma, sum_mu, sum_sigast, sum2_drude_omega_p, sum2_drude_nu, sum_cp3_a[3], sum_cp3_phi[3], sum2_cp3_omega[3], sum2_cp3_gamma[3];

    cnt_drude = cnt_cp3 = cnt_cp = cnt_ade = cnt_plrc = cnt_pec = 0;
    sum_epsilon = sum_sigma = sum_mu = sum_sigast = sum2_drude_omega_p = sum2_drude_nu = 0.0;
    sum_cp3_a[0] = sum_cp3_a[1] = sum_cp3_a[2] = 0.0;
    sum_cp3_phi[0] = sum_cp3_phi[1] = sum_cp3_phi[2] = 0.0;
    sum2_cp3_omega[0] = sum2_cp3_omega[1] = sum2_cp3_omega[2] = 0.0;
    sum2_cp3_gamma[0] = sum2_cp3_gamma[1] = sum2_cp3_gamma[2] = 0.0;

    mat = mat_vacuum_linear; // prefill with vacuum parameters

    if (data->mp->nmat > 0) // use tabulated materials if possible
        mat.type = SV_MAT_LINTAB;

    for (l = 0; l < 3; l++) {
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                smmat = (sm->material)[l][m][n];

                switch (smmat.type) {
                    case SV_MAT_DRUDE:
                        cnt_drude++; break;

                    case SV_MAT_CP3:
                        cnt_cp3++; break;

                    case SV_MAT_CP:
                        cnt_cp++; break;

                    case SV_MAT_ADE:
                        cnt_ade++; break;

                    case SV_MAT_PLRC:
                        cnt_plrc++; break;

                    case SV_MAT_PEC:
                        cnt_pec++; break;
                }

                sum_epsilon += smmat.epsilon;
                sum_sigma += smmat.sigma;
                sum_mu += smmat.mu;
                sum_sigast += smmat.sigast;
                sum2_drude_omega_p += smmat.drude_omega_p * smmat.drude_omega_p;
                sum2_drude_nu += smmat.drude_nu * smmat.drude_nu;
                for (i = 0; i < 3; i++) {
                    sum_cp3_a[i] += smmat.cp3_a[i];
                    sum_cp3_phi[i] += smmat.cp3_phi[i];
                    sum2_cp3_omega[i] += smmat.cp3_omega[i] * smmat.cp3_omega[i];
                    sum2_cp3_gamma[i] += smmat.cp3_gamma[i] * smmat.cp3_gamma[i];
                }
            } // n
        } // m
    } // l

    // if there are any metals present, pick the most frequent one
    mat_max = mat.type;
    cnt_max = 0;
    if (cnt_drude > cnt_max) {
        mat_max = SV_MAT_DRUDE;
        cnt_max = cnt_drude;
    }
    if (cnt_cp3 > cnt_max) {
        mat_max = SV_MAT_CP3;
        cnt_max = cnt_cp3;
    }
    if (cnt_cp > cnt_max) {
        mat_max = SV_MAT_CP;
        cnt_max = cnt_cp;
    }
    if (cnt_ade > cnt_max) {
        mat_max = SV_MAT_ADE;
        cnt_max = cnt_ade;
    }
    if (cnt_plrc > cnt_max) {
        mat_max = SV_MAT_PLRC;
        cnt_max = cnt_plrc;
    }
    if (cnt_pec > cnt_max) {
        mat_max = SV_MAT_PEC;
        cnt_max = cnt_pec;
    }
    mat.type = mat_max;

    // suppose that these materials only appear next to vacuum; for arbitrary mixtures, implement a more complex comparison;
    // at present, the winner is given by the order below, not by the actual count of subcells
    /*
    if (cnt_cp >= 1 && cnt_cp > cnt_drude && cnt_cp > cnt_cp3 && cnt_cp > cnt_ade && cnt_cp > cnt_plrc && cnt_cp > cnt_pec))
    mat.type = SV_MAT_CP;
    else if (cnt_ade >= 1)
    mat.type = SV_MAT_ADE;
    else if (cnt_plrc >= 1)
    mat.type = SV_MAT_PLRC;
    else if (cnt_pec >= 1)
    mat.type = SV_MAT_PEC;
    // else: mat.type remains LINEAR or LINTAB
    */

    switch (mat.type) {
        case SV_MAT_LINEAR:
            mat.epsilon = sum_epsilon / 27;
            mat.sigma = sum_sigma / 27;
            mat.mu = sum_mu / 27;
            mat.sigast = sum_sigast / 27;

            sv_pool_set_data(data->mp, mat, data->set, coords->i, coords->j, coords->k);
            break;

        case SV_MAT_LINTAB:
            mat.epsilon = sum_epsilon / 27;
            mat.sigma = sum_sigma / 27;
            mat.mu = sum_mu / 27;
            mat.sigast = sum_sigast / 27;

            register_and_set_material(data->mp, data->set, mat, coords->i, coords->j, coords->k);
            break;

        case SV_MAT_DRUDE:
            mat.epsilon = sum_epsilon / 27;
            mat.drude_omega_p = sqrt(sum2_drude_omega_p / 27);
            mat.drude_nu = sqrt(sum2_drude_nu / 27);

            register_and_set_material(data->mp, data->set, mat, coords->i, coords->j, coords->k);
            break;

        case SV_MAT_CP3:
            mat.epsilon = sum_epsilon / 27;
            mat.sigma = sum_sigma / 27;
            for (i = 0; i < 3; i++) {
                mat.cp3_a[i] = sum_cp3_a[i] / 27;
                mat.cp3_phi[i] = sum_cp3_phi[i] / 27;
                mat.cp3_omega[i] = sqrt(sum2_cp3_omega[i] / 27);
                mat.cp3_gamma[i] = sqrt(sum2_cp3_gamma[i] / 27);
            }

            register_and_set_material(data->mp, data->set, mat, coords->i, coords->j, coords->k);
            break;

        case SV_MAT_CP:
            mat.epsilon = sum_epsilon / 27;
            mat.drude_omega_p = sqrt(sum2_drude_omega_p / 27);
            mat.drude_nu = sqrt(sum2_drude_nu / 27);
            for (i = 0; i < 2; i++) {
                mat.cp3_a[i] = sum_cp3_a[i] / 27;
                mat.cp3_phi[i] = sum_cp3_phi[i] / 27;
                mat.cp3_omega[i] = sqrt(sum2_cp3_omega[i] / 27);
                mat.cp3_gamma[i] = sqrt(sum2_cp3_gamma[i] / 27);
            }

            register_and_set_material(data->mp, data->set, mat, coords->i, coords->j, coords->k);
            break;

        case SV_MAT_ADE:
            mat.epsilon = sum_epsilon / 27;
            mat.sigma = sum_sigma / 27;
            mat.drude_omega_p = sqrt(sum2_drude_omega_p / 27);
            mat.drude_nu = sqrt(sum2_drude_nu / 27);
            for (i = 0; i < 2; i++) {
                mat.cp3_a[i] = sum_cp3_a[i] / 27;
                mat.cp3_phi[i] = sum_cp3_phi[i] / 27;
                mat.cp3_omega[i] = sqrt(sum2_cp3_omega[i] / 27);
                mat.cp3_gamma[i] = sqrt(sum2_cp3_gamma[i] / 27);
            }

            calculate_ade_a(&mat, mat.drude_omega_p, mat.drude_nu, data->set->plan.dt);
            calculate_ade_bp(0, &mat, mat.cp3_a[0], mat.cp3_phi[0], mat.cp3_omega[0], mat.cp3_gamma[0], data->set->plan.dt);
            calculate_ade_bp(1, &mat, mat.cp3_a[1], mat.cp3_phi[1], mat.cp3_omega[1], mat.cp3_gamma[1], data->set->plan.dt);
            calculate_ade_c(&mat, mat.epsilon, mat.sigma, data->set->plan.dt);

            register_and_set_material(data->mp, data->set, mat, coords->i, coords->j, coords->k);
            break;

        case SV_MAT_PLRC:
            mat.epsilon = sum_epsilon / 27;
            mat.sigma = sum_sigma / 27;
            mat.drude_omega_p = sqrt(sum2_drude_omega_p / 27);
            mat.drude_nu = sqrt(sum2_drude_nu / 27);
            for (i = 0; i < 2; i++) {
                mat.cp3_a[i] = sum_cp3_a[i] / 27;
                mat.cp3_phi[i] = sum_cp3_phi[i] / 27;
                mat.cp3_omega[i] = sqrt(sum2_cp3_omega[i] / 27);
                mat.cp3_gamma[i] = sqrt(sum2_cp3_gamma[i] / 27);
            }

            calculate_plrc_d(&mat, mat.drude_omega_p, mat.drude_nu, data->set->plan.dt);
            calculate_plrc_p(0, &mat, data->set->plan.dt);
            calculate_plrc_p(1, &mat, data->set->plan.dt);

            /*
            get_p_values(mat.cp3_a[0], mat.cp3_phi[0], mat.cp3_omega[0], mat.cp3_gamma[0], data->set->plan.dt,
                     &(mat.plrc_p_chi[0]), &(mat.plrc_p_xi[0]), &(mat.plrc_p_dchir[0]), &(mat.plrc_p_dchii[0]), &(mat.plrc_p_dxir[0]), &(mat.plrc_p_dxii[0]));
            get_p_values(mat.cp3_a[1], mat.cp3_phi[1], mat.cp3_omega[1], mat.cp3_gamma[1], data->set->plan.dt,
                     &(mat.plrc_p_chi[1]), &(mat.plrc_p_xi[1]), &(mat.plrc_p_dchir[1]), &(mat.plrc_p_dchii[1]), &(mat.plrc_p_dxir[1]), &(mat.plrc_p_dxii[1]));
            */

            register_and_set_material(data->mp, data->set, mat, coords->i, coords->j, coords->k);
            break;

        case SV_MAT_PEC: // impossible to mix
            mat.pos = sv_pool_get_matindex(data->mp, mat);
            sv_pool_set_data(data->mp, data->mp->mats[mat.pos], data->set, coords->i, coords->j, coords->k);
            break;

        default: // DATABASE -- TODO; now simply pick material from the center
            sv_pool_set_data(data->mp, (sm->material)[1][1][1], data->set, coords->i, coords->j, coords->k);
    }

    g_free(sm);
    return 0;

    /*
    // eps_eff = (mat.epsilon + (1-p)*1) * EPSILON_0; // primitive convex combination
    //  mp->d->sigma->data[i][j][k] = p * mat.sigma; // simple fraction
    // Looyenga power-law mixing formula (exponent 1/3): eps_eff^(1/3) = p*eps_1^(1/3) + (1-p)*eps_2^(1/3)
    // - empiricaly proved to yield usable results
    // Is the power formula for conductivities better than simple fraction?
    mp->d->epsilon->data[i][j][k] = exp( 3 *log(p*exp(1.0/3.0 * log(mat.epsilon)) + (1-p))) * EPSILON_0;
    mp->d->sigma->data[i][j][k] = exp(3 * log(p)) * mat.sigma;
    // mu_eff = (mu + (1-p)*1) * MU_0; // primitive convex combination
    // mp->d->sigast->data[i][j][k] = p * mat.sigast; // simple fraction
    mp->d->mu->data[i][j][k] = exp(3 * log(p*exp(1.0/3.0 * log(mat.mu)) + (1-p))) * MU_0;
    mp->d->sigast->data[i][j][k] = exp(3 * log(p)) * mat.sigast;
    if (mp->d->mat) mp->d->mat->data[i][j][k] = 0;
    }
    else mp->d->mat->data[i][j][k] = mat.pos; // !!! poresit doplnkove hodnoty v tabulce
    */
} // pruchod submeshemi

gboolean
scan_material(FILE *fr, SvMatProp *mat, gdouble dt)
{
    float epsilon, mu, sigma, sigast, omegap, nu, a, omega, phi, gamma;
    int type;

    // initialize all variables that might not be used
    mat->epsilon = mat->mu = 1;
    mat->sigma = mat->sigast = 0.0;
    mat->drude_omega_p = mat->drude_nu = 0.0;
    mat->cp3_a[0] = mat->cp3_a[1] = mat->cp3_a[2] = 0.0;
    mat->cp3_phi[0] = mat->cp3_phi[1] = mat->cp3_phi[2] = 0.0;
    mat->cp3_omega[0] = mat->cp3_omega[1] = mat->cp3_omega[2] = 0.0;
    mat->cp3_gamma[0] = mat->cp3_gamma[1] = mat->cp3_gamma[2] = 0.0;

    if (fscanf(fr, "%d", &type) != 1)
        return 1;
    if (type == SV_MAT_LINEAR || type == SV_MAT_LINTAB) {
        mat->type = type;

        if (fscanf(fr, "%f", &epsilon) != 1)
            return 1; //TODO this should be locale dependent
        if (fscanf(fr, "%f", &mu) != 1)
            return 1;
        if (fscanf(fr, "%f", &sigma) != 1)
            return 1;
        if (fscanf(fr, "%f", &sigast) != 1)
            return 1;
        mat->epsilon = epsilon;
        mat->mu = mu;
        mat->sigma = sigma;
        mat->sigast = sigast;
    } else if (type == SV_MAT_PEC) {
        mat->type = type;
    } else if (type == SV_MAT_DRUDE) {
        mat->type = type;

        if (fscanf(fr, "%g", &epsilon) != 1)
            return 1; //TODO this should be locale dependent
        if (fscanf(fr, "%g", &omegap) != 1)
            return 1;
        if (fscanf(fr, "%g", &nu) != 1)
            return 1;
        mat->epsilon = epsilon;
        mat->drude_omega_p = omegap;
        mat->drude_nu = nu;
    } else if (type == SV_MAT_CP) {
        mat->type = type;

        if (fscanf(fr, "%g", &epsilon) != 1)
            return 1;
        if (fscanf(fr, "%g", &omegap) != 1)
            return 1;
        if (fscanf(fr, "%g", &nu) != 1)
            return 1;
        mat->epsilon = epsilon;
        mat->drude_omega_p = omegap;
        mat->drude_nu = nu;

        if (fscanf(fr, "%g", &a) != 1)
            return 1;
        if (fscanf(fr, "%g", &phi) != 1)
            return 1;
        if (fscanf(fr, "%g", &omega) != 1)
            return 1;
        if (fscanf(fr, "%g", &gamma) != 1)
            return 1;
        mat->cp3_a[0] = a;
        mat->cp3_phi[0] = phi;
        mat->cp3_omega[0] = omega;
        mat->cp3_gamma[0] = gamma;
        //printf("%g %g %g %g\n", a, phi, omega, gamma);

        if (fscanf(fr, "%g", &a) != 1)
            return 1;
        if (fscanf(fr, "%g", &phi) != 1)
            return 1;
        if (fscanf(fr, "%g", &omega) != 1)
            return 1;
        if (fscanf(fr, "%g", &gamma) != 1)
            return 1;
        mat->cp3_a[1] = a;
        mat->cp3_phi[1] = phi;
        mat->cp3_omega[1] = omega;
        mat->cp3_gamma[1] = gamma;

        //printf("loaded a/gamma[0]: %g %g\n", mat->cp3_a[0], mat->cp3_gamma[0]);
    } else if (type == SV_MAT_ADE) {
        mat->type = type;

        if (fscanf(fr, "%g", &epsilon) != 1)
            return 1;
        if (fscanf(fr, "%g", &sigma) != 1)
            return 1;
        if (fscanf(fr, "%g", &omegap) != 1)
            return 1;
        if (fscanf(fr, "%g", &nu) != 1)
            return 1;
        mat->epsilon = epsilon;
        mat->sigma = sigma;
        mat->drude_omega_p = omegap;
        mat->drude_nu = nu;

        calculate_ade_a(mat, omegap, nu, dt);

        /*
            mat->ade_a0 = (nu*dt-2) / (nu*dt+2);
            mat->ade_a1 = 4 / (nu*dt+2);
            mat->ade_a2 = dt*dt*EPSILON_0*omegap*omegap / (2*(nu*dt+2));
        */

        if (fscanf(fr, "%g", &a) != 1)
            return 1;
        if (fscanf(fr, "%g", &phi) != 1)
            return 1;
        if (fscanf(fr, "%g", &omega) != 1)
            return 1;
        if (fscanf(fr, "%g", &gamma) != 1)
            return 1;
        mat->cp3_a[0] = a;
        mat->cp3_phi[0] = phi;
        mat->cp3_omega[0] = omega;
        mat->cp3_gamma[0] = gamma;

        calculate_ade_bp(0, mat, a, phi, omega, gamma, dt);

        /*
            mat->ade_bp0[0] = -(omega*omega*dt*dt + (2-gamma*dt)*(2-gamma*dt)) / ((omega*omega*dt*dt + (2+gamma*dt)*(2+gamma*dt)));
            mat->ade_bp1[0] = 2*(4 - (gamma*gamma + omega*omega)*dt*dt) / ((omega*omega*dt*dt + (2+gamma*dt)*(2+gamma*dt)));
            mat->ade_bp2[0] = 2*EPSILON_0*a*omega*dt*(omega*dt*cos(phi) + (2-gamma*dt)*sin(phi)) / ((omega*omega*dt*dt + (2+gamma*dt)*(2+gamma*dt)));
            mat->ade_bp3[0] = 4*EPSILON_0*a*omega*dt*(omega*dt*cos(phi) - gamma*dt*sin(phi)) / ((omega*omega*dt*dt + (2+gamma*dt)*(2+gamma*dt)));
            mat->ade_bp4[0] = 2*EPSILON_0*a*omega*dt*(omega*dt*cos(phi) - (2+gamma*dt)*sin(phi)) / ((omega*omega*dt*dt + (2+gamma*dt)*(2+gamma*dt)));
        */

        //printf("%g %g %g %g\n", a, phi, omega, gamma);

        if (fscanf(fr, "%g", &a) != 1)
            return 1;
        if (fscanf(fr, "%g", &phi) != 1)
            return 1;
        if (fscanf(fr, "%g", &omega) != 1)
            return 1;
        if (fscanf(fr, "%g", &gamma) != 1)
            return 1;
        mat->cp3_a[1] = a;
        mat->cp3_phi[1] = phi;
        mat->cp3_omega[1] = omega;
        mat->cp3_gamma[1] = gamma;

        calculate_ade_bp(1, mat, a, phi, omega, gamma, dt);

        /*
            mat->ade_bp0[1] = -(omega*omega*dt*dt + (2-gamma*dt)*(2-gamma*dt)) / ((omega*omega*dt*dt + (2+gamma*dt)*(2+gamma*dt)));
            mat->ade_bp1[1] = 2*(4 - (gamma*gamma + omega*omega)*dt*dt) / ((omega*omega*dt*dt + (2+gamma*dt)*(2+gamma*dt)));
            mat->ade_bp2[1] = 2*EPSILON_0*a*omega*dt*(omega*dt*cos(phi) + (2-gamma*dt)*sin(phi)) / ((omega*omega*dt*dt + (2+gamma*dt)*(2+gamma*dt)));
            mat->ade_bp3[1] = 4*EPSILON_0*a*omega*dt*(omega*dt*cos(phi) - gamma*dt*sin(phi)) / ((omega*omega*dt*dt + (2+gamma*dt)*(2+gamma*dt)));
            mat->ade_bp4[1] = 2*EPSILON_0*a*omega*dt*(omega*dt*cos(phi) - (2+gamma*dt)*sin(phi)) / ((omega*omega*dt*dt + (2+gamma*dt)*(2+gamma*dt)));
        */

        calculate_ade_c(mat, epsilon, sigma, dt);

        /*
            mat->ade_c0 = dt / (sigma*dt/2 + mat->ade_a2 + mat->ade_bp4[0] + mat->ade_bp4[1] + EPSILON_0*epsilon);
            mat->ade_c1 = 1 / (sigma*dt/2 + mat->ade_a2 + mat->ade_bp4[0] + mat->ade_bp4[1] + EPSILON_0*epsilon);
            mat->ade_c2 = -(mat->ade_a2 + mat->ade_bp2[0] + mat->ade_bp2[1]) / (sigma*dt/2 + mat->ade_a2 + mat->ade_bp4[0] + mat->ade_bp4[1] + EPSILON_0*epsilon);
            mat->ade_c3 = -(sigma*dt/2 + 2*mat->ade_a2 + mat->ade_bp3[0] + mat->ade_bp3[1] - EPSILON_0*epsilon)
            / (sigma*dt/2 + mat->ade_a2 + mat->ade_bp4[0] + mat->ade_bp4[1] + EPSILON_0*epsilon);
        */

        /* printf("ADE material parameters: a: %g %g %g    b: %g %g %g %g %g  %g %g %g %g %g,    c: %g %g %g %g\n", mat->ade_a0, mat->ade_a1, mat->ade_a2, mat->ade_bp0[0], mat->ade_bp1[0], mat->ade_bp2[0], mat->ade_bp3[0], mat->ade_bp4[0], mat->ade_bp0[1], mat->ade_bp1[1], mat->ade_bp2[1], mat->ade_bp3[1], mat->ade_bp4[1], mat->ade_c0, mat->ade_c1, mat->ade_c2, mat->ade_c3); */

        // printf("loaded a/gamma[0]: %g %g\n", mat->cp3_a[0], mat->cp3_gamma[0]);
    } else if (type == SV_MAT_PLRC) {
        mat->type = type;

        if (fscanf(fr, "%g", &epsilon) != 1)
            return 1;
        if (fscanf(fr, "%g", &sigma) != 1)
            return 1;
        if (fscanf(fr, "%g", &omegap) != 1)
            return 1;
        if (fscanf(fr, "%g", &nu) != 1)
            return 1;
        mat->epsilon = epsilon;
        mat->sigma = sigma;
        mat->drude_omega_p = omegap;
        mat->drude_nu = nu;

        calculate_plrc_d(mat, omegap, nu, dt);

        /*
            mat->plrc_d_chi = (omegap*omegap)/(nu*nu)*(exp(-nu*dt) + nu*dt - 1);
            mat->plrc_d_xi = mat->plrc_d_chi*(1/(1-exp(nu*dt)) + 1/(nu*dt)) - (omegap*omegap)/(nu*nu)*(1 - nu*dt/2*(1/tanh(nu*dt/2)));
            mat->plrc_d_dchi = -(omegap*omegap)/(nu*nu)*(1 - exp(-nu*dt))*(1 - exp(-nu*dt));
            mat->plrc_d_dxi = mat->plrc_d_dchi*(1/(1-exp(nu*dt)) + 1/(nu*dt));
        */

        if (fscanf(fr, "%g", &a) != 1)
            return 1;
        if (fscanf(fr, "%g", &phi) != 1)
            return 1;
        if (fscanf(fr, "%g", &omega) != 1)
            return 1;
        if (fscanf(fr, "%g", &gamma) != 1)
            return 1;
        mat->cp3_a[0] = a;
        mat->cp3_phi[0] = phi;
        mat->cp3_omega[0] = omega;
        mat->cp3_gamma[0] = gamma;

        calculate_plrc_p(0, mat, dt);

        // get_p_values replaced by calculate_plrc_p
            // get_p_values(a, phi, omega, gamma, dt, &(mat->plrc_p_chi[0]), &(mat->plrc_p_xi[0]), &(mat->plrc_p_dchir[0]), &(mat->plrc_p_dchii[0]), &(mat->plrc_p_dxir[0]), &(mat->plrc_p_dxii[0]));

            /* printf("plrc dterm %g %g %g %g p (%g %g) (%g %g) (%g %g), (%g %g)\n", mat->plrc_d_chi, mat->plrc_d_xi, mat->plrc_d_dchi, mat->plrc_d_dxi, mat->plrc_p_chi[0], 0.0, mat->plrc_p_xi[0], 0.0, mat->plrc_p_dchir[0], mat->plrc_p_dchii[0], mat->plrc_p_dxir[0], mat->plrc_p_dxii[0]); */

        if (fscanf(fr, "%g", &a) != 1)
            return 1;
        if (fscanf(fr, "%g", &phi) != 1)
            return 1;
        if (fscanf(fr, "%g", &omega) != 1)
            return 1;
        if (fscanf(fr, "%g", &gamma) != 1)
            return 1;
        mat->cp3_a[1] = a;
        mat->cp3_phi[1] = phi;
        mat->cp3_omega[1] = omega;
        mat->cp3_gamma[1] = gamma;

        calculate_plrc_p(1, mat, dt);
        // get_p_values(a, phi, omega, gamma, dt, &(mat->plrc_p_chi[1]), &(mat->plrc_p_xi[1]), &(mat->plrc_p_dchir[1]), &(mat->plrc_p_dchii[1]), &(mat->plrc_p_dxir[1]), &(mat->plrc_p_dxii[1]));

        /* printf("plrc dterm %g %g %g %g p (%g %g) (%g %g) (%g %g), (%g %g)\n", mat->plrc_d_chi, mat->plrc_d_xi, mat->plrc_d_dchi, mat->plrc_d_dxi, mat->plrc_p_chi[1], 0.0, mat->plrc_p_xi[1], 0.0, mat->plrc_p_dchir[1], mat->plrc_p_dchii[1], mat->plrc_p_dxir[1], mat->plrc_p_dxii[1]); */

    } else if (type == SV_MAT_CP3) {
        mat->type = type;

        if (fscanf(fr, "%g", &epsilon) != 1)
            return 1;
        if (fscanf(fr, "%g", &sigma) != 1)
            return 1;
        if (fscanf(fr, "%g", &a) != 1)
            return 1;
        if (fscanf(fr, "%g", &phi) != 1)
            return 1;
        if (fscanf(fr, "%g", &omega) != 1)
            return 1;
        if (fscanf(fr, "%g", &gamma) != 1)
            return 1;
        mat->epsilon = epsilon;
        mat->sigma = sigma;
        mat->cp3_a[0] = a;
        mat->cp3_phi[0] = phi;
        mat->cp3_omega[0] = omega;
        mat->cp3_gamma[0] = gamma;

        if (fscanf(fr, "%g", &a) != 1)
            return 1;
        if (fscanf(fr, "%g", &phi) != 1)
            return 1;
        if (fscanf(fr, "%g", &omega) != 1)
            return 1;
        if (fscanf(fr, "%g", &gamma) != 1)
            return 1;
        mat->cp3_a[1] = a;
        mat->cp3_phi[1] = phi;
        mat->cp3_omega[1] = omega;
        mat->cp3_gamma[1] = gamma;

        if (fscanf(fr, "%g", &a) != 1)
            return 1;
        if (fscanf(fr, "%g", &phi) != 1)
            return 1;
        if (fscanf(fr, "%g", &omega) != 1)
            return 1;
        if (fscanf(fr, "%g", &gamma) != 1)
            return 1;
        mat->cp3_a[2] = a;
        mat->cp3_phi[2] = phi;
        mat->cp3_omega[2] = omega;
        mat->cp3_gamma[2] = gamma;

        printf("loaded a/gamma[0]: %g %g\n", mat->cp3_a[0], mat->cp3_gamma[0]);
    }

    /*
    if (mat->type != SV_MAT_LINEAR) {
    if (sv_pool_get_matindex(mp, *mat) == -1)
        sv_pool_add_material(mp, mat, FALSE);
    }
    */

    //printf("scan material finished\n");
    return 0;
}

gint
load_tets(gchar *filebase, GArray *tetrahedrons, gint attribute_pos, gint attribute_val,
          gdouble xshift, gdouble yshift, gdouble zshift, gdouble xmult, gdouble ymult, gdouble zmult,
          gint verbose, SvMatProp *mat, gint *ngtot)
{
    gchar filenodes[256];
    gchar filetets[256];
    FILE *fr;
    gint i, j, nn, natt, nt, npnts, buf, nloaded, nodecorr;
    gint tp[4];
    gfloat *nx, *ny, *nz;
    gchar c, line[100];
    SvTetrahedron tx;
    gint attval[100];
    gdouble xmin, xmax, ymin, ymax, zmin, zmax;

    g_snprintf(filenodes, sizeof(filenodes), "%s.node", filebase);
    g_snprintf(filetets, sizeof(filetets), "%s.ele", filebase);

    if (verbose)
        printf("Tetrahedrons will be loaded from %s and %s, at position %d, value %d\n", filenodes, filetets, attribute_pos, attribute_val);

    if (!(fr = fopen(filenodes, "r"))) {
        fprintf(stderr, "Node file %s cannot be opened for reading\n", filenodes);
        return 1;
    }

    do {
        c = getc(fr);
        ungetc(c, fr);
        if (c == '#' || c == '\n') {
            fgets(line, 100, fr);
        } else
            break;
    } while (1);

    if (fscanf(fr, "%d", &nn) != 1) {
        fclose(fr);
        return 1;
    }  //number of nodes
    if (fscanf(fr, "%d", &buf) != 1) {
        fclose(fr);
        return 1;
    }   //always 3
    if (fscanf(fr, "%d", &natt) != 1) {
        fclose(fr);
        return 1;
    }  //number of attributes
    if (fscanf(fr, "%d", &buf) != 1) {
        fclose(fr);
        return 1;
    }    //boundary markers - do not care

    nx = (gfloat *)g_malloc(nn * sizeof(gfloat));
    ny = (gfloat *)g_malloc(nn * sizeof(gfloat));
    nz = (gfloat *)g_malloc(nn * sizeof(gfloat));

    do {
        c = getc(fr);
        ungetc(c, fr);
        if (c == '#' || c == '\n') {
            fgets(line, 100, fr);
        } else
            break;
    } while (1);

    xmin = ymin = zmin = G_MAXDOUBLE;
    xmax = ymax = zmax = -G_MAXDOUBLE;

    nodecorr = 0; //sometimes nodes are indexed from 1

    for (i = 0; i < nn; i++) {
        fscanf(fr, "%d", &buf);
        if (i == 0 && buf == 1)
            nodecorr = 1;

        fscanf(fr, "%f", nx + i);
        fscanf(fr, "%f", ny + i);
        fscanf(fr, "%f", nz + i);
        nx[i] = (gfloat)(nx[i] * xmult + xshift);
        ny[i] = (gfloat)(ny[i] * ymult + yshift);
        nz[i] = (gfloat)(nz[i] * zmult + zshift);
        //printf("Loading %d: %g %g %g\n", i, nx[i], ny[i], nz[i]);

        if (nx[i] > xmax)
            xmax = nx[i];
        if (nx[i] < xmin)
            xmin = nx[i];
        if (ny[i] > ymax)
            ymax = ny[i];
        if (ny[i] < ymin)
            ymin = ny[i];
        if (nz[i] > zmax)
            zmax = nz[i];
        if (nz[i] < zmin)
            zmin = nz[i];

        for (j = 0; j < natt; j++)
            fscanf(fr, "%d", &buf);
    }
    fclose(fr);

    if (verbose)
        printf("%d nodes loaded from %s file\n", nn, filenodes);

    if (!(fr = fopen(filetets, "r"))) {
        fprintf(stderr, "Element file %s cannot be opened for reading\n", filetets);
        return 1;
    }

    do {
        c = getc(fr);
        ungetc(c, fr);
        if (c == '#' || c == '\n') {
            fgets(line, 100, fr);
        } else
            break;
    } while (1);

    if (fscanf(fr, "%d", &nt) != 1) {
        fclose(fr);
        return 1;
    }    //number of tetrahedrons
    if (fscanf(fr, "%d", &npnts) != 1) {
        fclose(fr);
        return 1;
    }    //points in tetrahedron (4 or 10)
    if (fscanf(fr, "%d", &natt) != 1) {
        fclose(fr);
        return 1;
    }  //number of attributes

    do {
        c = getc(fr);
        ungetc(c, fr);
        if (c == '#' || c == '\n') {
            fgets(line, 100, fr);
        } else
            break;
    } while (1);

    nloaded = 0;

    if (verbose) {
        if (natt && attribute_pos != -1 && attribute_val != -1)
            printf("File %s has %d attributes, nodes with attribute %d value %d will be used only\n", filetets, natt, attribute_pos, attribute_val);
        else
            printf("File attributes will be ignored, all tetrahedrons loaded.\n");
    }
    for (i = 0; i < nt; i++) {

        if (fscanf(fr, "%d", &buf) != 1) {
            fclose(fr);
            return 1;
        }
        if (fscanf(fr, "%d", tp) != 1) {
            fclose(fr);
            return 1;
        }
        if (fscanf(fr, "%d", tp + 1) != 1) {
            fclose(fr);
            return 1;
        }
        if (fscanf(fr, "%d", tp + 2) != 1) {
            fclose(fr);
            return 1;
        }
        if (fscanf(fr, "%d", tp + 3) != 1) {
            fclose(fr);
            return 1;
        }
        if (npnts == 10)
            for (j = 0; j < 6; j++)
                if (fscanf(fr, "%d", &buf) != 1) {
                    fclose(fr);
                    return 1;
                }

        for (j = 0; j < natt; j++) {
            if (fscanf(fr, "%d", attval + j) != 1) {
                fclose(fr);
                return 1;
            }
        }

        //printf("%d vs %d, %d\n", attval[attribute_pos], attribute_val, attribute_pos);

        if (natt == 0 || (attribute_pos == -1) || (attribute_val == -1) || (natt > 0 && (natt <= attribute_pos) && (attval[attribute_pos] == attribute_val))) {

            nloaded++;

            //printf("tetrahedron %d: %d %d %d %d\n", i, tp[0], tp[1], tp[2], tp[3]);

            tx.pnt1[0] = nx[tp[0] - nodecorr]; tx.pnt1[1] = ny[tp[0] - nodecorr]; tx.pnt1[2] = nz[tp[0] - nodecorr];
            tx.pnt2[0] = nx[tp[1] - nodecorr]; tx.pnt2[1] = ny[tp[1] - nodecorr]; tx.pnt2[2] = nz[tp[1] - nodecorr];
            tx.pnt3[0] = nx[tp[2] - nodecorr]; tx.pnt3[1] = ny[tp[2] - nodecorr]; tx.pnt3[2] = nz[tp[2] - nodecorr];
            tx.pnt4[0] = nx[tp[3] - nodecorr]; tx.pnt4[1] = ny[tp[3] - nodecorr]; tx.pnt4[2] = nz[tp[3] - nodecorr];
            tx.setpart = 1;

            //printf("tetrahedron %d really: %g %g %g   %g %g %g   %g %g %g   %g %g %g\n", i, tx.pnt1[0], tx.pnt1[1], tx.pnt1[2], tx.pnt2[0], tx.pnt2[1], tx.pnt2[2], tx.pnt3[0], tx.pnt3[1], tx.pnt3[2], tx.pnt4[0], tx.pnt4[1], tx.pnt4[2]);

            tx.mat.pos = mat->pos;
            tx.mat.type = mat->type; tx.mat.epsilon = mat->epsilon, tx.mat.mu = mat->mu;
            tx.mat.sigma = mat->sigma; tx.mat.sigast = mat->sigast, tx.mat.drude_omega_p = mat->drude_omega_p;
            tx.mat.drude_nu = mat->drude_nu;
            for (j = 0; j < 3; j++) {
                tx.mat.cp3_a[j] = mat->cp3_a[j];
                tx.mat.cp3_phi[j] = mat->cp3_phi[j];
                tx.mat.cp3_omega[j] = mat->cp3_omega[j];
                tx.mat.cp3_gamma[j] = mat->cp3_gamma[j];
            }
            tx.n = (*ngtot)++;

            g_array_append_val(tetrahedrons, tx);
        }
    }
    fclose(fr);

    if (verbose)
        printf("%d tetrahedrons loaded from %s file, span (%d %d %d) ... (%d %d %d) voxels\n", nloaded, filetets, (gint)xmin, (gint)ymin, (gint)zmin, (gint)xmax, (gint)ymax, (gint)zmax);

    return 0;
}

#define MAGIC "GWYO"
#define MAGIC2 "GWYP"
#define MAGIC_SIZE (sizeof(MAGIC)-1)

gint
load_gwydd(gchar *filename, SvGwydd *gx, gint verbose, gint channel)
{
    gsize size = 0;
    gsize pos = 0;
    guchar *buffer = NULL;
    GError *err = NULL;
    GObject *object;
    GwyContainer *container;
    gchar data_key[30];
    gchar mask_key[30];
    const gchar **keys;
    guint i;

    if (verbose)
        printf("Loading Gwyddion file %s\n", filename);

    if (!gwy_file_get_contents(filename, &buffer, &size, &err)) {
        fprintf(stderr, "Error: cannot load file\n");
        gx->dfield = NULL;
        return 1;
    }
    if (size < MAGIC_SIZE || (memcmp(buffer, MAGIC, MAGIC_SIZE) && memcmp(buffer, MAGIC2, MAGIC_SIZE))) {
        fprintf(stderr, "Error: wrong file type.\n");
        gwy_file_abandon_contents(buffer, size, &err);
        gx->dfield = NULL;
        return 1;
    }
    if (!memcmp(buffer, MAGIC, MAGIC_SIZE)) {
        fprintf(stderr, "Error: File in old format, now unsupported.\n");
        gwy_file_abandon_contents(buffer, size, &err);
        gx->dfield = NULL;
        return 1;
    } else
        object = gwy_serializable_deserialize(buffer + MAGIC_SIZE,
                                              size - MAGIC_SIZE, &pos);

    gwy_file_abandon_contents(buffer, size, &err);
    if (!object) {
        fprintf(stderr, "Error: Deserialization failed.\n");
        gx->dfield = NULL;
        return 1;
    }
    if (!GWY_IS_CONTAINER(object)) {
        fprintf(stderr, "Error: Deserialization resulted in something unexpected.\n");
        g_object_unref(object);
        gx->dfield = NULL;
        return 1;
    }

    container = GWY_CONTAINER(object);
    gwy_app_data_browser_add(container);

    g_snprintf(data_key, sizeof(data_key), "/%i/data", channel);
    g_snprintf(mask_key, sizeof(mask_key), "/%i/mask", channel);

    gx->dfield = gwy_container_get_object(container, g_quark_from_string(data_key));
    if (!gx->dfield) {
        fprintf(stderr, "Error: datafield not found, here follows the list of keys in file\n");
        keys = gwy_container_keys_by_name(container);
        for (i = 0; i < gwy_container_get_n_items(container); i++) {
            printf("%s\n", keys[i]);
        }
        gx->dfield = NULL;
        return 1;
    }

    if (gx->mask != 0)
    {
        gx->mfield = gwy_container_get_object(container, g_quark_from_string(mask_key));

        if (!gx->mfield) {
            fprintf(stderr, "Error: mask not found but requested to be applied\n");
        }
    }

    if (!GWY_IS_DATA_FIELD(gx->dfield)) {
        fprintf(stderr, "Error: datafield not loaded properly\n");
        gx->dfield = NULL;
        return 1;
    }

    if (GWY_IS_DATA_FIELD(gx->dfield))
        printf("Datafield loaded\n");
    if (gx->mask != 0 && GWY_IS_DATA_FIELD(gx->mfield))
        printf("Mask loaded\n");

    //gwy_object_unref(container);

    return 0;
}

/*static inline
gboolean ggwy_data_field_inside(GwyDataField *data_field, gint i, gint j)
{
    if (i >= 0 && j >= 0 && i < data_field->xres && j < data_field->yres)
        return TRUE;
    else
        return FALSE;
}*/

static void
add_sphere(gdouble ***buffer, gint pi, gint pj, gint pk, gdouble factor, gdouble radius, gint xres, gint yres, gint zres)
{
    gint i, j, k;
    gdouble rdist;

    for (i = (gint)MAX(0, pi - radius); i < MIN(xres, pi + radius); i++) {
        for (j = (gint)MAX(0, pj - radius); j < MIN(xres, pj + radius); j++) {
            for (k = (gint)MAX(0, pk - radius); k < MIN(xres, pk + radius); k++) {
                rdist = MAX(0, (radius - sqrt((i - pi)*(i - pi) + (j - pj)*(j - pj) + (k - pk)*(k - pk))) / radius);
                buffer[i][j][k] += factor*rdist;
                //                printf("%d %d %d %g\n", i, j, k, factor*rdist);
            }
        }
    }
}

static gboolean
closetomat(SvICube *mat, gint pi, gint pj, gint pk, gint matindex, gdouble radius, gint xres, gint yres, gint zres)
{
    gint i, j, k;

    for (i = (gint)MAX(0, pi - radius); i < MIN(xres, pi + radius); i++) {
        for (j = (gint)MAX(0, pj - radius); j < MIN(xres, pj + radius); j++) {
            for (k = (gint)MAX(0, pk - radius); k < MIN(xres, pk + radius); k++) {
                if (mat->data[i][j][k] == matindex) return 1;
            }
        }
    }
    return 0;
}

/*
static gdouble
disttomat(SvICube *mat, gint pi, gint pj, gint pk, gint matindex, gdouble radius, gint xres, gint yres, gint zres)
{
    gint i, j, k;
    gdouble dist, mindist = G_MAXDOUBLE;

    for (i = MAX(0, pi-radius); i<MIN(xres, pi+radius); i++) {
        for (j = MAX(0, pj-radius); j<MIN(xres, pj+radius); j++) {
            for (k = MAX(0, pk-radius); k<MIN(xres, pk+radius); k++) {
                if (mat->data[i][j][k] == matindex) {
            dist = ((i-pi)*(i-pi) + (j-pj)*(j-pj) + (k-pk)*(k-pk));
            if (dist < mindist) mindist = dist;
                }
            }
        }
    }
    return sqrt(mindist);
}
*/

static void
addmat(gdouble ***buffer, gint pi, gint pj, gint pk, gdouble radius, gint xres, gint yres, gint zres)
{
    gint i, j, k;
    for (i = (gint)MAX(0, pi - radius); i < MIN(xres, pi + radius); i++) {
        for (j = (gint)MAX(0, pj - radius); j < MIN(xres, pj + radius); j++) {
            for (k = (gint)MAX(0, pk - radius); k < MIN(xres, pk + radius); k++) {
                if (buffer[i][j][k] < 1)
                    buffer[i][j][k] = CLAMP(MAX((radius - sqrt((i - pi)*(i - pi) + (j - pj)*(j - pj) + (k - pk)*(k - pk))), 0) / radius, buffer[i][j][k], 1);
            }
        }
    }
    buffer[pi][pj][pk] = 1;
}

static gboolean
closetovoid(SvICube *mat, gint pi, gint pj, gint pk, gint matindex, gdouble radius, gint xres, gint yres, gint zres)
{
    gint i, j, k;

    for (i = (gint)MAX(0, pi - radius); i < MIN(xres, pi + radius); i++) {
        for (j = (gint)MAX(0, pj - radius); j < MIN(xres, pj + radius); j++) {
            for (k = (gint)MAX(0, pk - radius); k < MIN(xres, pk + radius); k++)
                if (mat->data[i][j][k] != matindex)
                    return 1;
        }
    }
    return 0;
}

/*
gint evalnbs(gdouble ***buffer, gint posx, gint posy, gint posz, gint matindex, gint addindex)
{
    gint i, j, k;
    gint count = 0;

    for (i=(posx-1); i<=(posx+1); i++) {
        for (j=(posy-1); j<=(posy+1); j++) {
            for (k=(posz-1); k<=(posz+1); k++)
            if (((gint)buffer[i][j][k])==matindex || ((gint)buffer[i][j][k])==addindex) count++;
        }
    }
    return count;
}
*/

gint
evalnbs(gint ***buffer, gint posx, gint posy, gint posz, gint matindex, gint addindex, gint maxsofar)
{
    gint i, j, k;
    gint count = 0;

    for (i = (posx - 1); i <= (posx + 1); i++) {
        for (j = (posy - 1); j <= (posy + 1); j++) {
            for (k = (posz - 1); k <= (posz + 1); k++) {
                if (!(((gint)buffer[i][j][k]) == matindex || ((gint)buffer[i][j][k]) == addindex))
                    count++;
                if (count > maxsofar)
                    return 27 - count;
            }
        }
    }

    return 27 - count;
}

gboolean eval_closecheck(SvPool *mp, SvSet *set, gint ipos, gint jpos, gint kpos, gint skipdepth)
{
    gint i, j, k;

    for (i = MAX(0, ipos - skipdepth); i < MIN(ipos + skipdepth, set->sp.xres); i++) {
        for (j = MAX(0, jpos - skipdepth); j < MIN(jpos + skipdepth, set->sp.yres); j++) {
            for (k = MAX(0, kpos - skipdepth); k < MIN(kpos + skipdepth, set->sp.zres); k++) {
                if (sv_yee_data_get_tsf_mattype(mp->d, set, mp->mats, mp->nmat, i, j, k) != 0) {
                    //                       printf("skip %d %d %d\n", ipos, jpos, kpos);
                    return TRUE;
                }
            }
        }
    }

    return FALSE;
}

/*load vtk file, of exres, eyres, ezres is not zero, check if the size is the same.
Data must be allready allocated, otherwise it only reports size*/
gint load_vtk(gchar *filename, gint exres, gint eyres, gint ezres, gint *rxres, gint *ryres, gint *rzres, gint ***data)
{
    gint i, j, k, xres, yres, zres, nexpected;
    gchar buffer[100];
    FILE *fr;

    /*load sample data*/
    fr = fopen(filename, "r");
    if (!fr) {
        fprintf(stderr, "Error: can not open file %s\n", filename);
        return 1;
    }

    fgets(buffer, sizeof(buffer), fr);
    fgets(buffer, sizeof(buffer), fr);
    fscanf(fr, "%s", buffer);
    if (strcmp(buffer, "ASCII") != 0) {
        fprintf(stderr, "Error in VTK file, it should be ascii\n");
        fclose(fr);
        return 1;
    }
    fscanf(fr, "%s", buffer);
    fscanf(fr, "%s", buffer);
    if (strcmp(buffer, "STRUCTURED_POINTS") != 0) {
        fprintf(stderr, "Error in VTK file, it should be structured points dataset\n");
        fclose(fr);
        return 1;
    }

    fscanf(fr, "%s", buffer);
    if (strcmp(buffer, "DIMENSIONS") != 0) {
        fprintf(stderr, "Error in VTK file, DIMENSIONS keyword is expected here\n");
        fclose(fr);
        return 1;
    }
    fscanf(fr, "%d", &xres);
    fscanf(fr, "%d", &yres);
    fscanf(fr, "%d", &zres);
    if ((exres + eyres + ezres) != 0) {
        if (xres != exres || yres != eyres || zres != ezres) {
            fprintf(stderr, "Error in VTK file, Dimensions do not match expected values\n");
            fclose(fr);
            return 1;
        }
    }
    fgets(buffer, sizeof(buffer), fr);
    fgets(buffer, sizeof(buffer), fr);
    fgets(buffer, sizeof(buffer), fr);
    fscanf(fr, "%s", buffer);
    if (strcmp(buffer, "POINT_DATA") != 0) {
        fprintf(stderr, "Error in VTK file, POINT_DATA keyword is expected here\n");
        fclose(fr);
        return 1;
    }
    fscanf(fr, "%d", &nexpected);
    if (nexpected != (xres*yres*zres)) {
        fprintf(stderr, "Error in VTK file, POINT_DATA does not match dimensions\n");
        fclose(fr);
        return 1;
    }
    fgets(buffer, sizeof(buffer), fr);
    fgets(buffer, sizeof(buffer), fr);
    fgets(buffer, sizeof(buffer), fr);
    if (data != NULL) {
        for (k = 0; k < zres; k++) {
            for (j = 0; j < yres; j++) {
                for (i = 0; i < xres; i++) {
                    if (fscanf(fr, "%d", &(data[i][j][k])) != 1) {
                        fprintf(stderr, "Error in VTK file, Error while reading datapoints\n");
                        fclose(fr);
                        return 1;
                    }
                }
            }
        }
    }
    *rxres = xres;
    *ryres = yres;
    *rzres = zres;

    return 0;
}

static void
allocate_dplrc(SvYeeData *d, gint np)
{
    gint j;

    d->dplrcx = (SvDCube**)g_malloc(np * sizeof(SvDCube*));
    d->dplrcy = (SvDCube**)g_malloc(np * sizeof(SvDCube*));
    d->dplrcz = (SvDCube**)g_malloc(np * sizeof(SvDCube*));

    for (j = 0; j < np; j++) {
       d->dplrcx[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
       d->dplrcy[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
       d->dplrcz[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
    }
}

static void
allocate_plrc(SvYeeData *d, gint np) 
{
    gint j;

    d->plrcx = (SvDCube**)g_malloc(np * sizeof(SvDCube*));
    d->plrcy = (SvDCube**)g_malloc(np * sizeof(SvDCube*));
    d->plrcz = (SvDCube**)g_malloc(np * sizeof(SvDCube*));

    for (j = 0; j < np; j++) {
       d->plrcx[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
       d->plrcy[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
       d->plrcz[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
    }
    d->exp = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
    d->eyp = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
    d->ezp = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);

}

static void
allocate_iplrc(SvYeeData *d, gint np)
{
    gint j;

    d->iplrcx = (SvDCube**)g_malloc(np * sizeof(SvDCube*));
    d->iplrcy = (SvDCube**)g_malloc(np * sizeof(SvDCube*));
    d->iplrcz = (SvDCube**)g_malloc(np * sizeof(SvDCube*));

    for (j = 0; j < np; j++) {
       d->iplrcx[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
       d->iplrcy[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
       d->iplrcz[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
    }
}

static void 
allocate_ade(SvYeeData *d)
{
    gint j, np;

    d->exp = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
    d->eyp = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
    d->ezp = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);

    np = 1;

    d->dpxp = (SvDCube**)g_malloc(np * sizeof(SvDCube*));
    d->dpyp = (SvDCube**)g_malloc(np * sizeof(SvDCube*));
    d->dpzp = (SvDCube**)g_malloc(np * sizeof(SvDCube*));

    for (j = 0; j < np; j++) {
       d->dpxp[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
       d->dpyp[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
       d->dpzp[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
    }

    d->dpx = (SvDCube**)g_malloc(np * sizeof(SvDCube*));
    d->dpy = (SvDCube**)g_malloc(np * sizeof(SvDCube*));
    d->dpz = (SvDCube**)g_malloc(np * sizeof(SvDCube*));

    for (j = 0; j < np; j++) {
       d->dpx[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
       d->dpy[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
       d->dpz[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
    }

    np = 2;

    d->pxp = (SvDCube**)g_malloc(np * sizeof(SvDCube*));
    d->pyp = (SvDCube**)g_malloc(np * sizeof(SvDCube*));
    d->pzp = (SvDCube**)g_malloc(np * sizeof(SvDCube*));

    for (j = 0; j < np; j++) {
       d->pxp[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
       d->pyp[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
       d->pzp[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
    }

    d->px = (SvDCube**)g_malloc(np * sizeof(SvDCube*));
    d->py = (SvDCube**)g_malloc(np * sizeof(SvDCube*));
    d->pz = (SvDCube**)g_malloc(np * sizeof(SvDCube*));

    for (j = 0; j < np; j++) {
       d->px[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
       d->py[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
       d->pz[j] = sv_dcube_new(d->xres, d->yres, d->zres, d->xres * d->dx, d->yres * d->dy, d->zres * d->dz, 1);
    }



}


gint
load_data(SvPool *mp, SvSet *set)
{
    gint i, j, k, n, ii, jj, kk, it, xres, yres, zres, np, matindex, nmatindex, mi;
    guint mui;
    gint type;
    SvSphere sx, ssx;
    SvVoxel vx, svx;
    SvCylinder cx;
    SvCone cnx;
    SvRCone rcnx;
    SvTetrahedron tx;
    SvGwydd gx;
    GArray *spheres;
    GArray *voxels;
    GArray *cylinders;
    GArray *cones;
    GArray *rcones;
    GArray *tetrahedrons;
    GArray *gwydds;
    gint bnd_xmin, bnd_xmax, bnd_ymin, bnd_ymax, bnd_zmin, bnd_zmax;
    gint npos, ngtot = 0;
    FILE *fr;
    gboolean plrc, plrci, dplrc, ade;
    SvMatProp mat;
    GLfloat xshift, yshift, zshift, xmult, ymult, zmult;
    gint gwydd_channel = 0, mask;
    gdouble ***buffer;
    gint ***ibuffer;
    gchar filebase[256];
    gint attribute_pos = 0, attribute_val = 0;
    gdouble xreal, yreal, zreal, val;
    GRand *rnd;
    gint gxres, gyres, gzres, sub, face;
    gdouble dx, dy, dz, dd, posx, posy, posz, maxdd;
    gint nbdist, nbscore, nbbest, ni, nj, nk, iposx, iposy, iposz;
    gdouble probability, theta, phi;
    GwyExpr *expr;
    const gchar *const var_names[] = {"x", "y", "z"};
    guint var_positions[G_N_ELEMENTS(var_names)];
    gdouble vars[G_N_ELEMENTS(var_names) + 1];
    gint domat, jt;
    gint landx[100], landy[100], landz[100], land[100];
    const gboolean entity_antialiasing = set->sm.smooth; // experimental; default FALSE, now set to MEDIUM_SMOOTH directive
    Submesh *sm;
    GTree *submeshes;
    MpSetData data_mp_set;
    gdouble epsilon, sigma, nepsilon, nsigma;
    gint vxres, vyres, vzres;
    SvSg *sg;

    probability = theta = phi = 0;
    dx = dy = dz = dd = posx = posy = posz = maxdd = 0;
    nbdist = nbscore = nbbest = ni = nj = nk = iposx = iposy = iposz = 0;
    gxres = gyres = gzres = sub = face = 0;

#ifndef G_OS_WIN32
    struct timeval time;
#endif

    if (set->sm.in_voxel_filename != NULL) {
        fr = fopen(set->sm.in_voxel_filename, "rb");
        if (fr == NULL) {
            fprintf(stderr, "Error: cannot open material file %s for reading\n", set->sm.in_voxel_filename);
        }

        if (set->plan.matmode != SV_MATMODE_NONE && fr) {
            if (set->sc.verbose > 1)
                printf("Loading MEDIUM_LINEAR data...\n");

            fread(&xres, sizeof(gint), 1, fr);
            fread(&yres, sizeof(gint), 1, fr);
            fread(&zres, sizeof(gint), 1, fr);
            if (xres != set->sp.xres || yres != set->sp.yres || zres != set->sp.zres) {
                fprintf(stderr, "Error: material dimensions differ from pool settings\n");
                return 1;
            }

            if (set->plan.matmode == SV_MATMODE_MAGNETIC) {
                fseek(fr, 2 * xres*yres*zres * sizeof(gfloat), SEEK_SET);
            } else {
                for (i = 0; i < xres; i++) {
                    for (j = 0; j < yres; j++)
                        fread((mp->d->epsilon->data[i][j]), sizeof(gfloat), zres, fr);
                }
                sv_fcube_multiply(mp->d->epsilon, (gfloat)EPSILON_0);

                for (i = 0; i < xres; i++) {
                    for (j = 0; j < yres; j++)
                        fread(mp->d->sigma->data[i][j], sizeof(gfloat), zres, fr);
                }
            }

            if (set->plan.matmode == SV_MATMODE_ELECTRIC) {
                fclose(fr);
                return 0;
            }

            for (i = 0; i < xres; i++) {
                for (j = 0; j < yres; j++)
                    fread(mp->d->mu->data[i][j], sizeof(gfloat), zres, fr);
            }
            sv_fcube_multiply(mp->d->mu, (gfloat)MU_0);

            for (i = 0; i < xres; i++) {
                for (j = 0; j < yres; j++)
                    fread(mp->d->sigast->data[i][j], sizeof(gfloat), zres, fr);
            }

            fclose(fr);
        }
    }

    if (set->sm.in_vector_filename != NULL) {
        fr = fopen(set->sm.in_vector_filename, "r");
        if (fr == NULL)
            fprintf(stderr, "Error: cannot open material vector file for loading.\n");
        else {
            if (set->sc.verbose > 1)
                printf("Loading MEDIUM_VECTOR data...\n");
            xres = set->sp.xres;
            yres = set->sp.yres;
            zres = set->sp.zres;

            spheres = g_array_new(FALSE, FALSE, sizeof(SvSphere));
            cones = g_array_new(FALSE, FALSE, sizeof(SvCone));
            rcones = g_array_new(FALSE, FALSE, sizeof(SvRCone));
            cylinders = g_array_new(FALSE, FALSE, sizeof(SvCylinder));
            voxels = g_array_new(FALSE, FALSE, sizeof(SvVoxel));
            tetrahedrons = g_array_new(FALSE, FALSE, sizeof(SvTetrahedron));
            gwydds = g_array_new(FALSE, FALSE, sizeof(SvGwydd));

            while (fscanf(fr, "%d", &type) != EOF) { //TODO use proper locale handling functions here
                switch (type) {
                    case SV_ENTITY_SPHERE:
                        scan_point(fr, sx.pnt1, sx.pnt1 + 1, sx.pnt1 + 2);
                        scan_radius(fr, &sx.radius);
                        scan_material(fr, &sx.mat, set->plan.dt);
                        if (sx.mat.type != SV_MAT_LINEAR) {
                            if ((matindex = sv_pool_get_matindex(mp, sx.mat)) == -1)
                                sv_pool_add_material(mp, &sx.mat, FALSE);
                            else
                                sx.mat.pos = matindex;
                        }
                        sx.n = ngtot++;
                        g_array_append_val(spheres, sx);
                        break;

                    case SV_ENTITY_CYLINDER:
                        scan_point(fr, cx.pnt1, cx.pnt1 + 1, cx.pnt1 + 2);
                        scan_point(fr, cx.pnt2, cx.pnt2 + 1, cx.pnt2 + 2);
                        scan_radius(fr, &cx.radius);
                        scan_material(fr, &cx.mat, set->plan.dt);
                        if (cx.mat.type != SV_MAT_LINEAR) {
                            if ((matindex = sv_pool_get_matindex(mp, cx.mat)) == -1)
                                sv_pool_add_material(mp, &cx.mat, FALSE);
                            else cx.mat.pos = matindex;
                        }
                        cx.n = ngtot++;
                        g_array_append_val(cylinders, cx);
                        break;

                    case SV_ENTITY_CONE:
                        scan_point(fr, cnx.pnt1, cnx.pnt1 + 1, cnx.pnt1 + 2);
                        scan_point(fr, cnx.pnt2, cnx.pnt2 + 1, cnx.pnt2 + 2);
                        scan_radius(fr, &cnx.radius);
                        scan_material(fr, &cnx.mat, set->plan.dt);
                        if (cnx.mat.type != SV_MAT_LINEAR) {
                            if ((matindex = sv_pool_get_matindex(mp, cnx.mat)) == -1)
                                sv_pool_add_material(mp, &cnx.mat, FALSE);
                            else
                                cnx.mat.pos = matindex;
                        }
                        cnx.n = ngtot++;
                        g_array_append_val(cones, cnx);
                        break;

                    case SV_ENTITY_RCONE:
                        scan_point(fr, rcnx.pnt1, rcnx.pnt1 + 1, rcnx.pnt1 + 2);
                        scan_point(fr, rcnx.pnt2, rcnx.pnt2 + 1, rcnx.pnt2 + 2);
                        scan_radius(fr, &rcnx.radius1);
                        scan_radius(fr, &rcnx.radius2);
                        scan_material(fr, &rcnx.mat, set->plan.dt);
                        if (rcnx.mat.type != SV_MAT_LINEAR) {
                            if ((matindex = sv_pool_get_matindex(mp, rcnx.mat)) == -1)
                                sv_pool_add_material(mp, &rcnx.mat, FALSE);
                            else
                                rcnx.mat.pos = matindex;
                        }
                        rcnx.n = ngtot++;
                        g_array_append_val(rcones, rcnx);
                        break;

                    case SV_ENTITY_VOXEL:
                        scan_point(fr, vx.pnt1, vx.pnt1 + 1, vx.pnt1 + 2);
                        scan_point(fr, vx.pnt2, vx.pnt2 + 1, vx.pnt2 + 2);
                        scan_material(fr, &vx.mat, set->plan.dt);
                        if (vx.mat.type != SV_MAT_LINEAR) {
                            if ((matindex = sv_pool_get_matindex(mp, vx.mat)) == -1)
                                sv_pool_add_material(mp, &vx.mat, FALSE);
                            else
                                vx.mat.pos = matindex;
                        }
                        vx.n = ngtot++;
                        g_array_append_val(voxels, vx);
                        break;

                    case SV_ENTITY_TETRAHEDRON:
                        scan_point(fr, tx.pnt1, tx.pnt1 + 1, tx.pnt1 + 2);
                        scan_point(fr, tx.pnt2, tx.pnt2 + 1, tx.pnt2 + 2);
                        scan_point(fr, tx.pnt3, tx.pnt3 + 1, tx.pnt3 + 2);
                        scan_point(fr, tx.pnt4, tx.pnt4 + 1, tx.pnt4 + 2);
                        scan_material(fr, &tx.mat, set->plan.dt);
                        if (tx.mat.type != SV_MAT_LINEAR) {
                            if ((matindex = sv_pool_get_matindex(mp, tx.mat)) == -1)
                                sv_pool_add_material(mp, &tx.mat, FALSE);
                            else
                                tx.mat.pos = matindex;
                        }
                        tx.n = ngtot++;
                        g_array_append_val(tetrahedrons, tx);
                        break;

                    case SV_ENTITY_TETGEN:
                        fscanf(fr, "%s", filebase);             //TODO this should be locale dependent
                        scan_int(fr, &attribute_pos);
                        scan_int(fr, &attribute_val);
                        scan_point(fr, &xshift, &yshift, &zshift);
                        scan_point(fr, &xmult, &ymult, &zmult);
                        scan_material(fr, &mat, set->plan.dt);
                        if (mat.type != SV_MAT_LINEAR) {
                            if ((matindex = sv_pool_get_matindex(mp, mat)) == -1)
                                sv_pool_add_material(mp, &mat, FALSE);
                            else
                                mat.pos = matindex;
                        }
                        load_tets(filebase, tetrahedrons, attribute_pos, attribute_val, xshift, yshift, zshift, xmult, ymult, zmult, set->sc.verbose > 1, &mat, &ngtot);
                        break;

                    case SV_ENTITY_GWYDDION:
                        fscanf(fr, "%s", filebase);             //TODO this should be locale dependent
                        scan_int(fr, &gwydd_channel);
                        scan_int(fr, &gx.mask);
                        scan_int(fr, &gx.i);
                        scan_int(fr, &gx.j);
                        scan_int(fr, &gx.k);
                        scan_point(fr, &gx.xoffset, &gx.yoffset, &gx.zoffset);
                        scan_int(fr, &gx.depth);
                        scan_material(fr, &gx.mat, set->plan.dt);
                        if (gx.mat.type != SV_MAT_LINEAR) {
                            if ((matindex = sv_pool_get_matindex(mp, gx.mat)) == -1)
                                sv_pool_add_material(mp, &gx.mat, FALSE);
                            else
                                gx.mat.pos = matindex;
                        }
                        load_gwydd(filebase, &gx, set->sc.verbose > 1, gwydd_channel);

                        gx.n = ngtot++;
                        g_array_append_val(gwydds, gx);
                        break;

                    default:
                        break;
                }
            }

            /*alloc tabulated fields if necessary*/
            if (mp->nmat != 0) {
                mp->d->mat = sv_icube_new(set->sp.xres, set->sp.yres, set->sp.zres, set->sp.xres*set->sp.dx, set->sp.yres*set->sp.dy, set->sp.zres*set->sp.dz, 1);
                sv_icube_fill(mp->d->mat, 0); // init with vacuum


                if (mp->sg && mp->sg->nsubgrid > 0) {
                   for (n=0; n<mp->sg->nsubgrid; n++)
                   {
                        sg = mp->sg->sg[n];
                        sg->d->mat = sv_icube_new(sg->d->xres, sg->d->yres, sg->d->zres, sg->d->xres*sg->d->dx, sg->d->yres*sg->d->dy, sg->d->zres*sg->d->dz, 1);
                        sv_icube_fill(sg->d->mat, 0);
                   }
                }

            }

            //if vtk data are available, load it now 
            if (set->sm.in_vtk_data_filename != NULL) {
                if (load_vtk(set->sm.in_vtk_data_filename, set->sp.xres, set->sp.yres, set->sp.zres, &vxres, &vyres, &vzres, NULL))
                    return 1;
                else
                    load_vtk(set->sm.in_vtk_data_filename, set->sp.xres, set->sp.yres, set->sp.zres, &vxres, &vyres, &vzres, mp->d->mat->data);
            }


            /*now process the arrays*/
            npos = 0;

            if (entity_antialiasing)
                printf("Geometry antialiasing on, processing can take some time...\n");
            submeshes = g_tree_new((GCompareFunc)submeshcmp);

            do {
                for (mui = 0; mui < spheres->len; mui++) {
                    sx = g_array_index(spheres, SvSphere, mui);

                    if (npos != sx.n)
                        continue;
                    else
                        npos++;

                    for (i = (gint)MAX(0, sx.pnt1[0] - sx.radius - 1); i < MIN(xres, sx.pnt1[0] + sx.radius + 1); i++) { //x
                        for (j = (gint)MAX(0, sx.pnt1[1] - sx.radius - 1); j < MIN(yres, sx.pnt1[1] + sx.radius + 1); j++) { //y
                            for (k = (gint)MAX(0, sx.pnt1[2] - sx.radius - 1); k < MIN(zres, sx.pnt1[2] + sx.radius + 1); k++) { //z
                                if (entity_antialiasing) {
                                    sm = init_submesh(i, j, k, submeshes, mp, set);
                                    add_mat_to_submesh_sphere(i, j, k, sx.pnt1, sx.radius, sm, sx.mat);
                                    check_and_note_mixed(submeshes, sm, mp, set);
                                } else { // neni antialiasing: je-li stred ve sfere, prdnu to tam
                                    if (is_in_sphere(i + 0.5f, j + 0.5f, k + 0.5f, sx.pnt1, sx.radius))
                                        sv_pool_set_data(mp, sx.mat, set, i, j, k);
                                }
                            } // k
                        } // j
                    } // i

                  if (mp->sg && mp->sg->nsubgrid > 0) {

                        for (n=0; n<mp->sg->nsubgrid; n++) 
                        {
                            sg = mp->sg->sg[n];
                            printf("adjusting material to subgrid %d\n", n);

                            //check if the object is inside subgrid
                            if (!((sx.pnt1[0]+sx.radius)>=sg->ifrom && (sx.pnt1[1]+sx.radius)>=sg->jfrom && (sx.pnt1[2]+sx.radius)>=sg->kfrom 
                                  && (sx.pnt1[0]-sx.radius)<sg->ito &&  (sx.pnt1[1]-sx.radius)<sg->jto && (sx.pnt1[2]-sx.radius)<sg->kto)) continue;

                            ssx.pnt1[0] = (vx.pnt1[0] - sg->ifrom)*sg->division;
                            ssx.pnt1[1] = (vx.pnt1[1] - sg->jfrom)*sg->division;
                            ssx.pnt1[2] = (vx.pnt1[2] - sg->kfrom)*sg->division;
                            ssx.radius = sx.radius*sg->division;

                            for (i = (gint)MAX(0, ssx.pnt1[0] - ssx.radius - 1); i < MIN(sg->d->xres, ssx.pnt1[0] + ssx.radius + 1); i++) { //x
                               for (j = (gint)MAX(0, ssx.pnt1[1] - ssx.radius - 1); j < MIN(sg->d->yres, ssx.pnt1[1] + ssx.radius + 1); j++) { //y
                                  for (k = (gint)MAX(0, ssx.pnt1[2] - ssx.radius - 1); k < MIN(sg->d->zres, ssx.pnt1[2] + ssx.radius + 1); k++) { //z
                                    if (is_in_sphere(i + 0.5f, j + 0.5f, k + 0.5f, ssx.pnt1, ssx.radius))
                                         sv_yee_set_data(sg->d, sx.mat, set, i, j, k);
                                  } // k
                               } // j
                            } // i
                        } 

                    }


                } //spheres

                //    printf("we have %d tetrahedrons\n", tetrahedrons->len);
                for (mui = 0; mui < tetrahedrons->len; mui++) {
                    tx = g_array_index(tetrahedrons, SvTetrahedron, mui);

                    if (npos != tx.n)
                        continue;
                    else
                        npos++;

                    //      printf("th: (%g %g %g) (%g %g %g) (%g %g %g) (%g %g %g) %d\n", tx.pnt1[0], tx.pnt1[1], tx.pnt1[2], tx.pnt2[0], tx.pnt2[1], tx.pnt2[2],
                    //           tx.pnt3[0], tx.pnt3[1], tx.pnt3[2], tx.pnt4[0], tx.pnt4[1], tx.pnt4[2], tx.mat.type);

                    bnd_xmin = (gint)MIN(MIN(tx.pnt1[0], tx.pnt2[0]), MIN(tx.pnt3[0], tx.pnt4[0]));
                    bnd_xmin = CLAMP(bnd_xmin, 0, xres - 1);
                    bnd_ymin = (gint)MIN(MIN(tx.pnt1[1], tx.pnt2[1]), MIN(tx.pnt3[1], tx.pnt4[1]));
                    bnd_ymin = CLAMP(bnd_ymin, 0, yres - 1);
                    bnd_zmin = (gint)MIN(MIN(tx.pnt1[2], tx.pnt2[2]), MIN(tx.pnt3[2], tx.pnt4[2]));
                    bnd_zmin = CLAMP(bnd_zmin, 0, zres - 1);

                    bnd_xmax = (gint)MAX(MAX(tx.pnt1[0], tx.pnt2[0]), MAX(tx.pnt3[0], tx.pnt4[0]));
                    bnd_xmax = CLAMP(bnd_xmax, 0, xres - 1);
                    bnd_ymax = (gint)MAX(MAX(tx.pnt1[1], tx.pnt2[1]), MAX(tx.pnt3[1], tx.pnt4[1]));
                    bnd_ymax = CLAMP(bnd_ymax, 0, yres - 1);
                    bnd_zmax = (gint)MAX(MAX(tx.pnt1[2], tx.pnt2[2]), MAX(tx.pnt3[2], tx.pnt4[2]));
                    bnd_zmax = CLAMP(bnd_zmax, 0, zres - 1);

                    for (i = bnd_xmin; i <= bnd_xmax; i++) { //x
                        for (j = bnd_ymin; j <= bnd_ymax; j++) { //y
                            for (k = bnd_zmin; k <= bnd_zmax; k++) { //z
                                if (entity_antialiasing) {
                                    sm = init_submesh(i, j, k, submeshes, mp, set);
                                    add_mat_to_submesh_tetrahedron(i, j, k, tx.pnt1, tx.pnt2, tx.pnt3, tx.pnt4, sm, tx.mat);
                                    check_and_note_mixed(submeshes, sm, mp, set);
                                } else {
                                    if (is_in_tetrahedron(i + 0.5f, j + 0.5f, k + 0.5f, tx.pnt1, tx.pnt2, tx.pnt3, tx.pnt4))
                                        sv_pool_set_data(mp, tx.mat, set, i, j, k);
                                }
                            } // k
                        } // j
                    } // i
                }

                for (mui = 0; mui < voxels->len; mui++) {
                    vx = g_array_index(voxels, SvVoxel, mui);

                    if (npos != vx.n)
                        continue;
                    else
                        npos++;

                    for (i = (gint)MAX(0, vx.pnt1[0]); i < MIN(xres, vx.pnt2[0]); i++) { //x
                        for (j = (gint)MAX(0, vx.pnt1[1]); j < MIN(yres, vx.pnt2[1]); j++) { //y
                            for (k = (gint)MAX(0, vx.pnt1[2]); k < MIN(zres, vx.pnt2[2]); k++)
                                sv_pool_set_data(mp, vx.mat, set, i, j, k);
                        }
                    }

                    if (mp->sg && mp->sg->nsubgrid > 0) {

                        for (n=0; n<mp->sg->nsubgrid; n++) 
                        {
                            sg = mp->sg->sg[n];
                            printf("adjusting material to subgrid %d\n", n);

                            //check if the object is inside subgrid
                            if (!(vx.pnt1[0]<=sg->ito && vx.pnt2[0]>=sg->ifrom && vx.pnt1[1]<=sg->jto && vx.pnt2[1]>=sg->jfrom && vx.pnt1[2]<=sg->kto && vx.pnt2[2]>=sg->kfrom)) continue;

                            svx.pnt1[0] = (vx.pnt1[0] - sg->ifrom)*sg->division;
                            svx.pnt2[0] = (vx.pnt2[0] - sg->ifrom)*sg->division;

                            svx.pnt1[1] = (vx.pnt1[1] - sg->jfrom)*sg->division;
                            svx.pnt2[1] = (vx.pnt2[1] - sg->jfrom)*sg->division;

                            svx.pnt1[2] = (vx.pnt1[2] - sg->kfrom)*sg->division;
                            svx.pnt2[2] = (vx.pnt2[2] - sg->kfrom)*sg->division;

                            printf("new fill: %g %g %g    %g %g %g\n", svx.pnt1[0], svx.pnt1[1], svx.pnt1[2], svx.pnt2[0], svx.pnt2[1], svx.pnt2[2]);


                            for (i = (gint)MAX(0, svx.pnt1[0]); i < MIN(sg->d->xres, svx.pnt2[0]); i++) { //x
                                for (j = (gint)MAX(0, svx.pnt1[1]); j < MIN(sg->d->yres, svx.pnt2[1]); j++) { //y
                                    for (k = (gint)MAX(0, svx.pnt1[2]); k < MIN(sg->d->zres, svx.pnt2[2]); k++)

                                        sv_yee_set_data(sg->d, vx.mat, set, i, j, k);
                                }
                            }

                        } 

                    }

                }

                for (mui = 0; mui < cylinders->len; mui++) {
                    cx = g_array_index(cylinders, SvCylinder, mui);

                    if (npos != cx.n)
                        continue;
                    else
                        npos++;

                    for (i = (gint)MAX(0, MIN(cx.pnt1[0], cx.pnt2[0]) - cx.radius - 1); i < MIN(xres, MAX(cx.pnt1[0], cx.pnt2[0]) + cx.radius + 1); i++) { //x
                        for (j = (gint)MAX(0, MIN(cx.pnt1[1], cx.pnt2[1]) - cx.radius - 1); j < MIN(yres, MAX(cx.pnt1[1], cx.pnt2[1]) + cx.radius + 1); j++) { //y
                            for (k = (gint)MAX(0, MIN(cx.pnt1[2], cx.pnt2[2]) - cx.radius - 1); k < MIN(zres, MAX(cx.pnt1[2], cx.pnt2[2]) + cx.radius + 1); k++) { //z
                                if (entity_antialiasing) {
                                    sm = init_submesh(i, j, k, submeshes, mp, set);
                                    add_mat_to_submesh_cylinder(i, j, k, cx.pnt1, cx.pnt2, cx.radius, sm, cx.mat);
                                    check_and_note_mixed(submeshes, sm, mp, set);
                                } else {
                                    if (is_in_cylinder(i + 0.5f, j + 0.5f, k + 0.5f, cx.pnt1, cx.pnt2, cx.radius)) {
                                        sv_pool_set_data(mp, cx.mat, set, i, j, k);
                                        //if (i==64 && j==116 && k==48) printf("XXXXXXXXXXXXXXXXXXXX\n");
                                    }
                                }
                            } // k
                        }  // j
                    }  // i
                }

                for (mui = 0; mui < cones->len; mui++) {
                    cnx = g_array_index(cones, SvCone, mui);

                    if (npos != cnx.n)
                        continue;
                    else
                        npos++;

                    for (i = (gint)MAX(0, MIN(cnx.pnt1[0], cnx.pnt2[0]) - cnx.radius - 1); i < MIN(xres, MAX(cnx.pnt1[0], cnx.pnt2[0]) + cnx.radius + 1); i++) { //x
                        for (j = (gint)MAX(0, MIN(cnx.pnt1[1], cnx.pnt2[1]) - cnx.radius - 1); j < MIN(yres, MAX(cnx.pnt1[1], cnx.pnt2[1]) + cnx.radius + 1); j++) { //y
                            for (k = (gint)MAX(0, MIN(cnx.pnt1[2], cnx.pnt2[2]) - cnx.radius - 1); k < MIN(zres, MAX(cnx.pnt1[2], cnx.pnt2[2]) + cnx.radius + 1); k++) { //z
                                if (entity_antialiasing) {
                                    sm = init_submesh(i, j, k, submeshes, mp, set);
                                    add_mat_to_submesh_cone(i, j, k, cnx.pnt1, cnx.pnt2, cnx.radius, sm, cnx.mat);
                                    check_and_note_mixed(submeshes, sm, mp, set);
                                } else {
                                    if (is_in_cone(i + 0.5f, j + 0.5f, k + 0.5f, cnx.pnt1, cnx.pnt2, cnx.radius))
                                        sv_pool_set_data(mp, cnx.mat, set, i, j, k);
                                }
                            } // k
                        } // j
                    } // i
                }

                for (mui = 0; mui < rcones->len; mui++) {
                    rcnx = g_array_index(rcones, SvRCone, mui);

                    if (npos != rcnx.n)
                        continue;
                    else
                        npos++;

                    for (i = (gint)MAX(0, MIN(rcnx.pnt1[0], rcnx.pnt2[0]) - MAX(rcnx.radius1, rcnx.radius2) - 1);
                         i < MIN(xres, MAX(rcnx.pnt1[0], rcnx.pnt2[0]) + MAX(rcnx.radius1, rcnx.radius2) + 1); i++) { //x
                        for (j = (gint)MAX(0, MIN(rcnx.pnt1[1], rcnx.pnt2[1]) - MAX(rcnx.radius1, rcnx.radius2) - 1);
                             j < MIN(yres, MAX(rcnx.pnt1[1], rcnx.pnt2[1]) + MAX(rcnx.radius1, rcnx.radius2) + 1); j++) { //y
                            for (k = (gint)MAX(0, MIN(rcnx.pnt1[2], rcnx.pnt2[2]) - MAX(rcnx.radius1, rcnx.radius2) - 1);
                                 k < MIN(zres, MAX(rcnx.pnt1[2], rcnx.pnt2[2]) + MAX(rcnx.radius1, rcnx.radius2) + 1); k++) { //z
                                if (entity_antialiasing) {
                                    sm = init_submesh(i, j, k, submeshes, mp, set);
                                    add_mat_to_submesh_rcone(i, j, k, rcnx.pnt1, rcnx.pnt2, rcnx.radius1, rcnx.radius2, sm, rcnx.mat);
                                    check_and_note_mixed(submeshes, sm, mp, set);
                                } else {
                                    if (is_in_rcone(i + 0.5f, j + 0.5f, k + 0.5f, rcnx.pnt1, rcnx.pnt2, rcnx.radius1, rcnx.radius2))
                                        sv_pool_set_data(mp, rcnx.mat, set, i, j, k);
                                }
                            } // k
                        } // j
                    } // i
                }

                for (mui = 0; mui < gwydds->len; mui++) {
                    gx = g_array_index(gwydds, SvGwydd, mui);
                    if (npos != gx.n)
                        continue;
                    else
                        npos++;

                    if (!gx.dfield)
                        printf("gwyddion field not found, skipped\n");
                    else {
                        //printf("gwyddion field in pos %d %d %d, depth %d, mask %d\n", gx.i, gx.j, gx.k, gx.depth, gx.mask);

                        xres = set->sp.xres;
                        xreal = set->sp.dx*xres;

                        yres = set->sp.yres;
                        yreal = set->sp.dy*yres;

                        zres = set->sp.zres;
                        zreal = set->sp.dz*zres;

                        //printf("cube: %d %d %d, %g %g %g, dfield %d %d, %g %g, minmax %g %g\n", xres, yres, zres, xreal, yreal, zreal, gxres, gyres, gxreal, gyreal, min, max);

                        for (i = 0; i < xres; i++) { //x
                            for (j = 0; j < yres; j++) { //y
                                for (k = 0; k < zres; k++) { //z
                                    mask = 0;

                                    if (gx.i != -1 && ggwy_data_field_inside(gx.dfield, (gint)(j - gx.yoffset / set->sp.dy), (gint)(k - gx.zoffset / set->sp.dz))) {
                                        val = gx.i + (gdouble)xres / xreal * (gwy_data_field_get_dval_real(gx.dfield, j*set->sp.dy - gx.yoffset, k*set->sp.dz - gx.zoffset, GWY_INTERPOLATION_BILINEAR) - gx.xoffset);
                                        if (gx.mask == 1 && gx.mfield)
                                            mask = (gint)gwy_data_field_get_dval_real(gx.mfield, (j*set->sp.dy - gx.xoffset), (k*set->sp.dz - gx.yoffset), GWY_INTERPOLATION_ROUND);
                                        else if (gx.mask == -1 && gx.mfield)
                                            mask = (gint)fabs(gwy_data_field_get_dval_real(gx.mfield, j*set->sp.dy - gx.xoffset, k*set->sp.dz - gx.yoffset, GWY_INTERPOLATION_ROUND) - 1);
                                        if (gx.depth > 0 && mask == 0 && i > val && i < (val + gx.depth))
                                            sv_pool_set_data(mp, gx.mat, set, i, j, k);
                                    }

                                    if (gx.j != -1 && ggwy_data_field_inside(gx.dfield, (gint)(i - gx.xoffset / set->sp.dx), (gint)(k - gx.zoffset / set->sp.dz))) {
                                        val = gx.j + (gdouble)yres / yreal * (gwy_data_field_get_dval_real(gx.dfield, i*set->sp.dx - gx.xoffset, k*set->sp.dz - gx.zoffset, GWY_INTERPOLATION_BILINEAR) - gx.yoffset);
                                        if (gx.mask == 1 && gx.mfield)
                                            mask = (gint)gwy_data_field_get_dval_real(gx.mfield, (gint)(i*set->sp.dx - gx.xoffset), (gint)(k*set->sp.dz - gx.yoffset), GWY_INTERPOLATION_ROUND);
                                        else if (gx.mask == -1 && gx.mfield)
                                            mask = (gint)fabs(gwy_data_field_get_dval_real(gx.mfield, (gint)(i*set->sp.dx - gx.xoffset), (gint)(k*set->sp.dz - gx.yoffset), GWY_INTERPOLATION_ROUND) - 1);
                                        if (gx.depth > 0 && mask == 0 && j > val && j < (val + gx.depth))
                                            sv_pool_set_data(mp, gx.mat, set, i, j, k);
                                    }

                                    if (gx.k != -1 && ggwy_data_field_inside(gx.dfield, (gint)(i - gx.xoffset / set->sp.dx), (gint)(j - gx.yoffset / set->sp.dy))) {
                                        val = gx.k + (gdouble)zres / zreal * (gwy_data_field_get_dval_real(gx.dfield, i*set->sp.dx - gx.xoffset, j*set->sp.dy - gx.yoffset, GWY_INTERPOLATION_BILINEAR) - gx.zoffset);
                                        if (gx.mask == 1 && gx.mfield)
                                            mask = (gint)gwy_data_field_get_dval_real(gx.mfield, i*set->sp.dx - gx.xoffset, j*set->sp.dy - gx.yoffset, GWY_INTERPOLATION_ROUND);
                                        else if (gx.mask == -1 && gx.mfield)
                                            mask = (gint)fabs(gwy_data_field_get_dval_real(gx.mfield, i*set->sp.dx - gx.xoffset, j*set->sp.dy - gx.yoffset, GWY_INTERPOLATION_ROUND) - 1);
                                        if (gx.depth > 0 && mask == 0 && k > val && k < (val + gx.depth))
                                            sv_pool_set_data(mp, gx.mat, set, i, j, k);
                                    }
                                } //z
                            } //y
                        } //x
                    }
                } // for gwydds
            } while (npos < ngtot);

            // finally resolve the mixed submeshes
            if (entity_antialiasing) {
                printf("Resolving aliased cells...\n");
                data_mp_set.mp = mp;
                data_mp_set.set = set;
                g_tree_foreach(submeshes, (GTraverseFunc)resolve_mixed, &data_mp_set);
            }
            g_tree_destroy(submeshes);

            //g_array_free(spheres, TRUE);

            if (mp->nmat != 0) {
                /*alloc structures for recursive accumulator, if necessary*/
                plrc = plrci = dplrc = ade = FALSE;
                for (i = 0; i < mp->nmat; i++) {
                    if (mp->mats[i].type == SV_MAT_DRUDE) {
                        dplrc = TRUE;
                        mp->mats[i].drude_omega_p *= EV_TO_J;// / LIGHT_SPEED;
                        mp->mats[i].drude_nu *= EV_TO_J;// / LIGHT_SPEED;
                        mp->mats[i].epsilon *= EPSILON_0;
                        mp->mats[i].mu *= MU_0;
                     } else if (mp->mats[i].type == SV_MAT_LINTAB) {
                        mp->mats[i].epsilon *= EPSILON_0;
                        mp->mats[i].mu *= MU_0;
                    } else if (mp->mats[i].type == SV_MAT_CP3) {
                        plrc = plrci = TRUE;
                        mp->mats[i].epsilon *= EPSILON_0;
                        mp->mats[i].mu *= MU_0;
                    } else if (mp->mats[i].type == SV_MAT_CP) {
                        plrc = dplrc = plrci = TRUE;
                        mp->mats[i].epsilon *= EPSILON_0;
                        mp->mats[i].mu *= MU_0;
                    } else if (mp->mats[i].type == SV_MAT_ADE) {
                        ade = TRUE;
                        mp->mats[i].epsilon *= EPSILON_0;
                        mp->mats[i].mu *= MU_0;
                    } else if (mp->mats[i].type == SV_MAT_PLRC) {
                        plrc = dplrc = plrci = TRUE;
                        mp->mats[i].epsilon *= EPSILON_0;
                        mp->mats[i].mu *= MU_0;
                    }
                }

                if (dplrc && set->sc.usegpu == 0) {
                    printf("allocating recursive accumulator\n");
                    np = 1;

                    allocate_dplrc(mp->d, np);

                    if (mp->sg && mp->sg->nsubgrid > 0) {
                       for (n=0; n<mp->sg->nsubgrid; n++) allocate_dplrc(mp->sg->sg[n]->d, np);
                    }
                } // if dplrc

                if (plrc && set->sc.usegpu == 0) {
                    printf("allocating recursive accumulator\n");
                    if (plrci)
                        np = 3;
                    else
                        np = 1;

                    allocate_plrc(mp->d, np);

                    if (mp->sg && mp->sg->nsubgrid > 0) {
                       for (n=0; n<mp->sg->nsubgrid; n++) allocate_plrc(mp->sg->sg[n]->d, np);
                    }
                } // if plrc

                if (plrci && set->sc.usegpu == 0) {
                    printf("allocating recursive accumulator imaginary part\n");
                    np = 3;

                    allocate_iplrc(mp->d, np);

                    if (mp->sg && mp->sg->nsubgrid > 0) {
                       for (n=0; n<mp->sg->nsubgrid; n++) allocate_iplrc(mp->sg->sg[n]->d, np);

                    }
                } // if plrci

                if (ade && set->sc.usegpu == 0) {
                    printf("allocating ADE storage\n");
                    allocate_ade(mp->d);

                    if (mp->sg && mp->sg->nsubgrid > 0) {
                       for (n=0; n<mp->sg->nsubgrid; n++) allocate_ade(mp->sg->sg[n]->d);

                    }
                } // if ade	    		
            } // if (mp->nmat != 0)
        } // if (fr == NULL) .... else
    }  // if ((set->sm.in_vector != NULL)

    /*roughness applied on tabulated material if requested*/
    if (set->sm.nroughens) {
        rnd = g_rand_new();

        buffer = (gdouble ***)g_malloc(xres * sizeof(gdouble**));
        for (i = 0; i < xres; i++) {
            buffer[i] = (gdouble **)g_malloc(yres * sizeof(gdouble*));
            for (j = 0; j < yres; j++)
                buffer[i][j] = (gdouble *)g_malloc(zres * sizeof(gdouble));
        }

        for (mi = 0; mi < set->sm.nroughens; mi++) {
            if (set->sc.verbose) {
                if (set->sm.rough_seed[mi] == -1)
                    printf("Roughening material %d with random seed: ", set->sm.rough_matindex[mi]);
                else
                    printf("Roughening material %d with seed %d: ", set->sm.rough_matindex[mi], set->sm.rough_seed[mi]);
                fflush(stdout);
            }

#ifdef G_OS_WIN32
            if (set->sm.rough_seed[mi] == -1)
                g_rand_set_seed(rnd, g_random_int() & 0x7fffffff);
            else
                g_rand_set_seed(rnd, set->sm.rough_seed[mi]);
#else
            if (set->sm.rough_seed[mi] == -1) {
                gettimeofday(&time, NULL);
                g_rand_set_seed(rnd, (time.tv_sec * 1000) + (time.tv_usec / 1000));
            } else g_rand_set_seed(rnd, set->sm.rough_seed[mi]);
#endif

            for (it = 0; it < set->sm.rough_iterations[mi]; it++) {
                for (i = 0; i < xres; i++) {
                    for (j = 0; j < yres; j++) {
                        for (k = 0; k < zres; k++) {
                            if (mp->d->mat->data[i][j][k] == set->sm.rough_matindex[mi])
                                buffer[i][j][k] = 1;
                            else
                                buffer[i][j][k] = -1;
                        }
                    }
                }
                for (i = 0; i < xres; i++) {
                    for (j = 0; j < yres; j++) {
                        for (k = 0; k < zres; k++) {
                            /*add some material if we are inside*/
                            if (mp->d->mat->data[i][j][k] == set->sm.rough_matindex[mi]) {
                                if (g_rand_double_range(rnd, 0, 1) < set->sm.rough_probability[mi]) {
                                    if (closetovoid(mp->d->mat, i, j, k, set->sm.rough_matindex[mi], set->sm.rough_radius_peak[mi] + set->sm.rough_radius_span[mi], xres, yres, zres)) {
                                        add_sphere(buffer, i, j, k, 1, g_rand_double_range(rnd, set->sm.rough_radius_peak[mi] - set->sm.rough_radius_span[mi], set->sm.rough_radius_peak[mi] + set->sm.rough_radius_span[mi]), xres, yres, zres);
                                    }
                                }
                            } else if (mp->d->mat->data[i][j][k] != set->sm.rough_matindex[mi]) {
                                if (g_rand_double_range(rnd, 0, 1) < set->sm.rough_probability[mi]) {
                                    if (closetomat(mp->d->mat, i, j, k, set->sm.rough_matindex[mi], set->sm.rough_radius_peak[mi] + set->sm.rough_radius_span[mi], xres, yres, zres)) {
                                        add_sphere(buffer, i, j, k, -1, g_rand_double_range(rnd, set->sm.rough_radius_peak[mi] - set->sm.rough_radius_span[mi], set->sm.rough_radius_peak[mi] + set->sm.rough_radius_span[mi]), xres, yres, zres);
                                    }
                                }
                            }
                        } // k
                    } // j
                } // i
                for (i = 0; i < xres; i++) {
                    for (j = 0; j < yres; j++) {
                        for (k = 0; k < zres; k++) {
                            if (buffer[i][j][k] > 0)
                                mp->d->mat->data[i][j][k] = set->sm.rough_matindex[mi];
                            else
                                if (mp->d->mat->data[i][j][k] == set->sm.rough_matindex[mi])
                                    mp->d->mat->data[i][j][k] = set->sm.rough_voidindex[mi];
                        }
                    }
                }
                if (set->sc.verbose) {
                    printf(".");
                    fflush(stdout);
                }
            } // it
            if (set->sc.verbose)
                printf("\n");
        } // m

        for (i = 0; i < xres; i++) {
            for (j = 0; j < yres; j++)
                g_free((void **)buffer[i][j]);
            g_free((void *)buffer[i]);
        }
        g_free((void *)buffer);
    } // if (set->sm.nroughens)

    /*roughness applied on tabulated material if requested*/
    if (set->sm.ngrowths) {
        rnd = g_rand_new();
        xres = set->sp.xres;
        yres = set->sp.yres;
        zres = set->sp.zres;

        for (mi = 0; mi < set->sm.ngrowths; mi++) {
            if (set->sm.grow_addindex[mi] >= mp->nmat || set->sm.grow_attachindex[mi] >= mp->nmat) {
                fprintf(stderr, "Error: requesting growth of material above listed material index (in total you have %d materials, want to grow No. %d on No. %d)\n",
                        mp->nmat, set->sm.grow_attachindex[mi], set->sm.grow_addindex[mi]);
                continue;
            }

            sub = set->sm.grow_subsampling[mi];
            gxres = sub*(set->sm.grow_in[mi] - set->sm.grow_i0[mi]);
            gyres = sub*(set->sm.grow_jn[mi] - set->sm.grow_j0[mi]);
            gzres = sub*(set->sm.grow_kn[mi] - set->sm.grow_k0[mi]);

            /*buffer is scaled by subsampling factor*/
            ibuffer = (gint ***)g_malloc(gxres * sizeof(gint**));
            for (i = 0; i < gxres; i++) {
                ibuffer[i] = (gint **)g_malloc(gyres * sizeof(gint*));
                for (j = 0; j < gyres; j++)
                    ibuffer[i][j] = (gint *)g_malloc(gzres * sizeof(gint));
            }
            /*copy scaled data to buffer*/
            for (i = 0; i < gxres; i++) {
                for (j = 0; j < gyres; j++) {
                    for (k = 0; k < gzres; k++) {

                        if ((i / sub + set->sm.grow_i0[mi]) >= 0 && (i / sub + set->sm.grow_i0[mi]) < xres &&
                            (j / sub + set->sm.grow_j0[mi]) >= 0 && (j / sub + set->sm.grow_j0[mi]) < yres &&
                            (k / sub + set->sm.grow_k0[mi]) >= 0 && (k / sub + set->sm.grow_k0[mi]) < zres &&
                            mp->d->mat->data[i / sub + set->sm.grow_i0[mi]][j / sub + set->sm.grow_j0[mi]][k / sub + set->sm.grow_k0[mi]] == set->sm.grow_addindex[mi])
                            ibuffer[i][j][k] = set->sm.grow_addindex[mi];
                        else
                            ibuffer[i][j][k] = -1;
                    }
                }
            }

            if (set->sc.verbose) {
                if (set->sm.grow_seed[mi] == -1)
                    printf("Growing material %d on %d with random seed: ", set->sm.grow_attachindex[mi], set->sm.grow_addindex[mi]);
                else
                    printf("Growing material %d on %d with seed %d: ", set->sm.grow_attachindex[mi], set->sm.grow_addindex[mi], set->sm.grow_seed[mi]);
                fflush(stdout);
            }

#ifdef G_OS_WIN32
            if (set->sm.grow_seed[mi] == -1)
                g_rand_set_seed(rnd, g_random_int() & 0x7fffffff);
            else
                g_rand_set_seed(rnd, set->sm.grow_seed[mi]);
#else
            if (set->sm.grow_seed[mi] == -1) {
                gettimeofday(&time, NULL);
                g_rand_set_seed(rnd, (time.tv_sec * 1000) + (time.tv_usec / 1000));
            } else
                g_rand_set_seed(rnd, set->sm.grow_seed[mi]);
#endif

            maxdd = sqrt(gxres*gxres + gyres*gyres + gzres*gzres);
            nbdist = (gint)set->sm.grow_mobility[mi];
            probability = set->sm.grow_probability[mi]; //fixed value 20

            it = 0;

            do {

#pragma omp parallel default(shared) private(posx, posy, posz, face, theta, phi, dx, dy, dz, dd, nbscore, iposx, iposy, iposz, nbbest, ni, nj, nk)
#pragma omp for nowait
                for (jt = 0; jt < 100; jt++) {
                    land[jt] = 0;
                    maxdd = sqrt(gxres*gxres + gyres*gyres + gzres*gzres);
                    nbdist = (gint)set->sm.grow_mobility[mi];
                    probability = set->sm.grow_probability[mi]; //fixed value 20

                    /*select random position of start, choose face and then position it (not yet)*/
                    face = g_rand_int_range(rnd, 0, 6);

                    if (face == 0 && set->sm.grow_skipi0[mi] == 0) {
                        posx = nbdist + 1;
                        posy = g_rand_int_range(rnd, 1, gyres - 1);
                        posz = g_rand_int_range(rnd, 1, gzres - 1);

                    } else if (face == 1 && set->sm.grow_skipin[mi] == 0) {
                        posx = gxres - nbdist - 2;
                        posy = g_rand_int_range(rnd, 1, gyres - 1);
                        posz = g_rand_int_range(rnd, 1, gzres - 1);
                    } else if (face == 2 && set->sm.grow_skipj0[mi] == 0) {
                        posy = nbdist + 1;
                        posx = g_rand_int_range(rnd, 1, gxres - 1);
                        posz = g_rand_int_range(rnd, 1, gzres - 1);
                    } else if (face == 3 && set->sm.grow_skipjn[mi] == 0) {
                        posy = gyres - nbdist - 2;
                        posx = g_rand_int_range(rnd, 1, gxres - 1);
                        posz = g_rand_int_range(rnd, 1, gzres - 1);
                    } else if (face == 4 && set->sm.grow_skipk0[mi] == 0) {
                        posz = nbdist + 1;
                        posx = g_rand_int_range(rnd, 1, gxres - 1);
                        posy = g_rand_int_range(rnd, 1, gyres - 1);
                    } else if (face == 5 && set->sm.grow_skipkn[mi] == 0) {
                        posz = gzres - nbdist - 2;
                        posx = g_rand_int_range(rnd, 1, gxres - 1);
                        posy = g_rand_int_range(rnd, 1, gyres - 1);
                    } else
                        continue;

                    /*select random flight direction*/
                    //dx = g_rand_double_range(rnd, -1, 1);
                    //dy = g_rand_double_range(rnd, -1, 1);
                    //dz = g_rand_double_range(rnd, -1, 1);

                    //theta = g_rand_double_range(rnd, 0, G_PI);
                    theta = acos(2.0*g_rand_double_range(rnd, 0, 1) - 1);
                    phi = g_rand_double_range(rnd, 0, 2 * G_PI);

                    dx = sin(theta)*cos(phi);
                    dy = sin(theta)*sin(phi);
                    dz = cos(theta);

                    //printf("0 0 0 %g %g %g\n", dx, dy, dz);
                    //printf("\n\n%g %g %g\n", posx, posy, posz);

                    /*flight*/
                    for (dd = 0; dd < maxdd; dd++) {
                        posx += dx;
                        posy += dy;
                        posz += dz;

                        //if (posx < -100 || posy < -100 || posz < -100 || posx >= (gxres + 100) || posy >= (gyres + 100) || posz >= (gzres + 100))
                        //    break;

                        if (posx < 1 || posy < 1 || posz < 1 || posx >= (gxres - 1) || posy >= (gyres - 1) || posz >= (gzres - 1))
                            continue;

                        if (set->sm.grow_skipin[mi] != 0) {
                            if (posx >= (gxres - 10))
                                break;
                        }

                        iposx = (gint)posx;
                        iposy = (gint)posy;
                        iposz = (gint)posz;

                        if (ibuffer[iposx][iposy][iposz] == set->sm.grow_addindex[mi] || ibuffer[iposx][iposy][iposz] == set->sm.grow_attachindex[mi])
                            break;

                       //printf("tak co bude %d %d %d   %d\n", iposx, iposy, iposz, ibuffer[iposx][iposy][iposz]);
                        //printf("%d %d %d\n", iposx, iposy, iposz);

                        if (ibuffer[iposx + 1][iposy][iposz] == set->sm.grow_addindex[mi] || ibuffer[iposx + 1][iposy][iposz] == set->sm.grow_attachindex[mi] ||
                            ibuffer[iposx - 1][iposy][iposz] == set->sm.grow_addindex[mi] || ibuffer[iposx - 1][iposy][iposz] == set->sm.grow_attachindex[mi] ||
                            ibuffer[iposx][iposy + 1][iposz] == set->sm.grow_addindex[mi] || ibuffer[iposx][iposy + 1][iposz] == set->sm.grow_attachindex[mi] ||
                            ibuffer[iposx][iposy - 1][iposz] == set->sm.grow_addindex[mi] || ibuffer[iposx][iposy - 1][iposz] == set->sm.grow_attachindex[mi] ||
                            ibuffer[iposx][iposy][iposz + 1] == set->sm.grow_addindex[mi] || ibuffer[iposx][iposy][iposz + 1] == set->sm.grow_attachindex[mi] ||
                            ibuffer[iposx][iposy][iposz - 1] == set->sm.grow_addindex[mi] || ibuffer[iposx][iposy][iposz - 1] == set->sm.grow_attachindex[mi]) {

                            /*relax if necessary*/

                            if (nbdist > 0) {
                                landx[jt] = iposx;
                                landy[jt] = iposy;
                                landz[jt] = iposz;
                                land[jt] = 1;
                                if (g_rand_double_range(rnd, 0, 1) < probability) {
                                    nbscore = nbbest = evalnbs(ibuffer, iposx, iposy, iposz, set->sm.grow_addindex[mi], set->sm.grow_attachindex[mi], 27);
                                    if (iposx > nbdist && iposy > nbdist && iposz > nbdist && iposx < (gxres - nbdist - 1) && iposy < (gyres - nbdist - 1) && iposz < (gzres - nbdist - 1)) {
                                        for (ni = (iposx - nbdist); ni < (iposx + nbdist); ni++) {
                                            for (nj = (iposy - nbdist); nj < (iposy + nbdist); nj++) {
                                                for (nk = (iposz - nbdist); nk < (iposz + nbdist); nk++) {
                                                    //evaluate only points at surface. i.e. not occupied but in touch with material
                                                    if ((ibuffer[ni][nj][nk] != set->sm.grow_addindex[mi] && ibuffer[ni][nj][nk] != set->sm.grow_attachindex[mi]) &&
                                                        (ibuffer[ni + 1][nj][nk] == set->sm.grow_addindex[mi] || ibuffer[ni + 1][nj][nk] == set->sm.grow_attachindex[mi] ||
                                                         ibuffer[ni - 1][nj][nk] == set->sm.grow_addindex[mi] || ibuffer[ni - 1][nj][nk] == set->sm.grow_attachindex[mi] ||
                                                         ibuffer[ni][nj + 1][nk] == set->sm.grow_addindex[mi] || ibuffer[ni][nj + 1][nk] == set->sm.grow_attachindex[mi] ||
                                                         ibuffer[ni][nj - 1][nk] == set->sm.grow_addindex[mi] || ibuffer[ni][nj - 1][nk] == set->sm.grow_attachindex[mi] ||
                                                         ibuffer[ni][nj][nk + 1] == set->sm.grow_addindex[mi] || ibuffer[ni][nj][nk + 1] == set->sm.grow_attachindex[mi] ||
                                                         ibuffer[ni][nj][nk - 1] == set->sm.grow_addindex[mi] || ibuffer[ni][nj][nk - 1] == set->sm.grow_attachindex[mi])) {
                                                        nbscore = evalnbs(ibuffer, ni, nj, nk, set->sm.grow_addindex[mi], set->sm.grow_attachindex[mi], 27 - nbbest);
                                                        if (nbscore > nbbest) {
                                                            nbbest = nbscore;
                                                            landx[jt] = ni;
                                                            landy[jt] = nj;
                                                            landz[jt] = nk;
                                                        }
                                                    }
                                                } // nk
                                            } // nj
                                        } // ni
                                    } // if
                                } // if        
                            } else { // purely ballistic deposition
                                landx[jt] = iposx;
                                landy[jt] = iposy;
                                landz[jt] = iposz;
                                land[jt] = 1;
                            }

                            break;
                        } // if ibuffer...
                    } // for dd
                } // for jt
                /*apply data*/
                if (set->sc.verbose && (it % 100000) == 0) {
                    printf(".");
                    fflush(stdout);
                }

                for (jt = 0; jt < 100; jt++) {
                    if (land[jt] == 1) {
                        // printf("yes at %d %d %d\n", landx[jt], landy[jt], landz[jt]);
                        ibuffer[landx[jt]][landy[jt]][landz[jt]] = set->sm.grow_attachindex[mi];
                    }
                }
                it += 100;

                /*debugging data output*/
                /*
                            if ((it%5000000) == 0) {
                                gchar *out = g_strdup_printf("xgrowth_%s_%04dM.vtk", set->so.cubs[0].filebase, (int)(it/1000000));

                                for (i = 0; i < gxres; i++) {
                                    for (j = 0; j < gyres; j++) {
                                        for (k = 0; k < gzres; k++) {
                                            if ((i/sub + set->sm.grow_i0[m])>=0 && (i/sub + set->sm.grow_i0[m])<xres &&
                                                (j/sub + set->sm.grow_j0[m])>=0 && (j/sub + set->sm.grow_j0[m])<yres &&
                                                (k/sub + set->sm.grow_k0[m])>=0 && (k/sub + set->sm.grow_k0[m])<zres &&
                                                ibuffer[i][j][k] == set->sm.grow_addindex[m]) {
                                                mp->d->mat->data[i/sub + set->sm.grow_i0[m]][j/sub + set->sm.grow_j0[m]][k/sub + set->sm.grow_k0[m]] = set->sm.grow_addindex[m];
                                            }
                                        } // k
                                    } // j
                                } // i

                                output_ivtk(mp->d->mat->data, set->sp.xres, set->sp.yres, set->sp.zres, out, "material");
                            }
                */
                /*end of debugging data output*/
            } while (it < set->sm.grow_nsteps[mi]);

            /*copy scaled data back from ibuffer*/
            for (i = 0; i < gxres; i++) {
                for (j = 0; j < gyres; j++) {
                    for (k = 0; k < gzres; k++) {
                        if ((i / sub + set->sm.grow_i0[mi]) >= 0 && (i / sub + set->sm.grow_i0[mi]) < xres &&
                            (j / sub + set->sm.grow_j0[mi]) >= 0 && (j / sub + set->sm.grow_j0[mi]) < yres &&
                            (k / sub + set->sm.grow_k0[mi]) >= 0 && (k / sub + set->sm.grow_k0[mi]) < zres &&
                            ibuffer[i][j][k] == set->sm.grow_attachindex[mi]) {
                            mp->d->mat->data[i / sub + set->sm.grow_i0[mi]][j / sub + set->sm.grow_j0[mi]][k / sub + set->sm.grow_k0[mi]] = set->sm.grow_attachindex[mi];
                        }
                    } // k
                } // j
            } // i
            /*free ibuffer*/
            for (i = 0; i < gxres; i++) {
                for (j = 0; j < gyres; j++)
                    g_free((void **)ibuffer[i][j]);
                g_free((void *)ibuffer[i]);
            }
            g_free((void *)ibuffer);
        } // for m
        if (set->sc.verbose)
            printf("\n");
    } // is set->ngrowths



    /*expression applied*/
    if (set->sm.nexprs) {
        buffer = (gdouble ***)g_malloc(xres * sizeof(gdouble**));
        for (i = 0; i < xres; i++) {
            buffer[i] = (gdouble **)g_malloc(yres * sizeof(gdouble*));
            for (j = 0; j < yres; j++)
                buffer[i][j] = (gdouble *)g_malloc(zres * sizeof(gdouble));
        }

        for (mi = 0; mi < set->sm.nexprs; mi++) {
            if (set->sc.verbose) {
                printf("Applying modifier expression %s to material %d: ", set->sm.expr_expr[mi], set->sm.expr_matindex[mi]);
                fflush(stdout);
            }

            expr = gwy_expr_new();
            gwy_expr_define_constant(expr, "pi", G_PI, NULL);
            if (!gwy_expr_compile(expr, set->sm.expr_expr[mi], NULL))
                fprintf(stderr, "Error: cannot compile expression for surface modifier\n");
            else {
                if (gwy_expr_resolve_variables(expr, G_N_ELEMENTS(var_names), var_names, var_positions)) {
                    fprintf(stderr, "Error: cannot evaluate expression for surface modifier: strange variables\n");
                } else {
                    printf("Variable order is %d %d %d\n", var_positions[0], var_positions[1], var_positions[2]);

                    for (i = 0; i < xres; i++) {
                        for (j = 0; j < yres; j++) {
                            for (k = 0; k < zres; k++) {
                                if (mp->d->mat->data[i][j][k] == set->sm.expr_matindex[mi])
                                    buffer[i][j][k] = 1;
                                else
                                    buffer[i][j][k] = 0;
                            } // k
                        } // j
                    } // i
                    for (i = 1; i < (xres - 1); i++) {
                        for (j = 1; j < (yres - 1); j++) {
                            for (k = 1; k < (zres - 1); k++) {
                                if (mp->d->mat->data[i][j][k] == set->sm.expr_matindex[mi]
                                    && (mp->d->mat->data[i + 1][j][k] != set->sm.expr_matindex[mi] || mp->d->mat->data[i - 1][j][k] != set->sm.expr_matindex[mi]
                                        || mp->d->mat->data[i][j + 1][k] != set->sm.expr_matindex[mi] || mp->d->mat->data[i][j - 1][k] != set->sm.expr_matindex[mi]
                                        || mp->d->mat->data[i][j][k + 1] != set->sm.expr_matindex[mi] || mp->d->mat->data[i][j][k - 1] != set->sm.expr_matindex[mi])) //do this only at boundary
                                    addmat(buffer, i, j, k, (gdouble)set->sm.expr_maxdist[mi], xres, yres, zres);
                            } // k
                        } // j
                        if (set->sc.verbose && (i % 10) == 0) {
                            printf(".");
                            fflush(stdout);
                        }
                    } // i
                    printf(":"); fflush(stdout);
                    for (i = 0; i < xres; i++) {
                        for (j = 0; j < yres; j++) {
                            for (k = 0; k < zres; k++) {
                                /*run the expression*/
                                vars[var_positions[0]] = (gdouble)i;
                                vars[var_positions[1]] = (gdouble)j;
                                vars[var_positions[2]] = (gdouble)k;
                                buffer[i][j][k] += gwy_expr_execute(expr, vars);
                            } // k
                        } // j
                        if (set->sc.verbose && (i % 10) == 0) {
                            printf(".");
                            fflush(stdout);
                        }
                    } // i
                    for (i = MAX(0, set->sm.expr_i0[mi]); i < MIN(xres, set->sm.expr_in[mi]); i++) {
                        for (j = MAX(0, set->sm.expr_j0[mi]); j < MIN(yres, set->sm.expr_jn[mi]); j++) {
                            for (k = MAX(0, set->sm.expr_k0[mi]); k < MIN(zres, set->sm.expr_kn[mi]); k++) {
                                if (buffer[i][j][k] >= 1)
                                    mp->d->mat->data[i][j][k] = set->sm.expr_matindex[mi];
                                else
                                    if (mp->d->mat->data[i][j][k] == set->sm.expr_matindex[mi])
                                        mp->d->mat->data[i][j][k] = set->sm.expr_voidindex[mi];
                            } // k
                        } // j
                    } // i
                } // else
            } // else

            gwy_expr_free(expr);

            if (set->sc.verbose)
                printf("Done\n");
        }

        for (i = 0; i < xres; i++) {
            for (j = 0; j < yres; j++)
                g_free((void **)buffer[i][j]);
            g_free((void *)buffer[i]);
        }
        g_free((void *)buffer);
    }

    /*soft blur if requested, use only one additional field for that*/
    /* SHOULD HAVE BEEN DISABLED (RS)
    if (set->sm.smooth) {
        buffer = (gdouble ***) g_malloc(xres*sizeof(gdouble**));
        for (i = 0; i < xres; i++) {
            buffer[i]= (gdouble **) g_malloc(yres*sizeof(gdouble*));
            for (j = 0; j < yres; j++)
                buffer[i][j]= (gdouble *) g_malloc(zres*sizeof(gdouble));
        }

        for (m = 0; m < set->sm.smooth; m++) {
            for (i = 1; i < (xres-1); i++) { //x
                for (j = 1; j < (yres-1); j++) {//y
                    for (k = 1; k < (zres-1); k++) {//z
                        buffer[i][j][k] = mp->d->epsilon->data[i][j][k]/2
                            + mp->d->epsilon->data[i-1][j][k]/12 + mp->d->epsilon->data[i+1][j][k]/12
                            + mp->d->epsilon->data[i][j-1][k]/12 + mp->d->epsilon->data[i][j+1][k]/12
                            + mp->d->epsilon->data[i][j][k-1]/12 + mp->d->epsilon->data[i][j][k+1]/12;
                    }
                }
            }
            for (i = 1; i < (xres-1); i++) { //x
                for (j = 1; j < (yres-1); j++) { //y
                    for (k = 1; k < (zres-1); k++) //z
                        mp->d->epsilon->data[i][j][k] = buffer[i][j][k];
                }
            }
        } // m

        for (i=0; i<xres; i++) {
            for (j=0; j<yres; j++)
                g_free((void **) buffer[i][j]);
            g_free((void *) buffer[i]);
        }
    g_free((void *)buffer);
    }
    */

    if (set->sm.crop) {
        xres = set->sp.xres;
        yres = set->sp.yres;
        zres = set->sp.zres;

        domat = 1;
        for (i = 0; i < xres; i++) {//x
            for (j = 0; j < yres; j++) { //y
                for (k = 0; k < zres; k++) { //z
                    if (i <= set->sm.crop_i0 || i >= set->sm.crop_in || j <= set->sm.crop_j0 || j >= set->sm.crop_jn || k <= set->sm.crop_k0 || k >= set->sm.crop_kn) {
                        // TODO - check whether sv_pool_set_data can be used
                        // mat.type = SV_MAT_LINEAR; mat.epsilon = set->sm.crop_epsilon; mat.sigma = set->sm.crop_sigma; mat.mu = set->sm.crop_mu; mat.sigast = set->sm.crop_sigast;
                        // sv_pool_set_data(mp, set, mat)
                        if (mp->nmat == 0) { //voxel-by voxel material
                            if (mp->set->plan.matmode == SV_MATMODE_FULL || mp->set->plan.matmode == SV_MATMODE_ELECTRIC) {
                                mp->d->epsilon->data[i][j][k] = (gfloat)(set->sm.crop_epsilon*EPSILON_0);
                                mp->d->sigma->data[i][j][k] = (gfloat)(set->sm.crop_sigma);
                            }
                            if (mp->set->plan.matmode == SV_MATMODE_FULL || mp->set->plan.matmode == SV_MATMODE_MAGNETIC) {
                                mp->d->mu->data[i][j][k] = (gfloat)(set->sm.crop_mu*MU_0);
                                mp->d->sigast->data[i][j][k] = (gfloat)(set->sm.crop_sigast);
                            }
                        } else { //tabulated material
                            if (domat) {
                                mat.type = SV_MAT_LINTAB;
                                mat.epsilon = set->sm.crop_epsilon;
                                mat.mu = set->sm.crop_mu;
                                mat.sigma = 0;
                                mat.sigast = 0;

                                printf("Crop: adding material for cropping medium.\n");
                                if (sv_pool_get_matindex(mp, mat) == -1)
                                    sv_pool_add_material(mp, &mat, TRUE);
                                domat = 0;
                            }
                            mp->d->mat->data[i][j][k] = mat.pos;
                        }
                    } // if
                } // k
            } // j
        } // i
    } // if

    if (set->sm.nspectrals) 
    {
        for (mi = 0; mi < set->sm.nspectrals; mi++) {

            modifier_fftshift(mp->d->mat, set->sm.spectral_matindex[mi], set->sm.spectral_sigma[mi], set->sm.spectral_t[mi], 
                           set->sm.spectral_sigma[mi], set->sm.spectral_t[mi], set->sm.spectral_sigma[mi], set->sm.spectral_t[mi], set->sm.spectral_seed[mi]);

//            mp->mats = modifier_fftshift_antialiased(mp->d->mat, mp->mats, &(mp->nmat), set->sm.spectral_matindex[mi], set->sm.spectral_sigma[mi], set->sm.spectral_t[mi], 
//                           set->sm.spectral_sigma[mi], set->sm.spectral_t[mi], set->sm.spectral_sigma[mi], set->sm.spectral_t[mi], set->sm.spectral_seed[mi], 2);
         }
    }

/*create source skipping field if requested, working for TSF source*/
    if (set->ss.tsf.box_boundary_skipdepth >= 0) {
        if (set->sc.verbose) {
            printf("Precomputing TSF source skip...");
            fflush(stdout);
        }


        mp->d->sourceskip = sv_icube_new(set->sp.xres, set->sp.yres, set->sp.zres, set->sp.xres*set->sp.dx, set->sp.yres*set->sp.dy, set->sp.zres*set->sp.dz, 1);
        sv_icube_fill(mp->d->sourceskip, 0); // init with vacuum

        for (i = 0; i < set->sp.xres; i++) {
            for (j = 0; j < set->sp.yres; j++) {
                for (k = 0; k < set->sp.zres; k++) {
                    /*do not evaluate internal part*/
                    if (i > (mp->src->tsf->box_i0 + set->ss.tsf.box_boundary_skipdepth)
                        && j > (mp->src->tsf->box_j0 + set->ss.tsf.box_boundary_skipdepth)
                        && k > (mp->src->tsf->box_k0 + set->ss.tsf.box_boundary_skipdepth)
                        && i < (mp->src->tsf->box_in - set->ss.tsf.box_boundary_skipdepth)
                        && j < (mp->src->tsf->box_jn - set->ss.tsf.box_boundary_skipdepth)
                        && k < (mp->src->tsf->box_kn - set->ss.tsf.box_boundary_skipdepth))
                        continue;

                    if (eval_closecheck(mp, set, i, j, k, set->ss.tsf.box_boundary_skipdepth))
                        mp->d->sourceskip->data[i][j][k] = 1;
                }
            }
        }
        if (set->sc.verbose)
            printf("Done.\n");
    }
    if (set->ss.ltsf.box_boundary_skipdepth >= 0) {
        if (set->sc.verbose) {
            printf("Precomputing LTSF source skip...");
            fflush(stdout);
        }


        mp->d->sourceskip = sv_icube_new(set->sp.xres, set->sp.yres, set->sp.zres, set->sp.xres*set->sp.dx, set->sp.yres*set->sp.dy, set->sp.zres*set->sp.dz, 1);
        sv_icube_fill(mp->d->sourceskip, 0); // init with vacuum

        for (i = 0; i < set->sp.xres; i++) {
            for (j = 0; j < set->sp.yres; j++) {
                for (k = 0; k < set->sp.zres; k++) {
                    /*do not evaluate internal part*/
                    if (i > (mp->src->ltsf->box_i0 + set->ss.ltsf.box_boundary_skipdepth)
                        && j > (mp->src->ltsf->box_j0 + set->ss.ltsf.box_boundary_skipdepth)
                        && k > (mp->src->ltsf->box_k0 + set->ss.ltsf.box_boundary_skipdepth)
                        && i < (mp->src->ltsf->box_in - set->ss.ltsf.box_boundary_skipdepth)
                        && j < (mp->src->ltsf->box_jn - set->ss.ltsf.box_boundary_skipdepth)
                        && k < (mp->src->ltsf->box_kn - set->ss.ltsf.box_boundary_skipdepth))
                        continue;

                    if (eval_closecheck(mp, set, i, j, k, set->ss.ltsf.box_boundary_skipdepth))
                        mp->d->sourceskip->data[i][j][k] = 1;
                }
            }
        }
        if (set->sc.verbose)
            printf("Done.\n");
    }

    //evaluate volume, for debug purposes only

    /*
    mattosum = 1;
    matsum = 0;

    for (i = 0; i < xres; i++) {//x
        for (j = 0; j < yres; j++) { //y
            for (k = 0; k < zres; k++) { //z
                 if (mp->d->mat->data[i][j][k] == mattosum) matsum += 1;
            }
        }
    }
    printf("volume: %d voxels %g m^3  effective radius %g nm\n", matsum, matsum*set->sp.dx*set->sp.dy*set->sp.dz, pow(3.0/4.0/G_PI*((gdouble)matsum*set->sp.dx*set->sp.dy*set->sp.dz), 1/3.0));
    */
    //end of evaluate volume




    /*fill precomputed internal boundaries*/
    /*there are different possible boundary variants, but now treat everything the same way and base it on epsilon and mat only*/

    if (set->sc.verbose) {
        printf("Precomputing internal boundaries...");
        fflush(stdout);
    }

    sv_icube_fill(mp->d->bnds, 0); // init with no boundaries

    for (i = 0; i < (set->sp.xres); i++) {
        for (j = 0; j < (set->sp.yres); j++) {
            for (k = 0; k < (set->sp.zres); k++) {
                if (i == 0 || j == 0 || k == 0 || i == (set->sp.xres - 1) || j == (set->sp.yres - 1) || k == (set->sp.zres - 1))
                    mp->d->bnds->data[i][j][k] = 1;
                else {
                    matindex = sv_yee_data_get_epsilon_sigma(mp->d, set, mp->mats, mp->nmat, i, j, k, &epsilon, &sigma);

                    for (ii = -1; ii <= 1; ii++) {
                        for (jj = -1; jj <= 1; jj++) {
                            for (kk = -1; kk <= 1; kk++) {
                                if (mp->d->bnds->data[i][j][k] == 1)
                                    continue;

                                nmatindex = sv_yee_data_get_epsilon_sigma(mp->d, set, mp->mats, mp->nmat, i + ii, j + jj, k + kk, &nepsilon, &nsigma);
                                if (nmatindex != matindex || nepsilon != epsilon || nsigma != sigma) {
                                    mp->d->bnds->data[i][j][k] = 1;
                                }
                            }
                        }
                    }
                }
                //              printf("bnd %d %d %d %d\n", i, j, k, mp->d->bnds->data[i][j][k]);

            }
        }
    }
    if (set->sc.verbose)
        printf("Done.\n");




    //printf("%d %d\n", mp->mats[1].type, mp->mats[2].type);

    /*gint m1=0, m2=0, m3=0;
    for (i=0; i<xres; i++) { //x
        for (j=0; j<yres; j++) { //y
            for (k=0; k<zres; k++) { //z
                if (mp->d->mat->data[i][j][k]==1) m1++;
                if (mp->d->mat->data[i][j][k]==2) m2++;
                if (mp->d->mat->data[i][j][k]==3) m3++;
            }
        }
    }

    printf("sphere outputs........:   %d %d %d\n", m1, m2, m3);*/

    return 0;
}

SvPool* init_and_make_plan(SvSet *set)
{
    SvPool *mp;
    gdouble alloc, ram, shared;
    gchar matfile[256];
    gint testout = SV_MATMODE_FULL, i;
    gint matmode_mat, matmode_script;

    set->plan.gpumode = SV_GPUMODE_NONE;

    set->plan.dt = 1.0 / LIGHT_SPEED / sqrt(1.0 / set->sp.dx / set->sp.dx + 1.0 / set->sp.dy / set->sp.dy + 1.0 / set->sp.dz / set->sp.dz) * set->sc.dtmult;
    mp = sv_pool_new(set);

    /*load sources and allocate what is necessary*/
    if (set->sc.verbose > 1)
        printf("Loading sources...\n");
    load_sources(mp, set);

    /*check if we are expected to use a database optical properties. If yes, write all the
     * vector material properties to a new file, already using the "right" material types and
     * parameters */
    if (set->sc.verbose > 1)
        printf("Testing material files...\n");
    if (set->sm.in_vector_filename != NULL && (testout = test_script_file(set->sm.in_vector_filename)) == SV_MATMODE_DATABASE) {
        if (override_matfile_from_database(set->sm.in_vector_filename, set->ss.lambda_min, set->ss.lambda_center, set->ss.lambda_max, set->sc.suffix, 0, NULL, NULL, set->sm.localfiles))
            return NULL;
        g_snprintf(matfile, sizeof(matfile), "tmp_matfile_%06d", set->sc.suffix);
        set->sm.in_vector_filename = g_strdup(matfile);
    }

    if (testout == SV_MATMODE_ERROR)
        fprintf(stderr, "Error loading material file\n");

    /*change string based summing output to real one*/
    for (i = 0; i < set->so.nsums; i++) {
        if (set->so.sums[i].stringbased) {
            if (subst_mat(set->so.sums[i].string,
                          &(set->so.sums[i].layered_epsilon), &(set->so.sums[i].layered_mu),
                          &(set->so.sums[i].layered_sigma), &(set->so.sums[i].layered_sigast),
                          set->ss.lambda_min, set->ss.lambda_center, set->ss.lambda_max, set->sm.localfiles))
                return NULL;
            //  printf("ssssssssss %g %g\n", set->so.sums[i].epsilon, set->so.sums[i].sigma);
            set->so.sums[i].layered_epsilon *= EPSILON_0;
            set->so.sums[i].layered_mu *= MU_0;
        }
    }

    /*determine matmode - whether to alloc all material fields*/
    /*first test user options - none at present*/
    /*then read material and script file if any and check for materials*/
    if (set->sm.matmode_check) { //if not, material was set by user already
        if (set->sc.verbose > 1)
            printf("Checking material mode...\n");
        matmode_mat = SV_MATMODE_NONE;
        matmode_script = SV_MATMODE_NONE;
        //printf("matmodes: %d %d\n", matmode_mat, matmode_script);

        if (set->sm.in_voxel_filename != NULL)
            matmode_mat = test_material_file(set->sm.in_voxel_filename);
        if (set->sm.in_vector_filename != NULL)
            matmode_script = test_script_file(set->sm.in_vector_filename);
        //printf("matmodes: %d %d\n", matmode_mat, matmode_script);

        if (matmode_mat == SV_MATMODE_NONE && matmode_script == SV_MATMODE_NONE)
            set->plan.matmode = SV_MATMODE_NONE;
        else if (matmode_mat == SV_MATMODE_ELECTRIC && (matmode_script == SV_MATMODE_NONE || matmode_script == SV_MATMODE_ELECTRIC))
            set->plan.matmode = SV_MATMODE_ELECTRIC;
        else if (matmode_mat == SV_MATMODE_MAGNETIC && (matmode_script == SV_MATMODE_NONE || matmode_script == SV_MATMODE_MAGNETIC))
            set->plan.matmode = SV_MATMODE_MAGNETIC;
        else if (matmode_script == SV_MATMODE_ELECTRIC && (matmode_mat == SV_MATMODE_NONE || matmode_mat == SV_MATMODE_ELECTRIC))
            set->plan.matmode = SV_MATMODE_ELECTRIC;
        else if (matmode_script == SV_MATMODE_MAGNETIC && (matmode_mat == SV_MATMODE_NONE || matmode_mat == SV_MATMODE_MAGNETIC))
            set->plan.matmode = SV_MATMODE_MAGNETIC;
        else
            set->plan.matmode = SV_MATMODE_FULL;
    } else
        set->plan.matmode = SV_MATMODE_FULL;

    /*assume what memory will be allocated*/
    alloc = set->sp.xres * set->sp.yres * set->sp.zres;
    if (set->plan.matmode == SV_MATMODE_FULL)
        ram = alloc * 10 * sizeof(gdouble);
    else if (set->plan.matmode == SV_MATMODE_ELECTRIC || set->plan.matmode == SV_MATMODE_MAGNETIC)
        ram = alloc * 8 * sizeof(gdouble);
    else
        ram = alloc * 6 * sizeof(gdouble);
    if (set->sc.verbose > 1)
        printf("RAM allocated memory: %g MB\n", ram / 1e6);

    if (!set->plan.gpumode == SV_GPUMODE_NONE) {
        if (set->plan.matmode == SV_MATMODE_FULL)
            shared = alloc * 10 * sizeof(float);
        else if (set->plan.matmode == SV_MATMODE_ELECTRIC || set->plan.matmode == SV_MATMODE_MAGNETIC)
            shared = alloc * 8 * sizeof(float);
        else
            shared = alloc * 6 * sizeof(float);
        if (set->sc.verbose > 1)
            printf("GPU allocated memory: %g MB\n", shared / 1e6);
    }

    /*alloc main fields*/

    if (set->sc.verbose > 1)
        printf("Initializing main fields...\n");
    sv_pool_allocate(set, mp);

    /*load sources and allocate what is necessary*/
    //if (set->sc.verbose>1) printf("Loading sources...\n");
    //load_sources(mp, set);

   /*determine how many threads will be used if auto mode is requested*/
    if (set->sc.nthreads == -1)
        set->sc.nthreads = omp_get_max_threads(); //omp does what it wants now, ignoring this and running all the available threads
    else
        omp_set_num_threads(set->sc.nthreads);

    /*load main fields from data files and allocate what is necessary*/
    if (set->sm.in_voxel_filename != NULL || set->sm.in_vector_filename != NULL) {
        if (set->sc.verbose > 1)
            printf("Loading material data...\n");
        load_data(mp, set);
    }

    /*allocate GPUs depending on what is in sources and data*/
    if (set->sc.verbose > 1)
        printf("Allocate GPU fields...\n");
    sv_pool_allocate_gpus(mp, set);

    /*allocate CPU outputs (gpu output data allocated in previous step)*/
    if (set->sc.verbose > 1)
        printf("Allocate outputs...\n");
    sv_pool_allocate_output(mp, set);

    if (set->plan.gpumode == SV_GPUMODE_NONE) {
        if (set->sc.verbose > 1)
            printf("Allocate farfield storage...\n");
        sv_pool_farfield_allocate_storage(mp, set);
    }

    if (set->sc.verbose > 1) {
        printf("Computation plan:\n");
        if (set->plan.matmode == SV_MATMODE_FULL)
            printf("Material properties: full\n");
        else if (set->plan.matmode == SV_MATMODE_ELECTRIC)
            printf("Material properties: electric material only\n");
        else if (set->plan.matmode == SV_MATMODE_MAGNETIC)
            printf("Material properties: mangetic material only\n");
        else if (set->plan.matmode == SV_MATMODE_NONE && mp->nmat != 0)
            printf("Material properties: tabulated\n");
        else
            printf("Material properties: none\n");

        printf("System has %d cores, up to %d thread(s) will be used for computation on CPU\n", get_n_cores(), set->sc.nthreads);

        if (set->plan.gpumode == SV_GPUMODE_NONE)
            printf("GPU use: none\n");
        else
            printf("GPU use: full computation\n");
        printf("Timestep: %g s\n", set->plan.dt);
    }

    return mp;
}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
