
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

 /*  output.c :
  *  all the output functions
  */

#include <stdlib.h>
#include "source.h"
#include "pool.h"
#include "constants.h"
#include "settings.h"
#include <math.h>
#include <libgwymodule/gwymodule-file.h>
#include <libgwydgets/gwydgets.h>
  //#include <app/gwyapp.h>

#define MAGIC "GWYO"
#define MAGIC2 "GWYP"
#define MAGIC_SIZE (sizeof(MAGIC)-1)

SvOutput*
sv_output_new()
{
    SvOutput* out = (SvOutput*)g_malloc(sizeof(SvOutput));

    out->data = gwy_container_new();
    out->last = 0;
    out->glast = 0;

    return out;
}

int
is_output(SvPool *mp, SvSet *set)
{
    gint i;

    for (i = 0; i < set->so.nimgs; i++) {
        if (set->so.imgs[i].step == 0)
            return 0;
        if (set->sc.step_act % set->so.imgs[i].step == 0)
            return 1;
    }

    for (i = 0; i < set->so.npnts; i++) {
        if (set->so.pnts[i].step == 0)
            return 0;
        if (set->sc.step_act % set->so.pnts[i].step == 0)
            return 1;
    }

    for (i = 0; i < set->so.nsums; i++) {
        if (set->so.sums[i].step == 0)
            return 0;
        if (set->sc.step_act % set->so.sums[i].step == 0)
            return 1;
    }

    for (i = 0; i < set->so.nforces; i++) {
        if (set->so.forces[i].step == 0)
            return 0;
        if (set->sc.step_act % set->so.forces[i].step == 0)
            return 1;
    }

    for (i = 0; i < set->so.ncubs; i++) {
        if (set->so.cubs[i].step == 0)
            return 0;
        if (set->sc.step_act % set->so.cubs[i].step == 0 && !(set->sc.step_act<set->so.cubs[i].start || set->sc.step_act>set->so.cubs[i].stop))
            return 1;
    }

    for (i = 0; i < set->so.nplns; i++) {
        if (set->so.plns[i].step == 0)
            return 0;
        if (set->sc.step_act % set->so.plns[i].step == 0 && !(set->sc.step_act<set->so.plns[i].start || set->sc.step_act>set->so.plns[i].stop))
            return 1;
    }

    return 0;
}

void
sv_pool_allocate_output(SvPool *mp, SvSet *set)
{
    gint i, n;
    gint nabss;

    mp->out->outsumdata = (gdouble **)g_malloc(set->so.nsums * sizeof(gdouble *));
    for (n = 0; n < set->so.nsums; n++) {
        mp->out->outsumdata[n] = (gdouble *)g_malloc(set->sc.nsteps * sizeof(gdouble));
        for (i = 0; i < set->sc.nsteps; i++)
            mp->out->outsumdata[n][i] = 0;
    }

    mp->out->outpointdata = (gdouble **)g_malloc(set->so.npnts * sizeof(gdouble *));
    for (n = 0; n < set->so.npnts; n++) {
        mp->out->outpointdata[n] = (gdouble *)g_malloc(set->sc.nsteps * 6 * sizeof(gdouble));
        for (i = 0; i < (6 * set->sc.nsteps); i++)
            mp->out->outpointdata[n][i] = 0;
    }

    mp->out->outforcedata = (gdouble **)g_malloc(set->so.nforces * sizeof(gdouble *));
    for (n = 0; n < set->so.nforces; n++) {
        mp->out->outforcedata[n] = (gdouble *)g_malloc(set->sc.nsteps * 3 * sizeof(gdouble));
        for (i = 0; i < (3 * set->sc.nsteps); i++)
            mp->out->outforcedata[n][i] = 0;
    }

    nabss = 0;
    for (n = 0; n < set->so.ncubs; n++)
        if (set->so.cubs[n].component == SV_OVOLUME_ABS || set->so.cubs[n].component == SV_OVOLUME_SUMALL || set->so.cubs[n].component == SV_OVOLUME_MAXALL) nabss++;

    printf("Allocating output: %d %d %d, %d\n", set->sp.xres, set->sp.yres, set->sp.zres, nabss);
    mp->out->outabsdata = (SvDCube **)g_malloc(nabss * sizeof(SvDCube *));
    for (n = 0; n < nabss; n++)
        mp->out->outabsdata[n] = sv_dcube_new(set->sp.xres, set->sp.yres, set->sp.zres, set->sp.xres, set->sp.yres, set->sp.zres, 1);
}


static void output_vtk_header(FILE *fw, gint xres, gint yres, gint zres, gchar *description)
{
    fprintf(fw, "# vtk DataFile Version 2.0\n");
    fprintf(fw, "%s\n", description);
    fprintf(fw, "ASCII\n");
    fprintf(fw, "DATASET STRUCTURED_POINTS\n");
    fprintf(fw, "DIMENSIONS %d %d %d\n", xres, yres, zres);
    fprintf(fw, "ASPECT_RATIO 1 1 1\n");
    fprintf(fw, "ORIGIN 0 0 0\n");
    fprintf(fw, "POINT_DATA %d\n", xres*yres*zres);
    fprintf(fw, "SCALARS volume_scalars double 1\n");
    fprintf(fw, "LOOKUP_TABLE default\n");
}

/*outputs global_material or another global field for debugging*/
void output_ivtk(int ***data, int xres, int yres, int zres, char *filename, char *description)
{
    FILE *fw;
    int i, j, k;

    fw = fopen(filename, "w");
    if (!fw) {
        fprintf(stderr, "Error: VTK output: cannot open file %s for writing\n", filename);
    }
    output_vtk_header(fw, xres, yres, zres, description);

    for (k = 0; k < zres; k++) {
        for (j = 0; j < yres; j++) {
            for (i = 0; i < xres; i++) {
                fprintf(fw, "%d\n", data[i][j][k]);
            }
        }
    }
    fclose(fw);
}

static void output_vtk(double ***data, int xres, int yres, int zres, char *filename, char *description)
{
    FILE *fw;
    int i, j, k;

    fw = fopen(filename, "w");
    if (!fw) {
        fprintf(stderr, "Error: VTK output: cannot open file %s for writing\n", filename);
    }
    output_vtk_header(fw, xres, yres, zres, description);

    for (k = 0; k < zres; k++) {
        for (j = 0; j < yres; j++) {
            for (i = 0; i < xres; i++) {
                fprintf(fw, "%g\n", data[i][j][k]);
            }
        }
    }

    fclose(fw);
}
static void output_power_vtk(double ***adata, double ***bdata, double ***cdata, int xres, int yres, int zres, char *filename, char *description)
{
    FILE *fw;
    int i, j, k;

    fw = fopen(filename, "w");
    if (!fw) {
        fprintf(stderr, "Error: VTK output: cannot open file %s for writing\n", filename);
    }
    output_vtk_header(fw, xres, yres, zres, description);

    for (k = 0; k < zres; k++) {
        for (j = 0; j < yres; j++) {
            for (i = 0; i < xres; i++) {
                fprintf(fw, "%g\n", sqrt(adata[i][j][k] * adata[i][j][k] + bdata[i][j][k] * bdata[i][j][k] + cdata[i][j][k] * cdata[i][j][k]));
            }
        }
    }

    fclose(fw);
}

static void output_material_vtk(SvPool *mp, SvSet *set, int xres, int yres, int zres, char *filename, char *description, SvOutputVolumeType component)
{
    FILE *fw;
    int i, j, k;
    double epsilon, sigma, mu, sigast;
    int matindex;

    fw = fopen(filename, "w");
    if (!fw) {
        fprintf(stderr, "Error: VTK output: cannot open file %s for writing\n", filename);
    }
    output_vtk_header(fw, xres, yres, zres, description);

    for (k = 0; k < zres; k++) {
        for (j = 0; j < yres; j++) {
            for (i = 0; i < xres; i++) {

                sv_yee_data_get_epsilon_sigma(mp->d, set, mp->mats, mp->nmat, i, j, k, &epsilon, &sigma);
                matindex = sv_yee_data_get_mu_sigast(mp->d, set, mp->mats, mp->nmat, i, j, k, &mu, &sigast);

                epsilon /= EPSILON_0;
                mu /= MU_0;

                if (component == SV_OVOLUME_EPSILON)
                    fprintf(fw, "%g\n", epsilon);
                else if (component == SV_OVOLUME_SIGMA)
                    fprintf(fw, "%g\n", sigma);
                else if (component == SV_OVOLUME_MU)
                    fprintf(fw, "%g\n", mu);
                else if (component == SV_OVOLUME_SIGAST)
                    fprintf(fw, "%g\n", sigast);
                else if (component == SV_OVOLUME_MATTYPE)
                    fprintf(fw, "%d\n", mp->mats[matindex].type);

            }
        }
    }

    fclose(fw);
}


void
sv_pool_getsum(SvPool *mp, SvSet *set)
{
    gint i, j, k, n;
    gdouble val, epsilon, sigma;

    if (set->so.nsums == 0)
        return;

    if (set->sc.verbose > 1) {
        printf("Collecting sums...   ");
        fflush(stdout);
    }

    for (n = 0; n < set->so.nsums; n++) {
        if (set->plan.gpumode == SV_GPUMODE_NONE) {
            for (i = set->so.sums[n].box_i0; i < set->so.sums[n].box_in; i++) {
                for (j = set->so.sums[n].box_j0; j < set->so.sums[n].box_jn; j++) {
                    for (k = set->so.sums[n].box_k0; k < set->so.sums[n].box_kn; k++) {

                        sv_yee_data_get_epsilon_sigma(mp->d, set, mp->mats, mp->nmat, i, j, k, &epsilon, &sigma);

                        if (set->so.sums[n].layered_epsilon == -1 || (fabs(set->so.sums[n].layered_epsilon - epsilon) < (EPSILON_0 / 100.0))) {
                            val = 0;
                            if (set->so.sums[n].component == SV_SUM_ALL || set->so.sums[n].component == SV_SUM_MAX || set->so.sums[n].component == SV_SUM_ABS || set->so.sums[n].component == SV_SUM_EX)
                                val += mp->d->ex->data[i][j][k] * mp->d->ex->data[i][j][k];
                            if (set->so.sums[n].component == SV_SUM_ALL || set->so.sums[n].component == SV_SUM_MAX || set->so.sums[n].component == SV_SUM_ABS || set->so.sums[n].component == SV_SUM_EY)
                                val += mp->d->ey->data[i][j][k] * mp->d->ey->data[i][j][k];
                            if (set->so.sums[n].component == SV_SUM_ALL || set->so.sums[n].component == SV_SUM_MAX || set->so.sums[n].component == SV_SUM_ABS || set->so.sums[n].component == SV_SUM_EZ)
                                val += mp->d->ez->data[i][j][k] * mp->d->ez->data[i][j][k];

                            if (set->so.sums[n].component == SV_SUM_MAX)
                                mp->out->outsumdata[n][set->sc.step_act] = MAX(mp->out->outsumdata[n][set->sc.step_act], val);
                            else if (set->so.sums[n].component == SV_SUM_ABS)
                                mp->out->outsumdata[n][set->sc.step_act] += sigma*val;
                            else
                                mp->out->outsumdata[n][set->sc.step_act] += val;
                        } // if
                    } // k
                } // j
            } // i
        } // if
    } // n

    if (set->sc.verbose > 1)
        printf("done.\n");
}

static void
get_all(SvDCube *ex, SvDCube *ey, SvDCube *ez, GwyDataField *df, gint i, gint j, gint k)
{
    gint ii, xres, yres;
    GwyDataField *dfy, *dfz;
    gdouble *xdata, *ydata, *zdata;
 
    dfy = gwy_data_field_new(5, 5, 5, 5, 1);
    dfz = gwy_data_field_new(5, 5, 5, 5, 1);

    sv_dcube_get_datafield(ex, df, i, j, k);
    sv_dcube_get_datafield(ey, dfy, i, j, k);
    sv_dcube_get_datafield(ez, dfz, i, j, k);
  
    xres = gwy_data_field_get_xres(df);
    yres = gwy_data_field_get_yres(df);

    xdata = gwy_data_field_get_data(df);
    ydata = gwy_data_field_get_data(dfy);
    zdata = gwy_data_field_get_data(dfz);


    for (ii=0; ii<(xres*yres); ii++) {
        xdata[ii] = xdata[ii]*xdata[ii] + ydata[ii]*ydata[ii] + zdata[ii]*zdata[ii];    
    }
    gwy_data_field_data_changed(df);

    g_object_unref(dfy);
    g_object_unref(dfz);

}
 
void
sv_pool_getpoints(SvPool *mp, SvSet *set)
{
    gint i;
    for (i = 0; i < set->so.npnts; i++) {
        mp->out->outpointdata[i][set->sc.step_act * 6] = mp->d->ex->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k];
        mp->out->outpointdata[i][set->sc.step_act * 6 + 1] = mp->d->ey->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k];
        mp->out->outpointdata[i][set->sc.step_act * 6 + 2] = mp->d->ez->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k];
        mp->out->outpointdata[i][set->sc.step_act * 6 + 3] = mp->d->hx->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k];
        mp->out->outpointdata[i][set->sc.step_act * 6 + 4] = mp->d->hy->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k];
        mp->out->outpointdata[i][set->sc.step_act * 6 + 5] = mp->d->hz->data[set->so.pnts[i].i][set->so.pnts[i].j][set->so.pnts[i].k];
    }
}

void
prod(double t[3][3], double n[3], double f[3])
{
    f[0] = t[0][0] * n[0] + t[0][1] * n[1] + t[0][2] * n[2];
    f[1] = t[1][0] * n[0] + t[1][1] * n[1] + t[1][2] * n[2];
    f[2] = t[2][0] * n[0] + t[2][1] * n[1] + t[2][2] * n[2];
}

static void
get_matprops(SvPool *mp, SvSet *set, gint i, gint j, gint k, gdouble *epsilon, gdouble *mu)
{
    *epsilon = sv_yee_data_get_epsilon(mp->d, set, mp->mats, mp->nmat, i, j, k);
    *mu = sv_yee_data_get_mu(mp->d, set, mp->mats, mp->nmat, i, j, k);
}

static void
assemble_matrix(gdouble t[3][3], gdouble ex, gdouble ey, gdouble ez, gdouble hx, gdouble hy, gdouble hz, gdouble epsilon, gdouble mu)
{
    gdouble val;

    val = (epsilon*(ex*ex + ey*ey + ez*ez) + mu*(hx*hx + hy*hy + hz*hz)) / 2.0;
    t[0][0] = epsilon*ex*ex + mu*hx*hx - val;
    t[1][1] = epsilon*ey*ey + mu*hy*hy - val;
    t[2][2] = epsilon*ez*ez + mu*hz*hz - val;
    t[0][1] = t[1][0] = epsilon*ex*ey + mu*hx*hy;
    t[0][2] = t[2][0] = epsilon*ex*ez + mu*hx*hz;
    t[1][2] = t[2][1] = epsilon*ey*ez + mu*hy*hz;
}

void
sv_pool_getforce(SvPool *mp, SvSet *set)
{
    gint n, i, j, k;
    gdouble faccx, faccy, faccz;
    gdouble t[3][3];
    gdouble normal[3], f[3];
    gdouble epsilon, mu, ex, ey, ez, hx, hy, hz;

    for (n = 0; n < set->so.nforces; n++) {
        faccx = faccy = faccz = 0;
        for (j = set->so.forces[n].box_j0; j < set->so.forces[n].box_jn; j++) {
            for (k = set->so.forces[n].box_k0; k < set->so.forces[n].box_kn; k++) {
                /*Local values of ex, ey, ez at x0. Normal of this facet is 1 0 0*/
                /*T_ij*/

                get_matprops(mp, set, set->so.forces[n].box_i0, j, k, &epsilon, &mu);

                ex = mp->d->ex->data[set->so.forces[n].box_i0][j][k];
                ey = mp->d->ey->data[set->so.forces[n].box_i0][j][k];
                ez = mp->d->ez->data[set->so.forces[n].box_i0][j][k];
                hx = mp->d->hx->data[set->so.forces[n].box_i0][j][k];
                hy = mp->d->hy->data[set->so.forces[n].box_i0][j][k];
                hz = mp->d->hz->data[set->so.forces[n].box_i0][j][k];

                assemble_matrix(t, ex, ey, ez, hx, hy, hz, epsilon, mu);

                normal[0] = 1; normal[1] = 0; normal[2] = 0;
                prod(t, normal, f);

                faccx += f[0];
                faccy += f[1];
                faccz += f[2];

                get_matprops(mp, set, set->so.forces[n].box_in, j, k, &epsilon, &mu);

                ex = mp->d->ex->data[set->so.forces[n].box_in][j][k];
                ey = mp->d->ey->data[set->so.forces[n].box_in][j][k];
                ez = mp->d->ez->data[set->so.forces[n].box_in][j][k];
                hx = mp->d->hx->data[set->so.forces[n].box_in][j][k];
                hy = mp->d->hy->data[set->so.forces[n].box_in][j][k];
                hz = mp->d->hz->data[set->so.forces[n].box_in][j][k];

                assemble_matrix(t, ex, ey, ez, hx, hy, hz, epsilon, mu);
                normal[0] = -1; normal[1] = 0; normal[2] = 0;
                prod(t, normal, f);
                faccx += f[0];
                faccy += f[1];
                faccz += f[2];
            } // k
        } // j
        for (i = set->so.forces[n].box_i0; i < set->so.forces[n].box_in; i++) {
            for (k = set->so.forces[n].box_k0; k < set->so.forces[n].box_kn; k++) {
                /*Local values of ex, ey, ez at x0. Normal of this facet is 1 0 0*/
                /*T_ij*/

                get_matprops(mp, set, i, set->so.forces[n].box_j0, k, &epsilon, &mu);

                ex = mp->d->ex->data[i][set->so.forces[n].box_j0][k];
                ey = mp->d->ey->data[i][set->so.forces[n].box_j0][k];
                ez = mp->d->ez->data[i][set->so.forces[n].box_j0][k];
                hx = mp->d->hx->data[i][set->so.forces[n].box_j0][k];
                hy = mp->d->hy->data[i][set->so.forces[n].box_j0][k];
                hz = mp->d->hz->data[i][set->so.forces[n].box_j0][k];

                assemble_matrix(t, ex, ey, ez, hx, hy, hz, epsilon, mu);
                normal[0] = 0; normal[1] = 1; normal[2] = 0;
                prod(t, normal, f);
                faccx += f[0];
                faccy += f[1];
                faccz += f[2];

                get_matprops(mp, set, i, set->so.forces[n].box_jn, k, &epsilon, &mu);

                ex = mp->d->ex->data[i][set->so.forces[n].box_jn][k];
                ey = mp->d->ey->data[i][set->so.forces[n].box_jn][k];
                ez = mp->d->ez->data[i][set->so.forces[n].box_jn][k];
                hx = mp->d->hx->data[i][set->so.forces[n].box_jn][k];
                hy = mp->d->hy->data[i][set->so.forces[n].box_jn][k];
                hz = mp->d->hz->data[i][set->so.forces[n].box_jn][k];

                assemble_matrix(t, ex, ey, ez, hx, hy, hz, epsilon, mu);
                normal[0] = 0; normal[1] = -1; normal[2] = 0;
                prod(t, normal, f);
                faccx += f[0];
                faccy += f[1];
                faccz += f[2];
            } // k
        } // i
        for (i = set->so.forces[n].box_i0; i < set->so.forces[n].box_in; i++) {
            for (j = set->so.forces[n].box_j0; j < set->so.forces[n].box_jn; j++) {

                get_matprops(mp, set, i, j, set->so.forces[n].box_k0, &epsilon, &mu);

                ex = mp->d->ex->data[i][j][set->so.forces[n].box_k0];
                ey = mp->d->ey->data[i][j][set->so.forces[n].box_k0];
                ez = mp->d->ez->data[i][j][set->so.forces[n].box_k0];
                hx = mp->d->hx->data[i][j][set->so.forces[n].box_k0];
                hy = mp->d->hy->data[i][j][set->so.forces[n].box_k0];
                hz = mp->d->hz->data[i][j][set->so.forces[n].box_k0];

                assemble_matrix(t, ex, ey, ez, hx, hy, hz, epsilon, mu);
                normal[0] = 0; normal[1] = 0; normal[2] = 1;
                prod(t, normal, f);
                faccx += f[0];
                faccy += f[1];
                faccz += f[2];

                get_matprops(mp, set, i, j, set->so.forces[n].box_kn, &epsilon, &mu);

                ex = mp->d->ex->data[i][j][set->so.forces[n].box_kn];
                ey = mp->d->ey->data[i][j][set->so.forces[n].box_kn];
                ez = mp->d->ez->data[i][j][set->so.forces[n].box_kn];
                hx = mp->d->hx->data[i][j][set->so.forces[n].box_kn];
                hy = mp->d->hy->data[i][j][set->so.forces[n].box_kn];
                hz = mp->d->hz->data[i][j][set->so.forces[n].box_kn];

                assemble_matrix(t, ex, ey, ez, hx, hy, hz, epsilon, mu);
                normal[0] = 0; normal[1] = 0; normal[2] = -1;
                prod(t, normal, f);
                faccx += f[0];
                faccy += f[1];
                faccz += f[2];
            } // j
        } // i

        mp->out->outforcedata[n][set->sc.step_act * 3] = faccx;
        mp->out->outforcedata[n][set->sc.step_act * 3 + 1] = faccy;
        mp->out->outforcedata[n][set->sc.step_act * 3 + 2] = faccz;
    } // n
}


void
dataline_get_spectrum(GwyDataLine *dline, GwyDataLine *spectrum)
{

    gint i, res;
    GwyDataLine *iin, *rout, *iout;
    gdouble *linedata, *rdata, *idata;

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

}

void
sv_pool_output(SvPool *mp, SvSet *set)
{
    gint i, j, k, n, col, row, indval, nav, nsg;
    gboolean outim = FALSE;
    gchar key[32], titlekey[32];
    gchar description[100];
    GByteArray *buffer;
    GwyDataField *df = NULL;
    gdouble *xdata = NULL, *ydata = NULL, avx, avy, avz;
    GwyGraphModel *gm = NULL;
    GwyGraphCurveModel *gc = NULL;
    gdouble val = 0;
    gint matindex;
    gchar filename[100];
    gdouble epsilon, sigma, mu, sigast;
    GwyDataLine *dline, *spectrumx, *spectrumy, *spectrumz;
    gdouble *dlinedata, *sdatax, *sdatay, *sdataz;

    //    gdouble f[3], t[3][3], tav[3][3], normal[3], faccx, faccy, faccz;

    FILE *fh;
    gchar *component_string[] = {"Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "All", "Cur", "Epsilon", "Sigma", "Mu", "Sigast", "Material"};

    if (set->sc.verbose > 1) {
        printf("Running output...   ");
        fflush(stdout);
    }

    //    a = sv_dcube_get_sum(mp->d->hx);
    //    b = sv_dcube_get_sum(mp->d->hy);
    //    c = sv_dcube_get_sum(mp->d->hz);

    //   fprintf(stderr, "Set diff: (%g %g %g)\n", a, b, c); 

    //    fprintf(stderr, "Set diff: (%g %g %g)\n", (gdouble)(mp->d->hx->xres-1)*(mp->d->hx->yres-1)*(mp->d->hx->zres-1)*1.0 - a, 
    //                               (gdouble)(mp->d->hx->xres-1)*(mp->d->hx->yres-1)*(mp->d->hx->zres-1)*2.0 - b,
    //                               (gdouble)(mp->d->hx->xres-1)*(mp->d->hx->yres-1)*(mp->d->hx->zres-1)*3.0 - c);

    for (i = 0; i < set->so.npnts; i++) {
        if (set->so.pnts[i].step == 0)
            continue;

        if (set->sc.step_act % set->so.pnts[i].step != 0 && set->sc.step_act < (set->sc.nsteps - 1))
            continue;

        fh = fopen(set->so.pnts[i].filebase, "w");
        if (fh == NULL) {
            fprintf(stderr, "Error: point output: cannot open file %s for writing\n", set->so.pnts[i].filebase);
            continue;
        }

        if (set->so.pnts[i].component == SV_COMP_ALLFFT) {
            if (set->sc.step_act == (set->sc.nsteps - 1)) {
                dline = gwy_data_line_new(set->so.pnts[i].stop - set->so.pnts[i].start,
                    (set->so.pnts[i].stop - set->so.pnts[i].start)*mp->set->plan.dt, FALSE);
                dlinedata = gwy_data_line_get_data(dline);

                for (j = set->so.pnts[i].start; j < set->so.pnts[i].stop; j++)
                    dlinedata[j - set->so.pnts[i].start] = mp->out->outpointdata[i][j * 6];
                gwy_data_line_data_changed(dline);
                spectrumx = gwy_data_line_new_alike(dline, TRUE);
                dataline_get_spectrum(dline, spectrumx);


                for (j = set->so.pnts[i].start; j < set->so.pnts[i].stop; j++)
                    dlinedata[j - set->so.pnts[i].start] = mp->out->outpointdata[i][j * 6 + 1];
                gwy_data_line_data_changed(dline);
                spectrumy = gwy_data_line_new_alike(dline, TRUE);
                dataline_get_spectrum(dline, spectrumy);

                for (j = set->so.pnts[i].start; j < set->so.pnts[i].stop; j++)
                    dlinedata[j - set->so.pnts[i].start] = mp->out->outpointdata[i][j * 6 + 2];
                gwy_data_line_data_changed(dline);
                spectrumz = gwy_data_line_new_alike(dline, TRUE);
                dataline_get_spectrum(dline, spectrumz);

                sdatax = gwy_data_line_get_data(spectrumx);
                sdatay = gwy_data_line_get_data(spectrumy);
                sdataz = gwy_data_line_get_data(spectrumz);

                for (j = 1; j < gwy_data_line_get_res(spectrumx); j++)
                    fprintf(fh, "%g %g %g %g %g\n", LIGHT_SPEED / ((gdouble)j / (((gdouble)(set->so.pnts[i].stop - set->so.pnts[i].start) / 2.0)*mp->set->plan.dt*2.0)),
                            sdatax[j] + sdatay[j] + sdataz[j],
                            sqrt(sdatax[j]),
                            sqrt(sdatay[j]),
                            sqrt(sdataz[j]));

            }
        } else {
            for (j = 0; j < set->sc.step_act; j++) {
                if (set->so.pnts[i].component == SV_COMP_ALL) {
                    fprintf(fh, "%d %g %g %g %g %g %g\n", j,
                            mp->out->outpointdata[i][j * 6], mp->out->outpointdata[i][j * 6 + 1],
                            mp->out->outpointdata[i][j * 6 + 2], mp->out->outpointdata[i][j * 6 + 3],
                            mp->out->outpointdata[i][j * 6 + 4], mp->out->outpointdata[i][j * 6 + 5]);
                } else {
                    switch (set->so.pnts[i].component) {
                        case SV_COMP_EX:
                            val = mp->out->outpointdata[i][j * 6];
                            break;

                        case SV_COMP_EY:
                            val = mp->out->outpointdata[i][j * 6 + 1];
                            break;

                        case SV_COMP_EZ:
                            val = mp->out->outpointdata[i][j * 6 + 2];
                            break;

                        case SV_COMP_HX:
                            val = mp->out->outpointdata[i][j * 6 + 3];
                            break;

                        case SV_COMP_HY:
                            val = mp->out->outpointdata[i][j * 6 + 4];
                            break;

                        case SV_COMP_HZ:
                            val = mp->out->outpointdata[i][j * 6 + 5];
                            break;

                        default:
                            break;
                    }

                    fprintf(fh, "%d %g\n", j, val);
                } // else
            } // j
        } // else

        fclose(fh);

    } // i

    for (n = 0; n < set->so.nsums; n++) {
        if (set->so.sums[i].step == 0)
            continue;

        if (set->sc.step_act % set->so.sums[n].step != 0 && set->sc.step_act < (set->sc.nsteps - 1))
            continue;

        fh = fopen(set->so.sums[n].filename, "w");
        if (fh == NULL) {
            fprintf(stderr, "Error: sum: cannot open file %s for writing\n", set->so.sums[i].filename);
            continue;
        }

        for (i = 0; i < set->sc.nsteps; i++) {
            if (set->so.sums[n].component == SV_SUM_MAX)
                fprintf(fh, "%d %g\n", i, mp->out->outsumdata[n][i]);
            else
                fprintf(fh, "%d %g\n", i, mp->out->outsumdata[n][i] * set->sp.dx*set->sp.dy*set->sp.dz);
        }
        fclose(fh);
    }

    for (n = 0; n < set->so.ncubs; n++) {
        if (set->so.cubs[i].step == 0)
            continue;
        if (set->sc.step_act % set->so.cubs[n].step != 0 && set->sc.step_act < (set->sc.nsteps - 1))
            continue;
        if (set->sc.step_act < set->so.cubs[n].start || set->sc.step_act > set->so.cubs[n].stop)
            continue;

        if (set->so.cubs[n].component == SV_OVOLUME_ABS || set->so.cubs[n].component == SV_OVOLUME_SUMALL || set->so.cubs[n].component == SV_OVOLUME_MAXALL) {

#pragma omp parallel default(shared) private(i, j, k, sigma)
#pragma omp for nowait

            for (i = 0; i < mp->d->ex->xres; i++) {
                for (j = 0; j < mp->d->ex->yres; j++) {
                    for (k = 0; k < mp->d->ex->zres; k++) {
                        if (set->so.cubs[n].component == SV_OVOLUME_ABS)
                            sigma = sv_yee_data_get_sigma(mp->d, set, mp->mats, mp->nmat, i, j, k);
                        else
                            sigma = 1;

                        if (set->so.cubs[n].component == SV_OVOLUME_MAXALL)
                            mp->out->outabsdata[set->so.cubs[n].k]->data[i][j][k] = MAX(mp->out->outabsdata[set->so.cubs[n].k]->data[i][j][k], (mp->d->ex->data[i][j][k] * mp->d->ex->data[i][j][k] + mp->d->ey->data[i][j][k] * mp->d->ey->data[i][j][k] + mp->d->ez->data[i][j][k] * mp->d->ez->data[i][j][k]));
                        else
                            mp->out->outabsdata[set->so.cubs[n].k]->data[i][j][k] += sigma * (mp->d->ex->data[i][j][k] * mp->d->ex->data[i][j][k] + mp->d->ey->data[i][j][k] * mp->d->ey->data[i][j][k] + mp->d->ez->data[i][j][k] * mp->d->ez->data[i][j][k]);
                    } // k
                } // j
            } // i
        } // if

        if (set->so.cubs[n].format == 1) { //this means ascii mode here
            if (set->so.cubs[n].component == SV_OVOLUME_ABS || set->so.cubs[n].component == SV_OVOLUME_SUMALL || set->so.cubs[n].component == SV_OVOLUME_MAXALL)
                g_snprintf(filename, 100, "%s", set->so.cubs[n].filebase);
            else
                g_snprintf(filename, 100, "%s_%.4d", set->so.cubs[n].filebase, set->sc.step_act);

            fh = fopen(filename, "w");
            if (fh == NULL) {
                fprintf(stderr, "Error: volume: cannot open file %s for writing\n", filename);
                continue;
            }

            fprintf(fh, "%d %d %d %d\n", mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, set->so.cubs[n].component);

            if ((set->so.cubs[n].component == SV_OVOLUME_MAT || set->so.cubs[n].component == SV_OVOLUME_MATTYPE) && mp->nmat == 0)
                fprintf(fh, "No data available\n");
            else {
                for (k = 0; k < mp->d->ex->zres; k++) {
                    for (j = 0; j < mp->d->ex->yres; j++) {
                        for (i = 0; i < mp->d->ex->xres; i++) {
                            sv_yee_data_get_epsilon_sigma(mp->d, set, mp->mats, mp->nmat, i, j, k, &epsilon, &sigma);
                            sv_yee_data_get_mu_sigast(mp->d, set, mp->mats, mp->nmat, i, j, k, &mu, &sigast);

                            epsilon /= EPSILON_0;
                            mu /= MU_0;

                            if (set->so.cubs[n].component == SV_OVOLUME_EX)
                                fprintf(fh, "%g ", mp->d->ex->data[i][j][k]);
                            else if (set->so.cubs[n].component == SV_OVOLUME_EY)
                                fprintf(fh, "%g ", mp->d->ey->data[i][j][k]);
                            else if (set->so.cubs[n].component == SV_OVOLUME_EZ)
                                fprintf(fh, "%g ", mp->d->ez->data[i][j][k]);
                            else if (set->so.cubs[n].component == SV_OVOLUME_HX)
                                fprintf(fh, "%g ", mp->d->ex->data[i][j][k]);

                            else if (set->so.cubs[n].component == SV_OVOLUME_HY)
                                fprintf(fh, "%g ", mp->d->ey->data[i][j][k]);
                            else if (set->so.cubs[n].component == SV_OVOLUME_HZ)
                                fprintf(fh, "%g ", mp->d->ez->data[i][j][k]);
                            else if (set->so.cubs[n].component == SV_OVOLUME_ALL)
                                fprintf(fh, "%g ", sqrt(mp->d->ex->data[i][j][k] * mp->d->ex->data[i][j][k] + mp->d->ey->data[i][j][k] * mp->d->ey->data[i][j][k] + mp->d->ez->data[i][j][k] * mp->d->ez->data[i][j][k]));
                            else if (set->so.cubs[n].component == SV_OVOLUME_EPSILON)
                                fprintf(fh, "%g ", epsilon);
                            else if (set->so.cubs[n].component == SV_OVOLUME_SIGMA)
                                fprintf(fh, "%g ", sigma);
                            else if (set->so.cubs[n].component == SV_OVOLUME_MU)
                                fprintf(fh, "%g ", mu);
                            else if (set->so.cubs[n].component == SV_OVOLUME_SIGAST)
                                fprintf(fh, "%g ", sigast);
                            else if (set->so.cubs[n].component == SV_OVOLUME_MAT)
                                fprintf(fh, "%d ", mp->d->mat->data[i][j][k]);
                            else if (set->so.cubs[n].component == SV_OVOLUME_MATTYPE)
                                fprintf(fh, "%d ", mp->mats[mp->d->mat->data[i][j][k]].type);
                            else if (set->so.cubs[n].component == SV_OVOLUME_ABS || set->so.cubs[n].component == SV_OVOLUME_SUMALL || set->so.cubs[n].component == SV_OVOLUME_MAXALL)
                                fprintf(fh, "%g ", mp->out->outabsdata[set->so.cubs[n].k]->data[i][j][k]);
                        } // i
                        fprintf(fh, "\n");
                    } // j
                    fprintf(fh, "\n");
                } // k
            } // else

            fclose(fh);
        } else if (set->so.cubs[n].format == 2) { //this means vtk mode here
            if (set->so.cubs[n].component == SV_OVOLUME_ABS || set->so.cubs[n].component == SV_OVOLUME_SUMALL || set->so.cubs[n].component == SV_OVOLUME_MAXALL)
                g_snprintf(filename, 100, "%s.vtk", set->so.cubs[n].filebase);
            else
                g_snprintf(filename, 100, "%s_%.4d.vtk", set->so.cubs[n].filebase, set->sc.step_act);

            if (set->so.cubs[n].component == SV_OVOLUME_EX)
                output_vtk(mp->d->ex->data, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "Ex");
            else if (set->so.cubs[n].component == SV_OVOLUME_EY)
                output_vtk(mp->d->ey->data, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "Ey");
            else if (set->so.cubs[n].component == SV_OVOLUME_EZ)
                output_vtk(mp->d->ez->data, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "Ez");
            else if (set->so.cubs[n].component == SV_OVOLUME_HX)
                output_vtk(mp->d->hx->data, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "Hx");
            else if (set->so.cubs[n].component == SV_OVOLUME_HY)
                output_vtk(mp->d->hy->data, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "Hy");
            else if (set->so.cubs[n].component == SV_OVOLUME_HZ)
                output_vtk(mp->d->hz->data, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "Hz");
            else if (set->so.cubs[n].component == SV_OVOLUME_ALL)
                output_power_vtk(mp->d->ex->data, mp->d->ey->data, mp->d->ez->data, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "All");
            else if (set->so.cubs[n].component == SV_OVOLUME_EPSILON)
                output_material_vtk(mp, set, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "Epsilon", set->so.cubs[n].component);
            else if (set->so.cubs[n].component == SV_OVOLUME_SIGMA)
                output_material_vtk(mp, set, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "Sigma", set->so.cubs[n].component);
            else if (set->so.cubs[n].component == SV_OVOLUME_MU)
                output_material_vtk(mp, set, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "Mu", set->so.cubs[n].component);
            else if (set->so.cubs[n].component == SV_OVOLUME_SIGAST)
                output_material_vtk(mp, set, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "Sigast", set->so.cubs[n].component);
            else if (set->so.cubs[n].component == SV_OVOLUME_MAT)
                output_ivtk(mp->d->mat->data, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "Material");
            else if (set->so.cubs[n].component == SV_OVOLUME_MATTYPE)
                output_material_vtk(mp, set, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "Mattype", set->so.cubs[n].component);
            else if (set->so.cubs[n].component == SV_OVOLUME_ABS)
                output_vtk(mp->out->outabsdata[set->so.cubs[n].k]->data, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "Abs");
            else if (set->so.cubs[n].component == SV_OVOLUME_SUMALL)
                output_vtk(mp->out->outabsdata[set->so.cubs[n].k]->data, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "Sumall");
            else if (set->so.cubs[n].component == SV_OVOLUME_MAXALL)
                output_vtk(mp->out->outabsdata[set->so.cubs[n].k]->data, mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres, filename, "Maxall");
        } else { //this means binary mode
            if ((set->so.cubs[n].component == SV_OVOLUME_ABS || set->so.cubs[n].component == SV_OVOLUME_SUMALL || set->so.cubs[n].component == SV_OVOLUME_MAXALL) && !(set->sc.step_act > (set->so.cubs[n].stop - set->so.cubs[n].step)))
                continue;

            if (set->so.cubs[n].component == SV_OVOLUME_ABS || set->so.cubs[n].component == SV_OVOLUME_SUMALL || set->so.cubs[n].component == SV_OVOLUME_MAXALL)
                g_snprintf(filename, 100, "%s.raw", set->so.cubs[n].filebase);
            else
                g_snprintf(filename, 100, "%s_%.4d.raw", set->so.cubs[n].filebase, set->sc.step_act);

            fh = fopen(filename, "wb");
            if (fh == NULL) {
                fprintf(stderr, "Error: volume: cannot open file %s for writing\n", filename);
                continue;
            }

            if ((set->so.cubs[n].component == SV_OVOLUME_MAT || set->so.cubs[n].component == SV_OVOLUME_MATTYPE) && mp->nmat == 0)
                fprintf(stderr, "No data available for binary output.\n");
            else {
                for (k = 0; k < mp->d->ex->zres; k++) {
                    for (j = 0; j < mp->d->ex->yres; j++) {
                        for (i = 0; i < mp->d->ex->xres; i++) {
                            sv_yee_data_get_epsilon_sigma(mp->d, set, mp->mats, mp->nmat, i, j, k, &epsilon, &sigma);
                            matindex = sv_yee_data_get_mu_sigast(mp->d, set, mp->mats, mp->nmat, i, j, k, &mu, &sigast);

                            epsilon /= EPSILON_0;
                            mu /= MU_0;

                            if (set->so.cubs[n].component == SV_OVOLUME_EX)
                                val = mp->d->ex->data[i][j][k];
                            else if (set->so.cubs[n].component == SV_OVOLUME_EY)
                                val = mp->d->ey->data[i][j][k];
                            else if (set->so.cubs[n].component == SV_OVOLUME_EZ)
                                val = mp->d->ez->data[i][j][k];
                            else if (set->so.cubs[n].component == SV_OVOLUME_HX)
                                val = mp->d->hx->data[i][j][k];
                            else if (set->so.cubs[n].component == SV_OVOLUME_HY)
                                val = mp->d->hy->data[i][j][k];
                            else if (set->so.cubs[n].component == SV_OVOLUME_HZ)
                                val = mp->d->hz->data[i][j][k];
                            else if (set->so.cubs[n].component == SV_OVOLUME_ALL)
                                val = sqrt(mp->d->ex->data[i][j][k] * mp->d->ex->data[i][j][k] + mp->d->ey->data[i][j][k] * mp->d->ey->data[i][j][k] + mp->d->ez->data[i][j][k] * mp->d->ez->data[i][j][k]);
                            else if (set->so.cubs[n].component == SV_OVOLUME_EPSILON)
                                val = epsilon;
                            else if (set->so.cubs[n].component == SV_OVOLUME_SIGMA)
                                val = sigma;
                            else if (set->so.cubs[n].component == SV_OVOLUME_MU)
                                val = mu;
                            else if (set->so.cubs[n].component == SV_OVOLUME_SIGAST)
                                val = sigast;
                            else if (set->so.cubs[n].component == SV_OVOLUME_MAT)
                                val = (gdouble)matindex;
                            else if (set->so.cubs[n].component == SV_OVOLUME_MATTYPE)
                                val = (gdouble)mp->mats[matindex].type;
                            else if (set->so.cubs[n].component == SV_OVOLUME_ABS || set->so.cubs[n].component == SV_OVOLUME_SUMALL || set->so.cubs[n].component == SV_OVOLUME_MAXALL)
                                val = mp->out->outabsdata[set->so.cubs[n].k]->data[i][j][k];

                            fwrite(&val, sizeof(gdouble), 1, fh);
                        } // i
                    } // j
                } // k
            } // else
            fclose(fh);
        } // else
    } // n

    for (n = 0; n < set->so.nforces; n++) {
        //FIXME: average over last set->ss.lambda_center/set->sc.dx values;

        if (set->so.forces[i].step == 0)
            continue;

        if (set->sc.step_act % set->so.forces[n].step != 0 && set->sc.step_act < (set->sc.nsteps - 1))
            continue;

        fh = fopen(set->so.forces[n].filename, "w");
        if (fh == NULL) {
            fprintf(stderr, "Error: force: cannot open file %s for writing\n", set->so.forces[i].filename);
            continue;
        }

        for (i = 1; i < set->sc.nsteps; i++) {
            avx = avy = avz = nav = 0;
            for (j = (gint)MAX(0, i - (5.0*set->ss.lambda_center / mp->set->plan.dt / LIGHT_SPEED)); j < i; j++) {
                avx += mp->out->outforcedata[n][3 * j];
                avy += mp->out->outforcedata[n][3 * j + 1];
                avz += mp->out->outforcedata[n][3 * j + 2];
                nav++;
            }

            fprintf(fh, "%d %g %g %g %g %g %g\n", i, mp->out->outforcedata[n][3 * i] * set->sp.dx*set->sp.dy,
                    mp->out->outforcedata[n][3 * i + 1] * set->sp.dx*set->sp.dy,
                    mp->out->outforcedata[n][3 * i + 2] * set->sp.dx*set->sp.dy,
                    avx / nav*set->sp.dx*set->sp.dy,
                    avy / nav*set->sp.dx*set->sp.dy,
                    avz / nav*set->sp.dx*set->sp.dy); //FIXME: here we mix dx and dy - this should be done at accumulation level properly for every side
        }
        fclose(fh);
    }

    for (i = 0; i < set->so.nimgs; i++) {
        if (set->so.imgs[i].step == 0)
            continue;

        if (set->sc.step_act % set->so.imgs[i].step != 0)
            continue;
        else
            outim = TRUE;

        df = gwy_data_field_new(5, 5, 5, 5, 1); /*FIXME, is the duplication necessary?*/
        gwy_data_field_set_si_unit_xy(df, gwy_si_unit_new("pixels"));
        gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("V"));


        if (set->so.imgs[i].component == SV_COMP_EX)
            sv_dcube_get_datafield(mp->d->ex, df, set->so.imgs[i].i, set->so.imgs[i].j, set->so.imgs[i].k);
        else if (set->so.imgs[i].component == SV_COMP_EY)
            sv_dcube_get_datafield(mp->d->ey, df, set->so.imgs[i].i, set->so.imgs[i].j, set->so.imgs[i].k);
        else if (set->so.imgs[i].component == SV_COMP_EZ)
            sv_dcube_get_datafield(mp->d->ez, df, set->so.imgs[i].i, set->so.imgs[i].j, set->so.imgs[i].k);
        else if (set->so.imgs[i].component == SV_COMP_HX)
            sv_dcube_get_datafield(mp->d->hx, df, set->so.imgs[i].i, set->so.imgs[i].j, set->so.imgs[i].k);
        else if (set->so.imgs[i].component == SV_COMP_HY)
            sv_dcube_get_datafield(mp->d->hy, df, set->so.imgs[i].i, set->so.imgs[i].j, set->so.imgs[i].k);
        else if (set->so.imgs[i].component == SV_COMP_HZ)
            sv_dcube_get_datafield(mp->d->hz, df, set->so.imgs[i].i, set->so.imgs[i].j, set->so.imgs[i].k);
        else if (set->so.imgs[i].component == SV_COMP_ALL)
            get_all(mp->d->ex, mp->d->ey, mp->d->ez, df, set->so.imgs[i].i, set->so.imgs[i].j, set->so.imgs[i].k);
         else if (set->so.imgs[i].component == SV_COMP_MAT) {
            if (!mp->d->mat)
                fprintf(stderr, "Array \"material\" does not exist in your present material treatment regime, skipping its output\n");
            else
                sv_icube_get_datafield(mp->d->mat, df, set->so.imgs[i].i, set->so.imgs[i].j, set->so.imgs[i].k);
        } else if (set->so.imgs[i].component == SV_COMP_EPSILON || set->so.imgs[i].component == SV_COMP_SIGMA || set->so.imgs[i].component == SV_COMP_MU || set->so.imgs[i].component == SV_COMP_SIGAST) {
            //TODO, here interpret also tabulated and Drude materials
            if (set->so.imgs[i].component == SV_COMP_EPSILON) {
                if (!mp->d->epsilon)
                    fprintf(stderr, "Array \"epsilon\" does not exist in your present material treatment regime, skipping its output\n");
                else {
                    sv_fcube_get_datafield(mp->d->epsilon, df, set->so.imgs[i].i, set->so.imgs[i].j, set->so.imgs[i].k);
                    gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("F/m"));
                }
            }
            if (set->so.imgs[i].component == SV_COMP_SIGMA) {
                if (!mp->d->sigma)
                    fprintf(stderr, "Array \"sigma\" does not exist in your present material treatment regime, skipping its output\n");
                else {
                    sv_fcube_get_datafield(mp->d->sigma, df, set->so.imgs[i].i, set->so.imgs[i].j, set->so.imgs[i].k);
                    gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("S/m"));
                }
            }
            if (set->so.imgs[i].component == SV_COMP_MU) {
                if (!mp->d->mu)
                    fprintf(stderr, "Array \"mu\" does not exist in your present material treatment regime, skipping its output\n");
                else {
                    sv_fcube_get_datafield(mp->d->mu, df, set->so.imgs[i].i, set->so.imgs[i].j, set->so.imgs[i].k);
                    gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("H/m"));
                }
            }
            if (set->so.imgs[i].component == SV_COMP_SIGAST) {
                if (!mp->d->sigast)
                    fprintf(stderr, "Array \"sigast\" does not exist in your present material treatment regime, skipping its output\n");
                else {
                    sv_fcube_get_datafield(mp->d->sigast, df, set->so.imgs[i].i, set->so.imgs[i].j, set->so.imgs[i].k);
                    gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("H"));
                }
            }
        }

        g_snprintf(key, sizeof(key), "/%u/data", mp->out->last);
        g_snprintf(titlekey, sizeof(titlekey), "/%u/data/title", mp->out->last++);
        g_snprintf(description, sizeof(description), "%s step %d, %s", set->so.imgs[i].filebase, set->sc.step_act, component_string[set->so.imgs[i].component]);

        gwy_container_set_object_by_name(mp->out->data, g_strdup(key), df);
        gwy_container_set_string_by_name(mp->out->data, g_strdup(titlekey), (guchar *)g_strdup(description));

        g_object_unref(df);
    }

    for (i = 0; i < set->so.nsubgrid_imgs; i++) {
        if (set->so.subgrid_imgs[i].step == 0)
            continue;

        if (set->sc.step_act % set->so.subgrid_imgs[i].step != 0)
            continue;
        else
            outim = TRUE;

        df = gwy_data_field_new(5, 5, 5, 5, 1); /*FIXME, is the duplication necessary?*/
        gwy_data_field_set_si_unit_xy(df, gwy_si_unit_new("pixels"));
        gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("V"));

        nsg = set->so.subgrid_imgs[i].format; //subgrid number

        if (set->so.subgrid_imgs[i].component == SV_COMP_EX)
            sv_dcube_get_datafield(mp->sg->sg[nsg]->d->ex, df, set->so.subgrid_imgs[i].i, set->so.subgrid_imgs[i].j, set->so.subgrid_imgs[i].k);
        else if (set->so.subgrid_imgs[i].component == SV_COMP_EY)
            sv_dcube_get_datafield(mp->sg->sg[nsg]->d->ey, df, set->so.subgrid_imgs[i].i, set->so.subgrid_imgs[i].j, set->so.subgrid_imgs[i].k);
        else if (set->so.subgrid_imgs[i].component == SV_COMP_EZ)
            sv_dcube_get_datafield(mp->sg->sg[nsg]->d->ez, df, set->so.subgrid_imgs[i].i, set->so.subgrid_imgs[i].j, set->so.subgrid_imgs[i].k);
        else if (set->so.subgrid_imgs[i].component == SV_COMP_HX)
            sv_dcube_get_datafield(mp->sg->sg[nsg]->d->hx, df, set->so.subgrid_imgs[i].i, set->so.subgrid_imgs[i].j, set->so.subgrid_imgs[i].k);
        else if (set->so.subgrid_imgs[i].component == SV_COMP_HY)
            sv_dcube_get_datafield(mp->sg->sg[nsg]->d->hy, df, set->so.subgrid_imgs[i].i, set->so.subgrid_imgs[i].j, set->so.subgrid_imgs[i].k);
        else if (set->so.subgrid_imgs[i].component == SV_COMP_HZ)
            sv_dcube_get_datafield(mp->sg->sg[nsg]->d->hz, df, set->so.subgrid_imgs[i].i, set->so.subgrid_imgs[i].j, set->so.subgrid_imgs[i].k);
        else if (set->so.subgrid_imgs[i].component == SV_COMP_MAT) {
            if (!mp->d->mat)
                fprintf(stderr, "Array \"material\" does not exist in your present material treatment regime, skipping its output\n");
            else
                sv_icube_get_datafield(mp->sg->sg[nsg]->d->mat, df, set->so.subgrid_imgs[i].i, set->so.subgrid_imgs[i].j, set->so.subgrid_imgs[i].k);
        } else if (set->so.subgrid_imgs[i].component == SV_COMP_EPSILON || set->so.subgrid_imgs[i].component == SV_COMP_SIGMA || set->so.subgrid_imgs[i].component == SV_COMP_MU || set->so.subgrid_imgs[i].component == SV_COMP_SIGAST) {
            //TODO, here interpret also tabulated and Drude materials
            if (set->so.subgrid_imgs[i].component == SV_COMP_EPSILON && mp->d->epsilon) {
                if (!mp->d->epsilon)
                    fprintf(stderr, "Array \"epsilon\" does not exist in your present material treatment regime, skipping its output\n");
                else {
                    sv_fcube_get_datafield(mp->sg->sg[nsg]->d->epsilon, df, set->so.subgrid_imgs[i].i, set->so.subgrid_imgs[i].j, set->so.subgrid_imgs[i].k);
                    gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("F/m"));
                }
            }
            if (set->so.subgrid_imgs[i].component == SV_COMP_SIGMA && mp->d->sigma) {
                if (!mp->d->sigma)
                    fprintf(stderr, "Array \"sigma\" does not exist in your present material treatment regime, skipping its output\n");
                else {
                    sv_fcube_get_datafield(mp->sg->sg[nsg]->d->sigma, df, set->so.subgrid_imgs[i].i, set->so.subgrid_imgs[i].j, set->so.subgrid_imgs[i].k);
                    gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("S/m"));
                }
            }
            if (set->so.subgrid_imgs[i].component == SV_COMP_MU && mp->d->mu) {
                if (!mp->d->mu)
                    fprintf(stderr, "Array \"mu\" does not exist in your present material treatment regime, skipping its output\n");
                else {
                    sv_fcube_get_datafield(mp->sg->sg[nsg]->d->mu, df, set->so.subgrid_imgs[i].i, set->so.subgrid_imgs[i].j, set->so.subgrid_imgs[i].k);
                    gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("H/m"));
                }
            }
            if (set->so.subgrid_imgs[i].component == SV_COMP_SIGAST && mp->d->sigast) {
                if (!mp->d->sigast)
                    fprintf(stderr, "Array \"sigast\" does not exist in your present material treatment regime, skipping its output\n");
                else {
                    sv_fcube_get_datafield(mp->sg->sg[nsg]->d->sigast, df, set->so.subgrid_imgs[i].i, set->so.subgrid_imgs[i].j, set->so.subgrid_imgs[i].k);
                    gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("H"));
                }
            }
        }

        g_snprintf(key, sizeof(key), "/%u/data", mp->out->last);
        g_snprintf(titlekey, sizeof(titlekey), "/%u/data/title", mp->out->last++);
        g_snprintf(description, sizeof(description), "%s step %d, %s", set->so.subgrid_imgs[i].filebase, set->sc.step_act, component_string[set->so.subgrid_imgs[i].component]);

        gwy_container_set_object_by_name(mp->out->data, g_strdup(key), df);
        gwy_container_set_string_by_name(mp->out->data, g_strdup(titlekey), (guchar *)g_strdup(description));

        g_object_unref(df);
    }


    for (n = 0; n < set->so.nplns; n++) {
        if (set->so.plns[i].step == 0)
            continue;

        if (set->sc.step_act % set->so.plns[n].step != 0 && set->sc.step_act < (set->sc.nsteps - 1))
            continue;
        if (set->sc.step_act <= set->so.plns[n].start || set->sc.step_act > set->so.plns[n].stop)
            continue;

        if (set->so.plns[n].format) { //this means ascii mode here
            g_snprintf(filename, 100, "%s_%.4d", set->so.plns[n].filebase, set->sc.step_act);

            fh = fopen(filename, "w");
            if (fh == NULL) {
                fprintf(stderr, "Error: plane: cannot open file %s for writing\n", filename);
                continue;
            }
        } else {
            g_snprintf(filename, 100, "%s_%.4d.raw", set->so.plns[n].filebase, set->sc.step_act);

            fh = fopen(filename, "wb");
            if (fh == NULL) {
                fprintf(stderr, "Error: plane: cannot open file %s for writing\n", filename);
                continue;
            }
        }

        if (set->so.plns[n].j == -1 && set->so.plns[n].k == -1) {
            for (k = 0; k < mp->d->ex->zres; k++) {
                for (j = 0; j < mp->d->ex->yres; j++) {
                    switch (set->so.plns[n].component) {
                        case SV_COMP_EX:
                            val = mp->d->ex->data[set->so.plns[n].i][j][k];
                            break;

                        case SV_COMP_EY:
                            val = mp->d->ey->data[set->so.plns[n].i][j][k];
                            break;

                        case SV_COMP_EZ:
                            val = mp->d->ez->data[set->so.plns[n].i][j][k];
                            break;

                        case SV_COMP_HX:
                            val = mp->d->hx->data[set->so.plns[n].i][j][k];
                            break;

                        case SV_COMP_HY:
                            val = mp->d->hy->data[set->so.plns[n].i][j][k];
                            break;

                        case SV_COMP_HZ:
                            val = mp->d->hz->data[set->so.plns[n].i][j][k];
                            break;

                        default:
                            fprintf(stderr, "Error: Unsupported output type for plane output\n");
                            break;
                    }

                    if (set->so.plns[n].format)
                        fprintf(fh, "%10.10g ", val);
                    else
                        fwrite(&val, sizeof(gdouble), 1, fh);
                } // j
                if (set->so.plns[n].format)
                    fprintf(fh, "\n");
            } // k
        } // if

        if (set->so.plns[n].i == -1 && set->so.plns[n].k == -1) {
            for (k = 0; k < mp->d->ex->zres; k++) {
                for (i = 0; i < mp->d->ex->xres; i++) {
                    switch (set->so.plns[n].component) {
                        case SV_COMP_EX:
                            val = mp->d->ex->data[i][set->so.plns[n].j][k];
                            break;

                        case SV_COMP_EY:
                            val = mp->d->ey->data[i][set->so.plns[n].j][k];
                            break;

                        case SV_COMP_EZ:
                            val = mp->d->ez->data[i][set->so.plns[n].j][k];
                            break;

                        case SV_COMP_HX:
                            val = mp->d->hx->data[i][set->so.plns[n].j][k];
                            break;

                        case SV_COMP_HY:
                            val = mp->d->hy->data[i][set->so.plns[n].j][k];
                            break;

                        case SV_COMP_HZ:
                            val = mp->d->hz->data[i][set->so.plns[n].j][k];
                            break;

                        default:
                            fprintf(stderr, "Error: Unsupported output type for plane output\n");
                            break;
                    }

                    if (set->so.plns[n].format)
                        fprintf(fh, "%10.10g ", val);
                    else
                        fwrite(&val, sizeof(gdouble), 1, fh);
                } // i
                if (set->so.plns[n].format)
                    fprintf(fh, "\n");
            } // k
        } // if

        if (set->so.plns[n].i == -1 && set->so.plns[n].j == -1) {
            for (j = 0; j < mp->d->ex->yres; j++) {
                for (i = 0; i < mp->d->ex->xres; i++) {
                    switch (set->so.plns[n].component) {
                        case SV_COMP_EX:
                            val = mp->d->ex->data[i][j][set->so.plns[n].k];
                            break;

                        case SV_COMP_EY:
                            val = mp->d->ey->data[i][j][set->so.plns[n].k];
                            break;

                        case SV_COMP_EZ:
                            val = mp->d->ez->data[i][j][set->so.plns[n].k];
                            break;

                        case SV_COMP_HX:
                            val = mp->d->hx->data[i][j][set->so.plns[n].k];
                            break;

                        case SV_COMP_HY:
                            val = mp->d->hy->data[i][j][set->so.plns[n].k];
                            break;

                        case SV_COMP_HZ:
                            val = mp->d->hz->data[i][j][set->so.plns[n].k];
                            break;

                        default:
                            fprintf(stderr, "Error: Unsupported output type for plane output\n");
                            break;
                    }

                    if (set->so.plns[n].format)
                        fprintf(fh, "%10.10g ", val);
                    else
                        fwrite(&val, sizeof(gdouble), 1, fh);
                } // i
                if (set->so.plns[n].format)
                    fprintf(fh, "\n");
            } // j
        } // if

        if (fh)
            fclose(fh);
    }

    if (set->sf.nrs > 0) {
        if (set->sc.step_act % 1 == 0) {
            for (n = 0; n < set->sf.nrs; n++) {
                if (!set->sf.source_filename[n])
                    continue;
                fh = fopen(set->sf.source_filename[n], "w");
                if (fh == NULL) { fprintf(stderr, "Error: farfield: cannot open file %s for writing\n", set->sf.source_filename[n]); continue; }

                for (i = 0; i < gwy_data_line_get_res(mp->farfield->rpoints[n].ex); i++) {
                    fprintf(fh, "%d %g %g %g\n", i, -mp->farfield->rpoints[n].ex->data[i],
                            -mp->farfield->rpoints[n].ey->data[i],
                            -mp->farfield->rpoints[n].ez->data[i]);
                }
                fclose(fh);
            } // n
        } // if
        if (set->sc.step_act == (set->sc.nsteps - 1)) {
            outim = TRUE;
            for (indval = 1; indval <= set->sf.nsets; indval++) { /*sets are counted from 1*/
                if (set->sf.setxres[indval] == 1 || set->sf.setyres[indval] == 1) {
                    xdata = (gdouble *)g_malloc(MAX(set->sf.setxres[indval], set->sf.setyres[indval]) * sizeof(gdouble));
                    ydata = (gdouble *)g_malloc(MAX(set->sf.setxres[indval], set->sf.setyres[indval]) * sizeof(gdouble));
                } else {
                    df = gwy_data_field_new(set->sf.setxres[indval], set->sf.setyres[indval], set->sf.setxres[indval], set->sf.setyres[indval], 1);
                    gwy_data_field_set_si_unit_xy(df, gwy_si_unit_new("pixels"));
                    gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("V"));
                }
                n = 0;
                while (set->sf.individual[n] != indval)
                    n++;

                for (col = 0; col < set->sf.setxres[indval]; col++) {
                    for (row = 0; row < set->sf.setyres[indval]; row++) {
                        val = 0;   /*output only total power at the farfield point*/
                        for (i = 0; i < gwy_data_line_get_res(mp->farfield->rpoints[n].ex); i++) {
                            val += (mp->farfield->rpoints[n].ex->data[i] * mp->farfield->rpoints[n].ex->data[i])
                                + (mp->farfield->rpoints[n].ey->data[i] * mp->farfield->rpoints[n].ey->data[i])
                                + (mp->farfield->rpoints[n].ez->data[i] * mp->farfield->rpoints[n].ez->data[i]);
                        }
                        if (set->sf.setxres[indval] == 1 || set->sf.setyres[indval] == 1) {
                            xdata[n] = n;
                            ydata[n] = val;
                        } else
                            gwy_data_field_set_val(df, col, row, val);
                        n++;
                    } // row
                } // col

                if (set->sf.setxres[indval] == 1 || set->sf.setyres[indval] == 1) {
                    gm = gwy_graph_model_new();
                    gc = gwy_graph_curve_model_new();
                    gwy_graph_curve_model_set_data(gc, xdata, ydata, n);
                    gwy_graph_model_add_curve(gm, gc);

                    g_snprintf(key, sizeof(key), "/0/graph/graph/%u", mp->out->glast);
                    g_snprintf(titlekey, sizeof(titlekey), "/0/graph/graph/%u/title", mp->out->glast++);
                    g_snprintf(description, sizeof(description), "far field set %d", indval);

                    gwy_container_set_object_by_name(mp->out->data, g_strdup(key), gm);
                    gwy_container_set_string_by_name(mp->out->data, g_strdup(titlekey), (guchar *)g_strdup(description));
                    g_object_unref(gm);
                } else {
 
                    g_snprintf(key, sizeof(key), "/%u/data", mp->out->last);
                    g_snprintf(titlekey, sizeof(titlekey), "/%u/data/title", mp->out->last++);
                    g_snprintf(description, sizeof(description), "far field set %d", indval);

                    gwy_container_set_object_by_name(mp->out->data, g_strdup(key), df);
                    gwy_container_set_string_by_name(mp->out->data, g_strdup(titlekey), (guchar *)g_strdup(description));
                    g_object_unref(df);
                 }
           } // indval
        } // if
    } // if


    if (set->spf.nrs > 0) {
        if (set->sc.step_act % 1 == 0) {
            for (n = 0; n < set->spf.nrs; n++) {
                if (!set->spf.source_filename[n])
                    continue;

                fh = fopen(set->spf.source_filename[n], "w");
                if (fh == NULL) {
                    fprintf(stderr, "Error: periodic farfield: cannot open file %s for writing\n", set->spf.source_filename[n]);
                    continue;
                }

                for (i = 0; i < gwy_data_line_get_res(mp->pfarfield->prpoints[n].ex); i++) {
                    fprintf(fh, "%d %g %g %g\n", i, -mp->pfarfield->prpoints[n].ex->data[i],
                            -mp->pfarfield->prpoints[n].ey->data[i],
                            -mp->pfarfield->prpoints[n].ez->data[i]);
                }
                fclose(fh);
            }
        }
        /*at the end, output also area result if there was any requested*/
        if (set->sc.step_act == (set->sc.nsteps - 1)) {
            printf("spf nrs output\n");
            outim = TRUE;
            for (indval = 1; indval <= set->spf.nsets; indval++) { /*sets are counted from 1*/
                if (set->spf.setxres[indval] == 1 || set->spf.setyres[indval] == 1) {
                    xdata = (gdouble *)g_malloc(MAX(set->spf.setxres[indval], set->spf.setyres[indval]) * sizeof(gdouble));
                    ydata = (gdouble *)g_malloc(MAX(set->spf.setxres[indval], set->spf.setyres[indval]) * sizeof(gdouble));
                } else {
                    df = gwy_data_field_new(set->spf.setxres[indval], set->spf.setyres[indval], set->spf.setxres[indval], set->spf.setyres[indval], 1);
                    gwy_data_field_set_si_unit_xy(df, gwy_si_unit_new("pixels"));
                    gwy_data_field_set_si_unit_z(df, gwy_si_unit_new("V"));
                }

                n = 0;
                while (set->spf.individual[n] != indval)
                    n++;

                for (col = 0; col < set->spf.setxres[indval]; col++) {
                    for (row = 0; row < set->spf.setyres[indval]; row++) {
                        val = 0;   /*output only total power at the farfield point*/
                        for (i = 0; i < gwy_data_line_get_res(mp->pfarfield->prpoints[n].ex); i++) {
                            val += (mp->pfarfield->prpoints[n].ex->data[i] * mp->pfarfield->prpoints[n].ex->data[i])
                                + (mp->pfarfield->prpoints[n].ey->data[i] * mp->pfarfield->prpoints[n].ey->data[i])
                                + (mp->pfarfield->prpoints[n].ez->data[i] * mp->pfarfield->prpoints[n].ez->data[i]);
                        }
                        if (set->spf.setxres[indval] == 1 || set->spf.setyres[indval] == 1) {
                            xdata[n] = n;
                            ydata[n] = val;
                        } else
                            gwy_data_field_set_val(df, col, row, val);
                        n++;
                    } // row
                } // col

                if (set->spf.setxres[indval] == 1 || set->spf.setyres[indval] == 1) {
                    gm = gwy_graph_model_new();
                    gc = gwy_graph_curve_model_new();
                    gwy_graph_curve_model_set_data(gc, xdata, ydata, n);
                    gwy_graph_model_add_curve(gm, gc);

                    g_snprintf(key, sizeof(key), "/0/graph/graph/%u", mp->out->glast);
                    g_snprintf(titlekey, sizeof(titlekey), "/0/graph/graph/%u/title", mp->out->glast++);
                    g_snprintf(description, sizeof(description), "far field set %d", indval);

                    gwy_container_set_object_by_name(mp->out->data, g_strdup(key), gm);
                    gwy_container_set_string_by_name(mp->out->data, g_strdup(titlekey), (guchar *)g_strdup(description));
                    g_object_unref(gm);
                } else {
                    g_snprintf(key, sizeof(key), "/%u/data", mp->out->last);
                    g_snprintf(titlekey, sizeof(titlekey), "/%u/data/title", mp->out->last++);
                    g_snprintf(description, sizeof(description), "far field set %d", indval);

                    gwy_container_set_object_by_name(mp->out->data, g_strdup(key), df);
                    gwy_container_set_string_by_name(mp->out->data, g_strdup(titlekey), (guchar *)g_strdup(description));

                    g_object_unref(df);
                }
            } // indval
        } // if
    } // if

    if (outim) {
        /*save container*/
        fh = fopen(set->so.outfile, "wb");
        if (fh == NULL)
            fprintf(stderr, "Error opening file %s for writing.\n", set->so.outfile);
        else {
            buffer = gwy_serializable_serialize(G_OBJECT(mp->out->data), NULL);

            if (fwrite(MAGIC2, 1, MAGIC_SIZE, fh) != MAGIC_SIZE || fwrite(buffer->data, 1, buffer->len, fh) != buffer->len)
                fprintf(stderr, "Error saving file\n");

            fclose(fh);
            g_byte_array_free(buffer, TRUE);
        }
    }

    /*for debug purposes only, comment this out*/
    /*if (set->sc.step_act == 300)
      {
      fh = fopen("dump.txt", "w");

      fprintf(fh, "%d %d %d\n", mp->d->ex->xres, mp->d->ex->yres, mp->d->ex->zres);
      for (i=0; i<mp->d->ex->xres; i++)
      {
      for (j=0; j<mp->d->ex->yres; j++)
      {
      for (k=0; k<mp->d->ex->zres; k++)
      {
      fprintf(fh, "%g ", mp->d->ex->data[i][j][k]);

      }
      }
      }

      fclose(fh);
      }
    */
    /*end of debug output*/

    if (set->sc.verbose > 1)
        printf("done.\n");
}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
