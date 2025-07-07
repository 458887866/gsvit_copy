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


#include <stdio.h>
#include <glib.h>
#include <libprocess/gwyprocess.h>
#include "plan.h"
#include "pool.h"
#include "settings.h"
#include "constants.h"
#include "omp.h"
#include <math.h>
#include <glib.h>

void sv_pool_allocate_output(SvPool *mp, SvSet *set);

/*check whether there is electric or magnetic material in material file*/
SvMatMode 
test_material_file(gchar *filename)
{
    gboolean is_e = 0;
    gboolean is_m = 0;
    int i, j, xres, yres;
    gfloat value;
    FILE *fr = fopen(filename, "rb");

    if (!fr) return SV_MATMODE_NONE;

    fread(&xres, sizeof(gint), 1, fr);
    fread(&yres, sizeof(gint), 1, fr);
    for (i=0; i<xres; i++)
    {
        for (j=0; j<yres; j++)
        {
                   fread(&value, sizeof(gfloat), 1, fr);
                   if (value!=1) {
                       is_e = 1;
                       break;
                   }
            if (is_e) break;
        }
        if (is_e) break;
    }
    if (!is_e)
    {
        for (i=0; i<xres; i++)
        {
            for (j=0; j<yres; j++)
            {
                    fread(&value, sizeof(gfloat), 1, fr);
                    if (value!=0) {
                        is_e = 1;
                        break;
                    }
                if (is_e) break;
            }
            if (is_e) break;
        }
    }

    for (i=0; i<xres; i++)
    {
        for (j=0; j<yres; j++)
        {
                   fread(&value, sizeof(gfloat), 1, fr);
                   if (value!=1) {
                       is_m = 1;
                       break;
                   }
            if (is_m) break;
        }
        if (is_m) break;
    }
    if (!is_m)
    {
        for (i=0; i<xres; i++)
        {
            for (j=0; j<yres; j++)
            {
                    fread(&value, sizeof(gfloat), 1, fr);
                    if (value!=0) {
                        is_m = 1;
                        break;
                    }
                if (is_m) break;
            }
            if (is_m) break;
        }
    }
    fclose(fr);

    if ((is_e+is_m)==0) return SV_MATMODE_NONE;
    if (is_e==0 && is_m!=0) return SV_MATMODE_MAGNETIC;
    if (is_e!=0 && is_m==0) return SV_MATMODE_ELECTRIC;
    
    return SV_MATMODE_FULL;
}

/*check whether there is electric or magnetic material in material script file*/
SvMatMode 
test_script_file(gchar *filename)
{
    gint type;
    gint mattype;
    float epsilon, mu, sigma, sigast, buf, omegap, nu;
    FILE *fr;
    gboolean is_e = 0;
    gboolean is_m = 0;

    fr = fopen(filename, "r");
    while (fscanf(fr, "%d", &type)!=EOF) //TODO use proper locale handling functions here
    {
        switch (type) {
            case 4:
            fscanf(fr, "%f", &buf);
            fscanf(fr, "%f", &buf);
            fscanf(fr, "%f", &buf);
            fscanf(fr, "%f", &buf);
            fscanf(fr, "%d", &mattype);
            printf("type %d\n", mattype);
            if (mattype==SV_MAT_LINEAR || mattype==SV_MAT_LINTAB) {
                fscanf(fr, "%f", &epsilon);
                fscanf(fr, "%f", &mu);
                fscanf(fr, "%f", &sigma);
                fscanf(fr, "%f", &sigast);
                if (mattype==SV_MAT_LINEAR) { //for lintab don't alloc material fields
                    if (epsilon!=1.0 || sigma>0) is_e = 1;
                    if (mu!=1.0 || sigast>0) is_m = 1;
                }
            } else if (mattype==SV_MAT_DRUDE) {
                fscanf(fr, "%f", &epsilon);
                fscanf(fr, "%f", &omegap);
                fscanf(fr, "%f", &nu);
                /*do nothing, if only drude is here we dont need material fields*/
            } else if (type == SV_MAT_CP3)
            {

            } else {
                fprintf(stderr, "Error: Unknown material type?\n");
                return 0;
            }
            break;
        }
    }

    fclose(fr);

    if ((is_e+is_m)==0) return SV_MATMODE_NONE;
    if (is_e==0 && is_m!=0) return SV_MATMODE_MAGNETIC;
    if (is_e!=0 && is_m==0) return SV_MATMODE_ELECTRIC;

    return SV_MATMODE_NONE;
}

gboolean
load_sources(SvPool *mp, SvSet *set)
{
    FILE *fr;
    gint i, k, ndata, pos;
    gdouble ex, ey, ez, hx, hy, hz;

    if (set->ss.npnts)
        mp->src->sp = (SvSourcePoint *)g_malloc(set->ss.npnts*sizeof(SvSourcePoint));

    for (i=0; i<set->ss.npnts; i++)
    {
        fr = fopen(set->ss.pnts[i].filename, "r");
        fscanf(fr, "%d", &ndata);
        if (ndata<=0) return 1;
        if (ndata<set->sc.nsteps) 
            fprintf(stderr, "Warning: less data (%d) in \"%s\" than scheduled computation steps (%d)\n",
                    ndata,
                    set->ss.pnts[i].filename,
                    set->sc.nsteps);
        mp->src->sp[i].i = set->ss.pnts[i].i;
        mp->src->sp[i].j = set->ss.pnts[i].j;
        mp->src->sp[i].sdata.ndata = ndata;
        mp->src->sp[i].sdata.pos = (gint *)g_malloc(ndata*sizeof(gint));
        mp->src->sp[i].sdata.ex = (gdouble *)g_malloc(ndata*sizeof(gdouble));
        mp->src->sp[i].sdata.ey = (gdouble *)g_malloc(ndata*sizeof(gdouble));
        mp->src->sp[i].sdata.ez = (gdouble *)g_malloc(ndata*sizeof(gdouble));
        mp->src->sp[i].sdata.hx = (gdouble *)g_malloc(ndata*sizeof(gdouble));
        mp->src->sp[i].sdata.hy = (gdouble *)g_malloc(ndata*sizeof(gdouble));
        mp->src->sp[i].sdata.hz = (gdouble *)g_malloc(ndata*sizeof(gdouble));

        for (k=0; k<ndata; k++)
        {
            if (get_int(fr, &pos, "PointSource")) return 1;
            if (get_double(fr, &ex, "PointSource")) return 1;
            if (get_double(fr, &ey, "PointSource")) return 1;
            if (get_double(fr, &ez, "PointSource")) return 1;
            if (get_double(fr, &hx, "PointSource")) return 1;
            if (get_double(fr, &hy, "PointSource")) return 1;
            if (get_double(fr, &hz, "PointSource")) return 1;

            mp->src->sp[i].sdata.pos[k] = pos;
            mp->src->sp[i].sdata.ex[k] = ex;
            mp->src->sp[i].sdata.ey[k] = ey;
            mp->src->sp[i].sdata.ez[k] = ez;
            mp->src->sp[i].sdata.hx[k] = hx;
            mp->src->sp[i].sdata.hy[k] = hy;
            mp->src->sp[i].sdata.hz[k] = hz;
        }

        fclose(fr);
    }

    if ((set->ss.tsf.i0 + set->ss.tsf.j0  
         + set->ss.tsf.i1 + set->ss.tsf.j1)!=0)
    {
        mp->src->tsf = (SvSourceTSF *)g_malloc(sizeof(SvSourceTSF));
        mp->src->tsf->i0 = set->ss.tsf.i0;
        mp->src->tsf->j0 = set->ss.tsf.j0;
        mp->src->tsf->i1 = set->ss.tsf.i1;
        mp->src->tsf->j1 = set->ss.tsf.j1;
        mp->src->tsf->theta = set->ss.tsf.theta;
        mp->src->tsf->phi = set->ss.tsf.phi;
        mp->src->tsf->psi = set->ss.tsf.psi;
        mp->src->tsf->epsilon = set->ss.tsf.epsilon;
        mp->src->tsf->mu = set->ss.tsf.mu;
        mp->src->tsf->sigma = set->ss.tsf.sigma;
        mp->src->tsf->sigast = set->ss.tsf.sigast;

        fr = fopen(set->ss.tsf.filename, "r");
        fscanf(fr, "%d", &ndata);
        if (ndata<=0) return 1;
        if (ndata<set->sc.nsteps) 
            fprintf(stderr, "Warning: less data (%d) in \"%s\" than scheduled computation steps (%d) in TSF source\n",
                    ndata,
                    set->ss.pnts[i].filename,
                    set->sc.nsteps);
        mp->src->tsf->ndata = ndata;
        mp->src->tsf->pos = (gint *)g_malloc(ndata*sizeof(gint));
        mp->src->tsf->e = (gdouble *)g_malloc(ndata*sizeof(gdouble));

        for (k=0; k<ndata; k++)
        {
            if (get_int(fr, &pos, "TSF Source")) return 1;
            if (get_double(fr, &ex, "TSF Source")) return 1;

            mp->src->tsf->pos[k] = pos;
            mp->src->tsf->e[k] = ex;
        }
/*
       mp->src->tsf->corr = 1.0/sqrt(pow(sin(mp->src->tsf->theta), 4)
                                             *(pow(cos(mp->src->tsf->phi), 4)+pow(sin(mp->src->tsf->phi), 4))
                                             + pow(cos(mp->src->tsf->theta), 4));
        

        if (set->sc.verbose>1) printf("Angular correction factor for TSF: %g\n", mp->src->tsf->corr);
*/

        mp->src->tsf->corr = 1; /*FIXME this is strange, needs to be solved*/

        mp->src->tsf->jp = sv_1dpool_new((gint)sqrt(set->sp.xres*set->sp.xres+set->sp.yres*set->sp.yres)+10,
                                       set->sp.dx/mp->src->tsf->corr, set->plan.dt);

        sv_1dpool_set_source_pos(mp->src->tsf->jp, 0); 
        sv_1dpool_set_material(mp->src->tsf->jp, mp->src->tsf->epsilon, mp->src->tsf->mu);

        fclose(fr);
    }
return 0;
}

static gboolean is_in_circle(gdouble i, gdouble j, gdouble pnt1[2], gdouble radius)
{
    if (((i-pnt1[0])*(i-pnt1[0])+(j-pnt1[1])*(j-pnt1[1]))<=(radius*radius))
        return 1;
    else return 0;
}

#if(0)
static gdouble dist(gdouble *a, gdouble *b)
{
    return sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]));
}

static gdouble dot(gdouble *a, gdouble *b, gdouble *c)
{
    return (b[0]-a[0])*(c[0]-b[0]) + (b[1]-a[1])*(c[1]-b[1]);
}

static gdouble ltop(gdouble *a, gdouble *b, gdouble *c)
{
    gdouble a1 = a[0] - b[0];
    gdouble a2 = a[1] - b[1];
    gdouble a3 = a[2] - b[2];
    gdouble b1 = a[0] - c[0];
    gdouble b2 = a[1] - c[1];
    gdouble b3 = a[2] - c[2];

    return mag(a2*b3 - a3*b2, a3*b1 - a1*b3, a1*b2 - a2*b1) / mag(c[0] - b[0], c[1] - b[1], c[2] - b[2]);
}

static gdouble mag(gdouble aa, gdouble ab, gdouble ac)
{
    return sqrt(aa*aa + ab*ab + ac*ac);
}

static gdouble linedist(gdouble *a, gdouble *b, gdouble *c)
{
    gdouble dist, dot1, dot2;

    dist = ltop(a,b,c);

    dot1 = dot(b,c,a);
    if (dot1 > 0) return 1e6;

    dot2 = dot(c,b,a);
    if (dot2 > 0) return 1e6;

    return fabs(dist);
}
#endif

static gdouble dett(gdouble a11, gdouble a12, gdouble a13,
            gdouble a21, gdouble a22, gdouble a23,
            gdouble a31, gdouble a32, gdouble a33)
{
    return a11*a22*a33 - a11*a23*a32 + a12*a23*a31 - a12*a21*a33 + a13*a21*a32 - a13*a22*a31;
}

static gdouble determinant(gdouble *a, gdouble *b, gdouble *c)
{
  return dett(a[0], b[0], c[0], a[1], b[1], c[1], a[2], b[2], c[2]);
}

static gboolean is_in_triangle(gdouble i, gdouble j, gdouble *p1, gdouble *p2, gdouble *p3)
{
    gdouble d0[3], d1[3], d2[3], dn[3];
    gdouble rdet[4];

    /*det 0*/
    d0[0] = p1[0]; d0[1] = p1[1]; d0[2] = 1;
    d1[0] = p2[0]; d1[1] = p2[1]; d1[2] = 1;
    d2[0] = p3[0]; d2[1] = p3[1]; d2[2] = 1;
    dn[0] = i; dn[1] = j; dn[2] = 1;

    rdet[0] = determinant(d0, d1, d2);
    rdet[1] = determinant(dn, d1, d2);
    rdet[2] = determinant(d0, dn, d2);
    rdet[3] = determinant(d0, d1, dn);

    if (rdet[0]==0) {
        fprintf(stderr, "Warning: degenerated triangle\n");
    }
    if ((rdet[0]>=0 && rdet[1]>=0 && rdet[2]>=0 && rdet[3]>=0) ||
        (rdet[0]<=0 && rdet[1]<=0 && rdet[2]<=0 && rdet[3]<=0)) {
        return 1;
    }
    return 0;
}


void scan_point(FILE *fr, gdouble *px, gdouble *py)
{
    float x, y;
    fscanf(fr, "%f", &x);  //TODO this should be locale dependent
    fscanf(fr, "%f", &y);
    *px = x;
    *py = y;
}
void scan_radius(FILE *fr, gdouble *pradius)
{
    float radius;
    fscanf(fr, "%f", &radius); //TODO this should be locale dependent
    *pradius = radius;
}


void sv_pool_add_material(SvPool *mp, SvMatProp *mat)
{
    gint i;
    gboolean new = TRUE;
    gdouble diff = 1e-10;

    /*zero index material it is used to determine whether use tabulated material at given point or no*/
    if (mp->nmat==0) {
        mp->mats = (SvMatProp *)g_realloc(mp->mats, sizeof(SvMatProp));
        mp->mats[mp->nmat].epsilon = 1;
        mp->mats[mp->nmat].mu = 1;
        mp->mats[mp->nmat].sigma = 0;
        mp->mats[mp->nmat].sigast = 0;
        mp->nmat++;
    }

    /*check whether we already have this material*/
    for (i=1; i<mp->nmat; i++)
    {
        if (((mat->type==SV_MAT_LINEAR || mat->type==SV_MAT_LINTAB) 
            && (fabs(mat->epsilon-mp->mats[i].epsilon)<diff))
            ||
            ((mat->type==SV_MAT_DRUDE || mat->type==SV_MAT_CP3)
             && (fabs(mat->epsilon-mp->mats[i].epsilon)<diff))
           ) new=FALSE;
    }

    /*add it if necessary*/
    if (new) {
        mp->mats = (SvMatProp *)g_realloc(mp->mats, (mp->nmat+1)*sizeof(SvMatProp));
        mp->mats[mp->nmat].epsilon = mat->epsilon;
        mp->mats[mp->nmat].mu = mat->mu;
        mp->mats[mp->nmat].sigma = mat->sigma;
        mp->mats[mp->nmat].sigast = mat->sigast;
        mp->mats[mp->nmat].drude_omega_p = mat->drude_omega_p;
        mp->mats[mp->nmat].drude_nu = mat->drude_nu;
        mp->mats[mp->nmat].type = mat->type;
        for (i=0; i<3; i++) 
        {
            mp->mats[mp->nmat].cp3_a[i] = mat->cp3_a[i];
            mp->mats[mp->nmat].cp3_phi[i] = mat->cp3_phi[i];
            mp->mats[mp->nmat].cp3_omega[i] = mat->cp3_omega[i];
            mp->mats[mp->nmat].cp3_gamma[i] = mat->cp3_gamma[i];
        }
        mat->pos = mp->nmat;
        mp->nmat++;

        printf("material added\n");
    }
}

void scan_material(FILE *fr, SvMatProp *mat)
{
    float epsilon, mu, sigma, sigast, omegap, nu, a, omega, phi, gamma;
    int type;

    fscanf(fr, "%d", &type);
    if (type == SV_MAT_LINEAR || type == SV_MAT_LINTAB) {
        fscanf(fr, "%f", &epsilon); //TODO this should be locale dependent
        fscanf(fr, "%f", &mu);
        fscanf(fr, "%f", &sigma);
        fscanf(fr, "%f", &sigast);
        mat->epsilon = epsilon;
        mat->mu = mu;
        mat->sigma = sigma;
        mat->sigast = sigast;
        mat->type = type;
    } else if (type == SV_MAT_DRUDE) {
        fscanf(fr, "%f", &epsilon); //TODO this should be locale dependent
        fscanf(fr, "%f", &omegap);
        fscanf(fr, "%f", &nu);
        mat->epsilon = epsilon;
        mat->drude_omega_p = omegap;
        mat->drude_nu = nu;
        mat->type = type;
    } else if (type == SV_MAT_CP3) {
        mat->type = type;
        fscanf(fr, "%f", &epsilon);
        fscanf(fr, "%f", &sigma);
        fscanf(fr, "%f", &a);
        fscanf(fr, "%f", &phi);
        fscanf(fr, "%f", &omega);
        fscanf(fr, "%f", &gamma);
        mat->epsilon = epsilon;
        mat->sigma = sigma;
        mat->cp3_a[0] = a;
        mat->cp3_phi[0] = phi;
        mat->cp3_omega[0] = omega;
        mat->cp3_gamma[0] = gamma;

        fscanf(fr, "%f", &a);
        fscanf(fr, "%f", &phi);
        fscanf(fr, "%f", &omega);
        fscanf(fr, "%f", &gamma);
        mat->cp3_a[1] = a;
        mat->cp3_phi[1] = phi;
        mat->cp3_omega[1] = omega;
        mat->cp3_gamma[1] = gamma;

        fscanf(fr, "%f", &a);
        fscanf(fr, "%f", &phi);
        fscanf(fr, "%f", &omega);
        fscanf(fr, "%f", &gamma);
        mat->cp3_a[2] = a;
        mat->cp3_phi[2] = phi;
        mat->cp3_omega[2] = omega;
        mat->cp3_gamma[2] = gamma;
        printf("loaded a/gamma[0]: %g %g\n", mat->cp3_a[0], mat->cp3_gamma[0]);
    }
}



gint load_data(SvPool *mp, SvSet *set)
{
    gint i, j, xres, yres;
    guint m;
    gint type;
    SvCircle cx;
    SvRectangle rx;
    SvTriangle tx;
    GArray *circles = NULL;
    GArray *rectangles;
    GArray *triangles;
    gint npos, ngtot = 0;
    FILE *fr = NULL;   

    if(set->sm.in_material != NULL)
        fr = fopen(set->sm.in_material, "rb");

    if (set->plan.matmode != SV_MATMODE_NONE && fr) {
        if (set->sc.verbose) 
            printf("Loading MEDIUM_LINEAR data...\n");

        fread(&xres, sizeof(gint), 1, fr);
        fread(&yres, sizeof(gint), 1, fr);
        if (xres!=set->sp.xres || yres!=set->sp.yres)
        {
            fprintf(stderr, "Error: material dimensions differ from pool settings\n");
            return 1;
        }

        if (set->plan.matmode == SV_MATMODE_MAGNETIC) {
            fseek(fr, 2*xres*yres*sizeof(gfloat), SEEK_SET);
        }
        else {
            fread(mp->epsilon->data, sizeof(gfloat), xres*yres, fr);
            gwy_data_field_multiply(mp->epsilon, EPSILON_0);

            fread(mp->sigma->data, sizeof(gfloat), xres*yres, fr);
        }

        if (set->plan.matmode == SV_MATMODE_ELECTRIC) {
            fclose(fr);
            return 0;
        }
        fread(mp->mu->data, sizeof(gfloat), xres*yres, fr);
        gwy_data_field_multiply(mp->mu, MU_0);

        fread(mp->sigast->data, sizeof(gfloat), xres*yres, fr);

        fclose(fr);
    }

    fr = NULL;
    if (set->sm.in_vector != NULL)
        fr = fopen(set->sm.in_vector, "r");
    if (fr) {
        if (set->sc.verbose) 
            printf("Loading MEDIUM_VECTOR data...\n");
        xres = set->sp.xres;
        yres = set->sp.yres;

        circles = g_array_new (FALSE, FALSE, sizeof (SvCircle));
        rectangles = g_array_new (FALSE, FALSE, sizeof (SvRectangle));
        triangles = g_array_new (FALSE, FALSE, sizeof (SvTriangle));

        while (fscanf(fr, "%d", &type)!=EOF) {
            switch (type) {
                case 2:
                scan_point(fr, cx.pnt1, cx.pnt1+1);
                scan_radius(fr, &cx.radius);
                scan_material(fr, &cx.mat);
                printf("circle %g %g mat type %d\n", cx.pnt1[0], cx.pnt1[1], cx.mat.type);
                if (cx.mat.type!=0) 
                    sv_pool_add_material(mp, &cx.mat);
                cx.n = ngtot++;
                g_array_append_val(circles, cx);
                break;

                case 3:
                scan_point(fr, tx.pnt1, tx.pnt1+1);
                scan_point(fr, tx.pnt2, tx.pnt2+1);
                scan_point(fr, tx.pnt3, tx.pnt3+1);
                scan_material(fr, &tx.mat);
                printf("tx mat type %d\n", tx.mat.type);
                if (tx.mat.type!=0) 
                    sv_pool_add_material(mp, &tx.mat);
                tx.n = ngtot++;
                g_array_append_val(triangles, tx);
                break;

                case 4:
                scan_point(fr, rx.pnt1, rx.pnt1+1);
                scan_point(fr, rx.pnt2, rx.pnt2+1);
                scan_material(fr, &rx.mat);
                printf("rx mat type %d\n", rx.mat.type);
                if (rx.mat.type!=0) 
                    sv_pool_add_material(mp, &rx.mat);
                rx.n = ngtot++;
                g_array_append_val(rectangles, rx);
                break;

                default:
                break;
            }
        }

        printf("done, now process it\n");
        /*now process the arrays*/

        /*alloc tabulated fields if necessary*/
        if (mp->nmat!=0) 
        {
            mp->mat = gwy_data_field_new(set->sp.xres, set->sp.yres,
                                         set->sp.xres*set->sp.dx, 
                                         set->sp.yres*set->sp.dy, 
                                         1); 

            //no recursive accumulator here now
            for (i=0; i<mp->nmat; i++) {
                if (mp->mats[i].type==SV_MAT_DRUDE) {
                    mp->mats[i].drude_omega_p *= EV_TO_J/LIGHT_SPEED;
                    mp->mats[i].drude_nu *= EV_TO_J/LIGHT_SPEED;
                } else if (mp->mats[i].type==SV_MAT_LINTAB) {
                    mp->mats[i].epsilon *= EPSILON_0;
                    mp->mats[i].mu *= MU_0; 
                } else if (mp->mats[i].type==SV_MAT_CP3) {
                }
            }
        }

        npos = 0;

        do {
            for (m=0; m<circles->len; m++) {        
                cx = g_array_index(circles, SvCircle, m);

                if (npos != cx.n) 
                    continue;
                else 
                    npos++;

                for (i=(gint)MAX(0, cx.pnt1[0]-cx.radius-1); i<(gint)MIN(xres, cx.pnt1[0]+cx.radius+1); i++) { //x
                    for (j= (gint)MAX(0, cx.pnt1[1]-cx.radius-1); j<(gint)MIN(yres, cx.pnt1[1]+cx.radius+1); j++) { //y
                        if (is_in_circle(i, j, cx.pnt1, cx.radius)) {
                            if (cx.mat.type==0) {
                                if (set->plan.matmode == SV_MATMODE_ELECTRIC || set->plan.matmode == SV_MATMODE_FULL) {
                                    mp->epsilon->data[j*xres+i] = cx.mat.epsilon*EPSILON_0;
                                    mp->sigma->data[j*xres+i] = cx.mat.sigma;
                                } else if (set->plan.matmode == SV_MATMODE_MAGNETIC || set->plan.matmode == SV_MATMODE_FULL) {
                                    mp->mu->data[j*xres+i] = cx.mat.mu*MU_0;
                                    mp->sigast->data[j*xres+i] = cx.mat.sigast;
                                }
                            } else {
                                mp->mat->data[j*xres + i] = cx.mat.pos;
                            }
                        }
                    }
                }
            }

            for (m=0; m<triangles->len; m++) {        
                tx = g_array_index(triangles, SvTriangle, m);

                if (npos != tx.n) 
                    continue;
                else 
                    npos++;

                for (i=0; i<xres; i++) { //x
                    for (j=0; j<yres; j++) { //y
                        if (is_in_triangle(i, j, tx.pnt1, tx.pnt2, tx.pnt3)) {
                            if (tx.mat.type==0) {
                                if (set->plan.matmode == SV_MATMODE_ELECTRIC || set->plan.matmode == SV_MATMODE_FULL) {
                                    mp->epsilon->data[j*xres + i] = tx.mat.epsilon*EPSILON_0;
                                    mp->sigma->data[j*xres + i] = tx.mat.sigma;
                                } else if (set->plan.matmode == SV_MATMODE_MAGNETIC || set->plan.matmode == SV_MATMODE_FULL) {
                                    mp->mu->data[j*xres + i] = tx.mat.mu*MU_0;
                                    mp->sigast->data[j*xres + i] = tx.mat.sigast;
                                }
                            } else {
                                mp->mat->data[j*xres + i] = tx.mat.pos;
                            }

                        }
                    }
                }
            }
            for (m=0; m<rectangles->len; m++) {        
                rx = g_array_index(rectangles, SvRectangle, m);

                if (npos != rx.n) 
                    continue;
                else 
                    npos++;

                for (i= (gint)rx.pnt1[0]; i<rx.pnt2[0]; i++) { //x
                    for (j= (gint)rx.pnt1[1]; j<rx.pnt2[1]; j++) { //y
                        if (rx.mat.type==0) {
                            if (set->plan.matmode == SV_MATMODE_ELECTRIC || set->plan.matmode == SV_MATMODE_FULL) {
                                mp->epsilon->data[j*xres + i] = rx.mat.epsilon*EPSILON_0;
                                mp->sigma->data[j*xres + i] = rx.mat.sigma;
                            } else if (set->plan.matmode == SV_MATMODE_MAGNETIC || set->plan.matmode == SV_MATMODE_FULL) {
                                mp->mu->data[j*xres + i] = rx.mat.mu*MU_0;
                                mp->sigast->data[j*xres + i] = rx.mat.sigast;
                            }
                        } else {
                            mp->mat->data[j*xres + i] = rx.mat.pos;
                        }
                    }
                }
            }
        } while (npos<ngtot);
    }
    g_array_free(circles, TRUE);

    return 0;
}

SvPool* init_and_make_plan(SvSet *set)
{
    SvPool *mp;
    gdouble alloc, ram, shared;

    gint matmode_mat, matmode_script;

    set->plan.gpumode = SV_GPUMODE_NONE;

    /*determine matmode - whether to alloc all material fields*/
    /*first test user options - none at present*/
    /*then read material and script file if any and check for materials*/
    if (set->sm.matmode_check) { //if not, material was set by user already
        matmode_mat = SV_MATMODE_NONE;
        matmode_script = SV_MATMODE_NONE;
        printf("matmodes: %d %d\n", matmode_mat, matmode_script);

        if (set->sm.in_material != NULL) 
            matmode_mat = test_material_file(set->sm.in_material);
        if (set->sm.in_vector != NULL) 
            matmode_script = test_script_file(set->sm.in_vector);
        printf("matmodes: %d %d\n", matmode_mat, matmode_script);

        if (matmode_mat==SV_MATMODE_NONE && matmode_script==SV_MATMODE_NONE) 
            set->plan.matmode = SV_MATMODE_NONE;
        else if (matmode_mat==SV_MATMODE_ELECTRIC && (matmode_script==SV_MATMODE_NONE || matmode_script==SV_MATMODE_ELECTRIC))
            set->plan.matmode = SV_MATMODE_ELECTRIC;
        else if (matmode_mat==SV_MATMODE_MAGNETIC && (matmode_script==SV_MATMODE_NONE || matmode_script==SV_MATMODE_MAGNETIC))
            set->plan.matmode = SV_MATMODE_MAGNETIC;
        else if (matmode_script==SV_MATMODE_ELECTRIC && (matmode_mat==SV_MATMODE_NONE || matmode_mat==SV_MATMODE_ELECTRIC))
            set->plan.matmode = SV_MATMODE_ELECTRIC;
        else if (matmode_script==SV_MATMODE_MAGNETIC && (matmode_mat==SV_MATMODE_NONE || matmode_mat==SV_MATMODE_MAGNETIC))
            set->plan.matmode = SV_MATMODE_MAGNETIC;
         else set->plan.matmode = SV_MATMODE_FULL;
    } else 
        set->plan.matmode = SV_MATMODE_FULL;

    /*assume what memory will be allocated*/
    alloc = set->sp.xres*set->sp.yres;
    if (set->plan.matmode == SV_MATMODE_FULL) 
        ram = alloc*10*sizeof(gdouble);
    else if (set->plan.matmode == SV_MATMODE_ELECTRIC || set->plan.matmode == SV_MATMODE_MAGNETIC) 
        ram = alloc*8*sizeof(gdouble);
    else 
        ram = alloc*6*sizeof(gdouble);
    if (set->sc.verbose) 
        printf("RAM allocated memory: %g MB\n", ram/1e6);

    if (!set->plan.gpumode == SV_GPUMODE_NONE) {
        if (set->plan.matmode == SV_MATMODE_FULL) 
            shared = alloc*10*sizeof(float);
        else if (set->plan.matmode == SV_MATMODE_ELECTRIC || set->plan.matmode == SV_MATMODE_MAGNETIC) 
            shared = alloc*8*sizeof(float);
        else 
            shared = alloc*6*sizeof(float);
        if (set->sc.verbose)
            printf("GPU allocated memory: %g MB\n", shared/1e6);
    }

    /*alloc main fields*/
    set->plan.dt = 1/LIGHT_SPEED/sqrt(1/set->sp.dx/set->sp.dx+1/set->sp.dy/set->sp.dy);
    
    if (set->sc.verbose) 
        printf("Initializing main fields...\n");
    mp = sv_pool_new(set);

    /*load main fields from data files and allocate what is necessary*/
    if (set->sm.in_material != NULL || set->sm.in_vector != NULL) {
       if (set->sc.verbose) 
           printf("Loading material data...\n");
       load_data(mp, set); 
    }

    /*load sources and allocate what is necessary*/
    if (set->sc.verbose) 
        printf("Loading sources...\n");
    load_sources(mp, set);

    /*allocate GPUs depending on what is in sources and data*/
    if (set->sc.verbose) 
        printf("Allocate GPU fields...\n");
    sv_pool_allocate_gpus(mp, set);

    /*allocate CPU outputs (gpu output data allocated in previous step)*/
    if (set->sc.verbose) 
        printf("Allocate outputs...\n");
    sv_pool_allocate_output(mp, set);


    if (set->sc.nthreads == -1) 
        set->sc.nthreads = omp_get_max_threads(); //omp does what it wants now, ignoring this and running all the available threads
    else 
        omp_set_num_threads(set->sc.nthreads);


    if (set->sc.verbose) {
        printf("Computation plan:\n");
        if (set->plan.matmode == SV_MATMODE_FULL) 
            printf("Material properties: full\n");
        else if (set->plan.matmode == SV_MATMODE_ELECTRIC) 
            printf("Material properties: electric material only\n");
        else if (set->plan.matmode == SV_MATMODE_MAGNETIC) 
            printf("Material properties: mangetic material only\n");
        else if (set->plan.matmode == SV_MATMODE_NONE && mp->nmat!=0) 
            printf("Material properties: tabulated\n");
        else 
            printf("Material properties: none\n");

        if (set->plan.gpumode == SV_GPUMODE_NONE) 
            printf("GPU use: none\n");
        else 
            printf("GPU use: full computation\n");
        printf("Timestep: %g s\n", set->plan.dt);
    }

    return mp;
}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */


