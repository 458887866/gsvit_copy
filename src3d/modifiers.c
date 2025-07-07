
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


 /*  modifier.c :
  *  algorithms for shape modifiers
  */

#include "pool.h"
#include "settings.h"
#include "constants.h"
#include "modifiers.h"
#include <math.h>
#include <omp.h>
#include <libgwyddion/gwyddion.h>
#include <fftw3.h>

#if defined (_WIN32) && defined (_MSC_VER)
#include "../config_win.h"
#else
#include "../config.h"
#include <sys/time.h>
#endif


//static void
//brugemann(double epsilon1, double sigma1, double epsilon2, double sigma2, double p, double *eepsilon, double *esigma);


void output_vtk(float ***data, int xres, int yres, int zres, char *filename)
{
   FILE *fw;
   int i, j, k;

   fw = fopen(filename, "w");
   if (!fw) {
       fprintf(stderr, "Error: cannot open file %s for writing\n", filename);
   }
   fprintf(fw, "# vtk DataFile Version 2.0\n");
   fprintf(fw, "Volume example\n");
   fprintf(fw, "ASCII\n");
   fprintf(fw, "DATASET STRUCTURED_POINTS\n");
   fprintf(fw, "DIMENSIONS %d %d %d\n", xres, yres, zres);
   fprintf(fw, "ASPECT_RATIO 1 1 1\n");
   fprintf(fw, "ORIGIN 0 0 0\n");
   fprintf(fw, "POINT_DATA %d\n", xres*yres*zres);
   fprintf(fw, "SCALARS volume_scalars double 1\n");
   fprintf(fw, "LOOKUP_TABLE default\n");

   for (k=0; k<zres; k++) {
           for (j=0; j<yres; j++) {
               for (i=0; i<xres; i++) {
                   fprintf(fw, "%g\n", (gdouble)data[i][j][k]);
               }
           }
   }

   fclose(fw);
}





static void 
set_seed(GwyRandGenSet *rnd, gint seed) 
{

#ifdef G_OS_WIN32
    if (seed == -1)
        gwy_rand_gen_set_init(rnd, g_random_int() & 0x7fffffff);
    else
        gwy_rand_gen_set_init(rnd, seed);
#else
    struct timeval time;

    if (seed == -1) {
        gettimeofday(&time, NULL);
        gwy_rand_gen_set_init(rnd, (time.tv_sec * 1000) + (time.tv_usec / 1000));
    } else gwy_rand_gen_set_init(rnd, seed);
#endif

}


void
shift_fields_fill(SvFCube *shift, gint xres, gint yres, gint zres, gdouble sigma, gdouble tau, GwyRandGenSet *rnd)
{
    gint i, j, k, n;
    gdouble x, y, z, r, modulus, phi, sum;
    fftw_plan plan;

    fftw_complex *in;
    fftw_complex *out;

    in = (fftw_complex *)fftw_malloc(xres*yres*zres*sizeof(fftw_complex));
    out = (fftw_complex *)fftw_malloc(xres*yres*zres*sizeof(fftw_complex));


    plan = fftw_plan_dft_3d(xres, yres, zres,
                           in, out,     
                           FFTW_BACKWARD, FFTW_ESTIMATE);

    n = 0;
    for (i=0; i<xres; i++)
    {
        x = (i <= xres/2 ? i : xres-i)/(0.5*xres);
        for (j=0; j<yres; j++)
        {
            y = (j <= yres/2 ? j : yres-j)/(0.5*yres);
            for (k=0; k<zres; k++)
            {
                z = (k <= zres/2 ? k : zres-k)/(0.5*zres);
                r = tau*tau*(x*x + y*y + z*z);  //delit 4?

                modulus = exp(-r);     //t kazdeho posunu
                phi = 2*G_PI*gwy_rand_gen_set_double(rnd, 0);

                in[n][0] = modulus*cos(phi);
                in[n++][1] = modulus*sin(phi);
             }
        }
    }
    fftw_execute(plan);


    n = 0;
    sum = 0;
    for (i=0; i<xres; i++)
    {
        for (j=0; j<yres; j++)
        {
            for (k=0; k<zres; k++)
            {
                shift->data[i][j][k] = (gfloat)out[n++][0];
                sum += shift->data[i][j][k]*shift->data[i][j][k];
             }
        }
    }
    sum = sqrt(sum/(xres*yres*zres));
    sv_fcube_multiply(shift, (gfloat)(sigma/sum));  //sigma kazdeho posunu

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
}


void
modifier_fftshift(SvICube *mat, gint material, gdouble xsigma, gdouble xt, gdouble ysigma, gdouble yt, gdouble zsigma, gdouble zt, gint seed)
{
    SvFCube *xshift, *yshift, *zshift;
    SvICube *matbuf;
    gint i, j, k, ni, nj, nk;

    GwyRandGenSet *rnd = gwy_rand_gen_set_new(1);

    gint xres = mat->xres;
    gint yres = mat->yres;
    gint zres = mat->zres;

    xshift = sv_fcube_new(xres, yres, zres, mat->xreal, mat->yreal, mat->zreal, 0);
    yshift = sv_fcube_new_alike(xshift, 0);
    zshift = sv_fcube_new_alike(xshift, 0);

    matbuf = sv_icube_new_alike(mat, 1);

    set_seed(rnd, seed);

    //fill shift fields
    shift_fields_fill(xshift, xres, yres, zres, xsigma, xt, rnd);
    shift_fields_fill(yshift, xres, yres, zres, ysigma, yt, rnd);
    shift_fields_fill(zshift, xres, yres, zres, zsigma, zt, rnd);


    //perform the shift
    for (i=0; i<xres; i++)
    {
        for (j=0; j<yres; j++)
        {
            for (k=0; k<zres; k++)
            {
                ni = i + (gint)xshift->data[i][j][k];
                nj = j + (gint)yshift->data[i][j][k];
                nk = k + (gint)zshift->data[i][j][k];

                ni = CLAMP(ni, 0, xres-1);
                nj = CLAMP(nj, 0, yres-1);
                nk = CLAMP(nk, 0, zres-1);

                matbuf->data[i][j][k] = mat->data[ni][nj][nk];
            }
        }
    }

    //copy the result
    for (i=0; i<xres; i++)
    {
        for (j=0; j<yres; j++)
        {
            for (k=0; k<zres; k++)
            {
                mat->data[i][j][k] = matbuf->data[i][j][k];

            }
        }
    }

    sv_fcube_free(xshift);
    sv_fcube_free(yshift);
    sv_fcube_free(zshift);

    sv_icube_free(matbuf);
}

SvMatProp *
modifier_fftshift_antialiased(SvICube *mat, SvMatProp *mats, gint *nmat,
   gint material, gdouble xsigma, gdouble xt, gdouble ysigma, gdouble yt, gdouble zsigma, gdouble zt, gint seed, gint aa)
{
    gint i, j, k, bi, bj, bk;
    gint xres, yres, zres;
    gint bxres, byres, bzres;
    gint nmstart;
    gdouble fraction, p1/*, eepsilon, esigma*/;

    xres = mat->xres; 
    yres = mat->yres; 
    zres = mat->zres; 

    bxres = mat->xres*aa;
    byres = mat->yres*aa;
    bzres = mat->zres*aa;

    printf("antialiased modifier. nmat = %d\n", *nmat);

    // create additional 27  materials in the table, asumming that our object is only material 0 and vacuum, no general case so far
    nmstart = *nmat - 1;
    printf("mat nmstart %d: %g %g %g %g\n", nmstart, mats[nmstart].epsilon/EPSILON_0, mats[nmstart].mu/MU_0, mats[nmstart].sigma, mats[nmstart].sigast);
    mats = (SvMatProp *)g_realloc(mats, (*nmat + 27) * sizeof(SvMatProp));
    for (i=0; i<27; i++) {
        mats[nmstart + i + 1].type = mats[nmstart].type;

        p1 = (26-(double)i)/26.0;

        //brugemann(mats[nmstart].epsilon/EPSILON_0, mats[nmstart].sigma, 1, 0, p1, &eepsilon, &esigma);
        //mats[nmstart + i + 1].epsilon = eepsilon*EPSILON_0;
        //mats[nmstart + i + 1].sigma = esigma;

        //mats[nmstart + i + 1].epsilon = 1 + (gdouble)i*(mats[nmstart].epsilon/EPSILON_0 - 1)/26.0;
        //mats[nmstart + i + 1].sigma = (gdouble)i*mats[nmstart].sigma/26.0;
        mats[nmstart + i + 1].mu = 1 + (gdouble)i*(mats[nmstart].mu/MU_0 - 1)/26.0;
        mats[nmstart + i + 1].sigast = (gdouble)i*mats[nmstart].sigast/26.0;

        printf("mat %d  p %g : %g %g %g %g\n", nmstart + i + 1, p1,  mats[nmstart + i + 1].epsilon/EPSILON_0, mats[nmstart + i + 1].mu, mats[nmstart + i + 1].sigma, mats[nmstart + i + 1].sigast);
 
        mats[nmstart + i + 1].mu *= MU_0;

    }
    *nmat += 27;


    // set up upscaled buffer and copy data into it
    SvICube *bmat = sv_icube_new(bxres, byres, bzres, 1, 1, 1, 1);
    for (i=0; i<bxres; i++)
    {
        for (j=0; j<byres; j++)
        {
            for (k=0; k<bzres; k++)
            {
                if (mat->data[i/aa][j/aa][k/aa] != 0) bmat->data[i][j][k] = 1;
            }
        }
    }


    modifier_fftshift(bmat, 0, aa*xsigma, aa*xt, aa*ysigma, aa*yt, aa*zsigma, aa*zt, seed);

    for (i=0; i<xres; i++)
    {
        for (j=0; j<yres; j++)
        {
            for (k=0; k<zres; k++)
            {
                fraction = 0;
                for (bi = 0; bi<aa; bi++)
                {
                   for (bj = 0; bj<aa; bj++)
                   {
                      for (bk = 0; bk<aa; bk++)
                      {
                          fraction += (bmat->data[aa*i+bi][aa*j+bj][aa*k+bk]);
                      }
                   }
                }
                fraction /= aa*aa*aa;
                mat->data[i][j][k] = (gint)((double)(nmstart + 1 + 26*fraction));
            }
        }
    }

    sv_icube_free(bmat);

    return mats;
}
 
/*alters boundaries position by random shifts in 3d*/
void 
modifier_shift(SvICube *mat, gint material, gdouble xsigma, gdouble xt, gdouble ysigma, gdouble yt, gdouble zsigma, gdouble zt, gint seed)
{
    SvFCube *xshift, *yshift, *zshift;
    SvFCube *xnoise, *ynoise, *znoise;
    SvFCube *xkernel, *ykernel, *zkernel;
    SvICube *matbuf;

    GwyRandGenSet *rnd = gwy_rand_gen_set_new(1);

    gint xres = mat->xres;
    gint yres = mat->yres;
    gint zres = mat->zres;

    gint xgres = 25;
    gint ygres = 25;
    gint zgres = 25;

    gdouble rsq, xsum, ysum, zsum, xval, yval, zval;

    gint i, j, k, ki, kj, kk, ni, nj, nk;

    xnoise = sv_fcube_new(xres, yres, zres, mat->xreal, mat->yreal, mat->zreal, 0);
    ynoise = sv_fcube_new_alike(xnoise, 0);
    znoise = sv_fcube_new_alike(xnoise, 0);

    xshift = sv_fcube_new(xres, yres, zres, mat->xreal, mat->yreal, mat->zreal, 0);
    yshift = sv_fcube_new_alike(xshift, 0);
    zshift = sv_fcube_new_alike(xshift, 0);

    xkernel = sv_fcube_new(xgres, ygres, zgres, xgres, ygres, zgres, 0);
    ykernel = sv_fcube_new_alike(xkernel, 0);
    zkernel = sv_fcube_new_alike(xkernel, 0);

    matbuf = sv_icube_new_alike(mat, 1);

    set_seed(rnd, seed);    

    printf("random fill, sigma %g %g %g\n", xsigma, ysigma, zsigma);

    //fill by random numbers
    for (i=0; i<xres; i++) 
    {
        for (j=0; j<yres; j++)
        {
            for (k=0; k<zres; k++)
            {
                xnoise->data[i][j][k] = (gfloat)gwy_rand_gen_set_gaussian(rnd, 0, xsigma); 
                ynoise->data[i][j][k] = (gfloat)gwy_rand_gen_set_gaussian(rnd, 0, ysigma); 
                znoise->data[i][j][k] = (gfloat)gwy_rand_gen_set_gaussian(rnd, 0, zsigma); 
             }
        }
    }

    printf("gaussian kernel T %g %g %g\n", xt, yt, zt);

    //create gaussian arrays
    xsum = ysum = zsum = 0;
    for (i=0; i<xgres; i++) 
    {
        for (j=0; j<ygres; j++)
        {
            for (k=0; k<zgres; k++)
            {
                rsq = ((i-xgres/2)*(i-xgres/2) +  (j-ygres/2)*(j-ygres/2) +  (k-zgres/2)*(k-zgres/2));
                xval = exp(-rsq/xt);
                yval = exp(-rsq/yt);
                zval = exp(-rsq/zt);

                xsum += xval;
                ysum += yval;
                zsum += zval;

                xkernel->data[i][j][k] = (gfloat)xval;
                ykernel->data[i][j][k] = (gfloat)yval;
                zkernel->data[i][j][k] = (gfloat)zval;
            }
        }
    }

  //  output_vtk(xkernel->data, xgres, ygres, zgres, "xkernel.vtk");
  //  output_vtk(ykernel->data, xgres, ygres, zgres, "ykernel.vtk");
  //  output_vtk(zkernel->data, xgres, ygres, zgres, "zkernel.vtk");
  
    printf("convolution\n");

    //convolve manually
    for (i=0; i<xres; i++)
    {
        for (j=0; j<yres; j++)
        {
            for (k=0; k<zres; k++)
            {
                xval = yval = zval = 0;

                for (ki=0; ki<xgres; ki++)
                {
                    for (kj=0; kj<ygres; kj++)
                    {
                        for (kk=0; kk<zgres; kk++)
                        {
                            ni = i + ki - xgres/2;
                            if (ni<0) ni += xres;
                            if (ni>(xres-1)) ni -= xres;

                            nj = j + kj - ygres/2;
                            if (nj<0) nj += yres;
                            if (nj>(yres-1)) nj -= yres;

                            nk = k + kk - zgres/2;
                            if (nk<0) nk += zres;
                            if (nk>(zres-1)) nk -= zres;

                            xval += xnoise->data[ni][nj][nk]*xkernel->data[ki][kj][kk]; 
                            yval += ynoise->data[ni][nj][nk]*ykernel->data[ki][kj][kk]; 
                            zval += znoise->data[ni][nj][nk]*zkernel->data[ki][kj][kk];

                           //printf("xx %g %g   %g\n", xnoise->data[ni][nj][nk], xkernel->data[ki][kj][kk], xnoise->data[ni][nj][nk]*xkernel->data[ki][kj][kk]); 
                        }
                    }
                }
                xshift->data[i][j][k] = (gfloat)(xval/xsum);
                yshift->data[i][j][k] = (gfloat)(yval/ysum);
                zshift->data[i][j][k] = (gfloat)(zval/zsum);

                //printf("%g  %g\n", xval, xnoise->data[i][j][k]);//xshift->data[i][j][k]);
            }
        }
        printf(".");
        fflush(stdout);
    }

    printf("\nout files\n");

 //   output_vtk(xshift->data, xres, yres, zres, "xshift.vtk");
 //   output_vtk(yshift->data, xres, yres, zres, "yshift.vtk");
 //   output_vtk(zshift->data, xres, yres, zres, "zshift.vtk");

    //perform the shift
    for (i=0; i<xres; i++) 
    {
        for (j=0; j<yres; j++)
        {
            for (k=0; k<zres; k++)
            {
                ni = i + (gint)xshift->data[i][j][k];
                nj = j + (gint)yshift->data[i][j][k];
                nk = k + (gint)zshift->data[i][j][k];

                if (ni<0 || nj<0 || nk<0 || ni>=xres || nj>=yres || nk>=zres) {
                    continue;
                }

                matbuf->data[i][j][k] = mat->data[ni][nj][nk];

            }
        }
    }

    //copu the result
    for (i=0; i<xres; i++) 
    {
        for (j=0; j<yres; j++)
        {
            for (k=0; k<zres; k++)
            {
                mat->data[i][j][k] = matbuf->data[i][j][k];

            }
        }
    }

  
    sv_fcube_free(xnoise);
    sv_fcube_free(ynoise);
    sv_fcube_free(znoise);

    sv_fcube_free(xshift);
    sv_fcube_free(yshift);
    sv_fcube_free(zshift);

    sv_fcube_free(xkernel);
    sv_fcube_free(ykernel);
    sv_fcube_free(zkernel); 


}

///////////////// experimental part
/*
#include <complex.h>

static void
brugemann(double epsilon1, double sigma1, double epsilon2, double sigma2, double p, double *eepsilon, double *esigma)
{
    double dp, dm, omega = 1;
    double complex epsA, epsB, eresult, d, s, em, ep, e, pp;
  
    pp = p;
    epsA = epsilon1 + I*sigma1;
    epsB = epsilon2 + I*sigma2;

    d = epsA*(2 - 3*pp) - epsB*(1 - 3*pp);
    s = csqrt(d*d + 8.0*epsA*epsB);

    em = (d - s)*0.25;
    ep = (d + s)*0.25;
    if (cimag(ep) > 0.0) eresult = ep;
    else if (cimag(em) > 0.0) eresult = em;
    else {
        e = epsA*(1 - p) + epsB*p;
        dp = cabs(ep - e);
        dm = cabs(em - e);
        if  (dp<=dm) eresult = ep;
        else eresult = em;
    }

    *eepsilon = creal(eresult);
    *esigma = cimag(eresult);
}
*/


/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
