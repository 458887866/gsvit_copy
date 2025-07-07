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


/*  settings.h : 
 *  all the main structures that are used for driving computation, 
 *  generally data loaded from different input files
 */

#ifndef SV_SET
#define SV_SET

#include <glib.h>
#include <stdio.h>

#define MAX_GPUS 16

typedef enum {
    SV_BOUNDARY_NONE = 0,
    SV_BOUNDARY_PEC = 1,
    SV_BOUNDARY_LIAO = 2,
    SV_BOUNDARY_CPML = 3,
    SV_BOUNDARY_PERIODIC = 4
} SvBoundaryType;

typedef enum {
    SV_COMP_EX = 0,
    SV_COMP_EY = 1,
    SV_COMP_EZ = 2,
    SV_COMP_HX = 3,
    SV_COMP_HY = 4,
    SV_COMP_HZ = 5,
    SV_COMP_ALL = 6,
    SV_COMP_CUR = 7,
    SV_COMP_EPSILON = 8,
    SV_COMP_SIGMA = 9,
    SV_COMP_MU = 10,
    SV_COMP_SIGAST = 11
} SvCompType;

typedef enum {
    SV_SUM_EX = 0,
    SV_SUM_EY = 1,
    SV_SUM_EZ = 2,
    SV_SUM_ALL = 3,
    SV_SUM_ABS = 4
} SvSumType;

typedef enum {
    SV_MAT_LINEAR = 0,
    SV_MAT_LINTAB = 1,
    SV_MAT_DRUDE = 2,
    SV_MAT_CP3 = 3
} SvMatType;

typedef enum {
    SV_MATMODE_NONE = 0,
    SV_MATMODE_ELECTRIC = 1,
    SV_MATMODE_MAGNETIC = 2,
    SV_MATMODE_FULL = 3,
    SV_MATMODE_LIST = 4
} SvMatMode;

typedef enum {
    SV_GPUMODE_NONE = 0,
    SV_GPUMODE_FULL = 1
} SvGpuMode;

typedef struct
{
    gint i;
    gint j;
    gint component;
    gint skip;
    gint start;
    gint stop;
    gint logscale;
    gchar filebase[50];
} SvOutputPar;

typedef struct
{
    gint i;
    gint j;
    gint k;
    gchar *filename;

    gdouble theta;    /*input wave parameters*/
    gdouble phi;
    gint mode;        /*only for xsvit use*/
    gdouble source_wl;       /*only for xsvit use*/
    gdouble length;   /*only for xsvit use*/
    gdouble source_amplitude;/*only for xsvit use*/

} SvSrcPoint;

typedef struct
{
    gint i0;
    gint j0;
    gint i1;
    gint j1;
    gdouble theta;
    gdouble phi;
    gdouble psi;
    gdouble epsilon;
    gdouble mu;
    gdouble sigma;
    gdouble sigast;
    gchar *filename;
    gint mode;        /*only for xsvit use*/
    gdouble source_wl;       /*only for xsvit use*/
    gdouble length;   /*only for xsvit use*/
    gdouble source_amplitude;/*only for xsvit use*/

} SvSrcTSF;

typedef struct
{
    SvSrcPoint *pnts;
    SvSrcTSF tsf;
    gint npnts;    
} SvSetSource;

typedef struct
{
    SvBoundaryType bx0;
    SvBoundaryType bxn;
    SvBoundaryType by0;
    SvBoundaryType byn;
} SvSetBoundary;

typedef struct
{
    gint bx0pos;
    gint bxnpos;
    gint by0pos;
    gint bynpos;
    gint bz0pos;
    gint bznpos;
    SvBoundaryType bx0;
    SvBoundaryType bxn;
    SvBoundaryType by0;
    SvBoundaryType byn;
} SvSetMBoundary;

typedef struct
{
    gint component;
    gint i0;
    gint j0;
    gint i1;
    gint j1;
    gint skip;
    gdouble epsilon;
    gdouble mu;
    gdouble sigma;
    gdouble sigast;
    gchar filename[50];
} SvOutputSum;

typedef struct
{
    SvOutputPar *pnts;
    SvOutputPar *plns;
    SvOutputPar *imgs;
    SvOutputPar *cubs;
    SvOutputSum *sums;
    gint npnts;
    gint nplns;
    gint nimgs;
    gint ncubs;
    gint nsums;
    gchar *outfile;
} SvSetOutput;

typedef struct {
    gdouble dx;        /*diferences at real space*/
    gdouble dy;
    gint xres;         /*pool resolution in pixels*/
    gint yres;
} SvSetPool;

typedef struct {
    int nsteps;                  /*total number of steps*/ 
    int step_act;                /*actual time step*/
    int verbose;                 /*verbosity level*/
    int ugpu[MAX_GPUS];          /*use concrete GPUs*/
    int ngpu;                    /*maximal number of GPUs to use*/
    int nthreads;
    int suffix;                  /*suffix for temporary files*/
    int devicequery;             /*run devicequery for cards*/
} SvSetComp;

typedef struct {
    gchar* in_material;           /*complete material information*/
    gchar* in_vector;             /*script material information*/
    gboolean matmode_check;       /*whether to check material file actively*/
} SvSetMat;

typedef struct {
    gchar* savefile;              /*file to save or load nfff data*/
    gint i0;                      /*boundaries*/
    gint i1;
    gint j0;
    gint j1;
    gint save;                   /*1 save, 0 load (boundaries are then loaded)*/
    gint ramahi;                 /*1 saving ezp, ezm as hxp, hxm, otherwise all as expected*/
    gint nrs;
} SvSetFarfield;

typedef struct {
    gdouble dt;                  /*time difference*/
    SvMatMode matmode;           /*material mode*/
    SvGpuMode gpumode;           /*gpu mode*/
} SvPlan;

typedef struct _SvSet {
    SvSetPool sp;       /*data pool settings*/
    SvSetComp sc;       /*computation settings*/
    SvSetOutput so;     /*output settings*/
    SvSetSource ss;     /*source settings*/
    SvSetBoundary sb;   /*boundary settings*/
    SvSetMBoundary smb; /*in medium boundary settings*/
    SvSetMat sm;        /*material settings*/
    SvSetFarfield sf;   /*nfff settings*/
    SvPlan plan;        /*computation plan*/
    gint tmmode;        /*tm mode (1), te mode (0)*/
} SvSet;



int clear_settings(SvSet *set);

int parse_settings(gchar *filename, SvSet *set);

void remove_temporary_files(SvSet *set);

gint get_double(FILE *F, gdouble *val, gchar *key);

gint get_int(FILE *F, gint *val, gchar *key);

gchar* get_self_dir();

#endif

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */


