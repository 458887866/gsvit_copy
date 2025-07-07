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

#ifndef SETTINGS_H
#define SETTINGS_H

#include <glib.h>
#include <stdio.h>
#include "defaults.h"
#include "dcube.h"

#define MAX_GPUS    16
#define NUM_LSMP    50   /* layered source material properties */

/*boundary types*/
typedef enum {
    SV_BOUNDARY_NONE = 0,
    SV_BOUNDARY_PEC = 1,
    SV_BOUNDARY_LIAO = 2,
    SV_BOUNDARY_CPML = 3,
    SV_BOUNDARY_PERIODIC = 4
} SvBoundaryType;

/*field component for output*/
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
    SV_COMP_SIGAST = 11,
    SV_COMP_MAT = 12,
    SV_COMP_EMAX = 13,
    SV_COMP_EXMAX = 14,
    SV_COMP_EYMAX = 15,
    SV_COMP_EZMAX = 16,
    SV_COMP_HMAX = 17,
    SV_COMP_HXMAX = 18,
    SV_COMP_HYMAX = 19,
    SV_COMP_HZMAX = 20,
    SV_COMP_ALLFFT = 21
} SvCompType;

/*volume output*/
typedef enum {
    SV_OVOLUME_EX = 0,
    SV_OVOLUME_EY = 1,
    SV_OVOLUME_EZ = 2,
    SV_OVOLUME_HX = 3,
    SV_OVOLUME_HY = 4,
    SV_OVOLUME_HZ = 5,
    SV_OVOLUME_ALL = 6,
    SV_OVOLUME_EPSILON = 7,
    SV_OVOLUME_SIGMA = 8,
    SV_OVOLUME_MU = 9,
    SV_OVOLUME_SIGAST = 10,
    SV_OVOLUME_MAT = 11,
    SV_OVOLUME_MATTYPE = 12,
    SV_OVOLUME_ABS = 13,
    SV_OVOLUME_SUMALL = 14,
    SV_OVOLUME_MAXALL = 15
} SvOutputVolumeType;

/*summation for output*/
typedef enum {
    SV_SUM_EX = 0,
    SV_SUM_EY = 1,
    SV_SUM_EZ = 2,
    SV_SUM_ALL = 3,
    SV_SUM_ABS = 4,
    SV_SUM_MAX = 5
} SvSumType;

/*material type for input and further treatment*/
typedef enum {
    SV_MAT_LINEAR = 0, /*epsilon, sigma, mu, sigast*/
    SV_MAT_LINTAB = 1, /*the same, but given by table and index*/
    SV_MAT_DRUDE = 2,  /*rc drude model*/
    SV_MAT_CP3 = 3,    /*rc cp3 model*/
    SV_MAT_CP = 4,     /*rc cp model*/
    SV_MAT_ADE = 5,    /*drude-cp model for Auxiliary Differential Equations*/
    SV_MAT_PLRC = 6,   /*plrc drude-cp model*/
    SV_MAT_PEC = 10,   /*pec*/
    SV_MAT_DATABASE = 99 /*linear data loaded from database prior to computation*/
} SvMatType;

/*allocation and computation mode*/
typedef enum {
    SV_MATMODE_NONE = 0,
    SV_MATMODE_ELECTRIC = 1,
    SV_MATMODE_MAGNETIC = 2,
    SV_MATMODE_FULL = 3,
    SV_MATMODE_LIST = 4,
    SV_MATMODE_DATABASE = 99,
    SV_MATMODE_ERROR = 999
} SvMatMode;

/*use of GPU*/
typedef enum {
    SV_GPUMODE_NONE = 0,
    SV_GPUMODE_FULL = 1
} SvGpuMode;

/*leaf type - paramaters*/
typedef enum {    
    SV_TYPE_POOL = 0,
    SV_TYPE_BASIC = 1,
    SV_SRCTYPE_SF = 2,
    SV_SRCTYPE_TSF = 3,
    SV_SRCTYPE_TSFF = 4,
    SV_SRCTYPE_LTSF = 5,
    SV_SRCTYPE_LTSFF = 6,    
    SV_SRCTYPE_POINT = 7,
    SV_SRCTYPE_EXT = 8,
    SV_OUTTYPE_GENERAL = 9,
    SV_OUTTYPE_POINT = 10,
    SV_OUTTYPE_IMAGE = 11,
    SV_OUTTYPE_PLANE = 12,
    SV_OUTTYPE_VOLUME = 13,
    SV_OUTTYPE_SUM = 14,
    SV_OUTTYPE_SUMTAB = 15,
    SV_OUTTYPE_FORCE = 16,
    SV_TYPE_BND = 17,
    SV_TYPE_MATERIAL_PROP = 18,
    SV_TYPE_MATERIAL_GROW = 19,
    SV_TYPE_MATERIAL_ROUGHNESS = 20,
    SV_TYPE_MATERIAL_SPECTRAL = 21,
    SV_TYPE_MATERIAL_EXPRESSION = 22,
    SV_TYPE_NFFF_BOX = 23,
    SV_TYPE_NFFF_POINT = 24,
    SV_TYPE_NFFF_AREA = 25,
    SV_TYPE_PNFFF_BOX = 26,
    SV_TYPE_PNFFF_POINT = 27,
    SV_TYPE_PNFFF_AREA = 28,
    SV_OUTTYPE_SUBGRID_IMAGE = 29,
} SvType;

/*leaf type-material object*/
typedef enum {
    SV_TYPE_MAT_SPHERE = 0,
    SV_TYPE_MAT_VOXEL = 1,
    SV_TYPE_MAT_CYLINDER = 2,
    SV_TYPE_MAT_CONE = 3,
    SV_TYPE_MAT_RCONE = 4,
    SV_TYPE_MAT_GWYDD = 5,
    SV_TYPE_MAT_MESH = 6,
    SV_TYPE_MAT_TETRAHEDRON = 7,
} SvTypeMat;

/*general output parameters, not all of them used for every case and some have various meaning*/
typedef struct
{
    gint i;
    gint j;
    gint k;
    gint component; /*see SvCompType*/
    gint step;      /*output skip*/
    gint start;     /*output start*/
    gint stop;      /*output stop*/
    gint format;    /*used in plane and volume output it selects between binary and ascii mode. In subgrid it idendifies the subgrid number*/
    gchar filebase[100];
} SvOutputPar;

/*point source*/
typedef struct
{
    /* 0b. Point source origin properties */
    gint        point_origin_position_i;
    gint        point_origin_position_j;
    gint        point_origin_position_k;
    gdouble     point_origin_theta;          /*input origin parameters*/
    gdouble     point_origin_phi;

    /* 1. Source properties */
    gint        source_mode;        /*only for xsvit use*/
    gchar       *source_filename;
    gdouble     source_wl;          /*only for xsvit use*/
    gdouble     source_wlspan;      /*only for xsvit use*/
    gdouble     source_amplitude;   /*only for xsvit use*/
    gdouble     source_pulsewidth;  /*only for xsvit use*/
} SvSrcPoint;

/*externally prepared source (defined by E- and H-field dependencies on a plane*/
typedef struct
{
    gint i;           /*plane orientation and position (same as for plane output)*/
    gint j;
    gint k;
    gint ijstart;     /*plane position offset within the main grid*/
    gint jkstart;
    gint iextfrom;    /*subset selection from available data in external source*/
    gint jextfrom;
    gint iextto;
    gint jextto;
    gint shift;       /*time shift with respect to source data (in steps) */
    gint extxres;     /*x/y resolution of external source data*/
    gint extyres;
    gchar *filebase_ex;
    gchar *filebase_ey;
    gchar *filebase_ez;
    gchar *filebase_hx;
    gchar *filebase_hy;
    gchar *filebase_hz;

    gboolean valid;   /*TODO: implement in gsvit*/
} SvSrcExt;

/* Scattered field source (SF) */
typedef struct
{
    /* 1. Source properties */
    gint        source_mode;            /*only for xsvit use*/
    gchar       *source_filename;
    gdouble     source_wl;              /*only for xsvit use*/
    gdouble     source_wlspan;          /*only for xsvit use*/
    gdouble     source_amplitude;       /*only for xsvit use*/
    gdouble     source_pulsewidth;      /*only for xsvit use*/

    /* 2. Incident angle */
    gdouble     ia_theta;               /*input wave parameters [rad]*/
    gdouble     ia_phi;
    gdouble     ia_psi;

    gdouble     layered_epsilon;        /*material properties for whole TSF*/
    gdouble     layered_mu;
    gdouble     layered_sigma;
    gdouble     layered_sigast;

    gboolean    valid;   /*TODO: implement in gsvit*/
} SvSrcSF;

/* Total/scattered field source (TSF) */
typedef struct
{
    /* 0a. Box properties */
    gint        box_i0;          
    gint        box_j0;
    gint        box_k0;
    gint        box_in;
    gint        box_jn;
    gint        box_kn;
    gint        box_boundary_skipi0;    /*skipping of any of the volume faces*/
    gint        box_boundary_skipj0;
    gint        box_boundary_skipk0;
    gint        box_boundary_skipin;
    gint        box_boundary_skipjn;
    gint        box_boundary_skipkn;
    gint        box_boundary_skipdepth; /*skipping by material neighbourhood*/

    /* 1. Source properties */
    gint        source_mode;        /*only for xsvit use*/
    gchar       *source_filename;
    gdouble     source_wl;          /*only for xsvit use*/
    gdouble     source_wlspan;      /*only for xsvit use*/
    gdouble     source_amplitude;   /*only for xsvit use*/
    gint        source_pulsewidth;  /*only for xsvit use*/

    /* 2. Incident angle */
    gdouble     ia_theta;           /*input wave parameters [rad]*/
    gdouble     ia_phi;
    gdouble     ia_psi;
        
    gdouble     layered_epsilon;    /*material parameters for whole TSF*/
    gdouble     layered_mu;
    gdouble     layered_sigma;
    gdouble     layered_sigast;
       
    gint        gaussian;
    gdouble     gaussian_fxpos;
    gdouble     gaussian_fypos;
    gdouble     gaussian_rx;
    gdouble     gaussian_ry;

    gint        radial;
    gdouble     radial_fxpos;
    gdouble     radial_fypos;
    gdouble     radial_rx;
    gdouble     radial_ry;

    gint        fiber;    
    gdouble     fiber_fxpos;
    gdouble     fiber_fypos;
    gdouble     fiber_radius;
    gdouble     fiber_cutoff;
    gdouble     fiber_epsilon_core;
    gdouble     fiber_epsilon_cladding;

    gboolean    valid;              /*TODO: implement in gsvit*/
} SvSrcTSF;

/* Focused Total/Scattered field source (TSFF) */
typedef struct
{
    /* 0a. Box properties */
    gint box_i0;
    gint box_j0;
    gint box_k0;
    gint box_in;
    gint box_jn;
    gint box_kn;
    gint box_boundary_skipi0;
    gint box_boundary_skipin;
    gint box_boundary_skipj0;
    gint box_boundary_skipjn;
    gint box_boundary_skipk0;
    gint box_boundary_skipkn;
    gint box_boundary_skipdepth;    /*skipping by material neighbourhood*/

    /* 1. Source properties */
    gint    source_mode;            /*only for xsvit use*/
    gchar   *source_filename;
    gdouble source_wl;              /*only for xsvit use*/
    gdouble source_wlspan;          /*only for xsvit use*/
    gdouble source_amplitude;       /*only for xsvit use*/
    gdouble source_pulsewidth;      /*only for xsvit use*/

    /* 3. Focused source properties */
    gdouble focused_thetamax;       /* aperture max angle [rad] */
    gdouble focused_fdist;          /* focal distance */
    gdouble focused_pol;            /* polarization angle [rad]*/
    gint focused_nip;               /* N integration points */
    gint focused_mip;               /* M integration points */ 

    gdouble xshift;
    gdouble yshift;
    gdouble zshift;

    gboolean valid;                 /*TODO: implement in gsvit*/
} SvSrcTSFF;


typedef struct
{
    /* 0a. Box properties */
    gint    box_i0;                 /*LTSF volume boundaries (layered TSF)*/
    gint    box_j0;
    gint    box_k0;
    gint    box_in;
    gint    box_jn;
    gint    box_kn;
    gint    box_boundary_skipi0;    /*skipping of any of the volume faces*/
    gint    box_boundary_skipj0;
    gint    box_boundary_skipk0;
    gint    box_boundary_skipin;
    gint    box_boundary_skipjn;
    gint    box_boundary_skipkn;
    gint    box_boundary_skipdepth; /*skipping by material neighborhood*/

    /* 1. Source properties */
    gint    source_mode;            /*only for xsvit use*/
    gchar   *source_filename;
    gdouble source_wl;              /*only for xsvit use*/
    gdouble source_wlspan;          /*only for xsvit use*/    
    gdouble source_amplitude;       /*only for xsvit use*/
    gdouble source_pulsewidth;      /*only for xsvit use*/
    //gdouble length;               /*only for xsvit use*/

    /* 2. Incident angle */
    gdouble ia_theta;               /*input wave parameters [rad]*/
    gdouble ia_phi;
    gdouble ia_psi;

    /* 4. Layered source properties */
    gint    layered_count;              /* number of layers [0..NUM_LSMP-1]*/
    gint    layered_zpos[NUM_LSMP];
    gdouble layered_epsilon[NUM_LSMP];  /*material parameters for whole TSF, all in z direction*/
    gdouble layered_mu[NUM_LSMP];
    gdouble layered_sigma[NUM_LSMP];
    gdouble layered_sigast[NUM_LSMP];    
    gchar*  layered_material[NUM_LSMP];

    gint    gaussian;
    gdouble gaussian_fxpos;
    gdouble gaussian_fypos;
    gdouble gaussian_rx;
    gdouble gaussian_ry;

    gint    radial;
    gdouble radial_fxpos;
    gdouble radial_fypos;
    gdouble radial_rx;
    gdouble radial_ry;

    gint    fiber;
    gdouble fiber_radius;
    gdouble fiber_fxpos;
    gdouble fiber_fypos;
    gdouble fiber_cutoff;
    gdouble fiber_epsilon_core;
    gdouble fiber_epsilon_cladding;

    gboolean valid;   /*TODO: implement in gsvit*/
} SvSrcLTSF;

/* Layered Focused Total/Scattered field source (LTSFF) */
typedef struct
{
    /* 0a. Box properties */
    gint    box_i0;                 /*LTSFF volume boundaries (focused layered TSF)*/
    gint    box_j0;
    gint    box_k0;
    gint    box_in;
    gint    box_jn;
    gint    box_kn;
    gint    box_boundary_skipi0;    /*skipping of any of the volume faces*/
    gint    box_boundary_skipj0;
    gint    box_boundary_skipk0;
    gint    box_boundary_skipin;
    gint    box_boundary_skipjn;
    gint    box_boundary_skipkn;
    gint    box_boundary_skipdepth; /*skipping by material neighbourhood*/

    /* 1. Source properties */
    gint    source_mode;        /*only for xsvit use*/
    gchar   *source_filename;
    gdouble source_wl;          /*only for xsvit use*/
    gdouble source_wlspan;      /*only for xsvit use*/
    gdouble source_amplitude;   /*only for xsvit use*/
    gdouble source_pulsewidth;  /*only for xsvit use*/

    /* 3. Focused source properties */
    gdouble focused_thetamax;   /* aperture max angle [rad] */
    gdouble focused_fdist;      /* focal distance */
    gdouble focused_pol;        /* polarization angle [rad]*/
    gint focused_nip;           /* N integration points */
    gint focused_mip;           /* M integration points */

    /* 4. Layered source properties */
    gint    layered_count;              /* number of layers [0..NUM_LSMP-1]*/
    gint    layered_zpos[NUM_LSMP];
    gdouble layered_epsilon[NUM_LSMP];  /*material parameters for whole TSF, all in z direction*/
    gdouble layered_mu[NUM_LSMP];
    gdouble layered_sigma[NUM_LSMP];
    gdouble layered_sigast[NUM_LSMP];
    gchar*  layered_material[NUM_LSMP];

    gint fast;        

    gboolean valid;             /*TODO: implement in gsvit*/
} SvSrcLTSFF;

typedef struct
{
    /* 0a. Box properties */
    gint    box_i0;
    gint    box_j0;
    gint    box_k0;
    gint    box_in;
    gint    box_jn;
    gint    box_kn;

    gdouble layered_epsilon;    /*bounding material property to select where to apply*/
    gdouble lambda_peak;        /*generated signal peak wavelength*/
    gdouble lambda_width;       /*generated signal wavelengt span width*/
    gdouble density;            /*spatial density of generating dipoles*/
    gdouble strength;           /*multiplicative factor*/
    gint    source_mode;           /*0: uncontrolled, 1: intensity based, 2: absorption based, 3: intensity and direction based*/
    gint    startfrom;             /*when to start; for mode 1 and 2 this also means when to store actual local intensity to prevent feedback effects*/  
} SvSrcLocal;


/*all the source settings together*/
typedef struct
{
    SvSrcPoint  *pnts;
    gint        npnts;
    gint        npnts_allocated;

    SvSrcExt    ext;
    SvSrcTSF    tsf;
    SvSrcTSFF   tsff;
    SvSrcLTSF   ltsf;
    SvSrcLTSFF  ltsff;
    SvSrcSF     sf;
    SvSrcLocal  *locals;    
    gint        nlocals;
    gdouble     lambda_min;
    gdouble     lambda_max;
    gdouble     lambda_center;
} SvSetSource;


typedef struct {
    gint box_i0;
    gint box_j0;
    gint box_k0;
    gint box_in;
    gint box_jn;
    gint box_kn;
    gint division;
} SvSetSg;


typedef struct {
    SvSetSg sg[20];
    gint nsg;
} SvSetSubgrid;


/*all the computational volume boundary settings together*/
typedef struct
{
    SvBoundaryType bx0;
    SvBoundaryType bxn;
    SvBoundaryType by0;
    SvBoundaryType byn;
    SvBoundaryType bz0;
    SvBoundaryType bzn;

    gint depth_bx0;         /*all the following used for CPML case only*/
    gint m_bx0;
    gdouble sigma_bx0;
    gdouble a_bx0;
    gdouble kappa_bx0;

    gint depth_bxn;
    gint m_bxn;
    gdouble sigma_bxn;
    gdouble a_bxn;
    gdouble kappa_bxn;

    gint depth_by0;
    gint m_by0;
    gdouble sigma_by0;
    gdouble a_by0;
    gdouble kappa_by0;
    
    gint depth_byn;
    gint m_byn;
    gdouble sigma_byn;
    gdouble a_byn;
    gdouble kappa_byn;
    
    gint depth_bz0;
    gint m_bz0;
    gdouble sigma_bz0;
    gdouble a_bz0;
    gdouble kappa_bz0;
    
    gint depth_bzn;
    gint m_bzn;
    gdouble sigma_bzn;
    gdouble a_bzn;
    gdouble kappa_bzn;
} SvSetBoundary;

/*medium boundary settings (boundaries inside computational volume, e.g. periodic)*/
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
    SvBoundaryType bz0;
    SvBoundaryType bzn;
} SvSetMBoundary;

/*space summation output, e.g. for absorption*/
typedef struct
{
    gint component;
    gint box_i0;
    gint box_j0;
    gint box_k0;
    gint box_in;
    gint box_jn;
    gint box_kn;
    gint step;
    gdouble layered_epsilon;
    gdouble layered_mu;
    gdouble layered_sigma;
    gdouble layered_sigast;
    gchar filename[100];
    gboolean stringbased;
    gchar string[100];
} SvOutputSum;

typedef struct
{
    gint box_i0;   /*integration boundary definition*/
    gint box_j0;
    gint box_k0;
    gint box_in;
    gint box_jn;
    gint box_kn;
    gint step;
    gchar filename[100];
} SvOutputForce;

/*all the outputs together*/
typedef struct
{
    SvOutputPar     *pnts;          /*points*/
    SvOutputPar     *plns;          /*same as images but to a ascii/binary file*/
    SvOutputPar     *imgs;          /*images*/
    SvOutputPar     *cubs;          /*volumes*/
    SvOutputSum     *sums;          /*summation results*/
    SvOutputForce   *forces;        /*forces from maxwell stress tensor*/
    SvOutputPar     *subgrid_imgs;  /*images through subgrids*/
    gint            npnts;
    gint            npnts_allocated;
    gint            nplns;
    gint            nplns_allocated;
    gint            nimgs;
    gint            nimgs_allocated;
    gint            ncubs;
    gint            ncubs_allocated;
    gint            nsums;
    gint            nsums_allocated;
    gint            nforces;
    gint            nforces_allocated;
    gint            nsubgrid_imgs;
    gint            nsubgrid_imgs_allocated;
    gchar           *outfile;   /*Gwyddion filename used for most of the outputs*/
    gboolean        savespectrum;
} SvSetOutput;

/*computational volume settings*/
typedef struct {
    gdouble dx;        /*differences at real space*/
    gdouble dy;
    gdouble dz;
    gint xres;         /*computational volume resolution in pixels*/
    gint yres;
    gint zres;
} SvSetPool;

/*computation settings*/
typedef struct {
    gint nsteps;                    /*total number of steps*/ 
    gint step_act;                  /*actual time step*/
    gint verbose;                   /*verbosity level*/
    gint ugpu[MAX_GPUS];            /*use concrete GPUs*/
    gboolean usegpu;                /*use GPU or not*/
    gint nthreads;
    gint suffix;                    /*suffix for temporary files*/
    gint devicequery;               /*run devicequery for cards*/
    gdouble dtmult;                 /*timestep multiplication factor*/
    gchar *killer;                  /*external shutdown control filename*/
} SvSetComp;

/*material input*/
typedef struct {
    gchar* in_voxel_filename;       /*complete material information*/
    gchar* in_vector_filename;      /*script material information*/
    gchar* in_vtk_data_filename;    /*vtk data (integers)*/
    gboolean matmode_check;         /*whether to check material file actively*/
    gint smooth;                    /*if and how many times to run soft blur material average (only for mat type=0)*/    

    //sines are unused now
    gint nsines;                    /*if and how many different sine addition processes should be applied for mat type=1*/
    gdouble sine_theta[20];         /*direction parameters*/
    gdouble sine_phi[20];           /*direction parameters*/
    gdouble sine_wavelength[20];    /*wavelength in the direction, pixel coordinates*/
    gdouble sine_amplitude[20];     /*amplitude of addition strength*/
    gdouble sine_phase[20];         /*sine phase*/
    gboolean sine_randomphase[20];  /*generate random phase or not*/    

    /* grow modifier */
    gint ngrowths;
    gint grow_i0[MATERIAL_COUNT];   /*grow bounding box*/
    gint grow_j0[MATERIAL_COUNT];
    gint grow_k0[MATERIAL_COUNT];
    gint grow_in[MATERIAL_COUNT];
    gint grow_jn[MATERIAL_COUNT];
    gint grow_kn[MATERIAL_COUNT];
    gboolean grow_skipi0[MATERIAL_COUNT];       /*skip some boundaries*/
    gboolean grow_skipin[MATERIAL_COUNT];
    gboolean grow_skipj0[MATERIAL_COUNT];
    gboolean grow_skipjn[MATERIAL_COUNT];
    gboolean grow_skipk0[MATERIAL_COUNT];
    gboolean grow_skipkn[MATERIAL_COUNT];
    gint grow_attachindex[MATERIAL_COUNT];      /* material index where to attach */
    gint grow_addindex[MATERIAL_COUNT];         /* material index to grow (to add) */
    gint grow_subsampling[MATERIAL_COUNT];      /*subsampling level to increase sensitivity*/
    gint grow_nsteps[MATERIAL_COUNT];           /*number of particles to be set*/
    gdouble grow_mobility[MATERIAL_COUNT];      /*mobility length control parameter*/
    gdouble grow_probability[MATERIAL_COUNT];   /*mobility frequency control parameter*/
    gint grow_seed[MATERIAL_COUNT];             /*random seed or generate it (==0)*/

    /* roughness modifier */
    gint nroughens;                             /*if and how many different roughening processes should be applied for mat type=1*/
    gint rough_radius_peak[MATERIAL_COUNT];     /*roughening parameters follow*/
    gint rough_radius_span[MATERIAL_COUNT];
    gint rough_iterations[MATERIAL_COUNT];
    gdouble rough_probability[MATERIAL_COUNT];
    gint rough_matindex[MATERIAL_COUNT];
    gint rough_voidindex[MATERIAL_COUNT];
    gint rough_seed[MATERIAL_COUNT];

    /* spectral modifier */
    gint nspectrals;
    gdouble spectral_sigma[MATERIAL_COUNT];
    gdouble spectral_t[MATERIAL_COUNT];
    gint spectral_seed[MATERIAL_COUNT];
    gint spectral_matindex[MATERIAL_COUNT];

    /* expression modifier */
    gint nexprs;
    gint expr_i0[MATERIAL_COUNT];               /*expression bounding box*/
    gint expr_j0[MATERIAL_COUNT];
    gint expr_k0[MATERIAL_COUNT];
    gint expr_in[MATERIAL_COUNT];
    gint expr_jn[MATERIAL_COUNT];
    gint expr_kn[MATERIAL_COUNT];
    gint expr_matindex[MATERIAL_COUNT];
    gint expr_voidindex[MATERIAL_COUNT];
    gint expr_maxdist[MATERIAL_COUNT];
    gint expr_distmode[MATERIAL_COUNT];
    gchar expr_expr[MATERIAL_COUNT][MATERIAL_EXPR_EXPR_CHARS];

    gboolean crop;
    gint crop_i0;
    gint crop_in;
    gint crop_j0;
    gint crop_jn;
    gint crop_k0;
    gint crop_kn;
    gdouble crop_epsilon;
    gdouble crop_sigma;
    gdouble crop_mu;
    gdouble crop_sigast;    


    gboolean localfiles;        /*use local material files if there are some*/

} SvSetMaterial;

/*computational plan*/
typedef struct {
    gdouble dt;                 /*calculated time difference for all the computations*/
    SvMatMode matmode;          /*material mode to be used*/
    SvGpuMode gpumode;          /*gpu mode to be used*/
} SvPlan;

/*far field computation settings*/
typedef struct {
    gint box_i0;                /*integration surface boundaries*/
    gint box_in;
    gint box_j0;
    gint box_jn;
    gint box_k0;
    gint box_kn;
    gint *ri;                   /*calculation positions in voxel values*/
    gint *rj;
    gint *rk;
    gint *individual;           /*individually output, now always true*/
    gchar **source_filename;
    gint sumfrom;               /*unused now*/
    gint sumto;                 /*unused now*/
    gint nrs;                   /*number of calculation positions*/
    gint nsets;                 /*total number of far field point sets for collective output*/
    gint *setxres;              /*x resolution of the set*/
    gint *setyres;              /*y resolution of the set*/
    gint box_boundary_skipi0;   /*skipping boundaries or parts of them*/
    gint skipi0_jmin, skipi0_kmin, skipi0_jmax, skipi0_kmax;
    gint box_boundary_skipin;
    gint skipin_jmin, skipin_kmin, skipin_jmax, skipin_kmax;

    gint box_boundary_skipj0;
    gint skipj0_imin, skipj0_kmin, skipj0_imax, skipj0_kmax;
    gint box_boundary_skipjn;
    gint skipjn_imin, skipjn_kmin, skipjn_imax, skipjn_kmax;

    gint box_boundary_skipk0;
    gint skipk0_imin, skipk0_jmin, skipk0_imax, skipk0_jmax;
    gint box_boundary_skipkn;
    gint skipkn_imin, skipkn_jmin, skipkn_imax, skipkn_jmax;

    gint nareas;                /*for xsvit only*/
    gint area_thetares[20];     /*for xsvit only*/
    gint area_phires[20];       /*for xsvit only*/
    gint area_radius[20];       /*for xsvit only*/
    gdouble area_thetafrom[20]; /*for xsvit only*/
    gdouble area_phifrom[20];   /*for xsvit only*/
    gdouble area_thetato[20];   /*for xsvit only*/
    gdouble area_phito[20];     /*for xsvit only*/
    gint area_savefile[20];     /*for xsvit only*/

    gint nsquares;              /*all the next for xsvit only as well*/
    gint square_ijres[20];      /*first axis resolution, regardless absolute orientation (in fact x, or y)*/
    gint square_jkres[20];      /*second axis resolution, regardless absolute orientation (in fact y or z)*/
    gdouble square_ijfrom[20];
    gdouble square_jkfrom[20];
    gdouble square_ijto[20];
    gdouble square_jkto[20];
    gint square_orientation[20];    /*0, 1, 2 means x, y, z*/
    gdouble square_distance[20];    /*distance from zero in appropriate direction*/
    gdouble square_savefile[20];
} SvSetFarfield;

/*periodic far field computation settings*/
typedef struct {
    gint box_i0;        /*integration surface boundaries*/
    gint box_in;
    gint box_j0;
    gint box_jn;
    gint box_k0;
    gint box_kn;
    gint pimin;         /*periodicity range in i (x)*/
    gint pimax;    
    gint pjmin;         /*periodicity range in j (y)*/
    gint pjmax;
    gint *ri;           /*calculation positions in voxel values*/
    gint *rj;
    gint *rk;
    gint *individual;   /*individually output (0), or use it as a set (any positive)*/
    gchar **source_filename;
    gint nrs;           /*total number of far field points*/
    gint nsets;         /*total number of far field point sets for collective output*/
    gint *setxres;      /*x resolution of the set*/
    gint *setyres;      /*y resolution of the set*/
    gint box_boundary_skipk0;
    gint skipk0_imin, skipk0_jmin, skipk0_imax, skipk0_jmax;
    gint box_boundary_skipkn;
    gint skipkn_imin, skipkn_jmin, skipkn_imax, skipkn_jmax;

    gint nareas;                /*for xsvit only*/
    gint area_thetares[20];     /*for xsvit only*/
    gint area_phires[20];       /*for xsvit only*/
    gint area_radius[20];       /*for xsvit only*/
    gdouble area_thetafrom[20]; /*for xsvit only*/
    gdouble area_phifrom[20];   /*for xsvit only*/
    gdouble area_thetato[20];   /*for xsvit only*/
    gdouble area_phito[20];     /*for xsvit only*/
    gint area_savefile[20];     /*for xsvit only*/
    
    gint postprocess;
    gint ppstart;
} SvSetPFarfield;

typedef struct _SvSet {
    SvSetPool sp;         /*data pool settings*/
    SvSetComp sc;         /*computation settings*/
    SvSetOutput so;       /*output settings*/
    SvSetSource ss;       /*source settings*/
    SvSetBoundary sb;     /*boundary settings*/
    SvSetMBoundary smb;   /*inside medium boundary settings*/
    SvSetMaterial sm;     /*material settings*/
    SvSetFarfield sf;     /*far field computation*/
    SvSetPFarfield spf;   /*periodic far field computation*/
    SvSetSubgrid sg;      /*local mesh refinement*/
    SvPlan plan;          /*computation plan*/
} SvSet;


typedef struct _SvSetMat {
    GArray *spheres;
    GArray *voxels;
    GArray *cylinders;
    GArray *cones;
    GArray *rcones;
    GArray *tetrahedrons;
    GArray *gwydds;

    gint   ngtot;

    gint        nmats;
    gchar       **materials;
    gint        *overpos;
    //gboolean    *valid;

    SvDCube *gwyddata;
    gdouble gwydd_xmin[100];
    gdouble gwydd_xmax[100];
    gdouble gwydd_ymin[100];
    gdouble gwydd_ymax[100];
    gdouble gwydd_zmin[100];
    gdouble gwydd_zmax[100];
    gint gwydd_nvx[100];

    gint ntetgens;
    gint tetgen_start[100];
    gint tetgen_end[100];
    gdouble tetgen_xshift[100];
    gdouble tetgen_yshift[100];
    gdouble tetgen_zshift[100];
    gdouble tetgen_xmult[100];
    gdouble tetgen_ymult[100];
    gdouble tetgen_zmult[100];
    gint tetgen_attribute_pos[100];
    gint tetgen_attribute_val[100];
    gchar tetgen_filebase[100][100];
} SvSetMat;

/*material properties*/
typedef struct {
    // input data
    gint type;
    gdouble epsilon;
    gdouble mu;
    gdouble sigma;
    gdouble sigast;
    gdouble drude_omega_p;
    gdouble drude_nu;
    gdouble cp3_a[3];
    gdouble cp3_phi[3];
    gdouble cp3_omega[3];
    gdouble cp3_gamma[3];

    // calculated from input data:
    gdouble ade_a0;
    gdouble ade_a1;
    gdouble ade_a2;
    gdouble ade_bp0[2];
    gdouble ade_bp1[2];
    gdouble ade_bp2[2];
    gdouble ade_bp3[2];
    gdouble ade_bp4[2];
    gdouble ade_c0;
    gdouble ade_c1;
    gdouble ade_c2;
    gdouble ade_c3;
    gdouble plrc_d_chi;
    gdouble plrc_d_xi;
    gdouble plrc_d_dchi;
    gdouble plrc_d_dxi;
    gdouble plrc_p_chi[2];
    gdouble plrc_p_xi[2];
    gdouble plrc_p_dchir[2];
    gdouble plrc_p_dxir[2];
    gdouble plrc_p_dchii[2];
    gdouble plrc_p_dxii[2];

    // loaded from database
    //gchar *material;

    // asigned:
    gint pos;
} SvMatProp;


int parse_settings(gchar *filename, SvSet *set, gboolean called_from_gsvit);
int parse_settings_mat(gchar *filename, SvSet *set_par, SvSetMat *set_mat, gboolean called_from_gsvit);
//gint parse_settings_from_string(const gchar *str, SvSet *set, gboolean called_from_gsvit);
//gint parse_settings_from_file(gchar *filename, SvSet *set, gboolean called_from_gsvit);

void delete_point_sources(SvSet *set);

void clear_settings(SvSet *set, gboolean remove);
void clear_settings_mat(SvSetMat *set, gboolean remove);

void init_settings(SvSet *set, SvType type);
void init_add_settings_mat(SvSetMat *set, SvTypeMat type);
void init_settings_mat_prop(SvMatProp *mat);

gboolean is_source_valid(SvSet *set, SvType type);

int write_settings(gchar *filename, SvSet *set);

GString *sv_check(SvSet *set, GString *message, gboolean *ok);

void remove_temporary_files(SvSet *set);

gint get_double(FILE *F, gdouble *val, gchar *key);

gint get_int(FILE *F, gint *val, gchar *key);

gchar* get_self_dir();

gint get_self_dir_ex(gchar* buf, gint size);

void alloc_set_par_mat(SvSet *set_par, SvSetMat *set_mat, gboolean init_data);

void alloc_set_mat(SvSetMat *set);

#endif /* SETTINGS_H */

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */


