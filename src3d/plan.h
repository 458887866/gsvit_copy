
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


/*  plan.h : 
 *  create plan for computation, load data files
 */

#ifndef PLAN_H
#define PLAN_H

#include "settings.h"
#include "pool.h"

#ifdef G_OS_WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <GL/gl.h>


typedef enum {
   SV_ENTITY_SPHERE = 4,
   SV_ENTITY_CYLINDER = 7,
   SV_ENTITY_VOXEL = 8,
   SV_ENTITY_TETRAHEDRON = 12,
   SV_ENTITY_CONE = 107,
   SV_ENTITY_RCONE = 108,
   SV_ENTITY_HEIGHTFIELD = 20,
   SV_ENTITY_TETGEN = 21,
   SV_ENTITY_GWYDDION = 22
} SvEntityType;


typedef struct {
    GLfloat pnt1[3];
    GLfloat pnt2[3];
    GLfloat pnt3[3];
    GLfloat pnt4[3];
    SvMatProp mat;
    int n;
    int setpart;
} SvTetrahedron;

typedef struct {
    GLfloat pnt1[3];
    GLfloat pnt2[3];
    SvMatProp mat;
    int n;
} SvVoxel;

typedef struct {
    GLfloat pnt1[3];
    GLfloat pnt2[3];
    GLfloat radius;
    SvMatProp mat;
    int n;
} SvCylinder;

typedef struct {
    GLfloat pnt1[3];
    GLfloat pnt2[3];
    GLfloat radius;
    SvMatProp mat;
    int n;
} SvCone;

typedef struct {
    GLfloat pnt1[3];
    GLfloat pnt2[3];
    GLfloat radius1;
    GLfloat radius2;
    SvMatProp mat;
    int n;
} SvRCone;

typedef struct {
    GLfloat pnt1[3];
    GLfloat radius;
    SvMatProp mat;
    int n;
} SvSphere;

typedef struct {
    GwyDataField *dfield;
    GwyDataField *mfield;
    GLfloat xoffset;
    GLfloat yoffset;
    GLfloat zoffset;
    gint i;
    gint j;
    gint k;
    gint depth;
    gint mask;
    SvMatProp mat;
    int n;
    gint channel;        //for xsvit use only
    gchar filebase[256]; //for xsvit use only
} SvGwydd;



SvPool* init_and_make_plan(SvSet *set);

gboolean scan_point(FILE *fr, GLfloat *px, GLfloat *py, GLfloat *pz);
gboolean scan_radius(FILE *fr, GLfloat *pradius);
gboolean scan_int(FILE *fr, gint *pint);
gboolean scan_material(FILE *fr, SvMatProp *mat, gdouble dt);
SvMatMode test_script_file(gchar *filename);
gboolean override_matfile_from_database(gchar *filename, gdouble wf, gdouble wc, gdouble wt, gint suffix, gint *nmaterials, gchar **materials, gint *overpos, gboolean localfiles);
gint load_tets(gchar *filebase, GArray *tetrahedrons, gint attribute_pos, gint attribute_val, gdouble xshift, gdouble yshift, gdouble zshift, gdouble xmult, gdouble ymult, gdouble zmult,
               gint verbose, SvMatProp *mat, gint *ngtot);
gint load_gwydd(gchar *filename, SvGwydd *gx, gint verbose, gint channel);

#endif /* PLAN_H */
/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
