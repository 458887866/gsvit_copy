/*\
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


/*  scene.c :
*  3D scene functions
*/

#include "xs_scene.h"
#include <gtk/gtk.h>
#include <gtk/gtkgl.h>
#include <math.h>
#include <GL/glu.h>


#define ZOOM_STEP   0.2f
#define SCALE_NEW

//FILE *fr = NULL;


void ShowMatrix()
{
    GLdouble modelMatrix[16];
    GLdouble projMatrix[16];
    GLint viewport[4];
    glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
    glGetIntegerv(GL_VIEWPORT, viewport);
}

/*******************************    all the OpenGL stuff ********************************************/

//////////////////////////////////////////////////////////////////////////
/* draw functions*/

static void
draw_cone(GLfloat radius, GLfloat length, GLfloat x, GLfloat y, GLfloat z, GLfloat ax, GLfloat ay, GLfloat az)
{
    gint i, vcount = 30;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glTranslatef(x, y, z);
    glRotatef(ax, 1, 0, 0);
    glRotatef(ay, 0, 1, 0);
    glRotatef(az, 0, 0, 1);

    glBegin(GL_TRIANGLE_FAN);
    glVertex3f(0, 0, 0);
    for (i = 0; i < vcount; ++i) {
        glVertex3f((GLfloat)(radius*sin(i * 2 * G_PI / (vcount - 1))), (GLfloat)(radius*cos(i * 2 * G_PI / (vcount - 1))), length);
    }
    glEnd();

    glPopMatrix();
}

/* drawing sphere taken from http://ozark.hendrix.edu/~burch/cs/490/sched/feb8/ */

static void
draw_sphere(GLfloat radius, GLfloat x, GLfloat y, GLfloat z, GLfloat cr, GLfloat cg, GLfloat cb, GLfloat alpha)
{
    gint i, j, n = 36, m = 36;
    GLfloat lat0, z0, zr0, lat1, z1, zr1, lng, px, py;


    glPushMatrix();
    glTranslatef(x, y, z);
    glScalef(radius, radius, radius);
    for (i = 0; i <= n; i++) {
        lat0 = (GLfloat)(G_PI * (-0.5 + (GLfloat)(i - 1) / n));
        z0 = (GLfloat)sin(lat0);
        zr0 = (GLfloat)cos(lat0);

        lat1 = (GLfloat)(G_PI * (-0.5 + (GLfloat)i / n));
        z1 = (GLfloat)sin(lat1);
        zr1 = (GLfloat)cos(lat1);

        glBegin(GL_QUAD_STRIP);
        glColor4f(cr, cg, cb, alpha);
        for (j = 0; j <= m; j++) {
            lng = (GLfloat)(2 * G_PI * (GLfloat)(j - 1) / m);
            px = (GLfloat)cos(lng);
            py = (GLfloat)sin(lng);

            glNormal3f(px * zr0, py * zr0, z0);
            glVertex3f(px * zr0, py * zr0, z0);
            glNormal3f(px * zr1, py * zr1, z1);
            glVertex3f(px * zr1, py * zr1, z1);
        }
        glEnd();
    }
    glPopMatrix();
}

static void
draw_pointsource(GLfloat x, GLfloat y, GLfloat z, GLfloat cr, GLfloat cg, GLfloat cb, GLfloat alpha)
{
    gint ms = 20, ls = 13;
    draw_sphere(10, x, y, z, cr, cg, cb, alpha);

    glBegin(GL_LINES);
    glColor3f(cr, cg, cb);

    glVertex3f(x, y, z);
    glVertex3f(x + ms, y, z);
    glVertex3f(x, y, z);
    glVertex3f(x - ms, y, z);
    glVertex3f(x, y, z);
    glVertex3f(x, y + ms, z);
    glVertex3f(x, y, z);
    glVertex3f(x, y - ms, z);
    glVertex3f(x, y, z);
    glVertex3f(x, y, z + ms);
    glVertex3f(x, y, z);
    glVertex3f(x, y, z - ms);
    glVertex3f(x, y, z);
    glVertex3f(x + ls, y + ls, z + ls);
    glVertex3f(x, y, z);
    glVertex3f(x + ls, y + ls, z - ls);
    glVertex3f(x, y, z);
    glVertex3f(x + ls, y - ls, z + ls);
    glVertex3f(x, y, z);
    glVertex3f(x + ls, y - ls, z - ls);
    glVertex3f(x, y, z);
    glVertex3f(x - ls, y + ls, z + ls);
    glVertex3f(x, y, z);
    glVertex3f(x - ls, y + ls, z - ls);
    glVertex3f(x, y, z);
    glVertex3f(x - ls, y - ls, z + ls);
    glVertex3f(x, y, z);
    glVertex3f(x - ls, y - ls, z - ls);
    glEnd();
}

static void
draw_outpoint(GLfloat x, GLfloat y, GLfloat z, GLfloat cr, GLfloat cg, GLfloat cb, GLfloat alpha)
{
    gint ms = 10;
    draw_sphere(2, x, y, z, cr, cg, cb, alpha);

    glBegin(GL_LINES);
    glColor3f(cr, cg, cb);

    glVertex3f(x, y, z);
    glVertex3f(x + ms, y, z);
    glVertex3f(x, y, z);
    glVertex3f(x - ms, y, z);
    glVertex3f(x, y, z);
    glVertex3f(x, y + ms, z);
    glVertex3f(x, y, z);
    glVertex3f(x, y - ms, z);
    glVertex3f(x, y, z);
    glVertex3f(x, y, z + ms);
    glVertex3f(x, y, z);
    glVertex3f(x, y, z - ms);
    glVertex3f(x, y, z);
    glEnd();
}

static void
draw_cylinder(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat radius, GLfloat cr, GLfloat cg, GLfloat cb, GLfloat alpha)
{
    gint i;
    GLfloat x, y;
    GLfloat length = (GLfloat)sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) + (z2 - z1)*(z2 - z1));
    GLfloat a = 0;
    GLfloat step = (GLfloat)(2.*G_PI / 36.0);
    GLfloat dx, dy, dz, ax, ay, az, nx, ny, nz, omega;

    nx = 0;
    ny = 0;
    nz = 1;
    dx = x2 - x1;
    dy = y2 - y1;
    dz = z2 - z1;

    ax = ny*dz - nz*dy;
    ay = nz*dx - nx*dz;
    az = nx*dy - ny*dx;
    omega = (GLfloat)(180.0f / G_PI*acos((nx*dx + ny*dy + nz*dz) / sqrt(dx*dx + dy*dy + dz*dz)));

    glPushMatrix();
    glTranslatef(x1, y1, z1);
    glRotatef(omega, ax, ay, az);

    glBegin(GL_TRIANGLE_STRIP);
    glColor4f(cr, cg, cb, alpha);
    for (i = 0; i <= 36; ++i) {
        x = (GLfloat)(cos(a)*radius);
        y = (GLfloat)(sin(a)*radius);
        glVertex3f(x, y, 0);
        glVertex3f(x, y, length);

        a += step;
    }
    glEnd();
    glBegin(GL_POLYGON);
    for (i = 0; i <= 36; ++i) {
        glVertex3f((GLfloat)(radius*sin(i*step)), (GLfloat)(radius*cos(i*step)), 0);
    }
    glEnd();
    glBegin(GL_POLYGON);
    for (i = 0; i <= 36; ++i) {
        glVertex3f((GLfloat)(radius*sin(i*step)), (GLfloat)(radius*cos(i*step)), length);
    }
    glEnd();
    glBegin(GL_LINE_LOOP);
    glColor4f(cr, cg, cb, 1);
    for (i = 0; i <= 36; ++i) {
        glVertex3f((GLfloat)(radius*sin(i*step)), (GLfloat)(radius*cos(i*step)), 0);
    }
    glEnd();
    glBegin(GL_LINE_LOOP);
    for (i = 0; i <= 36; ++i) {
        glVertex3f((GLfloat)(radius*sin(i*step)), (GLfloat)(radius*cos(i*step)), length);
    }
    glEnd();
    glPopMatrix();
}

static void
draw_tetrahedron(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat x3, GLfloat y3, GLfloat z3, GLfloat x4, GLfloat y4, GLfloat z4,
                 GLfloat cr, GLfloat cg, GLfloat cb, GLfloat alpha)
{
    glBegin(GL_TRIANGLES);
    glColor4f(cr, cg, cb, alpha);
    glVertex3f(x1, y1, z1);
    glVertex3f(x2, y2, z2);
    glVertex3f(x3, y3, z3);

    glVertex3f(x1, y1, z1);
    glVertex3f(x3, y3, z3);
    glVertex3f(x4, y4, z4);

    glVertex3f(x1, y1, z1);
    glVertex3f(x2, y2, z2);
    glVertex3f(x4, y4, z4);

    glVertex3f(x2, y2, z2);
    glVertex3f(x3, y3, z3);
    glVertex3f(x4, y4, z4);
    glEnd();
}


static void
draw_rcone(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat radius1, GLfloat radius2, GLfloat cr, GLfloat cg, GLfloat cb, GLfloat alpha)
{
    gint i;
    GLfloat length = (GLfloat)sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) + (z2 - z1)*(z2 - z1));
    GLfloat a = 0;
    GLfloat step = (GLfloat)(2.0f*G_PI / 36.0);
    GLfloat dx, dy, dz, ax, ay, az, nx, ny, nz, omega;

    nx = 0;
    ny = 0;
    nz = 1;
    dx = x2 - x1;
    dy = y2 - y1;
    dz = z2 - z1;

    ax = ny*dz - nz*dy;
    ay = nz*dx - nx*dz;
    az = nx*dy - ny*dx;
    omega = (GLfloat)(180.0f / G_PI*acos((nx*dx + ny*dy + nz*dz) / sqrt(dx*dx + dy*dy + dz*dz)));

    glPushMatrix();
    glTranslatef(x1, y1, z1);
    glRotatef(omega, ax, ay, az);

    glBegin(GL_TRIANGLE_STRIP);
    glColor4f(cr, cg, cb, alpha);
    for (i = 0; i <= 36; ++i) {
        x1 = (GLfloat)(cos(a)*radius1);
        y1 = (GLfloat)(sin(a)*radius1);
        x2 = (GLfloat)(cos(a)*radius2);
        y2 = (GLfloat)(sin(a)*radius2);
        glVertex3f(x1, y1, 0);
        glVertex3f(x2, y2, length);

        a += step;
    }
    glEnd();
    glBegin(GL_POLYGON);
    for (i = 0; i <= 36; ++i) {
        glVertex3f((GLfloat)(radius1*sin(i*step)), (GLfloat)(radius1*cos(i*step)), 0);
    }
    glEnd();
    glBegin(GL_POLYGON);
    for (i = 0; i <= 36; ++i) {
        glVertex3f((GLfloat)(radius2*sin(i*step)), (GLfloat)(radius2*cos(i*step)), length);
    }
    glEnd();
    glBegin(GL_LINE_LOOP);
    glColor4f(cr, cg, cb, 1);
    for (i = 0; i <= 36; ++i) {
        glVertex3f((GLfloat)(radius1*sin(i*step)), (GLfloat)(radius1*cos(i*step)), 0);
    }
    glEnd();
    glBegin(GL_LINE_LOOP);
    for (i = 0; i <= 36; ++i) {
        glVertex3f((GLfloat)(radius2*sin(i*step)), (GLfloat)(radius2*cos(i*step)), length);
    }
    glEnd();
    glPopMatrix();
}


static void
draw_conem(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat radius, GLfloat cr, GLfloat cg, GLfloat cb, GLfloat alpha)
{
    gint i;
    GLfloat x, y;
    GLfloat length = (GLfloat)sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) + (z2 - z1)*(z2 - z1));
    GLfloat a = 0;
    GLfloat step = (GLfloat)(2.0f*G_PI / 36.0);
    GLfloat dx, dy, dz, ax, ay, az, nx, ny, nz, omega;

    nx = 0;
    ny = 0;
    nz = 1;
    dx = x2 - x1;
    dy = y2 - y1;
    dz = z2 - z1;

    ax = ny*dz - nz*dy;
    ay = nz*dx - nx*dz;
    az = nx*dy - ny*dx;
    omega = (GLfloat)(180.0f / G_PI*acos((nx*dx + ny*dy + nz*dz) / sqrt(dx*dx + dy*dy + dz*dz)));

    glPushMatrix();
    glTranslatef(x1, y1, z1);
    glRotatef(omega, ax, ay, az);

    glBegin(GL_TRIANGLE_STRIP);
    glColor4f(cr, cg, cb, alpha);
    for (i = 0; i <= 36; ++i) {
        x = (GLfloat)(cos(a)*radius);
        y = (GLfloat)(sin(a)*radius);
        glVertex3f(x, y, 0);
        glVertex3f(0, 0, length);

        a += step;
    }
    glEnd();
    glBegin(GL_POLYGON);
    for (i = 0; i <= 36; ++i) {
        glVertex3f((GLfloat)(radius*sin(i*step)), (GLfloat)(radius*cos(i*step)), 0);
    }
    glEnd();
    glBegin(GL_LINE_LOOP);
    glColor4f(cr, cg, cb, 1);
    for (i = 0; i <= 36; ++i) {
        glVertex3f((GLfloat)(radius*sin(i*step)), (GLfloat)(radius*cos(i*step)), 0);
    }
    glEnd();
    glPopMatrix();
}

static void
draw_line(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat cr, GLfloat cg, GLfloat cb, GLfloat alpha)
{
    glBegin(GL_LINES);
    glColor3f(cr, cg, cb);
    glVertex3f(z1, y1, z1);
    glVertex3f(x2, y2, z2);
    glEnd();
}

/*
static void
draw_arrow(gdouble x1, gdouble y1, gdouble z1, gdouble x2, gdouble y2, gdouble z2, gdouble cr, gdouble cg, gdouble cb, gdouble alpha)
{
double ax, ay, az;

draw_line(x1, y1, z1, x2, y2, z2, cr, cg, cb, alpha);
draw_cone(5, 20, x2, y2, z2, ax, ay, az);

//draw_cone(5, -20, 100.0, 0, 0, 0, 90, 0);


}*/

static void
fill_pp(GLfloat i0, GLfloat j0, GLfloat k0, GLfloat i1, GLfloat j1, GLfloat k1,
        GLfloat acr, GLfloat acg, GLfloat acb,
        GLfloat aalpha)
{
    glBegin(GL_QUADS);
    glColor4f(acr, acg, acb, aalpha);
    glVertex3f(i0, j0, k0);
    glVertex3f(i1, j0, k0);
    glVertex3f(i1, j1, k0);
    glVertex3f(i0, j1, k0);

    glVertex3f(i0, j0, k0);
    glVertex3f(i0, j0, k1);
    glVertex3f(i1, j0, k1);
    glVertex3f(i1, j0, k0);

    glVertex3f(i1, j1, k0);
    glVertex3f(i1, j1, k1);
    glVertex3f(i0, j1, k1);
    glVertex3f(i0, j1, k0);

    glVertex3f(i0, j0, k1);
    glVertex3f(i1, j0, k1);
    glVertex3f(i1, j1, k1);
    glVertex3f(i0, j1, k1);

    glVertex3f(i0, j0, k0);
    glVertex3f(i0, j1, k0);
    glVertex3f(i0, j1, k1);
    glVertex3f(i0, j0, k1);

    glVertex3f(i1, j0, k0);
    glVertex3f(i1, j0, k1);
    glVertex3f(i1, j1, k1);
    glVertex3f(i1, j1, k0);
    glEnd();
}


static void
draw_pp(GLfloat i0, GLfloat j0, GLfloat k0, GLfloat i1, GLfloat j1, GLfloat k1,
        GLfloat lcr, GLfloat lcg, GLfloat lcb, GLfloat acr, GLfloat acg, GLfloat acb,
        GLfloat aalpha)
{
    glBegin(GL_QUADS);
    glColor4f(acr, acg, acb, aalpha);
    glVertex3f(i0, j0, k0);
    glVertex3f(i1, j0, k0);
    glVertex3f(i1, j1, k0);
    glVertex3f(i0, j1, k0);

    glVertex3f(i0, j0, k0);
    glVertex3f(i0, j0, k1);
    glVertex3f(i1, j0, k1);
    glVertex3f(i1, j0, k0);

    glVertex3f(i1, j1, k0);
    glVertex3f(i1, j1, k1);
    glVertex3f(i0, j1, k1);
    glVertex3f(i0, j1, k0);

    glVertex3f(i0, j0, k1);
    glVertex3f(i1, j0, k1);
    glVertex3f(i1, j1, k1);
    glVertex3f(i0, j1, k1);

    glVertex3f(i0, j0, k0);
    glVertex3f(i0, j1, k0);
    glVertex3f(i0, j1, k1);
    glVertex3f(i0, j0, k1);

    glVertex3f(i1, j0, k0);
    glVertex3f(i1, j0, k1);
    glVertex3f(i1, j1, k1);
    glVertex3f(i1, j1, k0);
    glEnd();

    glBegin(GL_LINES);
    glColor3f(lcr, lcg, lcb);
    glVertex3f(i0, j0, k0);
    glVertex3f(i1, j0, k0);

    glVertex3f(i1, j0, k0);
    glVertex3f(i1, j1, k0);

    glVertex3f(i1, j1, k0);
    glVertex3f(i0, j1, k0);

    glVertex3f(i0, j1, k0);
    glVertex3f(i0, j0, k0);

    glVertex3f(i0, j0, k1);
    glVertex3f(i1, j0, k1);

    glVertex3f(i1, j0, k1);
    glVertex3f(i1, j1, k1);

    glVertex3f(i1, j1, k1);
    glVertex3f(i0, j1, k1);

    glVertex3f(i0, j1, k1);
    glVertex3f(i0, j0, k1);

    glVertex3f(i0, j0, k0);
    glVertex3f(i0, j0, k1);

    glVertex3f(i1, j0, k0);
    glVertex3f(i1, j0, k1);

    glVertex3f(i1, j1, k0);
    glVertex3f(i1, j1, k1);

    glVertex3f(i0, j1, k0);
    glVertex3f(i0, j1, k1);
    glEnd();
}

static void
draw_datafield(GwyDataField *datafield, gint i0, gint j0, gint k0, gint i1, gint j1, gint k1, gint axis, gboolean square)
{
    gint col, row, xres, yres;
    GLfloat shift, min, max, scale, val;
    GLfloat aalpha = 0.5;
    gdouble* data;

    if (!datafield)
        return;


    xres = gwy_data_field_get_xres(datafield);
    yres = gwy_data_field_get_yres(datafield);

    min = (GLfloat)gwy_data_field_get_min(datafield);
    max = (GLfloat)gwy_data_field_get_max(datafield);

    if (square) {
        shift = 0;
        max = MAX(min*min, max*max);
    } else {
        shift = min;
    }

    if (max > shift)
        scale = (GLfloat)(1.0 / (max - shift));
    else
        scale = 0;

    data = gwy_data_field_get_data(datafield);

    glBegin(GL_QUADS);
    if (axis == 0) {
        for (row = 0; row < yres; row++) {
            for (col = 0; col < xres; col++) {
                if (square)
                    val = (GLfloat)(scale*(data[col + xres*row] * data[col + xres*row]));
                else
                    val = (GLfloat)(scale*(data[col + xres*row] - shift));

                if (val < 0.5)
                    glColor4f(0, 0, (GLfloat)(2 * val), aalpha);
                else
                    glColor4f((GLfloat)(2 * (val - 0.5)), (GLfloat)(2 * (val - 0.5)), 1, aalpha);

                glVertex3f((GLfloat)i0, (GLfloat)(j0 + col), (GLfloat)(k0 + row));
                glVertex3f((GLfloat)i0, (GLfloat)(j0 + (col + 1)), (GLfloat)(k0 + row));
                glVertex3f((GLfloat)i0, (GLfloat)(j0 + (col + 1)), (GLfloat)(k0 + (row + 1)));
                glVertex3f((GLfloat)i0, (GLfloat)(j0 + (col)), (GLfloat)(k0 + (row + 1)));
            }
        }
    }
    if (axis == 1) {
        for (row = 0; row < yres; row++) {
            for (col = 0; col < xres; col++) {
                if (square)
                    val = (GLfloat)(scale*(data[col + xres*row] * data[col + xres*row]));
                else
                    val = (GLfloat)(scale*(data[col + xres*row] - shift));

                if (val < 0.5)
                    glColor4f(0, 0, 2 * val, aalpha);
                else
                    glColor4f((GLfloat)(2 * (val - 0.5)), (GLfloat)(2 * (val - 0.5)), 1, aalpha);

                glVertex3f((GLfloat)(i0 + col), (GLfloat)j0, (GLfloat)(k0 + row));
                glVertex3f((GLfloat)(i0 + (col + 1)), (GLfloat)j0, (GLfloat)(k0 + row));
                glVertex3f((GLfloat)(i0 + (col + 1)), (GLfloat)j0, (GLfloat)(k0 + (row + 1)));
                glVertex3f((GLfloat)(i0 + (col)), (GLfloat)j0, (GLfloat)(k0 + (row + 1)));
            }
        }
    }
    if (axis == 2) {
        for (row = 0; row < yres; row++) {
            for (col = 0; col < xres; col++) {
                if (square)
                    val = (GLfloat)(scale*(data[col + xres*row] * data[col + xres*row]));
                else val = (GLfloat)(scale*(data[col + xres*row] - shift));

                if (val < 0.5)
                    glColor4f(0, 0, 2 * val, aalpha);
                else
                    glColor4f((GLfloat)(2 * (val - 0.5)), (GLfloat)(2 * (val - 0.5)), 1, aalpha);
                glVertex3f((GLfloat)(i0 + col), (GLfloat)(j0 + row), (GLfloat)(k0));
                glVertex3f((GLfloat)(i0 + (col + 1)), (GLfloat)(j0 + row), (GLfloat)(k0));
                glVertex3f((GLfloat)(i0 + (col + 1)), (GLfloat)(j0 + (row + 1)), (GLfloat)(k0));
                glVertex3f((GLfloat)(i0 + (col)), (GLfloat)(j0 + (row + 1)), (GLfloat)(k0));
            }
        }
    }
    glEnd();
}

/* end of draw functions */
//////////////////////////////////////////////////////////////////////////


gboolean
scene_setup(GtkWidget *view_scene)
{
    GdkGLConfig *glconfig;

	glconfig = gdk_gl_config_new_by_mode(GDK_GL_MODE_RGB |
	                                     GDK_GL_MODE_DEPTH |
	                                     GDK_GL_MODE_DOUBLE);
	if (!glconfig) {
	    printf("Error: cannot start OpenGL rendering\n");
	    return FALSE;
	}
	if (!gtk_widget_set_gl_capability(view_scene, glconfig, NULL, TRUE, GDK_GL_RGBA_TYPE)) {
	    printf("Error: cannot start OpenGL rendering\n");
	    return FALSE;
	}

	return TRUE;
}


void
scene_reset_cb(GtkWidget *widget, XGControls *xgc)
{
    xgc->zoomfactor = 1.0f;
    xgc->viewposx = 0.0f;
    xgc->viewposy = 0.0f;
    xgc->viewposx_start = 0.0f;
    xgc->viewposy_start = 0.0f;
    xgc->viewposx_prev = 0.0f;
    xgc->viewposy_prev = 0.0f;
    xgc->xrot = 0.0f;
    xgc->yrot = 0.0f;

    //////////////////////////////////////////////////////////////////////////
    glLoadIdentity();
    glViewport(0, 0, xgc->view_scene->allocation.width, xgc->view_scene->allocation.height);
    gdouble aspect = (gdouble)xgc->view_scene->allocation.width / (gdouble)xgc->view_scene->allocation.height;

    if (xgc->view_scene->allocation.width < xgc->view_scene->allocation.height)
        glOrtho(-150, 150, -150 / aspect, 150 / aspect, -20050, 10000);
    else
        glOrtho(-150 * aspect, 150 * aspect, -150, 150, -20050, 10000);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //////////////////////////////////////////////////////////////////////////

    gtk_widget_queue_draw(xgc->view_scene);
}

//////////////////////////////////////////////////////////////////////////
/* scene callback functions */

gboolean 
scene_configure_cb(GtkWidget *da, GdkEventConfigure *event, XGControls *xgc)
{
    GdkGLContext *glcontext = gtk_widget_get_gl_context(da);
    GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(da);

    if (!gdk_gl_drawable_gl_begin(gldrawable, glcontext)) {
        g_assert_not_reached();
    }

    glLoadIdentity();
    glViewport(0, 0, da->allocation.width, da->allocation.height);

    gdouble aspect = (gdouble)da->allocation.width / (gdouble)da->allocation.height;

    if (da->allocation.width < da->allocation.height)
        glOrtho(-150, 150, -150 / aspect, 150 / aspect, -20050, 10000);
    else
        glOrtho(-150 * aspect, 150 * aspect, -150, 150, -20050, 10000);

    
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    gdk_gl_drawable_gl_end(gldrawable);

    return TRUE;
}

void
render_scene(XGControls *xgc)
{
    guint m;
    gint i, j, k;
    GLfloat xoff, yoff, zoff;
    SvSphere sx;
    SvVoxel vx;
    SvCylinder cx;
    SvCone cnx;
    SvRCone rcnx;
    SvTetrahedron tx;
    SvGwydd gx;


    if (xgc->data.set.sp.xres > 0 && xgc->data.set.sp.yres > 0 && xgc->data.set.sp.zres > 0) {
        if (xgc->data.set.sp.xres > 0)
            xoff = -xgc->data.set.sp.xres / 2.0f;
        if (xgc->data.set.sp.yres > 0)
            yoff = -xgc->data.set.sp.yres / 2.0f;
        if (xgc->data.set.sp.zres > 0)
            zoff = -xgc->data.set.sp.zres / 2.0f;
    } else {
        xoff = -50;
        yoff = -50;
        zoff = -50;
        //glScalef(0.08f, 0.08f, 0.08f);
    }

    glShadeModel(GL_FLAT);

    /*x axis*/
    glBegin(GL_LINES);
    glColor3f(1., 0., 0.);
    glVertex3f(0, 0, 0);
    glVertex3f(100.0, 0, 0);
    glEnd();
    draw_cone(5, -20, 100.0, 0, 0, 0, 90, 0);

    /*y axis*/
    glBegin(GL_LINES);
    glColor3f(0., 1., 0.);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 100.0, 0);
    glEnd();
    draw_cone(5, 20, 0, 100, 0, 90, 0, 0);

    /*z axis*/
    glBegin(GL_LINES);
    glColor3f(0., 0., 1.);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 100.0);
    glEnd();
    draw_cone(5, -20, 0, 0, 100, 0, 0, 0);

    if (xgc->data.set.sp.xres > 0 && xgc->data.set.sp.yres > 0 && xgc->data.set.sp.zres > 0) {
        /*pool*/
        if (xgc->data.is_pool)
            draw_pp(xoff, yoff, zoff, xgc->data.set.sp.xres + xoff, xgc->data.set.sp.yres + yoff, xgc->data.set.sp.zres + zoff,
                0.1f, 0.6f, 0.2f, 0.1f, 0.6f, 0.2f, 0.1f);
        /*tsf*/
        //if (xgc->data.is_sf && xgc->data.set.ss.sf.filename!=NULL) {
        if (xgc->data.is_sf && xgc->data.set.ss.sf.valid) {
            draw_pp(xoff, yoff, zoff, xgc->data.set.sp.xres + xoff, xgc->data.set.sp.yres + yoff, xgc->data.set.sp.zres + zoff,
                0.6f, 0.3f, 0.0f, 0.6f, 0.3f, 0.0f,
                0.1f);
        }
        //if (xgc->data.is_tsf && (xgc->data.set.ss.tsf.i0 + xgc->data.set.ss.tsf.j0 + xgc->data.set.ss.tsf.k0 + xgc->data.set.ss.tsf.i1 + xgc->data.set.ss.tsf.j1 + xgc->data.set.ss.tsf.k1)!=0) {
        if (xgc->data.is_tsf && xgc->data.set.ss.tsf.valid) {
            draw_pp(xgc->data.set.ss.tsf.box_i0 + xoff, xgc->data.set.ss.tsf.box_j0 + yoff, xgc->data.set.ss.tsf.box_k0 + zoff, xgc->data.set.ss.tsf.box_in + xoff, xgc->data.set.ss.tsf.box_jn + yoff, xgc->data.set.ss.tsf.box_kn + zoff,
                0.65f, 0.3f, 0.0f, 0.65f, 0.3f, 0.0f,
                0.15f);

            /*draw arrow*/
        }
        //if (xgc->data.is_tsff && (xgc->data.set.ss.tsff.i0 + xgc->data.set.ss.tsff.j0 + xgc->data.set.ss.tsff.k0 + xgc->data.set.ss.tsff.i1 + xgc->data.set.ss.tsff.j1 + xgc->data.set.ss.tsff.k1)!=0) {
        if (xgc->data.is_tsff && xgc->data.set.ss.tsff.valid) {
            draw_pp(xgc->data.set.ss.tsff.box_i0 + xoff, xgc->data.set.ss.tsff.box_j0 + yoff, xgc->data.set.ss.tsff.box_k0 + zoff, xgc->data.set.ss.tsff.box_in + xoff, xgc->data.set.ss.tsff.box_jn + yoff, xgc->data.set.ss.tsff.box_kn + zoff,
                0.7f, 0.3f, 0.0f, 0.7f, 0.3f, 0.0f,
                0.15f);
        }
        //if (xgc->data.is_ltsf && (xgc->data.set.ss.ltsf.i0 + xgc->data.set.ss.ltsf.j0 + xgc->data.set.ss.ltsf.k0 + xgc->data.set.ss.ltsf.i1 + xgc->data.set.ss.ltsf.j1 + xgc->data.set.ss.ltsf.k1)!=0) {
        if (xgc->data.is_ltsf && xgc->data.set.ss.ltsf.valid) {
            draw_pp(xgc->data.set.ss.ltsf.box_i0 + xoff, xgc->data.set.ss.ltsf.box_j0 + yoff, xgc->data.set.ss.ltsf.box_k0 + zoff, xgc->data.set.ss.ltsf.box_in + xoff, xgc->data.set.ss.ltsf.box_jn + yoff, xgc->data.set.ss.ltsf.box_kn + zoff,
                0.75f, 0.3f, 0.0f, 0.75f, 0.5f, 0.0f,
                0.15f);
        }
        //if (xgc->data.is_ltsff && (xgc->data.set.ss.ltsff.i0 + xgc->data.set.ss.ltsff.j0 + xgc->data.set.ss.ltsff.k0 + xgc->data.set.ss.ltsff.i1 + xgc->data.set.ss.ltsff.j1 + xgc->data.set.ss.ltsff.k1)!=0) {
        if (xgc->data.is_ltsff && xgc->data.set.ss.ltsff.valid) {
            draw_pp(xgc->data.set.ss.ltsff.box_i0 + xoff, xgc->data.set.ss.ltsff.box_j0 + yoff, xgc->data.set.ss.ltsff.box_k0 + zoff, xgc->data.set.ss.ltsff.box_in + xoff, xgc->data.set.ss.ltsff.box_jn + yoff, xgc->data.set.ss.ltsff.box_kn + zoff,
                0.8f, 0.3f, 0.0f, 0.8f, 0.3f, 0.0f,
                0.15f);
        }
        for (i = 0; i < xgc->data.set.ss.npnts; i++) {
            if (xgc->data.is_psrc[i]) {
                draw_pointsource(xgc->data.set.ss.pnts[i].point_origin_position_i + xoff, xgc->data.set.ss.pnts[i].point_origin_position_j + yoff, xgc->data.set.ss.pnts[i].point_origin_position_k + zoff, 1.0f, 0.8f, 0.0f, 0.5f);
            }
        }

        if (xgc->data.is_bx0)
            draw_pp(xoff, yoff, zoff, xoff, xgc->data.set.sp.yres + yoff, xgc->data.set.sp.zres + zoff, 0.6f, 0.1f, 0.1f, 0.6f, 0.1f, 0.1f, 0.5f);
        if (xgc->data.is_bxn)
            draw_pp(xgc->data.set.sp.xres + xoff, yoff, zoff, xgc->data.set.sp.xres + xoff, xgc->data.set.sp.yres + yoff, xgc->data.set.sp.zres + zoff, 0.6f, 0.1f, 0.1f, 0.6f, 0.1f, 0.1f, 0.5f);
        if (xgc->data.is_by0)
            draw_pp(xoff, yoff, zoff, xgc->data.set.sp.xres + xoff, yoff, xgc->data.set.sp.zres + zoff, 0.6f, 0.1f, 0.1f, 0.6f, 0.1f, 0.1f, 0.5f);
        if (xgc->data.is_byn)
            draw_pp(xoff, xgc->data.set.sp.yres + yoff, zoff, xgc->data.set.sp.yres + xoff, xgc->data.set.sp.yres + yoff, xgc->data.set.sp.zres + zoff, 0.6f, 0.1f, 0.1f, 0.6f, 0.1f, 0.1f, 0.5f);
        if (xgc->data.is_bz0)
            draw_pp(xoff, yoff, zoff, xgc->data.set.sp.xres + xoff, xgc->data.set.sp.yres + yoff, zoff, 0.6f, 0.1f, 0.1f, 0.6f, 0.1f, 0.1f, 0.5f);
        if (xgc->data.is_bzn)
            draw_pp(xoff, yoff, xgc->data.set.sp.zres + zoff, xgc->data.set.sp.xres + xoff, xgc->data.set.sp.yres + yoff, xgc->data.set.sp.zres + zoff, 0.6f, 0.1f, 0.1f, 0.6f, 0.1f, 0.1f, 0.5f);

        if (xgc->data.set.smb.bx0 == SV_BOUNDARY_PERIODIC && xgc->data.is_mbx0)
            draw_pp(xgc->data.set.smb.bx0pos + xoff, yoff, zoff, xgc->data.set.smb.bx0pos + xoff, xgc->data.set.sp.yres + yoff, xgc->data.set.sp.zres + zoff, 0.8f, 0.2f, 0.1f, 0.8f, 0.2f, 0.1f, 0.2f);
        if (xgc->data.set.smb.bxn == SV_BOUNDARY_PERIODIC && xgc->data.is_mbxn)
            draw_pp(xgc->data.set.smb.bxnpos + xoff, yoff, zoff, xgc->data.set.smb.bxnpos + xoff, xgc->data.set.sp.yres + yoff, xgc->data.set.sp.zres + zoff, 0.8f, 0.2f, 0.1f, 0.8f, 0.2f, 0.1f, 0.2f);

        if (xgc->data.set.smb.by0 == SV_BOUNDARY_PERIODIC && xgc->data.is_mby0)
            draw_pp(xoff, xgc->data.set.smb.by0pos + yoff, zoff, xgc->data.set.sp.xres + xoff, xgc->data.set.smb.by0pos + yoff, xgc->data.set.sp.zres + zoff, 0.8f, 0.2f, 0.1f, 0.8f, 0.2f, 0.1f, 0.2f);
        if (xgc->data.set.smb.byn == SV_BOUNDARY_PERIODIC && xgc->data.is_mbyn)
            draw_pp(xoff, xgc->data.set.smb.bynpos + yoff, zoff, xgc->data.set.sp.xres + xoff, xgc->data.set.smb.bynpos + yoff, xgc->data.set.sp.zres + zoff, 0.8f, 0.2f, 0.1f, 0.8f, 0.2f, 0.1f, 0.2f);

        if (xgc->data.set.smb.bz0 == SV_BOUNDARY_PERIODIC && xgc->data.is_mbz0)
            draw_pp(xoff, yoff, xgc->data.set.smb.bz0pos + zoff, xgc->data.set.sp.xres + xoff, xgc->data.set.sp.yres + yoff, xgc->data.set.smb.bz0pos + zoff, 0.8f, 0.2f, 0.1f, 0.8f, 0.2f, 0.1f, 0.2f);
        if (xgc->data.set.smb.bzn == SV_BOUNDARY_PERIODIC && xgc->data.is_mbzn)
            draw_pp(xoff, yoff, xgc->data.set.smb.bznpos + zoff, xgc->data.set.sp.xres + xoff, xgc->data.set.sp.yres + yoff, xgc->data.set.smb.bznpos + zoff, 0.8f, 0.2f, 0.1f, 0.8f, 0.2f, 0.1f, 0.2f);



        for (i = 0; i < xgc->data.set.so.npnts; i++) {
            if (xgc->data.is_outpnt[i])
                draw_outpoint(xgc->data.set.so.pnts[i].i + xoff, xgc->data.set.so.pnts[i].j + yoff, xgc->data.set.so.pnts[i].k + zoff, 0.0f, 0.0f, 1.0f, 0.5f);
        }
        for (i = 0; i < xgc->data.set.so.nplns; i++) {
            if (xgc->data.is_outpln[i]) {
                if (xgc->data.set.so.plns[i].j == -1 && xgc->data.set.so.plns[i].k == -1)
                    draw_pp(xgc->data.set.so.plns[i].i + xoff, yoff, zoff, xgc->data.set.so.plns[i].i + xoff + 1, xgc->data.set.sp.yres + yoff, xgc->data.set.sp.zres + zoff,
                        0.1f, 0.2f, 0.6f, 0.1f, 0.2f, 0.6f, 0.1f);
                if (xgc->data.set.so.plns[i].i == -1 && xgc->data.set.so.plns[i].k == -1)
                    draw_pp(xoff, xgc->data.set.so.plns[i].j + yoff, zoff, xgc->data.set.sp.xres + xoff, xgc->data.set.so.plns[i].j + yoff + 1, xgc->data.set.sp.zres + zoff,
                        0.1f, 0.2f, 0.6f, 0.1f, 0.2f, 0.6f, 0.1f);
                if (xgc->data.set.so.plns[i].i == -1 && xgc->data.set.so.plns[i].j == -1)
                    draw_pp(xoff, yoff, xgc->data.set.so.plns[i].k + zoff, xgc->data.set.sp.xres + xoff, xgc->data.set.sp.yres + yoff, xgc->data.set.so.plns[i].k + zoff + 1,
                        0.1f, 0.2f, 0.6f, 0.1f, 0.2f, 0.6f, 0.1f);
            }
        }
        for (i = 0; i < xgc->data.set.so.nimgs; i++) {
            if (xgc->data.is_outimg[i]) {
                xgc->image_selected = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->image_combo));

                if (xgc->data.set.so.imgs[i].j == -1 && xgc->data.set.so.imgs[i].k == -1) {

                    if (i == xgc->image_selected && gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->image_check))) {
                        draw_pp(xgc->data.set.so.imgs[i].i + xoff, yoff, zoff, xgc->data.set.so.imgs[i].i + xoff, xgc->data.set.sp.yres + yoff, xgc->data.set.sp.zres + zoff,
                            0.1f, 0.3f, 0.7f, 0.1f, 0.3f, 0.7f, 0.01f);
                        draw_datafield(xgc->outputfield[i], (gint)(xgc->data.set.so.imgs[i].i + xoff), (gint)yoff, (gint)zoff, (gint)(xgc->data.set.so.imgs[i].i + xoff), (gint)(xgc->data.set.sp.yres + yoff),
                            (gint)(xgc->data.set.sp.zres + zoff), 0, gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->image_squarecheck)));
                    }
                    else
                        draw_pp(xgc->data.set.so.imgs[i].i + xoff, yoff, zoff, xgc->data.set.so.imgs[i].i + xoff, xgc->data.set.sp.yres + yoff, xgc->data.set.sp.zres + zoff,
                            0.1f, 0.3f, 0.7f, 0.1f, 0.3f, 0.7f, 0.1f);
                }
                if (xgc->data.set.so.imgs[i].i == -1 && xgc->data.set.so.imgs[i].k == -1) {

                    if (i == xgc->image_selected && gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->image_check))) {
                        draw_pp(xoff, xgc->data.set.so.imgs[i].j + yoff, zoff, xgc->data.set.sp.xres + xoff, xgc->data.set.so.imgs[i].j + yoff, xgc->data.set.sp.zres + zoff,
                            0.1f, 0.3f, 0.7f, 0.1f, 0.3f, 0.7f, 0.01f);
                        draw_datafield(xgc->outputfield[i], (gint)xoff, (gint)(xgc->data.set.so.imgs[i].j + yoff), (gint)zoff, (gint)(xgc->data.set.sp.xres + xoff), (gint)(xgc->data.set.so.imgs[i].j + yoff),
                            (gint)(xgc->data.set.sp.zres + zoff), 1, gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->image_squarecheck)));
                    }
                    else
                        draw_pp(xoff, xgc->data.set.so.imgs[i].j + yoff, zoff, xgc->data.set.sp.xres + xoff, xgc->data.set.so.imgs[i].j + yoff, xgc->data.set.sp.zres + zoff,
                            0.1f, 0.3f, 0.7f, 0.1f, 0.3f, 0.7f, 0.1f);
                }
                if (xgc->data.set.so.imgs[i].i == -1 && xgc->data.set.so.imgs[i].j == -1) {

                    if (i == xgc->image_selected && gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->image_check))) {
                        draw_pp(xoff, yoff, xgc->data.set.so.imgs[i].k + zoff, xgc->data.set.sp.xres + xoff, xgc->data.set.sp.yres + yoff, xgc->data.set.so.imgs[i].k + zoff,
                            0.1f, 0.3f, 0.7f, 0.1f, 0.3f, 0.7f, 0.01f);
                        draw_datafield(xgc->outputfield[i], (gint)xoff, (gint)yoff, (gint)(xgc->data.set.so.imgs[i].k + zoff), (gint)(xgc->data.set.sp.xres + xoff), (gint)(xgc->data.set.sp.yres + yoff),
                            (gint)(xgc->data.set.so.imgs[i].k + zoff), 2, gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->image_squarecheck)));
                    }
                    else
                        draw_pp(xoff, yoff, xgc->data.set.so.imgs[i].k + zoff, xgc->data.set.sp.xres + xoff, xgc->data.set.sp.yres + yoff, xgc->data.set.so.imgs[i].k + zoff,
                            0.1f, 0.3f, 0.7f, 0.1f, 0.3f, 0.7f, 0.1f);
                }
            }

        }
        for (i = 0; i < xgc->data.set.so.nsums; i++) {
            if (xgc->data.is_outsum[i]) {
                draw_pp(xgc->data.set.so.sums[i].box_i0 + xoff, xgc->data.set.so.sums[i].box_j0 + yoff, xgc->data.set.so.sums[i].box_k0 + zoff, xgc->data.set.so.sums[i].box_in + xoff, xgc->data.set.so.sums[i].box_jn + yoff, xgc->data.set.so.sums[i].box_kn + zoff,
                    0.1f, 0.2f, 0.8f, 0.1f, 0.2f, 0.8f, 0.1f);
            }
        }
        for (i = 0; i < xgc->data.set.so.nforces; i++) {
            if (xgc->data.is_outforce[i]) {
                draw_pp(xgc->data.set.so.forces[i].box_i0 + xoff, xgc->data.set.so.forces[i].box_j0 + yoff, xgc->data.set.so.forces[i].box_k0 + zoff, xgc->data.set.so.forces[i].box_in + xoff, xgc->data.set.so.forces[i].box_jn + yoff, xgc->data.set.so.forces[i].box_kn + zoff,
                    0.1f, 0.2f, 0.8f, 0.1f, 0.2f, 0.8f, 0.1f);
            }
        }
        if (xgc->data.set.sf.nrs) {
            if (xgc->data.is_nfff)
                draw_pp(xgc->data.set.sf.box_i0 + xoff, xgc->data.set.sf.box_j0 + yoff, xgc->data.set.sf.box_k0 + zoff, xgc->data.set.sf.box_in + xoff, xgc->data.set.sf.box_jn + yoff, xgc->data.set.sf.box_kn + zoff,
                    0.6f, 0.0f, 0.6f, 0.6f, 0.0f, 0.6f, 0.1f);

            for (i = 0; i < xgc->data.set.sf.nrs; i++) {
                if (xgc->data.is_nfff_point[i])
                    draw_line(0, 0, 0, xgc->data.set.sf.ri[i] + xoff, xgc->data.set.sf.rj[i] + yoff, xgc->data.set.sf.rk[i] + zoff, 1.0f, 0.0f, 1.0f, 0.5f);
            }
        }
        if (xgc->data.set.spf.nrs) {
            if (xgc->data.is_pnfff)
                draw_pp(xgc->data.set.spf.box_i0 + xoff, xgc->data.set.spf.box_j0 + yoff, xgc->data.set.spf.box_k0 + zoff, xgc->data.set.spf.box_in + xoff, xgc->data.set.spf.box_jn + yoff, xgc->data.set.spf.box_kn + zoff,
                    0.6f, 0.0f, 0.6f, 0.6f, 0.0f, 0.6f, 0.1f);

            for (i = 0; i < xgc->data.set.spf.nrs; i++) {
                if (xgc->data.is_pnfff_point[i])
                    draw_line(0, 0, 0, xgc->data.set.spf.ri[i] + xoff, xgc->data.set.spf.rj[i] + yoff, xgc->data.set.spf.rk[i] + zoff, 1.0f, 0.1f, 1.0f, 0.5f);
            }
        }


        /*********************************   material  ****************************************/

        for (m = 0; m < xgc->data.set_mat.spheres->len; m++) {
            if (xgc->data.is_sphere[m]) {
                sx = g_array_index(xgc->data.set_mat.spheres, SvSphere, m);
                draw_sphere(sx.radius, sx.pnt1[0] + xoff, sx.pnt1[1] + yoff, sx.pnt1[2] + zoff, 1.0f, 1.0f, 1.0f, 0.5f);
            }

        }
        for (m = 0; m < xgc->data.set_mat.voxels->len; m++) {
            if (xgc->data.is_voxel[m]) {
                vx = g_array_index(xgc->data.set_mat.voxels, SvVoxel, m);
                draw_pp(vx.pnt1[0] + xoff, vx.pnt1[1] + yoff, vx.pnt1[2] + zoff, vx.pnt2[0] + xoff, vx.pnt2[1] + yoff, vx.pnt2[2] + zoff, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.5f);
            }

        }
        for (m = 0; m < xgc->data.set_mat.cylinders->len; m++) {
            if (xgc->data.is_cylinder[m]) {
                cx = g_array_index(xgc->data.set_mat.cylinders, SvCylinder, m);
                draw_cylinder(cx.pnt1[0] + xoff, cx.pnt1[1] + yoff, cx.pnt1[2] + zoff, cx.pnt2[0] + xoff, cx.pnt2[1] + yoff, cx.pnt2[2] + zoff, cx.radius, 1.0f, 1.0f, 1.0f, 0.5f);
            }
        }
        for (m = 0; m < xgc->data.set_mat.cones->len; m++) {
            if (xgc->data.is_cone[m]) {
                cnx = g_array_index(xgc->data.set_mat.cones, SvCone, m);
                draw_conem(cnx.pnt2[0] + xoff, cnx.pnt2[1] + yoff, cnx.pnt2[2] + zoff, cnx.pnt1[0] + xoff, cnx.pnt1[1] + yoff, cnx.pnt1[2] + zoff, cnx.radius, 1.0f, 1.0f, 1.0f, 0.5f);
            }
        }
        for (m = 0; m < xgc->data.set_mat.rcones->len; m++) {
            if (xgc->data.is_rcone[m]) {
                rcnx = g_array_index(xgc->data.set_mat.rcones, SvRCone, m);
                draw_rcone(rcnx.pnt1[0] + xoff, rcnx.pnt1[1] + yoff, rcnx.pnt1[2] + zoff, rcnx.pnt2[0] + xoff, rcnx.pnt2[1] + yoff, rcnx.pnt2[2] + zoff, rcnx.radius1, rcnx.radius2, 1.0f, 1.0f, 1.0f, 0.5f);
            }
        }
        for (m = 0; m < xgc->data.set_mat.tetrahedrons->len; m++) {
            if (xgc->data.is_tetrahedron[m]) {
                tx = g_array_index(xgc->data.set_mat.tetrahedrons, SvTetrahedron, m);
                draw_tetrahedron(tx.pnt1[0] + xoff, tx.pnt1[1] + yoff, tx.pnt1[2] + zoff,
                    tx.pnt2[0] + xoff, tx.pnt2[1] + yoff, tx.pnt2[2] + zoff,
                    tx.pnt3[0] + xoff, tx.pnt3[1] + yoff, tx.pnt3[2] + zoff,
                    tx.pnt4[0] + xoff, tx.pnt4[1] + yoff, tx.pnt4[2] + zoff,
                    1.0f, 1.0f, 1.0f, 0.2f);
            }

        }
        for (m = 0; m < xgc->data.set_mat.gwydds->len; m++) {
            if (xgc->data.is_gwydd[m] && xgc->data.set_mat.gwyddata) {

                gx = g_array_index(xgc->data.set_mat.gwydds, SvGwydd, m);
                //gx.i, gx.j, gx.k, xoffset, yoffset, zoffset, depth


                if (!gx.dfield) printf("gwyddion field not found, skipped\n");
                else {
                    for (i = 0; i < xgc->data.set.sp.xres; i++) { //x
                        for (j = 0; j < xgc->data.set.sp.yres; j++) { //y
                            for (k = 0; k < xgc->data.set.sp.zres; k++) { //z
                                if (xgc->data.set_mat.gwyddata->data[i][j][k] != 0) {
                                    fill_pp(i + xoff, j + yoff, k + zoff, i + 1 + xoff, j + 1 + yoff, k + 1 + zoff,
                                        0.8f,
                                        (GLfloat)(1.0 - xgc->data.set_mat.gwyddata->data[i][j][k]),
                                        (GLfloat)(1.0 - xgc->data.set_mat.gwyddata->data[i][j][k]),
                                        1.0f);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/* main visualization function */
gboolean
scene_expose_cb(GtkWidget *da, GdkEventExpose *event, XGControls *xgc)
{    
    GdkGLContext *glcontext = gtk_widget_get_gl_context(da);
    GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(da);        
    // g_print (" :: expose\n");

    if (!gdk_gl_drawable_gl_begin(gldrawable, glcontext)) {
        g_assert_not_reached();
    }

    /* general drawing settings */
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);    
    glPushMatrix();       

    glScalef(xgc->zoomfactor, xgc->zoomfactor, xgc->zoomfactor);
    glTranslatef((GLfloat)xgc->viewposx, (GLfloat)xgc->viewposy, 0.0f);
    glRotatef(xgc->xrot, 1.0f, 0.0f, 0.0f);
    glRotatef(xgc->yrot, 0.0f, 1.0f, 0.0f);

    render_scene(xgc);    

    glPopMatrix();        

    if (gdk_gl_drawable_is_double_buffered(gldrawable))
        gdk_gl_drawable_swap_buffers(gldrawable);
    else
        glFlush();

    gdk_gl_drawable_gl_end(gldrawable);

    return 1;
}

gboolean
scene_button_press_cb(GtkWidget* widget, GdkEventButton * event, XGControls *xgc)
{
    gfloat x, y;

    if (event->type == GDK_BUTTON_PRESS) {  
        if (event->button == 3) {           /* rotate */
            x = (gfloat)event->x;
            y = (gfloat)event->y;
            xgc->xdiff = x - xgc->yrot;
            xgc->ydiff = -y + xgc->xrot;

            return TRUE;
        }

        if (event->button == 2) {           /* pan */
            glPushMatrix();

            GLdouble modelMatrix[16];
            GLdouble projMatrix[16];
            GLint viewport[4];
            GLdouble objz, winy;
            glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
            glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
            glGetIntegerv(GL_VIEWPORT, viewport);

            winy = viewport[3] - event->y;
            gluUnProject(event->x, winy, 0, modelMatrix, projMatrix, viewport, &xgc->viewposx_start, &xgc->viewposy_start, &objz);

            glPopMatrix();

            return TRUE;
        }        
    } else if (event->type == GDK_2BUTTON_PRESS) {
        if (event->button == 3) {           /* fit */
            scene_reset_cb(widget, xgc);
            return TRUE;
        }
    }


    return FALSE;
}

gboolean
scene_button_release_cb(GtkWidget* widget, GdkEventButton * event, XGControls *xgc)
{
    if (event->type == GDK_BUTTON_RELEASE) {
        if (event->button == 2) {           /* pan */
            xgc->viewposx_prev = xgc->viewposx;
            xgc->viewposy_prev = xgc->viewposy;

            return TRUE;
        }
    }

    return FALSE;
}

get_move_to_cursor_offset(gdouble mouse_pos_x, gdouble mouse_pos_y, 
                          XGControls *xgc, gdouble* offset_x, gdouble* offset_y)
{
    glPushMatrix();

    GLdouble modelMatrix[16];
    GLdouble projMatrix[16];
    GLint viewport[4];
    GLdouble objx, objy, objz;
    glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
    glGetIntegerv(GL_VIEWPORT, viewport);

    GLdouble winy = viewport[3] - mouse_pos_y;

    gluUnProject(mouse_pos_x, winy, 0,
                 modelMatrix, projMatrix, viewport,
                 &objx, &objy, &objz);

    glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
    glGetIntegerv(GL_VIEWPORT, viewport);

    *offset_x = xgc->viewposx_prev + (objx - xgc->viewposx_start)/(xgc->zoomfactor);
    *offset_y = xgc->viewposy_prev + (objy - xgc->viewposy_start)/(xgc->zoomfactor);

    glPopMatrix();  
}

void
get_zoom_to_cursor_offset(gdouble mouse_pos_x, gdouble mouse_pos_y, XGControls *xgc, 
                          gdouble* offset_x, gdouble* offset_y)
{
    GLdouble modelMatrix[16];
    GLdouble projMatrix[16];
    GLint viewport[4];
    GLdouble objx, objy, objz, objx1, objy1, objz1, winy;

    /* scale to previous zoom factor and translate to actual view position (x, y, z) */
    glPushMatrix();

    glScalef(xgc->zoomfactor_prev, xgc->zoomfactor_prev, xgc->zoomfactor_prev);
    glTranslatef((GLfloat)xgc->viewposx, (GLfloat)xgc->viewposy, 0.0f);

    glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
    glGetIntegerv(GL_VIEWPORT, viewport);

    winy = viewport[3] - mouse_pos_y;

    gluUnProject(mouse_pos_x, winy, 0,
                 modelMatrix, projMatrix, viewport,
                 &objx, &objy, &objz);

    glPopMatrix();


    /* scale to actual zoom factor and translate to actual position */
    glPushMatrix();

    glScalef(xgc->zoomfactor, xgc->zoomfactor, xgc->zoomfactor);
    glTranslatef((GLfloat)xgc->viewposx, (GLfloat)xgc->viewposy, 0.0f);

    glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
    glGetIntegerv(GL_VIEWPORT, viewport);    

    gluUnProject(mouse_pos_x, winy, 0,
                 modelMatrix, projMatrix, viewport,
                 &objx1, &objy1, &objz1);
    
    glPopMatrix();
    

    *offset_x = objx1 - objx;
    *offset_y = objy1 - objy;
}

gboolean
scene_scroll_cb(GtkWidget* widget, GdkEventScroll * event, XGControls *xgc)
{
    gdouble offset_x, offset_y;

    offset_x = offset_y = 0.0f;

    if (event->type == GDK_SCROLL) {
        if (event->direction == GDK_SCROLL_UP) {
            xgc->zoomfactor_prev = xgc->zoomfactor;
            xgc->zoomfactor *= 1.0f + ZOOM_STEP;                        

            get_zoom_to_cursor_offset(event->x, event->y, xgc, &offset_x, &offset_y);
            xgc->viewposx = xgc->viewposx_prev + offset_x;
            xgc->viewposy = xgc->viewposy_prev + offset_y;            

            xgc->viewposx_prev = xgc->viewposx;
            xgc->viewposy_prev = xgc->viewposy;

            gtk_widget_queue_draw(xgc->view_scene);

            return TRUE;
        } else if (event->direction == GDK_SCROLL_DOWN) {
            xgc->zoomfactor_prev = xgc->zoomfactor;
            xgc->zoomfactor *= 1.0f - ZOOM_STEP;

            get_zoom_to_cursor_offset(event->x, event->y, xgc, &offset_x, &offset_y);
            xgc->viewposx = xgc->viewposx_prev + offset_x;
            xgc->viewposy = xgc->viewposy_prev + offset_y;            

            xgc->viewposx_prev = xgc->viewposx;
            xgc->viewposy_prev = xgc->viewposy;

            gtk_widget_queue_draw(xgc->view_scene);

            return TRUE;
        }
    } 

    return FALSE;
}


gboolean
scene_motion_notify_cb(GtkWidget *widget, GdkEventMotion *event, XGControls *xgc)
{
    gint x, y;
    gdouble offset_x, offset_y;
    GdkModifierType state;

    if (event->is_hint)
        gdk_window_get_pointer(event->window, &x, &y, &state);
    else {
        x = (gint)event->x;
        y = (gint)event->y;
        state = event->state;
    }

    if (state & GDK_BUTTON3_MASK) {                 /* rotate*/
        xgc->yrot = x - xgc->xdiff;
        xgc->xrot = y + xgc->ydiff;
        gtk_widget_queue_draw(xgc->view_scene);
    } else if (state & GDK_BUTTON2_MASK) {          /* pan */        
        offset_x = offset_y = 0.0f;
        get_move_to_cursor_offset(x, y, xgc, &offset_x, &offset_y);
        xgc->viewposx = offset_x;
        xgc->viewposy = offset_y;

        gtk_widget_queue_draw(xgc->view_scene);
    }

    return TRUE;
}

void
scene_zoom_in_cb(GtkWidget *widget, XGControls *xgc)
{
    xgc->zoomfactor *= 1.0f + ZOOM_STEP;
    gtk_widget_queue_draw(xgc->view_scene);
}

void
scene_zoom_out_cb(GtkWidget *widget, XGControls *xgc)
{
    xgc->zoomfactor *= 1.0f - ZOOM_STEP;
    gtk_widget_queue_draw(xgc->view_scene);
}

void
scene_zoom_fit_cb(GtkWidget *widget, XGControls *xgc)
{
    //xgc->zoomfactor = 5.0f;
    xgc->zoomfactor = 1.0f;
    xgc->viewposx = 0.0f;
    xgc->viewposy = 0.0f;
    xgc->xrot = 0.0f;
    xgc->yrot = 0.0f;
    gtk_widget_queue_draw(xgc->view_scene);
}

void
scene_move_view_left_cb(GtkWidget *widget, XGControls *xgc)
{
    //    xgc->timeout = g_timeout_add(50, move_right_timeout, xgc);
    xgc->viewposx += 30;

    xgc->viewposx_prev = xgc->viewposx;
    xgc->viewposy_prev = xgc->viewposy;

    gtk_widget_queue_draw(xgc->view_scene);

    //fclose(fr);
}

void
scene_move_view_right_cb(GtkWidget *widget, XGControls *xgc)
{
    xgc->viewposx -= 30;

    xgc->viewposx_prev = xgc->viewposx;
    xgc->viewposy_prev = xgc->viewposy;

    gtk_widget_queue_draw(xgc->view_scene);
}

void
scene_move_view_up_cb(GtkWidget *widget, XGControls *xgc)
{
    xgc->viewposy += 30;

    xgc->viewposx_prev = xgc->viewposx;
    xgc->viewposy_prev = xgc->viewposy;

    gtk_widget_queue_draw(xgc->view_scene);
}

void
scene_move_view_down_cb(GtkWidget *widget, XGControls *xgc)
{
    xgc->viewposy -= 30;

    xgc->viewposx_prev = xgc->viewposx;
    xgc->viewposy_prev = xgc->viewposy;

    gtk_widget_queue_draw(xgc->view_scene);
}

/* end of scene callback functions */
//////////////////////////////////////////////////////////////////////////

