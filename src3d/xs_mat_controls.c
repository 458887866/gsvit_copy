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


/* xs_mat_controls.c :
*  Material object controls functions
*/

#include "xs_mat_controls.h"
#include "xs_mat_tree.h"
#include "write.h"
#include <string.h>

void
mat_controls_sensitive(XGControls *xgc, gint id)
{
    SvSphere        *sphere;
    SvVoxel         *voxel;
    SvCylinder      *cyl;
    SvCone          *cone;
    SvRCone         *rcone;
    SvTetrahedron   *tthn;
    SvMatProp       mat;    

    memset(&mat, 0, sizeof(SvMatProp));

    if (id >= SET_MAT_SPHERE && id < SET_MAT_VOXEL) {
        sphere = &g_array_index(xgc->data.set_mat.spheres, SvSphere, id - SET_MAT_SPHERE);
        memcpy(&mat, &sphere->mat, sizeof(SvMatProp));
    } else if (id >= SET_MAT_VOXEL && id < SET_MAT_CYLINDER) {
        voxel = &g_array_index(xgc->data.set_mat.voxels, SvVoxel, id - SET_MAT_VOXEL);
        memcpy(&mat, &voxel->mat, sizeof(SvMatProp));
    } else if (id >= SET_MAT_CYLINDER && id < SET_MAT_CONE) {
        cyl = &g_array_index(xgc->data.set_mat.cylinders, SvCylinder, id - SET_MAT_CYLINDER);
        memcpy(&mat, &cyl->mat, sizeof(SvMatProp));
    } else if (id >= SET_MAT_CONE && id < SET_MAT_RCONE) {
        cone = &g_array_index(xgc->data.set_mat.cones, SvCone, id - SET_MAT_CONE);
        memcpy(&mat, &cone->mat, sizeof(SvMatProp));
    } else if (id >= SET_MAT_RCONE && id < SET_MAT_GWYDD) {
        rcone = &g_array_index(xgc->data.set_mat.rcones, SvRCone, id - SET_MAT_RCONE);
        memcpy(&mat, &rcone->mat, sizeof(SvMatProp));
    } else if (id >= SET_MAT_GWYDD && id < SET_MAT_MESH) {
    } else if (id >= SET_MAT_MESH && id < SET_MAT_TETRAHEDRON) {
    } else if (id >= SET_MAT_TETRAHEDRON) {
        tthn = &g_array_index(xgc->data.set_mat.tetrahedrons, SvTetrahedron, id - SET_MAT_TETRAHEDRON);
        memcpy(&mat, &tthn->mat, sizeof(SvMatProp));
    }

      

    gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs_mat.oc.epsilon_spin), (mat.type == 0));
    gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs_mat.oc.mu_spin), (mat.type == 0));
    gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs_mat.oc.sigma_spin), (mat.type == 0));
    gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs_mat.oc.sigast_spin), (mat.type == 0));

    gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs_mat.oc.tepsilon_spin), (mat.type == 1));
    gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs_mat.oc.tmu_spin), (mat.type == 1));
    gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs_mat.oc.tsigma_spin), (mat.type == 1));
    gtk_widget_set_sensitive(GTK_WIDGET(xgc->cs_mat.oc.tsigast_spin), (mat.type == 1));
}

void
mat_controls_changed(XGControls *xgc)
{
    if (xgc->mat_lock_controls_update) {
        return;
    }

//    gint i;
//    guint id, ugval;
    gchar buff[256];
    GtkWidget *focus;
    gchar size_um_buff[256] = {0};
    SvMatProp *pmat = NULL;
    const gchar *overriden = NULL;
    gboolean is_overriden = FALSE;
    gint posid = 0;
    gint i, index;
    gboolean getmat = FALSE;
    SvSphere* sphere;
    SvVoxel* voxel;
    SvCylinder* cyl;
    SvCone* cone;
    SvRCone* rcone;
    SvTetrahedron* tthn;    

    // remember focus at the beginning and restore at the end of ss_controls_changed()
    // Note: this works only for GtkEntry not for GtkSpinButton
    focus = gtk_window_get_focus(GTK_WINDOW(xgc->toplevel));

    if (xgc->mat_selid >= SET_MAT_SPHERE && xgc->mat_selid < SET_MAT_VOXEL) {
        /* I. 1. Sphere properties */
        index = xgc->mat_selid - SET_MAT_SPHERE;
        sphere = &g_array_index(xgc->data.set_mat.spheres, SvSphere, index);

        sphere->pnt1[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.sphere_center_i_spin));
        sphere->pnt1[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.sphere_center_j_spin));
        sphere->pnt1[2] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.sphere_center_k_spin));
        sphere->radius = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.sphere_radius_spin));        

        //sphere->mat.material = g_strdup(gtk_entry_get_text(GTK_ENTRY(xgc->cs_mat.oc.overriden)));
        //overriden = gtk_entry_get_text(GTK_ENTRY(xgc->cs_mat.oc.overriden));
        overriden = gtk_combo_box_get_active_text(GTK_COMBO_BOX(xgc->cs_mat.oc.overriden));

        pmat = &(sphere->mat);
        posid = sphere->n;

        /*if (overriden) { //remove overriden entry if there was any (point it nowhere)
            for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                if (xgc->data.set_mat.overpos[i] == posid)
                    xgc->data.set_mat.overpos[i] = -1;
                    //xgc->data.set_mat.valid[i] = FALSE;                    
            }
        }*/

        //overriden = xgc->data.set_mat.materials[posid];
    } else if (xgc->mat_selid >= SET_MAT_VOXEL && xgc->mat_selid < SET_MAT_CYLINDER) {
        /* I. 2. Box properties */
        index = /*xgc->data.set_mat.spheres->len +*/ (xgc->mat_selid - SET_MAT_VOXEL);
        voxel = &g_array_index(xgc->data.set_mat.voxels, SvVoxel, index);

        voxel->pnt1[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.box_i0_spin));
        voxel->pnt1[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.box_j0_spin));
        voxel->pnt1[2] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.box_k0_spin));
        voxel->pnt2[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.box_in_spin));
        voxel->pnt2[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.box_jn_spin));
        voxel->pnt2[2] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.box_kn_spin));        

        pmat = &(voxel->mat);
        posid = xgc->data.set_mat.spheres->len + voxel->n;
    } else if (xgc->mat_selid >= SET_MAT_CYLINDER && xgc->mat_selid < SET_MAT_CONE) {
        /* I. 3. Cylinder properties */
        index = xgc->mat_selid - SET_MAT_CYLINDER;
        cyl = &g_array_index(xgc->data.set_mat.cylinders, SvCylinder, index);

        cyl->pnt1[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cyl_i0_spin));
        cyl->pnt1[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cyl_j0_spin));
        cyl->pnt1[2] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cyl_k0_spin));
        cyl->pnt2[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cyl_in_spin));
        cyl->pnt2[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cyl_jn_spin));
        cyl->pnt2[2] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cyl_kn_spin));
        cyl->radius = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cyl_radius_spin));

        pmat = &(cyl->mat);
        posid = xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + cyl->n;
    } else if (xgc->mat_selid >= SET_MAT_CONE && xgc->mat_selid < SET_MAT_RCONE) {
        /* I. 3. Cone properties */
        index = xgc->mat_selid - SET_MAT_CONE;
        cone = &g_array_index(xgc->data.set_mat.cones, SvCone, index);

        cone->pnt1[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cone_i0_spin));
        cone->pnt1[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cone_j0_spin));
        cone->pnt1[2] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cone_k0_spin));
        cone->pnt2[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cone_in_spin));
        cone->pnt2[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cone_jn_spin));
        cone->pnt2[2] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cone_kn_spin));
        cone->radius = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cone_radius_spin));

        pmat = &(cone->mat);
        posid = xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + xgc->data.set_mat.cylinders->len + cone->n;
    } else if (xgc->mat_selid >= SET_MAT_RCONE && xgc->mat_selid < SET_MAT_GWYDD) {
        /* I. 4. Cut cone properties */
        index = xgc->mat_selid - SET_MAT_RCONE;
        rcone = &g_array_index(xgc->data.set_mat.rcones, SvRCone, index);

        rcone->pnt1[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_i0_spin));
        rcone->pnt1[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_j0_spin));
        rcone->pnt1[2] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_k0_spin));
        rcone->pnt2[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_in_spin));
        rcone->pnt2[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_jn_spin));
        rcone->pnt2[2] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_kn_spin));
        rcone->radius1 = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_radius1_spin));
        rcone->radius2 = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_radius2_spin));

        pmat = &(rcone->mat);
        posid = xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + xgc->data.set_mat.cylinders->len + xgc->data.set_mat.cones->len + rcone->n;
    } else if (xgc->mat_selid >= SET_MAT_GWYDD && xgc->mat_selid < SET_MAT_MESH) {
    } else if (xgc->mat_selid >= SET_MAT_MESH && xgc->mat_selid < SET_MAT_TETRAHEDRON) {
    } else if (xgc->mat_selid >= SET_MAT_TETRAHEDRON) {
        /* I. 5. Tetrahedron properties */
        index = /*xgc->data.set_mat.spheres->len + xgc->data.set_mat.cylinders->len + xgc->data.set_mat.cones->len + xgc->data.set_mat.rcones->len + xgc->data.set_mat.tetrahedrons->len +*/ (xgc->mat_selid - SET_MAT_TETRAHEDRON);
        tthn = &g_array_index(xgc->data.set_mat.tetrahedrons, SvTetrahedron, index);

        tthn->pnt1[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_i0_spin));
        tthn->pnt1[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_j0_spin));
        tthn->pnt1[2] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_k0_spin));
        tthn->pnt2[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_i1_spin));
        tthn->pnt2[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_j1_spin));
        tthn->pnt2[2] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_k1_spin));
        tthn->pnt3[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_i2_spin));
        tthn->pnt3[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_j2_spin));
        tthn->pnt3[2] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_k2_spin));
        tthn->pnt4[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_i3_spin));
        tthn->pnt4[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_j3_spin));
        tthn->pnt4[2] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_k3_spin));
        
        pmat = &(tthn->mat);
        posid = xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + xgc->data.set_mat.cylinders->len + xgc->data.set_mat.cones->len + xgc->data.set_mat.rcones->len + tthn->n;
    }

    //////////////////////////////////////////////////////////////////////////
    /*material handling common to all the objects*/    
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm0))) {
        pmat->type = 0;
        pmat->epsilon = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.epsilon_spin));
        pmat->mu = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.mu_spin));
        pmat->sigma = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.sigma_spin));
        pmat->sigast = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.sigast_spin));

        if (overriden) { //remove overriden entry if there was any (point it nowhere)
            for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                if (xgc->data.set_mat.overpos[i] == posid)
                    xgc->data.set_mat.overpos[i] = -1;
                    //xgc->data.set_mat.valid[i] = FALSE;
            }
        }
    } else if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm1))) {
        pmat->type = 1;
        pmat->epsilon = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tepsilon_spin));
        pmat->mu = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tmu_spin));
        pmat->sigma = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tsigma_spin));
        pmat->sigast = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tsigast_spin));

        /*if (overriden) { //remove overriden entry if there was any (point it nowhere)
            for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                if (xgc->data.set_mat.overpos[i] == posid)
                    xgc->data.set_mat.overpos[i] = -1;
            }
        }*/
    } else if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm4))) {
        getmat = TRUE;
        pmat->epsilon = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cepsilon_spin));
        pmat->drude_omega_p = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.omegap_spin));
        pmat->drude_nu = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.nu_spin));

        pmat->cp3_a[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.a_spin[0]));
        pmat->cp3_phi[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.ia_phi_spin[0]));
        pmat->cp3_omega[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.omega_spin[0]));
        pmat->cp3_gamma[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.gamma_spin[0]));

        pmat->cp3_a[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.a_spin[1]));
        pmat->cp3_phi[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.ia_phi_spin[1]));
        pmat->cp3_omega[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.omega_spin[1]));
        pmat->cp3_gamma[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.gamma_spin[1]));


        /*if (overriden) { //remove overriden entry if there was any (point it nowhere)
            for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                if (xgc->data.set_mat.overpos[i] == posid)
                    xgc->data.set_mat.overpos[i] = -1;
            }
        }*/
    } else if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm10))) {
        pmat->type = 10;

        /*if (overriden) { //remove overriden entry if there was any (point it nowhere)
            for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                if (xgc->data.set_mat.overpos[i] == posid)
                    xgc->data.set_mat.overpos[i] = -1;
            }
        }*/
    } else if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm99))) {
        pmat->type = 1;

        //overriden = gtk_entry_get_text(GTK_ENTRY(xgc->cs_mat.oc.overriden));
        overriden = gtk_combo_box_get_active_text(GTK_COMBO_BOX(xgc->cs_mat.oc.overriden));
        is_overriden = FALSE;

        /*for (i = 0; i < xgc->data.set_mat.nmats; i++) {
            if (xgc->data.set_mat.overpos[i] == posid) {                
                g_free(xgc->data.set_mat.materials[i]);
                if (strcmp(overriden, MAT_PROP_MATERIAL_NONE) != 0) {
                    xgc->data.set_mat.materials[i] = g_strdup(overriden);
                    pmat->type = 99;
                    is_overriden = TRUE;
                } else {
                    xgc->data.set_mat.overpos[i] = -1;
                }
            }
        }
        if (!is_overriden) {
            xgc->data.set_mat.materials[xgc->data.set_mat.nmats] = g_strdup(overriden);
            xgc->data.set_mat.overpos[i] = posid;
            //xgc->data.set_mat.valid[i] = TRUE;
            xgc->data.set_mat.nmats++;
            pmat->type = 99;
        }*/

        
        if (strcmp(overriden, MAT_PROP_MATERIAL_NONE) != 0) {
            for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                if (xgc->data.set_mat.overpos[i] == posid) {
                    g_free(xgc->data.set_mat.materials[i]);
                    xgc->data.set_mat.materials[i] = g_strdup(overriden);
                    pmat->type = 99;
                } else {
                    xgc->data.set_mat.overpos[i] = -1;
                }
            }
            if (xgc->data.set_mat.nmats == 0) {
                xgc->data.set_mat.materials[xgc->data.set_mat.nmats] = g_strdup(overriden);
                xgc->data.set_mat.overpos[i] = posid;
                xgc->data.set_mat.nmats++;
                pmat->type = 99;
            }
        }

        /*if (0) {
            xgc->data.set_mat.materials[xgc->data.set_mat.nmats] = g_strdup(gtk_entry_get_text(GTK_ENTRY(xgc->cs_mat.oc.overriden)));
            xgc->data.set_mat.overpos[xgc->data.set_mat.nmats] = xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + xgc->data.set_mat.cylinders->len
                + xgc->data.set_mat.cones->len + xgc->data.set_mat.rcones->len + xgc->data.set_mat.tetrahedrons->len + xgc->data.set_mat.gwydds->len - 1;
            xgc->data.set_mat.nmats++;
        } else {
            if (overriden) { //there was already an overriden material entry for this object
                for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                    if (xgc->data.set_mat.overpos[i] == posid) {
                        xgc->data.set_mat.materials[i] = g_strdup(gtk_entry_get_text(GTK_ENTRY(xgc->cs_mat.oc.overriden)));
                    }
                }
                pmat->type = 99;
            } else { //add new overriden material entry
                xgc->data.set_mat.materials[xgc->data.set_mat.nmats] = g_strdup(gtk_entry_get_text(GTK_ENTRY(xgc->cs_mat.oc.overriden)));
                xgc->data.set_mat.overpos[xgc->data.set_mat.nmats] = posid;
                xgc->data.set_mat.nmats++;
            }
        }*/
    }

    if (getmat) {
        if (strcmp(gtk_combo_box_get_active_text(GTK_COMBO_BOX(xgc->cs_mat.oc.metal_model)), "PLRC") == 0)
            pmat->type = 6;
        else if (strcmp(gtk_combo_box_get_active_text(GTK_COMBO_BOX(xgc->cs_mat.oc.metal_model)), "ADE") == 0)
            pmat->type = 5;
        else
            pmat->type = 4;
    }
    //////////////////////////////////////////////////////////////////////////

    mat_controls_sensitive(xgc, xgc->mat_selid);

    GtkTreeIter iter_parent, iter;
    GtkTreePath *path;
    gint leaf_pos = 0;
    gint node_id = SET_MAT_SPHERE;

    if (xgc->mat_selid >= SET_MAT_SPHERE && xgc->mat_selid < SET_MAT_VOXEL) {
        node_id = SET_MAT_SPHERE;
        leaf_pos = xgc->mat_selid - SET_MAT_SPHERE;
    } else if (xgc->mat_selid >= SET_MAT_VOXEL && xgc->mat_selid < SET_MAT_CYLINDER) {
        node_id = SET_MAT_VOXEL;
        leaf_pos = xgc->data.set_mat.spheres->len + (xgc->mat_selid - SET_MAT_VOXEL);
    } else if (xgc->mat_selid >= SET_MAT_CYLINDER && xgc->mat_selid < SET_MAT_CONE) {
        node_id = SET_MAT_CYLINDER;
        leaf_pos = xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + (xgc->mat_selid - SET_MAT_CYLINDER);
    } else if (xgc->mat_selid >= SET_MAT_CONE && xgc->mat_selid < SET_MAT_RCONE) {
        node_id = SET_MAT_CONE;
        leaf_pos = xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + xgc->data.set_mat.cylinders->len + (xgc->mat_selid - SET_MAT_CONE);
    } else if (xgc->mat_selid >= SET_MAT_RCONE && xgc->mat_selid < SET_MAT_GWYDD) {
        node_id = SET_MAT_RCONE;
        leaf_pos = xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + xgc->data.set_mat.cylinders->len + xgc->data.set_mat.cones->len + (xgc->mat_selid - SET_MAT_RCONE);
    } else if (xgc->mat_selid >= SET_MAT_GWYDD && xgc->mat_selid < SET_MAT_MESH) {
    } else if (xgc->mat_selid >= SET_MAT_MESH && xgc->mat_selid < SET_MAT_TETRAHEDRON) {
    } else if (xgc->mat_selid >= SET_MAT_TETRAHEDRON) {
        node_id = SET_MAT_TETRAHEDRON;
        leaf_pos = xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + xgc->data.set_mat.cylinders->len + xgc->data.set_mat.cones->len + xgc->data.set_mat.rcones->len + (xgc->mat_selid - SET_MAT_TETRAHEDRON);
    }

    path = gtk_tree_path_new_from_indices(ROW_MAT_OBJECT_INDEX, -1);
    gtk_tree_model_get_iter(GTK_TREE_MODEL(xgc->ts_mat), &iter_parent, path);

    mt_tree_assemble_row_text(xgc, buff, sizeof(buff), node_id, index);
    gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(xgc->ts_mat), &iter, &iter_parent, leaf_pos);
    gtk_tree_store_set(xgc->ts_mat, &iter, COLUMN_PARAMETER, buff, -1);


    gtk_widget_queue_draw(xgc->view_scene);

    // remember focus at the beginning and restore at the end of ss_controls_changed()
    // Note: this works only for GtkEntry not for GtkSpinButton
    if (NULL != focus)
        gtk_window_set_focus(GTK_WINDOW(xgc->toplevel), focus);
} /* mat_controls_changed() */

static gint
gtk_combo_box_get_text_index(GtkComboBox* combo, gchar* text)
{
    GtkTreeModel *model = gtk_combo_box_get_model(combo);
    GtkTreeIter iter;
    gchar *buf;
    gint i = 0;

    if (gtk_tree_model_get_iter_first(model, &iter)) {
        do {
            gtk_tree_model_get(model, &iter, 0, &buf, -1);
            if(strcmp(buf, text) == 0) {                
                return i;
            }
            i++;
        } while (gtk_tree_model_iter_next(model, &iter));
    }

    return -1;
}

static gboolean
set_values_object(XGControls* xgc, SvSetMat *set, SvTypeMat type, gint id)
{
    SvSphere sphere;
    SvVoxel voxel;
    SvCylinder cyl;
    SvCone cone;
    SvRCone rcone;
    SvTetrahedron tthn;
    SvMatProp mat;
    gint metalmodelval = MAT_PROP_METAL_MODEL;
    gint mattype = MAT_PROP_TYPE;
    gint i;
    gchar *overriden = NULL;
    gint posid = 0;

    memset(&mat, 0, sizeof(SvMatProp));
    mat.type = MAT_PROP_TYPE;
    mat.epsilon = MAT_PROP_EPS;
    mat.mu = MAT_PROP_MU;
    mat.sigma = MAT_PROP_SIGMA;
    mat.sigast = MAT_PROP_SIGAST;
    mat.drude_omega_p = MAT_PROP_OMEGA;
    mat.drude_nu = MAT_PROP_NU;
    mat.cp3_a[3] = MAT_PROP_A;
    mat.cp3_phi[3] = MAT_PROP_PHI;
    mat.cp3_omega[3] = MAT_PROP_OMEGA;
    mat.cp3_gamma[3] = MAT_PROP_GAMMA;
//    mat.material = MAT_PROP_MATERIAL_NONE;
    mat.pos = MAT_PROP_POS;    

    //sphere.pnt1[0] = MAT_SPHERE_CENTER_I;
    //sphere.pnt1[1] = MAT_SPHERE_CENTER_J;
    //sphere.pnt1[2] = MAT_SPHERE_CENTER_K;
    //sphere.radius = MAT_SPHERE_RADIUS;
    //memcpy(&sphere.mat, &mat, sizeof(SvMatProp));

    if (type == SV_TYPE_MAT_SPHERE) {
        /* I. 1. Sphere properties */
        sphere = g_array_index(xgc->data.set_mat.spheres, SvSphere, id - SET_MAT_SPHERE);
    } else if (type == SV_TYPE_MAT_VOXEL) {
        /* I. 2. Box properties */
        voxel = g_array_index(xgc->data.set_mat.voxels, SvVoxel, id - SET_MAT_VOXEL);
    } else if (type == SV_TYPE_MAT_CYLINDER) {
        /* I. 3. Cylinder properties */
        cyl = g_array_index(xgc->data.set_mat.cylinders, SvCylinder, id - SET_MAT_CYLINDER);
    } else if (type == SV_TYPE_MAT_CONE) {
        /* I. 4. Cone properties */
        cone = g_array_index(xgc->data.set_mat.cones, SvCone, id - SET_MAT_CONE);
    } else if (type == SV_TYPE_MAT_RCONE) {
        /* I. 5. Cut cone properties */
        rcone = g_array_index(xgc->data.set_mat.rcones, SvRCone, id - SET_MAT_RCONE);
    } else if (type == SV_TYPE_MAT_TETRAHEDRON) {
        /* I. 6. Tetrahedron properties */
        tthn = g_array_index(xgc->data.set_mat.tetrahedrons, SvTetrahedron, id - SET_MAT_TETRAHEDRON);
    }

    xgc->mat_lock_controls_update = TRUE;

    gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs_mat.oc.metal_model), metalmodelval);

    if (type == SV_TYPE_MAT_SPHERE) {
        /* I. 1. Sphere properties */
        memcpy(&mat, &sphere.mat.type, sizeof(SvMatProp));
        posid = sphere.n;

        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.sphere_center_i_spin), sphere.pnt1[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.sphere_center_j_spin), sphere.pnt1[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.sphere_center_k_spin), sphere.pnt1[2]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.sphere_radius_spin), sphere.radius);

        //gtk_entry_set_text(GTK_ENTRY(xgc->cs_mat.oc.overriden), MAT_PROP_MATERIAL_NONE);
        //gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs_mat.oc.overriden), MAT_PROP_MATERIAL_NONE_INDEX);
        if (xgc->data.set_mat.overpos[posid] != -1) {
            overriden = xgc->data.set_mat.materials[posid];
            //gtk_entry_set_text(GTK_ENTRY(xgc->cs_mat.oc.overriden), overriden);
            i = gtk_combo_box_get_text_index(GTK_COMBO_BOX(xgc->cs_mat.oc.overriden), overriden);
            if (i != -1)
                gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs_mat.oc.overriden), i);
            else
                gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs_mat.oc.overriden), MAT_PROP_MATERIAL_NONE_INDEX);
        }
    } else if (type == SV_TYPE_MAT_VOXEL) {
        /* I. 2. Box properties */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.box_i0_spin), voxel.pnt1[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.box_j0_spin), voxel.pnt1[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.box_k0_spin), voxel.pnt1[2]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.box_in_spin), voxel.pnt2[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.box_jn_spin), voxel.pnt2[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.box_kn_spin), voxel.pnt2[2]);

        memcpy(&mat, &voxel.mat.type, sizeof(SvMatProp));
        posid = xgc->data.set_mat.spheres->len + voxel.n;
    } else if (type == SV_TYPE_MAT_CYLINDER) {
        /* I. 3. Cylinder properties */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cyl_i0_spin), cyl.pnt1[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cyl_j0_spin), cyl.pnt1[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cyl_k0_spin), cyl.pnt1[2]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cyl_in_spin), cyl.pnt2[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cyl_jn_spin), cyl.pnt2[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cyl_kn_spin), cyl.pnt2[2]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cyl_radius_spin), cyl.radius);

        memcpy(&mat, &cyl.mat.type, sizeof(SvMatProp));
        posid = xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + cyl.n;
    } else if (type == SV_TYPE_MAT_CONE) {
        /* I. 4. Cone properties */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cone_i0_spin), cone.pnt1[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cone_j0_spin), cone.pnt1[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cone_k0_spin), cone.pnt1[2]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cone_in_spin), cone.pnt2[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cone_jn_spin), cone.pnt2[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cone_kn_spin), cone.pnt2[2]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cone_radius_spin), cone.radius);

        memcpy(&mat, &cone.mat.type, sizeof(SvMatProp));
        posid = xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + xgc->data.set_mat.cylinders->len + cone.n;
    } else if (type == SV_TYPE_MAT_RCONE) {
        /* I. 4. Cut cone properties */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_i0_spin), rcone.pnt1[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_j0_spin), rcone.pnt1[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_k0_spin), rcone.pnt1[2]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_in_spin), rcone.pnt2[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_jn_spin), rcone.pnt2[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_kn_spin), rcone.pnt2[2]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_radius1_spin), rcone.radius1);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.rcone_radius2_spin), rcone.radius2);

        memcpy(&mat, &rcone.mat.type, sizeof(SvMatProp));
        posid = xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + xgc->data.set_mat.cylinders->len + xgc->data.set_mat.cones->len + rcone.n;
    } else if (type == SV_TYPE_MAT_TETRAHEDRON) {
        /* I. 5. Tetrahedron properties */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_i0_spin), tthn.pnt1[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_j0_spin), tthn.pnt1[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_k0_spin), tthn.pnt1[2]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_i1_spin), tthn.pnt2[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_j1_spin), tthn.pnt2[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_k1_spin), tthn.pnt2[2]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_i2_spin), tthn.pnt3[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_j2_spin), tthn.pnt3[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_k2_spin), tthn.pnt3[2]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_i3_spin), tthn.pnt4[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_j3_spin), tthn.pnt4[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tthn_k3_spin), tthn.pnt4[2]);

        memcpy(&mat, &rcone.mat.type, sizeof(SvMatProp));
        posid = xgc->data.set_mat.spheres->len + xgc->data.set_mat.voxels->len + xgc->data.set_mat.cylinders->len + xgc->data.set_mat.cones->len + rcone.n;
    }

    mattype = mat.type;

    //////////////////////////////////////////////////////////////////////////
    /*material handling common to all the objects*/
    //gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm0), TRUE);
    //gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm1), TRUE);
    //gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm4b), TRUE);
    //gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm4), TRUE);
    //gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm10), TRUE);
    //gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm99), TRUE);

    for (i = 0; i < xgc->data.set_mat.nmats; i++) {
        if (xgc->data.set_mat.overpos[i] == posid) {
            overriden = g_strdup(xgc->data.set_mat.materials[i]);
            mat.type = 99;
        }
    }

    if (mat.type == 0) { /* 1) voxel-by voxel */
        mattype = 0;

        //gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm0), TRUE);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.epsilon_spin), mat.epsilon);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.mu_spin), mat.mu);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.sigma_spin), mat.sigma);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.sigast_spin), mat.sigast);

        if (overriden) { //remove overriden entry if there was any (point it nowhere)
            for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                if (xgc->data.set_mat.overpos[i] == posid)
                    xgc->data.set_mat.overpos[i] = -1;
                    //xgc->data.set_mat.valid[i] = FALSE;
            }
        }
    } else if (mat.type == 1) { /* 2) tabulated */
        mattype = 1;
        //gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm1), TRUE);

        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tepsilon_spin), mat.epsilon);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tmu_spin), mat.mu);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tsigma_spin), mat.sigma);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tsigast_spin), mat.sigast);

        if (overriden) { //remove overriden entry if there was any (point it nowhere)
            for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                if (xgc->data.set_mat.overpos[i] == posid)
                    xgc->data.set_mat.overpos[i] = -1;
            }
        }

        //pmat->type = 1;
        //pmat->epsilon = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tepsilon_spin));
        //pmat->mu = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tmu_spin));
        //pmat->sigma = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tsigma_spin));
        //pmat->sigast = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.tsigast_spin));
    } else if (mat.type == 10) { /* 3) perfect electric conductor */
        mattype = 10;

        /*pmat->type = 10;

        if (overriden) { //remove overriden entry if there was any (point it nowhere)
        for (i = 0; i < xgc->data.set_mat.nmats; i++) {
        if (xgc->data.set_mat.overpos[i] == posid)
        xgc->data.set_mat.overpos[i] = -1;
        }
        }*/
    } else if (mat.type == 4 || mat.type == 5 || mat.type == 6) {  /* 4) CP model based */
        //posid... id within 0..n objects
        //nmats... 0 to n of overriden materials
        //overpos[i] shows which object posid string [i] points to
        //materials[i] is the string

        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cepsilon_spin), mat.epsilon);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.omegap_spin), mat.drude_omega_p);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.nu_spin), mat.drude_nu);

        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.a_spin[0]), mat.cp3_a[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.ia_phi_spin[0]), mat.cp3_phi[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.omega_spin[0]), mat.cp3_omega[0]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.gamma_spin[0]), mat.cp3_gamma[0]);

        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.a_spin[1]), mat.cp3_a[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.ia_phi_spin[1]), mat.cp3_phi[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.omega_spin[1]), mat.cp3_omega[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.gamma_spin[1]), mat.cp3_gamma[1]);
               
        metalmodelval = mat.type;

        if (metalmodelval >= 4 && metalmodelval <= 6)
            gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs_mat.oc.metal_model), metalmodelval - 4);
        else
            gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs_mat.oc.metal_model), 1);        
    } else if (mat.type == 99) { /* Loaded from database */
        mattype = 99;
        //gtk_entry_set_text(GTK_ENTRY(xgc->cs_mat.oc.overriden), overriden);
        
        i = gtk_combo_box_get_text_index(GTK_COMBO_BOX(xgc->cs_mat.oc.overriden), overriden);
        if (i != -1)
            gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs_mat.oc.overriden), i);
        else
            gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->cs_mat.oc.overriden), MAT_PROP_MATERIAL_NONE_INDEX);

        /*getmat = TRUE;
        pmat->epsilon = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.cepsilon_spin));
        pmat->drude_omega_p = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.omegap_spin));
        pmat->drude_nu = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.nu_spin));

        pmat->cp3_a[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.a_spin[0]));
        pmat->cp3_phi[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.ia_phi_spin[0]));
        pmat->cp3_omega[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.omega_spin[0]));
        pmat->cp3_gamma[0] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.gamma_spin[0]));

        pmat->cp3_a[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.a_spin[1]));
        pmat->cp3_phi[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.ia_phi_spin[1]));
        pmat->cp3_omega[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.omega_spin[1]));
        pmat->cp3_gamma[1] = (GLfloat)gtk_spin_button_get_value(GTK_SPIN_BUTTON(xgc->cs_mat.oc.gamma_spin[1]));


        if (overriden) { //remove overriden entry if there was any (point it nowhere)
            for (i = 0; i < xgc->data.set_mat.nmats; i++) {
                if (xgc->data.set_mat.overpos[i] == posid)
                    xgc->data.set_mat.overpos[i] = -1;
            }
        }*/
    } 

    /*if (getmat) {
        if (strcmp(gtk_combo_box_get_active_text(GTK_COMBO_BOX(xgc->cs_mat.oc.combo_metal)), "PLRC") == 0)
            pmat->type = 6;
        else if (strcmp(gtk_combo_box_get_active_text(GTK_COMBO_BOX(xgc->cs_mat.oc.combo_metal)), "ADE") == 0)
            pmat->type = 5;
        else
            pmat->type = 4;
    }*/

    if (mattype == 0)       /* 1) voxel-by voxel */
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm0), TRUE);
    else if (mattype == 1)  /* 2) tabulated */
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm1), TRUE);    
    else if (mattype == 10) /* 3) perfect electric conductor */
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm10), TRUE);
    else if (mattype == 4 || mattype == 5 || mattype == 6)  /* 4) CP model based */
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm4), TRUE);
    else if (mattype == 99) /* 5) loaded from database */
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(xgc->cs_mat.oc.mm99), TRUE);
    //////////////////////////////////////////////////////////////////////////
    
    xgc->mat_lock_controls_update = FALSE;

    return FALSE;
} /* set_values_object */

void
mat_controls_row_select_process(XGControls *xgc, gint id, GtkTreePath *path)
{
    gboolean writeme = FALSE;
    gint index;

    if (id >= SET_MAT_SPHERE && id < SET_MAT_VOXEL) {
        /* I. 1. Sphere properties */
        writeme = set_values_object(xgc, &(xgc->data.set_mat), SV_TYPE_MAT_SPHERE, id);
    } else if (id >= SET_MAT_VOXEL && id < SET_MAT_CYLINDER) {
        /* I. 2. Box properties */
        writeme = set_values_object(xgc, &(xgc->data.set_mat), SV_TYPE_MAT_VOXEL, id);
    } else if (id >= SET_MAT_CYLINDER && id < SET_MAT_CONE) {
        /* I. 3. Cylinder properties */
        writeme = set_values_object(xgc, &(xgc->data.set_mat), SV_TYPE_MAT_CYLINDER, id);
    } else if (id >= SET_MAT_CONE && id < SET_MAT_RCONE) {
        /* I. 4. Cone properties */
        writeme = set_values_object(xgc, &(xgc->data.set_mat), SV_TYPE_MAT_CONE, id);
    } else if (id >= SET_MAT_RCONE && id < SET_MAT_GWYDD) {
        /* I. 5. Cut cone properties */
        writeme = set_values_object(xgc, &(xgc->data.set_mat), SV_TYPE_MAT_RCONE, id);
    } else if (id >= SET_MAT_GWYDD && id < SET_MAT_MESH) {
    } else if (id >= SET_MAT_MESH && id < SET_MAT_TETRAHEDRON) {
    } else if (id >= SET_MAT_TETRAHEDRON) {
        /* I. 8. Tetrahedron properties */
        writeme = set_values_object(xgc, &(xgc->data.set_mat), SV_TYPE_MAT_TETRAHEDRON, id);
    }


    /*rewrite the temporary matfile and then reload it and interpret completely in order to assure consistency between its text and gui interpretation*/
    if (writeme) {
        if (write_matfile(xgc)) {
            index = gtk_tree_path_get_indices(path)[0];
            path = gtk_tree_path_new_from_indices(0, -1);
            gtk_tree_view_expand_row(GTK_TREE_VIEW(xgc->tv_mat), path, TRUE);
        }
    }
}
