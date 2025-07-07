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


/* main_menu.c :
*  main menu functions
*/

#include "global.h"
#include "xs_main_menu.h"
#include "xsvit.h"
#include "xs_scene.h"
#include "xs_par_tree.h"
#include "write.h"

#include <gtk/gtk.h>
#include <glib/gstdio.h>

/*static void
activate_action (GtkAction *action)
{
    printf("activate_action\n");
}

static void
about_cb (GtkAction *action, GtkWidget *window)
{
printf("about_cb\n");
}*/

typedef struct
{
    GtkAction action;
} ToolMenuAction;

typedef struct
{
    GtkActionClass parent_class;
} ToolMenuActionClass;

G_DEFINE_TYPE(ToolMenuAction, tool_menu_action, GTK_TYPE_ACTION)

static void
tool_menu_action_class_init (ToolMenuActionClass *class)
{
    GTK_ACTION_CLASS (class)->toolbar_item_type = GTK_TYPE_MENU_TOOL_BUTTON;
}

static void
tool_menu_action_init (ToolMenuAction *action)
{
}

static GtkActionEntry entries[] = {
    {"FileMenu", NULL, "_File"},                /* name, stock id, label */
    {"EditParametersMenu", NULL, "Edit _parameters"}, /* name, stock id, label */
    {"EditMaterialsMenu", NULL, "Edit _material objects"}, /* name, stock id, label */
    {"ExecuteMenu", NULL, "E_xecute"},          /* name, stock id, label */
    {"ViewMenu", NULL, "_View"},                /* name, stock id, label */
    {"HelpMenu", NULL, "_Help"},                /* name, stock id, label */

    {"New", GTK_STOCK_NEW,                      /* name, stock id */
    "_New", "<control>N",                       /* label, accelerator */
    "Create new parameter file (*.par)",        /* tooltip */
    G_CALLBACK (file_new_cb)},
    {"Open", GTK_STOCK_OPEN,                    /* name, stock id */
    "_Open", "<control>O",                      /* label, accelerator */
    "Open parameter file (*.par)..",            /* tooltip */
    G_CALLBACK (file_open_cb)},
    {"Save", GTK_STOCK_SAVE,                    /* name, stock id */
    "_Save","<control>S",                       /* label, accelerator */
    "Save current file",                        /* tooltip */
    G_CALLBACK (file_save_cb)},
    {"SaveAs", GTK_STOCK_SAVE_AS,               /* name, stock id */
    "Save _As..", "",                           /* label, accelerator */
    "Save parameter and material file as (*.par, *.mat)..",                          /* tooltip */
    G_CALLBACK (file_save_as_cb)},
    {"SaveMatAs", GTK_STOCK_SAVE_AS,            /* name, stock id */
    "Save Material File As..", NULL,            /* label, accelerator */
    "Save material file under different name than stated in parameter file..",       /* tooltip */
    G_CALLBACK (file_save_mat_as_cb)},
    {"RecentFiles", NULL,                       /* name, stock id */
    "_Recent files", "",                        /* label, accelerator */
    "Open recently used file",                  /* tooltip */
    NULL/*G_CALLBACK (recent_chooser_item_activated_cb)*/},
    {"Preferences", NULL,                       /* name, stock id */
    "_Preferences", "",                         /* label, accelerator */
    "Set user preferences for XSvit",           /* tooltip */
    G_CALLBACK (preferences_cb)},
    {"Quit", GTK_STOCK_QUIT,                    /* name, stock id */
    "_Quit", "<control>Q",                      /* label, accelerator */
    "Quit program",                             /* tooltip */
    G_CALLBACK (quit_cb)},

    {"AddPointSource", NULL,                    /* name, stock id */
    "Add point source", "",                     /* label, accelerator */
    "Add point source",                         /* tooltip */
    G_CALLBACK (par_add_src_point_cb)},
    {"AddSFSource", NULL,                       /* name, stock id */
    "Add SF source", "",                        /* label, accelerator */
    "Add scattered field source",               /* tooltip */
    G_CALLBACK (par_add_src_sf_cb)},
    {"AddTSFSource", NULL,                      /* name, stock id */
    "Add TSF source", "",                       /* label, accelerator */
    "Add total/scattered field source",         /* tooltip */
    G_CALLBACK (par_add_src_tsf_cb)},
    {"AddTSFFSource", NULL,                     /* name, stock id */
    "Add TSFF source..", "",                    /* label, accelerator */
    "Add focused total/scattered field source", /* tooltip */
    G_CALLBACK (par_add_src_tsff_cb)},
    {"AddLTSFSource", NULL,                     /* name, stock id */
    "Add LTSF source", "",                      /* label, accelerator */
    "Add layered total/scattered field source", /* tooltip */
    G_CALLBACK (par_add_src_ltsf_cb)},
    {"AddLTSFFSource", NULL,                    /* name, stock id */
    "Add LTSFF source", "",                     /* label, accelerator */
    "Add layered focused total/scattered field source", /* tooltip */
    G_CALLBACK (par_add_src_ltsff_cb)},
    {"AddPointOutput", NULL,                    /* name, stock id */
    "Add point output", "",                     /* label, accelerator */
    "Add point output",                         /* tooltip */
    G_CALLBACK (par_add_out_point_cb)},
    {"AddImageOutput", NULL,                    /* name, stock id */
    "Add image output..", "",                   /* label, accelerator */
    "Add image output",                         /* tooltip */
    G_CALLBACK (par_add_out_image_cb)},
    {"AddPlaneOutput", NULL,                    /* name, stock id */
    "Add plane output", "",                     /* label, accelerator */
    "Add plane output",                         /* tooltip */
    G_CALLBACK (par_add_out_plane_cb)},
    {"AddVolumeOutput", NULL,                   /* name, stock id */
    "Add volume output", "",                    /* label, accelerator */
    "Add volume output",                        /* tooltip */
    G_CALLBACK (par_add_out_volume_cb)},
    {"AddSumOutput", NULL,                      /* name, stock id */
    "Add sum output", "",                       /* label, accelerator */
    "Add sum output",                           /* tooltip */
    G_CALLBACK (par_add_out_sum_cb)},
    {"AddForceOutput", NULL,                    /* name, stock id */
    "Add force output", "",                     /* label, accelerator */
    "Add force output",                         /* tooltip */
    G_CALLBACK (par_add_out_force_cb)},
    {"AddGrowthModifier", NULL,                 /* name, stock id */
    "Add growth modifier", "",                  /* label, accelerator */
    "Add roughness via ballistic deposition",   /* tooltip */
    G_CALLBACK (par_add_grow_cb)},
    {"AddRoughnessModifier", NULL,              /* name, stock id */
    "Add roughness modifier", "",               /* label, accelerator */
    "Add roughness via random Gaussians",       /* tooltip */
    G_CALLBACK (par_add_roughness_cb)},
    {"AddSpectralModifier", NULL,               /* name, stock id */
    "Add spectral modifier", "",                /* label, accelerator */
    "Add roughness via spectral synthesis",     /* tooltip */
    G_CALLBACK (par_add_spectral_cb)},
    {"AddExpressionModifier", NULL,             /* name, stock id */
    "Add expression modifier", "",              /* label, accelerator */
    "Add roughness via analytical expression",  /* tooltip */
    G_CALLBACK (par_add_expression_cb)},
    {"AddNFFFPoint", NULL,                      /* name, stock id */
    "Add NFFF point", "",                       /* label, accelerator */
    "Add near field to far field transform point",      /* tooltip */
    G_CALLBACK (par_add_nfffp_cb)},
    {"AddNFFFArea", NULL,                       /* name, stock id */
    "Add NFFF area", "",                        /* label, accelerator */
    "Add near field to far field area",         /* tooltip */
    G_CALLBACK (par_add_nfffa_cb)},
    {"AddPeriodicNFFFPoint", NULL,              /* name, stock id */
    "Add periodic NFFF point", "",              /* label, accelerator */
    "Add periodic near field to far field transform point", /* tooltip */
    G_CALLBACK (par_add_pnfffp_cb)},
    {"AddPeriodicNFFFArea", NULL,               /* name, stock id */
    "Add periodic NFFF area", "",               /* label, accelerator */
    "Add periodic near field to far field area",        /* tooltip */
    G_CALLBACK (par_add_pnfffa_cb)},
    {"RemoveParameter", NULL,                   /* name, stock id */
    "Remove parameter", "",                     /* label, accelerator */
    "Remove parameter",                         /* tooltip */
    G_CALLBACK (par_row_remove_cb)},

    {"AddSphere", NULL,                         /* name, stock id */
    "Add _sphere", "",                          /* label, accelerator */
    "Add sphere",                               /* tooltip */
    G_CALLBACK (mat_add_sphere_cb)},
    {"AddBox", NULL,                            /* name, stock id */
    "Add _box", "",                             /* label, accelerator */
    "Add box (parallelpiped)",                  /* tooltip */
    G_CALLBACK (mat_add_box_cb)},
    {"AddCylinder", NULL,                       /* name, stock id */
    "Add _cylinder", "",                        /* label, accelerator */
    "Add cylinder",                             /* tooltip */
    G_CALLBACK (mat_add_cylinder_cb)},
    {"AddCone", NULL,                           /* name, stock id */
    "Add c_one", "",                            /* label, accelerator */
    "Add cone",                                 /* tooltip */
    G_CALLBACK (mat_add_cone_cb)},
    {"AddCutCone", NULL,                        /* name, stock id */
    "Add c_ut cone", "",                        /* label, accelerator */
    "Add cut cone",                             /* tooltip */
    G_CALLBACK (mat_add_rcone_cb)},
    {"AddTetrahedron", NULL,                    /* name, stock id */
    "Add _tetrahedron", "",                     /* label, accelerator */
    "Add single tetrahedron",                   /* tooltip */
    G_CALLBACK (mat_add_tetrahedron_cb)},
    {"RemoveMatObject", NULL,                   /* name, stock id */
    "Remove material object", "",               /* label, accelerator */
    "Remove material object",                   /* tooltip */
    G_CALLBACK (mat_row_remove_cb)},

    {"Run", GTK_STOCK_MEDIA_PLAY,               /* name, stock id */
    "_Run", "F5",                               /* label, accelerator */
    "Start GSvit solver with the current parameter file",   /* tooltip */
    G_CALLBACK (gsvit_run_cb)},
    {"Stop", GTK_STOCK_MEDIA_STOP   ,           /* name, stock id */
    "_Stop", "<shift>F5",                       /* label, accelerator */
    "Kill GSvit solver",                        /* tooltip */
    G_CALLBACK (gsvit_stop_cb)},
    {"ViewResults", GTK_STOCK_FIND,             /* name, stock id */
    "_View results", "",                        /* label, accelerator */
    "Start GSvit solver with the current parameter file",   /* tooltip */
    G_CALLBACK (gwydd_cb)},

    {"ZoomToFit", GTK_STOCK_ZOOM_FIT,           /* name, stock id */
    "Zoom To _Fit", "",                         /* label, accelerator */
    "Zoom objects in 3D scene to fit the view", /* tooltip */
    G_CALLBACK (scene_reset_cb)},
    {"ZoomIn", GTK_STOCK_ZOOM_IN,               /* name, stock id */
    "Zoom _In", "",                             /* label, accelerator */
    "Zoom in objects in 3D scene",              /* tooltip */
    G_CALLBACK (scene_zoom_in_cb)},
    {"ZoomOut", GTK_STOCK_ZOOM_OUT,             /* name, stock id */
    "Zoom _Out", "",                            /* label, accelerator */
    "Zoom out objects in 3D scene",             /* tooltip */
    G_CALLBACK (scene_zoom_out_cb)},
    {"MoveCameraLeft", GTK_STOCK_GO_BACK,       /* name, stock id */
    "Move view _left", "",                    /* label, accelerator */
    "Move view left",                         /* tooltip */
    G_CALLBACK (scene_move_view_right_cb)},
    {"MoveCameraRight", GTK_STOCK_GO_FORWARD,   /* name, stock id */
    "Move view _right", "",                   /* label, accelerator */
    "Move view right",                        /* tooltip */
    G_CALLBACK (scene_move_view_left_cb)},
    {"MoveCameraUp", GTK_STOCK_GO_UP,           /* name, stock id */
    "Move view _up", "",                      /* label, accelerator */
    "Move view up",                           /* tooltip */
    G_CALLBACK (scene_move_view_up_cb)},
    {"MoveCameraDown", GTK_STOCK_GO_DOWN,       /* name, stock id */
    "Move view _down", "",                    /* label, accelerator */
    "Move view down",                         /* tooltip */
    G_CALLBACK (scene_move_view_down_cb)},

    {"CheckFilesConsitency", GTK_STOCK_PROPERTIES,       /* name, stock id */
    "Check consitency", "",                    /* label, accelerator */
    "Check parameter file and material file consitency",                         /* tooltip */
    G_CALLBACK (check_files_consistency_cb)},

    {"XSvitHelp", NULL,                         /* name, stock id */
    "_XSvit online help", "",                   /* label, accelerator */
    "Open online help for XSvit interface in web browser",                                    /* tooltip */
    G_CALLBACK (help_x_cb)},
    {"GSvitHelp", NULL,                         /* name, stock id */
    "_GSvit online help", "",                   /* label, accelerator */
    "Open online help for GSvit interface in web browser",                                    /* tooltip */
    G_CALLBACK (help_g_cb)},
    {"About", NULL,                             /* name, stock id */
    "_About", "",                               /* label, accelerator */
    "Authors and license",                      /* tooltip */
    G_CALLBACK (about_cb)},
};
static guint n_entries = G_N_ELEMENTS (entries);


static const gchar *ui_info =
"<ui>"
"  <menubar name='MenuBar'>"
"    <menu action='FileMenu'>"
"      <menuitem action='New'/>"
"      <menuitem action='Open'/>"
"      <menuitem action='Save'/>"
"      <menuitem action='SaveAs'/>"
"      <menuitem action='SaveMatAs'/>"
"      <separator/>"
"      <menuitem action='RecentFiles'/>"
"      <separator/>"
"      <menuitem action='Preferences'/>"
"      <separator/>"
"      <menuitem action='Quit'/>"
"    </menu>"
"    <menu action='EditParametersMenu'>"
"      <menuitem action='AddPointSource'/>"
"      <menuitem action='AddSFSource'/>"
"      <menuitem action='AddTSFSource'/>"
"      <menuitem action='AddTSFFSource'/>"
"      <menuitem action='AddLTSFSource'/>"
"      <menuitem action='AddLTSFFSource'/>"
"      <separator/>"
"      <menuitem action='AddPointOutput'/>"
"      <menuitem action='AddImageOutput'/>"
"      <menuitem action='AddPlaneOutput'/>"
"      <menuitem action='AddVolumeOutput'/>"
"      <menuitem action='AddSumOutput'/>"
"      <menuitem action='AddForceOutput'/>"
"      <separator/>"
"      <menuitem action='AddGrowthModifier'/>"
"      <menuitem action='AddRoughnessModifier'/>"
"      <menuitem action='AddSpectralModifier'/>"
"      <menuitem action='AddExpressionModifier'/>"
"      <separator/>"
"      <menuitem action='AddNFFFPoint'/>"
"      <menuitem action='AddNFFFArea'/>"
"      <menuitem action='AddPeriodicNFFFPoint'/>"
"      <menuitem action='AddPeriodicNFFFArea'/>"
"      <separator/>"
"      <menuitem action='RemoveParameter'/>"
"    </menu>"
"    <menu action='EditMaterialsMenu'>"
"      <menuitem action='AddSphere'/>"
"      <menuitem action='AddBox'/>"
"      <menuitem action='AddCylinder'/>"
"      <menuitem action='AddCone'/>"
"      <menuitem action='AddCutCone'/>"
"      <menuitem action='AddTetrahedron'/>"
"      <separator/>"
"      <menuitem action='RemoveMatObject'/>"
"    </menu>"
"    <menu action='ExecuteMenu'>"
"      <menuitem action='Run'/>"
"      <menuitem action='Stop'/>"
"      <menuitem action='ViewResults'/>"
"    </menu>"
"    <menu action='HelpMenu'>"
"      <menuitem action='XSvitHelp'/>"
"      <menuitem action='GSvitHelp'/>"
"      <menuitem action='About'/>"
"    </menu>"
"  </menubar>"
"  <toolbar name='ToolBar'>"
"    <toolitem action='New'/>"
"    <toolitem action='Open'/>"
"    <toolitem action='Save'/>"
"    <toolitem action='SaveAs'/>"
"    <separator action='Sep1'/>"
"    <toolitem action='ZoomToFit'/>"
"    <toolitem action='ZoomIn'/>"
"    <toolitem action='ZoomOut'/>"
"    <toolitem action='MoveCameraLeft'/>"
"    <toolitem action='MoveCameraRight'/>"
"    <toolitem action='MoveCameraUp'/>"
"    <toolitem action='MoveCameraDown'/>"
"    <separator action='Sep2'/>"
"    <toolitem action='CheckFilesConsitency'/>"
"    <toolitem action='Run'/>"
"    <toolitem action='Stop'/>"
"    <toolitem action='ViewResults'/>"
"    <separator action='Sep3'/>"
"    <toolitem action='Quit'/>"
"  </toolbar>"
"</ui>";


void create_main_menubar_and_toolbar(XGControls *xgc, GtkWidget** menubar, GtkWidget** toolbar)
{
    GtkUIManager *merge;
    GtkActionGroup *action_group;
    //GtkAction *open_action;
    //GtkWidget *bar;
    GtkWidget *window = xgc->toplevel;
    GError *error = NULL;

    action_group = gtk_action_group_new ("AppWindowActions");
    /*open_action = g_object_new (tool_menu_action_get_type(),
                                "name", "Open",
                                "label", "_Open",
                                "tooltip", "Open a file",
                                "stock-id", GTK_STOCK_OPEN,
                                NULL);
    gtk_action_group_add_action (action_group, open_action);
    g_object_unref (open_action);*/

    //gtk_action_group_add_actions (action_group, entries, n_entries, window);
    gtk_action_group_add_actions (action_group, entries, n_entries, xgc);

//    gtk_action_group_add_toggle_actions (action_group, toggle_entries, n_toggle_entries, NULL);
//    gtk_action_group_add_radio_actions (action_group, color_entries, n_color_entries, COLOR_RED, G_CALLBACK (activate_radio_action), NULL);
//    gtk_action_group_add_radio_actions (action_group, shape_entries, n_shape_entries, SHAPE_SQUARE, G_CALLBACK (activate_radio_action), NULL);

    merge = gtk_ui_manager_new();
    g_object_set_data_full(G_OBJECT (window), "ui-manager", merge, g_object_unref);
    gtk_ui_manager_insert_action_group(merge, action_group, 0);
    gtk_window_add_accel_group(GTK_WINDOW (window),
        gtk_ui_manager_get_accel_group(merge));

    if (!gtk_ui_manager_add_ui_from_string(merge, ui_info, -1, &error)) {
        g_message ("building menus failed: %s", error->message);
        g_error_free (error);
    }    

    *menubar = gtk_ui_manager_get_widget(merge, "/MenuBar");
    gtk_widget_show(*menubar);

    //////////////////////////////////////////////////////////////////////////
    //GtkWidget *recent_files;
    //recent_files = gtk_ui_manager_get_widget(merge, "/MenuBar/FileMenu/RecentFiles");

    xgc->recentfiles_menu_item = gtk_ui_manager_get_widget(merge, "/MenuBar/FileMenu/RecentFiles");

    //GList *list;
    //list = gtk_container_get_children(*menubar);
    //fm = g_list_nth(list, 0);

    //xgc->recentfiles_menu_item = gtk_menu_item_new_with_mnemonic("_Recent files");
    //gtk_menu_shell_append(GTK_MENU_SHELL(*menubar), xgc->recentfiles_menu_item);
    //gtk_menu_shell_append(GTK_MENU_SHELL(node->data), xgc->recentfiles_menu_item);

    //gtk_menu_append(fm, xgc->recentfiles_menu_item);

    //gtk_menu_shell_insert(GTK_MENU_SHELL(*menubar), xgc->recentfiles_menu_item, 1);

    //gtk_menu_shell_insert(GTK_MENU_SHELL(fm), xgc->recentfiles_menu_item, 1);
    //gtk_menu_item_set_submenu(GTK_MENU_ITEM(fm), xgc->recentfiles_menu_item);

    //GtkWidget *w1 =	gtk_menu_item_get_submenu(GTK_MENU_ITEM(*menubar));
    //gtk_menu_shell_append(GTK_MENU_SHELL(w1), xgc->recentfiles_menu_item);



    GtkWidget *recentchooser_menuitem = gtk_recent_chooser_menu_new_for_manager(gtk_recent_manager_get_default());

    gtk_recent_chooser_set_sort_type(GTK_RECENT_CHOOSER(recentchooser_menuitem), GTK_RECENT_SORT_MRU);

    GtkRecentFilter * recent_filter = gtk_recent_filter_new();
    gtk_recent_filter_add_pattern(recent_filter, "*.par");
    gtk_recent_chooser_set_filter(GTK_RECENT_CHOOSER(recentchooser_menuitem), recent_filter);
    gtk_recent_chooser_menu_set_show_numbers(GTK_RECENT_CHOOSER_MENU(recentchooser_menuitem), TRUE);
    gtk_recent_chooser_set_show_not_found(GTK_RECENT_CHOOSER(recentchooser_menuitem), FALSE);
    gtk_menu_item_set_submenu(GTK_MENU_ITEM(xgc->recentfiles_menu_item), recentchooser_menuitem);
    g_signal_connect(G_OBJECT(recentchooser_menuitem), "item-activated", G_CALLBACK(recent_chooser_item_activated_cb), xgc);
    //////////////////////////////////////////////////////////////////////////    

    xgc->new_menu_item = gtk_ui_manager_get_widget(merge, "/MenuBar/FileMenu/New");
    xgc->open_menu_item = gtk_ui_manager_get_widget(merge, "/MenuBar/FileMenu/Open");
    xgc->save_menu_item = gtk_ui_manager_get_widget(merge, "/MenuBar/FileMenu/Save");
    xgc->saveas_menu_item = gtk_ui_manager_get_widget(merge, "/MenuBar/FileMenu/SaveAs");
    xgc->savematas_menu_item = gtk_ui_manager_get_widget(merge, "/MenuBar/FileMenu/SaveMatAs");

    xgc->run_menu_item = gtk_ui_manager_get_widget(merge, "/MenuBar/ExecuteMenu/Run");
    xgc->stop_menu_item = gtk_ui_manager_get_widget(merge, "/MenuBar/ExecuteMenu/Stop");
    xgc->gwydd_menu_item = gtk_ui_manager_get_widget(merge, "/MenuBar/ExecuteMenu/ViewResults");

    *toolbar = gtk_ui_manager_get_widget (merge, "/ToolBar");
    gtk_toolbar_set_style(GTK_TOOLBAR(*toolbar), GTK_TOOLBAR_ICONS);
    gtk_widget_show (*toolbar);
    

    xgc->run_tool_item = GTK_TOOL_ITEM(gtk_ui_manager_get_widget(merge, "/ToolBar/Run"));
    xgc->stop_tool_item = GTK_TOOL_ITEM(gtk_ui_manager_get_widget(merge, "/ToolBar/Stop"));
    xgc->gwydd_tool_item = GTK_TOOL_ITEM(gtk_ui_manager_get_widget(merge, "/ToolBar/ViewResults"));
}

