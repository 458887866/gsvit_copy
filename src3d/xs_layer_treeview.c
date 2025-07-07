/*
*  Copyright (C) 2017 Petr Klapetek
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

/*  layer_treeview.c :
*  Tree view control for source layers
*/

#include <gtk/gtk.h>
#include <glib-object.h>
#include <gdk/gdkkeysyms.h>
#include <string.h>
#include <stdlib.h>
#include "xs_layer_treeview.h"
#include "xs_dialogs.h"


typedef struct
{
    gint		number;
    gint		zpos;
    gdouble     epsilon;
    gdouble     mu;
    gdouble     sigma;
    gdouble     sigast;
    gchar*      material;
} Item;

enum
{
    COLUMN_LAYER_NUMBER,
    COLUMN_LAYER_Z,
    COLUMN_LAYER_EPSILON,
    COLUMN_LAYER_MU,
    COLUMN_LAYER_SIGMA,
    COLUMN_LAYER_SIGAST,
    COLUMN_LAYER_MATERIAL,
    NUM_LAYER_COLUMNS
};

//static gint row_to_edit_after_tab = -1;
//static gint column_to_edit_after_tab = -1;

static GArray *articles = NULL;


#define GSVIT_TRREEVIEW_GET_PRIVATE(obj)  \
   (G_TYPE_INSTANCE_GET_PRIVATE((obj), GSVIT_TYPE_TREEVIEW, TreeViewWrapperPrivate))


static void tree_view_wrapper_finalize(GObject *object);

G_DEFINE_TYPE(TreeViewWrapper, tree_view_wrapper, G_TYPE_OBJECT)


TreeViewWrapper *treeviewwrapper;


static void
tree_view_wrapper_class_init(TreeViewWrapperClass *klass)
{
    GObjectClass *gobject_class = G_OBJECT_CLASS(klass);

    gobject_class->finalize = tree_view_wrapper_finalize;
}

static void
tree_view_wrapper_init(G_GNUC_UNUSED TreeViewWrapper *object)
{
}

static void
tree_view_wrapper_finalize(GObject *object)
{
    g_object_unref(((TreeViewWrapper *)object)->items_model);
}



/*
static void
add_items (void)
{
    Item item;

    g_return_if_fail (articles != NULL);

    item.number = 1;
    item.zpos = 0;
    item.epsilon = 1.0;
    item.mu = 1.0;
    item.sigma = 0.0;
    item.sigast = 0.0;
    g_array_append_vals (articles, &item, 1);
}*/

void
set_ltv_items(TreeViewWrapper *treeview_wrapper, gpointer* data)
{
    Item item;
    GtkTreeIter iter;
    guint i;
    LTV_Data *ltvd = (LTV_Data*)data;
    //GtkTreeModel *model = gtk_tree_view_get_model(treeview);
    GtkTreeModel *model = treeview_wrapper->items_model;

    if (articles != NULL)
        g_array_unref(articles);

    /* create array */
    articles = g_array_sized_new(FALSE, FALSE, sizeof(Item), ltvd->nlayers);

    for (i = 0; i < ltvd->nlayers; i++) {
        item.number = i + 1;
        item.zpos = ltvd->zpos[i];
        item.epsilon = ltvd->epsilon[i];
        item.mu = ltvd->mu[i];
        item.sigma = ltvd->sigma[i];
        item.sigast = ltvd->sigast[i];
        if (NULL != ltvd->material[i])
            item.material = g_strdup(ltvd->material[i]);
        g_array_append_vals(articles, &item, 1);
    }

    /* clear list store */
    gtk_list_store_clear(GTK_LIST_STORE(model));

    /* add items */
    for (i = 0; i < articles->len; i++) {
        gtk_list_store_append(GTK_LIST_STORE(model), &iter);

        gtk_list_store_set(GTK_LIST_STORE(model), &iter,
                           COLUMN_LAYER_NUMBER, g_array_index(articles, Item, i).number,
                           COLUMN_LAYER_Z, g_array_index(articles, Item, i).zpos,
                           COLUMN_LAYER_EPSILON, g_array_index(articles, Item, i).epsilon,
                           COLUMN_LAYER_MU, g_array_index(articles, Item, i).mu,
                           COLUMN_LAYER_SIGMA, g_array_index(articles, Item, i).sigma,
                           COLUMN_LAYER_SIGAST, g_array_index(articles, Item, i).sigast,
                           COLUMN_LAYER_MATERIAL, g_array_index(articles, Item, i).material,
                           -1);
    }
}


/* enumerate all items and set index */

static void
index_items(GtkTreeModel* model)
{
    int i = 0;
    GtkTreeIter iter;

    if (gtk_tree_model_get_iter_first(model, &iter)) {
        do {
            gtk_list_store_set(GTK_LIST_STORE(model), &iter, COLUMN_LAYER_NUMBER, i, -1);
            i++;
        } while (gtk_tree_model_iter_next(model, &iter));
    }
}

static void
insert_item_cb(GtkWidget *button, gpointer data)
{
    //GtkTreeModel *model = (GtkTreeModel *)data;
    //TreeViewData *treeviewdata = (TreeViewData *)data;
    TreeViewWrapper *tvw = (TreeViewWrapper *)data;
    GtkTreeModel *model = tvw->items_model;
    GtkTreeView * treeview = GTK_TREE_VIEW(tvw->parent_treeview);

    Item item;
    GtkTreeIter iter, iter_sibling;
    GtkTreeSelection *selection = gtk_tree_view_get_selection(treeview);

    g_return_if_fail(articles != NULL);

    item.number = gtk_tree_model_iter_n_children(GTK_TREE_MODEL(model), NULL) + 1;
    item.zpos = SOURCE_LAYERED_ZPOS;
    item.epsilon = SOURCE_LAYERED_EPS;
    item.mu = SOURCE_LAYERED_MU;
    item.sigma = SOURCE_LAYERED_SIGMA;
    item.sigast = SOURCE_LAYERED_SIGAST;
    item.material = g_strdup(SOURCE_LAYERED_MATERIAL);
   
    g_array_append_vals(articles, &item, 1);    

    if (gtk_tree_selection_get_selected(selection, NULL, &iter_sibling)) {
        gtk_list_store_insert_after(GTK_LIST_STORE(model), &iter, &iter_sibling);
    } else
        gtk_list_store_append(GTK_LIST_STORE(model), &iter);

    gtk_list_store_set(GTK_LIST_STORE(model), &iter,
                       COLUMN_LAYER_NUMBER, item.number,
                       COLUMN_LAYER_Z, item.zpos,
                       COLUMN_LAYER_EPSILON, item.epsilon,
                       COLUMN_LAYER_MU, item.mu,
                       COLUMN_LAYER_SIGMA, item.sigma,
                       COLUMN_LAYER_SIGAST, item.sigast,
                       COLUMN_LAYER_MATERIAL, item.material,
                       -1);

    index_items(model);

    g_free(item.material);

    g_signal_emit_by_name(tvw->parent_treeview, "cell-edited");
}

static GtkTreeModel *
create_items_model(void)
{
    guint i = 0;
    GtkListStore *model;
    GtkTreeIter iter;
    Item item;

    /* create array */
    articles = g_array_sized_new(FALSE, FALSE, sizeof(Item), 1);

    /* create list store */
    model = gtk_list_store_new(NUM_LAYER_COLUMNS, G_TYPE_INT, G_TYPE_INT, G_TYPE_DOUBLE, G_TYPE_DOUBLE, G_TYPE_DOUBLE, G_TYPE_DOUBLE, G_TYPE_STRING);

    item.number = gtk_tree_model_iter_n_children(GTK_TREE_MODEL(model), NULL) + 1;
    item.zpos = SOURCE_LAYERED_ZPOS;
    item.epsilon = SOURCE_LAYERED_EPS;
    item.mu = SOURCE_LAYERED_MU;
    item.sigma = SOURCE_LAYERED_SIGMA;
    item.sigast = SOURCE_LAYERED_SIGAST;
    item.material = g_strdup(SOURCE_LAYERED_MATERIAL);
    
    g_array_append_vals(articles, &item, 1);

    /* add items */
    for (i = 0; i < articles->len; i++) {
        gtk_list_store_append(model, &iter);

        gtk_list_store_set(model, &iter,
                           COLUMN_LAYER_NUMBER, g_array_index(articles, Item, i).number,
                           COLUMN_LAYER_Z, g_array_index(articles, Item, i).zpos,
                           COLUMN_LAYER_EPSILON, g_array_index(articles, Item, i).epsilon,
                           COLUMN_LAYER_MU, g_array_index(articles, Item, i).mu,
                           COLUMN_LAYER_SIGMA, g_array_index(articles, Item, i).sigma,
                           COLUMN_LAYER_SIGAST, g_array_index(articles, Item, i).sigast,
                           COLUMN_LAYER_MATERIAL, g_array_index(articles, Item, i).material,
                           -1);
    }

    //g_free(item.material);

    return GTK_TREE_MODEL(model);
}

static void
remove_item_cb(GtkWidget *button, gpointer data)
{
    GtkTreeIter iter;
    GtkTreeView *treeview = (GtkTreeView *)data;
    GtkTreeModel *model = gtk_tree_view_get_model(treeview);
    GtkTreeSelection *selection = gtk_tree_view_get_selection(treeview);

    if (gtk_tree_selection_get_selected(selection, NULL, &iter)) {
        gint i;
        GtkTreePath *path;

        path = gtk_tree_model_get_path(model, &iter);
        i = gtk_tree_path_get_indices(path)[0];
        gtk_list_store_remove(GTK_LIST_STORE(model), &iter);

        g_free(g_array_index(articles, Item, i).material);
        g_array_remove_index(articles, i);

        gtk_tree_selection_select_iter(selection, &iter);

        gtk_tree_path_free(path);

        index_items(model);

        g_signal_emit_by_name(treeview, "cell-edited");
    }
}

/*
static void
focus_next_cell_cb(GtkWidget *button, gpointer data)
{
    GtkTreeIter iter;
    GtkTreeView *treeview = (GtkTreeView *)data;
    GtkTreeModel *model = gtk_tree_view_get_model(treeview);
    GtkTreeSelection *selection = gtk_tree_view_get_selection(treeview);
    GtkTreeViewColumn *column = gtk_tree_view_get_column(treeview, COLUMN_LAYER_EPSILON);
    GtkTreePath *path;
    GtkTreeViewColumn *focus_column;

    gtk_tree_view_get_cursor (treeview, &path, &focus_column);
    gtk_tree_model_get_iter(model, &iter, path);
    //gtk_tree_model_iter_next(model, &iter);
    path = gtk_tree_model_get_path(model, &iter);
    gint row = gtk_tree_path_get_indices(path)[0];


    if (column_to_edit_after_tab == -1) {
        row_to_edit_after_tab = row;
        if (focus_column == gtk_tree_view_get_column(treeview, COLUMN_LAYER_Z))
            column_to_edit_after_tab = COLUMN_LAYER_EPSILON;
        else if (focus_column == gtk_tree_view_get_column(treeview, COLUMN_LAYER_EPSILON))
            column_to_edit_after_tab = COLUMN_LAYER_MU;
        else if (focus_column == gtk_tree_view_get_column(treeview, COLUMN_LAYER_MU))
            column_to_edit_after_tab = COLUMN_LAYER_SIGMA;
        else if (focus_column == gtk_tree_view_get_column(treeview, COLUMN_LAYER_SIGMA))
            column_to_edit_after_tab = COLUMN_LAYER_SIGAST;
        else if (focus_column == gtk_tree_view_get_column(treeview, COLUMN_LAYER_SIGAST))
            column_to_edit_after_tab = COLUMN_LAYER_MATERIAL;
        else if (focus_column == gtk_tree_view_get_column(treeview, COLUMN_LAYER_MATERIAL))
            column_to_edit_after_tab = COLUMN_LAYER_Z;
    }

    //gtk_tree_path_next(path);
    //gtk_tree_path_down(path);

    //gtk_tree_view_set_cursor_on_cell(treeview, path, NULL, NULL, TRUE);

    focus_column == gtk_tree_view_get_column(treeview, COLUMN_LAYER_EPSILON);

    GList *cell_layout, *l, *cl1, *cl2;
    cl1 = gtk_cell_layout_get_cells(GTK_CELL_LAYOUT(focus_column));
    cl2 = gtk_tree_view_column_get_cell_renderers(column);
    cell_layout = column->cell_list;
    gint nnn = g_list_length (cl1);

    for (l = cl1; l != NULL; l = l->next) {
        GtkCellRenderer* element_data = l->data;

    }
}
*/

static void
text_editing_started (GtkCellRenderer *cell,
    GtkCellEditable *editable,
    const gchar     *path,
    gpointer         data)
{
    if (GTK_IS_ENTRY (editable))

    {
        GtkEntry *entry = GTK_ENTRY (editable);

        /* ... create a GtkEntryCompletion */

        //gtk_entry_set_completion (entry, completion);
    }
}

gboolean
foreach_func_get_data(GtkTreeModel *model,
                      GtkTreePath  *path,
                      GtkTreeIter  *iter,
                      gpointer      data)
{
    LTV_Data *ltvd = (LTV_Data*)data;

    if (ltvd->nlayers >= NUM_LSMP)
        return TRUE;

    gtk_tree_model_get(model, iter,
                       COLUMN_LAYER_Z, &ltvd->zpos[ltvd->nlayers],
                       COLUMN_LAYER_EPSILON, &ltvd->epsilon[ltvd->nlayers],
                       COLUMN_LAYER_MU, &ltvd->mu[ltvd->nlayers],
                       COLUMN_LAYER_SIGMA, &ltvd->sigma[ltvd->nlayers],
                       COLUMN_LAYER_SIGAST, &ltvd->sigast[ltvd->nlayers],
                       COLUMN_LAYER_MATERIAL, &ltvd->material[ltvd->nlayers],
                       -1);

    ltvd->nlayers++;

    return FALSE; /* do not stop walking the store, call us with next row */
}

void
get_ltv_items(TreeViewWrapper *treeview_wrapper, gpointer *data)
{
    //GtkTreeModel *model = gtk_tree_view_get_model(treeview);
    GtkTreeModel *model = treeview_wrapper->items_model;
    LTV_Data *ltvd = (LTV_Data*)data;

    memset(ltvd, 0, sizeof(LTV_Data));
    gtk_tree_model_foreach(model, foreach_func_get_data, data);
}

void
set_ltv_item_data(TreeViewWrapper *treeview_wrapper, gpointer *data)
{
    //GtkTreeModel *model = gtk_tree_view_get_model(treeview);
    GtkTreeModel *model = treeview_wrapper->items_model;
    LTV_Data *ltvd = (LTV_Data*)data;

    memset(ltvd, 0, sizeof(LTV_Data));
    gtk_tree_model_foreach(model, foreach_func_get_data, data);
}

gboolean
foreach_func_set_data(GtkTreeModel *model,
                      GtkTreePath  *path,
                      GtkTreeIter  *iter,
                      gpointer      data)
{
    LTV_Data *ltvd = (LTV_Data*)data;

    if (ltvd->nlayers >= NUM_LSMP)
        return TRUE;

    gtk_tree_model_get(model, iter,
                       COLUMN_LAYER_Z, &ltvd->zpos[ltvd->nlayers],
                       COLUMN_LAYER_EPSILON, &ltvd->epsilon[ltvd->nlayers],
                       COLUMN_LAYER_MU, &ltvd->mu[ltvd->nlayers],
                       COLUMN_LAYER_SIGMA, &ltvd->sigma[ltvd->nlayers],
                       COLUMN_LAYER_SIGAST, &ltvd->sigast[ltvd->nlayers],
                       COLUMN_LAYER_MATERIAL, &ltvd->material[ltvd->nlayers],
                       -1);

    ltvd->nlayers++;

    return FALSE; /* do not stop walking the store, call us with next row */
}


static void
cell_edited(GtkCellRendererText *cell,
            const gchar         *path_string,
            const gchar         *new_text,
            gpointer             user_data)
{
    //GtkTreeModel *model = (GtkTreeModel *)data;
    //TreeViewData *treeviewdata = (TreeViewData *)data;
    TreeViewWrapper *treeviewdata = (TreeViewWrapper *)user_data;
    GtkTreeModel *model = treeviewdata->items_model;
    GtkTreePath *path = gtk_tree_path_new_from_string(path_string);
    GtkTreeIter iter;
    gint i = 0;

    gint column = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(cell), "column"));

    gtk_tree_model_get_iter(model, &iter, path);

    i = gtk_tree_path_get_indices(path)[0];

    switch (column) {
        case COLUMN_LAYER_NUMBER: {
            g_array_index(articles, Item, i).number = atoi(new_text);

            gtk_list_store_set(GTK_LIST_STORE(model), &iter, column, g_array_index(articles, Item, i).number, -1);
        }
        break;

        case COLUMN_LAYER_Z: {
            g_array_index(articles, Item, i).zpos = atoi(new_text);

            gtk_list_store_set(GTK_LIST_STORE(model), &iter, column, g_array_index(articles, Item, i).zpos, -1);
        }
        break;

        case COLUMN_LAYER_EPSILON: {
            g_array_index(articles, Item, i).epsilon = atof(new_text);

            gtk_list_store_set(GTK_LIST_STORE(model), &iter, column, g_array_index(articles, Item, i).epsilon, -1);
        }
        break;

        case COLUMN_LAYER_MU: {
            g_array_index(articles, Item, i).mu = atof(new_text);

            gtk_list_store_set(GTK_LIST_STORE(model), &iter, column, g_array_index(articles, Item, i).mu, -1);
        }
        break;

        case COLUMN_LAYER_SIGMA: {
            g_array_index(articles, Item, i).sigma = atof(new_text);

            gtk_list_store_set(GTK_LIST_STORE(model), &iter, column, g_array_index(articles, Item, i).sigma, -1);
        }
        break;

        case COLUMN_LAYER_SIGAST: {
            g_array_index(articles, Item, i).sigast = atof(new_text);

            gtk_list_store_set(GTK_LIST_STORE(model), &iter, column, g_array_index(articles, Item, i).sigast, -1);
        }
        break;

        case COLUMN_LAYER_MATERIAL: {
            //strncpy(g_array_index(articles, Item, i).material, new_text, sizeof(g_array_index(articles, Item, i).material));
            strcpy(g_array_index(articles, Item, i).material, new_text);

            gtk_list_store_set(GTK_LIST_STORE(model), &iter, column, g_array_index(articles, Item, i).material, -1);
        }
        break;
    }

    gtk_tree_path_free(path);


    g_signal_emit_by_name(treeviewdata->parent_treeview, "cell-edited");
}

/*
static void
cell_changed(GtkCellRendererCombo *cell,
             const gchar          *path_string,
             GtkTreeIter          *new_iter,
             CellRendererData*    user_data)
{
    //GtkTreeModel *model = (GtkTreeModel *)data;
    //TreeViewData *treeviewdata = (TreeViewData *)data;
    TreeViewWrapper *treeviewdata = (TreeViewWrapper *)user_data;
    GtkTreeModel *model = treeviewdata->items_model;
    GtkTreePath *path = gtk_tree_path_new_from_string(path_string);
    GtkTreeIter iter;
    gint i = 0;

    gint column = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(cell), "column"));

    gtk_tree_model_get_iter(model, &iter, path);

    i = gtk_tree_path_get_indices(path)[0];


    gchar *buf;

    gtk_tree_model_get(user_data->model_combo, new_iter, 0, &buf, -1);

    strncpy(g_array_index(articles, Item, i).material, buf, sizeof(g_array_index(articles, Item, i).material));

    gtk_list_store_set(GTK_LIST_STORE(user_data->model_treeview), &iter, column, g_array_index(articles, Item, i).material, -1);

    gtk_tree_path_free(path);


    g_signal_emit_by_name(user_data->model_treeview, "cell-edited");
}
*/

void
material_cell_data_func(GtkTreeViewColumn *col,
                        GtkCellRenderer   *renderer,
                        GtkTreeModel      *model,
                        GtkTreeIter       *iter,
                        gpointer           user_data)
{   
    //
    //gint column = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(renderer), "column"));
    glong column = (glong)user_data;
    GtkTreePath *path = gtk_tree_model_get_path(model, iter);
    gint row = gtk_tree_path_get_indices(path)[0];
    //

    gchar *buf;

    if (column == COLUMN_LAYER_MATERIAL) {
        gtk_tree_model_get(model, iter, column, &buf, -1);

        g_object_set(renderer, "text", buf, NULL);
    } else {
        gdouble         val;
        gchar           text[20];
        gchar           format_text[10];
        gint            digits;
        gboolean        sensitive;
        GtkAdjustment*  adjustment;        

        gtk_tree_model_get(model, iter, column, &val, -1);        

        g_object_get(renderer,
                     "digits", &digits,
                     "adjustment", &adjustment,
                     NULL);

        // create format string
        g_snprintf(format_text, sizeof(format_text), """%%.%df", digits);

        // check value range
        val = MAX(val, gtk_adjustment_get_lower(adjustment));
        val = MIN(val, gtk_adjustment_get_upper(adjustment));
        g_snprintf(text, sizeof(text), format_text, val);

        // set text to renderer
        //g_object_set(renderer, "text", text, NULL);

        gtk_tree_model_get(model, iter, COLUMN_LAYER_MATERIAL, &buf, -1);
        sensitive = (buf != NULL && strcmp(buf, SOURCE_LAYERED_MATERIAL) == 0);

        g_object_set(renderer, 
                     "text", text,
                     "sensitive", sensitive, 
                     "editable", sensitive, NULL);
    }       

    /*if (row == row_to_edit_after_tab && column == column_to_edit_after_tab) {
        g_object_set(renderer, "editing", TRUE, NULL);
        row_to_edit_after_tab = -1;
        column_to_edit_after_tab = -1;
    }*/
}

/*
static void 
fill_material_combo_list(GtkListStore **combo_list)
{
    GDir *dir;
    GError *error;
    const gchar *filename;
    gchar *path;

    path = get_spectra_path(NULL);

    dir = g_dir_open(path, 0, &error);
    gtk_list_store_insert_with_values(*combo_list, NULL, -1, 0, SOURCE_LAYERED_MATERIAL, -1);
    while ((filename = g_dir_read_name(dir))) {
        //printf("%s\n", filename);
        gtk_list_store_insert_with_values(*combo_list, NULL, -1, 0, g_strdup(filename), -1);
    }
}
*/

static void
add_columns(GtkTreeView  *treeview,
            TreeViewWrapper *treeviewdata
            /*TreeViewData *treeviewdata*/
            /*GtkTreeModel *items_model*/)
{
    GtkCellRenderer *renderer;
    GtkObject *adjustment;
    GtkTreeViewColumn *column;

    /* number column */
    renderer = gtk_cell_renderer_text_new();
    g_object_set(renderer, "width", 20, NULL);
    g_object_set_data(G_OBJECT(renderer), "column", GINT_TO_POINTER(COLUMN_LAYER_NUMBER));

    gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(treeview), -1, "No.", renderer, "text", COLUMN_LAYER_NUMBER, NULL);

    column = gtk_tree_view_get_column(GTK_TREE_VIEW(treeview), COLUMN_LAYER_NUMBER);
    gtk_tree_view_column_set_sizing(GTK_TREE_VIEW_COLUMN(column), GTK_TREE_VIEW_COLUMN_FIXED);
    gtk_tree_view_column_set_fixed_width(GTK_TREE_VIEW_COLUMN(column), 25);
    gtk_tree_view_column_set_alignment(GTK_TREE_VIEW_COLUMN(column), 0.0);
    //gtk_tree_view_column_set_clickable(GTK_TREE_VIEW_COLUMN(column), TRUE);

    /* Z column */
    renderer = gtk_cell_renderer_spin_new();
    g_object_set(renderer, "editable", TRUE, NULL);
    g_signal_connect(renderer, "edited", G_CALLBACK(cell_edited), treeviewdata);
    g_object_set_data(G_OBJECT(renderer), "column", GINT_TO_POINTER(COLUMN_LAYER_Z));

    adjustment = gtk_adjustment_new(0, 0, 1000, 1, 10, 0);
    g_object_set(renderer,
                 "editable", TRUE,
                 //"width", 50,
                 "adjustment", adjustment,
                 NULL);

    gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(treeview), -1, "Z [vx]", renderer, "text", COLUMN_LAYER_Z, NULL);

    column = gtk_tree_view_get_column(GTK_TREE_VIEW(treeview), COLUMN_LAYER_Z);
    gtk_tree_view_column_set_sizing(GTK_TREE_VIEW_COLUMN(column), GTK_TREE_VIEW_COLUMN_FIXED);
    gtk_tree_view_column_set_fixed_width(GTK_TREE_VIEW_COLUMN(column), 40);
    gtk_tree_view_column_set_alignment(GTK_TREE_VIEW_COLUMN(column), 0.0);
    //gtk_tree_view_column_set_clickable(GTK_TREE_VIEW_COLUMN(column), TRUE);

    /* Material properties (ε, μ, σ, σ*) */
    /* ε column */
    renderer = gtk_cell_renderer_spin_new();
    g_object_set(renderer, "editable", TRUE, NULL);
    g_signal_connect(renderer, "edited", G_CALLBACK(cell_edited), treeviewdata);
    g_object_set_data(G_OBJECT(renderer), "column", GINT_TO_POINTER(COLUMN_LAYER_EPSILON));

    adjustment = gtk_adjustment_new(0, 1, 100, 0.1, 10, 0);
    g_object_set(renderer,
                 "editable", TRUE,
                 "digits", 2,
                 "climb-rate", 0.1,
                 "adjustment", adjustment,
                 //"width", 50,
                 "alignment", PANGO_ALIGN_CENTER,
                 NULL);

    gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(treeview), -1, "ε", renderer, "text", COLUMN_LAYER_EPSILON, NULL);

    column = gtk_tree_view_get_column(GTK_TREE_VIEW(treeview), COLUMN_LAYER_EPSILON);
    gtk_tree_view_column_set_cell_data_func(column, renderer, material_cell_data_func, (gpointer)COLUMN_LAYER_EPSILON, NULL);
    gtk_tree_view_column_set_sizing(GTK_TREE_VIEW_COLUMN(column), GTK_TREE_VIEW_COLUMN_FIXED);
    gtk_tree_view_column_set_fixed_width(GTK_TREE_VIEW_COLUMN(column), 50);
    gtk_tree_view_column_set_alignment(GTK_TREE_VIEW_COLUMN(column), 0.0);

    /* μ column */
    renderer = gtk_cell_renderer_spin_new();
    g_object_set(renderer, "editable", TRUE, NULL);
    g_signal_connect(renderer, "edited", G_CALLBACK(cell_edited), treeviewdata);
    g_object_set_data(G_OBJECT(renderer), "column", GINT_TO_POINTER(COLUMN_LAYER_MU));

    adjustment = gtk_adjustment_new(0, 1, 100, 0.1, 10, 0);
    g_object_set(renderer,
                 "editable", TRUE,
                 "digits", 2,
                 "climb-rate", 0.1,
                 "adjustment", adjustment,
                 //"width", 50,
                 NULL);

    gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(treeview), -1, "μ", renderer, "text", COLUMN_LAYER_MU, NULL);

    column = gtk_tree_view_get_column(GTK_TREE_VIEW(treeview), COLUMN_LAYER_MU);
    gtk_tree_view_column_set_cell_data_func(column, renderer, material_cell_data_func, (gpointer)COLUMN_LAYER_MU, NULL);
    gtk_tree_view_column_set_sizing(GTK_TREE_VIEW_COLUMN(column), GTK_TREE_VIEW_COLUMN_FIXED);
    gtk_tree_view_column_set_fixed_width(GTK_TREE_VIEW_COLUMN(column), 50);
    gtk_tree_view_column_set_alignment(GTK_TREE_VIEW_COLUMN(column), 0.0);

    /* σ column */
    renderer = gtk_cell_renderer_spin_new();
    g_object_set(renderer, "editable", TRUE, NULL);
    g_signal_connect(renderer, "edited", G_CALLBACK(cell_edited), treeviewdata);
    g_object_set_data(G_OBJECT(renderer), "column", GINT_TO_POINTER(COLUMN_LAYER_SIGMA));

    adjustment = gtk_adjustment_new(0, 0, 9999999, 0.1, 10, 0);
    g_object_set(renderer,
                 "editable", TRUE,
                 "digits", 1,
                 "climb-rate", 0.1,
                 "adjustment", adjustment,
                 //"width", 50,
                 NULL);

    gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(treeview), -1, "σ", renderer, "text", COLUMN_LAYER_SIGMA, NULL);

    column = gtk_tree_view_get_column(GTK_TREE_VIEW(treeview), COLUMN_LAYER_SIGMA);
    gtk_tree_view_column_set_cell_data_func(column, renderer, material_cell_data_func, (gpointer)COLUMN_LAYER_SIGMA, NULL);
    gtk_tree_view_column_set_sizing(GTK_TREE_VIEW_COLUMN(column), GTK_TREE_VIEW_COLUMN_FIXED);
    gtk_tree_view_column_set_fixed_width(GTK_TREE_VIEW_COLUMN(column), 50);
    gtk_tree_view_column_set_alignment(GTK_TREE_VIEW_COLUMN(column), 0.0);

    /* σ* column */    
    renderer = gtk_cell_renderer_spin_new();
    g_object_set(renderer, "editable", TRUE, NULL);
    g_signal_connect(renderer, "edited", G_CALLBACK(cell_edited), treeviewdata);
    g_object_set_data(G_OBJECT(renderer), "column", GINT_TO_POINTER(COLUMN_LAYER_SIGAST));

    adjustment = gtk_adjustment_new(0, 0, 9999999, 0.1, 10, 0);
    g_object_set(renderer,
                 "editable", TRUE,
                 "digits", 1,
                 "climb-rate", 0.1,
                 "adjustment", adjustment,
                 //"width", 50,
                 NULL);

    gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(treeview), -1, "σ*", renderer, "text", COLUMN_LAYER_SIGAST, NULL);

    column = gtk_tree_view_get_column(GTK_TREE_VIEW(treeview), COLUMN_LAYER_SIGAST);
    gtk_tree_view_column_set_cell_data_func(column, renderer, material_cell_data_func, (gpointer)COLUMN_LAYER_SIGAST, NULL);
    gtk_tree_view_column_set_sizing(GTK_TREE_VIEW_COLUMN(column), GTK_TREE_VIEW_COLUMN_FIXED);
    gtk_tree_view_column_set_fixed_width(GTK_TREE_VIEW_COLUMN(column), 50);
    gtk_tree_view_column_set_alignment(GTK_TREE_VIEW_COLUMN(column), 0.0);

    /* Material column */
    GtkListStore *combo_list = gtk_list_store_new(1, G_TYPE_STRING);   // list for combo box
    xsv_fill_material_combo_list(&combo_list);

    treeviewdata->combo_model = GTK_TREE_MODEL(combo_list);

    renderer = gtk_cell_renderer_combo_new();
    g_object_set(renderer, "model", combo_list, "text-column", 0, "editable", TRUE, "has-entry", FALSE, NULL);      
    g_object_set_data(G_OBJECT(renderer), "column", GINT_TO_POINTER(COLUMN_LAYER_MATERIAL));
    g_signal_connect(renderer, "edited", G_CALLBACK(cell_edited), treeviewdata);

    gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(treeview), -1, "Material", renderer, "text", COLUMN_LAYER_MATERIAL, NULL);

    column = gtk_tree_view_get_column(GTK_TREE_VIEW(treeview), COLUMN_LAYER_MATERIAL);
    gtk_tree_view_column_set_cell_data_func(column, renderer, material_cell_data_func, (gpointer)COLUMN_LAYER_MATERIAL, NULL);
    gtk_tree_view_column_set_sizing(GTK_TREE_VIEW_COLUMN(column), GTK_TREE_VIEW_COLUMN_FIXED);
    gtk_tree_view_column_set_fixed_width(GTK_TREE_VIEW_COLUMN(column), 50);
    gtk_tree_view_column_set_alignment(GTK_TREE_VIEW_COLUMN(column), 0.0);    
}

gboolean
on_key_press(GtkWidget *widget, GdkEventKey *event, gpointer user_data)
{
    TreeViewWrapper *tvw = (TreeViewWrapper *)user_data;
    GtkTreeView *tv = GTK_TREE_VIEW(tvw->parent_treeview);

    switch (event->keyval)
    {
    case GDK_Insert:
        insert_item_cb(widget, tvw);
        break;
    case GDK_Delete:
        remove_item_cb(widget, tv);
        break;
//    case GDK_Tab:
//        focus_next_cell_cb(widget, tv);
//        break;
//    default:
        return FALSE;
    }

    return FALSE;
}

// treeview is created 
// returns treeview wrapper 

void
create_layer_treeview(TreeViewWrapper **treeview_wrapper, GtkWidget **packet)
{
    GtkWidget *vbox;
    GtkWidget *hbox;
    GtkWidget *sw;
    GtkWidget *tv;
    GtkWidget *button;
    GtkTreeModel *items_model;
    //GSvitTreeview *treeviewdata;
    //TreeViewData *treeviewdata;

    vbox = gtk_vbox_new(FALSE, 5);

    sw = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(sw), GTK_SHADOW_ETCHED_IN);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);    
    gtk_box_pack_start(GTK_BOX(vbox), sw, TRUE, TRUE, 0);

    /* create models */
    items_model = create_items_model();

    /* create tree view */
    tv = gtk_tree_view_new_with_model(items_model);
    gtk_tree_view_set_rules_hint(GTK_TREE_VIEW(tv), TRUE);
    gtk_tree_selection_set_mode(gtk_tree_view_get_selection(GTK_TREE_VIEW(tv)), GTK_SELECTION_SINGLE);

    //treeviewdata = treeviewdata_new();
    //treeviewdata = g_object_new(G_TYPE_OBJECT, NULL);

    treeviewwrapper = g_object_new(GSVIT_TYPE_TREEVIEW, NULL);

    g_object_ref(items_model);
    treeviewwrapper->items_model = items_model;
    treeviewwrapper->parent_treeview = tv;
    add_columns(GTK_TREE_VIEW(tv), treeviewwrapper/*, numbers_model*/);

    //g_object_unref(items_model);
    //g_object_unref(treeviewdata);

    gtk_container_add(GTK_CONTAINER(sw), tv);

    /* some buttons */
    hbox = gtk_hbox_new(TRUE, 4);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    button = gtk_button_new_with_label("Insert layer");
    g_signal_connect(button, "clicked", G_CALLBACK(insert_item_cb), treeviewwrapper /*items_model*/);
    gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, TRUE, 0);

    button = gtk_button_new_with_label("Remove layer");
    g_signal_connect(button, "clicked", G_CALLBACK(remove_item_cb), tv);
    gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, TRUE, 0);

    g_signal_connect (G_OBJECT(sw), "key_press_event", G_CALLBACK(on_key_press), treeviewwrapper);

    //*treeview = tv;
    *treeview_wrapper = treeviewwrapper;
    *packet = vbox;
}

void
destroy_layer_treeview()
{
    g_object_unref(treeviewwrapper);
}
