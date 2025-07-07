
#ifndef LAYER_TREEVIEW_H
#define LAYER_TREEVIEW_H

#include "settings.h"

typedef struct _LTV_Data {
    guint    nlayers;
    gint     zpos[NUM_LSMP];
    gdouble  epsilon[NUM_LSMP];
    gdouble  sigma[NUM_LSMP];
    gdouble  mu[NUM_LSMP];
    gdouble  sigast[NUM_LSMP];
    gchar*   material[NUM_LSMP];
} LTV_Data;


//typedef struct _TreeViewData       TreeViewData;

typedef struct
{
    GtkTreeModel    *items_model;
    GtkWidget       *parent_treeview;
} TreeViewData;



/* GSvit tree view object */
#define GSVIT_TYPE_TREEVIEW            (tree_view_wrapper_get_type())
#define GSVIT_TREEVIEW(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj), GSVIT_TYPE_TREEVIEW, TreeViewWrapper))
#define GSVIT_TREEVIEW_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass), GSVIT_TYPE_TREEVIEW, TreeViewWrapperClass))
#define GSVIT_IS_TREEVIEW(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj), GSVIT_TYPE_TREEVIEW))
#define GSVIT_IS_TREEVIEW_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass), GSVIT_TYPE_TREEVIEW))
#define GSVIT_TREEVIEW_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS((obj), GSVIT_TYPE_TREEVIEW, TreeViewWrapperClass))

typedef struct _TreeViewWrapper      TreeViewWrapper;
typedef struct _TreeViewWrapperClass TreeViewWrapperClass;

struct _TreeViewWrapper {
    GObject         parent_instance;

    GtkTreeModel    *items_model;
    GtkTreeModel    *combo_model;
    GtkWidget       *parent_treeview;    
};

struct _TreeViewWrapperClass {
    GObjectClass parent_class;
};

GType tree_view_wrapper_get_type (void) G_GNUC_CONST;
/* */


void create_layer_treeview(TreeViewWrapper **treeview_wrapper, GtkWidget **packet);
void destroy_layer_treeview();
void get_ltv_items(TreeViewWrapper *treeview_wrapper, gpointer* data);
void set_ltv_items(TreeViewWrapper *treeview_wrapper, gpointer* data);

#endif /* LAYER_TREEVIEW_H */