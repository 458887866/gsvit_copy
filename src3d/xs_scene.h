
#ifndef SCENE_H
#define SCENE_H

#include "global.h"

gboolean scene_setup(GtkWidget *view_scene);

void scene_reset_cb(GtkWidget *widget, XGControls *xgc);
gboolean scene_configure_cb(GtkWidget *da, GdkEventConfigure *event, XGControls *xgc);
gboolean scene_expose_cb(GtkWidget *da, GdkEventExpose *event, XGControls *xgc);
void scene_zoom_in_cb(GtkWidget *widget, XGControls *xgc);
void scene_zoom_out_cb(GtkWidget *widget, XGControls *xgc);
void scene_move_view_right_cb(GtkWidget *widget, XGControls *xgc);
void scene_move_view_left_cb(GtkWidget *widget, XGControls *xgc);
void scene_move_view_up_cb(GtkWidget *widget, XGControls *xgc);
void scene_move_view_down_cb(GtkWidget *widget, XGControls *xgc);
gboolean scene_configure_cb(GtkWidget *da, GdkEventConfigure *event, XGControls *xgc);
gboolean scene_expose_cb(GtkWidget *da, GdkEventExpose *event, XGControls *xgc);
gboolean scene_button_press_cb(GtkWidget* widget, GdkEventButton * event, XGControls *xgc);
gboolean scene_button_release_cb(GtkWidget* widget, GdkEventButton * event, XGControls *xgc);
gboolean scene_motion_notify_cb(GtkWidget *widget, GdkEventMotion *event, XGControls *xgc);
gboolean scene_scroll_cb(GtkWidget* widget, GdkEventScroll * event, XGControls *xgc);



#endif  /* SCENE_H */