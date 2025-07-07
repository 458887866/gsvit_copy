
#include <math.h>
#include "xs_dialogs.h"
#include "constants.h"
#include <string.h>


#define MAXPOS 10000

typedef struct {
    gchar *gsvitloc;
    GtkWidget *gsvitlab;
    gchar *gwyddloc;
    GtkWidget *gwyddlab;
    GtkWidget *toplevel;
} XGOptData;



void
xsv_fill_material_combo(GtkWidget **combo)
{
    GDir *dir;
    GError *error;
    const gchar *filename;
    gchar * path;

    path = get_spectra_path(NULL);

    dir = g_dir_open(path, 0, &error);
    //gtk_list_store_insert_with_values(*combo_list, NULL, -1, 0, SOURCE_LAYERED_MATERIAL, -1);
    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(*combo), MAT_PROP_MATERIAL_NONE);
    while ((filename = g_dir_read_name(dir))) {
        //printf("%s\n", filename);
        //gtk_list_store_insert_with_values(*combo_list, NULL, -1, 0, g_strdup(filename), -1);

        if (FALSE == (strcmp(filename, "Makefile.am") == 0 || strcmp(filename, "makefile.am") == 0))
            gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(*combo), filename);
        //gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(*combo), g_strdup(filename));
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(*combo), MAT_PROP_MATERIAL_NONE_INDEX);

    g_dir_close(dir);
}

void
xsv_fill_material_combo_list(GtkListStore **combo_list)
{
    GDir *dir;
    GError *error;
    const gchar *filename;
    gchar * path;

    path = get_spectra_path(NULL);

    dir = g_dir_open(path, 0, &error);
    gtk_list_store_insert_with_values(*combo_list, NULL, -1, 0, SOURCE_LAYERED_MATERIAL, -1);
    while ((filename = g_dir_read_name(dir))) {
        //printf("%s\n", filename);
        //gtk_list_store_insert_with_values(*combo_list, NULL, -1, 0, g_strdup(filename), -1);

        if (FALSE == (strcmp(filename, "Makefile.am") == 0 || strcmp(filename, "makefile.am") == 0))
            gtk_list_store_insert_with_values(*combo_list, NULL, -1, 0, filename, -1);
        //gtk_list_store_insert_with_values(*combo_list, NULL, -1, 0, g_strdup(filename), -1);
    }
    //gtk_combo_box_set_active(GTK_COMBO_BOX(*combo_list), SOURCE_LAYERED_MATERIAL_INDEX);

    g_dir_close(dir);
}


void choose_gsvit_cb(GtkWidget *widget, XGOptData *xo)
{
    GtkWidget *dialog;
    GtkFileFilter *filter = gtk_file_filter_new();
    gchar* dir = NULL;

    dialog = gtk_file_chooser_dialog_new("Choose GSvit location:",
                                         GTK_WINDOW(xo->toplevel),
                                         GTK_FILE_CHOOSER_ACTION_OPEN,
                                         GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                         GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                                         NULL);
    gtk_file_filter_set_name(filter, "GSvit executable");
    gtk_file_filter_add_pattern(filter, "gsvit*");
    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filter);
    if (xo->gsvitloc) {
        dir = g_path_get_dirname(xo->gsvitloc);
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), dir);
        g_free(dir);
    }

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        xo->gsvitloc = g_strdup(gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog)));
        gtk_entry_set_text(GTK_ENTRY(xo->gsvitlab), xo->gsvitloc);
    }
    gtk_widget_destroy(dialog);

}

void choose_gwydd_cb(GtkWidget *widget, XGOptData *xo)
{
    GtkWidget *dialog;
    GtkFileFilter *filter = gtk_file_filter_new();
    gchar* dir = NULL;

    dialog = gtk_file_chooser_dialog_new("Choose Gwyddion location:",
                                         GTK_WINDOW(xo->toplevel),
                                         GTK_FILE_CHOOSER_ACTION_OPEN,
                                         GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                         GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                                         NULL);
    gtk_file_filter_set_name(filter, "Gwyddion executable");
    gtk_file_filter_add_pattern(filter, "gwyddion*");
    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filter);
    if (xo->gwyddloc) {
        dir = g_path_get_dirname(xo->gwyddloc);
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), dir);
        g_free(dir);
    }

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        xo->gwyddloc = g_strdup(gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog)));
        gtk_entry_set_text(GTK_ENTRY(xo->gwyddlab), xo->gwyddloc);
    }
    gtk_widget_destroy(dialog);

}


void xsv_get_preferences(XGControls *xgc)
{
    GtkWidget *dialog, *label, *button;
    GtkWidget *table;
    gint response, row;
    XGOptData xo;
    gchar *data;
    gsize size;
    GError *error = NULL;
    gchar buff[256];

    if (xgc->data.gsvit_location)
        xo.gsvitloc = g_strdup(xgc->data.gsvit_location);
    else
        xo.gsvitloc = g_strdup(GENERAL_FILENAME);
    if (xgc->data.gwyddion_location)
        xo.gwyddloc = g_strdup(xgc->data.gwyddion_location);
    else
        xo.gwyddloc = g_strdup(GENERAL_FILENAME);


    /* Create the widgets */
    dialog = gtk_dialog_new_with_buttons("Preferences",
                                         GTK_WINDOW(xgc->toplevel),
                                         GTK_DIALOG_DESTROY_WITH_PARENT,
                                         GTK_STOCK_OK,
                                         GTK_RESPONSE_ACCEPT,
                                         GTK_STOCK_CANCEL,
                                         GTK_RESPONSE_REJECT,
                                         NULL);
    gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);

    xo.toplevel = dialog;

    table = gtk_table_new(4, 4, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(table), 2);
    gtk_table_set_col_spacings(GTK_TABLE(table), 6);
    gtk_container_set_border_width(GTK_CONTAINER(table), 4);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), table,
                       FALSE, FALSE, 4);



    row = 0;
    label = gtk_label_new("GSvit executable:");
    gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table), label, 0, 1, row, row + 1, GTK_FILL, 0, 0, 0);


    xo.gsvitlab = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xo.gsvitlab), 50);
    gtk_entry_set_text(GTK_ENTRY(xo.gsvitlab), xo.gsvitloc);
    gtk_entry_set_editable(GTK_ENTRY(xo.gsvitlab), FALSE);
    gtk_widget_set_can_focus(xo.gsvitlab, FALSE);
    gtk_table_attach(GTK_TABLE(table), xo.gsvitlab, 1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);


    button = gtk_button_new_with_label("Change...");
    gtk_table_attach(GTK_TABLE(table), button, 2, 3, row, row + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect(button, "clicked", G_CALLBACK(choose_gsvit_cb), &xo);
    row++;

    label = gtk_label_new("Gwyddion executable:");
    gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
    gtk_table_attach(GTK_TABLE(table), label, 0, 1, row, row + 1, GTK_FILL, 0, 0, 0);

    xo.gwyddlab = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(xo.gwyddlab), 50);
    gtk_entry_set_text(GTK_ENTRY(xo.gwyddlab), xo.gwyddloc);
    gtk_entry_set_editable(GTK_ENTRY(xo.gwyddlab), FALSE);
    gtk_widget_set_can_focus(xo.gwyddlab, FALSE);
    gtk_table_attach(GTK_TABLE(table), xo.gwyddlab, 1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);


    button = gtk_button_new_with_label("Change...");
    gtk_table_attach(GTK_TABLE(table), button, 2, 3, row, row + 1, GTK_FILL, 0, 0, 0);
    g_signal_connect(button, "clicked", G_CALLBACK(choose_gwydd_cb), &xo);



    gtk_widget_show_all(dialog);
    response = gtk_dialog_run(GTK_DIALOG(dialog));

    if (response == GTK_RESPONSE_ACCEPT) {

        xgc->data.gsvit_location = g_strdup(xo.gsvitloc);
        xgc->data.gwyddion_location = g_strdup(xo.gwyddloc);

        /*write all this to configuration file as well*/
        g_key_file_set_string(xgc->keyfile, "Locations", "GSvit", xgc->data.gsvit_location);
        g_key_file_set_string(xgc->keyfile, "Locations", "Gwyddion", xgc->data.gwyddion_location);

        data = g_key_file_to_data(xgc->keyfile, &size, &error);

        if (!data) {
            g_snprintf(buff, sizeof(buff), "Cannot save settings: %s", error->message);
            gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
            g_error_free(error);
        } else {
            if (!g_file_set_contents(xgc->config_path, data, size, &error)) {
                g_snprintf(buff, sizeof(buff), "Cannot save settings: %s", error->message);
                gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
                g_error_free(error);

            } else {
                g_snprintf(buff, sizeof(buff), "Settings saved");
                gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);
            }
        }
    }

    gtk_widget_destroy(dialog);
}





gboolean xsv_get_basic(XGControls *xgc)
{
    GtkWidget *dialog, *label;
    GtkObject *nsteps;
    GtkObject *verbose;
    GtkObject *nthreads;
    GtkObject *ugpu;
    GtkWidget *gpu;
    GtkWidget *table, *spin;
    gint response, row;
    gint ugval;

    /* Create the widgets */
    dialog = gtk_dialog_new_with_buttons("Basic computation parameters",
                                         GTK_WINDOW(xgc->toplevel),
                                         GTK_DIALOG_DESTROY_WITH_PARENT,
                                         GTK_STOCK_OK,
                                         GTK_RESPONSE_ACCEPT,
                                         GTK_STOCK_CANCEL,
                                         GTK_RESPONSE_REJECT,
                                         NULL);
    gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
    table = gtk_table_new(4, 4, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(table), 2);
    gtk_table_set_col_spacings(GTK_TABLE(table), 6);
    gtk_container_set_border_width(GTK_CONTAINER(table), 4);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), table,
                       FALSE, FALSE, 4);


    row = 0;
    label = gtk_label_new("Number of steps:");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table), label,
                     0, 1, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    nsteps = gtk_adjustment_new(xgc->data.set.sc.nsteps, 1, 1000000, 1, 10, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(nsteps), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(spin), 0);
    gtk_table_attach(GTK_TABLE(table), spin,
                     1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row++;

    label = gtk_label_new("Verbosity level:");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table), label,
                     0, 1, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    verbose = gtk_adjustment_new(xgc->data.set.sc.verbose, 0, 4, 1, 10, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(verbose), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(spin), 0);
    gtk_table_attach(GTK_TABLE(table), spin,
                     1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row++;

    label = gtk_label_new("Number of CPU cores (-1 for all):");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(GTK_TABLE(table), label,
                     0, 1, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    nthreads = gtk_adjustment_new(xgc->data.set.sc.nthreads, -1, 1000, 1, 10, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(nthreads), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(spin), 0);
    gtk_table_attach(GTK_TABLE(table), spin,
                     1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row++;

    gpu = gtk_check_button_new_with_mnemonic("_Use GPU No.:");
    gtk_table_attach(GTK_TABLE(table), gpu,
                     0, 1, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(gpu), xgc->data.set.sc.usegpu);

    ugval = 0;
    if (xgc->data.set.sc.ugpu[1] == 1) ugval = 1;
    else if (xgc->data.set.sc.ugpu[2] == 1) ugval = 2;
    else if (xgc->data.set.sc.ugpu[3] == 1) ugval = 3;

    ugpu = gtk_adjustment_new(ugval, 0, 4, 1, 10, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(ugpu), 1, 2);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(spin), 0);
    gtk_table_attach(GTK_TABLE(table), spin,
                     1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);

    row++;


    gtk_widget_show_all(dialog);
    response = gtk_dialog_run(GTK_DIALOG(dialog));

    if (response == GTK_RESPONSE_ACCEPT) {
        xgc->data.set.sc.nsteps = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(nsteps));
        xgc->data.set.sc.verbose = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(verbose));
        xgc->data.set.sc.nthreads = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(nthreads));
        ugval = (gint)gtk_adjustment_get_value(GTK_ADJUSTMENT(ugpu));
        xgc->data.set.sc.ugpu[0] = xgc->data.set.sc.ugpu[1] = xgc->data.set.sc.ugpu[2] = xgc->data.set.sc.ugpu[3] = 0;
        if (ugval == 0) xgc->data.set.sc.ugpu[0] = 1;
        else if (ugval == 1) xgc->data.set.sc.ugpu[1] = 1;
        else if (ugval == 2) xgc->data.set.sc.ugpu[2] = 1;
        else if (ugval == 3) xgc->data.set.sc.ugpu[3] = 1;

        xgc->data.set.sc.usegpu = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(gpu));

        gtk_widget_destroy(dialog);
        return TRUE;
    }
    gtk_widget_destroy(dialog);

    return FALSE;
}