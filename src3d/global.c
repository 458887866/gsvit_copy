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


/* global.c :
*  global functions
*/

#include "global.h"
#include <gtk/gtk.h>
#include <string.h>

void reset_par_elemets_visibility(XGControls *xgc)
{
    gint i;
    xgc->data.is_pool = TRUE;
    xgc->data.is_sf = TRUE;
    xgc->data.is_tsf = TRUE;
    xgc->data.is_tsff = TRUE;
    xgc->data.is_ltsf = TRUE;
    xgc->data.is_ltsff = TRUE;
    xgc->data.is_bx0 = xgc->data.is_bxn = xgc->data.is_by0 = xgc->data.is_byn = xgc->data.is_bz0 = xgc->data.is_bzn = FALSE;
    xgc->data.is_mbx0 = xgc->data.is_mbxn = xgc->data.is_mby0 = xgc->data.is_mbyn = xgc->data.is_mbz0 = xgc->data.is_mbzn = TRUE;

    for (i = 0; i < 1000; i++) {
        xgc->data.is_psrc[i] = TRUE;
        xgc->data.is_outpnt[i] = TRUE;
        xgc->data.is_outimg[i] = TRUE;
        xgc->data.is_outpln[i] = TRUE;
        /* xgc->data.is_outcub[i] = TRUE; */    /* volume output covers on the whole computational domain by now - visibility does not make sense */
        xgc->data.is_outsum[i] = TRUE;
        xgc->data.is_outforce[i] = TRUE;
    }

    xgc->data.is_nfff = TRUE;
    xgc->data.is_nfff_skipi0 = xgc->data.is_nfff_skipin = xgc->data.is_nfff_skipj0 = xgc->data.is_nfff_skipjn = xgc->data.is_nfff_skipk0 = xgc->data.is_nfff_skipkn = TRUE;

    xgc->data.is_pnfff = TRUE;
    xgc->data.is_pnfff_skipk0 = xgc->data.is_pnfff_skipkn = TRUE;

    for (i = 0; i < 100000; i++) {
        xgc->data.is_nfff_point[i] = TRUE;
        xgc->data.is_pnfff_point[i] = TRUE;
    }

    for (i = 0; i < 1000; i++) {
        xgc->data.is_sphere[i] = TRUE;
        xgc->data.is_voxel[i] = TRUE;
        xgc->data.is_cylinder[i] = TRUE;
        xgc->data.is_cone[i] = TRUE;
        xgc->data.is_rcone[i] = TRUE;
        xgc->data.is_gwydd[i] = TRUE;
    }
    for (i = 0; i < 1000000; i++) {
        xgc->data.is_tetrahedron[i] = TRUE;
    }
}

void set_main_window_title(XGControls *xgc)
{
    gchar title[256] = {0};

    if (NULL == xgc->data.parfilename)
        g_snprintf(title, sizeof(title), "Untitled - XSvit");
    else
        g_snprintf(title, sizeof(title), "%s - XSvit", xgc->data.parfilename);

    gtk_window_set_title(GTK_WINDOW(xgc->toplevel), title);
}

gboolean ggwy_data_field_inside(GwyDataField *data_field, gint i, gint j)
{
    return (i >= 0 && j >= 0 && i < data_field->xres && j < data_field->yres);
}

gchar* get_spectra_path(const gchar* string)
{
    gchar *path;
#ifdef G_OS_WIN32
    path = g_build_path("\\", get_self_dir(), "share", "gsvit", "data", "spectra", string, NULL);
#endif

#ifndef G_OS_WIN32
    path = g_build_path("/", g_path_get_dirname(g_path_get_dirname(get_self_dir())), "share", "gsvit", "data", "spectra", string, NULL);
#endif

    return path;
}

/**
* g_canonicalize_filename:
* @filename: (type filename): the name of the file
* @relative_to: (type filename) (nullable): the relative directory, or %NULL
* to use the current working directory
*
* Gets the canonical file name from @filename. All triple slashes are turned into
* single slashes, and all `..` and `.`s resolved against @relative_to.
*
* Symlinks are not followed, and the returned path is guaranteed to be absolute.
*
* If @filename is an absolute path, @relative_to is ignored. Otherwise,
* @relative_to will be prepended to @filename to make it absolute. @relative_to
* must be an absolute path, or %NULL. If @relative_to is %NULL, it'll fallback
* to g_get_current_dir().
*
* This function never fails, and will canonicalize file paths even if they don't
* exist.
*
* No file system I/O is done.
*
* Returns: (type filename) (transfer full): a newly allocated string with the
* canonical file path
* Since: 2.58
*/
gchar *
g_canonicalize_filename (const gchar *filename, const gchar *relative_to)
{
    gchar *canon, *start, *p, *q;
    guint i;

    g_return_val_if_fail (relative_to == NULL || g_path_is_absolute (relative_to), NULL);

    if (!g_path_is_absolute (filename))
    {
        gchar *cwd_allocated = NULL;
        const gchar  *cwd;

        if (relative_to != NULL)
            cwd = relative_to;
        else
            cwd = cwd_allocated = g_get_current_dir ();

        canon = g_build_filename (cwd, filename, NULL);
        g_free (cwd_allocated);
    }
    else
    {
        canon = g_strdup (filename);
    }

    start = (char *)g_path_skip_root (canon);

    if (start == NULL)
    {
        /* This shouldn't really happen, as g_get_current_dir() should
        return an absolute pathname, but bug 573843 shows this is
        not always happening */
        g_free (canon);
        return g_build_filename (G_DIR_SEPARATOR_S, filename, NULL);
    }

    /* POSIX allows double slashes at the start to
    * mean something special (as does windows too).
    * So, "//" != "/", but more than two slashes
    * is treated as "/".
    */
    i = 0;
    for (p = start - 1;
        (p >= canon) &&
        G_IS_DIR_SEPARATOR (*p);
        p--)
        i++;
    if (i > 2)
    {
        i -= 1;
        start -= i;
        memmove (start, start+i, strlen (start+i) + 1);
    }

    /* Make sure we're using the canonical dir separator */
    p++;
    while (p < start && G_IS_DIR_SEPARATOR (*p))
        *p++ = G_DIR_SEPARATOR;

    p = start;
    while (*p != 0)
    {
        if (p[0] == '.' && (p[1] == 0 || G_IS_DIR_SEPARATOR (p[1])))
        {
            memmove (p, p+1, strlen (p+1)+1);
        }
        else if (p[0] == '.' && p[1] == '.' && (p[2] == 0 || G_IS_DIR_SEPARATOR (p[2])))
        {
            q = p + 2;
            /* Skip previous separator */
            p = p - 2;
            if (p < start)
                p = start;
            while (p > start && !G_IS_DIR_SEPARATOR (*p))
                p--;
            if (G_IS_DIR_SEPARATOR (*p))
                *p++ = G_DIR_SEPARATOR;
            memmove (p, q, strlen (q)+1);
        }
        else
        {
            /* Skip until next separator */
            while (*p != 0 && !G_IS_DIR_SEPARATOR (*p))
                p++;

            if (*p != 0)
            {
                /* Canonicalize one separator */
                *p++ = G_DIR_SEPARATOR;
            }
        }

        /* Remove additional separators */
        q = p;
        while (*q && G_IS_DIR_SEPARATOR (*q))
            q++;

        if (p != q)
            memmove (p, q, strlen (q) + 1);
    }

    /* Remove trailing slashes */
    if (p > start && G_IS_DIR_SEPARATOR (*(p-1)))
        *(p-1) = 0;

    return canon;
}
