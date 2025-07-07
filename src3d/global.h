
#ifndef GLOBAL_H
#define GLOBAL_H

/*#include <gtk/gtk.h>
#include "dcube.h"
#include "settings.h"
#include <libgwymodule/gwymodule.h>*/

#include <gtk/gtk.h>
#include "settings.h"
#include "pool.h"
#include "plan.h"
#include "tests.h"
#include <libgwymodule/gwymodule.h>

#define TREE_VIEW_PAR_ROOT_ROWS 8
#define TREE_VIEW_MAT_ROOT_ROWS 1
#define FILES_TO_SHOW           100
#define IMAGES_TO_SHOW          100
#define TREE_MAXENT             20000  // maximum number of entities to put into tree store

enum {
    COLUMN_PARAMETER,
    COLUMN_CHECK,
    COLUMN_ID,

    COLUMN_SHOW_TOGGLE,
    NUM_COLUMNS
};


typedef struct {

    gchar *parfilename;
    gchar *matfilename;

    GPid pid;

    gboolean is_pool;
    gboolean is_psrc[1000];
    gboolean is_sf;
    gboolean is_tsf;
    gboolean is_tsff;
    gboolean is_ltsf;
    gboolean is_ltsff;
    gboolean is_outpnt[1000];
    gboolean is_outimg[1000];
    gboolean is_outpln[1000];
    /*gboolean is_outcub[1000];*/   /* volume output covers on the whole computational domain by now - visibility does not make sense */
    gboolean is_outsum[1000];
    gboolean is_outforce[1000];
    gboolean is_bx0;
    gboolean is_bxn;
    gboolean is_by0;
    gboolean is_byn;
    gboolean is_bz0;
    gboolean is_bzn;
    gboolean is_mbx0;
    gboolean is_mbxn;
    gboolean is_mby0;
    gboolean is_mbyn;
    gboolean is_mbz0;
    gboolean is_mbzn;
    gboolean is_nfff;
    gboolean is_nfff_skipi0;
    gboolean is_nfff_skipin;
    gboolean is_nfff_skipj0;
    gboolean is_nfff_skipjn;
    gboolean is_nfff_skipk0;
    gboolean is_nfff_skipkn;
    gboolean *is_nfff_point;
    gboolean is_pnfff;
    gboolean is_pnfff_skipk0;
    gboolean is_pnfff_skipkn;
    gboolean *is_pnfff_point;

    gboolean is_sphere[1000];
    gboolean is_voxel[1000];
    gboolean is_cylinder[1000];
    gboolean is_cone[1000];
    gboolean is_rcone[1000];
    gboolean *is_tetrahedron;
    gboolean is_gwydd[1000];

    /*GArray *spheres;
    GArray *voxels;
    GArray *cylinders;
    GArray *cones;
    GArray *rcones;
    GArray *tetrahedrons;
    GArray *gwydds;

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
    gchar tetgen_filebase[100][100];*/

    gchar *gsvit_location;
    gchar *gwyddion_location;

    SvSet set;
    SvSetMat set_mat;
    /*SourceFilenameBackup sfb;*/
} XGData;


typedef struct _XGControlsSetCompDomain {
    GtkWidget *vbg_outer_comp_domain;

    /* I. 1. Computational domain */
    GtkWidget *label_comp_domain_size;
    GtkObject *comp_domain_size_x;
    GtkWidget *comp_domain_size_x_spin;
    GtkObject *comp_domain_size_y;
    GtkWidget *comp_domain_size_y_spin;
    GtkObject *comp_domain_size_z;
    GtkWidget *comp_domain_size_z_spin;
    GtkWidget *label_comp_domain_spacing;
    GtkObject *comp_domain_spacing_x;
    GtkWidget *comp_domain_spacing_x_spin;
    GtkObject *comp_domain_spacing_y;
    GtkWidget *comp_domain_spacing_y_spin;
    GtkObject *comp_domain_spacing_z;
    GtkWidget *comp_domain_spacing_z_spin;
    GtkWidget *label_comp_domain_size_um;
    GtkWidget *comp_domain_size_x_um;
    GtkWidget *comp_domain_size_y_um;
    GtkWidget *comp_domain_size_z_um;
} XGControlsSetCompDomain;

typedef struct _XGControlsSetBasic {
    GtkWidget *vbg_outer_basic;

    /* II. 1. Basic parameters */
    GtkWidget *label_basic_nsteps;
    GtkObject *basic_nsteps;
    GtkWidget *basic_nsteps_spin;
    GtkWidget *label_basic_verbose;
    GtkObject *basic_verbose;
    GtkWidget *basic_verbose_spin;
    GtkWidget *label_basic_nthreads;
    GtkObject *basic_nthreads;
    GtkWidget *basic_nthreads_spin;
    GtkWidget *label_basic_ugpu;
    GtkWidget *basic_usegpu;
    GtkWidget *label_basic_ugpu_index;
    GtkObject *basic_ugpu_index;
    GtkWidget *basic_ugpu_index_spin;
} XGControlsSetBasic;

typedef struct _XGControlsSetSource {
    GtkWidget *vbg_outer_box;
    GtkWidget *vbg_outer_point_origin;
    GtkWidget *vbg_outer_source;
    GtkWidget *vbg_outer_ia;
    GtkWidget *vbg_outer_fs;
    GtkWidget *vbg_outer_ls;
    GtkWidget *vbg_outer_zmultiplier;

    /* III. A. Source box group */
    GtkWidget *label_box_start;
    GtkObject *box_i0;
    GtkWidget *box_i0_spin;
    GtkObject *box_j0;
    GtkWidget *box_j0_spin;
    GtkObject *box_k0;
    GtkWidget *box_k0_spin;
    GtkWidget *label_box_to;
    GtkObject *box_in;
    GtkWidget *box_in_spin;
    GtkObject *box_jn;
    GtkWidget *box_jn_spin;
    GtkObject *box_kn;
    GtkWidget *box_kn_spin;
    GtkWidget *label_box_size_um;
    GtkWidget *box_size_x_um;
    GtkWidget *box_size_y_um;
    GtkWidget *box_size_z_um;
    GtkWidget *label_box_boundary;
    GtkWidget *box_boundary_skipi0;
    GtkWidget *box_boundary_skipin;
    GtkWidget *box_boundary_skipj0;
    GtkWidget *box_boundary_skipjn;
    GtkWidget *box_boundary_skipk0;
    GtkWidget *box_boundary_skipkn;
    GtkWidget *box_enable_skipdepth;
    GtkWidget *label_box_skipdepth;
    GtkObject *box_skipdepth;
    GtkWidget *box_skipdepth_spin;

    /* III. 1. Point source properties */
    GtkWidget *label_point_origin_position;
    GtkObject *point_origin_position_i;
    GtkWidget *point_origin_position_i_spin;
    GtkObject *point_origin_position_j;
    GtkWidget *point_origin_position_j_spin;
    GtkObject *point_origin_position_k;
    GtkWidget *point_origin_position_k_spin;
    GtkWidget *label_point_origin_theta;
    GtkObject *point_origin_theta;
    GtkWidget *point_origin_theta_spin;
    GtkWidget *label_point_origin_phi;
    GtkObject *point_origin_phi;
    GtkWidget *point_origin_phi_spin;

    /* III. B. Source group */
    GtkWidget *label_source_mode;
    GtkWidget *source_mode;
    GtkWidget *label_source_filename;
    GtkWidget *source_filename;
    GtkWidget *source_button_change;
    GtkWidget *label_source_amplitude;
    GtkObject *source_amplitude;
    GtkWidget *source_amplitude_spin;
    GtkWidget *label_source_wavelength;
    GtkObject *source_wavelength;
    GtkWidget *source_wavelength_spin;
    GtkWidget *label_source_pulsewidth;
    GtkObject *source_pulsewidth;
    GtkWidget *source_pulsewidth_spin;

    /* III. C. Incident angle group */
    GtkWidget *label_incident_angle;
    GtkWidget *label_ia_theta;
    GtkObject *ia_theta;
    GtkWidget *ia_theta_spin;
    GtkWidget *label_ia_phi;
    GtkObject *ia_phi;
    GtkWidget *ia_phi_spin;
    GtkWidget *label_ia_psi;
    GtkObject *ia_psi;
    GtkWidget *ia_psi_spin;

    /* III. D. Focused source group */
    GtkWidget *label_fs_thetamax;
    GtkObject *fs_thetamax_deg;
    GtkWidget *fs_thetamax_deg_spin;
    GtkWidget *label_fs_fdist;
    GtkObject *fs_fdist;
    GtkWidget *fs_fdist_spin;
    GtkWidget *label_fs_polarisation;
    GtkObject *fs_pol_deg;
    GtkWidget *fs_pol_deg_spin;
    GtkWidget *label_fs_nip;
    GtkObject *fs_nip;
    GtkWidget *fs_nip_spin;
    GtkWidget *label_fs_mip;
    GtkObject *fs_mip;
    GtkWidget *fs_mip_spin;

    /* III. E. Layered source group */
    GtkWidget *hseparator_ls;
    GtkWidget *label_ls_prop;
    //GtkWidget *ls_prop;               /* treeview */
    GtkWidget *ls_prop_wrapper;         /* treeview wrapper */
    GtkWidget *layered_button_add;
    gint       layered_nlayers;

    /* III. F. Z multiplier group */
    GtkWidget *gaussian_mult_enable;
    GtkWidget *label_gaussian_mult_center;
    GtkObject *gaussian_mult_center_pos_i;
    GtkWidget *gaussian_mult_center_pos_i_spin;
    GtkObject *gaussian_mult_center_pos_j;
    GtkWidget *gaussian_mult_center_pos_j_spin;
    GtkWidget *label_gaussian_mult_radius;
    GtkObject *gaussian_mult_radius_i;
    GtkWidget *gaussian_mult_radius_i_spin;
    GtkObject *gaussian_mult_radius_j;
    GtkWidget *gaussian_mult_radius_j_spin;

    GtkWidget *radial_mult_enable;
    GtkWidget *label_radial_mult_center;
    GtkObject *radial_mult_center_pos_i;
    GtkWidget *radial_mult_center_pos_i_spin;
    GtkObject *radial_mult_center_pos_j;
    GtkWidget *radial_mult_center_pos_j_spin;
    GtkWidget *label_radial_mult_radius;
    GtkObject *radial_mult_radius_i;
    GtkWidget *radial_mult_radius_i_spin;
    GtkObject *radial_mult_radius_j;
    GtkWidget *radial_mult_radius_j_spin;

    GtkWidget *fiber_mult_enable;
    GtkWidget *label_fiber_mult_center;
    GtkObject *fiber_mult_center_pos_i;
    GtkWidget *fiber_mult_center_pos_i_spin;
    GtkObject *fiber_mult_center_pos_j;
    GtkWidget *fiber_mult_center_pos_j_spin;
    GtkWidget *label_fiber_mult_radius;
    GtkObject *fiber_mult_radius;
    GtkWidget *fiber_mult_radius_spin;
    GtkWidget *label_fiber_mult_cutoff;
    GtkObject *fiber_mult_cutoff;
    GtkWidget *fiber_mult_cutoff_spin;
    GtkWidget *label_fiber_mult_epsilon_core;
    GtkObject *fiber_mult_epsilon_core;
    GtkWidget *fiber_mult_epsilon_core_spin;
    GtkWidget *label_fiber_mult_epsilon_cladding;
    GtkObject *fiber_mult_epsilon_cladding;
    GtkWidget *fiber_mult_epsilon_cladding_spin;

} XGControlsSetSource;

typedef struct _XGControlsSetBoundary {
    /* IV. Boundary conditions */

    GtkWidget *vbg_outer_cvboundary;
    GtkWidget *vbg_outer_periodicboundary;

    /* IV. 1. Computational volume boundaries */

    GtkWidget *label_btype[6];
    GtkWidget *btype[6];

    GtkWidget *label_cvCPMLparams;
    GtkObject *bdepth[6];
    GtkWidget *bdepthspin[6];
    GtkObject *bm[6];
    GtkWidget *bmspin[6];
    GtkObject *bsigma[6];
    GtkWidget *bsigmaspin[6];
    GtkObject *bkappa[6];
    GtkWidget *bkappaspin[6];
    GtkObject *ba[6];
    GtkWidget *baspin[6];

    GtkWidget *mb[6];
    GtkObject *mbpos[6];
    GtkWidget *mbposspin[6];
} XGControlsSetBoundary;

typedef struct _XGControlsSetMedia {
    /* V. Media */

    /* V. 1. Material properties */

    GtkWidget *vbg_outer_material_prop;

    GtkWidget *label_material_voxel_filename;
    GtkWidget *material_voxel_filename;
    GtkWidget *material_voxel_button_change;

    GtkWidget *label_material_vector_filename;
    GtkWidget *material_vector_filename;
    GtkWidget *material_vector_button_change;

    GtkWidget *label_material_modecheck;
    GtkWidget *material_modecheck;

    GtkWidget *label_material_smoothsteps;
    GtkObject *material_smoothsteps;
    GtkWidget *material_smoothsteps_spin;

    /* V. 2. Add growth modifier */

    GtkWidget *vbg_outer_material_grow;

    GtkWidget *label_grow_span_from;
    GtkObject *grow_i0;
    GtkWidget *grow_i0_spin;
    GtkObject *grow_j0;
    GtkWidget *grow_j0_spin;
    GtkObject *grow_k0;
    GtkWidget *grow_k0_spin;
    GtkWidget *label_grow_span_to;
    GtkObject *grow_in;
    GtkWidget *grow_in_spin;
    GtkObject *grow_jn;
    GtkWidget *grow_jn_spin;
    GtkObject *grow_kn;
    GtkWidget *grow_kn_spin;
    GtkWidget *grow_skipi0;
    GtkWidget *grow_skipj0;
    GtkWidget *grow_skipk0;
    GtkWidget *grow_skipin;
    GtkWidget *grow_skipjn;
    GtkWidget *grow_skipkn;
    GtkWidget *label_grow_addindex;
    GtkObject *grow_addindex;
    GtkWidget *grow_addindex_spin;
    GtkWidget *label_grow_attachindex;
    GtkObject *grow_attachindex;
    GtkWidget *grow_attachindex_spin;
    GtkWidget *label_grow_subsampling;
    GtkObject *grow_subsampling;
    GtkWidget *grow_subsampling_spin;
    GtkWidget *label_grow_seed;
    GtkObject *grow_seed;
    GtkWidget *grow_seed_spin;
    GtkWidget *label_grow_nsteps;
    GtkObject *grow_nsteps;
    GtkWidget *grow_nsteps_spin;
    GtkWidget *label_grow_mobility;
    GtkObject *grow_mobility;
    GtkWidget *grow_mobility_spin;
    GtkWidget *label_grow_probability;
    GtkObject *grow_probability;
    GtkWidget *grow_probability_spin;

    /* V. 3. Add roughness modifier */

    GtkWidget *vbg_outer_material_rough;

    GtkWidget *label_rough_matindex;
    GtkObject *rough_matindex;
    GtkWidget *rough_matindex_spin;
    GtkWidget *label_rough_voidindex;
    GtkObject *rough_voidindex;
    GtkWidget *rough_voidindex_spin;
    GtkWidget *label_rough_iterations;
    GtkObject *rough_iterations;
    GtkWidget *rough_iterations_spin;
    GtkWidget *label_rough_seed;
    GtkObject *rough_seed;
    GtkWidget *rough_seed_spin;
    GtkWidget *label_rough_radiuspeak;
    GtkObject *rough_radiuspeak;
    GtkWidget *rough_radiuspeak_spin;
    GtkWidget *label_rough_radiusspan;
    GtkObject *rough_radiusspan;
    GtkWidget *rough_radiusspan_spin;
    GtkWidget *label_rough_probability;
    GtkObject *rough_probability;
    GtkWidget *rough_probability_spin;

    /* V. 4. Add spectral modifier */

    GtkWidget *vbg_outer_material_spectral;

    GtkWidget *label_spectral_sigma;
    GtkObject *spectral_sigma;
    GtkWidget *spectral_sigma_spin;
    GtkWidget *label_spectral_t;
    GtkObject *spectral_t;
    GtkWidget *spectral_t_spin;
    GtkWidget *label_spectral_matindex;
    GtkObject *spectral_matindex;
    GtkWidget *spectral_matindex_spin;
    GtkWidget *label_spectral_seed;
    GtkObject *spectral_seed;
    GtkWidget *spectral_seed_spin;

    /* V. 3. Add expression modifier */

    GtkWidget *vbg_outer_material_expr;

    GtkWidget *label_expr_span_from;
    GtkObject *expr_i0;
    GtkWidget *expr_i0_spin;
    GtkObject *expr_j0;
    GtkWidget *expr_j0_spin;
    GtkObject *expr_k0;
    GtkWidget *expr_k0_spin;
    GtkWidget *label_expr_span_to;
    GtkObject *expr_in;
    GtkWidget *expr_in_spin;
    GtkObject *expr_jn;
    GtkWidget *expr_jn_spin;
    GtkObject *expr_kn;
    GtkWidget *expr_kn_spin;
    GtkWidget *label_expr_matindex;
    GtkObject *expr_matindex;
    GtkWidget *expr_matindex_spin;
    GtkWidget *label_expr_voidindex;
    GtkObject *expr_voidindex;
    GtkWidget *expr_voidindex_spin;
    GtkWidget *label_expr_maxdist;
    GtkObject *expr_maxdist;
    GtkWidget *expr_maxdist_spin;
    GtkWidget *label_expr_distmode;
    GtkObject *expr_distmode;
    GtkWidget *expr_distmode_spin;
    GtkWidget *label_expr_expr;
    GtkWidget *expr_expr;
} XGControlsSetMedia;

typedef struct _XGControlsSetOutput {
    GtkWidget *vbg_outer_general_output;
    GtkWidget *vbg_outer_point_output;
    GtkWidget *vbg_outer_image_output;
    GtkWidget *vbg_outer_plane_output;
    GtkWidget *vbg_outer_volume_output;
    GtkWidget *vbg_outer_sum_output;
    GtkWidget *vbg_outer_force_output;
    GtkWidget *vbg_outer_box_output;

    /* VI. 1. General output */
    GtkWidget *label_general_output_filename;
    GtkWidget *general_output_filename;
    GtkWidget *general_output_filename_button_change;

    /* VI. 2. Point output properties */
    GtkWidget *label_point_output_origin_position;
    GtkObject *point_output_origin_position_i;
    GtkWidget *point_output_origin_position_i_spin;
    GtkObject *point_output_origin_position_j;
    GtkWidget *point_output_origin_position_j_spin;
    GtkObject *point_output_origin_position_k;
    GtkWidget *point_output_origin_position_k_spin;
    GtkWidget *label_point_output_step;
    GtkObject *point_output_step;
    GtkWidget *point_output_step_spin;
    GtkWidget *label_point_output_component;
    GtkWidget *point_output_component;
    GtkWidget *label_point_output_filename;
    GtkWidget *point_output_filename;
    GtkWidget *point_output_filename_button_change;

    /* VI. 3. Image output properties */
    GtkWidget *label_image_output_plane;
    GtkWidget *image_output_plane;
    GtkWidget *label_image_output_pos;
    GtkObject *image_output_pos;
    GtkWidget *image_output_pos_spin;
    GtkWidget *label_image_output_step;
    GtkObject *image_output_step;
    GtkWidget *image_output_step_spin;
    GtkWidget *label_image_output_component;
    GtkWidget *image_output_component;
    GtkWidget *label_image_output_label;
    GtkWidget *image_output_label;

    /* VI. 4. Plane output properties */
    GtkWidget *label_plane_output_plane;
    GtkWidget *plane_output_plane;
    GtkWidget *label_plane_output_pos;
    GtkObject *plane_output_pos;
    GtkWidget *plane_output_pos_spin;
    GtkWidget *label_plane_output_step;
    GtkObject *plane_output_step;
    GtkWidget *plane_output_step_spin;
    GtkWidget *label_plane_output_component;
    GtkWidget *plane_output_component;
    GtkWidget *label_plane_output_filename;
    GtkWidget *plane_output_filename;
    GtkWidget *plane_output_filename_button_change;
    GtkWidget *label_plane_output_format;
    GtkWidget *plane_output_format;
    GtkWidget *label_plane_output_start;
    GtkObject *plane_output_start;
    GtkWidget *plane_output_start_spin;
    GtkWidget *label_plane_output_stop;
    GtkObject *plane_output_stop;
    GtkWidget *plane_output_stop_spin;

    /* VI. 5. Volume output properties */
    GtkWidget *label_volume_output_step;
    GtkObject *volume_output_step;
    GtkWidget *volume_output_step_spin;
    GtkWidget *label_volume_output_component;
    GtkWidget *volume_output_component;
    GtkWidget *label_volume_output_filename;
    GtkWidget *volume_output_filename;
    GtkWidget *volume_output_filename_button_change;
    GtkWidget *label_volume_output_format;
    GtkWidget *volume_output_format;
    GtkWidget *label_volume_output_start;
    GtkObject *volume_output_start;
    GtkWidget *volume_output_start_spin;
    GtkWidget *label_volume_output_stop;
    GtkObject *volume_output_stop;
    GtkWidget *volume_output_stop_spin;

    /* VI. A. Box output group */
    GtkWidget *label_box_output_start;
    GtkObject *box_output_i0;
    GtkWidget *box_output_i0_spin;
    GtkObject *box_output_j0;
    GtkWidget *box_output_j0_spin;
    GtkObject *box_output_k0;
    GtkWidget *box_output_k0_spin;
    GtkWidget *label_box_output_to;
    GtkObject *box_output_in;
    GtkWidget *box_output_in_spin;
    GtkObject *box_output_jn;
    GtkWidget *box_output_jn_spin;
    GtkObject *box_output_kn;
    GtkWidget *box_output_kn_spin;
    GtkWidget *label_box_output_size_um;
    GtkWidget *box_output_size_x_um;
    GtkWidget *box_output_size_y_um;
    GtkWidget *box_output_size_z_um;

    /* VI. 5. Sum output properties */
    GtkWidget *label_sum_output_epsilon;
    GtkObject *sum_output_epsilon;
    GtkWidget *sum_output_epsilon_spin;
    GtkWidget *label_sum_output_mu;
    GtkObject *sum_output_mu;
    GtkWidget *sum_output_mu_spin;
    GtkWidget *label_sum_output_sigma;
    GtkObject *sum_output_sigma;
    GtkWidget *sum_output_sigma_spin;
    GtkWidget *label_sum_output_sigast;
    GtkObject *sum_output_sigast;
    GtkWidget *sum_output_sigast_spin;
    GtkWidget *label_sum_output_component;
    GtkWidget *sum_output_component;
    GtkWidget *label_sum_output_step;
    GtkObject *sum_output_step;
    GtkWidget *sum_output_step_spin;
    GtkWidget *label_sum_output_filename;
    GtkWidget *sum_output_filename;
    GtkWidget *sum_output_filename_button_change;

    /* VI. 6. Force output properties */
    GtkWidget *label_force_output_step;
    GtkObject *force_output_step;
    GtkWidget *force_output_step_spin;
    GtkWidget *label_force_output_filename;
    GtkWidget *force_output_filename;
    GtkWidget *force_output_filename_button_change;
} XGControlsSetOutput;

typedef struct _XGControlsSetNFFF {
    /* VII. NFFF */

    /* VII. 1. Near field to far field transform */

    GtkWidget *vbg_outer_nfff;

    GtkWidget *label_nfff_box_from;
    GtkObject *nfff_box_i0;
    GtkWidget *nfff_box_i0_spin;
    GtkObject *nfff_box_j0;
    GtkWidget *nfff_box_j0_spin;
    GtkObject *nfff_box_k0;
    GtkWidget *nfff_box_k0_spin;

    GtkWidget *label_nfff_box_to;
    GtkObject *nfff_box_in;
    GtkWidget *nfff_box_in_spin;
    GtkObject *nfff_box_jn;
    GtkWidget *nfff_box_jn_spin;
    GtkObject *nfff_box_kn;
    GtkWidget *nfff_box_kn_spin;

    GtkWidget *nfff_box_boundary_skipi0;
    GtkObject *nfff_skipi0_jmin;
    GtkWidget *nfff_skipi0_jmin_spin;
    GtkObject *nfff_skipi0_kmin;
    GtkWidget *nfff_skipi0_kmin_spin;
    GtkObject *nfff_skipi0_jmax;
    GtkWidget *nfff_skipi0_jmax_spin;
    GtkObject *nfff_skipi0_kmax;
    GtkWidget *nfff_skipi0_kmax_spin;

    GtkWidget *nfff_box_boundary_skipin;
    GtkObject *nfff_skipin_jmin;
    GtkWidget *nfff_skipin_jmin_spin;
    GtkObject *nfff_skipin_kmin;
    GtkWidget *nfff_skipin_kmin_spin;
    GtkObject *nfff_skipin_jmax;
    GtkWidget *nfff_skipin_jmax_spin;
    GtkObject *nfff_skipin_kmax;
    GtkWidget *nfff_skipin_kmax_spin;

    GtkWidget *nfff_box_boundary_skipj0;
    GtkObject *nfff_skipj0_imin;
    GtkWidget *nfff_skipj0_imin_spin;
    GtkObject *nfff_skipj0_kmin;
    GtkWidget *nfff_skipj0_kmin_spin;
    GtkObject *nfff_skipj0_imax;
    GtkWidget *nfff_skipj0_imax_spin;
    GtkObject *nfff_skipj0_kmax;
    GtkWidget *nfff_skipj0_kmax_spin;

    GtkWidget *nfff_box_boundary_skipjn;
    GtkObject *nfff_skipjn_imin;
    GtkWidget *nfff_skipjn_imin_spin;
    GtkObject *nfff_skipjn_kmin;
    GtkWidget *nfff_skipjn_kmin_spin;
    GtkObject *nfff_skipjn_imax;
    GtkWidget *nfff_skipjn_imax_spin;
    GtkObject *nfff_skipjn_kmax;
    GtkWidget *nfff_skipjn_kmax_spin;

    GtkWidget *nfff_box_boundary_skipk0;
    GtkObject *nfff_skipk0_imin;
    GtkWidget *nfff_skipk0_imin_spin;
    GtkObject *nfff_skipk0_jmin;
    GtkWidget *nfff_skipk0_jmin_spin;
    GtkObject *nfff_skipk0_imax;
    GtkWidget *nfff_skipk0_imax_spin;
    GtkObject *nfff_skipk0_jmax;
    GtkWidget *nfff_skipk0_jmax_spin;

    GtkWidget *nfff_box_boundary_skipkn;
    GtkObject *nfff_skipkn_imin;
    GtkWidget *nfff_skipkn_imin_spin;
    GtkObject *nfff_skipkn_jmin;
    GtkWidget *nfff_skipkn_jmin_spin;
    GtkObject *nfff_skipkn_imax;
    GtkWidget *nfff_skipkn_imax_spin;
    GtkObject *nfff_skipkn_jmax;
    GtkWidget *nfff_skipkn_jmax_spin;

    /* VII. 2. Near field to far field point */
    /* VIII. 2. Periodic near field to far field point */

    GtkWidget *vbg_outer_nfffp;

    GtkWidget *label_nfffp_pos;
    GtkObject *nfffp_pos_i;
    GtkWidget *nfffp_pos_i_spin;
    GtkObject *nfffp_pos_j;
    GtkWidget *nfffp_pos_j_spin;
    GtkObject *nfffp_pos_k;
    GtkWidget *nfffp_pos_k_spin;

    GtkWidget *label_nfffp_filename;
    GtkWidget *nfffp_filename;
    GtkWidget *nfffp_filename_button_change;

    /* VII. 3. Near field to far field area */
    /* VIII. 3. Periodic near field to far field area */

    GtkWidget *vbg_outer_nfffa;

    GtkWidget *label_nfffa_thetares;
    GtkObject *nfffa_thetares;
    GtkWidget *nfffa_thetares_spin;

    GtkWidget *label_nfffa_phires;
    GtkObject *nfffa_phires;
    GtkWidget *nfffa_phires_spin;

    GtkWidget *label_nfffa_radius;
    GtkObject *nfffa_radius;
    GtkWidget *nfffa_radius_spin;

    GtkWidget *label_nfffa_thetafrom;
    GtkObject *nfffa_thetafrom;
    GtkWidget *nfffa_thetafrom_spin;

    GtkWidget *label_nfffa_thetato;
    GtkObject *nfffa_thetato;
    GtkWidget *nfffa_thetato_spin;

    GtkWidget *label_nfffa_phifrom;
    GtkObject *nfffa_phifrom;
    GtkWidget *nfffa_phifrom_spin;

    GtkWidget *label_nfffa_phito;
    GtkObject *nfffa_phito;
    GtkWidget *nfffa_phito_spin;

    GtkWidget *nfffa_savefile;

    /* VIII. PNFFF */

    /* VIII. 1. Periodic near field to far field transform */

    GtkWidget *vbg_outer_pnfff;

    GtkWidget *label_pnfff_box_from;
    GtkObject *pnfff_box_i0;
    GtkWidget *pnfff_box_i0_spin;
    GtkObject *pnfff_box_j0;
    GtkWidget *pnfff_box_j0_spin;
    GtkObject *pnfff_box_k0;
    GtkWidget *pnfff_box_k0_spin;

    GtkWidget *label_pnfff_box_to;
    GtkObject *pnfff_box_in;
    GtkWidget *pnfff_box_in_spin;
    GtkObject *pnfff_box_jn;
    GtkWidget *pnfff_box_jn_spin;
    GtkObject *pnfff_box_kn;
    GtkWidget *pnfff_box_kn_spin;

    GtkWidget *label_nfff_box_size_um;
    GtkWidget *nfff_box_size_x_um;
    GtkWidget *nfff_box_size_y_um;
    GtkWidget *nfff_box_size_z_um;

    GtkWidget *label_pnfff_integrationx;
    GtkObject *pnfff_integration_xmin;
    GtkWidget *pnfff_integration_xmin_spin;
    GtkObject *pnfff_integration_xmax;
    GtkWidget *pnfff_integration_xmax_spin;

    GtkWidget *label_pnfff_integrationy;
    GtkObject *pnfff_integration_ymin;
    GtkWidget *pnfff_integration_ymin_spin;
    GtkObject *pnfff_integration_ymax;
    GtkWidget *pnfff_integration_ymax_spin;

    GtkWidget *pnfff_box_boundary_skipk0;
    GtkObject *pnfff_skipk0_imin;
    GtkWidget *pnfff_skipk0_imin_spin;
    GtkObject *pnfff_skipk0_jmin;
    GtkWidget *pnfff_skipk0_jmin_spin;
    GtkObject *pnfff_skipk0_imax;
    GtkWidget *pnfff_skipk0_imax_spin;
    GtkObject *pnfff_skipk0_jmax;
    GtkWidget *pnfff_skipk0_jmax_spin;

    GtkWidget *label_pnfff_box_size_um;
    GtkWidget *pnfff_box_size_x_um;
    GtkWidget *pnfff_box_size_y_um;
    GtkWidget *pnfff_box_size_z_um;

    GtkWidget *pnfff_box_boundary_skipkn;
    GtkObject *pnfff_skipkn_imin;
    GtkWidget *pnfff_skipkn_imin_spin;
    GtkObject *pnfff_skipkn_jmin;
    GtkWidget *pnfff_skipkn_jmin_spin;
    GtkObject *pnfff_skipkn_imax;
    GtkWidget *pnfff_skipkn_imax_spin;
    GtkObject *pnfff_skipkn_jmax;
    GtkWidget *pnfff_skipkn_jmax_spin;

} XGControlsSetNFFF;

typedef struct _XGControlsSet {
    GtkWidget               *vbox_page;
    GtkWidget               *sw_values;     /* scrolled window - parameter values */
    GtkWidget               *label_name;

    XGControlsSetCompDomain cc;         /* computational domain controls */
    XGControlsSetBasic      bac;        /* basic controls */
    XGControlsSetSource     sc;         /* source controls */
    XGControlsSetBoundary   boc;        /* boundary controls */
    XGControlsSetMedia      mc;         /* media controls */
    XGControlsSetOutput     oc;         /* output controls */
    XGControlsSetNFFF       fc;         /* field controls */
} XGControlsSet;

typedef struct _XGMatControlsSetObject {
    GtkWidget *vbg_outer_sphere;
    GtkWidget *vbg_outer_box;    
    GtkWidget *vbg_outer_cylinder;
    GtkWidget *vbg_outer_cone;
    GtkWidget *vbg_outer_rcone;
    GtkWidget *vbg_outer_tetrahedron;

    GtkWidget *vbg_outer_prop;

    /* I. 1. Sphere properties */
    GtkWidget *label_sphere_center;
    GtkObject *sphere_center_i;
    GtkWidget *sphere_center_i_spin;
    GtkObject *sphere_center_j;
    GtkWidget *sphere_center_j_spin;
    GtkObject *sphere_center_k;
    GtkWidget *sphere_center_k_spin;
    GtkWidget *label_sphere_radius;
    GtkObject *sphere_radius;
    GtkWidget *sphere_radius_spin;

    /* I. 2. Box properties */
    GtkWidget *label_box_from;
    GtkObject *box_i0;
    GtkWidget *box_i0_spin;
    GtkObject *box_j0;
    GtkWidget *box_j0_spin;
    GtkObject *box_k0;
    GtkWidget *box_k0_spin;
    GtkWidget *label_box_to;
    GtkObject *box_in;
    GtkWidget *box_in_spin;
    GtkObject *box_jn;
    GtkWidget *box_jn_spin;
    GtkObject *box_kn;
    GtkWidget *box_kn_spin;

    /* I. 2. Cylinder properties */
    GtkWidget *label_cyl_from;
    GtkObject *cyl_i0;
    GtkWidget *cyl_i0_spin;
    GtkObject *cyl_j0;
    GtkWidget *cyl_j0_spin;
    GtkObject *cyl_k0;
    GtkWidget *cyl_k0_spin;
    GtkWidget *label_cyl_to;
    GtkObject *cyl_in;
    GtkWidget *cyl_in_spin;
    GtkObject *cyl_jn;
    GtkWidget *cyl_jn_spin;
    GtkObject *cyl_kn;
    GtkWidget *cyl_kn_spin;
    GtkWidget *label_cyl_radius;
    GtkObject *cyl_radius;
    GtkWidget *cyl_radius_spin;

    /* I. 3. Cone properties */
    GtkWidget *label_cone_from;
    GtkObject *cone_i0;
    GtkWidget *cone_i0_spin;
    GtkObject *cone_j0;
    GtkWidget *cone_j0_spin;
    GtkObject *cone_k0;
    GtkWidget *cone_k0_spin;
    GtkWidget *label_cone_to;
    GtkObject *cone_in;
    GtkWidget *cone_in_spin;
    GtkObject *cone_jn;
    GtkWidget *cone_jn_spin;
    GtkObject *cone_kn;
    GtkWidget *cone_kn_spin;
    GtkWidget *label_cone_radius;
    GtkObject *cone_radius;
    GtkWidget *cone_radius_spin;

    /* I. 4. Cut cone properties */
    GtkWidget *label_rcone_from;
    GtkObject *rcone_i0;
    GtkWidget *rcone_i0_spin;
    GtkObject *rcone_j0;
    GtkWidget *rcone_j0_spin;
    GtkObject *rcone_k0;
    GtkWidget *rcone_k0_spin;
    GtkWidget *label_rcone_to;
    GtkObject *rcone_in;
    GtkWidget *rcone_in_spin;
    GtkObject *rcone_jn;
    GtkWidget *rcone_jn_spin;
    GtkObject *rcone_kn;
    GtkWidget *rcone_kn_spin;
    GtkWidget *label_rcone_radius1;
    GtkObject *rcone_radius1;
    GtkWidget *rcone_radius1_spin;
    GtkWidget *label_rcone_radius2;
    GtkObject *rcone_radius2;
    GtkWidget *rcone_radius2_spin;

    /* I. 5. Tetrahedron properties */
    GtkWidget *label_tthn_0;
    GtkObject *tthn_i0;
    GtkWidget *tthn_i0_spin;
    GtkObject *tthn_j0;
    GtkWidget *tthn_j0_spin;
    GtkObject *tthn_k0;
    GtkWidget *tthn_k0_spin;
    GtkWidget *label_tthn_1;
    GtkObject *tthn_i1;
    GtkWidget *tthn_i1_spin;
    GtkObject *tthn_j1;
    GtkWidget *tthn_j1_spin;
    GtkObject *tthn_k1;
    GtkWidget *tthn_k1_spin;
    GtkWidget *label_tthn_2;
    GtkObject *tthn_i2;
    GtkWidget *tthn_i2_spin;
    GtkObject *tthn_j2;
    GtkWidget *tthn_j2_spin;
    GtkObject *tthn_k2;
    GtkWidget *tthn_k2_spin;
    GtkWidget *label_tthn_3;
    GtkObject *tthn_i3;
    GtkWidget *tthn_i3_spin;
    GtkObject *tthn_j3;
    GtkWidget *tthn_j3_spin;
    GtkObject *tthn_k3;
    GtkWidget *tthn_k3_spin;

    /* I. A. Material properties group */
    GtkWidget *mm0;
    GtkWidget *mm1;
    GtkWidget *mm2;
    GtkWidget *mm4;
    GtkWidget *mm10;
    GtkWidget *mm99;

    GtkObject *epsilon;
    GtkWidget *epsilon_spin;
    GtkObject *mu;
    GtkWidget *mu_spin;
    GtkObject *sigma;
    GtkWidget *sigma_spin;    
    GtkObject *sigast;
    GtkWidget *sigast_spin;
        
    GtkObject *tepsilon;
    GtkWidget *tepsilon_spin;
    GtkObject *tmu;
    GtkWidget *tmu_spin;
    GtkObject *tsigma;
    GtkWidget *tsigma_spin;    
    GtkObject *tsigast;
    GtkWidget *tsigast_spin;

    GtkWidget *label_eps_omega_nu;
    GtkObject *cepsilon;
    GtkWidget *cepsilon_spin;
    GtkObject *omegap;
    GtkWidget *omegap_spin;
    GtkObject *nu;
    GtkWidget *nu_spin;

    GtkWidget *label_peak[2];
    GtkObject *a[2];
    GtkWidget *a_spin[2];
    GtkObject *ia_phi[2];
    GtkWidget *ia_phi_spin[2];
    GtkObject *omega[2];
    GtkWidget *omega_spin[2];
    GtkObject *gamma[2];
    GtkWidget *gamma_spin[2];

    GtkWidget *metal_model;
    GtkWidget *overriden;
} XGMatControlsSetObject;

typedef struct _XGMatControlsSet {
    GtkWidget               *vbox_page;
    GtkWidget               *sw_values; /* scrolled window - object values */
    GtkWidget               *label_name;

    XGMatControlsSetObject  oc;        /* object controls */
} XGMatControlsSet;

typedef struct {
    GtkWidget   *view_scene;

    GtkWidget   *toplevel;
    GtkWidget   *new_menu_item;
    GtkWidget   *open_menu_item;
    GtkWidget   *save_menu_item;
    GtkWidget   *saveas_menu_item;
    GtkWidget   *savematas_menu_item;
    GtkWidget   *recentfiles_menu_item;

    GtkWidget   *run_menu_item;
    GtkWidget   *stop_menu_item;

    GtkWidget   *gwydd_menu_item;
    GtkToolItem *check_tool_item;
    GtkToolItem *run_tool_item;
    GtkToolItem *stop_tool_item;
    GtkToolItem *gwydd_tool_item;

    GtkTextBuffer *tb_par;
    GtkTextBuffer *tb_mat;
    GtkTextBuffer *tb_out;
    GtkWidget     *tv_out;
    GtkWidget     *graph_combo;

    GtkWidget   *param_notebook;
    GtkWidget   *mat_notebook;

    GtkWidget   *statusbar;
    gint        statusbar_context_id;

    gboolean    par_file_success;
    gboolean    par_tree_success;

    gboolean    mat_file_success;
    gboolean    mat_tree_success;

    GwyGraphModel      *gmodel1;
    GwyGraphCurveModel *gcmodel1;
    GtkWidget          *graph_check;

    GtkWidget       *page_par_tree;
    GtkWidget       *page_mat_tree;

    GtkWidget       *tv_par;
    GtkTreeStore    *ts_par;
    gboolean        tv_par_expanded[TREE_VIEW_PAR_ROOT_ROWS];
    gchar           *tv_par_path_select_string;

    GtkWidget       *tv_mat;
    GtkTreeStore    *ts_mat;
    gboolean        tv_mat_expanded[TREE_VIEW_MAT_ROOT_ROWS];
    gchar           *tv_mat_path_select_string;

    gdouble         page_hscroll_val;
    gdouble         page_vscroll_val;

    gdouble         page_mat_hscroll_val;
    gdouble         page_mat_vscroll_val;

    gint        par_selid;
    GtkTreeIter par_seliter;
    gboolean    par_sel_removable;

    gint        mat_selid;
    GtkTreeIter mat_seliter;
    gboolean    mat_sel_removable;


    GtkWidget   *image_check;
    GtkWidget   *image_squarecheck;
    GtkWidget   *image_combo;
    gint        image_selected;

    XGControlsSet cs;
    XGMatControlsSet cs_mat;

    gint        nimagestoshow;
    gchar       *imagestoshow[100];

    gchar       *curdir;

    GLfloat     xrot;
    GLfloat     yrot;
    gfloat      xdiff;
    gfloat      ydiff;

    gint        undo;
    gint        logtimer;

    GLfloat     zoomfactor;
    GLfloat     zoomfactor_prev;
    GLdouble    viewposx;
    GLdouble    viewposy;
    GLdouble    viewposx_start;
    GLdouble    viewposy_start;
    GLdouble    viewposx_prev;
    GLdouble    viewposy_prev;
    //   guint     timeout;

    GKeyFile    *keyfile;
    gchar       *config_path;

    gchar       *tmpfilename;
    gchar       *tmpmatfilename;

    gchar       *filestoshow[FILES_TO_SHOW];
    gint        formattoshow[IMAGES_TO_SHOW];
    gint        nfilestoshow;

    gint         timeout_func;
    gint         watch;
    GFileMonitor *outputfile_monitor;
    GwyDataField *outputfield[100];
    GwyContainer *container;

    gchar        *loadfile;    

    gboolean    usekiller;

    gboolean    par_lock_controls_update;
    gboolean    mat_lock_controls_update;

    XGData      data;
} XGControls;


void reset_par_elemets_visibility(XGControls *xgc);

void set_main_window_title(XGControls *xgc);

gboolean ggwy_data_field_inside(GwyDataField *data_field, gint i, gint j);

gchar* get_spectra_path(const gchar* string);

gchar *g_canonicalize_filename (const gchar *filename, const gchar *relative_to);

#endif  /* GLOBAL_H */
