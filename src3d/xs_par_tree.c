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


/* par_tree.c :
*  Parameter file group tree functions
*/

#include "xs_par_tree.h"
#include "write.h"
#include "messages.h"
#include "string.h"

void
pt_assemble_row_text(XGControls    *xgc,
                     gchar         *buff,          /* output text */
                     gint          buff_size,      /* output text buffer size */
                     gint          node_id,        /* node group id (expandable node) */
                     gint          leaf_pos)       /* position in node */
{
    gint index = 0;

    if (node_id == SET_POOL) {
        if (leaf_pos == 0)
            g_snprintf(buff, buff_size, "Size [vx]: %d x %d x %d", xgc->data.set.sp.xres, xgc->data.set.sp.yres, xgc->data.set.sp.zres);
        else if (leaf_pos == 1)
            g_snprintf(buff, buff_size, "Spacing [μm]: %g x %g x %g", xgc->data.set.sp.dx*1e6, xgc->data.set.sp.dy*1e6, xgc->data.set.sp.dz*1e6);
    } else if (node_id == SET_BASIC) {
        if (leaf_pos == 0)
            g_snprintf(buff, buff_size, "Number of steps: %d", xgc->data.set.sc.nsteps);
        else if (leaf_pos == 1)
            g_snprintf(buff, buff_size, "Verbose level: %d", xgc->data.set.sc.verbose);
        else if (leaf_pos == 2)
            g_snprintf(buff, buff_size, "Number of threads: %d", xgc->data.set.sc.nthreads);
        else if (leaf_pos == 3) {
            if (xgc->data.set.sc.usegpu)
                g_snprintf(buff, buff_size, "Use GPU: yes");
            else
                g_snprintf(buff, buff_size, "Use GPU: no");
        } else if (leaf_pos == 4) {
            if (xgc->data.set.sc.ugpu[0] == 1)
                g_snprintf(buff, buff_size, "GPU No.: 0");
            else if (xgc->data.set.sc.ugpu[1] == 1)
                g_snprintf(buff, buff_size, "GPU No.: 1");
            else if (xgc->data.set.sc.ugpu[2] == 1)
                g_snprintf(buff, buff_size, "GPU No.: 2");
            else if (xgc->data.set.sc.ugpu[3] == 1)
                g_snprintf(buff, buff_size, "GPU No.: 3");
            else
                g_snprintf(buff, buff_size, "GPU No.: none");
        }
    } else if (node_id >= SET_PSOURCE && node_id < SET_POUT) {
        if (xgc->data.set.ss.pnts[leaf_pos].source_mode == 0)
            g_snprintf(buff, buff_size, "Point source: %d %d %d (%s)", xgc->data.set.ss.pnts[leaf_pos].point_origin_position_i, xgc->data.set.ss.pnts[leaf_pos].point_origin_position_j, xgc->data.set.ss.pnts[leaf_pos].point_origin_position_k, xgc->data.set.ss.pnts[leaf_pos].source_filename);
        else
            g_snprintf(buff, buff_size, "Point source: %d %d %d (mode %d)", xgc->data.set.ss.pnts[leaf_pos].point_origin_position_i, xgc->data.set.ss.pnts[leaf_pos].point_origin_position_j, xgc->data.set.ss.pnts[leaf_pos].point_origin_position_k, xgc->data.set.ss.pnts[leaf_pos].source_mode);
    } else if (node_id == SET_SF) {
        g_snprintf(buff, buff_size, "SF: angles %g %g %g deg", xgc->data.set.ss.sf.ia_theta * 180 / G_PI, xgc->data.set.ss.sf.ia_phi * 180 / G_PI, xgc->data.set.ss.sf.ia_psi * 180 / G_PI);
    } else if (node_id == SET_TSF) {
        g_snprintf(buff, buff_size, "TSF: %d %d %d ... %d %d %d, angles %g %g %g deg", xgc->data.set.ss.tsf.box_i0, xgc->data.set.ss.tsf.box_j0, xgc->data.set.ss.tsf.box_k0,
                   xgc->data.set.ss.tsf.box_in, xgc->data.set.ss.tsf.box_jn, xgc->data.set.ss.tsf.box_kn,
                   xgc->data.set.ss.tsf.ia_theta * 180 / G_PI, xgc->data.set.ss.tsf.ia_phi * 180 / G_PI, xgc->data.set.ss.tsf.ia_psi * 180 / G_PI);
    } else if (node_id == SET_TSFF) {
        g_snprintf(buff, buff_size, "TSFF: %d %d %d ... %d %d %d, tmax %g deg, %d x %d waves", xgc->data.set.ss.tsff.box_i0, xgc->data.set.ss.tsff.box_j0, xgc->data.set.ss.tsff.box_k0,
                   xgc->data.set.ss.tsff.box_in, xgc->data.set.ss.tsff.box_jn, xgc->data.set.ss.tsff.box_kn,
                   xgc->data.set.ss.tsff.focused_thetamax * 180 / G_PI, xgc->data.set.ss.tsff.focused_nip, xgc->data.set.ss.tsff.focused_mip);
    } else if (node_id == SET_LTSF) {
        g_snprintf(buff, buff_size, "LTSF: %d %d %d ... %d %d %d, angles %g %g %g deg", xgc->data.set.ss.ltsf.box_i0, xgc->data.set.ss.ltsf.box_j0, xgc->data.set.ss.ltsf.box_k0,
                   xgc->data.set.ss.ltsf.box_in, xgc->data.set.ss.ltsf.box_jn, xgc->data.set.ss.ltsf.box_kn,
                   xgc->data.set.ss.ltsf.ia_theta * 180 / G_PI, xgc->data.set.ss.ltsf.ia_phi * 180 / G_PI, xgc->data.set.ss.ltsf.ia_psi * 180 / G_PI);
    } else if (node_id == SET_LTSFF) {
        g_snprintf(buff, buff_size, "LTSFF: %d %d %d ... %d %d %d, tmax %g deg, %d x %d waves", xgc->data.set.ss.ltsff.box_i0, xgc->data.set.ss.ltsff.box_j0, xgc->data.set.ss.ltsff.box_k0,
                   xgc->data.set.ss.ltsff.box_in, xgc->data.set.ss.ltsff.box_jn, xgc->data.set.ss.ltsff.box_kn,
                   xgc->data.set.ss.ltsff.focused_thetamax * 180 / G_PI, xgc->data.set.ss.ltsff.focused_nip, xgc->data.set.ss.ltsff.focused_mip);
    } else if (node_id == SET_BND) {
        /* IV. 1. Boundary conditions */
        if (leaf_pos == 0) {
            if (xgc->data.set.sb.bx0 == SV_BOUNDARY_NONE)
                g_snprintf(buff, buff_size, "x0: none");
            else if (xgc->data.set.sb.bx0 == SV_BOUNDARY_PEC)
                g_snprintf(buff, buff_size, "x0: PEC");
            else if (xgc->data.set.sb.bx0 == SV_BOUNDARY_LIAO)
                g_snprintf(buff, buff_size, "x0: Liao");
            else if (xgc->data.set.sb.bx0 == SV_BOUNDARY_CPML)
                g_snprintf(buff, buff_size, "x0: CPML (depth=%d, m=%d, sigma=%g, a=%g, kappa=%g)", xgc->data.set.sb.depth_bx0, xgc->data.set.sb.m_bx0, xgc->data.set.sb.sigma_bx0, xgc->data.set.sb.a_bx0, xgc->data.set.sb.kappa_bx0);
            else if (xgc->data.set.sb.bx0 == SV_BOUNDARY_PERIODIC)
                g_snprintf(buff, buff_size, "x0: periodic");
        } else if (leaf_pos == 1) {
            if (xgc->data.set.sb.bxn == SV_BOUNDARY_NONE)
                g_snprintf(buff, buff_size, "xn: none");
            else if (xgc->data.set.sb.bxn == SV_BOUNDARY_PEC)
                g_snprintf(buff, buff_size, "xn: PEC");
            else if (xgc->data.set.sb.bxn == SV_BOUNDARY_LIAO)
                g_snprintf(buff, buff_size, "xn: Liao");
            else if (xgc->data.set.sb.bxn == SV_BOUNDARY_CPML)
                g_snprintf(buff, buff_size, "xn: CPML (depth=%d, m=%d, sigma=%g, a=%g, kappa=%g)", xgc->data.set.sb.depth_bxn, xgc->data.set.sb.m_bxn, xgc->data.set.sb.sigma_bxn, xgc->data.set.sb.a_bxn, xgc->data.set.sb.kappa_bxn);
            else if (xgc->data.set.sb.bxn == SV_BOUNDARY_PERIODIC)
                g_snprintf(buff, buff_size, "xn: periodic");
        } else if (leaf_pos == 2) {
            if (xgc->data.set.sb.by0 == SV_BOUNDARY_NONE)
                g_snprintf(buff, buff_size, "y0: none");
            else if (xgc->data.set.sb.by0 == SV_BOUNDARY_PEC)
                g_snprintf(buff, buff_size, "y0: PEC");
            else if (xgc->data.set.sb.by0 == SV_BOUNDARY_LIAO)
                g_snprintf(buff, buff_size, "y0: Liao");
            else if (xgc->data.set.sb.by0 == SV_BOUNDARY_CPML)
                g_snprintf(buff, buff_size, "y0: CPML (depth=%d, m=%d, sigma=%g, a=%g, kappa=%g)", xgc->data.set.sb.depth_by0, xgc->data.set.sb.m_by0, xgc->data.set.sb.sigma_by0, xgc->data.set.sb.a_by0, xgc->data.set.sb.kappa_by0);
            else if (xgc->data.set.sb.by0 == SV_BOUNDARY_PERIODIC)
                g_snprintf(buff, buff_size, "y0: periodic");
        } else if (leaf_pos == 3) {
            if (xgc->data.set.sb.byn == SV_BOUNDARY_NONE)
                g_snprintf(buff, buff_size, "yn: none");
            else if (xgc->data.set.sb.byn == SV_BOUNDARY_PEC)
                g_snprintf(buff, buff_size, "yn: PEC");
            else if (xgc->data.set.sb.byn == SV_BOUNDARY_LIAO)
                g_snprintf(buff, buff_size, "yn: Liao");
            else if (xgc->data.set.sb.byn == SV_BOUNDARY_CPML)
                g_snprintf(buff, buff_size, "yn: CPML (depth=%d, m=%d, sigma=%g, a=%g, kappa=%g)", xgc->data.set.sb.depth_byn, xgc->data.set.sb.m_byn, xgc->data.set.sb.sigma_byn, xgc->data.set.sb.a_byn, xgc->data.set.sb.kappa_byn);
            else if (xgc->data.set.sb.byn == SV_BOUNDARY_PERIODIC)
                g_snprintf(buff, buff_size, "yn: periodic");
        } else if (leaf_pos == 4) {
            if (xgc->data.set.sb.bz0 == SV_BOUNDARY_NONE)
                g_snprintf(buff, buff_size, "z0: none");
            else if (xgc->data.set.sb.bz0 == SV_BOUNDARY_PEC)
                g_snprintf(buff, buff_size, "z0: PEC");
            else if (xgc->data.set.sb.bz0 == SV_BOUNDARY_LIAO)
                g_snprintf(buff, buff_size, "z0: Liao");
            else if (xgc->data.set.sb.bz0 == SV_BOUNDARY_CPML)
                g_snprintf(buff, buff_size, "z0: CPML (depth=%d, m=%d, sigma=%g, a=%g, kappa=%g)", xgc->data.set.sb.depth_bz0, xgc->data.set.sb.m_bz0, xgc->data.set.sb.sigma_bz0, xgc->data.set.sb.a_bz0, xgc->data.set.sb.kappa_bz0);
            else if (xgc->data.set.sb.bz0 == SV_BOUNDARY_PERIODIC)
                g_snprintf(buff, buff_size, "z0: periodic");
        } else if (leaf_pos == 5) {
            if (xgc->data.set.sb.bzn == SV_BOUNDARY_NONE)
                g_snprintf(buff, buff_size, "zn: none");
            else if (xgc->data.set.sb.bzn == SV_BOUNDARY_PEC)
                g_snprintf(buff, buff_size, "zn: PEC");
            else if (xgc->data.set.sb.bzn == SV_BOUNDARY_LIAO)
                g_snprintf(buff, buff_size, "zn: Liao");
            else if (xgc->data.set.sb.bzn == SV_BOUNDARY_CPML)
                g_snprintf(buff, buff_size, "zn: CPML (depth=%d, m=%d, sigma=%g, a=%g, kappa=%g)", xgc->data.set.sb.depth_bzn, xgc->data.set.sb.m_bzn, xgc->data.set.sb.sigma_bzn, xgc->data.set.sb.a_bzn, xgc->data.set.sb.kappa_bzn);
            else if (xgc->data.set.sb.bzn == SV_BOUNDARY_PERIODIC)
                g_snprintf(buff, buff_size, "zn: periodic");
        } else if (leaf_pos == 6) {     /* optional periodic boundaries */
            g_snprintf(buff, buff_size, "x0 periodic boundary at %d", xgc->data.set.smb.bx0pos);
        } else if (leaf_pos == 7) {
            g_snprintf(buff, buff_size, "xn periodic boundary at %d", xgc->data.set.smb.bxnpos);
        } else if (leaf_pos == 8) {
            g_snprintf(buff, buff_size, "y0 periodic boundary at %d", xgc->data.set.smb.by0pos);
        } else if (leaf_pos == 9) {
            g_snprintf(buff, buff_size, "yn periodic boundary at %d", xgc->data.set.smb.bynpos);
        } else if (leaf_pos == 10) {
            g_snprintf(buff, buff_size, "z0 periodic boundary at %d", xgc->data.set.smb.bz0pos);
        } else if (leaf_pos == 11) {
            g_snprintf(buff, buff_size, "zn periodic boundary at %d", xgc->data.set.smb.bznpos);
        }
    } else if (node_id == SET_MEDIUM) {
        /* V. 1. Media */

        /* V. 1. Material properties */
        if (leaf_pos == 0) {
            if (xgc->data.set.sm.in_voxel_filename)
                g_snprintf(buff, buff_size, "Voxel-by-voxel data: %s", xgc->data.set.sm.in_voxel_filename);
            else
                g_snprintf(buff, buff_size, "Voxel-by-voxel data: undefined");
        } else if (leaf_pos == 1) {
            if (xgc->data.set.sm.in_vector_filename)
                g_snprintf(buff, buff_size, "Vector data: %s", xgc->data.set.sm.in_vector_filename);
            else
                g_snprintf(buff, buff_size, "Vector data: undefined");
        } else if (leaf_pos == 2) {
            if (xgc->data.set.sm.matmode_check == 1)
                g_snprintf(buff, buff_size, "Material mode checking: yes");
            else
                g_snprintf(buff, buff_size, "Material mode checking: no");
        } else if (leaf_pos == 3) {
            g_snprintf(buff, buff_size, "Material smoothing iterations: %d", xgc->data.set.sm.smooth);
        }
    } else if (node_id >= SET_GROW && node_id < SET_ROUGHNESS) {
        /* V. 2. Add growth modifier */

        index = leaf_pos - NOL_MEDIUM_MAT_PROPS;
        g_snprintf(buff, buff_size, "Grow material %d on  material %d", xgc->data.set.sm.grow_addindex[index], xgc->data.set.sm.grow_attachindex[index]);
    } else if (node_id >= SET_ROUGHNESS && node_id < SET_SPECTRAL) {
        /* V. 3. Add roughness modifier */

        index = leaf_pos - NOL_MEDIUM_MAT_PROPS - xgc->data.set.sm.ngrowths;
        g_snprintf(buff, buff_size, "Modify roughness %d via random Gaussians", xgc->data.set.sm.rough_matindex[index]);
    } else if (node_id >= SET_SPECTRAL && node_id < SET_EXPRESSION) {
        /* V. 4. Add spectral modifier */

        index = leaf_pos - NOL_MEDIUM_MAT_PROPS - xgc->data.set.sm.ngrowths - xgc->data.set.sm.nroughens;
        g_snprintf(buff, buff_size, "Modify roughness %d via spectral synthesis", xgc->data.set.sm.spectral_matindex[index]);
    } else if (node_id >= SET_EXPRESSION && node_id < SET_NFAREA) {
        /* V. 5. Add expression modifier */

        index = leaf_pos - NOL_MEDIUM_MAT_PROPS - xgc->data.set.sm.ngrowths - xgc->data.set.sm.nroughens - xgc->data.set.sm.nspectrals;
        g_snprintf(buff, buff_size, "Modify roughness %d via analytical expression", xgc->data.set.sm.expr_matindex[index]);
    } else if (node_id == SET_OUT) {
        g_snprintf(buff, buff_size, "General output file: %s", xgc->data.set.so.outfile);
    } else if (node_id >= SET_POUT && node_id < SET_IOUT) {
        /* VI. 2. Point output properties */

        index = leaf_pos - NOL_OUTPUT_PROPS;
        g_snprintf(buff, buff_size, "Point output (%s): %d %d %d every %d to %s",
                   cpstring(xgc->data.set.so.pnts[index].component),
                   xgc->data.set.so.pnts[index].i, xgc->data.set.so.pnts[index].j, xgc->data.set.so.pnts[index].k, xgc->data.set.so.pnts[index].step, xgc->data.set.so.pnts[index].filebase);
    } else if (node_id >= SET_IIOUT && node_id < SET_COUT) {
        /* VI. 3. Image output properties */

        index = leaf_pos - NOL_OUTPUT_PROPS - xgc->data.set.so.npnts;
        g_snprintf(buff, buff_size, "Image output (%s): %s every %d (%s)",
                   cpstring(xgc->data.set.so.imgs[index].component),
                   planestring(xgc->data.set.so.imgs[index].i, xgc->data.set.so.imgs[index].j, xgc->data.set.so.imgs[index].k),
                   xgc->data.set.so.imgs[index].step, xgc->data.set.so.imgs[index].filebase);
    } else if (node_id >= SET_IOUT && node_id < SET_IIOUT) {
        /* VI. 4. Plane output properties */

        index = leaf_pos - NOL_OUTPUT_PROPS - xgc->data.set.so.npnts - xgc->data.set.so.nimgs;
        g_snprintf(buff, buff_size, "Plane output (%s): %s every %d",
                   cpstring(xgc->data.set.so.plns[index].component),
                   planestring(xgc->data.set.so.plns[index].i,
                               xgc->data.set.so.plns[index].j, xgc->data.set.so.plns[index].k), xgc->data.set.so.plns[index].step);
    } else if (node_id >= SET_COUT && node_id < SET_SOUT) {
        /* VI. 5. Volume output properties */

        index = leaf_pos - NOL_OUTPUT_PROPS - xgc->data.set.so.npnts - xgc->data.set.so.nimgs - xgc->data.set.so.nplns;
        if (xgc->data.set.so.cubs[index].format)
            g_snprintf(buff, buff_size, "Volume output ASCII (%s) every %d",
                       vcpstring(xgc->data.set.so.cubs[index].component),
                       xgc->data.set.so.cubs[index].step);
        else
            g_snprintf(buff, buff_size, "Volume output binary (%s) every %d",
                       vcpstring(xgc->data.set.so.cubs[index].component),
                       xgc->data.set.so.cubs[index].step);
    } else if (node_id >= SET_SOUT && node_id < SET_FOUT) {
        /* VI. 5. Sum output properties */

        index = leaf_pos - NOL_OUTPUT_PROPS - xgc->data.set.so.npnts - xgc->data.set.so.nimgs - xgc->data.set.so.nplns - xgc->data.set.so.ncubs;
        if (!xgc->data.set.so.sums[index].stringbased) {
            if (xgc->data.set.so.sums[index].layered_epsilon == -1)
                g_snprintf(buff, buff_size, "Sum (%s) every %d for any material",
                           scpstring(xgc->data.set.so.sums[index].component), xgc->data.set.so.sums[index].step);
            else
                g_snprintf(buff, buff_size, "Sum (%s) every %d for mat (%g %g %g %g)",
                           scpstring(xgc->data.set.so.sums[index].component),
                           xgc->data.set.so.sums[index].step, xgc->data.set.so.sums[index].layered_epsilon, xgc->data.set.so.sums[index].layered_mu, xgc->data.set.so.sums[index].layered_sigma, xgc->data.set.so.sums[index].layered_sigast);
        } else
            g_snprintf(buff, buff_size, "Sum (%s) every %d for %s",
                       scpstring(xgc->data.set.so.sums[index].component),
                       xgc->data.set.so.sums[index].step, xgc->data.set.so.sums[index].string);
    } else if (node_id >= SET_FOUT && node_id < SET_GROW) {
        /* VI. 6. Force output properties */

        index = leaf_pos - NOL_OUTPUT_PROPS - xgc->data.set.so.npnts - xgc->data.set.so.nimgs - xgc->data.set.so.nplns - xgc->data.set.so.ncubs - xgc->data.set.so.nsums;
        g_snprintf(buff, buff_size, "Force every %d in box (%d %d %d ... %d %d %d) to %s", xgc->data.set.so.forces[index].step,
                   xgc->data.set.so.forces[index].box_i0, xgc->data.set.so.forces[index].box_j0, xgc->data.set.so.forces[index].box_k0, xgc->data.set.so.forces[index].box_in, xgc->data.set.so.forces[index].box_jn, xgc->data.set.so.forces[index].box_kn,
                   xgc->data.set.so.forces[index].filename);
    } else if (node_id == SET_NFFF) {
        /* VII. 1. Near field to far field transform box */

        if (leaf_pos == 0)
            g_snprintf(buff, buff_size, "NFFF box: %d %d %d ... %d %d %d", xgc->data.set.sf.box_i0, xgc->data.set.sf.box_j0, xgc->data.set.sf.box_k0, xgc->data.set.sf.box_in, xgc->data.set.sf.box_jn, xgc->data.set.sf.box_kn);
        else if (leaf_pos == 1) {
            if (xgc->data.set.sf.box_boundary_skipi0)
                g_snprintf(buff, buff_size, "NFFF skip i0: %d %d ... %d %d", xgc->data.set.sf.skipi0_jmin, xgc->data.set.sf.skipi0_kmin, xgc->data.set.sf.skipi0_jmax, xgc->data.set.sf.skipi0_kmax);
            else
                g_snprintf(buff, buff_size, "NFFF skip i0: none");
        } else if (leaf_pos == 2) {
            if (xgc->data.set.sf.box_boundary_skipin)
                g_snprintf(buff, buff_size, "NFFF skip in: %d %d ... %d %d", xgc->data.set.sf.skipin_jmin, xgc->data.set.sf.skipin_kmin, xgc->data.set.sf.skipin_jmax, xgc->data.set.sf.skipin_kmax);
            else
                g_snprintf(buff, buff_size, "NFFF skip in: none");
        } else if (leaf_pos == 3) {
            if (xgc->data.set.sf.box_boundary_skipj0)
                g_snprintf(buff, buff_size, "NFFF skip j0: %d %d ... %d %d", xgc->data.set.sf.skipj0_imin, xgc->data.set.sf.skipj0_kmin, xgc->data.set.sf.skipj0_imax, xgc->data.set.sf.skipj0_kmax);
            else
                g_snprintf(buff, buff_size, "NFFF skip j0: none");
        } else if (leaf_pos == 4) {
            if (xgc->data.set.sf.box_boundary_skipjn)
                g_snprintf(buff, buff_size, "NFFF skip jn: %d %d ... %d %d", xgc->data.set.sf.skipjn_imin, xgc->data.set.sf.skipjn_kmin, xgc->data.set.sf.skipjn_imax, xgc->data.set.sf.skipjn_kmax);
            else
                g_snprintf(buff, buff_size, "NFFF skip jn: none");
        } else if (leaf_pos == 5) {
            if (xgc->data.set.sf.box_boundary_skipk0)
                g_snprintf(buff, buff_size, "NFFF skip k0: %d %d ... %d %d", xgc->data.set.sf.skipk0_imin, xgc->data.set.sf.skipk0_jmin, xgc->data.set.sf.skipk0_imax, xgc->data.set.sf.skipk0_jmax);
            else
                g_snprintf(buff, buff_size, "NFFF skip k0: none");
        } else if (leaf_pos == 6) {
            if (xgc->data.set.sf.box_boundary_skipkn)
                g_snprintf(buff, buff_size, "NFFF skip kn: %d %d ... %d %d", xgc->data.set.sf.skipkn_imin, xgc->data.set.sf.skipkn_jmin, xgc->data.set.sf.skipkn_imax, xgc->data.set.sf.skipkn_jmax);
            else
                g_snprintf(buff, buff_size, "NFFF skip kn: none");
        }
    } else if (node_id >= SET_NFFFP && node_id < SET_PNFAREA) {
        /* VII. 2. Near field to far field point */

        index = leaf_pos - NOL_NFFF_PROPS;
        g_snprintf(buff, buff_size, "Far field point at (%d %d %d) to %s", xgc->data.set.sf.ri[index], xgc->data.set.sf.rj[index], xgc->data.set.sf.rk[index], xgc->data.set.sf.source_filename[index]);
    } else if (node_id >= SET_NFAREA && node_id < SET_NFFFP) {
        /* VII. 3. Near field to far field area */

        index = leaf_pos - NOL_NFFF_PROPS - xgc->data.set.sf.nrs;
        g_snprintf(buff, buff_size, "Far field area (%dx%d)", xgc->data.set.sf.area_thetares[index], xgc->data.set.sf.area_phires[index]);
    } else if (node_id == SET_PNFFF) {
        /* VIII. 1. Periodic near field to far field transform box */

        if (leaf_pos == 0)
            g_snprintf(buff, buff_size, "Periodic NFFF box: %d %d %d ... %d %d %d", xgc->data.set.spf.box_i0, xgc->data.set.spf.box_j0, xgc->data.set.spf.box_k0, xgc->data.set.spf.box_in, xgc->data.set.spf.box_jn, xgc->data.set.spf.box_kn);
        else if (leaf_pos == 1) {
            if (xgc->data.set.spf.box_boundary_skipk0)
                g_snprintf(buff, buff_size, "Periodic NFFF skip k0: %d %d ... %d %d", xgc->data.set.spf.skipk0_imin, xgc->data.set.spf.skipk0_jmin, xgc->data.set.spf.skipk0_imax, xgc->data.set.spf.skipk0_jmax);
            else
                g_snprintf(buff, buff_size, "Periodic NFFF skip k0: none");
        } else if (leaf_pos == 2) {
            if (xgc->data.set.spf.box_boundary_skipkn)
                g_snprintf(buff, buff_size, "Periodic NFFF skip kn: %d %d ... %d %d", xgc->data.set.spf.skipkn_imin, xgc->data.set.spf.skipkn_jmin, xgc->data.set.spf.skipkn_imax, xgc->data.set.spf.skipkn_jmax);
            else
                g_snprintf(buff, buff_size, "Periodic NFFF skip kn: none");
        }
    } else if (node_id >= SET_PNFFFP) {
        /* VIII. 2. Periodic near field to far field point */

        index = leaf_pos - NOL_PNFFF_PROPS;
        g_snprintf(buff, buff_size, "Periodic far field point at (%d %d %d) to %s", xgc->data.set.spf.ri[index], xgc->data.set.spf.rj[index], xgc->data.set.spf.rk[index], xgc->data.set.spf.source_filename[index]);
    } else if (node_id >= SET_PNFAREA && node_id < SET_PNFFFP) {
        /* VIII. 3. Periodic near field to far field area */

        index = leaf_pos - NOL_PNFFF_PROPS - xgc->data.set.spf.nrs;
        g_snprintf(buff, buff_size, "Periodic far field area (%dx%d)", xgc->data.set.spf.area_thetares[index], xgc->data.set.spf.area_phires[index]);
    }
}

void pt_remeber_expanded_rows_and_selection(XGControls *xgc)
{
    /* remember expanded rows and selection */

    GtkTreePath *path_expand = NULL, *path_select = NULL;
    gint i;

    for (i = 0; i < TREE_VIEW_PAR_ROOT_ROWS; i++) {
        path_expand = gtk_tree_path_new_from_indices(i, -1);
        xgc->tv_par_expanded[i] = gtk_tree_view_row_expanded(GTK_TREE_VIEW(xgc->tv_par), path_expand);
        gtk_tree_path_free(path_expand);
    }

    if (xgc->par_selid != -1) {
        path_select = gtk_tree_model_get_path(GTK_TREE_MODEL(xgc->ts_par), &(xgc->par_seliter));
        if (NULL != xgc->tv_par_path_select_string)
            g_free(xgc->tv_par_path_select_string);
        xgc->tv_par_path_select_string = gtk_tree_path_to_string(path_select);
    }

    xgc->page_hscroll_val = gtk_adjustment_get_value(gtk_scrolled_window_get_hadjustment(GTK_SCROLLED_WINDOW(xgc->cs.sw_values)));
    xgc->page_vscroll_val = gtk_adjustment_get_value(gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(xgc->cs.sw_values)));
}

void pt_restore_expanded_rows_and_selection(XGControls *xgc)
{
    /* restore expanded rows and selection */

    gint i;
    GtkTreePath *path_expand = NULL, *path_select = NULL;
    GtkTreeSelection *selection = NULL;

    for (i = 0; i < TREE_VIEW_PAR_ROOT_ROWS; i++) {
        if (xgc->tv_par_expanded[i]) {
            path_expand = gtk_tree_path_new_from_indices(i, -1);
            gtk_tree_view_expand_row(GTK_TREE_VIEW(xgc->tv_par), path_expand, TRUE);
            gtk_tree_path_free(path_expand);
        }
    }

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(xgc->tv_par));
    if (xgc->tv_par_path_select_string != NULL) {
        path_select = gtk_tree_path_new_from_string(xgc->tv_par_path_select_string);
        if (path_select != NULL)
            gtk_tree_selection_select_path(selection, path_select);
    }

    gtk_adjustment_set_value(gtk_scrolled_window_get_hadjustment(GTK_SCROLLED_WINDOW(xgc->cs.sw_values)), xgc->page_hscroll_val);
    gtk_adjustment_set_value(gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(xgc->cs.sw_values)), xgc->page_vscroll_val);
}

void pt_create_tree(XGControls *xgc)
{
    GtkTextIter     start;
    GtkTextIter     end;
    GtkWidget       *dialog;
    gchar           buff[256];
    GtkTreeIter     iter, child;
    gint            i, active;
    gchar           *contents;
    GError          *err = NULL;
    gchar           *textbuf;

    /* write parfile and read it exactly same way as in gsvit */

    if (NULL == xgc->tmpfilename)
        xgc->tmpfilename = get_temporary_par_filename(xgc);

    if (xgc->tmpfilename) {
        pt_remeber_expanded_rows_and_selection(xgc);

        gtk_text_buffer_get_start_iter(xgc->tb_par, &start);
        gtk_text_buffer_get_end_iter(xgc->tb_par, &end);
        textbuf = gtk_text_buffer_get_text(xgc->tb_par, &start, &end, TRUE);
        if (!g_file_set_contents(xgc->tmpfilename, textbuf, -1, &err)) {
            dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                                            GTK_DIALOG_DESTROY_WITH_PARENT,
                                            GTK_MESSAGE_ERROR,
                                            GTK_BUTTONS_CLOSE,
                                            "Error writing to temporary file: %s",
                                            g_strdup(err->message));
            gtk_dialog_run(GTK_DIALOG(dialog));
            gtk_widget_destroy(dialog);
        }
        //printf("dir: %s\n", g_get_current_dir());

        /*clear settings*/
        clear_settings(&(xgc->data.set), TRUE);        
        //////////////////////////////////////////////////////////////////////////

        if (parse_settings(xgc->tmpfilename, &(xgc->data.set), FALSE)) {
            /* no errors in parfile -> create tree*/

            xgc->nfilestoshow = 0;
            memset(xgc->filestoshow, 0, FILES_TO_SHOW * sizeof(gchar*));
            xgc->nimagestoshow = 0;
            memset(xgc->imagestoshow, 0, IMAGES_TO_SHOW * sizeof(gchar*));

            for (i = 0; i < 100; i++)
                xgc->outputfield[i] = NULL;

            gtk_tree_store_clear(xgc->ts_par);

            gtk_tree_store_append(xgc->ts_par, &iter, NULL);
            gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, "Computational domain", COLUMN_CHECK, FALSE, COLUMN_ID, -1, COLUMN_SHOW_TOGGLE, FALSE, -1);

            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_POOL, 0);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_pool, COLUMN_ID, SET_POOL, COLUMN_SHOW_TOGGLE, TRUE, -1);

            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_POOL, 1);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_POOL, COLUMN_SHOW_TOGGLE, FALSE, -1);


            gtk_tree_store_append(xgc->ts_par, &iter, NULL);
            gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, "Basic parameters", COLUMN_CHECK, FALSE, COLUMN_ID, -2, COLUMN_SHOW_TOGGLE, FALSE, -1);

            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BASIC, 0);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BASIC, COLUMN_SHOW_TOGGLE, FALSE, -1);

            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BASIC, 1);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BASIC, COLUMN_SHOW_TOGGLE, FALSE, -1);

            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BASIC, 2);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BASIC, COLUMN_SHOW_TOGGLE, FALSE, -1);

            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BASIC, 3);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BASIC, COLUMN_SHOW_TOGGLE, FALSE, -1);

            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BASIC, 4);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BASIC, COLUMN_SHOW_TOGGLE, FALSE, -1);

            gtk_tree_store_append(xgc->ts_par, &iter, NULL);
            gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, "Sources", COLUMN_CHECK, FALSE, COLUMN_ID, -3, COLUMN_SHOW_TOGGLE, FALSE, -1);

            /* Point source */
            for (i = 0; i < MIN(xgc->data.set.ss.npnts, TREE_MAXENT); i++) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_PSOURCE, i);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                //get_ps_store_text(xgc, i, buff, sizeof(buff));
                //                if (xgc->data.set.ss.pnts[i].mode == 0) 
                //                    g_snprintf(buff, sizeof(buff), "Point source: %d %d %d (%s)", xgc->data.set.ss.pnts[i].i, xgc->data.set.ss.pnts[i].j, xgc->data.set.ss.pnts[i].k, xgc->data.set.ss.pnts[i].filename); 
                //                else 
                //                    g_snprintf(buff, sizeof(buff), "Point source: %d %d %d (mode %d)", xgc->data.set.ss.pnts[i].i, xgc->data.set.ss.pnts[i].j, xgc->data.set.ss.pnts[i].k, xgc->data.set.ss.pnts[i].mode);

                //gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_PSOURCE + i, COLUMN_SHOW_TOGGLE, TRUE, -1);

                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_psrc[i], COLUMN_ID, SET_PSOURCE + i, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }

            /* Source SF */
            //if (xgc->data.set.ss.sf.filename) {
            if (xgc->data.set.ss.sf.valid) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_SF, -1);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                //get_sf_store_text(xgc, buff, sizeof(buff));
                //g_snprintf(buff, sizeof(buff), "SF: angles %g %g %g deg", xgc->data.set.ss.sf.theta*180/G_PI, xgc->data.set.ss.sf.phi*180/G_PI, xgc->data.set.ss.sf.psi*180/G_PI); 
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_sf, COLUMN_ID, SET_SF, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }

            /* Source TSF */
            //if ((xgc->data.set.ss.tsf.i0 + xgc->data.set.ss.tsf.j0 + xgc->data.set.ss.tsf.k0 + xgc->data.set.ss.tsf.i1 + xgc->data.set.ss.tsf.j1 + xgc->data.set.ss.tsf.k1) != 0) {
            if (xgc->data.set.ss.tsf.valid) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_TSF, -1);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                //get_tsf_store_text(xgc, buff, sizeof(buff));
                // 	            g_snprintf(buff, sizeof(buff), "TSF: %d %d %d ... %d %d %d, angles %g %g %g deg", xgc->data.set.ss.tsf.i0, xgc->data.set.ss.tsf.j0, xgc->data.set.ss.tsf.k0,
                //                                                                                                  xgc->data.set.ss.tsf.i1, xgc->data.set.ss.tsf.j1, xgc->data.set.ss.tsf.k1, 
                //                                                                                                  xgc->data.set.ss.tsf.theta*180/G_PI, xgc->data.set.ss.tsf.phi*180/G_PI, xgc->data.set.ss.tsf.psi*180/G_PI); 

                //gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_TSF, COLUMN_SHOW_TOGGLE, TRUE, -1);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_tsf, COLUMN_ID, SET_TSF, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }

            /* Source TSFF */
            //if ((xgc->data.set.ss.tsff.i0 + xgc->data.set.ss.tsff.j0 + xgc->data.set.ss.tsff.k0 + xgc->data.set.ss.tsff.i1 + xgc->data.set.ss.tsff.j1 + xgc->data.set.ss.tsff.k1) != 0) {
            if (xgc->data.set.ss.tsff.valid) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_TSFF, -1);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                //get_tsff_store_text(xgc, buff, sizeof(buff));
                // 	            g_snprintf(buff, sizeof(buff), "TSFF: %d %d %d ... %d %d %d, tmax %g deg, %d x %d waves", xgc->data.set.ss.tsff.i0, xgc->data.set.ss.tsff.j0, xgc->data.set.ss.tsff.k0,
                //                                                                                                          xgc->data.set.ss.tsff.i1, xgc->data.set.ss.tsff.j1, xgc->data.set.ss.tsff.k1, 
                //                                                                                                          xgc->data.set.ss.tsff.thetamax*180/G_PI, xgc->data.set.ss.tsff.nint, xgc->data.set.ss.tsff.mint); 

                //gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_TSFF, COLUMN_SHOW_TOGGLE, TRUE, -1);

                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_tsff, COLUMN_ID, SET_TSFF, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }

            /* Source LTSF */
            //if ((xgc->data.set.ss.ltsf.i0 + xgc->data.set.ss.ltsf.j0 + xgc->data.set.ss.ltsf.k0 + xgc->data.set.ss.ltsf.i1 + xgc->data.set.ss.ltsf.j1 + xgc->data.set.ss.ltsf.k1) != 0) {
            if (xgc->data.set.ss.ltsf.valid) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_LTSF, -1);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                //get_ltsf_store_text(xgc, buff, sizeof(buff));
                //                g_snprintf(buff, sizeof(buff), "LTSF: %d %d %d ... %d %d %d, angles %g %g %g deg", xgc->data.set.ss.ltsf.i0, xgc->data.set.ss.ltsf.j0, xgc->data.set.ss.ltsf.k0,
                //                                                                                                   xgc->data.set.ss.ltsf.i1, xgc->data.set.ss.ltsf.j1, xgc->data.set.ss.ltsf.k1, 
                //                                                                                                   xgc->data.set.ss.ltsf.theta*180/G_PI, xgc->data.set.ss.ltsf.phi*180/G_PI, xgc->data.set.ss.ltsf.psi*180/G_PI); 

                //gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_LTSF, COLUMN_SHOW_TOGGLE, TRUE, -1);

                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_ltsf, COLUMN_ID, SET_LTSF, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }

            /* Source LTSFF */
            //if ((xgc->data.set.ss.ltsff.i0 + xgc->data.set.ss.ltsff.j0 + xgc->data.set.ss.ltsff.k0 + xgc->data.set.ss.ltsff.i1 + xgc->data.set.ss.ltsff.j1 + xgc->data.set.ss.ltsff.k1) != 0) {
            if (xgc->data.set.ss.ltsff.valid) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_LTSFF, -1);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                //get_ltsff_store_text(xgc, buff, sizeof(buff));
                //                g_snprintf(buff, sizeof(buff), "LTSFF: %d %d %d ... %d %d %d, tmax %g deg, %d x %d waves", xgc->data.set.ss.ltsff.i0, xgc->data.set.ss.ltsff.j0, xgc->data.set.ss.ltsff.k0,
                //                                                                                                           xgc->data.set.ss.ltsff.i1, xgc->data.set.ss.ltsff.j1, xgc->data.set.ss.ltsff.k1, 
                //                                                                                                           xgc->data.set.ss.ltsff.thetamax*180/G_PI, xgc->data.set.ss.ltsff.nint, xgc->data.set.ss.ltsff.mint); 

                //gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_LTSFF, COLUMN_SHOW_TOGGLE, TRUE, -1);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_ltsff, COLUMN_ID, SET_LTSFF, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }


            gtk_tree_store_append(xgc->ts_par, &iter, NULL);
            gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, "Boundary conditions", COLUMN_CHECK, FALSE, COLUMN_ID, -3, COLUMN_SHOW_TOGGLE, FALSE, -1);

            /* Boundary x0*/
            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 0);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            /*if (xgc->data.set.sb.bx0==SV_BOUNDARY_NONE)
            g_snprintf(buff, sizeof(buff), "x0: none");
            else if (xgc->data.set.sb.bx0==SV_BOUNDARY_PEC)
            g_snprintf(buff, sizeof(buff), "x0: PEC");
            else if (xgc->data.set.sb.bx0==SV_BOUNDARY_LIAO)
            g_snprintf(buff, sizeof(buff), "x0: Liao");
            else if (xgc->data.set.sb.bx0==SV_BOUNDARY_CPML)
            g_snprintf(buff, sizeof(buff), "x0: CPML (depth=%d, m=%d, sigma=%g, a=%g, kappa=%g)", xgc->data.set.sb.depth_bx0, xgc->data.set.sb.m_bx0, xgc->data.set.sb.sigma_bx0, xgc->data.set.sb.a_bx0, xgc->data.set.sb.kappa_bx0);
            else if (xgc->data.set.sb.bx0==SV_BOUNDARY_PERIODIC)
            g_snprintf(buff, sizeof(buff), "x0: periodic");*/
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_bx0, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);

            /* Boundary xn */
            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 1);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            /*if (xgc->data.set.sb.bxn==SV_BOUNDARY_NONE)
            g_snprintf(buff, sizeof(buff), "xn: none");
            else if (xgc->data.set.sb.bxn==SV_BOUNDARY_PEC)
            g_snprintf(buff, sizeof(buff), "xn: PEC");
            else if (xgc->data.set.sb.bxn==SV_BOUNDARY_LIAO)
            g_snprintf(buff, sizeof(buff), "xn: Liao");
            else if (xgc->data.set.sb.bxn==SV_BOUNDARY_CPML)
            g_snprintf(buff, sizeof(buff), "xn: CPML (depth=%d, m=%d, sigma=%g, a=%g, kappa=%g)", xgc->data.set.sb.depth_bxn, xgc->data.set.sb.m_bxn, xgc->data.set.sb.sigma_bxn, xgc->data.set.sb.a_bxn, xgc->data.set.sb.kappa_bxn);
            else if (xgc->data.set.sb.bxn==SV_BOUNDARY_PERIODIC)
            g_snprintf(buff, sizeof(buff), "xn: periodic");*/
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_bxn, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);

            /* Boundary y0 */
            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 2);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            /*if (xgc->data.set.sb.by0==SV_BOUNDARY_NONE)
            g_snprintf(buff, sizeof(buff), "y0: none");
            else if (xgc->data.set.sb.by0==SV_BOUNDARY_PEC)
            g_snprintf(buff, sizeof(buff), "y0: PEC");
            else if (xgc->data.set.sb.by0==SV_BOUNDARY_LIAO)
            g_snprintf(buff, sizeof(buff), "y0: Liao");
            else if (xgc->data.set.sb.by0==SV_BOUNDARY_CPML)
            g_snprintf(buff, sizeof(buff), "y0: CPML (depth=%d, m=%d, sigma=%g, a=%g, kappa=%g)", xgc->data.set.sb.depth_by0, xgc->data.set.sb.m_by0, xgc->data.set.sb.sigma_by0, xgc->data.set.sb.a_by0, xgc->data.set.sb.kappa_by0);
            else if (xgc->data.set.sb.by0==SV_BOUNDARY_PERIODIC)
            g_snprintf(buff, sizeof(buff), "y0: periodic");*/
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_by0, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);

            /* Boundary yn */
            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 3);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            /*if (xgc->data.set.sb.byn==SV_BOUNDARY_NONE)
            g_snprintf(buff, sizeof(buff), "yn: none");
            else if (xgc->data.set.sb.byn==SV_BOUNDARY_PEC)
            g_snprintf(buff, sizeof(buff), "yn: PEC");
            else if (xgc->data.set.sb.byn==SV_BOUNDARY_LIAO)
            g_snprintf(buff, sizeof(buff), "yn: Liao");
            else if (xgc->data.set.sb.byn==SV_BOUNDARY_CPML)
            g_snprintf(buff, sizeof(buff), "yn: CPML (depth=%d, m=%d, sigma=%g, a=%g, kappa=%g)", xgc->data.set.sb.depth_byn, xgc->data.set.sb.m_byn, xgc->data.set.sb.sigma_byn, xgc->data.set.sb.a_byn, xgc->data.set.sb.kappa_byn);
            else if (xgc->data.set.sb.byn==SV_BOUNDARY_PERIODIC)
            g_snprintf(buff, sizeof(buff), "yn: periodic");*/
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_byn, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);

            /* Boundary z0 */
            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 4);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            /*if (xgc->data.set.sb.bz0==SV_BOUNDARY_NONE)
            g_snprintf(buff, sizeof(buff), "z0: none");
            else if (xgc->data.set.sb.bz0==SV_BOUNDARY_PEC)
            g_snprintf(buff, sizeof(buff), "z0: PEC");
            else if (xgc->data.set.sb.bz0==SV_BOUNDARY_LIAO)
            g_snprintf(buff, sizeof(buff), "z0: Liao");
            else if (xgc->data.set.sb.bz0==SV_BOUNDARY_CPML)
            g_snprintf(buff, sizeof(buff), "z0: CPML (depth=%d, m=%d, sigma=%g, a=%g, kappa=%g)", xgc->data.set.sb.depth_bz0, xgc->data.set.sb.m_bz0, xgc->data.set.sb.sigma_bz0, xgc->data.set.sb.a_bz0, xgc->data.set.sb.kappa_bz0);
            else if (xgc->data.set.sb.bz0==SV_BOUNDARY_PERIODIC)
            g_snprintf(buff, sizeof(buff), "z0: periodic");*/
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_bz0, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);

            /* Boundary zn */
            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 5);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            /*if (xgc->data.set.sb.bzn==SV_BOUNDARY_NONE)
            g_snprintf(buff, sizeof(buff), "zn: none");
            else if (xgc->data.set.sb.bzn==SV_BOUNDARY_PEC)
            g_snprintf(buff, sizeof(buff), "zn: PEC");
            else if (xgc->data.set.sb.bzn==SV_BOUNDARY_LIAO)
            g_snprintf(buff, sizeof(buff), "zn: Liao");
            else if (xgc->data.set.sb.bzn==SV_BOUNDARY_CPML)
            g_snprintf(buff, sizeof(buff), "zn: CPML (depth=%d, m=%d, sigma=%g, a=%g, kappa=%g)", xgc->data.set.sb.depth_bzn, xgc->data.set.sb.m_bzn, xgc->data.set.sb.sigma_bzn, xgc->data.set.sb.a_bzn, xgc->data.set.sb.kappa_bzn);
            else if (xgc->data.set.sb.bzn==SV_BOUNDARY_PERIODIC)
            g_snprintf(buff, sizeof(buff), "zn: periodic");*/
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_bzn, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);

            /* optional periodic boundaries */
            if (xgc->data.set.smb.bx0 == SV_BOUNDARY_PERIODIC) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 6);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_mbx0, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
            if (xgc->data.set.smb.bxn == SV_BOUNDARY_PERIODIC) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 7);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_mbxn, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
            if (xgc->data.set.smb.by0 == SV_BOUNDARY_PERIODIC) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 8);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_mby0, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
            if (xgc->data.set.smb.byn == SV_BOUNDARY_PERIODIC) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 9);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_mbyn, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
            if (xgc->data.set.smb.bz0 == SV_BOUNDARY_PERIODIC) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 10);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_mbz0, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
            if (xgc->data.set.smb.bzn == SV_BOUNDARY_PERIODIC) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_BND, 11);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_mbzn, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }

            /*if (xgc->data.set.smb.bx0 == SV_BOUNDARY_PERIODIC) {
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            g_snprintf(buff, sizeof(buff), "x0 periodic boundary at %d", xgc->data.set.smb.bx0pos);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
            if (xgc->data.set.smb.bxn == SV_BOUNDARY_PERIODIC) {
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            g_snprintf(buff, sizeof(buff), "xn periodic boundary at %d", xgc->data.set.smb.bxnpos);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
            if (xgc->data.set.smb.by0 == SV_BOUNDARY_PERIODIC) {
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            g_snprintf(buff, sizeof(buff), "y0 periodic boundary at %d", xgc->data.set.smb.by0pos);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
            if (xgc->data.set.smb.byn == SV_BOUNDARY_PERIODIC) {
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            g_snprintf(buff, sizeof(buff), "yn periodic boundary at %d", xgc->data.set.smb.bynpos);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
            if (xgc->data.set.smb.bz0 == SV_BOUNDARY_PERIODIC) {
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            g_snprintf(buff, sizeof(buff), "z0 periodic boundary at %d", xgc->data.set.smb.bz0pos);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }
            if (xgc->data.set.smb.bzn == SV_BOUNDARY_PERIODIC) {
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            g_snprintf(buff, sizeof(buff), "zn periodic boundary at %d", xgc->data.set.smb.bznpos);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_BND, COLUMN_SHOW_TOGGLE, TRUE, -1);
            }*/


            gtk_tree_store_append(xgc->ts_par, &iter, NULL);
            gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, "Media", COLUMN_CHECK, FALSE, COLUMN_ID, -4, -1);

            /* Voxel data*/
            /*if (xgc->data.set.sm.in_voxel_filename)*/ {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_MEDIUM, 0);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                //g_snprintf(buff, sizeof(buff), "Voxel-by-voxel data: %s", xgc->data.set.sm.in_voxel_filename); 
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MEDIUM, -1);
            }

            /* Vector data */
            /*if (xgc->data.set.sm.in_vector_filename)*/ {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_MEDIUM, 1);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                //g_snprintf(buff, sizeof(buff), "Vector data: %s", xgc->data.set.sm.in_vector_filename); 
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MEDIUM, -1);
            }

            /* Material mode */
            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_MEDIUM, 2);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            /*if (xgc->data.set.sm.matmode_check==1)
            g_snprintf(buff, sizeof(buff), "Material mode checking: yes");
            else
            g_snprintf(buff, sizeof(buff), "Material mode checking: no"); */
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MEDIUM, -1);

            /* Smoothing iterations */
            /* g_snprintf(buff, sizeof(buff), "Medium smoothing iterations: %d", xgc->data.set.sm.smooth); */
            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_MEDIUM, 3);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_MEDIUM, -1);

            /* V. 2. Add growth modifier */
            for (i = 0; i < MIN(xgc->data.set.sm.ngrowths, TREE_MAXENT); i++) {
                /* g_snprintf(buff, sizeof(buff), "Grow medium %d on  medium %d", xgc->data.set.sm.grow_maxindex[i], xgc->data.set.sm.grow_attachindex[i]); */
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_GROW, i + NOL_MEDIUM_MAT_PROPS);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_GROW + i, COLUMN_SHOW_TOGGLE, FALSE, -1);
            }

            /*for (i=0; i<MIN(xgc->data.set.sm.ngrowths, TREE_MAXENT); i++) {
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            g_snprintf(buff, sizeof(buff), "Grow medium %d on  medium %d", xgc->data.set.sm.grow_maxindex[i], xgc->data.set.sm.grow_attachindex[i]);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_GROW+i, COLUMN_SHOW_TOGGLE, FALSE, -1);
            }*/


            /* V. 3. Add roughness modifier */
            for (i = 0; i < MIN(xgc->data.set.sm.nroughens, TREE_MAXENT); i++) {
                /* g_snprintf(buff, sizeof(buff), "Add roughness to medium %d", xgc->data.set.sm.rough_matindex[i]); */
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_ROUGHNESS, i + NOL_MEDIUM_MAT_PROPS + xgc->data.set.sm.ngrowths);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_ROUGHNESS + i, COLUMN_SHOW_TOGGLE, FALSE, -1);
            }

            /*for (i = 0; i<MIN(xgc->data.set.sm.nroughens, TREE_MAXENT); i++) {
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            g_snprintf(buff, sizeof(buff), "Add roughness to medium %d", xgc->data.set.sm.rough_matindex[i]);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_ROUGHEN+i, COLUMN_SHOW_TOGGLE, FALSE, -1);
            }*/

            /* V. 4. Add spectral modifier */
            for (i = 0; i < MIN(xgc->data.set.sm.nspectrals, TREE_MAXENT); i++) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_SPECTRAL, i + NOL_MEDIUM_MAT_PROPS + xgc->data.set.sm.ngrowths + xgc->data.set.sm.nroughens);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_SPECTRAL + i, COLUMN_SHOW_TOGGLE, FALSE, -1);
            }

            /* V. 5. Add expression modifier */
            for (i = 0; i < MIN(xgc->data.set.sm.nexprs, TREE_MAXENT); i++) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_EXPRESSION, i + NOL_MEDIUM_MAT_PROPS + xgc->data.set.sm.ngrowths + xgc->data.set.sm.nroughens + xgc->data.set.sm.nspectrals);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_EXPRESSION + i, COLUMN_SHOW_TOGGLE, FALSE, -1);
            }

            gtk_tree_store_append(xgc->ts_par, &iter, NULL);
            gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, "Outputs", COLUMN_CHECK, FALSE, COLUMN_ID, -5, COLUMN_SHOW_TOGGLE, FALSE, -1);

            /* General output file */
            pt_assemble_row_text(xgc, buff, sizeof(buff), SET_OUT, 0);
            gtk_tree_store_append(xgc->ts_par, &child, &iter);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_OUT, COLUMN_SHOW_TOGGLE, FALSE, -1);

            /*gtk_tree_store_append(xgc->ts_par, &child, &iter);
            g_snprintf(buff, sizeof(buff), "General output file: %s", xgc->data.set.so.outfile);
            gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_OUT, COLUMN_SHOW_TOGGLE, FALSE, -1);*/

            /* Point output */
            for (i = 0; i < MIN(xgc->data.set.so.npnts, TREE_MAXENT); i++) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_POUT, NOL_OUTPUT_PROPS + i);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_outpnt[i], COLUMN_ID, SET_POUT + i, COLUMN_SHOW_TOGGLE, TRUE, -1);

                /*gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "Point output (%s): %d %d %d every %d to %s",
                cpstring(xgc->data.set.so.pnts[i].component),
                xgc->data.set.so.pnts[i].i, xgc->data.set.so.pnts[i].j, xgc->data.set.so.pnts[i].k, xgc->data.set.so.pnts[i].step, xgc->data.set.so.pnts[i].filebase);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_POUT+i, COLUMN_SHOW_TOGGLE, TRUE, -1);*/

                if (xgc->nfilestoshow < 99) {
                    xgc->filestoshow[xgc->nfilestoshow] = g_strdup(xgc->data.set.so.pnts[i].filebase);
                    xgc->formattoshow[xgc->nfilestoshow] = xgc->data.set.so.pnts[i].component;
                    xgc->nfilestoshow += 1;
                }
            }

            /* Image output */
            for (i = 0; i < MIN(xgc->data.set.so.nimgs, TREE_MAXENT); i++) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_IIOUT, NOL_OUTPUT_PROPS + xgc->data.set.so.npnts + i);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_outimg[i], COLUMN_ID, SET_IIOUT + i, COLUMN_SHOW_TOGGLE, TRUE, -1);

                /*gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "Image output (%s): %s every %d (%s)",
                cpstring(xgc->data.set.so.imgs[i].component),
                planestring(xgc->data.set.so.imgs[i].i, xgc->data.set.so.imgs[i].j, xgc->data.set.so.imgs[i].k),
                xgc->data.set.so.imgs[i].step, xgc->data.set.so.imgs[i].filebase);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_IIOUT+i, COLUMN_SHOW_TOGGLE, TRUE, -1);*/

                if (xgc->nimagestoshow < 99) {
                    g_snprintf(buff, sizeof(buff), "%s component, %s, %s",
                               cpstring(xgc->data.set.so.imgs[i].component),
                               planestring(xgc->data.set.so.imgs[i].i, xgc->data.set.so.imgs[i].j, xgc->data.set.so.imgs[i].k), xgc->data.set.so.imgs[i].filebase);
                    xgc->imagestoshow[xgc->nimagestoshow] = g_strdup(buff);
                    xgc->nimagestoshow += 1;
                }
            }

            /* Plane output */
            for (i = 0; i < MIN(xgc->data.set.so.nplns, TREE_MAXENT); i++) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_IOUT, NOL_OUTPUT_PROPS + xgc->data.set.so.npnts + xgc->data.set.so.nimgs + i);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_outpln[i], COLUMN_ID, SET_IOUT + i, COLUMN_SHOW_TOGGLE, TRUE, -1);

                /*gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "Plane output (%s): %s every %d",
                cpstring(xgc->data.set.so.plns[i].component),
                planestring(xgc->data.set.so.plns[i].i,
                xgc->data.set.so.plns[i].j, xgc->data.set.so.plns[i].k), xgc->data.set.so.plns[i].step);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_IOUT+i, COLUMN_SHOW_TOGGLE, TRUE, -1);*/                
            }

            /* Volume output */
            for (i = 0; i < MIN(xgc->data.set.so.ncubs, TREE_MAXENT); i++) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_COUT, NOL_OUTPUT_PROPS + xgc->data.set.so.npnts + xgc->data.set.so.nimgs + xgc->data.set.so.nplns + i);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_COUT + i, COLUMN_SHOW_TOGGLE, FALSE, -1);   /* volume output covers on the whole computational domain by now - visibility does not make sense */

                /*gtk_tree_store_append(xgc->ts_par, &child, &iter);
                if (xgc->data.set.so.cubs[i].format)
                g_snprintf(buff, sizeof(buff), "Volume output ASCII (%s) every %d",
                vcpstring(xgc->data.set.so.cubs[i].component),
                xgc->data.set.so.cubs[i].step);
                else
                g_snprintf(buff, sizeof(buff), "Volume output binary (%s) every %d",
                vcpstring(xgc->data.set.so.cubs[i].component),
                xgc->data.set.so.cubs[i].step);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_COUT+i, COLUMN_SHOW_TOGGLE, FALSE, -1);*/
            }

            /* Sum output */
            for (i = 0; i < MIN(xgc->data.set.so.nsums, TREE_MAXENT); i++) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_SOUT, NOL_OUTPUT_PROPS + xgc->data.set.so.npnts + xgc->data.set.so.nimgs + xgc->data.set.so.nplns + xgc->data.set.so.ncubs + i);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_outsum[i], COLUMN_ID, SET_SOUT + i, COLUMN_SHOW_TOGGLE, TRUE, -1);

                /*gtk_tree_store_append(xgc->ts_par, &child, &iter);
                if (!xgc->data.set.so.sums[i].stringbased) {
                if (xgc->data.set.so.sums[i].layered_epsilon == -1)
                g_snprintf(buff, sizeof(buff), "Sum (%s) every %d for any material",
                scpstring(xgc->data.set.so.sums[i].component), xgc->data.set.so.sums[i].step);
                else
                g_snprintf(buff, sizeof(buff), "Sum (%s) every %d for mat (%g %g %g %g)",
                scpstring(xgc->data.set.so.sums[i].component),
                xgc->data.set.so.sums[i].step, xgc->data.set.so.sums[i].layered_epsilon/EPSILON_0, xgc->data.set.so.sums[i].layered_mu/MU_0, xgc->data.set.so.sums[i].layered_sigma, xgc->data.set.so.sums[i].layered_sigast);
                } else
                g_snprintf(buff, sizeof(buff), "Sum (%s) every %d for %s",
                scpstring(xgc->data.set.so.sums[i].component),
                xgc->data.set.so.sums[i].step, xgc->data.set.so.sums[i].string);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_SOUT+i, COLUMN_SHOW_TOGGLE, TRUE, -1);*/

                if (xgc->nfilestoshow < 99) {
                    xgc->filestoshow[xgc->nfilestoshow] = g_strdup(xgc->data.set.so.sums[i].filename);
                    xgc->formattoshow[xgc->nfilestoshow] = xgc->data.set.so.sums[i].component;
                    xgc->nfilestoshow += 1;
                }
            }

            /* Force output */
            for (i = 0; i < MIN(xgc->data.set.so.nforces, TREE_MAXENT); i++) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_FOUT, NOL_OUTPUT_PROPS + xgc->data.set.so.npnts + xgc->data.set.so.nimgs + xgc->data.set.so.nplns + xgc->data.set.so.ncubs + xgc->data.set.so.nsums + i);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_outforce[i], COLUMN_ID, SET_FOUT + i, COLUMN_SHOW_TOGGLE, TRUE, -1);

                /*gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "Force every %d in box (%d %d %d ... %d %d %d) to %s", xgc->data.set.so.forces[i].step,
                xgc->data.set.so.forces[i].cuboid_i0, xgc->data.set.so.forces[i].cuboid_j0, xgc->data.set.so.forces[i].cuboid_k0, xgc->data.set.so.forces[i].cuboid_in, xgc->data.set.so.forces[i].cuboid_jn, xgc->data.set.so.forces[i].cuboid_kn,
                xgc->data.set.so.forces[i].filename);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_FOUT+i, COLUMN_SHOW_TOGGLE, TRUE, -1);*/

                if (xgc->nfilestoshow < 99) {
                    xgc->filestoshow[xgc->nfilestoshow] = g_strdup(xgc->data.set.so.forces[i].filename);
                    xgc->formattoshow[xgc->nfilestoshow] = 100;
                    xgc->nfilestoshow += 1;
                }
            }

            gtk_tree_store_append(xgc->ts_par, &iter, NULL);
            gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, "NFFF", COLUMN_CHECK, FALSE, COLUMN_ID, -6, COLUMN_SHOW_TOGGLE, FALSE, -1);

            if (xgc->data.set.sf.nrs || xgc->data.set.sf.nareas) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFF, 0);    /* NFFF box*/
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_nfff, COLUMN_ID, SET_NFFF, COLUMN_SHOW_TOGGLE, TRUE, -1);

                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFF, 1);    /* NFFF skips */
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_NFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);

                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFF, 2);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_NFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);

                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFF, 3);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_NFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);

                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFF, 4);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_NFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);

                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFF, 5);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_NFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);

                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFF, 6);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_NFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);

                /*gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "NFFF box: %d %d %d ... %d %d %d", xgc->data.set.sf.cuboid_i0, xgc->data.set.sf.cuboid_j0, xgc->data.set.sf.cuboid_k0, xgc->data.set.sf.cuboid_in, xgc->data.set.sf.cuboid_jn, xgc->data.set.sf.cuboid_kn);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_NFFF, COLUMN_SHOW_TOGGLE, TRUE, -1);

                if (xgc->data.set.sf.cuboid_boundary_skipi0) {

                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "NFFF skip i0: %d %d ... %d %d", xgc->data.set.sf.skipi0_jmin, xgc->data.set.sf.skipi0_kmin, xgc->data.set.sf.skipi0_jmax, xgc->data.set.sf.skipi0_kmax);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_NFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);
                }
                if (xgc->data.set.sf.cuboid_boundary_skipin) {
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "NFFF skip in: %d %d ... %d %d", xgc->data.set.sf.skipin_jmin, xgc->data.set.sf.skipin_kmin, xgc->data.set.sf.skipin_jmax, xgc->data.set.sf.skipin_kmax);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_NFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);
                }

                if (xgc->data.set.sf.cuboid_boundary_skipj0) {
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "NFFF skip j0: %d %d ... %d %d", xgc->data.set.sf.skipj0_imin, xgc->data.set.sf.skipj0_kmin, xgc->data.set.sf.skipj0_imax, xgc->data.set.sf.skipj0_kmax);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_NFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);
                }

                if (xgc->data.set.sf.cuboid_boundary_skipjn) {
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "NFFF skip jn: %d %d ... %d %d", xgc->data.set.sf.skipjn_imin, xgc->data.set.sf.skipjn_kmin, xgc->data.set.sf.skipjn_imax, xgc->data.set.sf.skipjn_kmax);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_NFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);
                }

                if (xgc->data.set.sf.cuboid_boundary_skipk0) {
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "NFFF skip k0: %d %d ... %d %d", xgc->data.set.sf.skipk0_imin, xgc->data.set.sf.skipk0_jmin, xgc->data.set.sf.skipk0_imax, xgc->data.set.sf.skipk0_jmax);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_NFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);
                }

                if (xgc->data.set.sf.cuboid_boundary_skipkn) {
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "NFFF skip kn: %d %d ... %d %d", xgc->data.set.sf.skipkn_imin, xgc->data.set.sf.skipkn_jmin, xgc->data.set.sf.skipkn_imax, xgc->data.set.sf.skipkn_jmax);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_NFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);
                }*/

                /* NFFF points */
                for (i = 0; i < xgc->data.set.sf.nrs; i++) {
                    pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFFFP, NOL_NFFF_PROPS + i);
                    gtk_tree_store_append(xgc->ts_par, &child, &iter);
                    gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_nfff_point[i], COLUMN_ID, SET_NFFFP + i, COLUMN_SHOW_TOGGLE, TRUE, -1);
                }

                /*for (i=0; i<xgc->data.set.sf.nrs; i++) {
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                if (xgc->data.set.sf.source_filename[i]) {
                g_snprintf(buff, sizeof(buff), "Far field point at (%d %d %d) to %s", xgc->data.set.sf.ri[i], xgc->data.set.sf.rj[i], xgc->data.set.sf.rk[i], xgc->data.set.sf.source_filename[i]);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_NFFFP+i, COLUMN_SHOW_TOGGLE, TRUE, -1);
                }
                }*/

                /* NFFF areas */

                for (i = 0; i < MIN(xgc->data.set.sf.nareas, TREE_MAXENT); i++) {
                    pt_assemble_row_text(xgc, buff, sizeof(buff), SET_NFAREA, NOL_NFFF_PROPS + xgc->data.set.sf.nrs + i);
                    gtk_tree_store_append(xgc->ts_par, &child, &iter);
                    gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_NFAREA + i, COLUMN_SHOW_TOGGLE, FALSE, -1);
                }

                /*for (i=0; i<MIN(xgc->data.set.sf.nareas, TREE_MAXENT); i++) {
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "Far field area (%dx%d)", xgc->data.set.sf.area_thetares[i], xgc->data.set.sf.area_phires[i]);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_NFAREA+i, -1);
                }*/
            }

            gtk_tree_store_append(xgc->ts_par, &iter, NULL);
            gtk_tree_store_set(xgc->ts_par, &iter, COLUMN_PARAMETER, "Periodic NFFF", COLUMN_CHECK, FALSE, COLUMN_ID, -7, COLUMN_SHOW_TOGGLE, FALSE, -1);

            if (xgc->data.set.spf.nrs || xgc->data.set.spf.nareas) {
                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_PNFFF, 0);    /* PNFFF box*/
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_pnfff, COLUMN_ID, SET_PNFFF, COLUMN_SHOW_TOGGLE, TRUE, -1);

                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_PNFFF, 1);    /* PNFFF skips */
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_PNFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);

                pt_assemble_row_text(xgc, buff, sizeof(buff), SET_PNFFF, 2);
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_PNFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);

                /*gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "Periodic NFFF box: %d %d %d ... %d %d %d", xgc->data.set.spf.cuboid_i0, xgc->data.set.spf.cuboid_j0, xgc->data.set.spf.cuboid_k0, xgc->data.set.spf.cuboid_in, xgc->data.set.spf.cuboid_jn, xgc->data.set.spf.cuboid_kn);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_PNFFF, COLUMN_SHOW_TOGGLE, TRUE, -1);

                if (xgc->data.set.spf.cuboid_boundary_skipk0) {
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "Periodic NFFF skip k0: %d %d ... %d %d", xgc->data.set.spf.skipk0_imin, xgc->data.set.spf.skipk0_jmin, xgc->data.set.spf.skipk0_imax, xgc->data.set.spf.skipk0_jmax);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_PNFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);
                }

                if (xgc->data.set.spf.cuboid_boundary_skipkn) {
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "Periodic NFFF skip kn: %d %d ... %d %d", xgc->data.set.spf.skipkn_imin, xgc->data.set.spf.skipkn_jmin, xgc->data.set.spf.skipkn_imax, xgc->data.set.spf.skipkn_jmax);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_PNFFF, COLUMN_SHOW_TOGGLE, FALSE, -1);
                }*/

                /* Periodic NFFF points */

                for (i = 0; i < MIN(xgc->data.set.spf.nrs, TREE_MAXENT); i++) {
                    pt_assemble_row_text(xgc, buff, sizeof(buff), SET_PNFFFP, NOL_PNFFF_PROPS + i);
                    gtk_tree_store_append(xgc->ts_par, &child, &iter);
                    gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, xgc->data.is_pnfff_point[i], COLUMN_ID, SET_PNFFFP + i, COLUMN_SHOW_TOGGLE, TRUE, -1);
                }

                /*for (i=0; i<MIN(xgc->data.set.spf.nrs, TREE_MAXENT); i++) {
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                if ( xgc->data.set.spf.source_filename[i]) {
                g_snprintf(buff, sizeof(buff), "Far field point at (%d %d %d) to %s", xgc->data.set.spf.ri[i], xgc->data.set.spf.rj[i], xgc->data.set.spf.rk[i], xgc->data.set.spf.source_filename[i]);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_PNFFFP+i, COLUMN_SHOW_TOGGLE, TRUE, -1);
                }
                }*/

                /* Periodic NFFF areas */

                for (i = 0; i < MIN(xgc->data.set.spf.nareas, TREE_MAXENT); i++) {
                    pt_assemble_row_text(xgc, buff, sizeof(buff), SET_PNFAREA, NOL_PNFFF_PROPS + xgc->data.set.spf.nrs + i);
                    gtk_tree_store_append(xgc->ts_par, &child, &iter);
                    gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, FALSE, COLUMN_ID, SET_PNFAREA + i, COLUMN_SHOW_TOGGLE, FALSE, -1);
                }

                /*for (i = 0; i<MIN(xgc->data.set.spf.nareas, TREE_MAXENT); i++) {
                gtk_tree_store_append(xgc->ts_par, &child, &iter);
                g_snprintf(buff, sizeof(buff), "Periodic far field area (%dx%d)", xgc->data.set.spf.area_thetares[i], xgc->data.set.spf.area_phires[i]);
                gtk_tree_store_set(xgc->ts_par, &child, COLUMN_PARAMETER, buff, COLUMN_CHECK, TRUE, COLUMN_ID, SET_PNFAREA+i, COLUMN_SHOW_TOGGLE, TRUE, -1);
                }*/
            }

            if (g_strcmp0(xgc->data.matfilename, xgc->data.set.sm.in_vector_filename) != 0) {  //new matfile loaded
                if (xgc->data.set.sm.in_vector_filename != NULL) {
                    //      printf("%s vs %s\n", xgc->data.matfilename, xgc->data.set.sm.in_vector);
                    xgc->data.matfilename = g_strdup(xgc->data.set.sm.in_vector_filename);
                }
                /*load material file*/
                if (xgc->data.set.sm.in_vector_filename != NULL) {
                    if (g_file_get_contents(xgc->data.matfilename, &contents, NULL, NULL)) {
                        g_snprintf(buff, sizeof(buff), "Material successfully loaded from %s", xgc->data.matfilename);
                        gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);

                        gtk_text_buffer_set_text(xgc->tb_mat, contents, -1);
                        /* 2018-05-31 - commented - gtk_text_buffer_set_text() already calls mat_file_changed() */
                        /* mat_file_changed(xgc, TRUE); */
                    }
                }
            }

            active = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->graph_combo));

            for (i = 0; i < 100; i++)
                gtk_combo_box_remove_text(GTK_COMBO_BOX(xgc->graph_combo), 0);

            for (i = 0; i < xgc->nfilestoshow; i++) {
                if (xgc->filestoshow[i] != NULL)
                    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc->graph_combo), xgc->filestoshow[i]);
            }

            if (active >= 0)
                gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->graph_combo), active);
            else if (xgc->nfilestoshow > 0)
                gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->graph_combo), 0);


            active = gtk_combo_box_get_active(GTK_COMBO_BOX(xgc->image_combo));

            for (i = 0; i < 100; i++)
                //gtk_combo_box_remove_text(GTK_COMBO_BOX(xgc->image_combo), 0);
                gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(xgc->image_combo), 0);
            for (i = 0; i < xgc->nimagestoshow; i++) {
                if (xgc->imagestoshow[i] != NULL)
                    gtk_combo_box_append_text(GTK_COMBO_BOX(xgc->image_combo), xgc->imagestoshow[i]);
                    //gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(xgc->image_combo), xgc->imagestoshow[i]);
            }
            if (active >= 0)
                gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->image_combo), active);
            else if (xgc->nimagestoshow > 0)
                gtk_combo_box_set_active(GTK_COMBO_BOX(xgc->image_combo), 0);

            g_snprintf(buff, sizeof(buff), "Ready");
            gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);

            gtk_widget_queue_draw(xgc->view_scene);

            xgc->par_tree_success = TRUE;
            xgc->par_file_success = TRUE;            

            gtk_widget_set_sensitive(xgc->page_par_tree, xgc->par_tree_success);            
        } else {
            g_snprintf(buff, sizeof(buff), MSG_SB_ERROR_PARSING_PARFILE);
            gtk_statusbar_push(GTK_STATUSBAR(xgc->statusbar), xgc->statusbar_context_id, buff);

            gtk_widget_queue_draw(xgc->view_scene);

            xgc->par_tree_success = FALSE;
            xgc->par_file_success = FALSE;

            gtk_widget_set_sensitive(xgc->page_par_tree, xgc->par_tree_success);
            
            dialog = gtk_message_dialog_new(GTK_WINDOW(xgc->toplevel),
                                            GTK_DIALOG_DESTROY_WITH_PARENT,
                                            GTK_MESSAGE_WARNING,
                                            GTK_BUTTONS_OK,
                                            MSG_MB_ERROR_PARSING_PARFILE);
            gtk_dialog_run(GTK_DIALOG(dialog));
            gtk_widget_destroy(dialog);

            gtk_notebook_set_current_page(GTK_NOTEBOOK(xgc->param_notebook), 1);
        }        
        //g_free(filename);

        pt_restore_expanded_rows_and_selection(xgc);
    }
} /* pt_create_tree */