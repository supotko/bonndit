#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
================================================================================
Compute fiber orientation distribution functions for diffusion weighted MRI data
================================================================================
"""

import argparse
import os

import numpy as np
from bonndit.michi import fields, dwmri
from bonndit.shore import ShoreModel, ShoreFit
# from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table


# import nibabel as nib


# To Do: Add option to first compute the diffustion tensors which are needed to estimate the response functions.

# To Do: Add option to automatically build masks for csf gm wm.

# To Do: Handle nrrd as well as nii input

# To Do: Enable saving in different output formats

def main():
    parser = argparse.ArgumentParser(
        description='This script computes fiber orientation distribution functions (fODFs) \
        as described in "Versatile, Robust and Efficient Tractography With Constrained Higher \
        Order Tensor fODFs" by Ankele et al. (2017)', add_help=False)

    parser.add_argument('indir',
                        help='Folder containing all required input files.')

    optional = parser.add_argument_group('optional')
    optional.add_argument("-h", "--help", action="help", help="show this help message and exit")
    optional.add_argument('-o', '--outdir',
                          help='folder in which the output will be saved')

    flags = parser.add_argument_group('flags (optional)', '')
    flags.add_argument('-v', '--verbose', action='store_true',
                       help='show a progress bars for calculation of the response function and the deconvolution')
    flags.add_argument('-R', '--responseonly', action='store_true',
                       help='calculate and save only the response functions')
    flags.add_argument('-V', '--volumes', action='store_true',
                       help='output the volume fractions (csf/gm/wm) after deconvolution')

    shoreopts = parser.add_argument_group('shore options (optional)', '')
    shoreopts.add_argument('-r', '--order', default=4, type=int,
                        help='Order of the shore basis')
    shoreopts.add_argument('-z', '--zeta', default=700, type=float,
                           help='Radial scaling factor')
    shoreopts.add_argument('-t', '--tau', default=1 / (4 * np.pi ** 2), type=float,
                           help='q-scaling')
    shoreopts.add_argument('-f', '--fawm', default=0.7, type=float,
                           help='The WM FA threshold')

    deconvopts = parser.add_argument_group('deconvolution options (optional)', '')
    deconvopts.add_argument('-c', '--constraint', choices=['hpsd', 'nonneg', 'none'], default='hpsd',
                            help='constraint for the fODFs')

    filenaming = parser.add_argument_group('file naming (optional)', 'Specify custom names for output files')
    filenaming.add_argument('-S', '--responseout', default='response.npz',
                            help='response function output name - filetype: .npz')
    filenaming.add_argument('-O', '--fodfout', default='fodf.nrrd',
                            help='fODF filename - filetype: .nrrd / .nii / .nii.gz')
    filenaming.add_argument('-W', '--whitematter', default='wmvolume.nrrd',
                            help='wm volume filename - filetype: .nrrd / .nii / .nii.gz')
    filenaming.add_argument('-G', '--graymatter', default='gmvolume.nrrd',
                            help='gm volume filename - filetype: .nrrd / .nii / .nii.gz')
    filenaming.add_argument('-C', '--cerebrospinalfluid', default='csfvolume.nrrd',
                            help='csf volume filename - filetype: .nrrd / .nii / .nii.gz')

    args = parser.parse_args()
    order = args.order
    zeta = args.zeta
    tau = args.tau
    fawm = args.fawm
    verbose = args.verbose
    indir = args.indir
    if not args.outdir:
        outdir = indir
    else:
        outdir = args.outdir

    # Load fractional anisotropy
    # dti_fa = nib.load(os.path.join(indir, "dti_FA.nii.gz"))
    dti_fa, meta = fields.load_scalar(os.path.join(indir, "dti_FA.nii.gz"))

    # Load DTI mask
    # dti_mask = nib.load(os.path.join(indir, "mask.nii.gz"))
    dti_mask, _ = fields.load_scalar(os.path.join(indir, "mask.nii.gz"))

    # Load and adjust tissue segmentation masks
    # csf_mask = nib.load(os.path.join(indir, "fast_pve_0.nii.gz"))
    csf_mask, _ = fields.load_scalar(os.path.join(indir, "fast_pve_0.nii.gz"))
    # gm_mask = nib.load(os.path.join(indir, "fast_pve_1.nii.gz"))
    gm_mask, _ = fields.load_scalar(os.path.join(indir, "fast_pve_1.nii.gz"))
    # wm_mask = nib.load(os.path.join(indir, "fast_pve_2.nii.gz"))
    wm_mask, _ = fields.load_scalar(os.path.join(indir, "fast_pve_2.nii.gz"))

    # dti_vecs = nib.load(os.path.join(indir, "dti_V1.nii.gz"))
    dti_vecs, _ = fields.load_vector(os.path.join(indir, "dti_V1.nii.gz"))

    # data = nib.load(os.path.join(indir, "data.nii.gz"))

    # bvals, bvecs = read_bvals_bvecs(os.path.join(indir, "bvals"),
    #                                os.path.join(indir, "bvecs"))
    # gtab = gradient_table(bvals, bvecs)
    data, gtabm, meta = dwmri.load(os.path.join(indir, "data.nii.gz"))

    gtab = gradient_table(gtabm.bvals, gtabm.bvecs)

    # If not only the response has to be computed
    if not args.responseonly:
        # Check if response functions already exist.
        if os.path.exists(os.path.join(outdir, args.responseout)):
            fit = ShoreFit.old_load(os.path.join(outdir, args.responseout))
            if verbose:
                print('Loaded existing response functions.')
    # Force recalculate the response if response only is specified
    else:
        model = ShoreModel(gtab, order, zeta, tau)
        fit = model.fit(data, wm_mask, gm_mask, csf_mask, dti_mask,
                        dti_fa, dti_vecs, fawm, verbose=verbose)
        fit.old_save(os.path.join(outdir, args.responseout))

    # Check if the response only flag is set.
    if not args.responseonly:
        out, wmout, gmout, csfout, mask = fit.fodf(data, verbose=verbose, pos=args.constraint)

        fields.save_tensor(os.path.join(args.outdir, args.fodfout), out, mask=mask, meta=meta)

        # Check if the volumes flag is set.
        if args.volumes:
            fields.save_scalar(os.path.join(args.outdir, args.whitematter), wmout, meta)
            fields.save_scalar(os.path.join(args.outdir, args.graymatter), gmout, meta)
            fields.save_scalar(os.path.join(args.outdir, args.cerebrospinalfluid), csfout, meta)




if __name__ == "__main__":
    main()
