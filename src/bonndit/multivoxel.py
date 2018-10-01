import errno
import multiprocessing as mp
import os

import nibabel as nib
import numpy as np
from dipy.core.gradients import gradient_table
from dipy.reconst.multi_voxel import MultiVoxelFit
from tqdm import tqdm


def fit_helper(args_kwargs):
    args = args_kwargs[1:-1]
    kwargs = args_kwargs[-1]
    func = args_kwargs[0]
    return func(*args, **kwargs)


class MultiVoxelFitter(object):
    def __init__(self, model, cpus=1, verbose=False, desc=""):
        self.model = model
        self.cpus = cpus
        self.verbose = verbose
        self.desc = desc

    def fit(self, fit_func, data, per_voxel_data, mask=None):

        space = data.shape[:-1]

        if mask is None:
            mask = np.ones(space)
        if mask.shape != space:
            raise ValueError("mask and data shape do not match")

        # Convert integer to boolean mask
        mask = np.ma.make_mask(mask)

        # 1000 chunks for the progressbar to run smoother
        chunksize = max(1, int(np.prod(data.shape[:-1]) / 1000))

        # collect kwargs which are the same for all voxels
        args_kwargs = []

        for ijk in np.ndindex(*data.shape[:-1]):
            if mask[ijk]:
                # collect kwargs per voxel if specified in the kwargs
                per_voxel_kwargs = {key: per_voxel_data[key][ijk]
                                    for key in per_voxel_data}

                args_kwargs.append((fit_func, data[ijk], per_voxel_kwargs))
                # print(new_kwargs)

        if self.cpus == 1:
            coeffs = list(tqdm(map(fit_helper, args_kwargs),
                               total=np.prod(data.shape[:-1]),
                               disable=not self.verbose,
                               desc=self.desc))
        else:
            with mp.Pool(self.cpus) as p:
                coeffs = list(tqdm(p.imap(fit_helper, args_kwargs,
                                          chunksize),
                                   total=np.prod(data.shape[:-1]),
                                   disable=not self.verbose,
                                   desc=self.desc))

        fit_array = np.empty(data.shape[:-1], dtype=object)
        fit_array[mask] = coeffs

        return MultiVoxel(self.model, fit_array, mask)

class MultiVoxel(MultiVoxelFit):
    def __init__(self, model, fit_array, mask):
        super().__init__(model, fit_array, mask)

    @classmethod
    def load(cls, filepath, model_class, fit_class):
        filecontent = np.load(filepath)

        gtab = gradient_table(filecontent['bvals'], filecontent['bvecs'])
        model_params = {key: filecontent[key] for key in filecontent.keys()
                        if key not in ['data', 'mask', 'bvals', 'bvecs']}
        model = model_class(gtab, **model_params)

        data = filecontent['data']
        mask = filecontent['mask']

        fit_array = np.empty(data.shape[:-1], dtype=object)
        for ijk in np.ndindex(*data.shape[:-1]):
            if mask[ijk]:
                fit_array[ijk] = fit_class(data[ijk])

        return cls(model, fit_array, mask)

    def save(self, filepath, affine=None, type='npz'):
        try:
            os.makedirs(os.path.dirname(filepath))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        if affine is None:
            affine = np.zeros((4, 4))

        data = self.fit_array.coeffs
        mask = self.mask
        if type == 'npz':
            np.savez(filepath, data=data, mask=mask, **self._model_params)


        elif type == 'nii':
            img = nib.Nifti1Image(data, affine=affine)
            nib.save(img, filepath)
