from __future__ import division

import errno
import logging
import multiprocessing as mp
import os
import sys

import numpy as np
import numpy.linalg as la
from dipy.core.geometry import cart2sphere
from dipy.core.gradients import gradient_table
from dipy.reconst.shm import real_sph_harm
from tqdm import tqdm

from bonndit.michi import esh
from .gradients import gtab_reorient

try:
    from itertools import imap
except ImportError:
    # For Python 3 imap was removed as gloabl map now returns an iterator
    imap = map


class SphericalHarmonicsModel(object):
    def __init__(self, gtab, order=4):
        """

        :param gtab:
        :param order:
        """
        self.gtab = gtab
        self.order = order

        # A gradient table without small bvalues,
        # depends on b0_threshold of gtab
        self._gtab = gradient_table(self.gtab.bvals[~self.gtab.b0s_mask],
                                    self.gtab.bvecs[~self.gtab.b0s_mask, :])

        # Ignore division by zero warning
        # dipy.core.geometry.cart2sphere -> theta = np.arccos(z / r)
        with np.errstate(divide='ignore', invalid='ignore'):
            self.sh_m = esh_matrix(self.order, self.gtab)

    def _fit_helper(self, data_vecs, rcond=None):
        """

        :param data_vecs:
        :param rcond:
        :return:
        """
        signal, vec = data_vecs[0], data_vecs[1]

        if vec is not None:
            with np.errstate(divide='ignore', invalid='ignore'):
                sh_m = esh_matrix(self.order, gtab_reorient(self._gtab, vec))

        else:
            sh_m = self.sh_m

        return la.lstsq(sh_m, signal, rcond)[0]

    def fit(self, data, vecs=None, verbose=False, cpus=1, desc=""):
        """

        :param data:
        :param vecs:
        :param verbose:
        :param cpus:
        :param desc:
        :return:
        """

        # 1000 chunks for the progressbar to run smoother
        chunksize = max(1, int(np.prod(data.shape[:-1]) / 1000))

        # If no vectors are specified create array of Nones for iteration.
        if type(vecs) != np.ndarray and vecs is None:
            vecs = np.empty(data.shape[:-1], dtype=object)

        # Calculate average b0 signal in white matter
        b0_avg = np.mean(data[..., self.gtab.b0s_mask])

        # Remove small bvalues, depends on b0_threshold of gtab
        data = data[..., ~self.gtab.b0s_mask]

        # Iterate over the data indices; show progress with tqdm
        # multiple processes for python > 3
        if sys.version_info[0] < 3 or cpus == 1:
            sh_coeff = list(tqdm(imap(self._fit_helper,
                                      zip(list(data), list(vecs))),
                                 total=np.prod(data.shape[:-1]),
                                 disable=not verbose,
                                 desc=desc))
        else:
            with mp.Pool(cpus) as p:
                sh_coeff = list(tqdm(p.imap(self._fit_helper,
                                            zip(list(data), list(vecs)),
                                            chunksize),
                                     total=np.prod(data.shape[:-1]),
                                     disable=not verbose,
                                     desc=desc))

        return SphericalHarmonicsFit(self, np.array(sh_coeff), b0_avg)


class SphericalHarmonicsFit(object):
    def __init__(self, model, coefs, b0_avg):
        """

        :param model:
        :param coefs:
        :param kernel:
        """
        self.model = model
        self.coefs = coefs
        self.b0_avg = b0_avg

        self.order = model.order
        self.gtab = model.gtab

    @classmethod
    def load(cls, filepath):
        """ Load a precalculated SphericalHarmonicsFit object from a file.

        :param filepath: path to the saved SphericalHarmonicsFit object
        :return: SphericalHarmonicsFit object which contains wm response
        function
        """
        response = np.load(filepath)

        gtab = gradient_table(response['bvals'], response['bvecs'])
        model = SphericalHarmonicsModel(gtab, response['order'])

        return cls(model, response['coefs'], response['b0_avg'])

    def save(self, filepath):
        """ Save a mtShoreFit object to a file.

        :param filepath: path to the file
        """
        try:
            os.makedirs(os.path.dirname(filepath))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        np.savez(filepath, coefs=self.coefs,
                 order=self.order, bvals=self.gtab.bvals,
                 bvecs=self.gtab.bvecs, b0_avg=self.b0_avg)


class ShResponseEstimator(object):
    def __init__(self, gtab, order=4):
        """

        :param gtab:
        :param order:
        """
        self.gtab = gtab
        self.order = order

    def fit(self, data, dti_vecs, wm_mask, verbose=False, cpus=1):
        """

        :param data:
        :param dti_vecs:
        :param wm_mask:
        :param verbose:
        :param cpus:
        :return:
        """
        # Check if tissue masks give at least a single voxel
        if np.sum(wm_mask.get_data()) < 1:
            raise ValueError('No white matter voxels specified by wm_mask. '
                             'A corresponding response can not be computed.')

        # Select white matter voxels
        wm_voxels = data.get_data()[wm_mask.get_data() == 1]
        wm_vecs = dti_vecs.get_data()[wm_mask.get_data() == 1]

        # Calculate white matter response
        wm_sh_coefs = SphericalHarmonicsModel(
            self.gtab, self.order).fit(wm_voxels, wm_vecs, verbose=verbose,
                                       cpus=cpus,
                                       desc='WM response').coefs
        wm_sh_coef = self.sh_accumulate(wm_sh_coefs)
        signal_wm = self.sh_compress(wm_sh_coef)

        return ShResponse(self, signal_wm)

    def sh_accumulate(self, sh_coefs):
        """

        :param sh_coefs:
        :return:
        """
        sh_accum = np.zeros_like(sh_coefs[0])
        accum_count = 0

        # Iterate over the data indices
        for i in np.ndindex(*sh_coefs.shape[:-1]):
            sh_accum += sh_coefs[i]
            accum_count += 1
        if accum_count == 0:
            return sh_accum

        # Do not show divide by zero warnings
        with np.errstate(divide='ignore', invalid='ignore'):
            return sh_accum / accum_count

    def sh_compress(self, coefs):
        """ Compress the spherical harmonics coefficients

        An axial symetric response function aligned to the z-axis can be
        described fully using only the z-rotational part of the spherical
        harmonics coefficients. This functions selects the zonal harmonics with
        even order from an array with spherical harmonics coefficients.

        :param coefs: spherical harmonics coefficients
        :return: z-rotational part of the spherical harmonics coefficients
        """
        zonal_coefs = np.zeros(esh.get_kernel_size(self.order))

        counter = 0
        for l in range(0, self.order + 1):
            counter = counter + l
            if l % 2 == 0:
                zonal_coefs[int(l / 2)] = coefs[counter]

        # This is what happens above
        # r[0] = coefs[0]
        # r[1] = coefs[3]
        # r[2] = coefs[10]
        # ...

        return zonal_coefs


class ShResponse(object):
    def __init__(self, model, sh_coef, kernel="rank1"):
        """

        :param model:
        :param sh_coef:
        :param kernel:
        """
        self.model = model
        self.gtab = model.gtab
        self.order = model.order

        self.wm_response = sh_coef

        # A gradient table without small bvalues,
        # depends on b0_threshold of gtab
        self._gtab = gradient_table(self.gtab.bvals[~self.gtab.b0s_mask],
                                    self.gtab.bvecs[~self.gtab.b0s_mask, :])

        # The deconvolution kernels are computed in set_kernel
        self.kernel_type = kernel
        self.kernel_wm = None
        self.set_kernel(kernel)

    def set_kernel(self, kernel):
        """

        :param kernel:
        :return:
        """
        # Get deconvolution kernel
        if kernel == "rank1":
            self.kernel_wm = esh.make_kernel_rank1(self.wm_response)
        elif kernel == "delta":
            self.kernel_wm = esh.make_kernel_delta(self.wm_response)
        else:
            msg = "{} is not a valid option for kernel. " \
                  "Use 'rank1' or 'delta'.".format(kernel)
            raise ValueError(msg)


    @classmethod
    def load(cls, filepath):
        """ Load a precalculated mtShoreFit object from a file.

        :param filepath: path to the saved mtShoreFit object
        :return: mtShoreFit object which contains response functions for white
        matter, gray matter and CSF
        """
        response = np.load(filepath)

        gtab = gradient_table(response['bvals'], response['bvecs'])
        model = ShResponseEstimator(gtab, response['order'])

        return cls(model, response['wm_resp'])

    def save(self, filepath):
        """ Save a mtShoreFit object to a file.

        :param filepath: path to the file
        """
        try:
            os.makedirs(os.path.dirname(filepath))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        np.savez(filepath, wm_resp=self.wm_response, order=self.order,
                 bvals=self.gtab.bvals, bvecs=self.gtab.bvecs)

    def odf(self, data, pos='hpsd', mask=None, kernel="rank1", verbose=False,
            cpus=1):
        """

        :param data:
        :param pos:
        :param mask:
        :param kernel:
        :param verbose:
        :param cpus:
        :return:
        """
        if self.kernel_type != kernel:
            self.set_kernel(kernel)

        data = data.get_data()
        space = data.shape[:-1]

        if not mask:
            mask = np.ones(space)
        else:
            mask = mask.get_data()
        # Convert integer to boolean mask
        mask = np.ma.make_mask(mask)

        # Create convolution matrix
        conv_mat = self.sh_convolution_matrix(kernel)
        with np.errstate(divide='ignore', invalid='ignore'):
            cond_number = la.cond(conv_mat)
            logging.info('Condition number of convolution matrtix:' +
                         str(cond_number))

        # 1000 chunks for the progressbar to run smoother
        chunksize = max(1, int(np.prod(data.shape[:-1]) / 1000))

        # Deconvolve the DWI signal
        deconv = {'none': self.deconvolve, 'hpsd': self.deconvolve_hpsd,
                  'nonneg': self.deconvolve_nonneg}
        data = data[mask, :]
        try:
            func = deconv[pos]
            if sys.version_info[0] < 3 or cpus == 1:
                result = list(
                    tqdm(imap(partial(func, conv_matrix=conv_mat),
                              data),
                         total=np.prod(data.shape[:-1]),
                         disable=not verbose,
                         desc='Optimization'))
            else:
                with mp.Pool(cpus) as p:
                    result = list(tqdm(p.imap(partial(func,
                                                      conv_matrix=conv_mat),
                                              data, chunksize=chunksize),
                                       total=np.prod(data.shape[:-1]),
                                       disable=not verbose,
                                       desc='Optimization'))
        except KeyError:
            raise ValueError(
                ('"{}" is not supported as a constraint, please' +
                 ' choose from [hpsd, nonneg, none]').format(pos))

        # Return fODFs and Volume fractions as separate numpy.ndarray objects
        NN = esh.LENGTH[self.order]
        out = np.zeros(space + (NN,))
        wmout = np.zeros(space)

        out[mask, :] = [esh.esh_to_sym(x[:NN]) for x in result]
        f = self.kernel_csf[0][0] / max(self.signal_csf[0], 1e-10)
        wmout[mask] = [x[0] * f for x in result]

        return out, wmout

    def deconvolve(self, data, conv_matrix):
        """

        :param data:
        :param conv_matrix:
        :return:
        """
        NN = esh.LENGTH[self.order]
        deconvolution_result = np.zeros(data.shape[:-1] + (NN + 2,))

        for i in np.ndindex(*data.shape[:-1]):
            signal = data[i]
            deconvolution_result[i] = la.lstsq(conv_matrix, signal,
                                               rcond=None)[0]

        return deconvolution_result

    def deconvolve_hpsd(self, data, conv_matrix):
        """

        :param data:
        :param conv_matrix:
        :return:
        """
        NN = esh.LENGTH[self.order]
        deconvolution_result = np.zeros(data.shape[:-1] + (NN + 2,))

        cvxopt.solvers.options['show_progress'] = False
        # set up QP problem from normal equations
        P = cvxopt.matrix(np.ascontiguousarray(np.dot(conv_matrix.T,
                                                      conv_matrix)))

        # positive definiteness constraint on ODF
        ind = tensor.H_index_matrix(self.order).reshape(-1)
        N = len(ind)

        # set up positive definiteness constraints
        G = np.zeros((N + 2, NN + 2))
        # constrain GM/CSF VFs to be non-negative: orthant constraints
        G[0, NN] = -1
        G[1, NN + 1] = -1
        esh2sym = esh.esh_to_sym_matrix(self.order)
        for i in range(N):
            G[i + 2, :NN] = -esh2sym[ind[i], :]
        h = np.zeros(N + 2)

        # initialize with partly GM, CSF, and isotropic ODF
        init = np.zeros(NN + 2)
        init[0] = 0.3
        init[1] = 0.3
        init[2] = 0.3 * self.signal_csf[0] / self.kernel_csf[0][0]
        init = cvxopt.matrix(np.ascontiguousarray(init))

        G = cvxopt.matrix(np.ascontiguousarray(G))
        h = cvxopt.matrix(np.ascontiguousarray(h))

        for i in np.ndindex(*data.shape[:-1]):
            signal = data[i]
            q = cvxopt.matrix(np.ascontiguousarray(-1 * np.dot(conv_matrix.T,
                                                               signal)))

            # NS = len(np.array(T{4,6,8}.TT).reshape(-1))
            NS = tensor.LENGTH[self.order // 2]

            # first two are orthant constraints, rest positive definiteness
            dims = {'l': 2, 'q': [], 's': [NS]}

            # This init stuff is a HACK.
            # It empirically removes some isolated failure cases
            # first, allow it to use its own initialization
            try:
                sol = cvxopt.solvers.coneqp(P, q, G, h, dims)
            except ValueError as e:
                logging.error("Error with cvxopt initialization: {}".format(e))
                return np.zeros(NN + 2)
            if sol['status'] != 'optimal':
                # try again with our initialization
                try:
                    sol = cvxopt.solvers.coneqp(P, q, G, h, dims,
                                                initvals={'x': init})
                except ValueError as e:
                    logging.error("Error with custom initialization: "
                                  "{}".format(e))
                    return np.zeros(NN + 2)
                if sol['status'] != 'optimal':
                    logging.debug('Optimization unsuccessful - '
                                  'Constraint: {}'.format('hpsd'))

            deconvolution_result[i] = np.array(sol['x'])[:, 0]

        return deconvolution_result

    def deconvolve_nonneg(self, data, conv_matrix):
        """

        :param data:
        :param conv_matrix:
        :return:
        """
        NN = esh.LENGTH[self.order]
        deconvolution_result = np.zeros(data.shape[:-1] + (NN + 2,))

        cvxopt.solvers.options['show_progress'] = False
        # set up QP problem from normal equations
        P = cvxopt.matrix(np.ascontiguousarray(np.dot(conv_matrix.T,
                                                      conv_matrix)))

        # set up non-negativity constraints
        NC = LOTS_OF_DIRECTIONS.shape[0]
        G = np.zeros((NC + 2, NN + 2))
        G[:NC, :NN] = -esh.matrix(self.order, LOTS_OF_DIRECTIONS)

        # also constrain GM/CSF VFs to be non-negative
        G[NC, NN] = -1
        G[NC + 1, NN + 1] = -1
        h = np.zeros(NC + 2)

        G = cvxopt.matrix(np.ascontiguousarray(G))
        h = cvxopt.matrix(np.ascontiguousarray(h))
        for i in np.ndindex(*data.shape[:-1]):
            signal = data[i]
            q = cvxopt.matrix(np.ascontiguousarray(-1 * np.dot(conv_matrix.T,
                                                               signal)))

            sol = cvxopt.solvers.qp(P, q, G, h)
            if sol['status'] != 'optimal':
                logging.debug('Optimization unsuccessful - '
                              'Voxel: {}, Constraint: {}'.format(i, 'nonneg'))

            deconvolution_result[i] = np.array(sol['x'])[:, 0]

        return deconvolution_result

    def sh_convolution_matrix(self, kernel="rank1"):
        """

        :param kernel:
        :return:
        """
        # TODO: Improve this by using esh_matrix
        if self.kernel_type != kernel:
            self.set_kernel(kernel)

        # Build matrix that maps ODF to signal
        M = np.zeros((self._gtab.bvals.shape[0], esh.LENGTH[self.order]))
        r, theta, phi = cart2sphere(self._gtabbvecs[:, 0],
                                    self._gtabbvecs[:, 1],
                                    self._gtabbvecs[:, 2])
        theta[np.isnan(theta)] = 0
        counter = 0
        for l in range(0, self.order + 1, 2):
            for m in range(-l, l + 1):
                M[:, counter] = (
                    real_sph_harm(m, l, theta, phi) * kernel_wm[l / 2])
                counter += 1

        return M



def esh_matrix(order, gtab):
    """ Matrix that evaluates SH coeffs in the given directions

    :param order:
    :param gtab:
    :return:
    """
    bvecs = gtab.bvecs
    r, theta, phi = cart2sphere(bvecs[:, 0], bvecs[:, 1], bvecs[:, 2])
    theta[np.isnan(theta)] = 0
    M = np.zeros((bvecs.shape[0], esh.LENGTH[order]))
    counter = 0
    for l in range(0, order + 1, 2):
        for m in range(-l, l + 1):
            M[:, counter] = real_sph_harm(m, l, theta, phi)
            counter += 1
    return M
