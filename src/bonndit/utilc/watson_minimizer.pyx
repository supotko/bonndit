#!python
#cython: language_level=3

cimport cython
import cython as cython
from cython.parallel import prange

import numpy as np
cimport numpy as np
np.import_array()

cdef extern from "<Ceres_Cython.h>":
        void minimize_watson_mult(double* parameters, double* signal_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* dj_p, double* loss_p, int amount, int lmax_p, int num_of_dir_p) nogil
        double* minimize_watson(double* parameters, double* signal_p, double* modif_x_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* dj_p, int num_of_dir_p, double* loss_p) nogil
        double* minimize_watson_ls(double* parameters, double* signal_p, double* modif_x_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* dj_p, int num_of_dir_p, double* loss_p) nogil
        #void minimize_watson_multiple(double* mult_paramters, double* signal_p, double* modif_x_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* dj_p, int num_of_dir_p, int num_of_voxels, double* loss_p) nogil
        #void getSHRotateRealCoef(double* pysh_v_c, double* rot_pysh_v_c, double* angles_v_c, double* dj_c, int lmax)
        int init_python()
        #void getdj(double* dj)

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef int init_py():
        return init_python()

# @cython.boundscheck(False)  # Deactivate bounds checking
# @cython.wraparound(False)   # Deactivate negative indexing.
# def mw(double[:] x, double[:] signal, double[:] modif_x, double[:] est_signal, double[:,:] dipy_v, double[:,:,:,:] pysh_v, double[:,:,:,:] rot_pysh_v, double[:] angles_v, double[:,:,:] dj, double[:] loss):
#         cdef double* result
#         cdef double[:] result_memview
#         with nogil:
#                 result = test.minimize_watson(&x[0],&signal[0],&modif_x[0],&est_signal[0],&dipy_v[0,0],&pysh_v[0,0,0,0],&rot_pysh_v[0,0,0,0],&angles_v[0],&dj[0,0,0],3, &loss[0])
#         result_memview = <double[:12]>result
#         return np.asarray(result_memview)

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def mw_openmp_mult(double[:,:] x, double[:,:] signal, double[:,:] est_signal, double[:,:] dipy_v, double[:,:,:,:] pysh_v, double[:,:,:,:] rot_pysh_v, double[:,:] angles_v, double[:,:,:] dj, double[:] loss, int amount, int lmax, int num_of_dir):
        with nogil:
                minimize_watson_mult(&x[0,0],&signal[0,0],&est_signal[0,0],&dipy_v[0,0],&pysh_v[0,0,0,0],&rot_pysh_v[0,0,0,0],&angles_v[0,0],&dj[0,0,0],&loss[0],amount,lmax,num_of_dir)

# @cython.boundscheck(False)  # Deactivate bounds checking
# @cython.wraparound(False)   # Deactivate negative indexing.
# def mw_ls(double[:] x, double[:] signal, double[:] modif_x, double[:] est_signal, double[:,:] dipy_v, double[:,:,:,:] pysh_v, double[:,:,:,:] rot_pysh_v, double[:] angles_v, double[:,:,:] dj, double[:] loss):
#         cdef double* result
#         cdef double[:] result_memview
#         with nogil:
#                 result = test.minimize_watson_ls(&x[0],&signal[0],&modif_x[0],&est_signal[0],&dipy_v[0,0],&pysh_v[0,0,0,0],&rot_pysh_v[0,0,0,0],&angles_v[0],&dj[0,0,0],3, &loss[0])
#         result_memview = <double[:12]>result
#         return np.asarray(result_memview)

# @cython.boundscheck(False)  # Deactivate bounds checking 
# @cython.wraparound(False)   # Deactivate negative indexing.
# cdef void c_mw_multiple(double[:,:] x, double[:,:] signal, double[:,:] modif_x, double[:,:] est_signal, double[:,:] dipy_v, double[:,:,:,:] pysh_v, double[:,:,:,:] rot_pysh_v, double[:,:] angles_v, double[:,:,:] dj, double[:,:] loss) nogil:
#         cdef int i, k = x.shape[0]

#         for i in range(k):
#                 test.minimize_watson(&x[i,0],&signal[i,0],&modif_x[i,0],&est_signal[i,0],&dipy_v[i,0],&pysh_v[i,0,0,0],&rot_pysh_v[i,0,0,0],&angles_v[i,0],&dj[0,0,0],3, &loss[i,0])

# @cython.boundscheck(False)  # Deactivate bounds checking
# @cython.wraparound(False)   # Deactivate negative indexing.
# def mw_multiple(double[:,:] x, double[:,:] signal, double[:,:] modif_x, double[:,:] est_signal, double[:,:] dipy_v, double[:,:,:,:] pysh_v, double[:,:,:,:] rot_pysh_v, double[:,:] angles_v, double[:,:,:] dj, double[:,:] loss):
#         c_mw_multiple(x, signal, modif_x, est_signal, dipy_v, pysh_v, rot_pysh_v, angles_v, dj, loss)

# @cython.boundscheck(False)  # Deactivate bounds checking
# @cython.wraparound(False)   # Deactivate negative indexing.
# cdef void c_mw_multiple_t(double[:,:] x, double[:,:] signal, double[:,:] modif_x, double[:,:] est_signal, double[:,:,:] dipy_v, double[:,:,:,:,:] pysh_v, double[:,:,:,:] dj, double[:,:] loss) nogil:
#         cdef int i, k = x.shape[0]

#         for i in prange(k, nogil=False):
#                 test.minimize_watson(&x[i,0],&signal[i,0],&modif_x[i,0],&est_signal[i,0],&dipy_v[i,0,0],&pysh_v[i,0,0,0,0],&dj[i,0,0,0],3, &loss[i,0])

# @cython.boundscheck(False)  # Deactivate bounds checking
# @cython.wraparound(False)   # Deactivate negative indexing.
# def mw_multiple_t(double[:,:] x, double[:,:] signal, double[:,:] modif_x, double[:,:] est_signal, double[:,:,:] dipy_v, double[:,:,:,:,:] pysh_v, double[:,:,:,:] dj, double[:,:] loss):
#         c_mw_multiple_t(x, signal, modif_x, est_signal, dipy_v, pysh_v, dj, loss)

# @cython.boundscheck(False)  # Deactivate bounds checking
# @cython.wraparound(False)   # Deactivate negative indexing.
# cdef void c_mw_multiple_gil(double[:,:] x, double[:,:] signal, double[:] modif_x, double[:] est_signal, double[:,:] dipy_v, double[:,:,:,:] pysh_v, double[:,:,:,:] rot_pysh_v, double[:] angles_v, double[:,:,:] dj, double[:,:] loss):
#         cdef int i, k = x.shape[0]

#         for i in range(k):
#                 test.minimize_watson(&x[i,0],&signal[i,0],&modif_x[0],&est_signal[0],&dipy_v[0,0],&pysh_v[0,0,0,0],&rot_pysh_v[0,0,0,0],&angles_v[0],&dj[0,0,0],3, &loss[i,0])

# @cython.boundscheck(False)  # Deactivate bounds checking
# @cython.wraparound(False)   # Deactivate negative indexing.
# def mw_multiple_gil(double[:,:] x, double[:,:] signal, double[:] modif_x, double[:] est_signal, double[:,:] dipy_v, double[:,:,:,:] pysh_v, double[:,:,:] dj, double[:,:] loss):
#         c_mw_multiple_gil(x, signal, modif_x, est_signal, dipy_v, pysh_v, dj, loss)




# @cython.boundscheck(False)  # Deactivate bounds checking
# @cython.wraparound(False)   # Deactivate negative indexing.
# cdef void c_native_mw_multiple(double[:,:] x, double[:,:] signal, double[:] modif_x, double[:] est_signal, double[:,:] dipy_v, double[:,:,:,:] pysh_v, double[:,:,:] dj, int num_of_voxels, double[:,:] loss) nogil:
#         test.minimize_watson_multiple(&x[0,0],&signal[0,0],&modif_x[0],&est_signal[0],&dipy_v[0,0],&pysh_v[0,0,0,0],&dj[0,0,0],3,num_of_voxels, &loss[0,0])

# @cython.boundscheck(False)  # Deactivate bounds checking
# @cython.wraparound(False)   # Deactivate negative indexing.
# def mw_multiple_cpp(double[:,:] x, double[:,:] signal, double[:] modif_x, double[:] est_signal, double[:,:] dipy_v, double[:,:,:,:] pysh_v, double[:,:,:] dj, int num_of_voxels, double[:,:] loss):
#         c_native_mw_multiple(x, signal, modif_x, est_signal, dipy_v, pysh_v, dj, num_of_voxels, loss)

# def getSHRotateRealCoef(double[:,:,:] pysh_v_c, double[:,:,:] rot_pysh_v_c, double[:] angles_v_c, double[:,:,:] dj_c, int lmax):
#         test.getSHRotateRealCoef(&pysh_v_c[0,0,0], &rot_pysh_v_c[0,0,0], &angles_v_c[0], &dj_c[0,0,0], lmax)

# def getdj(double[:,:,:] dj_c):
#         getdj(&dj_c[0,0,0])