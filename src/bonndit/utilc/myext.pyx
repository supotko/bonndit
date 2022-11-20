#!python
#cython: language_level=3, boundscheck=False, wraparound=False

cimport cython
import cython as cython
from cython.parallel import prange

import numpy as np
cimport numpy as np
np.import_array()

cdef extern from "Ceres_Cython.h":
        void minimize_watson_mult_o4(double* parameters, double* signal_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* loss_p, int amount, int num_of_dir_p, int no_spread) nogil
        void minimize_watson_mult_o6(double* parameters, double* signal_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* loss_p, int amount, int num_of_dir_p, int no_spread) nogil
        void minimize_watson_mult_o8(double* parameters, double* signal_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* loss_p, int amount, int num_of_dir_p, int no_spread) nogil

cdef void mw_openmp_mult(double[:,:] x, double[:,:] signal, double[:,:] est_signal, double[:,:] dipy_v, double[:,:,:,:] pysh_v, double[:,:,:,:] rot_pysh_v, double[:,:] angles_v, double[:] loss, int amount, int order, int num_of_dir, int no_spread) nogil:
        if order == 4:
                mw_openmp_mult_o4(x,signal,est_signal,dipy_v,pysh_v,rot_pysh_v,angles_v,loss,amount,num_of_dir,no_spread)
        elif order == 6:
                mw_openmp_mult_o8(x,signal,est_signal,dipy_v,pysh_v,rot_pysh_v,angles_v,loss,amount,num_of_dir,no_spread)
        elif order == 8:
                mw_openmp_mult_o8(x,signal,est_signal,dipy_v,pysh_v,rot_pysh_v,angles_v,loss,amount,num_of_dir,no_spread)

cdef void mw_openmp_mult_o4(double[:,:] x, double[:,:] signal, double[:,:] est_signal, double[:,:] dipy_v, double[:,:,:,:] pysh_v, double[:,:,:,:] rot_pysh_v, double[:,:] angles_v, double[:] loss, int amount, int num_of_dir, int no_spread) nogil:
        minimize_watson_mult_o4(&x[0,0],&signal[0,0],&est_signal[0,0],&dipy_v[0,0],&pysh_v[0,0,0,0],&rot_pysh_v[0,0,0,0],&angles_v[0,0],&loss[0],amount,num_of_dir,no_spread)

cdef void mw_openmp_mult_o6(double[:,:] x, double[:,:] signal, double[:,:] est_signal, double[:,:] dipy_v, double[:,:,:,:] pysh_v, double[:,:,:,:] rot_pysh_v, double[:,:] angles_v, double[:] loss, int amount, int num_of_dir, int no_spread) nogil:
        minimize_watson_mult_o6(&x[0,0],&signal[0,0],&est_signal[0,0],&dipy_v[0,0],&pysh_v[0,0,0,0],&rot_pysh_v[0,0,0,0],&angles_v[0,0],&loss[0],amount,num_of_dir,no_spread)

cdef void mw_openmp_mult_o8(double[:,:] x, double[:,:] signal, double[:,:] est_signal, double[:,:] dipy_v, double[:,:,:,:] pysh_v, double[:,:,:,:] rot_pysh_v, double[:,:] angles_v, double[:] loss, int amount, int num_of_dir, int no_spread) nogil:
        minimize_watson_mult_o8(&x[0,0],&signal[0,0],&est_signal[0,0],&dipy_v[0,0],&pysh_v[0,0,0,0],&rot_pysh_v[0,0,0,0],&angles_v[0,0],&loss[0],amount,num_of_dir,no_spread)