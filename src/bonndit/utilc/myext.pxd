#%%cython --annotate
#cython: language_level=3, boundscheck=False, wraparound=False, warn.unused=True, warn.unused_args=True,
# warn.unused_results=True

cdef extern from "Ceres_Cython.h":
        void minimize_watson_mult_o4(double* parameters, double* signal_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* loss_p, int amount, int num_of_dir_p, int no_spread) nogil
        void minimize_watson_mult_o6(double* parameters, double* signal_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* loss_p, int amount, int num_of_dir_p, int no_spread) nogil
        void minimize_watson_mult_o8(double* parameters, double* signal_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* loss_p, int amount, int num_of_dir_p, int no_spread) nogil

cdef void mw_openmp_mult_o4(double[:,:], double[:,:], double[:,:], double[:,:], double[:,:,:,:], double[:,:,:,:], double[:,:], double[:], int, int, int) nogil
cdef void mw_openmp_mult_o6(double[:,:], double[:,:], double[:,:], double[:,:], double[:,:,:,:], double[:,:,:,:], double[:,:], double[:], int, int, int) nogil
cdef void mw_openmp_mult_o8(double[:,:], double[:,:], double[:,:], double[:,:], double[:,:,:,:], double[:,:,:,:], double[:,:], double[:], int, int, int) nogil