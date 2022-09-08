cdef extern from "Ceres_Cython.h":
        void minimize_watson_mult(double* parameters, double* signal_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* dj_p, double* loss_p, int amount, int lmax_p, int num_of_dir_p) nogil
        double* minimize_watson(double* parameters, double* signal_p, double* modif_x_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* dj_p, int num_of_dir_p, double* loss_p) nogil
        double* minimize_watson_ls(double* parameters, double* signal_p, double* modif_x_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* dj_p, int num_of_dir_p, double* loss_p) nogil
        #void minimize_watson_multiple(double* mult_paramters, double* signal_p, double* modif_x_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* dj_p, int num_of_dir_p, int num_of_voxels, double* loss_p) nogil
        void getSHRotateRealCoef(double* pysh_v_c, double* rot_pysh_v_c, double* angles_v_c, double* dj_c, int lmax)
        int init_python()
        void getdj(double* dj)