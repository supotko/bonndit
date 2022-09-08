#ifndef CERES_CYTHON_H
#define CERES_CYTHON_H

void minimize_watson_mult(double* parameters, double* signal_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* dj_p, double* loss_p, int amount, int lmax_p, int num_of_dir_p);
double* minimize_watson(double parameters[12], double* signal_p, double* modif_x_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* dj_p, int num_of_dir_p, double* loss_p);
double* minimize_watson_ls(double parameters[12], double* signal_p, double* modif_x_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* dj_p, int num_of_dir_p, double* loss_p);
void minimize_watson_multiple(double* mult_paramters, double* signal_p, double* modif_x_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* dj_p, int num_of_dir_p, int num_of_voxels, double* loss_p);
int init_python();
void getdj(double* dj);
void minimize_watson_mult_o4(double* parameters, double* signal_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* dj_p, double* loss_p, int amount, int lmax_p, int num_of_dir_p);
void minimize_watson_single_o4(double* parameters, double* signal_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* loss_p, int num_of_dir_p);

#endif