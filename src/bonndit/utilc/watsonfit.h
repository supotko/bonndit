#ifndef WATSONFIT_H
#define WATSONFIT_H

void minimize_watson_mult_o4(double* parameters, double* signal_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* loss_p, int amount, int num_of_dir_p, int no_spread);
void minimize_watson_mult_o6(double* parameters, double* signal_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* loss_p, int amount, int num_of_dir_p, int no_spread);
void minimize_watson_mult_o8(double* parameters, double* signal_p, double* est_signal_p, double* dipy_v_p, double* pysh_v_p, double* rot_pysh_v_p, double* angles_v_p, double* loss_p, int amount, int num_of_dir_p, int no_spread);

#endif