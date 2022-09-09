#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_8859598183063454939);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5756722787744684903);
void car_H_mod_fun(double *state, double *out_1241046553396427050);
void car_f_fun(double *state, double dt, double *out_3265477864722710865);
void car_F_fun(double *state, double dt, double *out_1998648392208618460);
void car_h_25(double *state, double *unused, double *out_7967737535555518081);
void car_H_25(double *state, double *unused, double *out_4817880321206984679);
void car_h_24(double *state, double *unused, double *out_38287872680748422);
void car_H_24(double *state, double *unused, double *out_6990529920212484245);
void car_h_30(double *state, double *unused, double *out_3122496255442127692);
void car_H_30(double *state, double *unused, double *out_2299547362699736052);
void car_h_26(double *state, double *unused, double *out_2091745227628403381);
void car_H_26(double *state, double *unused, double *out_8559383640081040903);
void car_h_27(double *state, double *unused, double *out_6081054456872404523);
void car_H_27(double *state, double *unused, double *out_75953291515792835);
void car_h_29(double *state, double *unused, double *out_4472750121952623734);
void car_H_29(double *state, double *unused, double *out_1789316018385343868);
void car_h_28(double *state, double *unused, double *out_8893606560142257906);
void car_H_28(double *state, double *unused, double *out_6871715035454874442);
void car_h_31(double *state, double *unused, double *out_86464884801468896);
void car_H_31(double *state, double *unused, double *out_9185591742314392379);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}