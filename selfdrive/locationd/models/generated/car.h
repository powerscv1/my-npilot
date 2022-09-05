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
void car_err_fun(double *nom_x, double *delta_x, double *out_670054212114508416);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6195146868335722971);
void car_H_mod_fun(double *state, double *out_4468006615074037814);
void car_f_fun(double *state, double dt, double *out_2117109067336928460);
void car_F_fun(double *state, double dt, double *out_6711754359709666424);
void car_h_25(double *state, double *unused, double *out_3703869159224570191);
void car_H_25(double *state, double *unused, double *out_8740659271524641218);
void car_h_24(double *state, double *unused, double *out_6866048557873250388);
void car_H_24(double *state, double *unused, double *out_5182733672722318082);
void car_h_30(double *state, double *unused, double *out_66419630970547841);
void car_H_30(double *state, double *unused, double *out_5178388472057302200);
void car_h_26(double *state, double *unused, double *out_8915549541036292924);
void car_H_26(double *state, double *unused, double *out_5964581483310854174);
void car_h_27(double *state, double *unused, double *out_1297804108436275454);
void car_H_27(double *state, double *unused, double *out_3003625160256877289);
void car_h_29(double *state, double *unused, double *out_4935253636690297804);
void car_H_29(double *state, double *unused, double *out_5688619816371694384);
void car_h_28(double *state, double *unused, double *out_4703326159247582961);
void car_H_28(double *state, double *unused, double *out_7652250087937020635);
void car_h_31(double *state, double *unused, double *out_3756689605024318134);
void car_H_31(double *state, double *unused, double *out_5338373381077502698);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}