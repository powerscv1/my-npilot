#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_2337975499004216015);
void live_err_fun(double *nom_x, double *delta_x, double *out_7211975948903307689);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_1869955993172153366);
void live_H_mod_fun(double *state, double *out_7005854625022851394);
void live_f_fun(double *state, double dt, double *out_5267891945337088690);
void live_F_fun(double *state, double dt, double *out_2581783222534495222);
void live_h_4(double *state, double *unused, double *out_6899066197363664588);
void live_H_4(double *state, double *unused, double *out_7659050984387611860);
void live_h_9(double *state, double *unused, double *out_1671602116294396863);
void live_H_9(double *state, double *unused, double *out_371832049123164390);
void live_h_10(double *state, double *unused, double *out_119849606827224068);
void live_H_10(double *state, double *unused, double *out_8976141317441350647);
void live_h_12(double *state, double *unused, double *out_6955664201105629571);
void live_H_12(double *state, double *unused, double *out_2639594576355650065);
void live_h_31(double *state, double *unused, double *out_105963961587279508);
void live_H_31(double *state, double *unused, double *out_4292388927015004484);
void live_h_32(double *state, double *unused, double *out_5633413682939672987);
void live_H_32(double *state, double *unused, double *out_7178870349863458980);
void live_h_13(double *state, double *unused, double *out_102321889000685665);
void live_H_13(double *state, double *unused, double *out_2030752425708420177);
void live_h_14(double *state, double *unused, double *out_1671602116294396863);
void live_H_14(double *state, double *unused, double *out_371832049123164390);
void live_h_33(double *state, double *unused, double *out_1238518388196214819);
void live_H_33(double *state, double *unused, double *out_1141831922376146880);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}