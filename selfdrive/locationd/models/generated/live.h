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
void live_H(double *in_vec, double *out_1661342885356962876);
void live_err_fun(double *nom_x, double *delta_x, double *out_2279440072678078976);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_5203322416787072396);
void live_H_mod_fun(double *state, double *out_7445739159746980786);
void live_f_fun(double *state, double dt, double *out_8029348862025951388);
void live_F_fun(double *state, double dt, double *out_4044065417594102198);
void live_h_4(double *state, double *unused, double *out_4601813282528246286);
void live_H_4(double *state, double *unused, double *out_6822193976202631850);
void live_h_9(double *state, double *unused, double *out_2050477823565117837);
void live_H_9(double *state, double *unused, double *out_6581004329573041205);
void live_h_10(double *state, double *unused, double *out_2110006911353143150);
void live_H_10(double *state, double *unused, double *out_2584628456959966932);
void live_h_12(double *state, double *unused, double *out_2288494542137110932);
void live_H_12(double *state, double *unused, double *out_1802737568170670055);
void live_h_31(double *state, double *unused, double *out_3320005927060955619);
void live_H_31(double *state, double *unused, double *out_3455531918830024474);
void live_h_32(double *state, double *unused, double *out_3151615380295529888);
void live_H_32(double *state, double *unused, double *out_4490669804141502383);
void live_h_13(double *state, double *unused, double *out_1003991054242780984);
void live_H_13(double *state, double *unused, double *out_9086786425408962614);
void live_h_14(double *state, double *unused, double *out_2050477823565117837);
void live_H_14(double *state, double *unused, double *out_6581004329573041205);
void live_h_33(double *state, double *unused, double *out_8352926782010900739);
void live_H_33(double *state, double *unused, double *out_304974914191166870);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}