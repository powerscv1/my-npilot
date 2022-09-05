#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_670054212114508416) {
   out_670054212114508416[0] = delta_x[0] + nom_x[0];
   out_670054212114508416[1] = delta_x[1] + nom_x[1];
   out_670054212114508416[2] = delta_x[2] + nom_x[2];
   out_670054212114508416[3] = delta_x[3] + nom_x[3];
   out_670054212114508416[4] = delta_x[4] + nom_x[4];
   out_670054212114508416[5] = delta_x[5] + nom_x[5];
   out_670054212114508416[6] = delta_x[6] + nom_x[6];
   out_670054212114508416[7] = delta_x[7] + nom_x[7];
   out_670054212114508416[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_6195146868335722971) {
   out_6195146868335722971[0] = -nom_x[0] + true_x[0];
   out_6195146868335722971[1] = -nom_x[1] + true_x[1];
   out_6195146868335722971[2] = -nom_x[2] + true_x[2];
   out_6195146868335722971[3] = -nom_x[3] + true_x[3];
   out_6195146868335722971[4] = -nom_x[4] + true_x[4];
   out_6195146868335722971[5] = -nom_x[5] + true_x[5];
   out_6195146868335722971[6] = -nom_x[6] + true_x[6];
   out_6195146868335722971[7] = -nom_x[7] + true_x[7];
   out_6195146868335722971[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_4468006615074037814) {
   out_4468006615074037814[0] = 1.0;
   out_4468006615074037814[1] = 0;
   out_4468006615074037814[2] = 0;
   out_4468006615074037814[3] = 0;
   out_4468006615074037814[4] = 0;
   out_4468006615074037814[5] = 0;
   out_4468006615074037814[6] = 0;
   out_4468006615074037814[7] = 0;
   out_4468006615074037814[8] = 0;
   out_4468006615074037814[9] = 0;
   out_4468006615074037814[10] = 1.0;
   out_4468006615074037814[11] = 0;
   out_4468006615074037814[12] = 0;
   out_4468006615074037814[13] = 0;
   out_4468006615074037814[14] = 0;
   out_4468006615074037814[15] = 0;
   out_4468006615074037814[16] = 0;
   out_4468006615074037814[17] = 0;
   out_4468006615074037814[18] = 0;
   out_4468006615074037814[19] = 0;
   out_4468006615074037814[20] = 1.0;
   out_4468006615074037814[21] = 0;
   out_4468006615074037814[22] = 0;
   out_4468006615074037814[23] = 0;
   out_4468006615074037814[24] = 0;
   out_4468006615074037814[25] = 0;
   out_4468006615074037814[26] = 0;
   out_4468006615074037814[27] = 0;
   out_4468006615074037814[28] = 0;
   out_4468006615074037814[29] = 0;
   out_4468006615074037814[30] = 1.0;
   out_4468006615074037814[31] = 0;
   out_4468006615074037814[32] = 0;
   out_4468006615074037814[33] = 0;
   out_4468006615074037814[34] = 0;
   out_4468006615074037814[35] = 0;
   out_4468006615074037814[36] = 0;
   out_4468006615074037814[37] = 0;
   out_4468006615074037814[38] = 0;
   out_4468006615074037814[39] = 0;
   out_4468006615074037814[40] = 1.0;
   out_4468006615074037814[41] = 0;
   out_4468006615074037814[42] = 0;
   out_4468006615074037814[43] = 0;
   out_4468006615074037814[44] = 0;
   out_4468006615074037814[45] = 0;
   out_4468006615074037814[46] = 0;
   out_4468006615074037814[47] = 0;
   out_4468006615074037814[48] = 0;
   out_4468006615074037814[49] = 0;
   out_4468006615074037814[50] = 1.0;
   out_4468006615074037814[51] = 0;
   out_4468006615074037814[52] = 0;
   out_4468006615074037814[53] = 0;
   out_4468006615074037814[54] = 0;
   out_4468006615074037814[55] = 0;
   out_4468006615074037814[56] = 0;
   out_4468006615074037814[57] = 0;
   out_4468006615074037814[58] = 0;
   out_4468006615074037814[59] = 0;
   out_4468006615074037814[60] = 1.0;
   out_4468006615074037814[61] = 0;
   out_4468006615074037814[62] = 0;
   out_4468006615074037814[63] = 0;
   out_4468006615074037814[64] = 0;
   out_4468006615074037814[65] = 0;
   out_4468006615074037814[66] = 0;
   out_4468006615074037814[67] = 0;
   out_4468006615074037814[68] = 0;
   out_4468006615074037814[69] = 0;
   out_4468006615074037814[70] = 1.0;
   out_4468006615074037814[71] = 0;
   out_4468006615074037814[72] = 0;
   out_4468006615074037814[73] = 0;
   out_4468006615074037814[74] = 0;
   out_4468006615074037814[75] = 0;
   out_4468006615074037814[76] = 0;
   out_4468006615074037814[77] = 0;
   out_4468006615074037814[78] = 0;
   out_4468006615074037814[79] = 0;
   out_4468006615074037814[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_2117109067336928460) {
   out_2117109067336928460[0] = state[0];
   out_2117109067336928460[1] = state[1];
   out_2117109067336928460[2] = state[2];
   out_2117109067336928460[3] = state[3];
   out_2117109067336928460[4] = state[4];
   out_2117109067336928460[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2117109067336928460[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2117109067336928460[7] = state[7];
   out_2117109067336928460[8] = state[8];
}
void F_fun(double *state, double dt, double *out_6711754359709666424) {
   out_6711754359709666424[0] = 1;
   out_6711754359709666424[1] = 0;
   out_6711754359709666424[2] = 0;
   out_6711754359709666424[3] = 0;
   out_6711754359709666424[4] = 0;
   out_6711754359709666424[5] = 0;
   out_6711754359709666424[6] = 0;
   out_6711754359709666424[7] = 0;
   out_6711754359709666424[8] = 0;
   out_6711754359709666424[9] = 0;
   out_6711754359709666424[10] = 1;
   out_6711754359709666424[11] = 0;
   out_6711754359709666424[12] = 0;
   out_6711754359709666424[13] = 0;
   out_6711754359709666424[14] = 0;
   out_6711754359709666424[15] = 0;
   out_6711754359709666424[16] = 0;
   out_6711754359709666424[17] = 0;
   out_6711754359709666424[18] = 0;
   out_6711754359709666424[19] = 0;
   out_6711754359709666424[20] = 1;
   out_6711754359709666424[21] = 0;
   out_6711754359709666424[22] = 0;
   out_6711754359709666424[23] = 0;
   out_6711754359709666424[24] = 0;
   out_6711754359709666424[25] = 0;
   out_6711754359709666424[26] = 0;
   out_6711754359709666424[27] = 0;
   out_6711754359709666424[28] = 0;
   out_6711754359709666424[29] = 0;
   out_6711754359709666424[30] = 1;
   out_6711754359709666424[31] = 0;
   out_6711754359709666424[32] = 0;
   out_6711754359709666424[33] = 0;
   out_6711754359709666424[34] = 0;
   out_6711754359709666424[35] = 0;
   out_6711754359709666424[36] = 0;
   out_6711754359709666424[37] = 0;
   out_6711754359709666424[38] = 0;
   out_6711754359709666424[39] = 0;
   out_6711754359709666424[40] = 1;
   out_6711754359709666424[41] = 0;
   out_6711754359709666424[42] = 0;
   out_6711754359709666424[43] = 0;
   out_6711754359709666424[44] = 0;
   out_6711754359709666424[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_6711754359709666424[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_6711754359709666424[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6711754359709666424[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6711754359709666424[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_6711754359709666424[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_6711754359709666424[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_6711754359709666424[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_6711754359709666424[53] = -9.8000000000000007*dt;
   out_6711754359709666424[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_6711754359709666424[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_6711754359709666424[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6711754359709666424[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6711754359709666424[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_6711754359709666424[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_6711754359709666424[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_6711754359709666424[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6711754359709666424[62] = 0;
   out_6711754359709666424[63] = 0;
   out_6711754359709666424[64] = 0;
   out_6711754359709666424[65] = 0;
   out_6711754359709666424[66] = 0;
   out_6711754359709666424[67] = 0;
   out_6711754359709666424[68] = 0;
   out_6711754359709666424[69] = 0;
   out_6711754359709666424[70] = 1;
   out_6711754359709666424[71] = 0;
   out_6711754359709666424[72] = 0;
   out_6711754359709666424[73] = 0;
   out_6711754359709666424[74] = 0;
   out_6711754359709666424[75] = 0;
   out_6711754359709666424[76] = 0;
   out_6711754359709666424[77] = 0;
   out_6711754359709666424[78] = 0;
   out_6711754359709666424[79] = 0;
   out_6711754359709666424[80] = 1;
}
void h_25(double *state, double *unused, double *out_3703869159224570191) {
   out_3703869159224570191[0] = state[6];
}
void H_25(double *state, double *unused, double *out_8740659271524641218) {
   out_8740659271524641218[0] = 0;
   out_8740659271524641218[1] = 0;
   out_8740659271524641218[2] = 0;
   out_8740659271524641218[3] = 0;
   out_8740659271524641218[4] = 0;
   out_8740659271524641218[5] = 0;
   out_8740659271524641218[6] = 1;
   out_8740659271524641218[7] = 0;
   out_8740659271524641218[8] = 0;
}
void h_24(double *state, double *unused, double *out_6866048557873250388) {
   out_6866048557873250388[0] = state[4];
   out_6866048557873250388[1] = state[5];
}
void H_24(double *state, double *unused, double *out_5182733672722318082) {
   out_5182733672722318082[0] = 0;
   out_5182733672722318082[1] = 0;
   out_5182733672722318082[2] = 0;
   out_5182733672722318082[3] = 0;
   out_5182733672722318082[4] = 1;
   out_5182733672722318082[5] = 0;
   out_5182733672722318082[6] = 0;
   out_5182733672722318082[7] = 0;
   out_5182733672722318082[8] = 0;
   out_5182733672722318082[9] = 0;
   out_5182733672722318082[10] = 0;
   out_5182733672722318082[11] = 0;
   out_5182733672722318082[12] = 0;
   out_5182733672722318082[13] = 0;
   out_5182733672722318082[14] = 1;
   out_5182733672722318082[15] = 0;
   out_5182733672722318082[16] = 0;
   out_5182733672722318082[17] = 0;
}
void h_30(double *state, double *unused, double *out_66419630970547841) {
   out_66419630970547841[0] = state[4];
}
void H_30(double *state, double *unused, double *out_5178388472057302200) {
   out_5178388472057302200[0] = 0;
   out_5178388472057302200[1] = 0;
   out_5178388472057302200[2] = 0;
   out_5178388472057302200[3] = 0;
   out_5178388472057302200[4] = 1;
   out_5178388472057302200[5] = 0;
   out_5178388472057302200[6] = 0;
   out_5178388472057302200[7] = 0;
   out_5178388472057302200[8] = 0;
}
void h_26(double *state, double *unused, double *out_8915549541036292924) {
   out_8915549541036292924[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5964581483310854174) {
   out_5964581483310854174[0] = 0;
   out_5964581483310854174[1] = 0;
   out_5964581483310854174[2] = 0;
   out_5964581483310854174[3] = 0;
   out_5964581483310854174[4] = 0;
   out_5964581483310854174[5] = 0;
   out_5964581483310854174[6] = 0;
   out_5964581483310854174[7] = 1;
   out_5964581483310854174[8] = 0;
}
void h_27(double *state, double *unused, double *out_1297804108436275454) {
   out_1297804108436275454[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3003625160256877289) {
   out_3003625160256877289[0] = 0;
   out_3003625160256877289[1] = 0;
   out_3003625160256877289[2] = 0;
   out_3003625160256877289[3] = 1;
   out_3003625160256877289[4] = 0;
   out_3003625160256877289[5] = 0;
   out_3003625160256877289[6] = 0;
   out_3003625160256877289[7] = 0;
   out_3003625160256877289[8] = 0;
}
void h_29(double *state, double *unused, double *out_4935253636690297804) {
   out_4935253636690297804[0] = state[1];
}
void H_29(double *state, double *unused, double *out_5688619816371694384) {
   out_5688619816371694384[0] = 0;
   out_5688619816371694384[1] = 1;
   out_5688619816371694384[2] = 0;
   out_5688619816371694384[3] = 0;
   out_5688619816371694384[4] = 0;
   out_5688619816371694384[5] = 0;
   out_5688619816371694384[6] = 0;
   out_5688619816371694384[7] = 0;
   out_5688619816371694384[8] = 0;
}
void h_28(double *state, double *unused, double *out_4703326159247582961) {
   out_4703326159247582961[0] = state[0];
}
void H_28(double *state, double *unused, double *out_7652250087937020635) {
   out_7652250087937020635[0] = 1;
   out_7652250087937020635[1] = 0;
   out_7652250087937020635[2] = 0;
   out_7652250087937020635[3] = 0;
   out_7652250087937020635[4] = 0;
   out_7652250087937020635[5] = 0;
   out_7652250087937020635[6] = 0;
   out_7652250087937020635[7] = 0;
   out_7652250087937020635[8] = 0;
}
void h_31(double *state, double *unused, double *out_3756689605024318134) {
   out_3756689605024318134[0] = state[8];
}
void H_31(double *state, double *unused, double *out_5338373381077502698) {
   out_5338373381077502698[0] = 0;
   out_5338373381077502698[1] = 0;
   out_5338373381077502698[2] = 0;
   out_5338373381077502698[3] = 0;
   out_5338373381077502698[4] = 0;
   out_5338373381077502698[5] = 0;
   out_5338373381077502698[6] = 0;
   out_5338373381077502698[7] = 0;
   out_5338373381077502698[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_670054212114508416) {
  err_fun(nom_x, delta_x, out_670054212114508416);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6195146868335722971) {
  inv_err_fun(nom_x, true_x, out_6195146868335722971);
}
void car_H_mod_fun(double *state, double *out_4468006615074037814) {
  H_mod_fun(state, out_4468006615074037814);
}
void car_f_fun(double *state, double dt, double *out_2117109067336928460) {
  f_fun(state,  dt, out_2117109067336928460);
}
void car_F_fun(double *state, double dt, double *out_6711754359709666424) {
  F_fun(state,  dt, out_6711754359709666424);
}
void car_h_25(double *state, double *unused, double *out_3703869159224570191) {
  h_25(state, unused, out_3703869159224570191);
}
void car_H_25(double *state, double *unused, double *out_8740659271524641218) {
  H_25(state, unused, out_8740659271524641218);
}
void car_h_24(double *state, double *unused, double *out_6866048557873250388) {
  h_24(state, unused, out_6866048557873250388);
}
void car_H_24(double *state, double *unused, double *out_5182733672722318082) {
  H_24(state, unused, out_5182733672722318082);
}
void car_h_30(double *state, double *unused, double *out_66419630970547841) {
  h_30(state, unused, out_66419630970547841);
}
void car_H_30(double *state, double *unused, double *out_5178388472057302200) {
  H_30(state, unused, out_5178388472057302200);
}
void car_h_26(double *state, double *unused, double *out_8915549541036292924) {
  h_26(state, unused, out_8915549541036292924);
}
void car_H_26(double *state, double *unused, double *out_5964581483310854174) {
  H_26(state, unused, out_5964581483310854174);
}
void car_h_27(double *state, double *unused, double *out_1297804108436275454) {
  h_27(state, unused, out_1297804108436275454);
}
void car_H_27(double *state, double *unused, double *out_3003625160256877289) {
  H_27(state, unused, out_3003625160256877289);
}
void car_h_29(double *state, double *unused, double *out_4935253636690297804) {
  h_29(state, unused, out_4935253636690297804);
}
void car_H_29(double *state, double *unused, double *out_5688619816371694384) {
  H_29(state, unused, out_5688619816371694384);
}
void car_h_28(double *state, double *unused, double *out_4703326159247582961) {
  h_28(state, unused, out_4703326159247582961);
}
void car_H_28(double *state, double *unused, double *out_7652250087937020635) {
  H_28(state, unused, out_7652250087937020635);
}
void car_h_31(double *state, double *unused, double *out_3756689605024318134) {
  h_31(state, unused, out_3756689605024318134);
}
void car_H_31(double *state, double *unused, double *out_5338373381077502698) {
  H_31(state, unused, out_5338373381077502698);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
