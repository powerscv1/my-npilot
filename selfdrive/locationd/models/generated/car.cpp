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
void err_fun(double *nom_x, double *delta_x, double *out_8859598183063454939) {
   out_8859598183063454939[0] = delta_x[0] + nom_x[0];
   out_8859598183063454939[1] = delta_x[1] + nom_x[1];
   out_8859598183063454939[2] = delta_x[2] + nom_x[2];
   out_8859598183063454939[3] = delta_x[3] + nom_x[3];
   out_8859598183063454939[4] = delta_x[4] + nom_x[4];
   out_8859598183063454939[5] = delta_x[5] + nom_x[5];
   out_8859598183063454939[6] = delta_x[6] + nom_x[6];
   out_8859598183063454939[7] = delta_x[7] + nom_x[7];
   out_8859598183063454939[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5756722787744684903) {
   out_5756722787744684903[0] = -nom_x[0] + true_x[0];
   out_5756722787744684903[1] = -nom_x[1] + true_x[1];
   out_5756722787744684903[2] = -nom_x[2] + true_x[2];
   out_5756722787744684903[3] = -nom_x[3] + true_x[3];
   out_5756722787744684903[4] = -nom_x[4] + true_x[4];
   out_5756722787744684903[5] = -nom_x[5] + true_x[5];
   out_5756722787744684903[6] = -nom_x[6] + true_x[6];
   out_5756722787744684903[7] = -nom_x[7] + true_x[7];
   out_5756722787744684903[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_1241046553396427050) {
   out_1241046553396427050[0] = 1.0;
   out_1241046553396427050[1] = 0;
   out_1241046553396427050[2] = 0;
   out_1241046553396427050[3] = 0;
   out_1241046553396427050[4] = 0;
   out_1241046553396427050[5] = 0;
   out_1241046553396427050[6] = 0;
   out_1241046553396427050[7] = 0;
   out_1241046553396427050[8] = 0;
   out_1241046553396427050[9] = 0;
   out_1241046553396427050[10] = 1.0;
   out_1241046553396427050[11] = 0;
   out_1241046553396427050[12] = 0;
   out_1241046553396427050[13] = 0;
   out_1241046553396427050[14] = 0;
   out_1241046553396427050[15] = 0;
   out_1241046553396427050[16] = 0;
   out_1241046553396427050[17] = 0;
   out_1241046553396427050[18] = 0;
   out_1241046553396427050[19] = 0;
   out_1241046553396427050[20] = 1.0;
   out_1241046553396427050[21] = 0;
   out_1241046553396427050[22] = 0;
   out_1241046553396427050[23] = 0;
   out_1241046553396427050[24] = 0;
   out_1241046553396427050[25] = 0;
   out_1241046553396427050[26] = 0;
   out_1241046553396427050[27] = 0;
   out_1241046553396427050[28] = 0;
   out_1241046553396427050[29] = 0;
   out_1241046553396427050[30] = 1.0;
   out_1241046553396427050[31] = 0;
   out_1241046553396427050[32] = 0;
   out_1241046553396427050[33] = 0;
   out_1241046553396427050[34] = 0;
   out_1241046553396427050[35] = 0;
   out_1241046553396427050[36] = 0;
   out_1241046553396427050[37] = 0;
   out_1241046553396427050[38] = 0;
   out_1241046553396427050[39] = 0;
   out_1241046553396427050[40] = 1.0;
   out_1241046553396427050[41] = 0;
   out_1241046553396427050[42] = 0;
   out_1241046553396427050[43] = 0;
   out_1241046553396427050[44] = 0;
   out_1241046553396427050[45] = 0;
   out_1241046553396427050[46] = 0;
   out_1241046553396427050[47] = 0;
   out_1241046553396427050[48] = 0;
   out_1241046553396427050[49] = 0;
   out_1241046553396427050[50] = 1.0;
   out_1241046553396427050[51] = 0;
   out_1241046553396427050[52] = 0;
   out_1241046553396427050[53] = 0;
   out_1241046553396427050[54] = 0;
   out_1241046553396427050[55] = 0;
   out_1241046553396427050[56] = 0;
   out_1241046553396427050[57] = 0;
   out_1241046553396427050[58] = 0;
   out_1241046553396427050[59] = 0;
   out_1241046553396427050[60] = 1.0;
   out_1241046553396427050[61] = 0;
   out_1241046553396427050[62] = 0;
   out_1241046553396427050[63] = 0;
   out_1241046553396427050[64] = 0;
   out_1241046553396427050[65] = 0;
   out_1241046553396427050[66] = 0;
   out_1241046553396427050[67] = 0;
   out_1241046553396427050[68] = 0;
   out_1241046553396427050[69] = 0;
   out_1241046553396427050[70] = 1.0;
   out_1241046553396427050[71] = 0;
   out_1241046553396427050[72] = 0;
   out_1241046553396427050[73] = 0;
   out_1241046553396427050[74] = 0;
   out_1241046553396427050[75] = 0;
   out_1241046553396427050[76] = 0;
   out_1241046553396427050[77] = 0;
   out_1241046553396427050[78] = 0;
   out_1241046553396427050[79] = 0;
   out_1241046553396427050[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_3265477864722710865) {
   out_3265477864722710865[0] = state[0];
   out_3265477864722710865[1] = state[1];
   out_3265477864722710865[2] = state[2];
   out_3265477864722710865[3] = state[3];
   out_3265477864722710865[4] = state[4];
   out_3265477864722710865[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_3265477864722710865[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_3265477864722710865[7] = state[7];
   out_3265477864722710865[8] = state[8];
}
void F_fun(double *state, double dt, double *out_1998648392208618460) {
   out_1998648392208618460[0] = 1;
   out_1998648392208618460[1] = 0;
   out_1998648392208618460[2] = 0;
   out_1998648392208618460[3] = 0;
   out_1998648392208618460[4] = 0;
   out_1998648392208618460[5] = 0;
   out_1998648392208618460[6] = 0;
   out_1998648392208618460[7] = 0;
   out_1998648392208618460[8] = 0;
   out_1998648392208618460[9] = 0;
   out_1998648392208618460[10] = 1;
   out_1998648392208618460[11] = 0;
   out_1998648392208618460[12] = 0;
   out_1998648392208618460[13] = 0;
   out_1998648392208618460[14] = 0;
   out_1998648392208618460[15] = 0;
   out_1998648392208618460[16] = 0;
   out_1998648392208618460[17] = 0;
   out_1998648392208618460[18] = 0;
   out_1998648392208618460[19] = 0;
   out_1998648392208618460[20] = 1;
   out_1998648392208618460[21] = 0;
   out_1998648392208618460[22] = 0;
   out_1998648392208618460[23] = 0;
   out_1998648392208618460[24] = 0;
   out_1998648392208618460[25] = 0;
   out_1998648392208618460[26] = 0;
   out_1998648392208618460[27] = 0;
   out_1998648392208618460[28] = 0;
   out_1998648392208618460[29] = 0;
   out_1998648392208618460[30] = 1;
   out_1998648392208618460[31] = 0;
   out_1998648392208618460[32] = 0;
   out_1998648392208618460[33] = 0;
   out_1998648392208618460[34] = 0;
   out_1998648392208618460[35] = 0;
   out_1998648392208618460[36] = 0;
   out_1998648392208618460[37] = 0;
   out_1998648392208618460[38] = 0;
   out_1998648392208618460[39] = 0;
   out_1998648392208618460[40] = 1;
   out_1998648392208618460[41] = 0;
   out_1998648392208618460[42] = 0;
   out_1998648392208618460[43] = 0;
   out_1998648392208618460[44] = 0;
   out_1998648392208618460[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_1998648392208618460[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_1998648392208618460[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1998648392208618460[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1998648392208618460[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_1998648392208618460[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_1998648392208618460[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_1998648392208618460[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_1998648392208618460[53] = -9.8000000000000007*dt;
   out_1998648392208618460[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_1998648392208618460[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_1998648392208618460[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1998648392208618460[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1998648392208618460[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_1998648392208618460[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_1998648392208618460[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_1998648392208618460[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1998648392208618460[62] = 0;
   out_1998648392208618460[63] = 0;
   out_1998648392208618460[64] = 0;
   out_1998648392208618460[65] = 0;
   out_1998648392208618460[66] = 0;
   out_1998648392208618460[67] = 0;
   out_1998648392208618460[68] = 0;
   out_1998648392208618460[69] = 0;
   out_1998648392208618460[70] = 1;
   out_1998648392208618460[71] = 0;
   out_1998648392208618460[72] = 0;
   out_1998648392208618460[73] = 0;
   out_1998648392208618460[74] = 0;
   out_1998648392208618460[75] = 0;
   out_1998648392208618460[76] = 0;
   out_1998648392208618460[77] = 0;
   out_1998648392208618460[78] = 0;
   out_1998648392208618460[79] = 0;
   out_1998648392208618460[80] = 1;
}
void h_25(double *state, double *unused, double *out_7967737535555518081) {
   out_7967737535555518081[0] = state[6];
}
void H_25(double *state, double *unused, double *out_4817880321206984679) {
   out_4817880321206984679[0] = 0;
   out_4817880321206984679[1] = 0;
   out_4817880321206984679[2] = 0;
   out_4817880321206984679[3] = 0;
   out_4817880321206984679[4] = 0;
   out_4817880321206984679[5] = 0;
   out_4817880321206984679[6] = 1;
   out_4817880321206984679[7] = 0;
   out_4817880321206984679[8] = 0;
}
void h_24(double *state, double *unused, double *out_38287872680748422) {
   out_38287872680748422[0] = state[4];
   out_38287872680748422[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6990529920212484245) {
   out_6990529920212484245[0] = 0;
   out_6990529920212484245[1] = 0;
   out_6990529920212484245[2] = 0;
   out_6990529920212484245[3] = 0;
   out_6990529920212484245[4] = 1;
   out_6990529920212484245[5] = 0;
   out_6990529920212484245[6] = 0;
   out_6990529920212484245[7] = 0;
   out_6990529920212484245[8] = 0;
   out_6990529920212484245[9] = 0;
   out_6990529920212484245[10] = 0;
   out_6990529920212484245[11] = 0;
   out_6990529920212484245[12] = 0;
   out_6990529920212484245[13] = 0;
   out_6990529920212484245[14] = 1;
   out_6990529920212484245[15] = 0;
   out_6990529920212484245[16] = 0;
   out_6990529920212484245[17] = 0;
}
void h_30(double *state, double *unused, double *out_3122496255442127692) {
   out_3122496255442127692[0] = state[4];
}
void H_30(double *state, double *unused, double *out_2299547362699736052) {
   out_2299547362699736052[0] = 0;
   out_2299547362699736052[1] = 0;
   out_2299547362699736052[2] = 0;
   out_2299547362699736052[3] = 0;
   out_2299547362699736052[4] = 1;
   out_2299547362699736052[5] = 0;
   out_2299547362699736052[6] = 0;
   out_2299547362699736052[7] = 0;
   out_2299547362699736052[8] = 0;
}
void h_26(double *state, double *unused, double *out_2091745227628403381) {
   out_2091745227628403381[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8559383640081040903) {
   out_8559383640081040903[0] = 0;
   out_8559383640081040903[1] = 0;
   out_8559383640081040903[2] = 0;
   out_8559383640081040903[3] = 0;
   out_8559383640081040903[4] = 0;
   out_8559383640081040903[5] = 0;
   out_8559383640081040903[6] = 0;
   out_8559383640081040903[7] = 1;
   out_8559383640081040903[8] = 0;
}
void h_27(double *state, double *unused, double *out_6081054456872404523) {
   out_6081054456872404523[0] = state[3];
}
void H_27(double *state, double *unused, double *out_75953291515792835) {
   out_75953291515792835[0] = 0;
   out_75953291515792835[1] = 0;
   out_75953291515792835[2] = 0;
   out_75953291515792835[3] = 1;
   out_75953291515792835[4] = 0;
   out_75953291515792835[5] = 0;
   out_75953291515792835[6] = 0;
   out_75953291515792835[7] = 0;
   out_75953291515792835[8] = 0;
}
void h_29(double *state, double *unused, double *out_4472750121952623734) {
   out_4472750121952623734[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1789316018385343868) {
   out_1789316018385343868[0] = 0;
   out_1789316018385343868[1] = 1;
   out_1789316018385343868[2] = 0;
   out_1789316018385343868[3] = 0;
   out_1789316018385343868[4] = 0;
   out_1789316018385343868[5] = 0;
   out_1789316018385343868[6] = 0;
   out_1789316018385343868[7] = 0;
   out_1789316018385343868[8] = 0;
}
void h_28(double *state, double *unused, double *out_8893606560142257906) {
   out_8893606560142257906[0] = state[0];
}
void H_28(double *state, double *unused, double *out_6871715035454874442) {
   out_6871715035454874442[0] = 1;
   out_6871715035454874442[1] = 0;
   out_6871715035454874442[2] = 0;
   out_6871715035454874442[3] = 0;
   out_6871715035454874442[4] = 0;
   out_6871715035454874442[5] = 0;
   out_6871715035454874442[6] = 0;
   out_6871715035454874442[7] = 0;
   out_6871715035454874442[8] = 0;
}
void h_31(double *state, double *unused, double *out_86464884801468896) {
   out_86464884801468896[0] = state[8];
}
void H_31(double *state, double *unused, double *out_9185591742314392379) {
   out_9185591742314392379[0] = 0;
   out_9185591742314392379[1] = 0;
   out_9185591742314392379[2] = 0;
   out_9185591742314392379[3] = 0;
   out_9185591742314392379[4] = 0;
   out_9185591742314392379[5] = 0;
   out_9185591742314392379[6] = 0;
   out_9185591742314392379[7] = 0;
   out_9185591742314392379[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_8859598183063454939) {
  err_fun(nom_x, delta_x, out_8859598183063454939);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5756722787744684903) {
  inv_err_fun(nom_x, true_x, out_5756722787744684903);
}
void car_H_mod_fun(double *state, double *out_1241046553396427050) {
  H_mod_fun(state, out_1241046553396427050);
}
void car_f_fun(double *state, double dt, double *out_3265477864722710865) {
  f_fun(state,  dt, out_3265477864722710865);
}
void car_F_fun(double *state, double dt, double *out_1998648392208618460) {
  F_fun(state,  dt, out_1998648392208618460);
}
void car_h_25(double *state, double *unused, double *out_7967737535555518081) {
  h_25(state, unused, out_7967737535555518081);
}
void car_H_25(double *state, double *unused, double *out_4817880321206984679) {
  H_25(state, unused, out_4817880321206984679);
}
void car_h_24(double *state, double *unused, double *out_38287872680748422) {
  h_24(state, unused, out_38287872680748422);
}
void car_H_24(double *state, double *unused, double *out_6990529920212484245) {
  H_24(state, unused, out_6990529920212484245);
}
void car_h_30(double *state, double *unused, double *out_3122496255442127692) {
  h_30(state, unused, out_3122496255442127692);
}
void car_H_30(double *state, double *unused, double *out_2299547362699736052) {
  H_30(state, unused, out_2299547362699736052);
}
void car_h_26(double *state, double *unused, double *out_2091745227628403381) {
  h_26(state, unused, out_2091745227628403381);
}
void car_H_26(double *state, double *unused, double *out_8559383640081040903) {
  H_26(state, unused, out_8559383640081040903);
}
void car_h_27(double *state, double *unused, double *out_6081054456872404523) {
  h_27(state, unused, out_6081054456872404523);
}
void car_H_27(double *state, double *unused, double *out_75953291515792835) {
  H_27(state, unused, out_75953291515792835);
}
void car_h_29(double *state, double *unused, double *out_4472750121952623734) {
  h_29(state, unused, out_4472750121952623734);
}
void car_H_29(double *state, double *unused, double *out_1789316018385343868) {
  H_29(state, unused, out_1789316018385343868);
}
void car_h_28(double *state, double *unused, double *out_8893606560142257906) {
  h_28(state, unused, out_8893606560142257906);
}
void car_H_28(double *state, double *unused, double *out_6871715035454874442) {
  H_28(state, unused, out_6871715035454874442);
}
void car_h_31(double *state, double *unused, double *out_86464884801468896) {
  h_31(state, unused, out_86464884801468896);
}
void car_H_31(double *state, double *unused, double *out_9185591742314392379) {
  H_31(state, unused, out_9185591742314392379);
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
