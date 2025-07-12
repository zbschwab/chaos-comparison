#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <cblas.h>


# pragma once

// linearly-damped driven pendulum = lddp
// quadratically-damped driven pendulum = qddp
// double pendulum = dp

// ddp constants (same for linear and quadratic damping)
typedef struct cons_ddp {
    double omega0;  // natural frequency
    double omegad;  // driving frequency
    double beta;  // damping
} cons_ddp_t;

// ddp state: pendulum angle and velocity
typedef struct state_ddp {
    double phi;
    double omega;
} state_ddp_t;

// ddp velocity and acceleration
typedef struct deriv_ddp {
    double dphi;
    double d2phi;
} deriv_ddp_t;

// dp constants
typedef struct cons_dp {
    double m_1;  // mass at end of 1st pendulum in kg
    double m_2;  // mass at end of 2nd pendulum in kg
    double l_1;  // length of 1st pendulum in meters
    double l_2;  // length of 2nd pendulum in meters
} cons_dp_t;

// dp state: pendulum angles and velocities
typedef struct state_dp {
    double theta_1;
    double theta_2;
    double omega_1;
    double omega_2;
} state_dp_t;

// dp velocity and acceleration (of both pendulums)
typedef struct deriv_dp {
    double dtheta_1;
    double dtheta_2;
    double d2theta_1;
    double d2theta_2;
} deriv_dp_t;

// RK45 integrator step return struct (next state, normed error)
typedef struct step {
    double* s_next;
    double err;
} step_t;

// signum function
int sgn(double val);

// takes lddp state and returns derivative of each state variable
deriv_ddp_t deriv_lddp(state_ddp_t *s, cons_ddp_t *c, double t, double gamma);

// takes qddp state and returns derivative of each state variable
deriv_ddp_t deriv_qddp(state_ddp_t *s, cons_ddp_t *c, double t, double gamma);

// takes dp state and returns derivative of each state variable
deriv_dp_t deriv_dp(state_dp_t *s, cons_dp_t *c);

// calculate one step of lddp's phase space trajectory with RK45 integration algorithm
void rk45_lddp_step(step_t* next, state_ddp_t *s, cons_ddp_t *c, double *t, double *dt, double gamma);

// calculate one step of qddp's phase space trajectory with RK45 integration algorithm
int rk45_qddp_step(state_ddp_t *s, cons_ddp_t *c, double *t, double *dt, double gamma);

// calculate one step of dp's phase space trajectory with RK45 integration algorithm
int rk45_dp_step(state_dp_t *s, cons_dp_t *c, double *t, double *dt);

// takes lddp state and returns jacobian matrix (local linearization)
double* jac_lddp(state_ddp_t *s, cons_ddp_t *c);

// takes qddp state and returns jacobian matrix (local linearization)
double* jac_qddp(state_ddp_t *s, cons_ddp_t *c, double t);

// takes dp state and returns jacobian matrix (local linearization)
double* jac_dp(state_dp_t *s, cons_dp_t *c);

// calculate one step of lddp's deviation vector with a LTM (linearized tangent map) using RK45
int rk45_LTM_lddp_step(state_ddp_t *s, cons_ddp_t *c, double *dev_vec, double *dt);

// calculate one step of qddp's deviation vector with a LTM (linearized tangent map) using RK45
// double rk45_LTM_qddp_step(state_ddp_t *s, cons_ddp_t *c, double t, double dt, double gamma);

// // calculate one step of dp's deviation vector with a LTM using RK45
// double rk45_LTM_dp_step(state_dp_t *s, cons_dp_t *c, double t, double dt);