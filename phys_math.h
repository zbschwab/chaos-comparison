#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <linmath.h>

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

// signum function
int sgn(double val);

// takes ddp state and returns derivative of each state variable
deriv_ddp_t deriv_lddp(state_ddp_t *s, cons_ddp_t *c, double t, double gamma);

// takes ddp state and returns derivative of each state variable
deriv_ddp_t deriv_qddp(state_ddp_t *s, cons_ddp_t *c, double t, double gamma);

// takes dp state and returns derivative of each state variable
deriv_dp_t deriv_dp(state_dp_t *s, cons_dp_t *c);

// integration step w/ RK45 alg for ddp
double rk45_ddp_step(state_ddp_t *s, cons_ddp_t *c, double t, double dt, double gamma);

// integration step w/ RK45 alg for dp
double rk45_dp_step(state_dp_t *s, cons_dp_t *c, double dt);
