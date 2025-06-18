#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

# pragma once


// double pendulum constants
typedef struct cons_dp {
    double m_1;  // mass at end of 1st pendulum in kg
    double m_2;  // mass at end of 2nd pendulum in kg
    double l_1;  // length of 1st pendulum in meters
    double l_2;  // length of 2nd pendulum in meters
} cons_dp_t;

// double pendulum state: pendulum angles and their 1st derivatives
typedef struct state_dp {
    double theta_1;
    double theta_2;
    double omega_1;
    double omega_2;
} state_dp_t;

// double pendulum 1st and 2nd derivatives (of both angles)
typedef struct deriv_dp {
    double dtheta_1;
    double dtheta_2;
    double d2theta_1;
    double d2theta_2;
} deriv_dp_t;

// takes double pendulum state and returns derivative of each state variable
deriv_dp_t deriv(state_dp_t *s, cons_dp_t *c);