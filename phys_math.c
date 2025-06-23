#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <linmath.h>

#include "phys_math.h"

#define PI 3.1415927
#define G 9.81


// signum function (int --> double)
int sgn(double val) {
    return (0.0 < val) - (val < 0.0);
}

// takes lddp state and returns the derivative of each state variable
deriv_ddp_t deriv_lddp(state_ddp_t *s, cons_ddp_t *c, double t, double gamma) {
    deriv_ddp_t d = {0};
    d.dphi = s->omega;
    d.d2phi = (gamma * pow(c->omega0, 2) * cos(c->omegad * t)) - (2 * c->beta * s->omega) - (pow(c->omega0, 2) * sin(s->phi));

    return d;
}


// takes qddp state and returns the derivative of each state variable
deriv_ddp_t deriv_qddp(state_ddp_t *s, cons_ddp_t *c, double t, double gamma) {
    deriv_ddp_t d = {0};
    d.dphi = s->omega;
    d.d2phi = (gamma * pow(c->omega0, 2) * cos(c->omegad * t)) - (2 * c->beta * sgn(s->omega) * pow(s->omega, 2)) - (pow(c->omega0, 2) * sin(s->phi));

    return d;
}


// takes dp state and returns the derivative of each state variable
deriv_dp_t deriv_dp(state_dp_t *s, cons_dp_t *c) {
    // calculate f and alpha for coupled ODEs (see equation sheet)
    double f_1, f_2, alpha_1, alpha_2, denom, M, diff;

    M = (c->m_2 / (c->m_1 + c->m_2));
    diff = s->theta_1 - s->theta_2;

    f_1 = - (c->l_2 / c->l_1) * M * pow(s->omega_2, 2) * sin(diff) - (G / c->l_1 * sin(s->theta_1));
    f_2 = (c->l_1 / c->l_2) * pow(s->omega_1, 2) * sin(diff) - (G / c->l_2 * sin(s->theta_2));

    alpha_1 = (c->l_2 / c->l_1) * M * cos(diff);
    alpha_2 = (c->l_1 / c->l_2) * cos(diff);
    denom = 1 - (alpha_1 * alpha_2);

    deriv_dp_t d = {0};
    d.dtheta_1 = s->omega_1;
    d.dtheta_2 = s->omega_2;
    d.d2theta_1 = (f_1 - alpha_1 * f_2) / denom;
    d.d2theta_2 = (- alpha_2 * f_1 + f_2) / denom;

    return d;
}

