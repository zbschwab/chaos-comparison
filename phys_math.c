#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <cblas.h>

#include "phys_math.h"
#include "butcher_tableau.h"

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

// integration step w/ adaptive RK45 alg for lddp
double rk45_lddp_step(state_ddp_t *s, cons_ddp_t *c, double t, double dt, double gamma) {
    state_ddp_t temp = *s;
    deriv_ddp_t k1, k2, k3, k4, k5, k6;

    k1 = deriv_lddp(s, c, t + dt, gamma);

    temp.phi = s->phi + B12 * k1.dphi;
    temp.omega = s->omega + B12 * k1.d2phi;
    k2 = deriv_lddp(&temp, c, t + A2*dt, gamma);

    temp.phi = s->phi + B13 * k1.dphi + B23 * k2.dphi;
    temp.omega = s->omega + B13 * k1.d2phi + B23 * k2.d2phi;
    k3 = deriv_lddp(&temp, c, t + A3*dt, gamma);

    temp.phi = s->phi + B14 * k1.dphi + B24 * k2.dphi + B34 * k3.dphi;
    temp.omega = s->omega + B14 * k1.d2phi + B24 * k2.d2phi + B34 * k3.d2phi;
    k4 = deriv_lddp(&temp, c, t + A4*dt, gamma);

    temp.phi = s->phi + B15 * k1.dphi + B25 * k2.dphi + B35 * k3.dphi + B45 * k4.dphi;
    temp.omega = s->omega + B15 * k1.d2phi + B25 * k2.d2phi + B35 * k3.d2phi + B45 * k4.d2phi;
    k5 = deriv_lddp(&temp, c, t + A5*dt, gamma);

    temp.phi = s->phi + B16 * k1.dphi + B26 * k2.dphi + B36 * k3.dphi + B46 * k4.dphi + B56 * k5.dphi;
    temp.omega = s->omega + B16 * k1.d2phi + B26 * k2.d2phi + B36 * k3.d2phi + B46 * k4.d2phi + B56 * k5.d2phi;
    k6 = deriv_lddp(&temp, c, t + A6*dt, gamma);

    // calculate new state
    s->phi += dt * (CH1 * k1.dphi + CH3 * k3.dphi + CH4 * k4.dphi + CH5 * k5.dphi + CH6 * k6.dphi);
    s->omega += dt * (CH1 * k1.d2phi + CH3 * k3.d2phi + CH4 * k4.d2phi + CH5 * k5.d2phi + CH6 * k6.d2phi);

    // calculate truncation error
    double dphi_err, d2phi_err, trunc_err;
    dphi_err = fabs(CT1 * k1.dphi + CT3 * k3.dphi + CT4 * k4.dphi + CT5 * k5.dphi + CT6 * k6.dphi);
    d2phi_err = fabs(CT1 * k1.d2phi + CT3 * k3.d2phi + CT4 * k4.d2phi + CT5 * k5.d2phi + CT6 * k6.d2phi);
    trunc_err = fmax(dphi_err, d2phi_err);

    return trunc_err;
}

// integration step w/ adaptive RK45 alg for qddp
double rk45_qddp_step(state_ddp_t *s, cons_ddp_t *c, double t, double dt, double gamma) {
    state_ddp_t temp = *s;
    deriv_ddp_t k1, k2, k3, k4, k5, k6;

    k1 = deriv_qddp(s, c, t + dt, gamma);

    temp.phi = s->phi + B12 * k1.dphi;
    temp.omega = s->omega + B12 * k1.d2phi;
    k2 = deriv_qddp(&temp, c, t + A2*dt, gamma);

    temp.phi = s->phi + B13 * k1.dphi + B23 * k2.dphi;
    temp.omega = s->omega + B13 * k1.d2phi + B23 * k2.d2phi;
    k3 = deriv_qddp(&temp, c, t + A3*dt, gamma);

    temp.phi = s->phi + B14 * k1.dphi + B24 * k2.dphi + B34 * k3.dphi;
    temp.omega = s->omega + B14 * k1.d2phi + B24 * k2.d2phi + B34 * k3.d2phi;
    k4 = deriv_qddp(&temp, c, t + A4*dt, gamma);

    temp.phi = s->phi + B15 * k1.dphi + B25 * k2.dphi + B35 * k3.dphi + B45 * k4.dphi;
    temp.omega = s->omega + B15 * k1.d2phi + B25 * k2.d2phi + B35 * k3.d2phi + B45 * k4.d2phi;
    k5 = deriv_qddp(&temp, c, t + A5*dt, gamma);

    temp.phi = s->phi + B16 * k1.dphi + B26 * k2.dphi + B36 * k3.dphi + B46 * k4.dphi + B56 * k5.dphi;
    temp.omega = s->omega + B16 * k1.d2phi + B26 * k2.d2phi + B36 * k3.d2phi + B46 * k4.d2phi + B56 * k5.d2phi;
    k6 = deriv_qddp(&temp, c, t + A6*dt, gamma);

    // calculate new state
    s->phi += dt * (CH1 * k1.dphi + CH3 * k3.dphi + CH4 * k4.dphi + CH5 * k5.dphi + CH6 * k6.dphi);
    s->omega += dt * (CH1 * k1.d2phi + CH3 * k3.d2phi + CH4 * k4.d2phi + CH5 * k5.d2phi + CH6 * k6.d2phi);

    // calculate truncation error
    double dphi_err, d2phi_err, trunc_err;
    dphi_err = fabs(CT1 * k1.dphi + CT3 * k3.dphi + CT4 * k4.dphi + CT5 * k5.dphi + CT6 * k6.dphi);
    d2phi_err = fabs(CT1 * k1.d2phi + CT3 * k3.d2phi + CT4 * k4.d2phi + CT5 * k5.d2phi + CT6 * k6.d2phi);
    trunc_err = fmax(dphi_err, d2phi_err);

    return trunc_err;
}

// integration step w/ RK45 alg for dp
double rk45_dp_step(state_dp_t *s, cons_dp_t *c, double dt) {
    state_dp_t temp = *s;
    deriv_dp_t k1, k2, k3, k4, k5, k6;

    k1 = deriv_dp(s, c);

    temp.theta_1 = s->theta_1 + B12 * k1.dtheta_1;
    temp.theta_2 = s->theta_2 + B12 * k1.dtheta_2;
    temp.omega_1 = s->omega_1 + B12 * k1.d2theta_1;
    temp.omega_2 = s->omega_2 + B12 * k1.d2theta_2;
    k2 = deriv_dp(&temp, c);

    temp.theta_1 = s->theta_1 + B13 * k1.dtheta_1 + B23 * k2.dtheta_1;
    temp.theta_2 = s->theta_2 + B13 * k1.dtheta_2 + B23 * k2.dtheta_2;
    temp.omega_1= s->omega_1 + B13 * k1.d2theta_1 + B23 * k2.d2theta_1;
    temp.omega_2= s->omega_2 + B13 * k1.d2theta_2 + B23 * k2.d2theta_2;
    k3 = deriv_dp(&temp, c);

    temp.theta_1 = s->theta_1 + B14 * k1.dtheta_1 + B24 * k2.dtheta_1 + B34 * k3.dtheta_1;
    temp.theta_2 = s->theta_2 + B14 * k1.dtheta_2 + B24 * k2.dtheta_2 + B34 * k3.dtheta_2;
    temp.omega_1= s->omega_1 + B14 * k1.d2theta_1 + B24 * k2.d2theta_1 + B34 * k3.d2theta_1;
    temp.omega_2= s->omega_2 + B14 * k1.d2theta_2 + B24 * k2.d2theta_2 + B34 * k3.d2theta_2;
    k4 = deriv_dp(&temp, c);

    temp.theta_1 = s->theta_1 + B15 * k1.dtheta_1 + B25 * k2.dtheta_1 + B35 * k3.dtheta_1 + B45 * k4.dtheta_1;
    temp.theta_2 = s->theta_2 + B15 * k1.dtheta_2 + B25 * k2.dtheta_2 + B35 * k3.dtheta_2 + B45 * k4.dtheta_2;
    temp.omega_1 = s->omega_1 + B15 * k1.d2theta_1 + B25 * k2.d2theta_1 + B35 * k3.d2theta_1 + B45 * k4.d2theta_1;
    temp.omega_2 = s->omega_2 + B15 * k1.d2theta_2 + B25 * k2.d2theta_2 + B35 * k3.d2theta_2 + B45 * k4.d2theta_2;
    k5 = deriv_dp(&temp, c);

    temp.theta_1 = s->theta_1 + B16 * k1.dtheta_1 + B26 * k2.dtheta_1 + B36 * k3.dtheta_1 + B46 * k4.dtheta_1 + B56 * k5.dtheta_1;
    temp.theta_2 = s->theta_2 + B16 * k1.dtheta_2 + B26 * k2.dtheta_2 + B36 * k3.dtheta_2 + B46 * k4.dtheta_2 + B56 * k5.dtheta_2;
    temp.omega_1 = s->omega_1 + B16 * k1.d2theta_1 + B26 * k2.d2theta_1 + B36 * k3.d2theta_1 + B46 * k4.d2theta_1 + B56 * k5.d2theta_1;
    temp.omega_2 = s->omega_2 + B16 * k1.d2theta_2 + B26 * k2.d2theta_2 + B36 * k3.d2theta_2 + B46 * k4.d2theta_2 + B56 * k5.d2theta_2;
    k6 = deriv_dp(&temp, c);

    // calculate new state
    s->theta_1 += dt * (CH1 * k1.dtheta_1 + CH3 * k3.dtheta_1 + CH4 * k4.dtheta_1 + CH5 * k5.dtheta_1 + CH6 * k6.dtheta_1);
    s->theta_2 += dt * (CH1 * k1.dtheta_2 + CH3 * k3.dtheta_2 + CH4 * k4.dtheta_2 + CH5 * k5.dtheta_2 + CH6 * k6.dtheta_2);
    s->omega_1 += dt * (CH1 * k1.d2theta_1 + CH3 * k3.d2theta_1 + CH4 * k4.d2theta_1 + CH5 * k5.d2theta_1 + CH6 * k6.d2theta_1);
    s->omega_2 += dt * (CH1 * k1.d2theta_2 + CH3 * k3.d2theta_2 + CH4 * k4.d2theta_2 + CH5 * k5.d2theta_2 + CH6 * k6.d2theta_2);

    // calculate truncation error
    double dtheta_1_err, dtheta_2_err, d2theta_1_err, d2theta_2_err, trunc_err;
    dtheta_1_err = fabs(CT1 * k1.dtheta_1 + CT3 * k3.dtheta_1 + CT4 * k4.dtheta_1 + CT5 * k5.dtheta_1 + CT6 * k6.dtheta_1);
    dtheta_2_err = fabs(CT1 * k1.dtheta_2 + CT3 * k3.dtheta_2 + CT4 * k4.dtheta_2 + CT5 * k5.dtheta_2 + CT6 * k6.dtheta_2);
    d2theta_1_err = fabs(CT1 * k1.d2theta_1 + CT3 * k3.d2theta_1 + CT4 * k4.d2theta_1 + CT5 * k5.d2theta_1 + CT6 * k6.d2theta_1);
    d2theta_2_err = fabs(CT1 * k1.d2theta_2 + CT3 * k3.d2theta_2 + CT4 * k4.d2theta_2 + CT5 * k5.d2theta_2 + CT6 * k6.d2theta_2);
    trunc_err = fmax(fmax(dtheta_1_err, dtheta_2_err), fmax(d2theta_1_err, d2theta_2_err));

    return trunc_err;
}

// takes lddp state and returns jacobian matrix (local linearization)
double* jac_lddp(state_ddp_t *s, cons_ddp_t *c, double t) {
    double* jac = (double*)malloc(4*sizeof(double));
    jac[0] = 0;
    jac[1] = 1;
    jac[2] = - pow(c->omega0, 2) * cos(s->phi);
    jac[3] = - 2 * c->beta;

    return jac;
}
