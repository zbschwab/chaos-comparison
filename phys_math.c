#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <cblas.h>

#include "phys_math.h"
#include "butcher_tableau.h"

#define PI 3.1415927
#define G 9.81

#define ALPHA 1.0
#define BETA 1.0
#define ATOL 0.0001
#define RTOL 0.001
#define DT_MIN 0.00001


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
void rk45_lddp_step(traj_step_ddp_t* next, state_ddp_t *s, cons_ddp_t *c, double *t, double *dt, double gamma) {
    deriv_ddp_t k1, k2, k3, k4, k5, k6, k7;
    state_ddp_t temp = *s;
    state_ddp_t b4, b5;

    k1 = deriv_lddp(s, c, *t, gamma);

    temp.phi = s->phi + A12 * k1.dphi;
    temp.omega = s->omega + A12 * k1.d2phi;
    k2 = deriv_lddp(&temp, c, *t + C2*(*dt), gamma);

    temp.phi = s->phi + A13 * k1.dphi + A23 * k2.dphi;
    temp.omega = s->omega + A13 * k1.d2phi + A23 * k2.d2phi;
    k3 = deriv_lddp(&temp, c, *t + C3*(*dt), gamma);

    temp.phi = s->phi + A14 * k1.dphi + A24 * k2.dphi + A34 * k3.dphi;
    temp.omega = s->omega + A14 * k1.d2phi + A24 * k2.d2phi + A34 * k3.d2phi;
    k4 = deriv_lddp(&temp, c, *t + C4*(*dt), gamma);

    temp.phi = s->phi + A15 * k1.dphi + A25 * k2.dphi + A35 * k3.dphi + A45 * k4.dphi;
    temp.omega = s->omega + A15 * k1.d2phi + A25 * k2.d2phi + A35 * k3.d2phi + A45 * k4.d2phi;
    k5 = deriv_lddp(&temp, c, *t + C5*(*dt), gamma);

    temp.phi = s->phi + A16 * k1.dphi + A26 * k2.dphi + A36 * k3.dphi + A46 * k4.dphi + A56 * k5.dphi;
    temp.omega = s->omega + A16 * k1.d2phi + A26 * k2.d2phi + A36 * k3.d2phi + A46 * k4.d2phi + A56 * k5.d2phi;
    k6 = deriv_lddp(&temp, c, *t + C6*(*dt), gamma);

    temp.phi = s->phi + A17 * k1.dphi + A37 * k3.dphi + A47 * k4.dphi + A57 * k5.dphi + A67 * k6.dphi;
    temp.omega = s->omega + A17 * k1.d2phi + A37 * k3.d2phi + A47 * k4.d2phi + A57 * k5.d2phi + A67 * k6.d2phi;
    k7 = deriv_lddp(&temp, c, *t + C7*(*dt), gamma);

    // 4th-order estimate
    b4.phi = s->phi + (*dt) * (A17 * k1.dphi + A37 * k3.dphi + A47 * k4.dphi + A57 * k5.dphi + A67 * k6.dphi);
    b4.omega = s->omega + (*dt) * (A17 * k1.d2phi + A37 * k3.d2phi + A47 * k4.d2phi + A57 * k5.d2phi + A67 * k6.d2phi);

    // 5th-order estimate
    b5.phi = s->phi + (*dt) * (B1 * k1.dphi + B3 * k3.dphi + B4 * k4.dphi + B5 * k5.dphi + B6 * k6.dphi + B7 * k7.dphi);
    b5.omega = s->omega + (*dt) * (B1 * k1.d2phi + B3 * k3.d2phi + B4 * k4.d2phi + B5 * k5.d2phi + B6 * k6.d2phi + B7 * k7.d2phi);
    
    // compute + norm truncation error
    double scale_phi = ATOL + RTOL * fabs(s->phi);
    double scale_omega = ATOL + RTOL * fabs(s->omega);
    if (scale_phi < 1e-14) scale_phi = 1e-14;
    if (scale_omega < 1e-14) scale_omega = 1e-14;

    double e0 = (b5.phi - b4.phi) / scale_phi;
    double e1 = (b5.omega - b4.omega) / scale_omega;

    double err_norm = sqrt((e0*e0 + e1*e1) / 2.0);

    // update step
    next->s_next = b5;
    next->err = err_norm;

}

// integration step w/ adaptive RK45 alg for qddp
int rk45_qddp_step(state_ddp_t *s, cons_ddp_t *c, double *t, double *dt, double gamma) {
    deriv_ddp_t k1, k2, k3, k4, k5, k6, k7;
    state_ddp_t temp = *s;
    state_ddp_t b4, b5;

    k1 = deriv_qddp(s, c, *t, gamma);

    temp.phi = s->phi + A12 * k1.dphi;
    temp.omega = s->omega + A12 * k1.d2phi;
    k2 = deriv_qddp(&temp, c, *t + C2*(*dt), gamma);

    temp.phi = s->phi + A13 * k1.dphi + A23 * k2.dphi;
    temp.omega = s->omega + A13 * k1.d2phi + A23 * k2.d2phi;
    k3 = deriv_qddp(&temp, c, *t + C3*(*dt), gamma);

    temp.phi = s->phi + A14 * k1.dphi + A24 * k2.dphi + A34 * k3.dphi;
    temp.omega = s->omega + A14 * k1.d2phi + A24 * k2.d2phi + A34 * k3.d2phi;
    k4 = deriv_qddp(&temp, c, *t + C4*(*dt), gamma);

    temp.phi = s->phi + A15 * k1.dphi + A25 * k2.dphi + A35 * k3.dphi + A45 * k4.dphi;
    temp.omega = s->omega + A15 * k1.d2phi + A25 * k2.d2phi + A35 * k3.d2phi + A45 * k4.d2phi;
    k5 = deriv_qddp(&temp, c, *t + C5*(*dt), gamma);

    temp.phi = s->phi + A16 * k1.dphi + A26 * k2.dphi + A36 * k3.dphi + A46 * k4.dphi + A56 * k5.dphi;
    temp.omega = s->omega + A16 * k1.d2phi + A26 * k2.d2phi + A36 * k3.d2phi + A46 * k4.d2phi + A56 * k5.d2phi;
    k6 = deriv_qddp(&temp, c, *t + C6*(*dt), gamma);

    temp.phi = s->phi + A17 * k1.dphi + A37 * k3.dphi + A47 * k4.dphi + A57 * k5.dphi + A67 * k6.dphi;
    temp.omega = s->omega + A17 * k1.d2phi + A37 * k3.d2phi + A47 * k4.d2phi + A57 * k5.d2phi + A67 * k6.d2phi;
    k7 = deriv_qddp(&temp, c, *t + C7*(*dt), gamma);

    // 4th-order estimate
    b4.phi = s->phi + (*dt) * (A17 * k1.dphi + A37 * k3.dphi + A47 * k4.dphi + A57 * k5.dphi + A67 * k6.dphi);
    b4.omega = s->omega + (*dt) * (A17 * k1.d2phi + A37 * k3.d2phi + A47 * k4.d2phi + A57 * k5.d2phi + A67 * k6.d2phi);

    // 5th-order estimate
    b5.phi = s->phi + (*dt) * (B1 * k1.dphi + B3 * k3.dphi + B4 * k4.dphi + B5 * k5.dphi + B6 * k6.dphi + B7 * k7.dphi);
    b5.omega = s->omega + (*dt) * (B1 * k1.d2phi + B3 * k3.d2phi + B4 * k4.d2phi + B5 * k5.d2phi + B6 * k6.d2phi + B7 * k7.d2phi);
    
    // compute + norm truncation error
    double scale_phi = ATOL + RTOL * fabs(s->phi);
    double scale_omega = ATOL + RTOL * fabs(s->omega);
    if (scale_phi < 1e-14) scale_phi = 1e-14;
    if (scale_omega < 1e-14) scale_omega = 1e-14;

    double e0 = (b5.phi - b4.phi) / scale_phi;
    double e1 = (b5.omega - b4.omega) / scale_omega;

    double err_norm = sqrt((e0*e0 + e1*e1) / 2.0);

    // compute new dt (w/ safety clamp)
    double exp = 1.0 / 5.0;
    double factor = 0.9 * pow(1.0 / err_norm, exp);
    if (factor < 0.1) {factor = 0.1;}
    if (factor > 5.0) {factor = 5.0;}
    (*dt) *= factor;

    double step = *dt;
    printf("step complete. dt = %.12e, err_norm = %lf, rk4 = %lf, %lf, rk5 = %lf, %lf\n", step, err_norm, b4.phi, b4.omega, b5.phi, b5.omega);

    // update state if error norm < 1
    if (err_norm <= 1.0) {
        s->phi = b5.phi;
        s->omega = b5.omega;
        *t += (*dt);  // update t
        return EXIT_SUCCESS;
    }
    
    return EXIT_FAILURE;
}

// integration step w/ RK45 alg for dp
int rk45_dp_step(state_dp_t *s, cons_dp_t *c, double *t, double *dt) {
    state_dp_t temp = *s;
    deriv_dp_t k1, k2, k3, k4, k5, k6, k7;
    state_dp_t b4, b5;

    k1 = deriv_dp(s, c);

    temp.theta_1 = s->theta_1 + A12 * k1.dtheta_1;
    temp.theta_2 = s->theta_2 + A12 * k1.dtheta_2;
    temp.omega_1 = s->omega_1 + A12 * k1.d2theta_1;
    temp.omega_2 = s->omega_2 + A12 * k1.d2theta_2;
    k2 = deriv_dp(&temp, c);

    temp.theta_1 = s->theta_1 + A13 * k1.dtheta_1 + A23 * k2.dtheta_1;
    temp.theta_2 = s->theta_2 + A13 * k1.dtheta_2 + A23 * k2.dtheta_2;
    temp.omega_1= s->omega_1 + A13 * k1.d2theta_1 + A23 * k2.d2theta_1;
    temp.omega_2= s->omega_2 + A13 * k1.d2theta_2 + A23 * k2.d2theta_2;
    k3 = deriv_dp(&temp, c);

    temp.theta_1 = s->theta_1 + A14 * k1.dtheta_1 + A24 * k2.dtheta_1 + A34 * k3.dtheta_1;
    temp.theta_2 = s->theta_2 + A14 * k1.dtheta_2 + A24 * k2.dtheta_2 + A34 * k3.dtheta_2;
    temp.omega_1= s->omega_1 + A14 * k1.d2theta_1 + A24 * k2.d2theta_1 + A34 * k3.d2theta_1;
    temp.omega_2= s->omega_2 + A14 * k1.d2theta_2 + A24 * k2.d2theta_2 + A34 * k3.d2theta_2;
    k4 = deriv_dp(&temp, c);

    temp.theta_1 = s->theta_1 + A15 * k1.dtheta_1 + A25 * k2.dtheta_1 + A35 * k3.dtheta_1 + A45 * k4.dtheta_1;
    temp.theta_2 = s->theta_2 + A15 * k1.dtheta_2 + A25 * k2.dtheta_2 + A35 * k3.dtheta_2 + A45 * k4.dtheta_2;
    temp.omega_1 = s->omega_1 + A15 * k1.d2theta_1 + A25 * k2.d2theta_1 + A35 * k3.d2theta_1 + A45 * k4.d2theta_1;
    temp.omega_2 = s->omega_2 + A15 * k1.d2theta_2 + A25 * k2.d2theta_2 + A35 * k3.d2theta_2 + A45 * k4.d2theta_2;
    k5 = deriv_dp(&temp, c);

    temp.theta_1 = s->theta_1 + A16 * k1.dtheta_1 + A26 * k2.dtheta_1 + A36 * k3.dtheta_1 + A46 * k4.dtheta_1 + A56 * k5.dtheta_1;
    temp.theta_2 = s->theta_2 + A16 * k1.dtheta_2 + A26 * k2.dtheta_2 + A36 * k3.dtheta_2 + A46 * k4.dtheta_2 + A56 * k5.dtheta_2;
    temp.omega_1 = s->omega_1 + A16 * k1.d2theta_1 + A26 * k2.d2theta_1 + A36 * k3.d2theta_1 + A46 * k4.d2theta_1 + A56 * k5.d2theta_1;
    temp.omega_2 = s->omega_2 + A16 * k1.d2theta_2 + A26 * k2.d2theta_2 + A36 * k3.d2theta_2 + A46 * k4.d2theta_2 + A56 * k5.d2theta_2;
    k6 = deriv_dp(&temp, c);

    temp.theta_1 = s->theta_1 + A17 * k1.dtheta_1 + A37 * k3.dtheta_1 + A47 * k4.dtheta_1 + A57 * k5.dtheta_1 + A67 * k2.dtheta_1;
    temp.theta_2 = s->theta_2 + A17 * k1.dtheta_2 + A37 * k3.dtheta_2 + A47 * k4.dtheta_2 + A57 * k5.dtheta_2 + A67 * k2.dtheta_2;
    temp.omega_1 = s->omega_1 + A17 * k1.d2theta_1 + A37 * k3.d2theta_1 + A47 * k4.d2theta_1 + A57 * k5.d2theta_1 + A67 * k2.d2theta_1;
    temp.omega_2 = s->omega_2 + A17 * k1.d2theta_2 + A37 * k3.d2theta_2 + A47 * k4.d2theta_2 + A57 * k5.d2theta_2 + A67 * k2.d2theta_1;
    k7 = deriv_dp(&temp, c);

    // 4th-order estimate
    b4.theta_1 = s->theta_1 + (*dt) * (A17 * k1.dtheta_1 + A37 * k3.dtheta_1 + A47 * k4.dtheta_1 + A57 * k5.dtheta_1 + A67 * k6.dtheta_1);
    b4.theta_2 = s->theta_2 + (*dt) * (A17 * k1.dtheta_2 + A37 * k3.dtheta_2 + A47 * k4.dtheta_2 + A57 * k5.dtheta_2 + A67 * k6.dtheta_2);
    b4.omega_1 = s->omega_2 + (*dt) * (A17 * k1.d2theta_1 + A37 * k3.d2theta_1 + A47 * k4.d2theta_1 + A57 * k5.d2theta_1 + A67 * k6.d2theta_1);
    b4.omega_2 = s->omega_2 + (*dt) * (A17 * k1.d2theta_2 + A37 * k3.d2theta_2 + A47 * k4.d2theta_2 + A57 * k5.d2theta_2 + A67 * k6.d2theta_2);

    // 5th-order estimate
    b5.theta_1 = s->theta_1 + (*dt) * (B1 * k1.dtheta_1 + B3 * k3.dtheta_1 + B4 * k4.dtheta_1 + B5 * k5.dtheta_1 + B6 * k6.dtheta_1 + B7 * k7.dtheta_1);
    b5.theta_2 = s->theta_2 + (*dt) * (B1 * k1.dtheta_2 + B3 * k3.dtheta_2 + B4 * k4.dtheta_2 + B5 * k5.dtheta_2 + B6 * k6.dtheta_2 + B7 * k7.dtheta_2);
    b5.omega_1 = s->omega_1 + (*dt) * (B1 * k1.d2theta_1 + B3 * k3.d2theta_1 + B4 * k4.d2theta_1 + B5 * k5.d2theta_1 + B6 * k6.d2theta_1 + B7 * k7.d2theta_1);
    b5.omega_2 = s->omega_2 + (*dt) * (B1 * k1.d2theta_2 + B3 * k3.d2theta_2 + B4 * k4.d2theta_2 + B5 * k5.d2theta_2 + B6 * k6.d2theta_2 + B7 * k7.d2theta_1);

    // compute + norm truncation error
    double scale_theta_1 = ATOL + RTOL * fabs(s->theta_1);
    double scale_theta_2 = ATOL + RTOL * fabs(s->theta_2);
    double scale_omega_1 = ATOL + RTOL * fabs(s->omega_1);
    double scale_omega_2 = ATOL + RTOL * fabs(s->omega_2);
    if (scale_theta_1 < 1e-14) scale_theta_1 = 1e-14;
    if (scale_theta_2 < 1e-14) scale_theta_2 = 1e-14;
    if (scale_omega_1 < 1e-14) scale_omega_1 = 1e-14;
    if (scale_omega_2 < 1e-14) scale_omega_2 = 1e-14;

    double e0_1 = (b5.theta_1 - b4.theta_1) / scale_theta_1;
    double e0_2 = (b5.theta_2 - b4.theta_2) / scale_theta_2;
    double e1_1 = (b5.omega_1 - b4.omega_1) / scale_omega_1;
    double e1_2 = (b5.omega_2 - b4.omega_2) / scale_omega_2;

    double err_norm = sqrt((e0_1*e0_1 + e0_2*e0_2 + e1_1*e1_1 + e1_2*e1_2) / 4.0);

    // compute new dt (w/ safety clamp)
    double exp = 1.0 / 5.0;
    double factor = 0.9 * pow(1.0 / err_norm, exp);
    if (factor < 1.0) {factor = 1.0;}
    if (factor > 5.0) {factor = 5.0;}
    (*dt) *= factor;

    double step = *dt;
    printf("dt = %.12e, err_norm = %lf, rk4 = %lf, %lf, %lf, %lf, rk5 = %lf, %lf, %lf, %lf\n", step, err_norm, b4.theta_1, b4.theta_2, b4.omega_1, b4.omega_2, b5.theta_1, b5.theta_2, b5.omega_1, b5.omega_2);
    
    // update state if error norm < 1
    if (err_norm <= 1.0) {
        s->theta_1 = b5.theta_1;
        s->theta_2 = b5.theta_2;
        s->omega_1 = b5.omega_1;
        s->omega_2 = b5.omega_2;
        *t += (*dt);  // update t
        return EXIT_SUCCESS;
    } 

    return EXIT_FAILURE;

}

// takes lddp state and returns jacobian matrix (local linearization)
double* jac_lddp(state_ddp_t *s, cons_ddp_t *c) {
    double* jac = (double*)calloc(4, sizeof(double));
    jac[0] = 0;
    jac[1] = 1;
    jac[2] = - pow(c->omega0, 2) * cos(s->phi);
    jac[3] = - 2 * c->beta;

    return jac;
}

// takes qddp state and returns jacobian matrix (local linearization)
double* jac_qddp(state_ddp_t *s, cons_ddp_t *c, double t) {
    double* jac = (double*)calloc(4, sizeof(double));
    jac[0] = 0;
    jac[1] = 1;
    jac[2] = - pow(c->omega0, 2) * cos(s->phi);
    jac[3] = - 4 * c->beta * sin(s->omega) * s->omega;

    return jac;
}

// takes dp state and returns jacobian matrix (local linearization)
double* jac_dp(state_dp_t *s, cons_dp_t *c) {
    double* jac = (double*)calloc(16, sizeof(double));
    // row 1 (left to right)
    jac[0] = 0;
    jac[1] = 0;
    jac[2] = 1;
    jac[3] = 0;

    // row 2 (left to right)
    jac[4] = 0;
    jac[5] = 0;
    jac[6] = 0;
    jac[7] = 1;

    // row 3 (left to right)
    jac[8] = (2 * c->m_2 * (pow(s->omega_2, 2) * c->l_2 * c->m_2 * sin(s->theta_1) * sin(s->theta_2) * sin(s->theta_1 - s->theta_2) + G * (c->m_1 + c->m_2) * sin(s->theta_2) - c->m_2 * (pow(s->omega_1, 2) * c->l_1 * sin(s->theta_2) * sin(s->theta_1 - s->theta_2) + G) * sin(s->theta_1) * cos(s->theta_1 - s->theta_2)) * sin(s->theta_1) * sin(s->theta_1 - s->theta_2) * cos(s->theta_1 - s->theta_2) + (-c->m_1 + c->m_2 * pow(cos(s->theta_1 - s->theta_2), 2) - c->m_2) * (-pow(s->omega_1, 2) * c->l_1 * c->m_2 * pow(sin(s->theta_1), 2) * sin(s->theta_2) * pow(cos(s->theta_1 - s->theta_2), 2) + pow(s->omega_2, 2) * c->l_2 * c->m_2 * pow(sin(s->theta_1), 2) * sin(s->theta_2) * cos(s->theta_1 - s->theta_2) - G * (c->m_1 + c->m_2) * sin(s->theta_2) * cos(s->theta_1) + c->m_2 * (pow(s->omega_1, 2) * c->l_1 * sin(s->theta_2) * sin(s->theta_1 - s->theta_2) + G) * pow(sin(s->theta_1), 2) * sin(s->theta_1 - s->theta_2))) / (c->l_1 * pow(-c->m_1 + c->m_2 * pow(cos(s->theta_1 - s->theta_2), 2) - c->m_2, 2) * pow(sin(s->theta_1), 2) * sin(s->theta_2));
    jac[9] = c->m_2 * ((-c->m_1 + c->m_2 * pow(cos(s->theta_1 - s->theta_2), 2) - c->m_2) * (-pow(s->omega_2, 2) * c->l_2 * pow(sin(s->theta_2), 2) * cos(s->theta_1 - s->theta_2) - (pow(s->omega_1, 2) * c->l_1 * sin(s->theta_2) * sin(s->theta_1 - s->theta_2) + G) * sin(s->theta_2) * sin(s->theta_1 - s->theta_2) + (pow(s->omega_1, 2) * c->l_1 * pow(sin(s->theta_2), 2) * cos(s->theta_1 - s->theta_2) + G * cos(s->theta_2)) * cos(s->theta_1 - s->theta_2)) * sin(s->theta_1) - 2 * (pow(s->omega_2, 2) * c->l_2 * c->m_2 * sin(s->theta_1) * sin(s->theta_2) * sin(s->theta_1 - s->theta_2) + G * (c->m_1 + c->m_2) * sin(s->theta_2) - c->m_2 * (pow(s->omega_1, 2) * c->l_1 * sin(s->theta_2) * sin(s->theta_1 - s->theta_2) + G) * sin(s->theta_1) * cos(s->theta_1 - s->theta_2)) * sin(s->theta_2) * sin(s->theta_1 - s->theta_2) * cos(s->theta_1 - s->theta_2)) / (c->l_1 * pow(-c->m_1 + c->m_2 * pow(cos(s->theta_1 - s->theta_2), 2) - c->m_2, 2) * sin(s->theta_1) * pow(sin(s->theta_2), 2));
    jac[10] = 2 * s->omega_1 * c->m_2 * sin(2 * s->theta_1 - 2 * s->theta_2) / (2 * c->m_1 - c->m_2 * cos(2 * s->theta_1 - 2 * s->theta_2) + c->m_2);
    jac[11] = 2 * s->omega_2 * c->l_2 * c->m_2 * sin(s->theta_1 - s->theta_2) / (c->l_1 * (-c->m_1 + c->m_2 * pow(cos(s->theta_1 - s->theta_2), 2) - c->m_2));

    // row 4 (left to right)
    jac[12] = (2 * c->m_2 * (pow(s->omega_1, 2) * c->l_1 * (c->m_1 + c->m_2) * sin(s->theta_1) * sin(s->theta_2) * sin(s->theta_1 - s->theta_2) + G * (c->m_1 + c->m_2) * sin(s->theta_1) - (pow(s->omega_2, 2) * c->l_2 * c->m_2 * sin(s->theta_1) * sin(s->theta_1 - s->theta_2) + G * (c->m_1 + c->m_2)) * sin(s->theta_2) * cos(s->theta_1 - s->theta_2)) * sin(s->theta_1) * sin(s->theta_1 - s->theta_2) * cos(s->theta_1 - s->theta_2) + (-c->m_1 + c->m_2 * pow(cos(s->theta_1 - s->theta_2), 2) - c->m_2) * (pow(s->omega_1, 2) * c->l_1 * (c->m_1 + c->m_2) * pow(sin(s->theta_1), 2) * cos(s->theta_1 - s->theta_2) + (pow(s->omega_2, 2) * c->l_2 * c->m_2 * sin(s->theta_1) * sin(s->theta_1 - s->theta_2) + G * (c->m_1 + c->m_2)) * sin(s->theta_1) * sin(s->theta_1 - s->theta_2) - (pow(s->omega_2, 2) * c->l_2 * c->m_2 * pow(sin(s->theta_1), 2) * cos(s->theta_1 - s->theta_2) - G * (c->m_1 + c->m_2) * cos(s->theta_1)) * cos(s->theta_1 - s->theta_2)) * sin(s->theta_2)) / (c->l_2 * pow(-c->m_1 + c->m_2 * pow(cos(s->theta_1 - s->theta_2), 2) - c->m_2, 2) * pow(sin(s->theta_1), 2) * sin(s->theta_2));
    jac[13] = (-2 * c->m_2 * (pow(s->omega_1, 2) * c->l_1 * (c->m_1 + c->m_2) * sin(s->theta_1) * sin(s->theta_2) * sin(s->theta_1 - s->theta_2) + G * (c->m_1 + c->m_2) * sin(s->theta_1) - (pow(s->omega_2, 2) * c->l_2 * c->m_2 * sin(s->theta_1) * sin(s->theta_1 - s->theta_2) + G * (c->m_1 + c->m_2)) * sin(s->theta_2) * cos(s->theta_1 - s->theta_2)) * sin(s->theta_2) * sin(s->theta_1 - s->theta_2) * cos(s->theta_1 - s->theta_2) + (-c->m_1 + c->m_2 * pow(cos(s->theta_1 - s->theta_2), 2) - c->m_2) * (-pow(s->omega_1, 2) * c->l_1 * (c->m_1 + c->m_2) * sin(s->theta_1) * pow(sin(s->theta_2), 2) * cos(s->theta_1 - s->theta_2) + pow(s->omega_2, 2) * c->l_2 * c->m_2 * sin(s->theta_1) * pow(sin(s->theta_2), 2) * pow(cos(s->theta_1 - s->theta_2), 2) - G * (c->m_1 + c->m_2) * sin(s->theta_1) * cos(s->theta_2) - (pow(s->omega_2, 2) * c->l_2 * c->m_2 * sin(s->theta_1) * sin(s->theta_1 - s->theta_2) + G * (c->m_1 + c->m_2)) * pow(sin(s->theta_2), 2) * sin(s->theta_1 - s->theta_2))) / (c->l_2 * pow(-c->m_1 + c->m_2 * pow(cos(s->theta_1 - s->theta_2), 2) - c->m_2, 2) * sin(s->theta_1) * pow(sin(s->theta_2), 2));
    jac[14] = 2 * s->omega_1 * c->l_1 * (c->m_1 + c->m_2) * sin(s->theta_1 - s->theta_2) / (c->l_2 * (-c->m_1 + c->m_2 * pow(cos(s->theta_1 - s->theta_2), 2) - c->m_2));
    jac[15] = 2 * s->omega_2 * c->m_2 * sin(2 * s->theta_1 - 2 * s->theta_2) / (2 * c->m_1 - c->m_2 * cos(2 * s->theta_1 - 2 * s->theta_2) + c->m_2);

    return jac;
}

// calculate one step of lddp's deviation vector with a LTM (linearized tangent map) using RK45
void rk45_LTM_lddp_step(dev_step_ddp_t* next, state_ddp_t *s, cons_ddp_t *c, double *dev_vec, double *dt) {
    double* jac;
    double k1[2], k2[2], k3[2], k4[2], k5[2], k6[2], k7[2], b4[2], b5[2], temp[2];

    // pointer to 4-len array representing 2x2 matrix
    jac = jac_lddp(s, c);
    printf("jac = [%lf, %lf; %lf, %lf]\n", jac[0], jac[1], jac[2], jac[3]);


    // jac x dev_vec = k1 (matrix-matrix multiplication)
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                2, 2,
                ALPHA,  jac, 2,
                        dev_vec, 1,
                BETA,   k1, 1);

    // dev_vec + k1 * A12 * dt = temp (scalar-matrix multiplication)
    // jac x temp = k2
    memcpy(temp, dev_vec, 2*sizeof(double));
    cblas_daxpy(2, A12 * (*dt), k1, 1, temp, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                2, 2,
                ALPHA,  jac, 2,
                        temp, 1,
                BETA,   k2, 1);
    
    // dev_vec + dt * (k1 * A13 + k2 * A23) = temp
    // jac x temp = k3
    memcpy(temp, dev_vec, 2*sizeof(double));
    cblas_daxpy(2, A13 * (*dt), k1, 1, temp, 1);
    cblas_daxpy(2, A23 * (*dt), k2, 1, temp, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                2, 2,
                ALPHA,  jac, 2,
                        temp, 1,
                BETA,   k3, 1);

    // dev_vec + dt * (k1 * A14 + k2 * A24 + k3 * A34) = temp
    // jac x temp = k4
    memcpy(temp, dev_vec, 2*sizeof(double));
    cblas_daxpy(2, A14 * (*dt), k1, 1, temp, 1);
    cblas_daxpy(2, A24 * (*dt), k2, 1, temp, 1);
    cblas_daxpy(2, A34 * (*dt), k3, 1, temp, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                2, 2,
                ALPHA,  jac, 2,
                        temp, 1,
                BETA,   k4, 1);
    
    // dev_vec + dt * (k1 * A15 + k2 * A25 + k3 * A35 + k4 * A45) = temp
    // jac x temp = k5
    memcpy(temp, dev_vec, 2*sizeof(double));
    cblas_daxpy(2, A15 * (*dt), k1, 1, temp, 1);
    cblas_daxpy(2, A25 * (*dt), k2, 1, temp, 1);
    cblas_daxpy(2, A35 * (*dt), k3, 1, temp, 1);
    cblas_daxpy(2, A45 * (*dt), k4, 1, temp, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                2, 2,
                ALPHA,  jac, 2,
                        temp, 1,
                BETA,   k5, 1);

    // dev_vec + dt * (k1 * A16 + k2 * A26 + k3 * A36 + k4 * A46 + k5 * A56) = temp
    // jac x temp = k6
    memcpy(temp, dev_vec, 2*sizeof(double));
    cblas_daxpy(2, A16 * (*dt), k1, 1, temp, 1);
    cblas_daxpy(2, A26 * (*dt), k2, 1, temp, 1);
    cblas_daxpy(2, A36 * (*dt), k3, 1, temp, 1);
    cblas_daxpy(2, A46 * (*dt), k4, 1, temp, 1);
    cblas_daxpy(2, A56 * (*dt), k5, 1, temp, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                2, 2,
                ALPHA,  jac, 2,
                        temp, 1,
                BETA,   k6, 1);

    // dev_vec + dt * (k1 * A17 + k3 * A37 + k4 * A47 + k5 * A57 + k6 * A67) = temp
    // jac x temp = k7
    memcpy(temp, dev_vec, 2*sizeof(double));
    cblas_daxpy(2, A17 * (*dt), k1, 1, temp, 1);
    cblas_daxpy(2, A37 * (*dt), k3, 1, temp, 1);
    cblas_daxpy(2, A47 * (*dt), k4, 1, temp, 1);
    cblas_daxpy(2, A57 * (*dt), k5, 1, temp, 1);
    cblas_daxpy(2, A67 * (*dt), k6, 1, temp, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                2, 2,
                ALPHA,  jac, 2,
                        temp, 1,
                BETA,   k7, 1);

    // 4th and 5th-order estimate
    for (int i = 0; i < 2; i++) {
        b4[i] = dev_vec[i] + (*dt) * (A17 * k1[i] + A37 * k3[i] + A47 * k4[i] + A57 * k5[i] + A67 *k6[i]);
        b5[i] = dev_vec[i] + (*dt) * (B1 * k1[i] + B3 * k3[i] + B4 * k4[i] + B5 * k5[i] + B6 * k6[i] + B7 * k7[i]);
    }

    // compute + norm truncation error
    double scale[2], e[2];
    double err_norm = 0;

    for (int i = 0; i < 2; i++) {
        scale[i] = ATOL + RTOL * fabs(jac[i]);
        if (scale[i] < 1e-14) {
            scale[i] = 1e-14;
        }
        e[i] = (b5[i] - b4[i]) / scale[i];
        err_norm += pow(e[i], 2);
    }

    err_norm = sqrt(err_norm / 2.0);

    // renorm step
    double norm = sqrt(pow(b5[0], 2) + pow(b5[1], 2));
    for (int i = 0; i < 2; i++) {
        b5[i] /= norm;
    }

    
    
    // update step
    for (int i = 0; i < 2; i++) {
        next->d_next[i] = b5[i];
    }
    next->err = err_norm;
    next->norm = norm;

}

// // calculate one step of qddp's deviation vector with a LTM (linearized tangent map) using RK45
// double rk45_LTM_qddp_step(state_ddp_t *s, cons_ddp_t *c, double t, double dt, double gamma) {

// }

// // calculate one step of dp's deviation vector with a LTM using RK45
// double rk45_LTM_dp_step(state_dp_t *s, cons_dp_t *c, double t, double dt) {

// }





