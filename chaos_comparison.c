#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <cblas.h>

#include "phys_math.h"

#define COMPUTE_STEPS 1000000
#define GAMMA_STEPS 100

int main(void) {
    // set constants and initial conditions:
    double t_step = 0.01;

    // lddp
    cons_ddp_t* c_lddp = (cons_ddp_t*)malloc(sizeof(cons_ddp_t));
    state_ddp_t* s_lddp = (state_ddp_t*)malloc(sizeof(state_ddp_t));;
    double* lddp_jac = (double*)malloc(sizeof(double));

    c_lddp->beta = 1;
    c_lddp->omega0 = 1;
    c_lddp->omegad = 1;

    s_lddp->phi = 2;
    s_lddp->omega = 2;

    // qddp
    cons_ddp_t* c_qddp = (cons_ddp_t*)malloc(sizeof(cons_ddp_t));
    state_ddp_t* s_qddp = (state_ddp_t*)malloc(sizeof(state_ddp_t));;
    double* qddp_jac = (double*)malloc(sizeof(double));

    c_qddp->beta = 1;
    c_qddp->omega0 = 1;
    c_qddp->omegad = 1;

    s_qddp->phi = 2;
    s_qddp->omega = 2;

    // dp
    cons_dp_t* c_dp = (cons_dp_t*)malloc(sizeof(cons_dp_t));
    state_dp_t* s_dp = (state_dp_t*)malloc(sizeof(state_dp_t));;
    double* dp_jac = (double*)malloc(sizeof(double));

    c_dp->m_1 = 1;
    c_dp->m_2 = 1;
    c_dp->l_1 = 1;
    c_dp->l_2 = 1;

    s_dp->theta_1 = 2;
    s_dp->theta_2 = 2;
    s_dp->omega_1 = 2;
    s_dp->omega_2 = 2;
    
    // make pertubation parameter arrays:
    // lddp + qddp
    double gamma_arr[GAMMA_STEPS] = {0};
    double gamma_start = 1;
    double gamma_stop = 1.1;
    double gamma_step = (gamma_stop - gamma_start)/(GAMMA_STEPS - 1);

    int i = 0;
    for (double j = 0; j < gamma_step; j += gamma_step) {
        gamma_arr[i] = gamma_start + j * gamma_step;
    }

    // dp

    // test:

    lddp_jac = jac_lddp(s_lddp, c_lddp, t_step);


    printf(" %lf ", *lddp_jac);

    // // initialize array to store convergence steps
    // double lddp_lyp_con[6000];

    // for (int i = 0; i < COMPUTE_STEPS; i++) {
    //     lddp_lyp_con[i] = 
    //     t_step = jac_lddp(s_lddp, c_lddp, t_step);
    //     //printf(" %lf ", lddp_jac[i]);
    // }

    // qddp


    // dp


}