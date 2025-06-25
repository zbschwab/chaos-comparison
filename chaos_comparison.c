#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <cblas.h>

#include "phys_math.h"

#define COMPUTE_STEPS 1000000

int main(void) {
    // lddp
    cons_ddp_t* c_lddp = (cons_ddp_t*)malloc(sizeof(cons_ddp_t));
    state_ddp_t* s_lddp = (state_ddp_t*)malloc(sizeof(state_ddp_t));;
    double* lddp_jac = (double*)malloc(sizeof(double));

    c_lddp->beta = 1;
    c_lddp->omega0 = 1;
    c_lddp->omegad = 1;
    s_lddp->phi = 2;
    s_lddp->omega = 2;

    double t_step = 0.01;

    lddp_jac = jac_lddp(s_lddp, c_lddp, t_step);

    // test:
    printf(" %lf ", lddp_jac);

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