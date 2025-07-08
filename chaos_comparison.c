#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <cblas.h>

#include "phys_math.h"
#include "butcher_tableau.h"

#define COMPUTE_STEPS 1000000
#define GAMMA_STEPS 100
#define PI 3.1415927

int main(void) {
    // set constants and initial conditions:
    double t_init = 0.0;
    double dt_init = 0.0001;
    double* t = &t_init;
    double* dt = &dt_init;
    int status;

    // lddp
    cons_ddp_t* c_lddp = (cons_ddp_t*)malloc(sizeof(cons_ddp_t));
    state_ddp_t* s_lddp = (state_ddp_t*)malloc(sizeof(state_ddp_t));;
    double* lddp_jac = (double*)malloc(sizeof(double));

    c_lddp->beta = 3.0/4.0 * PI;
    c_lddp->omega0 = 3.0 * PI;
    c_lddp->omegad = 2.0 * PI;

    s_lddp->phi = 0.0;
    s_lddp->omega = 0.0;

    // qddp
    cons_ddp_t* c_qddp = (cons_ddp_t*)malloc(sizeof(cons_ddp_t));
    state_ddp_t* s_qddp = (state_ddp_t*)malloc(sizeof(state_ddp_t));;
    double* qddp_jac = (double*)malloc(sizeof(double));

    c_qddp->beta = 1.0;
    c_qddp->omega0 = 1.0;
    c_qddp->omegad = 1.0;

    s_qddp->phi = 2.0;
    s_qddp->omega = 2.0;

    // dp
    cons_dp_t* c_dp = (cons_dp_t*)malloc(sizeof(cons_dp_t));
    state_dp_t* s_dp = (state_dp_t*)malloc(sizeof(state_dp_t));;
    double* dp_jac = (double*)malloc(sizeof(double));

    c_dp->m_1 = 1.0;
    c_dp->m_2 = 1.0;
    c_dp->l_1 = 1.0;
    c_dp->l_2 = 1.0;

    s_dp->theta_1 = 2.0;
    s_dp->theta_2 = 2.0;
    s_dp->omega_1 = 2.0;
    s_dp->omega_2 = 2.0;
    
    // make pertubation parameter arrays:
    // lddp + qddp
    double gamma_arr[GAMMA_STEPS] = {0};
    double gamma_start = 1.0;
    double gamma_stop = 1.1;
    double gamma_step = (gamma_stop - gamma_start)/(GAMMA_STEPS - 1);

    int i = 0;
    for (double j = 0; j < gamma_step; j += gamma_step) {
        gamma_arr[i] = gamma_start + j * gamma_step;
        i++;
    }

    // dp

    // test:

    lddp_jac = jac_lddp(s_lddp, c_lddp);

    printf(" %lf ", *lddp_jac);

    double *lddp_traj = (double*)malloc(2*COMPUTE_STEPS*sizeof(double));

    i = 0;
    int j = 0;
    while (j < 20) {
        j++;
        status = rk45_lddp_step(s_lddp, c_lddp, t, dt, 1.1);
        if (status == EXIT_SUCCESS) {
            lddp_traj[i] = s_lddp->phi;
            lddp_traj[i+1] = s_lddp->omega;
            i += 2;
            printf("SUCCESS\n");
        } else {continue;}
    }

    // int i;
    // while (i < 2*COMPUTE_STEPS) {
    //     status = rk45_lddp_step(s_lddp, c_lddp, t, dt, 1.1);
    //     if (status == EXIT_SUCCESS) {
    //         lddp_traj[i] = s_lddp->phi;
    //         lddp_traj[i+1] = s_lddp->omega;
    //         i += 2;
    //     }
    // }

    // // initialize array to store convergence steps
    // double lddp_lyp_con[6000];

    // for (int i = 0; i < COMPUTE_STEPS; i++) {
    //     lddp_lyp_con[i] = 
    //     t_step = jac_lddp(s_lddp, c_lddp, t_step);
    //     //printf(" %lf ", lddp_jac[i]);
    // }

    // qddp


    // dp

    // // write data to csv
    // FILE *lddp_file;

    // lddp_file = fopen("lddp_lyp.csv", "w");
    // fprintf(lddp_file,"gamma, max_lyp\n");
    // for (int i = 0; i < gamma_step; i++) {
    //     fprintf(lddp_file, "%lf, %lf\n", gamma_arr[i], lddp_mlyp[i]);
    // }
    // fclose(lddp_file);

    // FILE *qddp_file;
    // qddp_file = fopen("qddp_lyp.csv", "w");
    // fprintf(qddp_file,"gamma, max_lyp\n");
    // for (int i = 0; i < gamma_step; i++) {
    //     fprintf(qddp_file, "%lf, %lf\n", gamma_arr[i], qddp_mlyp[i]);
    // }
    // fclose(qddp_file);

    free(s_lddp);
    free(s_qddp);
    free(s_dp);
    free(c_lddp);
    free(c_qddp);
    free(c_dp);

    free(lddp_traj);
}