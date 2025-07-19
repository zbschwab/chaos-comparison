#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <cblas.h>

#include "phys_math.h"
#include "butcher_tableau.h"

#define PI 3.1415927
#define EXP 0.2

#define COMPUTE_STEPS 1e6 
#define GAMMA_STEPS 100
#define THETA_STEPS 360  // range is 180, increment by 0.5 deg


int main(void) {
    // set initial system conditions: lddp
    cons_ddp_t* c_lddp = (cons_ddp_t*)malloc(sizeof(cons_ddp_t));
    state_ddp_t* s_lddp = (state_ddp_t*)malloc(sizeof(state_ddp_t));;
    double* lddp_jac = (double*)malloc(sizeof(double));

    c_lddp->beta = 1.0; //3.0/4.0 * PI;
    c_lddp->omega0 = 1.0; //3.0 * PI;
    c_lddp->omegad = 2.0 * PI;

    s_lddp->phi = 0.0;
    s_lddp->omega = 0.0;

    // set initial system conditions: qddp
    cons_ddp_t* c_qddp = (cons_ddp_t*)malloc(sizeof(cons_ddp_t));
    state_ddp_t* s_qddp = (state_ddp_t*)malloc(sizeof(state_ddp_t));;
    double* qddp_jac = (double*)malloc(sizeof(double));

    c_qddp->beta = 1.0;
    c_qddp->omega0 = 1.0;
    c_qddp->omegad = 1.0;

    s_qddp->phi = 2.0;
    s_qddp->omega = 2.0;

    // set initial system conditions: dp
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
    
    // make lddp + qddp perturbation param arrays (gamma = forcing constant)
    double gamma_arr[GAMMA_STEPS];
    double gamma_start = 1.0;
    double gamma_stop = 1.1;
    double gamma_step = (gamma_stop - gamma_start)/(GAMMA_STEPS - 1);

    printf("gamma_arr:\n");
    for (int i = 0; i < GAMMA_STEPS; i++) {
        gamma_arr[i] = gamma_start + i * gamma_step;
        //printf(" %lf ", gamma_arr[i]);
    }
    //printf("\n\ngamma_arr[73] = %lf\n\n", gamma_arr[73]);

    // make dp perturbation param arrays (initial angles of release)
    double theta_1[THETA_STEPS], theta_2[THETA_STEPS];
    float theta_step = 180.0 / THETA_STEPS;

    for (int i = 0; i < THETA_STEPS; i++) {
        theta_1[i] = i * theta_step;
        theta_2[i] = i * theta_step;
    }

    // lddp initial deviation vector (arbitrary normed vector, 2x1 vector)
    double* d_lddp = (double*)calloc(2, sizeof(double));
    d_lddp[0] = 1.0e-8; //1.0;

    // qddp initial deviation vector (2x1 vector)
    double* d_qddp = (double*)calloc(4, sizeof(double));
    d_qddp[0] = 1e-8;
    
    // dp initial deviation vector (4x1 vector)
    double* d_dp = (double*)calloc(16, sizeof(double));
    d_dp[0] = 1e-8;

    // initial integrator conditions + setup
    double t_init = 0.0;
    double dt_init = 3.8e-7; //0.0001;
    double* t = &t_init;
    double* dt = &dt_init;
    double traj_err, dev_err, err;

    // trajectory step arrays (position, velocity)
    double *lddp_traj = (double*)malloc(2*COMPUTE_STEPS*sizeof(double));
    double *qddp_traj = (double*)malloc(2*COMPUTE_STEPS*sizeof(double));
    double *dp_traj = (double*)malloc(4*COMPUTE_STEPS*sizeof(double));

    // calculate a singular max lyapunov exponent: lddp
    int i = 0;
    traj_step_ddp_t* traj_step = (traj_step_ddp_t*)malloc(sizeof(traj_step_ddp_t));
    dev_step_ddp_t* dev_step = (dev_step_ddp_t*)malloc(sizeof(dev_step_ddp_t));
    double* maxlyp_sum_lddp = (double*)calloc(1, sizeof(double));
    double t_forward = 20.0;
    *dt = 4.8e-05;  // test, delete later
    double gamma_stable = 0.0; //0.9;  // test, delete later
    double gamma_lo = 1.073;  // test, delete later
    double gamma_hi = 1.105;
    printf("\ngamma_arr[73] = %lf\n\n", gamma_arr[73]);
    //*dt = 1e-04;
    while (i < COMPUTE_STEPS) {
        //if (*dt < 1.0e-04) {*dt = 1.0e-04;}
        rk45_lddp_step(traj_step, dev_step, s_lddp, d_lddp, c_lddp, t, dt, gamma_stable); //gamma_arr[]
        err = fmax(traj_step->err, dev_step->err);
        
        if (err < 1.0) {
            // accept step
            *t += (*dt);
            *s_lddp = traj_step->next;
            memcpy(d_lddp, dev_step->next, 2*sizeof(double));

            // save accepted trajectory and increment loop
            lddp_traj[i] = s_lddp->phi;
            lddp_traj[i+1] = s_lddp->omega;
            i += 2;
            // printf("SUCCESS: phi = %lf, omega = %lf, dev = %lf, %lf, dt = %.10e\n", s_lddp->phi, s_lddp->omega, d_lddp[0], d_lddp[1], *dt);
            
            // update accumulated norm (max lyapunov sum), discarding transients (first 10 secs)
            if (*t > t_forward) {
                *maxlyp_sum_lddp += log(dev_step->norm);
            }
            if (i % 10000 == 0){
                printf("gamma = %lf\n\n", gamma_stable);
                printf("dev = %lf, %lf, dt = %.10e\n, t = %lf, dev_step->err = %lf, dev_step->norm = %.10e, maxlyp_sum: %.10e\n", d_lddp[0], d_lddp[1], *dt, *t, dev_step->err, dev_step->norm, *maxlyp_sum_lddp);
            }

        } else {
            printf("max err = %lf\n", err);
        }

        // compute new dt (w/ safety clamp)
        double factor = 0.9 * pow(1.0 / err, EXP);
        if (factor < 0.1) {factor = 0.1;}
        if (factor > 5.0) {factor = 5.0;}
        (*dt) *= factor;
    }

    double maxlyp_lddp = *maxlyp_sum_lddp / (*t - t_forward);
    printf("\nmaxlyp_sum / t = maxlyp_lddp:\n%lf / %lf = %lf\n", *maxlyp_sum_lddp, *t - t_forward, maxlyp_lddp);

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

    free(lddp_jac);
    free(qddp_jac);
    free(dp_jac);

    free(d_lddp);
    free(d_qddp);
    free(d_dp);

    free(lddp_traj);
    free(qddp_traj);
    free(dp_traj);

    free(traj_step);
    free(dev_step);
    free(maxlyp_sum_lddp);
}