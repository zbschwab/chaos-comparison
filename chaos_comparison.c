#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <math.h>
#include <cblas.h>

#include "phys_math.h"
#include "butcher_tableau.h"

#define PORT 8080

#define PI 3.1415927
#define EXP 0.2

#define COMPUTE_STEPS 1e6
#define GAMMA_STEPS 100
#define THETA_STEPS 360 // range is 180, increment by 0.5 deg

// Send a across a socket with a header that includes the message length.
int send_message(int fd, double* message, size_t num_doubles) {
    // If the message is NULL, set errno to EINVAL and return an error
    if (message == NULL) {
        errno = EINVAL;
        return -1;
    }

    // Calculate the total number of bytes to send for the double array
    size_t total_bytes_to_send = num_doubles * sizeof(double);

    // First, send the length of the message in a size_t
    if (write(fd, &total_bytes_to_send, sizeof(size_t)) != sizeof(size_t)) {
        // Writing failed, so return an error
        printf("Writing failed\n");
        return -1;
    }

    // Now we can send the message. Loop until the entire message has been written.
    size_t bytes_written = 0;
    while (bytes_written < total_bytes_to_send) {
        // Try to write the entire remaining message
        ssize_t rc = write(fd, (char*)message + bytes_written, total_bytes_to_send - bytes_written);
        // Did the write fail? If so, return an error
        if (rc <= 0) {
            perror("Writing message data failed");
            return -1;
        }
        // If there was no error, write returned the number of bytes written
        bytes_written += rc;
    }

    return 0;
}

int confirm_sent(int fd) {
    char buffer[16];
    ssize_t n = recv(fd, buffer, sizeof(buffer) - 1, 0);
    if (n <= 0) {
        perror("Confirmation failed");
        return -1;
    }

    buffer[n] = '\0';
    if (strcmp(buffer, "OK") == 0) {
        return 0;  // success
    } else {
        fprintf(stderr, "Unexpected confirmation: %s\n", buffer);
        return -1;
    }
}

int main(int argc, char const* argv[])
{
    // set initial system conditions: lddp
    cons_ddp_t *c_lddp = (cons_ddp_t *)malloc(sizeof(cons_ddp_t));
    state_ddp_t *s_lddp = (state_ddp_t *)malloc(sizeof(state_ddp_t));
    ;
    double *lddp_jac = (double *)malloc(sizeof(double));

    c_lddp->beta = 3.0/4.0;
    c_lddp->omega0 = 3.0 * PI;
    c_lddp->omegad = 2.0 * PI;

    s_lddp->phi = 0.0;
    s_lddp->omega = 0.0;

    // set initial system conditions: qddp
    cons_ddp_t *c_qddp = (cons_ddp_t *)malloc(sizeof(cons_ddp_t));
    state_ddp_t *s_qddp = (state_ddp_t *)malloc(sizeof(state_ddp_t));
    ;
    double *qddp_jac = (double *)malloc(sizeof(double));

    c_qddp->beta = 1.0;
    c_qddp->omega0 = 1.0;
    c_qddp->omegad = 1.0;

    s_qddp->phi = 2.0;
    s_qddp->omega = 2.0;

    // set initial system conditions: dp
    cons_dp_t *c_dp = (cons_dp_t *)malloc(sizeof(cons_dp_t));
    state_dp_t *s_dp = (state_dp_t *)malloc(sizeof(state_dp_t));
    ;
    double *dp_jac = (double *)malloc(sizeof(double));

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
    double gamma_step = (gamma_stop - gamma_start) / (GAMMA_STEPS - 1);

    printf("gamma_arr:\n");
    for (int i = 0; i < GAMMA_STEPS; i++)
    {
        gamma_arr[i] = gamma_start + i * gamma_step;
        // printf(" %lf ", gamma_arr[i]);
    }
    // printf("\n\ngamma_arr[73] = %lf\n\n", gamma_arr[73]);

    // make dp perturbation param arrays (initial angles of release)
    double theta_1[THETA_STEPS], theta_2[THETA_STEPS];
    float theta_step = 180.0 / THETA_STEPS;

    for (int i = 0; i < THETA_STEPS; i++)
    {
        theta_1[i] = i * theta_step;
        theta_2[i] = i * theta_step;
    }

    // lddp initial deviation vector (arbitrary normed vector, 2x1 vector)
    double *d_lddp = (double *)calloc(2, sizeof(double));
    d_lddp[0] = 1.0e-8; // 1.0;

    // qddp initial deviation vector (2x1 vector)
    double *d_qddp = (double *)calloc(4, sizeof(double));
    d_qddp[0] = 1e-8;

    // dp initial deviation vector (4x1 vector)
    double *d_dp = (double *)calloc(16, sizeof(double));
    d_dp[0] = 1e-8;

    // initial integrator conditions + setup
    double t_init = 0.0;
    double dt_init = 1.0e-5; // 3.8e-7; //0.0001;
    double *t = &t_init;
    double *dt = &dt_init;
    double traj_err, dev_err, err;

    // trajectory step arrays (position, velocity)
    double *lddp_traj = (double *)malloc(2 * COMPUTE_STEPS * sizeof(double));
    double *qddp_traj = (double *)malloc(2 * COMPUTE_STEPS * sizeof(double));
    double *dp_traj = (double *)malloc(4 * COMPUTE_STEPS * sizeof(double));

    // calculate a singular max lyapunov exponent: lddp
    int i = 0;
    traj_step_ddp_t *traj_step = (traj_step_ddp_t *)malloc(sizeof(traj_step_ddp_t));
    dev_step_ddp_t *dev_step = (dev_step_ddp_t *)malloc(sizeof(dev_step_ddp_t));
    double *maxlyp_sum_lddp = (double *)calloc(1, sizeof(double));
    double t_forward = 20.0;
    //*dt = 4.8e-05;  // test, delete later
    double gamma_stable = 0.9; // 0.9;  // test, delete later
    double gamma_lo = 1.073;   // test, delete later
    double gamma_hi = 1.105;
    printf("\ngamma_arr[73] = %lf\n\n", gamma_arr[73]);
    //*dt = 1e-04;
    while (i < COMPUTE_STEPS)
    {
        rk45_lddp_step(traj_step, dev_step, s_lddp, d_lddp, c_lddp, t, dt, gamma_lo); // gamma_arr[]
        err = fmax(traj_step->err, dev_step->err);

        if (err < 1.0)
        {
            // accept step
            *t += (*dt);
            *s_lddp = traj_step->next;
            memcpy(d_lddp, dev_step->next, 2 * sizeof(double));

            // save accepted trajectory and increment loop
            lddp_traj[i] = s_lddp->phi;
            lddp_traj[i + 1] = s_lddp->omega;
            i += 2;

            // update accumulated norm (max lyapunov sum), discarding transients (first 10 secs)
            if (*t > t_forward) {
                *maxlyp_sum_lddp += log(dev_step->norm);
            } if (i % 10000 == 0) {
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

    int client_fd, status = -1;
    struct sockaddr_in server_addr;

    // create socket
    if ((client_fd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
        perror("Socket creation failed");
        return 1;
    }

    server_addr.sin_family = AF_INET;
    server_addr.sin_port = htons(PORT);

    if (inet_pton(AF_INET, "127.0.0.1", &server_addr.sin_addr) <= 0) {
        perror("Invalid address/ Address not supported");
        goto cleanup;
    }

    // connect to server
    if (connect(client_fd, (struct sockaddr *)&server_addr, sizeof(server_addr)) < 0) {
        perror("Connection Failed");
        goto cleanup;
    }

    // Example double array to send
    double data_to_send[] = {1.23, 4.56, 7.89, 10.11, 12.13};
    size_t num_elements = sizeof(data_to_send) / sizeof(data_to_send[0]);

    printf("Sending %zu doubles (%zu bytes total)...\n", num_elements, num_elements * sizeof(double));

    // Send the message
    if (send_message(client_fd, data_to_send, num_elements) == 0) {
        printf("Data sent successfully.\n");
        // wait for server to confirm all data was successfully
        if (confirm_sent(client_fd) == 0) {
            // let the server know you're done writing
            shutdown(client_fd, SHUT_WR);
            status = 0;
        } else {
            printf("Failed to confirm data was recieved.\n");
            goto cleanup;
        }
    } else {
        printf("Failed to send data.\n");
    }

    cleanup:
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

        close(client_fd);  // closing the connected socket
        return status;
}