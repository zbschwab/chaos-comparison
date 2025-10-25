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
// #include <cblas.h>

#define PORT 8080

#define PI 3.1415927
#define EXP 0.2

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

int write_to_csv(double* test_arr1, double* test_arr2) {
    FILE *test_file;

    test_file = fopen("test_data.csv", "w"); // wx for no overwrite
    if (test_file == NULL) {
        if (errno == EEXIST) {
            perror("Error: file already exists");
            return 1;
        }
    }

    fprintf(test_file, "test_arr1, test_arr2\n");
    for (int i = 0; i < 1000; i++) {
        fprintf(test_file, "%lf, %lf\n", test_arr1[i], test_arr2[i]);
    }
    fclose(test_file);

    return 0;
}

int main(int argc, char const* argv[])
{
    // test writing data to csv
    double* test_arr1 = (double*)malloc(1000 * sizeof(double));
    double* test_arr2 = (double*)malloc(1000 * sizeof(double));
    for (int i = 0; i < 1000; i++) {
        test_arr1[i] = i;
        test_arr2[i] = 2*i;
    }

    if (write_to_csv(test_arr1, test_arr2) == 0) {
        printf("data written to file named `test_data.csv`\n");
    }
    
    // test socket connection to python notebook
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

    size_t num_doubles = sizeof(test_arr1);

    printf("Sending %zu doubles (%zu bytes total)...\n", num_doubles, num_doubles * sizeof(double));

    // Send the message
    if (send_message(client_fd, test_arr1, num_doubles) == 0) {
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
        free(test_arr1);
        free(test_arr2);

        close(client_fd);  // closing the connected socket
        return status;
}