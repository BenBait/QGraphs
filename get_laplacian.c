#include <stdlib.h>
#include <stdio.h>

typedef int T;

/* sum of a row is the degree of that node in the adjacency matrix, A */
/* A is just one row */
int *get_degree(T *A, int n) {

    // row and column counters
    int r = 0;
    int c = 0;
    // d contains the diagonal elements of the degree matrix
    // i.e. d[i] = D_ii
    int *d = calloc(n, sizeof(T));
    // accumulate the value of d for each row in this var
    int curr_d = 0;
    int idx = 0;

    for (r = 0; r < n; r++) {
        for (c = 0; c < n; c++) {
            // traverse through the array row major
            idx = c + r*n;
            curr_d += A[idx];
        }

        d[r] = curr_d;
        curr_d = 0;
    }

    return d;
}

/* Laplacian is the degree matrix minus the adjacency matrix */
T *get_laplacian(T *A, int n) {

    // get the diagonals of the degree matrix
    int *d = get_degree(A, n);
    
    // row and column indices
    int r = 0;
    int c = 0;

    // index in the adjacency matrix
    int idx = 0;

    // update the Laplacian in A, might be better to have it all in place
    // L = D - A, D is diagonal
    for (r = 0; r < n; r++) {
        for (c = 0; c < n; c++) {
            idx = c + r*n;
            A[idx] = (-1) * A[idx];
        }
        A[r+r*n] += d[r];
    }

    free(d);
    return A;
}

int main(int argc, char *argv[]) {

    // dimension of matrix
    int n = 3;

    // the data in the adjacency matrix
    T *A = calloc(n*n, sizeof(T));

    // get_laplacian turns the adjacency matrix in to the degree matrix,
    // in place
    A = get_laplacian(A, n);

    // row and column counters
    int r = 0;
    int c = 0;

    //counter through adjacency matrix
    int idx = 0;
    for (r = 0; r < n; r++) {
        for (c = 0; c < n; c++) {
            idx = c + r*n;
            printf("%d\t", A[idx]);
        }
        printf("\n");
    }

    free(A);
    return 0;
}
