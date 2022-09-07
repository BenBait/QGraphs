#include <stdlib.h>
#include <stdio.h>

typedef int T;

typedef struct SparseM_entry {
    int r,c;
    T val;
} SparseM_entry;

typedef struct SparseM {
    SparseM_entry *data;
} SparseM;

/* sum of a row is the degree of that node in the adjacency matrix, A */
/* A is just one row */
int *get_degree(SparseM *A, int num_entries, int n) {

    // loop counter
    int i = 0;
    // row and column counters
    int r = 0;
    int c = 0;
    // d contains the diagonal elements of the degree matrix
    // i.e. d[i] = D_ii
    int *d = calloc(n, sizeof(T));

    // accumulate the value of d for each row in this var
    int curr_d = 0;
    int idx = 0;

    for (i = 0; i < num_entries; i++) {
        // get nonzero entries
        r = (A->data[i]).r;
        // don't actually need the column here

        // accumulate the degree for each entry
        d[r] += (A->data[i]).val;
    }

    return d;
}

/* Laplacian is the degree matrix minus the adjacency matrix */
SparseM *get_laplacian(SparseM *A, int num_entries, int n) {

    // loop counter
    int i = 0;
    // get the diagonals of the degree matrix
    int *d = get_degree(A, num_entries, n);
    
    // row and column indices
    int r = 0;
    int c = 0;

    // index in the adjacency matrix
    int idx = 0;

    // update the Laplacian in A, might be better to have it all in place
    // L = D - A, D is diagonal
    // 
    // NEED TO ADD BOUNDS CHECKING: NOT MORE THAN N VALUES FOR ANY GIVEN R
    for (i = 0; i < num_entries; i++) {
        // get nonzero entries
        (A->data[i]).val = (-1) * (A->data[i]).val;

        c = (A->data[i]).c;
        r = (A->data[i]).r;

        if (r == c)
            (A->data[i]).val += d[r];
    }

    free(d);
    return A;
}

int main(int argc, char *argv[]) {

    // loop counter
    int i = 0;
    // dimension of matrix
    int n = 3;
    // number of nonzero entries in the identity matrix
    int num_entries = n;

    SparseM *A = malloc(sizeof(SparseM));
    A->data    = malloc(num_entries*sizeof(SparseM_entry));

    // make identity matrix
    for (int i = 0; i < n; i++) {
        // NEED BOUND CHECKING
        (A->data[i]).r = i;
        (A->data[i]).c = i;
        (A->data[i]).val = 1;
    }

    // get_laplacian turns the adjacency matrix in to the degree matrix,
    // in place
    A = get_laplacian(A, num_entries, n);
    // row and column counters
    int r = 0;
    int c = 0;
    // value at nonzero entry
    int val = 0;

    printf("(r, c, value)\n");

    for (i = 0; i < num_entries; i++) {
        // get nonzero entries
        c = (A->data[i]).c;
        r = (A->data[i]).r;
        val = (A->data[i]).val;
        printf("(%d, %d, %d)\n", r, c, val);
    }
    /*
     * WILL WORK ON PRINTING OUT FULL MATRIX

    //counter through adjacency matrix
    int idx = 0;
    for (r = 0; r < n; r++) {
        for (c = 0; c < n; c++) {
            idx = c + r*n;
            printf("%d\t", A[idx]);
        }
        printf("\n");
    }
    */

    free(A);
    return 0;
}
