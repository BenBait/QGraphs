#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>
// functions for generating matrices
//#include "LaplacianHelpers.h"
#include "SparseMethods.h"

#define EPSILON 0.00000000001

T check_DFT_eig(T *L, int n, int Pr, int id);
T check_simple_eig(T *L, int n, int Pr, int id);

int main (int argc, char *argv[]) 
{
    int i = 0; // iterator
    // n_complete is the number of nodes before segmentation
    if (argc < 3) {
        fprintf(stderr, "PASS 3 ARGUMENTS: n, Print(0 or 1), and DFT eig(0 or 1)\n");
        exit(-1);
    }

    // number of nodes before segmenting
    int n = atoi(argv[1]);
    // if Pr, print
    int Pr  = atoi(argv[2]);
    int DFT = atoi(argv[3]);

    // the n by n probabilistic laplacian
    T *L = calloc(n*n, sizeof(*L));

    // p is the transition probability for all nodes; note that the probability
    // of going to any other node is equal
    double p = 1 / ((double)n - 1);
    int row = 0;
    // if we are at a diagonal, then the value is 1, otherwise it is p
    // note that i is iterated before row
    for (i = 0; i < n*n; i++, row += (i % n == 0) ? 1 : 0)
        L[i] = (i % n == row) ? 1 : -1*p;

    if (Pr) {
        printf("L:\n");
        print_dense_from_dense(L, n);
    }

    T eig_val = 0;
    // now set the eigenvector depending on the command line args
    if (DFT) { // get the DFT modes as the eigenvector
        // there should be n eigenvectors given by check_DFT_eig
        for (i = 0; i < n; i++) {
            eig_val = check_DFT_eig(L, n, 0, i);
            if (Pr)
                printf("eigenvalue(%d): %f + i%f\n", i, creal(eig_val), cimag(eig_val));
        }
    } else { // simple eigenvectors with 1, -1 moving down the vector
        for (i = 0; i < n; i++) {
            eig_val = check_simple_eig(L, n, 0, i);
            if (Pr)
                printf("eigenvalue(%d): %f + i%f\n", i, creal(eig_val), cimag(eig_val));
        }
    }

    free(L);
    return 0;
}

// pass the id to enumerate the eigenvectors and also raise omega to that power
T check_DFT_eig(T *L, int n, int Pr, int id) 
{
    int i = 0; // iterator

    // allocate space for the eigenvector (use the heap just in case)
    T *v = calloc(n, sizeof(*v));
    // we will need the nth roots of unity for this eig
    T raise_e_to_this = 2*M_PI*img*id / (double)n; 
    T omega = cexp(raise_e_to_this);
    // the norm of the eigenvalue must be one, so normalize with this value
    T p = 1 / sqrt(n);
    // now get the normalized eigenvector

    for (i = 0; i < n; i++)
        v[i] = p*cpow(omega, i);
    // multiply the Laplacian with v to see if it is an eigenvector
    T *prod = multiply_mat_vec(L, v, n);

    int j = 0;
    // go to a nonzero entry or the end of the vector
    while ((creal(prod[j]) - 0.0 < EPSILON && cimag(prod[j]) - 0.0 < EPSILON)
           && j < n) {
        j++;
    }
    // if all the entries in prod are 0, 0 is the eigenvalue
    if (j == n) return 0.0;

    T eig_val = v[j] / prod[j];

    T curr_eig;

    for (i = j+1; i < n; i++) {
        curr_eig = v[i] / prod[i];
        if (creal(eig_val) - creal(curr_eig) > EPSILON ||
            cimag(eig_val) - cimag(curr_eig) > EPSILON) {
            fprintf(stderr, "DFT eigenvector failure\n");
            exit(-1);
        }
    }

    if (Pr) {
        printf("v\n(");
        for (i = 0; i < n; i++)
            printf("%f + i%f\t", creal(v[i]), cimag(v[i]));

        printf(")\nproduct\n(");
        for (i = 0; i < n; i++)
            printf("%f + i%f\t", creal(prod[i]), cimag(prod[i]));
        printf(")\n");
    }
    free(v);
    free(prod);
    
    return eig_val;
}

// pass the id to enumerate the eigenvectors and mark the location of the first
// 1 (followed by -1)
T check_simple_eig(T *L, int n, int Pr, int id) 
{
    int i = 0; // iterator

    // allocate space for the eigenvector (use the heap just in case)
    T *v = calloc(n, sizeof(*v));
    // normalize with p
    T p;
    // now get the normalized eigenvector
    if (id == 0) {
        p = 1 / sqrt(n);
        for (i = 0; i < n; i++)
            v[i] = p;
    } else {
        p = 1 / sqrt(2);
        v[id-1] = p;
        v[id] = -1*p;
    }
    // multiply the Laplacian with v to see if it is an eigenvector
    T *prod = multiply_mat_vec(L, v, n);

    int j = 0;
    // go to a nonzero entry or the end of the vector
    while ((creal(prod[j]) - 0.0 < EPSILON && cimag(prod[j]) - 0.0 < EPSILON)
           && j < n) {
        j++;
    }
    // if all the entries in prod are 0, 0 is the eigenvalue
    if (j == n) return 0.0;

    T eig_val = v[j] / prod[j];

    T curr_eig;

    for (i = j+1; i < n; i++) {
        curr_eig = v[i] / prod[i];
        if (creal(eig_val) - creal(curr_eig) > EPSILON ||
            cimag(eig_val) - cimag(curr_eig) > EPSILON) {
            fprintf(stderr, "simple eigenvector failure\n");
            exit(-1);
        }
    }

    if (Pr) {
        printf("v\n(");
        for (i = 0; i < n; i++)
            printf("%f + i%f\t", creal(v[i]), cimag(v[i]));

        printf(")\nproduct\n(");
        for (i = 0; i < n; i++)
            printf("%f + i%f\t", creal(prod[i]), cimag(prod[i]));
        printf(")\n");
    }
    free(v);
    free(prod);
    
    return eig_val;
}
