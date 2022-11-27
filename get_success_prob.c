#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
// functions for generating matrices
#include "LaplacianHelpers.h"
#include "SparseMethods.h"

T get_success_probability (cs *A, // adjacency matrix for segmented graph
                           int num_entries, // entries in segmented adj matrix
                           int n_complete,  // nodes in complete graph
                           unsigned long n, // nodes in segmented graph
                           int w,   // target node
                           T p,     // segmented edge node transition probability
                           T gamma, // variation term in Hamiltonian
                           T t,     // time for exp(iHt) calculation
                           int Pr);

// helper function to test with the identity matrix
cs *get_diag (int n, int diag);

int main (int argc, char *argv[]) {

    // n_complete is the number of nodes before segmentation
    if (argc < 5) {
        fprintf(stderr, "PASS 5 ARGUMENTS: n, p, gamma, w, and Print(0 or 1)\n"); 
        exit(1);
    }

    // number of nodes before segmenting
    int n_complete = atoi(argv[1]);
    // probability parameter p (was 1/2 initially)
    T p = atof(argv[2]);
    // scale factor of Laplacian for the Hamiltionian
    T gamma = atof(argv[3]);
    // index we are searching for
    int w = atoi(argv[4]);
    // if Pr, print
    int Pr = atoi(argv[5]);
    // segs in the number of segments in each edge after decomposing
    int segs = 3;

    // we have to get the number of paths through the inside of the graph
    // the first 2 nodes have n_complete - 2 connections
    // the next n_complete nodes have n_complete - m - 2
    int num_inside = 2*(n_complete - segs);
    // now there are n_complete - 2 more nodes, each with 1 fewer new path
    int added_edges = n_complete - segs - 1;
    // added_edges is 0 
    while (added_edges > 0) {
        num_inside += added_edges;
        added_edges--;
    }
    // n in the adjacency matrix for the segmented matrix
    unsigned long n = n_complete*segs + num_inside*(segs - 1);
    printf("n in semented graph: %ld\n", n);

    // original entries have (n_complete - 1) adjacencies, 
    // nodes in segments each have 2 adjacencies
    int num_entries = n_complete*(n_complete - 1) + (n-n_complete)*2;

    cs *A = get_segmented_complete(num_entries, n, n_complete, p, segs);

    if (Pr) {
        printf("A:\n");
        print_dense(A, num_entries, n);
    }

    T t = 1;
    get_success_probability(A, num_entries, n_complete, n, w, p, gamma, t, Pr);
    return 0;

    // factor that normalizes the dirac delta can be calculated as follows
    /*
    T *f = calloc(n, sizeof(T));
    T *g = calloc(n, sizeof(T));
    f[w] = 1.0;
    g[w] = 1.0;
    T norm_factor = herm_inner_prod(f, g, n, n_complete, p);
    free(f);
    free(g);
    */
}

T get_success_probability (cs *A,
                           int num_entries, // entries in segmented adj matrix
                           int n_complete,  // nodes in complete graph
                           unsigned long n, // nodes in segmented graph
                           int w,   // target node
                           T p,     // segmented edge node transition probability
                           T gamma, // variation term in Hamiltonian
                           T t,     // time for exp(iHt) calculation
                           int Pr) 
{
    cs *H = get_hamiltonian(A, n, w, gamma);

    // entries in the Hamiltionian is equal to entries in the adjacecncy matrix
    // plus the diagonal (which is all 0 in A)
    // don't worry, this is calculated in get_hamiltonian, we just need it for printing

    int num_H_entries = num_entries + n;
    if (Pr) {
        printf("\nH:\n");
        print_dense(H, num_H_entries, n);
    }

    // NOTE: Hamiltoian has 1 less entry once oracle matrix is subtracted
    // multiply every entry in the Hamiltonian by i and the time
    // before taking the exponential
    double complex a = 0.0+1.0*img;
    double complex b = 1.0+0.0*img;
    //double complex prodab = a*b;
    printf("mat_val: %f + i*%f\n", creal(H->x[0]), cimag(H->x[0]));
    printf("PRINTING MAT\n");
    for (int i = 0; i < num_H_entries; i++) {
        //H->x[i] *= just_i*t;
        double complex prod = a*H->x[i];
        printf("prod: %f + i*%f\n", creal(prod), cimag(prod));
        //H->x[i] = a*H->x[i];
        H->x[i] = prod;
        printf("mat_val: %f + i*%f\n", creal(H->x[i]), cimag(H->x[i]));
    }

    if (Pr) {
        printf("\nH:\n");
        print_dense(H, num_H_entries, n);
    }

    // we already know what norm_factor is, see bottom of this function to see
    // how it is calculated
    double M_high = (double)n_complete - 1;
    double M_low = 1 / (1.0 - p);

    T norm_factor = (w < n_complete) ? M_high : M_low;
    printf("Herm Inner Prod : %f + i*%f\n", creal(norm_factor), cimag(norm_factor));

    T k_fact = 2;
    T *prod = sq_scale_and_addI_dense (H, num_H_entries, n, k_fact);
    T *H_dense = get_dense(H, num_H_entries, n);
    T *exp_H = calloc(n*n, sizeof(T));
    add_mat(exp_H, prod, n);

    int accuracy = 3;
    int i;
    for (i = 0; i < accuracy; i++) {
        k_fact *= k_fact+1;
        T *next_prod = mult_dense_from_dense (prod, H_dense, n, k_fact);
        add_mat(exp_H, next_prod, n);
        prod = next_prod;
    }
    if (Pr) {
        printf("exp(H):\n");
        print_dense_from_dense(exp_H, n);
    }
    T *s   = calloc(n, sizeof(T));
    T *e_w = calloc(n, sizeof(T));
    e_w[w] = 1 / norm_factor;
    T vol  = M_high*n_complete + M_low*(n - n_complete);
    for (i = 0; i < n; i++) {
        s[i] = 1.0 / sqrt(vol);
    }
    T *exp_H_s = multiply_mat_vec(exp_H, s, n);
    T success_prob = herm_inner_prod(e_w, exp_H_s, n, n_complete, p);

    printf("\nSuccess Prob: %f + i*%f\n", creal(success_prob), cimag(success_prob));

    free(A->i);
    free(A->p);
    free(A->x);
    free(A);

    free(H->i);
    free(H->p);
    free(H->x);
    free(H);

    free(s);
    free(e_w);
}

cs *get_diag (int n, int diag)
{
    cs *A = malloc(sizeof(cs));
    assert(A != NULL);
    // square matrix; m = n = number of columns
    A->m  = n;
    A->n  = n;
    A->nz = n;
    A->nzmax = NZMAX;
    A->x = calloc(n, sizeof(T));
    A->p = calloc(n, sizeof(int));
    A->i = calloc(n, sizeof(int));
    assert(A->x != NULL);
    assert(A->p != NULL);
    assert(A->i != NULL);

    int id = 0;
    int j;
    for (j = 0; j < n; j++) {
        cs_push(A, j, j, diag, &id);
    }
    return A;
}
