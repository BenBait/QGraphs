#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
// functions for generating matrices
#include "LaplacianHelpers.h"
#include "SparseMethods.h"
#include "exp/matrix_exponential.h"

T get_success_probability (cs *A, // adjacency matrix for segmented graph
                           int num_entries, // entries in segmented adj matrix
                           int n_complete,  // nodes in complete graph
                           unsigned long n, // nodes in segmented graph
                           int w,   // target node
                           T p,     // segmented edge node transition probability
                           int t, // time value to checkout
                           double gamma); // gamma coefficient

// helper function to test with the identity matrix
cs *get_diag (int n, int diag);

// helper function computes exponential of Hamiltonian with a Taylor expansion
T *get_exp_H (cs *H, int num_H_entries, int n, int Pr);
T *diag_get_exp_H (cs *H, int num_H_entries, int n, int Pr);

int main (int argc, char *argv[])
{
    // n_complete is the number of nodes before segmentation
    if (argc < 5) {
        const char *err_message = "PASS 5 ARGUMENTS: n, p, w, t, gamma\n";
        fprintf(stderr, "%s", err_message); 
        exit(1);
    }

    // number of nodes before segmenting
    int n_complete = atoi(argv[1]);
    // probability parameter p (was 1/2 initially)
    T p = atof(argv[2]);
    // index we are searching for
    int w = atoi(argv[3]);
    // time value
    int t = atoi(argv[4]);
    // gamma scale factor of Laplacian for the Hamiltionian
    double gamma= atof(argv[5]);
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

    // original entries have (n_complete - 1) adjacencies, 
    // nodes in segments each have 2 adjacencies
    int num_entries = n_complete*(n_complete - 1) + (n-n_complete)*2;

    cs *A = get_segmented_complete(num_entries, n, n_complete, p, segs);

    get_success_probability(A, num_entries, n_complete, n, w, p, t, gamma);
    return 0;

}

T get_success_probability (cs *A,
                           int num_entries, // entries in segmented adj matrix
                           int n_complete,  // nodes in complete graph
                           unsigned long n, // nodes in segmented graph
                           int w,   // target node
                           T p,     // segmented edge node transition probability
                           int t, // maximum time value to checkout
                           double gamma) // gamma coefficient
{
    int i; // general iterator

    // entries in the Hamiltionian is equal to entries in the adjacecncy matrix
    // plus the diagonal (which is all 0 in A);
    // this is calculated in get_hamiltonian, we just need it for printing

    int num_H_entries = num_entries + n;

    // we need the following values for every calculation of the success prob

    // we already know what norm_factor is, see bottom of this function to see
    // how it is calculated
    double M_high = (double)n_complete - 1;
    double M_low = 1 / (1.0 - p);

    T norm_factor = (w < n_complete) ? M_high : M_low;

    // the basis vector for the desired entry
    T *e_w = calloc(n, sizeof(T));
    e_w[w] = 1 / norm_factor;
    // the volume of the graph
    T vol  = M_high*n_complete + M_low*(n - n_complete);

    // uniform superposition state
    T *s   = calloc(n, sizeof(T));
    for (i = 0; i < n; i++)
        s[i] = 1.0 / sqrt(vol);

    cs *L = get_laplacian(A, num_H_entries, n);
    fprint_dense("laplacian.matrix", L, num_H_entries, n);
    // iterate through the variation term in Hamiltonian
    cs *H = get_hamiltonian(A, n, w, (T)gamma);
    fprint_dense("hamiltonian.matrix", H, num_H_entries, n);
    //printf("eigenvectors of the Hamiltonian:\t\n"
    T scale = img*(T)t;

    for (i = 0; i < num_H_entries; i++)
        H->x[i] = H->x[i] * scale;

    // need to replace this with a diagonalization technique
    //T *exp_H = get_exp_H (H, num_H_entries, n, Pr);
    // use an external library first to get something to compare to
    // Gamal's results in python
    T *exp_H = c8mat_expm1(n, get_dense(H, num_H_entries, n));

    // multiply H by the uniform distribution state
    T *exp_H_s = multiply_mat_vec (exp_H, s, n);
    // std basis state inner prod exp(H_s)
    T std_prod_exp = herm_inner_prod (e_w, exp_H_s, n, n_complete, p);
    // no need for conjugate because we just get the square of the norm
    double real_part = creal(std_prod_exp) * creal(std_prod_exp);
    double imag_part = cimag(std_prod_exp) * cimag(std_prod_exp);
    double success_prob = real_part + imag_part;
    // instead of printing it, write it to a file with the time
    //FILE *F = fopen("single_success_prob.dat", "w");
    //fprintf(F, "%d, %f, %f\n", n_complete, w, p, gamma, t, success_prob);
    // print in Latex table format
    printf("%d & %d & %f & %f & %d & %f \\\\ \n", n_complete, w, creal(p), gamma, t, success_prob);
    //fclose(F);

    free(H->i);
    free(H->p);
    free(H->x);
    free(H);

    free(A->i);
    free(A->p);
    free(A->x);
    free(A);

    free(s);
    free(e_w);
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

/* this method of getting the exponential of H is with a diagonalization */
T *diag_get_exp_H (cs *H,
                  int num_H_entries,
                  int n,
                  int Pr)
{
    T *H_dense = get_dense (H, num_H_entries, n);
    T *exp_H = calloc (n*n, sizeof(T));

    if (Pr) {
        printf ("exp(H):\n");
        print_dense_from_dense (exp_H, n);
    }

    return exp_H;
}
/* this method of getting the exponential of H is with a Taylor expansion */
/*
T *get_exp_H (cs *H,
              int num_H_entries,
              int n,
              int Pr)
{
    // k_fact is denominator in taylor expansion
    T k_fact = 2;
    // prod is numerator in taylor expansion
    T *prod = sq_scale_and_addI_dense (H, num_H_entries, n, k_fact);
    T *H_dense = get_dense(H, num_H_entries, n);
    T *exp_H = calloc(n*n, sizeof(T));
    add_mat(exp_H, prod, n);

    int accuracy = 3;
    int i;
    for (i = 0; i < accuracy; i++) {
        k_fact *= k_fact+1;
        T *next_prod = mult_dense_from_dense (prod, H_dense, n, k_fact);
        add_mat (exp_H, next_prod, n);
        prod = next_prod;
    }
    if (Pr) {
        printf("exp(H):\n");
        print_dense_from_dense(exp_H, n);
    }

    return exp_H;
}
*/
