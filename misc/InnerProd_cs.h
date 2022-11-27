#ifndef INNER_P
#define INNER_P

#include<stdio.h>
#include"laplacianHelpers.h"

cs *multiply_diag(cs *A, cs*B, int n, int n_complete, int num_entries)
{
    int n = M->nz;
    for (int i = 0; i < num_entries; i++) {
        A->
    }

    return C;
}
cs *dagger(cs *L, unsigned long n, int n_complete, T p)
{
    cs *L_dag = malloc(sizeof(cs));
    assert(L_dag != NULL);
    // doing conjugation requires these 2 other diagonal matrices
    cs *M = malloc(sizeof(cs));
    assert(M != NULL);
    cs *M_inv = malloc(sizeof(cs));
    assert(M_inv != NULL);

    M->m = n;                        M_inv->m = n;
    M->n = n;                        M_inv->n = n;
    M->nz = n;                       M_inv->nz = n;
    M->nzmax = n;                    M_inv->nzmax = n;
    M->x = calloc(n, sizeof(T));     M_inv->x = calloc(n, sizeof(T));
    M->p = calloc(n, sizeof(int));   M_inv->p = calloc(n, sizeof(int));
    M->i = calloc(n, sizeof(int));   M_inv->i = calloc(n, sizeof(int));
    assert(M->x != NULL);            assert(M_inv->x != NULL);
    assert(M->p != NULL);            assert(M_inv->p != NULL);
    assert(M->i != NULL);            assert(M_inv->i != NULL);

    int i = 0;
    int j = 0;
    double val     = (double)n_complete - 1;
    double val_inv = 1 / ((double)n_complete - 1);
    while (i < n_complete) {
        cs_push(M, i, i, val, &i);
        cs_push(M_inv, j, j, val_inv, &j);
    }

    val     = 1 / ((double)1 - p);
    val_inv = (double)1 - p;
    while (i < n) {
        cs_push(M, i, i, val, &i);
        cs_push(M_inv, j, j, val_inv, &j);
    }
    printf("M:\n");
    cs_print(M, 0);
    printf("M_inv:\n");
    cs_print(M_inv, 0);

    L_dag = cs_multiply(L, M_inv);
    cs_print(L_dag, 0);
    printf("NEXT\n");
    //exit(1);
    L_dag = cs_multiply(M, L_dag);
    cs_print(L_dag, 0);
    printf("NEXT\n");

    return L_dag;
}
#endif
