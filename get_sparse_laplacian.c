#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
// functions for generating matrices
#include "LaplacianHelpers.h"
#include "SparseMethods.h"

int main (int argc, char *argv[]) {

    // n_complete is the number of nodes before segmentation
    if (argc < 3) {
        fprintf(stderr, "PASS 3 ARGUMENTS: n, p, and Print(0 or 1)\n");
        exit(-1);
    }

    // number of nodes before segmenting
    int n_complete = atoi(argv[1]);
    // probability parameter p (was 1/2 initially)
    T p = atof(argv[2]);
    // if Pr, print
    int Pr = atoi(argv[3]);
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
    printf("n in segmented graph: %ld\n", n);

    // original entries have (n_complete - 1) adjacencies, 
    // nodes in segments each have 2 adjacencies
    int num_entries = n_complete*(n_complete - 1) + (n-n_complete)*2;

    cs *A = get_segmented_complete(num_entries, n, n_complete, p, segs);

    if (Pr) {
        printf("A:\n");
        print_dense(A, num_entries, n);
    }
    //printf("Checking if rows all sum to one...\n");
    //dumb_row_sum_SparseM(A, num_entries, n, Pr);
    //printf("looks good\n");

    // get_laplacian turns the adjacency matrix in to the degree matrix
    cs *L = get_laplacian(A, num_entries, n);
    free(A->i);
    free(A->p);
    free(A->x);
    free(A);
    // entries in the Laplacian is equal to entries in the adjacecncy matrix
    // plus the diagonal (which is all 0 in A)
    // don't worry, this is calculated in get_laplacian, we just need it for printing
    int num_L_entries = num_entries + n;

    if (Pr) {
        printf("\nL:\n");
        print_dense(L, num_L_entries, n);
    }

    cs *L_dag = dagger(L, n, n_complete, p);

    if (Pr) {
        printf("\nL_dag:\n");
        print_dense(L_dag, num_L_entries, n);
    }

    free(L->i);
    free(L->p);
    free(L->x);
    free(L);
    free(L_dag->i);
    free(L_dag->p);
    free(L_dag->x);
    free(L_dag);
    return 0;
}
