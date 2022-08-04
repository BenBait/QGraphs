#include "stdio.h"
#include "stdlib.h"
#include "assert.h"

typedef struct Matrix {

    // width of matrix
    int N;
    // size of the entries in the matrix
    int val_width;
    // number of non-zero entries in the matrix
    int num_entries;
    // array of 2-entry-arrays 
    int *rc_entries;

} Matrix;

int main(int argc, char *argv[]) {

    if (argc != 4) {
        fprintf(stderr, "Please pass 3 arguments\n");
        return 1;
    }

    // create an empty matrix object
    Matrix M = {0, 0, 0, NULL};

    // initialize the dimension and width of the values for M
    M.N           = atoi(argv[1]);
    M.val_width   = atoi(argv[2]);
    M.num_entries = atoi(argv[3]);

    // 2 entries, row and column, for every entry in M data
    M.rc_entries = calloc(M.num_entries, 2*sizeof(int));
    assert(M.rc_entries != NULL);

    free(M.rc_entries);

    return 0;
}
