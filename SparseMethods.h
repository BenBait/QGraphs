#ifndef INNER_P
#define INNER_P

#include<stdio.h>
#include<complex.h>

#ifdef I
#undef I
#define img _Complex_I
#endif

#define NZMAX 10000
//typedef double complex T;
typedef double complex T;

typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
    int nzmax ;	    /* maximum number of entries */
    int nz ;	    /* # of entries */
    int m ;	        /* number of rows */
    int n ;	        /* number of columns */
    int *p ;	    /* column indices */
    int *i ;	    /* row indices */
    T *x ;	        /* numerical values, size nzmax */
} cs ;

void cs_push (cs *A, int r, int c, double val, int *idy)
{
    int id = *idy;
    A->i[id] = r;
    A->p[id] = c;
    A->x[id] = val;
    (*idy)++;
}

// multiplying a vector times M just scales the entries
T *multiply_M_vec (T *A,
                   unsigned long n,
                   int n_complete,
                   float val,
                   float val_low)
{
    // the number of entries scaled by the n_complete weight is the
    // top n_complete entries
    int num_top = n_complete;

    T *L = malloc(sizeof(T)*n);
    assert(L != NULL);

    int i = 0;
    for (i = 0; i < num_top; i++) {
        L[i] = A[i]*val;
    }

    while (i < n) {
        L[i] = A[i]*val_low;
        i++;
    }

    return L;
}

// we use the standard inner product after multiplying the left vector by the
// matrix M in the hermitian inner product
T std_inner_prod (T *f,
                  T *g,
                  unsigned long n)
{
    T L = 0;

    int i = 0;
    for (i = 0; i < n; i++) {
        L += f[i]*g[i];
    }
    
    return L;
}

// f and g are vectors that we want to find the hermitian inner product of
T herm_inner_prod (T *f,
                   T *g,
                   unsigned long n,
                   int n_complete,
                   T p)
{
    double M_high = (double)n_complete - 1;
    double M_low = 1 / ((double)1 - p);

    T *L_tmp = multiply_M_vec(g, n, n_complete, M_high, M_low);
    assert(L_tmp != NULL);
        
    double M_inv_high = 1 / ((double)n_complete - 1);
    double M_inv_low = (double)1 - p;

    T L = std_inner_prod(f, L_tmp, n);

    free(L_tmp);

    return L;
}

/*
 * The following functions are for finding the hermitian conjugate of a matrix
 * The hermitian conjugate is the "dagger" function, and it is used in checking
 * that our matrix M does indeed induce a self-adjoint Laplacian in this
 * Hilbert space
 */

// if multiplying with diagonal M on the left, scale entries in column
//cs *multiply_M_mat(cs *A, int n_complete, float val, float val_low)
void multiply_M_mat(cs *A,
                    int n_complete,
                    float val,
                    float val_low)
{
    // the number of entries in the top n_complete rows is n_complete
    // (n_complete-1 connections and the diagonal)
    int num_entries_top = n_complete*n_complete;
    int num_entries = A->nz;

    int i = 0;
    for (i = 0; i < num_entries_top; i++) {
        A->x[i] = 4*val;
    }
    while (i < num_entries) {
        A->x[i] *= val_low;
        i++;
    }
}

// helper transpose function used in calculating dagger
// this function copies the transpose of A in a different matrix
cs *transpose(cs *A) 
{
    cs *L = malloc(sizeof(cs));
    assert(L != NULL);
    L->m = A->m;
    L->n = A->n;
    int num_entries = A->nz;
    L->nz = num_entries;
    L->nzmax = NZMAX;
    L->x = calloc(num_entries, sizeof(T));
    L->p = calloc(num_entries, sizeof(int));
    L->i = calloc(num_entries, sizeof(int));
    assert(L->x != NULL);
    assert(L->p != NULL);
    assert(L->i != NULL);

    for (int i = 0; i < num_entries; i++) {
        L->x[i] = A->x[i];
        L->p[i] = A->i[i];
        L->i[i] = A->p[i];
    }

    return L;
}

cs *dagger(cs *L, unsigned long n, int n_complete, T p)
{
    double M_high = (double)n_complete - 1;
    double M_low = 1 / ((double)1 - p);
    multiply_M_mat(L, n_complete, M_high, M_low);
        
    double M_inv_high = 1 / ((double)n_complete - 1);
    double M_inv_low = (double)1 - p;

    cs *L_dag_tmp = transpose(L);

    multiply_M_mat(L_dag_tmp, n_complete, M_inv_high, M_inv_low);

    cs *L_dag = transpose(L_dag_tmp);

    // free what is inside the matrix first, then the pointer to the matrix
    free(L_dag_tmp->x);
    free(L_dag_tmp->p);
    free(L_dag_tmp->i);
    free(L_dag_tmp);

    return L_dag;
}

// turn a sparse matrix into a full matrix
T *get_dense(cs *A,
             int num_entries,
             unsigned long n)
{
    T *full = calloc(n*n, sizeof(T));

    int r, c;
    T val;
    int full_id;

    int i;
    for (i = 0; i < num_entries; i++) {
        r = A->i[i];
        c = A->p[i];
        val = A->x[i];
        full_id = c + r*n;
        assert(full_id < n*n);
        full[full_id] = val;
    }

    return full;
}

T *sq_scale_and_addI_dense (cs *A,
                           int num_entries,
                           unsigned long n,
                           T scale)
{
    T *out = calloc(n*n, sizeof(T));

    int r, c , i;
    T val;
    int full_id, lhs_id, rhs_id;

    T *full_A = get_dense(A, num_entries, n);

    for (r = 0; r < n; r++) {
        for (c = 0; c < n; c++) {
            full_id = r + c*n;
            // dot_product lhs row with rhs col
            val = 0;
            for (i = 0; i < n; i++) {
                // iterate through lhs row major
                lhs_id = i + r*n;
                // iterate through rhs col major
                rhs_id = c + i*n;
                val += full_A[lhs_id]*full_A[rhs_id];
            }
            //printf("ID: %d\n", full_id);
            if (r == c) val += 1.0;
            out[full_id] = val;
        }
    }

    free (full_A);
    return out;
}

T *mult_dense_from_dense (T *A,
                          T *B,
                          unsigned long n,
                          T scale)
{
    T *out = calloc(n*n, sizeof(T));

    int r, c, i;
    T val;
    int full_id, lhs_id, rhs_id;

    for (r = 0; r < n; r++) {
        for (c = 0; c < n; c++) {
            full_id = r + c*n;
            // dot_product lhs row with rhs col
            val = 0;
            for (i = 0; i < n; i++) {
                // iterate through lhs row major
                lhs_id = i + r*n;
                // iterate through rhs col major
                rhs_id = c + i*n;
                val += A[lhs_id]*B[rhs_id];
            }
            if (r == c) val += 1.0;
            out[full_id] = val / scale;
        }
    }

    return out;
}

T *multiply_mat_vec (T *A,
                     T *s,
                     unsigned long n)
{
    T *out = calloc(n, sizeof(T));

    int r, c;
    T val;
    int lhs_id;

    for (r = 0; r < n; r++) {
        val = 0;
        for (c = 0; c < n; c++) {
            // dot_product lhs row with rhs col
            // iterate through lhs row major
            lhs_id = c + r*n;
            // iterate through vector
            val += A[lhs_id]*s[c];
        }
        out[r] = val;
    }

    return out;
}

void add_mat(T *A, T* B, unsigned long n) {
    int i;
    int n_sqr = n*n;
    for (i = 0; i < n_sqr; i++) {
        A[i] += B[i];
    }
}

void print_dense_from_dense (T *full,
                             unsigned long n)
{
    int r, c;
    int full_id;

    for (r = 0; r < n; r++) {
        for (c = 0; c < n; c++) {
            full_id = c + r*n;
            if (cimag(full[full_id]) != 0)
                printf("%.2f+i*%0.2f\t", creal(full[full_id]), cimag(full[full_id]));
            else
                printf("%.2f\t", creal(full[full_id]));
        }
        printf("\n");
    }
}

void fprint_dense(const char file_loc[],
                  cs *A,
                  int num_entries,
                  unsigned long n)
{
    int r, c;
    T val;
    int full_id;

    T *full = get_dense(A, num_entries, n);

    FILE *F = fopen(file_loc, "w");
    for (r = 0; r < n; r++) {
        for (c = 0; c < n-1; c++) {
            full_id = c + r*n;
            if (cimag(full[full_id]) != 0)
                fprintf(F, "%.2f+i*%0.2f\t", creal(full[full_id]), cimag(full[full_id]));
            else
                fprintf(F, "%.2f\t", creal(full[full_id]));
        }
        // don't print the tab for the last character, just newline
        full_id = (n-1) + r*n;
        if (cimag(full[full_id]) != 0)
            fprintf(F, "%.2f+i*%0.2f\n", creal(full[full_id]), cimag(full[full_id]));
        else
            fprintf(F, "%.2f\n", creal(full[full_id]));
    }

    free(full);
}

void print_dense(cs *A,
                 int num_entries,
                 unsigned long n)
{
    int r, c;
    T val;
    int full_id;

    T *full = get_dense(A, num_entries, n);

    for (r = 0; r < n; r++) {
        for (c = 0; c < n; c++) {
            full_id = c + r*n;
            if (cimag(full[full_id]) != 0)
                printf("%.2f+i*%0.2f\t", creal(full[full_id]), cimag(full[full_id]));
            else
                printf("%.2f\t", creal(full[full_id]));
        }
        printf("\n");
    }

    free(full);
}

#endif
