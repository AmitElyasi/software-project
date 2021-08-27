#include <stdio.h> 


/*utils for python c api*/
int pyMat_to_C_array(PyObject* pyMat, float* mat, int dim);

/*return pyList ocject in the shape of (n,m) */
PyObject* c_array_to_pyMat(double* mat, int n, int m);

/*start of functions for Normalized Spectral Clustering implementation*/

/* given A real symmetric matrix updates A and a matrix V that A is diagonal and all of A eigenvalues
on the diagonal and V columns is the correspondence eigenvector
*/
void jacobi_algorithm_for_eigenvalues(float *A, float *V, int n);

void form_diagonal_mat(float *diagonal_mat, float *weighted_adj_mat , int n);

void form_weighted_adj_mat(float *mat, float *data_points, int dim, int n);

/** given D=diag(diagonal_mat) and W=weights_mat
 * calculates I-DWD into normalized_laplacian */
void calc_normalized_laplacian(float *normalized_laplacian, float *diagonal_mat, float *weights_mat, int dim);

/** given a square matrix (mat), sets row and col to be
 * the indexes of the off-diagonal largest absolute value */
void indexes_of_max_off_diag(float *mat, int *row, int *col, int dim);

/** given a square matrix (mat) ,
 * the indexes of its off-diagonal largest absolute value (i,j)
 * and it's dimension (dim)
 * returns t
 */
float get_t(float *mat, int i, int j, int dim);

/** given t, returns c
 */
float get_c(float t);

/** fill given matrix with the data from input file */
void read_data(FILE* fp, float *data_points, char *line, int n, int dim);

/*
* calculate distance between two vectors sqrt(sum_(i=1)^n (v1_i-v2_i)^2) (norm2)
*/
float distance(float *v1, float *v2, int dim,int row_v1, int row_v2);

float get(float* arr, int i, int j, int dim);

void set(float* arr, int i, int j, int dim, float item);

