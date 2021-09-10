#include <stdio.h> 

extern const int MAX_ITER;

/**
 * find the eigengap and return the number of cluster
 */
int calc_eigenvalue_gap(double *mat, double *sorted_eigenvalues, int n);

/**
 * given A real symmetric matrix A and a unit matrix V, updates A and V
 * with A to be an eigenvalues diagonal matrix,
 * and V's columns to be the corresponding eigenvectors.
 */
void jacobi_algorithm_for_eigenvalues(double *A, double *V, int n);

/**
 * form a unit matrix in mat
 * */
void form_unit_matrix(double *mat, int n);

/**
 * given the weighted adjacency matrix,
 * form the diagonal degree matrix into diagonal_mat
 */
void form_diagonal_mat(double *diagonal_mat, double *weighted_adj_mat , int n);

/**
 * given data points, form it's weighted adjacency matrix into mat
 */
void form_weighted_adj_mat(double *mat, double *data_points, int dim, int n);

/**
 * given D=diag(diagonal_mat) and W=weights_mat
 * calculates I-DWD into normalized_laplacian
 */
void calc_normalized_laplacian(double *normalized_laplacian, double *diagonal_mat, double *weights_mat, int dim);

/**
 * given a square matrix (mat), sets row and col to be
 * the indexes of the off-diagonal largest absolute value
 * */
void indexes_of_max_off_diag(double *mat, int *row, int *col, int dim);


/**
 * given a square matrix (mat) ,
 * the indexes of its off-diagonal largest absolute value (i,j)
 * and it's dimension (dim)
 * returns t
 */
double get_t(double *mat, int i, int j, int dim);

/**
 * given t, returns c
 */
double get_c(double t);

/**
 * given matrix mat calculates the off(mat)^2 = notm_f(mat) - sum_(i=1)^n a_ii
 */
double calc_off_square(double *mat,int n);

/**
 * given diagonal matrix [n*n] and an empty array [n], fills the array with the diagonal values
 */
void diag_to_array(double *diag, double *arr, int n);

/**
 * given array of length dim, and its dimention, fills the array with integers
 * such that arr[i]=i
 */
void fill_ordered_ints(int *arr, int dim);

/**
 * given an array- arr, and another array- indexes, both of length dim,
 * fills indexes with integers [0,dim-1] such that their order dictates a sorted order on arr values.
 * (using stable QuickSort algorithm)
 */
void quickSort_indexes(double *arr, int *indexes, int dim);

/** fill given matrix with the data from input file */
void read_data(FILE* fp, double *data_points, char *line, int n, int dim);

/**
 * Calculate distance between two vectors sqrt(sum_(i=1)^n (v1_i-v2_i)^2) (L2 norm)
 */
double distance(double *v1, double *v2, int dim,int row_v1, int row_v2);

/**
 * Given the eigenvectors matrix [n*n] (as columns),
 * an array [n] that dictates an order on the eigenvalues,
 * n (the mat dimention), k, and a result matrix [n*k],
 * fills the result matrix with the normalized first k eigenvectors as columns
 */
void normalized_k_eigenvectors(double *eigenvecs, int *indexes, int n, int k, double *result);

/**
 * returns arr[i][j]
 */
double get(double* arr, int i, int j, int dim);

/**
 * Sets arr[j][i]=item
 */
void set(double* arr, int i, int j, int dim, double item);

/**
 * Return target matrix for goal == "wam"
 */
double *wam(double *data_points, int n, int dim);

/**
 * Return an array that contains the diagonal of the target matrix for goal == "ddg"
 */
double *ddg(double *data_points,int n, int dim);

/**
 * Return target matrix for goal == "lnorm"
 */
double *lnorm(double *data_points, int n, int dim);

/**
 * Return target matrix in shape (n+1,n)
 * the first row of the returned matrix will be the eigenvalues,
 * the other n rows will be the eigenvectors matrix when the goal == "jacobi"
 * */
double *jacobi(double *data_points, int n);

/**
 * Return target matrix T in shape (n,k) when goal == "spk"
 */
double *spk(double *data_points, int n , int dim, int *k);

/**
 * K-means algorithm implementation
 */
int kmeans(int k, double* data_points, double* centroids, double* utl, int dim, int n);

/**
 * Assigns data points to their closest cluster (measure distance from the centroid)
 * updates the number of cluster for each data point
 */
void assign(double* data_points, double* clusters, int dim, int n, int k);


/**
 * Re-estimates a centroid for each cluster:
 * for each cluster calculate the average of the points assign to it,
 * updates centroids to be the average vector,
 * returns 1 if the old centroids are equal to the new ones.
 */
short re_estimate(double* data_points, double* clusters,double *utl, int dim, int n, int k) ;

/**
 * Adds vec2 to vec1 coordinate wise
 */
void vec_sum(double* vec1, double* vec2, int dim, int row_vec1, int row_vec2);

/**
 * Zeros a given matrix
 */
void zero_mat(double* clusters , int dim, int n);

/**
 * Prints matrix in the template required
 */
void print_matrix(double * mat, int row, int col);

typedef enum {WAM, DDG, LNORM, JACOBI, SPK} Goal;

/**
 * Prints matrix in the template required
 */
Goal get_enum(char* goal_string);

/**
 * Copies the fist last_row_to_copy rows of mat_to_copy to result_mat
 */
void copy_rows(double*, double*, int, int, int);