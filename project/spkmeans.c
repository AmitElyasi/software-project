
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Python.h>

/*utils for python c api*/
int pyMat_to_C_array(PyObject* pyMat, float* mat, int dim){
    int i,j,m,n;
    PyObject* pyVec;
    PyObject* pyItem;
    /* Is it a list? */
    if (!PyList_Check(pyMat))
        return 0;
    /* Get the size of it and build the output list */
    n = PyList_Size(pyMat);  /*  Same as in Python len(_list)  */
    /* Go over each item of the list and reduce it */
    for (i = 0; i < n; i++) {
        pyVec = PyList_GetItem(pyMat, i);
        if (!PyList_Check(pyVec)){  /* We only print lists */
            return 0;
        }
        m = PyList_Size(pyVec);
        for (j = 0; j < m; j++) {
            pyItem = PyList_GetItem(pyVec, j);
            set(mat, i, j , dim, PyFloat_AsDouble(pyItem));
    
            if (get(mat, i, j, dim) == -1 && PyErr_Occurred()){
                /* float too big to fit in a C double, bail out */
                puts("Error parsing a list to C matrix");
                return 0;
            }
        }
    }
    return 1;
}
/*return pyList ocject in the shape of (n,m) */
PyObject* c_array_to_pyMat(double* mat, int n, int m){
    int i, j;
    PyObject *pyItem, *pyVec, *pyMat;
    pyMat = PyList_New(0);
    for (i=0; i < n; i++){
        pyVec = PyList_New(0);
        for (j = 0; j < m; j++){
            pyItem = Py_BuildValue("d", get(mat, i, j, m));
            PyList_Append(pyVec, pyItem);
        }
        PyList_Append(pyMat, pyVec);
    }
    return pyMat;
}

/*start of Normalized Spectral Clustering implementation*/

/* given A real symmetric matrix updates A and a matrix V that A is diagonal and all of A eigenvalues
on the diagonal and V columns is the correspondence eigenvector
*/
void jacobi_algorithm_for_eigenvalues(float *A, float *V, int n){
    float off_of_A,off_of_Atag,epsilon,s,c,t,val;
    int row,col,i;
    /*if A diagonl metrix we skip the while loop with this starting settings*/
    off_of_A = 0;  
    off_of_Atag = calc_off_square(A, n);
    epsilon = 0.001;


    while(fabs(off_of_A - off_of_Atag) > epsilon){
        off_of_A = off_of_Atag;
        indexes_of_max_off_diag(A, &row, &col, n);
        t = get_t(A, row, col, n);
        c = get_c(t);
        s = t * c;
        
        /*update A */
        for(i=0;i<n;i++){
            if(i != row && i != col){
                val = c * get(A, i, row , n) - s * get(A, i, col, n);
                set(A, i, row, n, val);
                set(A, row, i, n, val);
                val = c * get(A, i, col , n) + s * get(A, i, row, n);
                set(A, i, col, n, val);
                set(A, col, i, n, val);
            }else{
                set(A, row, row, n, powf(c,2)* get(A, i, i, n) + powf(s,2) * get(A, col,  col, n) - 2 * s * c * get(A, row, col, n));
                set(A, col, col, n, powf(s,2)* get(A, i, i, n) + powf(c,2) * get(A, col, col, n) + 2 * s * c * get(A, row, col, n));
                set(A, row, col, n, (powf(c,2)*  - powf(s,2)) * get(A, row,  col, n) -  s * c * (get(A, row, row, n) - get(A, col, col, n)));
                set(A, col, row, n, powf(c,2)* get(A, i, i, n) + powf(s,2) * get(A, col,  col, n) - 2 * s * c * get(A, row, col, n));
            }
        }
        off_of_Atag = calc_off_square(A,n);
        /*update V */
        for(i = 0;i<n;i++){
            set(V, i, row, n, c * get(V, i, row, n) - s * get(V, i , col, n));
            set(V, i, row, n, c * get(V, i, col, n) + s * get(V, i , row, n));
        }
        
    }
    
}

void form_diagonal_mat(float *diagonal_mat, float *weighted_adj_mat , int n){
    int i,j;
    float d = 0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            d += weighted_adj_mat[0];
            weighted_adj_mat++;
        }
        diagonal_mat[i] = d;
    }
}


void form_weighted_adj_mat(float *mat, float *data_points, int dim, int n){
    float distance(float * , float *, int, int, int),get(float *,int, int, int),w;
    void set(float *, int, int, int, float);
    int i,j;
    
    for(i = 0;i<n;i++){
        for(j=0;j<n;j++){
            w = distance(data_points, data_points, dim, i, j)/2;
            w = exp(-w);
            set(mat, i, j, n , w);
        }
    }
}


/** given D=diag(diagonal_mat) and W=weights_mat
 * calculates I-DWD into normalized_laplacian */
void calc_normalized_laplacian(float *normalized_laplacian, float *diagonal_mat, float *weights_mat, int dim){
    void set(float *, int, int, int, float);
    float get(float*, int, int, int), result = 0;
    int i,j;

    for(i = 0;i<dim;i++){
        for(j=0;j<dim;j++){
            if(i==j) {
                result = 1 - (1/diagonal_mat[i]) * get(weights_mat, i, i, dim);
            }else{
                result = (-1) * (1/sqrt(diagonal_mat[i])) * get(weights_mat, i, j, dim) * (1/sqrt(diagonal_mat[j]));
            }
            set(normalized_laplacian, i, j, dim , result);
        }
    }
}


/** given a square matrix (mat), sets row and col to be
 * the indexes of the off-diagonal largest absolute value */
void indexes_of_max_off_diag(float *mat, int *row, int *col, int dim){
    float max=0,val;
    int i,j;

    /*set indexes to 0 by defualt, for the case where all off-diagonal values are 0 */
    row=0;
    col=0;

    for(i = 0;i<dim;i++){
        for(j=0;j<dim;j++){
            val = get(mat,i,j,dim);
            if( i!=j && abs(val) >= max ){
                max = val ;
                row = i;
                col = j;
            }
        }
    }
}


/** given a square matrix (mat) ,
 * the indexes of its off-diagonal largest absolute value (i,j)
 * and it's dimension (dim)
 * returns t
 */
float get_t(float *mat, int i, int j, int dim){
    float theta, t;

    theta = (get(mat,j,j,dim)-get(mat,i,i,dim))/(2*get(mat,i,j,dim));
    t=1/(abs(theta)+ sqrt(powf(theta,2)+1));
    if ( theta < 0 ){
        t = (-1) * t;
    }
    return t;
}


/** given t, returns c
 */
float get_c(float t){
    return  1/ (sqrt(powf(t,2)+1));
}
/* given matrix mat calculates the off(mat)^2 = notm_f(mat) - sum_(i=1)^n a_ii*/
float calc_off_square(float *mat,int n){
    float result = 0,get(float *, int ,int, int);
    int i,j;

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i != j){
                result += powf(get(mat, i, j, n),2);
            }
        }
    }
    return result;

}


/**
 * given an array [arr], and an array- indexes of length [dim] such that indexes[i]=i for any i,
 * changes the order of the indexes such that the new order is the sorted order of the array values.
 * (using QuickSort algorithm)
 */
void quickSort_indexes(float *arr, int *indexes, int dim){
    void quickSort_rec(float*, int*, int, int);
    quickSort_rec(arr, indexes, 0, dim-1);
}


void quickSort_rec(float *arr, int *indexes, int low, int high) {
    void quickSort_rec(float*, int *, int, int);
    int partition(float *, int *, int, int);
    int mid;
    if (low < high) {
        mid = partition(arr, indexes, low, high);
        quickSort_rec(arr, indexes, low, mid - 1);
        quickSort_rec(arr, indexes, mid + 1, high);
    }
}


int partition(float *arr, int *indexes, int low, int high){
    float pivot = arr[indexes[high]];
    int i = low-1, j, temp;

    for(j=low; j <= high-1; j++){
        if (arr[indexes[j]] < pivot){
            i++;
            temp = indexes[i];
            indexes[i] = indexes[j];
            indexes[j] = temp;
        }
    }
    temp = indexes[i+1];
    indexes[i+1] = indexes[high];
    indexes[high] = temp;

    return i+1;
}



/** fill given matrix with the data from input file */
void read_data(FILE* fp, float *data_points, char *line, int n, int dim){
    void set(float *, int, int, int, float);
    char *p;
    int size, i, j;

    size = 50 * dim;
    for(i = 0; i < n; i++){

        /* p is a pointer to a beging of a line in the file */
        p = fgets(line, size, fp);

        for(j = 0; j < dim; j++){

            /* extract float form the line */
            set(data_points, i, j, dim + 1, strtod(p, &p));

            p += 1; /* skip comma */
        }
    }
}

/*
* calculate distance between two vectors sqrt(sum_(i=1)^n (v1_i-v2_i)^2) (norm2)
*/
float distance(float *v1, float *v2, int dim,int row_v1, int row_v2){
    int i;
    float result = 0;
    float x;
    float get(float *, int, int, int);

    for(i = 0;i < dim;i++){
        x = (get(v1, row_v1, i, dim + 1)-get(v2, row_v2, i, dim));
        x *= x;
        result += x;
    }
    sqrt(result);
    return result;
}


float get(float* arr, int i, int j, int dim){
    int index;

    index = (i*dim + j);
    return arr[index];
}


void set(float* arr, int i, int j, int dim, float item){
    int index;

    index = (i*dim + j);
    arr[index] = item;
}


int main( int argc, char* argv[]) {
    void read_data(FILE*, float *, char *, int, int ), form_weighted_adj_mat(float *, float *, int , int);
    int max_iter, dim, k, n;
    long bOfFile;
    float *diagonal_mat;
    float *data_points;
    float *centroids;
    float *weighted_adj_mat;
    char *line;

    /* reading arguments */
    k = strtol(argv[1], NULL, 10);
    if(argc == 3){
        max_iter = strtol(argv[2], NULL, 10);

        if (max_iter<=0){
            printf("INPUT ERROR:\nmaximum iterations is invalid");
            return 1;
        }
    }

    else{
        max_iter = 200;
    }
    bOfFile = ftell(stdin);/*save the address of the beginning of the file */
    n = num_of_lines(stdin);

    if(k<=0){
        printf("INPUT ERROR:\nk is invalid");
        return 1;
    }
    if(n <= k){
        printf("INPUT ERROR:\nthere are less then k=%d data points",k);
        return 1;
    }

    fseek(stdin, bOfFile, SEEK_SET);/*set the file position back to the beginning */
    dim = num_of_columns(stdin);
    fseek(stdin, bOfFile, SEEK_SET);/*set the file position back to the beginning */

    line = malloc(sizeof(char) * (30*dim));
    /* build matrix that contins all the points */
    data_points = malloc(sizeof(float) * dim * n);
    weighted_adj_mat = malloc(sizeof(float) * n * n);
    /*vactor that represent the D matrix*/
    diagonal_mat = malloc(sizeof(float) * n);
    /* build matrix that contins all the centroids */
    centroids = malloc(sizeof(float) * (dim*k));
    read_data(stdin, data_points, line, n, dim);
    form_weighted_adj_mat(weighted_adj_mat,data_points, dim, n);
    form_diagonal_mat(diagonal_mat, weighted_adj_mat, n);
    print_centroids(centroids, k, dim);


    /* free the memory used */
    free(line);
    free(data_points);
    free(weighted_adj_mat);
    free(centroids);

    return 0;
}
