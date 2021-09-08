/**
 * Normalized Spectral Clustering implementation
 * Software Project
 *
 * Nizan Shemi
 * 206962912
 *
 * Amit Elyasi
 * 316291434
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


/*start of Normalized Spectral Clustering implementation*/

/*find the eigangap heuristic and return the number of cluster*/

int calc_eiganvalue_gap(double  *mat, int *sorted_eiganvalues_indexes, int n){
    double *deltas;
    double get(double  *, int, int, int);
    int i,index,next_index,half_n,max = 0,result = 0;
    
    deltas = malloc(sizeof(double ) * n);

    for(i=0;i<n-1;i++){
        index = sorted_eiganvalues_indexes[i];
        next_index = sorted_eiganvalues_indexes[i+1];
        deltas[i] = fabs(get(mat, index, index, n) - get(mat, next_index , next_index, n));
    }
    half_n = (int) (n/2);
    for(i=0;i<half_n;i++){
        if(max < deltas[i]){
            max = deltas[i];
            result = i;
        }
    }
    return result;
}

/** given A real symmetric matrix updates A and a matrix V that A is diagonal and all of A eigenvalues
 * on the diagonal and V columns is the correspondence eigenvector
 */
void jacobi_algorithm_for_eigenvalues(double  *A, double  *V, int n){
    double get(double *, int, int, int),get_c(double ),get_t(double  *,int, int, int),calc_off_square(double  *, int);
    double off_of_A,off_of_Atag,epsilon,s,c,t,val_row,val_col,a_ii,a_jj,a_ij;
    void indexes_of_max_off_diag(double  *, int *, int *, int),set(double  *, int, int, int ,double );
    int row,col,i,counter = 0;
    /*if A diagonl metrix we skip the while loop with this starting settings*/
    off_of_A = 0;
    off_of_Atag = calc_off_square(A, n);
    epsilon = pow(10, -15);


    while(fabs(off_of_A - off_of_Atag) > epsilon && counter < 100 ){
        counter++;
        off_of_A = off_of_Atag;
        indexes_of_max_off_diag(A, &row, &col, n);
        t = get_t(A, row, col, n);
        c = get_c(t);
        s = t * c;

        /*update A */
        for(i=0;i<n;i++){
            if(i != row && i != col){
                val_row = c * get(A, i, row , n) - s * get(A, i, col, n);
                val_col = c * get(A, i, col , n) + s * get(A, i, row, n);
                set(A, i, row, n, val_row);
                set(A, row, i, n, val_row);
                set(A, i, col, n, val_col);
                set(A, col, i, n, val_col);
            }
        }
        a_ii = get(A, row, row, n);
        a_jj = get(A, col, col, n);
        a_ij = get(A, row, col, n);
        set(A, row, row, n, (c * c) * a_ii + (s * s) * a_jj - 2 * s * c * a_ij);
        set(A, col, col, n, (s * s) * a_ii + (c * c) * a_jj + 2 * s * c * a_ij);
        set(A, row, col, n, ((c * c) - (s * s)) * a_ij + s * c * (a_ii - a_jj));
        set(A, col, row, n, ((c * c) - (s * s)) * a_ij + s * c * (a_ii - a_jj));
            
        off_of_Atag = calc_off_square(A,n);
        /*update V */
        for(i = 0;i<n;i++){
            val_row = get(V, i, row, n);
            val_col = get(V, i, col, n);
            
            set(V, i, row, n, c * val_row  - s * val_col);
            set(V, i, col, n, c * val_col + s * val_row);
        }

    }
}

/*initial V to the uint matrix*/
void form_V(double  *mat, int n){
    void set(double  *, int, int, int, double );
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){
                set(mat , i, j, n, 1.0);
            }else{
                set(mat , i, j, n, 0.0);
            }
        }
    }
}

void form_diagonal_mat(double  *diagonal_mat, double  *weighted_adj_mat , int n){
    int i,j;
    double  d = 0;
    for(i=0;i<n;i++){
        d = 0;
        for(j=0;j<n;j++){
            d += weighted_adj_mat[0];
            weighted_adj_mat++;
        }
        diagonal_mat[i] = d;
    }
}


void form_weighted_adj_mat(double  *mat, double  *data_points, int dim, int n){
    double  distance(double  * , double  *, int, int, int),get(double  *,int, int, int),w;
    void set(double  *, int, int, int, double );
    int i,j;

    for(i = 0;i<n;i++){
        for(j=0;j<n;j++){
            if(i == j){
                set(mat, i, j, n , 0);    
            }else{
                w = distance(data_points, data_points, dim, i, j)/2;
                w = exp(-w);
                set(mat, i, j, n , w);
            }
        }
    }
}


/** given D=diag(diagonal_mat) and W=weights_mat
 * calculates I-DWD into normalized_laplacian */
void calc_normalized_laplacian(double  *normalized_laplacian, double  *diagonal_mat, double  *weights_mat, int dim){
    void set(double  *, int, int, int, double );
    double  get(double *, int, int, int), result = 0;
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
void indexes_of_max_off_diag(double  *mat, int *row, int *col, int dim){
    double  max=0,val=0,get(double  *, int, int, int);
    int i,j;

    for(i = 0;i<dim;i++){
        for(j=i+1;j<dim;j++){
            val = get(mat,i,j,dim);
            if( i!=j && fabs(val) > max ){
                max = fabs(val);
                *row = i;
                *col = j;
            }
        }
    }
}


/** given a square matrix (mat) ,
 * the indexes of its off-diagonal largest absolute value (i,j)
 * and it's dimension (dim)
 * returns t
 */
double  get_t(double  *mat, int i, int j, int dim){
    double  theta,t,get(double  *, int, int, int);

    theta = (get(mat,j,j,dim)-get(mat,i,i,dim))/(2*get(mat,i,j,dim));
    t=1/(abs(theta)+ sqrt(pow(theta,2)+1));
    if ( theta < 0 ){
        t = (-1) * t;
    }
    return t;
}


/** given t, returns c
 */
double  get_c(double  t){
    return  1/ (sqrt(pow(t,2)+1));
}
/* given matrix mat calculates the off(mat)^2 = notm_f(mat) - sum_(i=1)^n a_ii*/
double calc_off_square(double  *mat,int n){
    double  result = 0;
    double  get(double  *, int ,int, int);
    int i,j;

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i != j){
                result += pow(get(mat, i, j, n),2);
            }
        }
    }
    return result;
}


/**
 * given diagonal matrix [n*n] and an empty array [n], fills the array with the diagonal values
 */
void diag_to_array(double  *diag, double  *arr, int n){
    double  get(double  *, int, int, int);
    int i;
    for (i=0;i<n;i++){
        arr[i]=get(diag, i, i, n);
    }
}


/** given array of length dim, and its dimention, fills the array with integers
 * such that arr[i]=i
 */
void fill_ordered_ints(int *arr, int dim){
    int i;
    for (i=0 ; i<dim ; i++){
        arr[i]=i;
    }
}


/**
 * given an array- arr, and another array- indexes, both of length dim,
 * fills indexes with integers [0,dim-1] such that their order dictates a sorted order on arr values.
 * (using QuickSort algorithm)
 */
void quickSort_indexes(double  *arr, int *indexes, int dim){
    void fill_ordered_ints(int *, int );
    void quickSort_rec(double *, int*, int, int);
    

    fill_ordered_ints(indexes, dim);
    quickSort_rec(arr, indexes, 0, dim-1);
}

/** recursive QuickSort
 */
void quickSort_rec(double  *arr, int *indexes, int low, int high) {
    void quickSort_rec(double *, int *, int, int);
    int partition(double  *, int *, int, int);
    int mid = 0;

    if (low < high) {
        mid = partition(arr, indexes, low, high);
        quickSort_rec(arr, indexes, low, mid - 1);
        quickSort_rec(arr, indexes, mid + 1, high);
    }
}

/** Implementation of the partition algorithm that is used for QS
 */
int partition(double  *arr, int *indexes, int low, int high){
    double  pivot = arr[indexes[high]];
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
void read_data(FILE* fp, double  *data_points, char *line, int n, int dim){
    void set(double  *, int, int, int, double );
    char *p;
    int size, i, j;

    size = 50 * dim;
    for(i = 0; i < n; i++){

        /* p is a pointer to a beging of a line in the file */
        p = fgets(line, size, fp);

        for(j = 0; j < dim; j++){

            /* extract double  form the line */
            set(data_points, i, j, dim, strtod(p, &p));

            p += 1; /* skip comma */
        }
    }
}

/*
 * counts the lines in input file
 */
int num_of_lines(FILE *fp){
    int ch;
    int lines = 1;

    while(!feof(fp))
    {
        ch = fgetc(fp);
        if(ch == '\n')
        {
            lines++;
        }
    }
    return lines;
}

/* count the columns in input file */
int num_of_columns(FILE *fp){
    int ch;
    int columns = 1;

    while(!feof(fp))
    {
        ch = fgetc(fp);
        if(ch == ','){
            columns++;
        }
        if(ch == '\n'){
            break;
        }
    }
    return columns;
}

/*
 * prints the centroids in the template requierd
 */
void print_matrix(double * mat, int row, int col){
    int i,j,dag = 0;
    double  get(double  *, int, int, int);
    if(row == -1){
        dag = 1;
        row = col;
    }

    for(i = 0;i<row;i++){
        for(j = 0;j<col-1;j++){
            if(dag == 1){
                if(i == j){
                    printf("%0.4f,", mat[i]);
                }else{
                    printf("%0.4f,", 0.0);
                }
            }else{
                printf("%0.4f,", get(mat, i, j, col));
            }
        }
        if(dag == 1){
            if(i == j){
                printf("%0.4f\n", mat[i]);
            }else{
                    printf("%0.4f\n", 0.0);
            }
        }else{
            printf("%0.4f\n", get(mat, i, col-1, col));
        }
    }
}


/** given the eigenvectors matrix [n*n] (as columns), an array [n] that dictates an order on the eigenvalues,
 * n (the mat dimention), k, and a result matrix [n*k],
 * fills the result matrix with the normalized first k eigenvectors as columns
 */
void normalized_k_eigenvectors(double  *eigenvecs, int *indexes, int n, int k, double  *result){
    double get(double  *, int, int, int);
    void set(double  *, int, int, int, double );
    int i, j, t;

    double sum = 0;

    for(i=0;i<k;i++){
        j=indexes[i];
        for (t=0; t < n; t++){
            sum += get(eigenvecs, t, j, n)*get(eigenvecs, t, j, n);
        }
        sum = sqrt(sum);
        for (t=0; t < n; t++){
            set(result, t, i, (k+1), get(eigenvecs, t, j, n)/sum);
        }
    }
}


/** calculate distance between two vectors sqrt(sum_(i=1)^n (v1_i-v2_i)^2) (L2 norm)
 */
double distance(double  *v1, double  *v2, int dim,int row_v1, int row_v2){
    int i;
    double  result = 0;
    double  x;
    double  get(double  *, int, int, int);

    for(i = 0;i < dim;i++){
        x = (get(v1, row_v1, i, dim)-get(v2, row_v2, i, dim));
        x *= x;
        result += x;
    }
    result = sqrt(result);
    return result;
}

/** calculate distance between two vectors sqrt(sum_(i=1)^n (v1_i-v2_i)^2) (L2 norm)
 */
double kmeans_distance(double  *v1, double  *v2, int dim,int row_v1, int row_v2){
    int i;
    double  result = 0;
    double  x;
    double  get(double  *, int, int, int);

    for(i = 0;i < dim;i++){
        x = (get(v1, row_v1, i, dim + 1)-get(v2, row_v2, i, dim));
        x *= x;
        result += x;
    }
    sqrt(result);
    return result;
}


/**
 * returns arr[i][j]
 */
double get(double * arr, int i, int j, int dim){
    int index;
    double x;

    index = (i*dim + j);
    x = arr[index];
    return arr[index];
}


/**
 * Sets arr[j][i]=item
 */
void set(double * arr, int i, int j, int dim, double  item){
    int index;

    index = (i*dim + j);
    arr[index] = item;
}

/* kmeans code*/

int kmeans(int k, double  *data_points, double  *centroids, double  *utl ,int max_iter, int dim, int n){
    void assign(double  *, double  *, int, int, int);
    short re_estimate(double  *,double  *, double  *, int , int, int);
    short convergece = 1;

    int i;

    for (i=0; i<max_iter; i++){
        assign(data_points, centroids, dim, n, k);
        convergence = re_estimate(data_points, centroids, utl, dim, n, k);
        if (convergence == 1) {
            return 0;
        }
    }
    return 0;
}


/*
 * assigns data points to their closest cluster (measure distance from the centroid)
 * updates the number of cluster for each data point
 */
void assign(double * data_points, double * clusters, int dim, int n, int k){
    int int_max = 2147483647;
    int cluster = 0;
    int v,c;
    void set(double  *, int, int, int, double );
    double  min_dis, dis,kmeans_distance(double  *,double  *, int , int , int);

    min_dis = int_max;
    for(v = 0; v < n; v++){
        for(c = 0;c < k; c++){
            dis = kmeans_distance(data_points, clusters, dim, v, c);
            if( dis <= min_dis){
                min_dis = dis;
                cluster = c;
            }
        }
        set(data_points, v, dim, dim + 1, cluster);
        min_dis = int_max;
    }
}


/*
 * re-estimates a centroid for each cluster:
 * for each cluster calculate the average of the points assign to it,
 * updates centroids to be the average vector,
 * returns 1 if the old centroids are equal to the new ones.
 */
short re_estimate(double * data_points, double * clusters,double  *utl, int dim, int n, int k){
    short isEqual = 1;
    int i, j;
    double  x,get(double  *, int, int, int);
    void zero_mat(double  *, int , int), set(double  *, int , int, int, double ),vec_sum(double  *, double  *, int, int, int);

    zero_mat(utl, dim + 1, k);

    /* sum all vectors for each cluster */
    for (i = 0; i < n; i++) {
        j = get(data_points, i, dim, dim+1);
        vec_sum(utl, data_points, dim, j, i);
        x = get(utl, j, dim, dim + 1) + 1;
        set(utl, j, dim, dim + 1, x);
    }

    /* Divides each sum by the number of vectors to get average */
    for (i = 0; i < k; i++) {
        for (j = 0; j < dim; j++) {
            x = get(utl, i, j, dim+1);
            set(utl, i, j, dim + 1, (x / get(utl, i, dim, dim+1)));
        }
    }

    /* Compare the old centroids to the new ones */
    for (i = 0; i < k; i++) {
        for (j = 0; j < dim; j++) {
            if (!(((get(clusters, i, j, dim) - get(utl, i, j, dim+1)) < 0.000001) && ((get(clusters, i, j, dim) - get(utl, i, j, dim+1)) > -0.000001))) {
                isEqual = 0;
                break;
            }
        }
        if (!isEqual){
            break;
        }
    }

    /* if there is no change in centroid- we have reach convergence */
    if (isEqual == 1) {
        return 1;
    }

    /* copy the new centroids to the old ones place */
    for (i = 0; i < k; i++) {
        for (j = 0; j < dim; j++) {
            x = get(utl, i, j, dim+1);
            set(clusters, i, j, dim, x);
        }
    }
    return isEqual;
}

/*
 * adds vec2 to vec1 coordinate wise
 */
void vec_sum(double * vec1, double * vec2, int dim, int row_vec1, int row_vec2){
    int i;
    void set(double  *, int , int, int, double );
    double sum,get(double  *, int, int, int);
    
    for(i = 0;i < dim;i++){
        sum = get(vec1, row_vec1, i, dim+1) + get(vec2, row_vec2, i, dim+1);
        set(vec1, row_vec1, i, dim + 1, sum);
    }
}


/*
 * zeros a given matrix from row start to row end
 */
void zero_mat(double * clusters , int dim, int n){
    int i,j;
    void set(double *, int , int, int, double );

    for(i = 0; i < n; i++){
        for(j=0; j < dim; j++){
            set(clusters, i, j, dim, 0);
        }
    }
}

/*return target matrix for goal == "wam" */
double *wam(double  *data_points, int n, int dim){
    void form_weighted_adj_mat(double  *, double  *, int, int);
    double  *target_matrix;

    target_matrix = malloc(sizeof(double ) * n * n);
    form_weighted_adj_mat(target_matrix, data_points, dim, n);

    return target_matrix;
}

/*return an array that contains the diagnoal of the target mateix for goal == "ddg" */
double *ddg(double  *data_points,int n, int dim){
    void form_diagonal_mat(double  *, double  *, int); 
    double  *target_diagnoal, *weighted_adj_mat, *wam(double  *, int, int);

    target_diagnoal = malloc(sizeof(double ) * n);
    weighted_adj_mat = wam(data_points, n, dim);
    form_diagonal_mat(target_diagnoal, weighted_adj_mat, n);

    free(weighted_adj_mat);

    return target_diagnoal;
}

/*return target matrix for goal == "lnorm" */
double *lnorm(double  *data_points, int n, int dim){
    void form_diagonal_mat(double  *, double  *, int),form_weighted_adj_mat(double  *, double  *, int, int);
    void calc_normalized_laplacian(double  *, double  *, double  * , int);
    double  *target_matrix,*weighted_adj_mat,*diagonal_mat;

    target_matrix = malloc(sizeof(double ) * n * n);
    weighted_adj_mat = malloc(sizeof(double ) * n * n);
    diagonal_mat = malloc(sizeof(double ) * n);

    form_weighted_adj_mat(weighted_adj_mat, data_points, dim, n);
    form_diagonal_mat(diagonal_mat, weighted_adj_mat, n);
    calc_normalized_laplacian(target_matrix, diagonal_mat, weighted_adj_mat, n);

    free(weighted_adj_mat);
    free(diagonal_mat);

    return target_matrix;
}

/*return target matrix in shape (n+1,n) the first row is the eigonvalue the other n row is the eigonvector matrix when the goal == "lnorm" */
double *jacobi(double  *data_points, int n){
    double  *traget_matrix,*V;
    void jacobi_algorithm_for_eigenvalues(double  *, double  *, int), form_V(double  *, int);
    int i,j;

    traget_matrix = malloc(sizeof(double ) * n * (n+1));
    V = malloc(sizeof(double ) * n * n);
    form_V(V, n);
    jacobi_algorithm_for_eigenvalues(data_points, V, n);

    for(i=0;i<=n;i++){
        for(j=0;j<n;j++){
            if(i == 0){
                set(traget_matrix, i, j, n, get(data_points, j, j, n));
            }else{
                set(traget_matrix, i, j, n, get(V, j, (i-1), n));
            }
        }
    }

    free(V);

    return traget_matrix;
}

/*return target matrix T in shape (n,k) when goal == "spk"*/
double *spk(double  *data_points, int n , int dim, int *k){
    double *traget_matrix,*lnorm(double  *, int, int),*V,*normalized_laplacian,*eigonvalues;
    void jacobi_algorithm_for_eigenvalues(double  *, double  *, int), form_V(double  *, int);
    void diag_to_array(double  *, double  *, int),quickSort_indexes(double  *, int *, int)
    ,normalized_k_eigenvectors(double  *, int *, int, int, double  *);
    int *indexes;
    
    V = malloc(sizeof(double ) * n * n);
    indexes = malloc(sizeof(int) * n);
    eigonvalues= malloc(sizeof(double ) * n);
    normalized_laplacian = lnorm(data_points, n, dim);
    print_matrix(normalized_laplacian, n, n);

    form_V(V, n);
    jacobi_algorithm_for_eigenvalues(normalized_laplacian, V, n);
    print_matrix(normalized_laplacian, n, n);
    print_matrix(V, n, n);
    diag_to_array(normalized_laplacian, eigonvalues, n);
    quickSort_indexes(eigonvalues, indexes, n);

    if(*k == 0){
        *k = calc_eiganvalue_gap(normalized_laplacian, indexes, n);
    }

    traget_matrix = malloc(sizeof(double ) * n * (*k+1));
    normalized_k_eigenvectors(V, indexes, n, *k, traget_matrix);

    free(normalized_laplacian);
    free(indexes);
    free(eigonvalues);
    free(V);

    return traget_matrix;
}

int main( int argc, char* argv[]) {

    void print_matrix(double *, int, int), read_data(FILE*, double  *, char *, int, int);
    int kmeans(int , double  *, double  *, double  *, int, int ,int); 
    int max_iter, dim, k, n,rows,cols,i,j;
    FILE *f;
    long bOfFile;
    double *target_matrix,*data_points,*centroids,*util,*wam(double  *, int , int),*lnorm(double  *, int , int),*ddg(double  *, int , int)
          ,*spk(double  *, int , int, int *),*jacobi(double  *, int);
    char *line,*goal;

    max_iter = 300;
    if(argc != 4){
        printf("invalid input");
        return 1;
    }
    /* reading data form arguments */
    k = strtol(argv[1], NULL, 10);
    goal = argv[2];
    f = fopen(argv[3],"r");
    bOfFile = ftell(f);/*save the address of the beginning of the file */
    n = num_of_lines(f);
    fseek(f, bOfFile, SEEK_SET);/*set the file position back to the beginning */
    dim = num_of_columns(f);
    fseek(f, bOfFile, SEEK_SET);/*set the file position back to the beginning */
    line = malloc(sizeof(char) * (30*dim));

    /* build matrix that contins all the points */
    data_points = malloc(sizeof(double ) * dim * n);
    read_data(f, data_points, line, n, dim);
    fclose(f);

    /*calculate the goal matrix*/
    if(!strcmp(goal, "wam")){
        target_matrix = wam(data_points, n, dim);
        rows = n;
        cols = n;
    }
    else if(!strcmp(goal, "ddg")){
        target_matrix = ddg(data_points, n, dim);
        rows = -1;
        cols = n;
    }
    else if(!strcmp(goal, "lnorm")){
        target_matrix = lnorm(data_points, n, dim);
        rows = n;
        cols = n;
    }
    else if(!strcmp(goal, "jacobi")){
        target_matrix = jacobi(data_points, n);
        rows = (n+1);
        cols = n;
    }else{
        target_matrix = spk(data_points, n, dim, &k);
        rows = k;
        cols = k;
    }
    if(strcmp(goal, "spk")){
        print_matrix(target_matrix, rows, cols);

        free(line);
        free(data_points);
        free(target_matrix);
        return 0;
    }

    centroids = malloc(sizeof(double ) * k * k);
    util = malloc(sizeof(double ) * k * (k+1));
    for(i = 0;i<k;i++){
        for(j=0;j<k;j++){
            set(centroids, i, j, k, get(target_matrix, i, j, (k+1)));
        }
    }
    kmeans(k, target_matrix, centroids, util , max_iter, k, n);
    print_matrix(centroids, k, k);
    /* free the memory used */
    free(line);
    free(data_points);
    free(target_matrix);
    free(centroids);
    free(util);

    return 0;
}