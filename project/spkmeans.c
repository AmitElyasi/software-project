#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


/*start of Normalized Spectral Clustering implementation*/

/*find the eigangap heuristic and return the number of cluster*/
int calc_eiganvalue_gap(float *mat, int *sorted_eiganvalues_indexes, int n){
    float *deltas;
    float get(float *, int, int, int);
    int i,index,half_n,max = 0,result = 0;
    
    deltas = malloc(sizeof(float) * n);
    for(i=0;i<n-1;i++){
        index = sorted_eiganvalues_indexes[i];
        deltas[i] = fabs(get(mat, index, index, n) - get(mat, index+1, index+1, n));
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
void jacobi_algorithm_for_eigenvalues(float *A, float *V, int n){
    float get(float*, int, int, int),get_c(float),get_t(float *,int, int, int),calc_off_square(float *, int);
    float off_of_A,off_of_Atag,epsilon,s,c,t,val;
    void indexes_of_max_off_diag(float *, int *, int *, int),set(float *, int, int, int ,float);
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
                val = c * get(A, i, row , n) - s * get(A, i, col, n);
                set(A, i, row, n, val);
                set(A, row, i, n, val);
                val = c * get(A, i, col , n) + s * get(A, i, row, n);
                set(A, i, col, n, val);
                set(A, col, i, n, val);
            }else{
                set(A, row, row, n, pow(c,2)* get(A, i, i, n) + pow(s,2) * get(A, col,  col, n) - 2 * s * c * get(A, row, col, n));
                set(A, col, col, n, pow(s,2)* get(A, i, i, n) + pow(c,2) * get(A, col, col, n) + 2 * s * c * get(A, row, col, n));
                set(A, row, col, n, (pow(c,2)*  - pow(s,2)) * get(A, row,  col, n) -  s * c * (get(A, row, row, n) - get(A, col, col, n)));
                set(A, col, row, n, pow(c,2)* get(A, i, i, n) + pow(s,2) * get(A, col,  col, n) - 2 * s * c * get(A, row, col, n));
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

/*initial V to the uint matrix*/
void form_V(float *mat, int n){
    void set(float *, int, int, int, float);
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){
                set(mat , i, j, n, 0.0);
            }else{
                set(mat , i, j, n, 1.0);
            }
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
    float max=0,val=0,get(float *, int, int, int);
    int i,j;

    for(i = 0;i<dim;i++){
        for(j=i;j<dim;j++){
            val = get(mat,i,j,dim);
            if( i!=j && abs(val) >= max ){
                if(max == val && *row >= i && *col >= j){
                max = val ;
                *row = i;
                *col = j;
            }
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
    float theta,t,get(float *, int, int, int);

    theta = (get(mat,j,j,dim)-get(mat,i,i,dim))/(2*get(mat,i,j,dim));
    t=1/(abs(theta)+ sqrt(pow(theta,2)+1));
    if ( theta < 0 ){
        t = (-1) * t;
    }
    return t;
}


/** given t, returns c
 */
float get_c(float t){
    return  1/ (sqrt(pow(t,2)+1));
}
/* given matrix mat calculates the off(mat)^2 = notm_f(mat) - sum_(i=1)^n a_ii*/
float calc_off_square(float *mat,int n){
    float result = 0;
    float get(float *, int ,int, int);
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
void diag_to_array(float *diag, float *arr, int n){
    float get(float *, int, int, int);
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
void quickSort_indexes(float *arr, int *indexes, int dim){
    void fill_ordered_ints(int *, int );
    void quickSort_rec(float*, int*, int, int);
    
    fill_ordered_ints(indexes, dim);
    quickSort_rec(arr, indexes, 0, dim-1);
}

/** recursive QuickSort
 */
void quickSort_rec(float *arr, int *indexes, int low, int high) {
    void quickSort_rec(float*, int *, int, int);
    int partition(float *, int *, int, int);
    int mid = 0;
    
    if (low < high) {
        mid = partition(arr, indexes, low, high);
        quickSort_rec(arr, indexes, low, mid - 1);
        quickSort_rec(arr, indexes, mid + 1, high);
    }
}

/** Implementation of the partition algorithm that is used for QS
 */
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
 * counts the lines in input file
 */
int num_of_lines(FILE *fp){
    int ch;
    int lines = 0;

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
void print_matrix(float* mat, int row, int col){
    int i,j,dag = 0;
    float get(float *, int, int, int);
    if(row == -1){
        dag = 1;
        row = col;
    }

    for(i = 0;i<row;i++){
        for(j = 0;j<col;j++){
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
void normalized_k_eigenvectors(float *eigenvecs, int *indexes, int n, int k, float *result){
    float get(float *, int, int, int);
    void set(float *, int, int, int, float);
    int i, j, t;
    float sum;

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


/**
 * returns arr[i][j]
 */
float get(float* arr, int i, int j, int dim){
    int index;

    index = (i*dim + j);
    return arr[index];
}


/**
 * Sets arr[j][i]=item
 */
void set(float* arr, int i, int j, int dim, float item){
    int index;

    index = (i*dim + j);
    arr[index] = item;
}

/* kmeans code*/
int kmeans(int k, float *data_points, float *centroids, float *utl ,int max_iter, int dim, int n){
    void assign(float *, float *, int, int, int);
    short re_estimate(float *,float *, float *, int , int, int);
    short convergece = 1;
    int i;

    for (i=0; i<max_iter; i++){
        assign(data_points, centroids, dim, n, k);
        convergece = re_estimate(data_points, centroids, utl, dim, n, k);
        if (convergece == 1) {
            return 0;
        }
    }
    return 0;
}


/*
 * assigns data points to their closest cluster (measure distance from the centroid)
 * updates the number of cluster for each data point
 */
void assign(float* data_points, float* clusters, int dim, int n, int k){
    int int_max = 2147483647;
    int cluster = 0;
    int v,c;
    void set(float *, int, int, int, float);
    float min_dis, dis,distance(float *,float *, int , int , int);

    min_dis = int_max;
    for(v = 0; v < n; v++){
        for(c = 0;c < k; c++){
            dis = distance(data_points, clusters, dim, v, c);
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
short re_estimate(float* data_points, float* clusters,float *utl, int dim, int n, int k){
    short isEqual = 1;
    int i, j;
    float x,get(float *, int, int, int);
    void zero_mat(float *, int , int), set(float *, int , int, int, float),vec_sum(float *, float *, int, int, int);

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
void vec_sum(float* vec1, float* vec2, int dim, int row_vec1, int row_vec2){
    int i;
    void set(float *, int , int, int, float);
    float sum,get(float *, int, int, int);
    
    for(i = 0;i < dim;i++){
        sum = get(vec1, row_vec1, i, dim+1) + get(vec2, row_vec2, i, dim+1);
        set(vec1, row_vec1, i, dim + 1, sum);
    }
}


/*
 * zeros a given matrix from row start to row end
 */
void zero_mat(float* clusters , int dim, int n){
    int i,j;
    void set(float *, int , int, int, float);
    
    for(i = 0; i < n; i++){
        for(j=0; j < dim; j++){
            set(clusters, i, j, dim, 0);
        }
    }
}

/*return target matrix for goal == "wam" */
float *wam(float *data_points, int n, int dim){
    void form_weighted_adj_mat(float *, float *, int, int);
    float *target_matrix;

    target_matrix = malloc(sizeof(float) * n * n);
    form_weighted_adj_mat(target_matrix, data_points, dim, n);

    return target_matrix;
}

/*return an array that contains the diagnoal of the target mateix for goal == "ddg" */
float *ddg(float *data_points,int n, int dim){
    void form_diagonal_mat(float *, float *, int); 
    float *target_diagnoal, *weighted_adj_mat, *wam(float *, int, int);


    target_diagnoal = malloc(sizeof(float) * n);
    weighted_adj_mat = wam(data_points, n, dim);
    form_diagonal_mat(target_diagnoal, weighted_adj_mat, n);
    
    free(weighted_adj_mat);

    return target_diagnoal;
}

/*return target matrix for goal == "lnorm" */
float *lnorm(float *data_points, int n, int dim){
    void form_diagonal_mat(float *, float *, int),form_weighted_adj_mat(float *, float *, int, int);
    void calc_normalized_laplacian(float *, float *, float * , int);
    float *target_matrix,*weighted_adj_mat,*diagonal_mat;

    target_matrix = malloc(sizeof(float) * n * n);
    weighted_adj_mat = malloc(sizeof(float) * n * n);
    diagonal_mat = malloc(sizeof(float) * n);

    form_weighted_adj_mat(weighted_adj_mat, data_points, dim, n);
    form_diagonal_mat(diagonal_mat, weighted_adj_mat, n);
    calc_normalized_laplacian(target_matrix, diagonal_mat, weighted_adj_mat, n);

    free(weighted_adj_mat);
    free(diagonal_mat);
    
    return target_matrix;
}

/*return target matrix in shape (n+1,n) the first row is the eigonvalue the other n row is the eigonvector matrix when the goal == "lnorm" */
float *jacobi(float *data_points, int n){
    float *traget_matrix,*V;
    void jacobi_algorithm_for_eigenvalues(float *, float *, int), form_V(float *, int);
    int i,j;

    traget_matrix = malloc(sizeof(float) * n * (n+1));
    V = malloc(sizeof(float) * n * n);
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
float *spk(float *data_points, int n , int dim, int *k){
    float *traget_matrix,*lnorm(float *, int, int),*V,*normalized_laplacian,*eigonvalues;
    void jacobi_algorithm_for_eigenvalues(float *, float *, int), form_V(float *, int);
    void diag_to_array(float *, float *, int),quickSort_indexes(float *, int *, int)
    ,normalized_k_eigenvectors(float *, int *, int, int, float *);
    int *indexes;

    
    V = malloc(sizeof(float) * n * n);
    indexes = malloc(sizeof(int) * n);
    eigonvalues= malloc(sizeof(float) * n);
    normalized_laplacian = lnorm(data_points, n, dim);

    form_V(V, n);
    jacobi_algorithm_for_eigenvalues(normalized_laplacian, V, n);
    diag_to_array(normalized_laplacian, eigonvalues, n);
    quickSort_indexes(eigonvalues, indexes, n);
    
    if(*k == 0){
        *k = calc_eiganvalue_gap(normalized_laplacian, indexes, n);
    }

    traget_matrix = malloc(sizeof(float) * n * (*k+1));
    normalized_k_eigenvectors(V, indexes, n, *k, traget_matrix);

    free(normalized_laplacian);
    free(indexes);
    free(eigonvalues);
    free(V);
    
    return traget_matrix;
}

int main( int argc, char* argv[]) {
    void print_matrix(float*, int, int), read_data(FILE*, float *, char *, int, int);
    int kmeans(int , float *, float *, float *, int, int ,int); 
    int max_iter, dim, k, n,rows,cols,i,j;
    FILE *f;
    long bOfFile;
    float *target_matrix,*data_points,*centroids,*util,*wam(float *, int , int),*lnorm(float *, int , int),*ddg(float *, int , int)
        ,*spk(float *, int , int, int *),*jacobi(float *, int);
    char *line,*goal;

    max_iter = 240;
    if(argc != 3){
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
    dim = num_of_columns(stdin);
    fseek(stdin, bOfFile, SEEK_SET);/*set the file position back to the beginning */
    line = malloc(sizeof(char) * (30*dim));
    
    /* build matrix that contins all the points */
    data_points = malloc(sizeof(float) * dim * n);
    read_data(f, data_points, line, n, dim);
    fclose(f);
    
    /*calculate the goal matrix*/
    if(strcmp(goal, "wam")){
        target_matrix = wam(data_points, n, dim);
        rows = n;
        cols = n;
    }
    else if(strcmp(goal, "ddg")){
        target_matrix = ddg(data_points, n, dim);
        rows = -1;
        cols = n;
    }
    else if(strcmp(goal, "lnorm")){
        target_matrix = lnorm(data_points, n, dim);
        rows = n;
        cols = n;
    }
    else if(strcmp(goal, "jacobi")){
        target_matrix = jacobi(data_points, n);
        rows = (n+1);
        cols = n;
    }else{
        target_matrix = spk(data_points, n, dim, &k);
        rows = k;
        cols = k;
    }
    if(!strcmp(goal, "spk")){
        print_matrix(target_matrix, rows, cols);
        
        free(line);
        free(data_points);
        free(target_matrix);
        return 0;
    }

    centroids = malloc(sizeof(float) * k * k);
    util = malloc(sizeof(float) * k * (k+1));
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
