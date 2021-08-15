
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void form_diagonal_mat(float *diagonal_mat, float *weighted_adj_mat , int n, int dim){
    int i,j;
    float d = 0;
    for(i=0;i<n;i++){
        for(j=0;j<dim;j++){
            d += weighted_adj_mat[0];
            weighted_adj_mat++;
        }
        d = sqrt(d);
        diagonal_mat[i] = 1/d;
    }
} 
void form_weighted_adj_mat(float *mat, float *data_points, int dim, int n){
    float distance(float * , float *, int, int, int),get(float *,int, int, int),w;
    void set(float *, int, int, int, float);
    int i,j;
    
    for(i = 0;i<n;i++){
        for(j=0;j<dim;j++){
            w = distance(data_points, data_points, dim, i, j)/2;
            w = exp(-w);
            set(mat, i, j, dim , w);
        }
    }
}
/* fill given matrix with the data from input file */
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
    weighted_adj_mat = malloc(sizeof(float) * dim * n);
    /*vactor that represent the D matrix*/
    diagonal_mat = malloc(sizeof(float) * n);
    /* build matrix that contins all the centroids */
    centroids = malloc(sizeof(float) * (dim*k));
    read_data(stdin, data_points, line, n, dim);
    form_weighted_adj_mat(weighted_adj_mat,data_points, dim, n);
    form_diagonal_mat(diagonal_mat, n);
    print_centroids(centroids, k, dim);


    /* free the memory used */
    free(line);
    free(data_points);
    free(weighted_adj_mat);
    free(centroids);

    return 0;
}