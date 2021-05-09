#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*
 * implementations of kmeans algorithm
 */
float** kmeans(int k, float** data_points ,int max_iter, int dim, int n){

    void assign(float**, float** ,int, int, int);
    void re_estimate(float**, float** ,int, int, int);
    float** build_clusters(int, int, float**);

    float** clusters = build_clusters(k, dim, data_points);
    for (int i=0; i<max_iter; i++){
        assign(data_points, clusters, dim, n, k);
        re_estimate(data_points, clusters, dim, n, k);
    }
    return clusters;
}


/*
 * assigns data points to their closest cluster (measure distance from the centroid)
 * updates the number of cluster for each data point
 */
void assign(float** data_points, float** clusters, int dim, int n, int k){
    int cluster;
    float distance(float *, float *, int);

    float min_dis = INT_MAX;
    for(int v = 0; v < n; v++){
        for(int c = 0;c < k; c++){
            float dis = distance(data_points[v], clusters[c], dim);
            if( dis <= min_dis){
                min_dis = dis;
                cluster = c;
            }
        }
        data_points[v][dim] = cluster;
        min_dis = INT_MAX;
    }
}


/*
* calculate distance between two vectors (v1-v2)^2
*/
float distance(float *v1, float *v2, int dim){
    float result = 0;
    for(int i = 0;i < dim;i++){
        result += pow(v1[i]-v2[i], 2);
    }

    return result;
}


/*
 * adds vec2 to vec1 coordinate wise
 */
void vec_sum(float* vec1, float* vec2, int dim){
    for(int i = 0;i < dim;i++){
        vec1[i] += vec2[i];
    }
}


/*
 * re-estimates a centroid for each cluster:
 * for each cluster calculate the average of the points assign to it,
 * updates centroids to be the average vector
 */
void re_estimate(float** data_points, float** clusters, int dim, int n, int k){
    int j;
    void vec_sum(float*, float*, int);
    void zero_mat(float**, int, int);

    zero_mat(clusters, dim, k);

    for(int i=0; i<n ; i++){
        j = data_points[i][dim];
        vec_sum(clusters[j], data_points[i], dim);
        clusters[j][dim]++;
    }

    for (int i=0; i < k ; i++){
        for(int j=0; j<dim; j++){
            clusters[i][j] = clusters[i][j]/clusters[i][dim];
        }
    }
}


/*
 * zeros the cluster matrix
 */
void zero_mat(float** clusters , int dim , int k){
    for(int i=0; i<k; i++){
        for(int j=0; j<=dim; j++){
            clusters[i][j]=0;
        }
    }
}



float** read_data(FILE* fp, int n, int m){
    long int bOfFile;
    float ** centroids;
    int num_of_columns(FILE *);
    int num_of_lines(FILE *);
    void fillVectors(FILE *, float **, int , int);

    /* build matrix that contins all the points */
    float **vectors = (float **) malloc(sizeof(float *) * n);
    for(int i = 0; i < n; i++){
        vectors[i] = (float *) malloc(sizeof(float) * (m+1));
    }
    fillVectors(stdin, vectors, m, n);

    return vectors;
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


/* fill given matrix with the data from input file */
void fillVectors(FILE *fp, float **vectors, int m, int n){
    char *p;
    int size = 50 * m ;
    char *line = malloc(sizeof(char) * size);

    for(int i = 0; i < n; i++){
        p = fgets(line, size, stdin);/* p is a pointer to a beging of a line in the file */
        for(int j = 0; j < m; j++){
            vectors[i][j] = strtod(p, &p);/* extract float form the line */
            p += 1;/* skip comma */
        }
        vectors[i][m] = 0;
    }
}


float ** build_clusters(int k, int m, float **vectors){
    float **centroid;
    centroid = (float **) malloc(sizeof(float *) * k);
    for(int i = 0; i < k; i++){
        centroid[i] = (float *) malloc(sizeof(float) * (m+1));
        for(int j = 0;j < m;j++){
            centroid[i][j] = vectors[i][j];
        }
        centroid[i][m] = 0;
    }
    return centroid;
}


/*
 * prints the centroids in the template requierd
 */
void print_centroids(float** clusters, int k, int dim){
    for(int i = 0; i < k;i++){
        for(int j = 0;j < dim-1;j++){
            printf("%0.4f,", clusters[i][j]);
        }
        printf("%.4f\n", clusters[i][dim-1]);
    }
}


int main( int argc, char* argv[]) {
    float** kmeans(int, float** , int, int, int);
    float** read_data(FILE*, int, int );
    void print_centroids(float**, int, int);
    int max_iter;

    /* reading arguments */
    int k = strtol(argv[1], NULL, 10);
    if(argc == 3){
        max_iter = strtol(argv[2], NULL, 10);
    }
    else{
        max_iter = 200;
    }
    long bOfFile = ftell(stdin);/*save the address of the beginning of the file */
    int n = num_of_lines(stdin);
    fseek(stdin, bOfFile, SEEK_SET);/*set the file position back to the beginning */
    int m = num_of_columns(stdin);
    fseek(stdin, bOfFile, SEEK_SET);/*set the file position back to the beginning */


    float** data_points = read_data(stdin, n, m);
    float** centroids = kmeans(k, data_points, max_iter, m, n);
    print_centroids(centroids, k, m);

}
