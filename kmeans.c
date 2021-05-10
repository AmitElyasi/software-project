/*
 * Kmeans implementation in C
 * Software Project ex. 1
 *
 * Amit Elyasi
 * 316291434
 *
 * Nizan Shemi
 * 206962912
 *
 */

#include <stdio.h>
#include <stdlib.h>


/*
 * implementations of kmeans algorithm:
 * given n data point (in R^d), group the data into k clusters,
 * each data point is assigned to exacly one cluster (note that k<n)
 */
float** kmeans(int k, float** data_points ,int max_iter, int dim, int n){
    short convergece = 1;
    void assign(float**, float** ,int, int, int);
    short re_estimate(float**, float** ,int, int, int);
    float** build_clusters(int, int, float**);

    float** clusters = build_clusters(k, dim, data_points);
    for (int i=0; i<max_iter; i++){
        assign(data_points, clusters, dim, n, k);
        convergece = re_estimate(data_points, clusters, dim, n, k);
        if (convergece == 1) {
            return clusters;
        }
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
 * re-estimates a centroid for each cluster:
 * for each cluster calculate the average of the points assign to it,
 * updates centroids to be the average vector,
 * returns 1 if the old centroids are equal to the new ones.
 */
short re_estimate(float** data_points, float** clusters, int dim, int n, int k) {
    int j;
    void vec_sum(float *, float *, int);
    void zero_mat(float **, int, int, int);
    short isEqual = 1;

    zero_mat(clusters, dim + 1, k, 2 * k);

    // sum all vectors for each cluster
    for (int i = 0; i < n; i++) {
        j = data_points[i][dim];
        vec_sum(clusters[j + k], data_points[i], dim);
        clusters[j + k][dim]++;
    }

    // Divides each sum by the number of vectors to get average
    for (int i = k; i < 2 * k; i++) {
        for (int j = 0; j < dim; j++) {
            clusters[i][j] = clusters[i][j] / clusters[i][dim];
        }
    }

    // Compare the old centroids to the new ones
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < dim; j++) {
            if (!(((clusters[i][j] - clusters[i + k][j]) < 0.000001) && ((clusters[i][j] - clusters[i + k][j]) > -0.000001))) {
                isEqual = 0;
                break;
            }
        }
        if (!isEqual){
            break;
        }
    }

    // if there is no change in centroid- we have reach convergence
    if (isEqual == 1) {
        return 1;
    }

    // copy the new centroids to the old ones place
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < dim; j++) {
            clusters[i][j]=clusters[i+k][j];
        }
    }
}


/*
* calculate distance between two vectors (v1-v2)^2
*/
float distance(float *v1, float *v2, int dim){
    float result = 0;
    for(int i = 0;i < dim;i++){
        result += (v1[i]-v2[i])*(v1[i]-v2[i]);
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
 * zeros a given matrix from row start to row end
 */
void zero_mat(float** clusters , int dim , int start, int end){
    for(int i=start; i<end; i++){
        for(int j=0; j<=dim; j++){
            clusters[i][j]=0;
        }
    }
}


float** read_data(FILE* fp, int n, int dim){
    long int bOfFile;
    float ** centroids;
    int num_of_columns(FILE *);
    int num_of_lines(FILE *);
    void fillVectors(FILE *, float **, int , int);

    /* build matrix that contins all the points */
    float **vectors = (float **) malloc(sizeof(float *) * n);
    for(int i = 0; i < n; i++){
        vectors[i] = (float *) malloc(sizeof(float) * (dim+1));
    }
    fillVectors(stdin, vectors, dim, n);

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
void fillVectors(FILE *fp, float **vectors, int dim, int n){
    char *p;
    int size = 50 * dim ;
    char *line = malloc(sizeof(char) * size);

    for(int i = 0; i < n; i++){

        //p is a pointer to a beging of a line in the file
        p = fgets(line, size, stdin);

        for(int j = 0; j < dim; j++){

            //extract float form the line
            vectors[i][j] = strtod(p, &p);

            p += 1; // skip comma
        }
        vectors[i][dim] = 0;
    }
    free(line);
}


float ** build_clusters(int k, int dim, float **vectors){
    float **centroid;
    centroid = (float **) malloc(sizeof(float *) * 2 * k);

    for(int i = 0; i < k; i++){
        centroid[i] = (float *) malloc(sizeof(float) * (dim));
        for(int j = 0;j < dim;j++){
            centroid[i][j] = vectors[i][j];
        }
    }
    for(int i = k; i < 2 * k; i++){
        centroid[i] = (float *) malloc(sizeof(float) * (dim+1));
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

        if (max_iter<=0){
            printf("INPUT ERROR:\nmaximum iterations is invalid");
            return 1;
        }
    }
    else{
        max_iter = 200;
    }
    long bOfFile = ftell(stdin);/*save the address of the beginning of the file */
    int n = num_of_lines(stdin);

    if(k<=0){
        printf("INPUT ERROR:\nk is invalid");
        return 1;
    }
    if(n < k){
        printf("INPUT ERROR:\nthere are less then k=%d data points",k);
        return 1;
    }

    fseek(stdin, bOfFile, SEEK_SET);/*set the file position back to the beginning */
    int dim = num_of_columns(stdin);
    fseek(stdin, bOfFile, SEEK_SET);/*set the file position back to the beginning */


    float** data_points = read_data(stdin, n, dim);
    float** centroids = kmeans(k, data_points, max_iter, dim, n);
    print_centroids(centroids, k, dim);


    // free the memory used
    for(int i = 0;i < n;i++) {
        free(data_points[i]);
    }
    free(data_points);

    for(int i = 0;i < 2*k ;i++){
        free(centroids[i]);
    }
    free(centroids);
}