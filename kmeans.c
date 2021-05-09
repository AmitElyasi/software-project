#include <stdio.h>
#include <stdlib.h>

/*
 * implementations of kmeans algorithm
 */
float** kmeans(int k, float** data_points ,int max_iter){
    void assign(float**, float** );
    void re_estimate(float**, float** );
    float** clusters = built_clusters(int, int, float**);
    
    //TODO
    }


/*
 * assigns data points to their closest cluster (measure distance from the centroid)
 * updates the number of cluster for each data point
 */
void assign(float** data_points, float** clusters){
    //TODO
}


/*
 * re-estimates a centroid for each cluster:
 * for each cluster calculate the average of the points assign to it,
 *  updates centroids to be the average vector
 */
void re_estimate(float** data_points, float** clusters){
    //TODO
}


float** read_data(FILE* fp){
    int n,m;
    long int bOfFile;
    float ** centroids;
    int num_of_columns(FILE *);
    int num_of_lines(FILE *);
    void fillVectors(FILE *, float **, int , int);

    bOfFile = ftell(stdin);/*save the address of the begining of the file */
    n = num_of_lines(stdin);
    fseek(stdin, bOfFile, SEEK_SET);/*set the file postion back to the begining */
    m = num_of_columns(stdin);
    fseek(stdin, bOfFile, SEEK_SET);/*set the file postion back to the begining */

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
    int size = 8 * m ;
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
void print_centroids(float** clusters){
    //TODO
}


int main() {
    float** kmeans(int, float ,int );
    float** read_data(FILE*);
    int max_iter;
    
    /*
     * reading arguments
     */
    int k = strtol(argv[1], NULL, 10);
    if(argc == 3){
        max_iter = strtol(argv[2], NULL, 10);
    }
    else{
        max_iter = 200; 
    } 
    float** data_points = read_data(stdin);
    
    float** centroids = kmeans(k, data_points, max_iter);
    print_centroids(centroids);            
    
}
