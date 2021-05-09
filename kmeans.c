#include <stdlib.h>
#include <stdio.h>

int main(int argc, char const *argv[])
{
    int n,m,max_iter;
    long int bOfFile;
    float ** centroids;
    int num_of_columns(FILE *);
    int num_of_lines(FILE *);
    void fillVectors(FILE *, float **, int , int);

    int k = strtol(argv[1], NULL, 10);
    if(argc == 3){
        max_iter = strtol(argv[2], NULL, 10);
    }else max_iter = 200;
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
    centroids = build_centroids(k,m ,vectors);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            float x = vectors[i][j];
            printf("%f ", x);
        }
        printf("claster %f \n", vectors[i][m]);
    }
    return 0;
}

/* count the lines in input file */
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
        if(ch == '\n')

        {
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

float ** build_centroids(int k, int m, float **vectors){
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
    