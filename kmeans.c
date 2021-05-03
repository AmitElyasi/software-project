#include <stdlib.h>
#include <stdio.h>

int main(int argc, char const *argv[])
{
    int i;
    int j;
    int n;
    int k = strtol(argv[1], NULL, 10);
    int max_iter;
    long int bOfFile;
    int num_of_lines(FILE *);
    void fillVectors(FILE *, float **, int , int);

    if(argc == 3){
        max_iter = strtol(argv[2], NULL, 10);
    }else max_iter = 200;
    bOfFile = ftell(stdin);
    n = num_of_lines(stdin);
    fseek(stdin, bOfFile, SEEK_SET);
    float **vectors = (float **) malloc(sizeof(float *) * n); 
    for(i = 0; i < n; i++){
        vectors[i] = (float *) malloc(sizeof(float) * k);
    }
    fillVectors(stdin, vectors, k, n);
    return 0;
}

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
void fillVectors(FILE *fp, float **vectors, int k, int n){
    char *line = malloc(sizeof(char) * (7 * k));
    int ch;
    int size = 7 * k;
    int i , j;
    char *p;
    
    for(i = 0; i < n; i++){
        p = fgets(line, size, stdin);
        for(j = 0; j < k; j++){
            vectors[i][j] = strtod(p, &p);
            p += 1;
        }
    }

    
}
    