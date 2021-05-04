#include <stdlib.h>
#include <stdio.h>

int main(int argc, char const *argv[])
{
    int i;
    int j;
    int n;
    long int bOfFile;
    int max_iter;
    int num_of_lines(FILE *);
    void fillVectors(FILE *, float **, int , int);

    int k = strtol(argv[1], NULL, 10);
    if(argc == 3){
        max_iter = strtol(argv[2], NULL, 10);
    }else max_iter = 200;
    bOfFile = ftell(stdin);/*save the addres of the begining of the file */
    n = num_of_lines(stdin);
    fseek(stdin, bOfFile, SEEK_SET);/*set the file postion back to the begining */
    /* build matrix that contins all the points */
    float **vectors = (float **) malloc(sizeof(float *) * n); 
    for(i = 0; i < n; i++){
        vectors[i] = (float *) malloc(sizeof(float) * k);
    }
    fillVectors(stdin, vectors, k, n);
    return 0;
    
}
/* count the line in input file */
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
/* fill given matrix with the data from input file */
void fillVectors(FILE *fp, float **vectors, int k, int n){
    int ch;
    int i , j;
    char *p;
    char *line = malloc(sizeof(char) * (8 * k));
    int size = 8 * k ;
    
    for(i = 0; i < n; i++){
        p = fgets(line, size, stdin);/* p is a pointer to a beging of a line in the file */
        for(j = 0; j < k; j++){
            vectors[i][j] = strtod(p, &p);/* extract float form the line */
            p += 1;/* skip comma */
        }
    }

    
}
    