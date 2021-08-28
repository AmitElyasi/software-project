#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/**
 * given an array [arr], and an array- indexes of length [dim] such that indexes[i]=i for any i,
 * changes the order of the indexes such that the new order is the sorted order of the array values.
 * (using QuickSort algorithm)
 */
void quickSort(float *arr, int *indexes, int dim){
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

void my_print(float *arr, int *indexes, int dim){
    int i;
    printf("[");
    for (i=0; i<dim; i++){
        printf("%f ,",arr[indexes[i]]);
    }
    printf("]");
}

void fill_ordered_ints(int *arr, int dim){
    int i;
    for (i=0 ; i<dim ; i++){
        arr[i]=i;
    }
}

void fill_random_floats(float *arr, int dim){
    int i;
    for (i=0 ; i<dim ; i++){
        arr[i]=((float)rand()/(float)(RAND_MAX)) * 100;
    }
}


int main( int argc, char* argv[]){
    srand ( time(NULL) );
    int dim = 10 ;
    float *arr;
    int *indexes;

    arr = malloc(sizeof(float) * (dim));
    fill_random_floats(arr, dim);

    indexes = malloc(sizeof(int) * (dim));
    fill_ordered_ints(indexes, dim);

    my_print(arr, indexes, dim);
    quickSort(arr, indexes, dim);
    printf("\n");
    my_print(arr, indexes, dim);

}

