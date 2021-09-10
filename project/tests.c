#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/**
 * given an array- arr, and another array- indexes, both of length dim,
 * fills indexes with integers [0,dim-1] such that their order dictates a sorted order on arr values.
 * (using stable QuickSort algorithm)
 */
void quickSort_indexes(double *arr, int *indexes, int dim){
    void fill_ordered_ints(int *, int );
    void quickSort_rec(double*, int*, int, int);
    void stable(double *, int *, int);

    fill_ordered_ints(indexes, dim);
    quickSort_rec(arr, indexes, 0, dim-1);
    stable(arr, indexes, dim);
}

/** recursive QuickSort
 */
void quickSort_rec(double *arr, int *indexes, int low, int high) {
    void quickSort_rec(double*, int *, int, int);
    int partition(double *, int *, int, int);
    int mid = 0;

    if (low < high) {
        mid = partition(arr, indexes, low, high);
        quickSort_rec(arr, indexes, low, mid - 1);
        quickSort_rec(arr, indexes, mid + 1, high);
    }
}

/** Implementation of the partition algorithm that is used for QS
 */
int partition(double *arr, int *indexes, int low, int high){
    double pivot = arr[indexes[high]];
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

/**
 * stables the sort
 */
void stable(double *arr, int *indexes, int dim){
    void arr_int_to_double(double*, int*, int);
    void sorted_double_to_int(double*, int*, int*, int);
    void fill_ordered_ints(int *, int );
    double *double_indexes;
    int low=0,high=1;
    int* ordered_ints;

    double_indexes = malloc(sizeof(double ) * (dim));
    ordered_ints = malloc(sizeof(int) * (dim));

    if(dim<=1){
        free(double_indexes);
        free(ordered_ints);
        return;
    }
    fill_ordered_ints(ordered_ints, dim);
    arr_int_to_double(double_indexes, indexes, dim);


    while(high < dim) {
        while (high + 1 < dim && arr[(int)double_indexes[high+1]] == arr[(int)double_indexes[low]]) {
            high += 1;
        }
        if (arr[(int)double_indexes[high]] != arr[(int)double_indexes[low]]){
            low += 1;
            high = low + 1;
            continue;
        }
        quickSort_rec(double_indexes, ordered_ints, low, high);
        low = high + 1;
        high = low + 1;
    }
    sorted_double_to_int(double_indexes,indexes,ordered_ints,dim);

    free(double_indexes);
    free(ordered_ints);
}


void my_print(double *arr, int *indexes, int dim){
    int i;
    printf("[");
    for (i=0; i<dim; i++){
        printf("%f ,",arr[indexes[i]]);
    }
    printf("]\n");
}

void index_print(int *indexes, int dim){
    int i;
    printf("[");
    for (i=0; i<dim; i++){
        printf("%i ,",indexes[i]);
    }
    printf("]\n");
}

void fill_ordered_ints(int *arr, int dim){
    int i;
    for (i=0 ; i<dim ; i++){
        arr[i]=i;
    }
}

void fill_random_doubles(double *arr, int dim){
    int i;
    for (i=0 ; i<dim ; i++){
        arr[i]=((double)rand()/(double)(RAND_MAX)) * 100;
    }
}


void arr_int_to_double(double* f_arr, int* i_arr, int dim){
    int i;
    for(i=0; i<dim; i++){
        f_arr[i] = (double) i_arr[i];
    }
}

void sorted_double_to_int(double* f_arr, int* i_arr,int* sorted_indexes, int dim){
    int i;
    for(i=0; i<dim; i++){
        i_arr[i] = (int) f_arr[sorted_indexes[i]];
    }
}


int main( int argc, char* argv[]){


    printf("%0.4f\n", 8.88886);


    //void index_print(int *indexes, int dim)
    srand ( time(NULL) );
    int dim = 10 ;
    double *arr;
    int *indexes;

    arr = malloc(sizeof(double) * (dim));
    //fill_random_doubles(arr, dim);

    arr[0]=7;
    arr[1]=7;
    arr[2]=5;
    arr[3]=1;
    arr[4]=7;
    arr[5]=5;
    arr[6]=1;
    arr[7]=2;
    arr[8]=1;
    arr[9]=7;


    indexes = malloc(sizeof(int) * (dim));
    fill_ordered_ints(indexes, dim);

    my_print(arr, indexes, dim);
    //quickSort(arr, indexes, dim);
    my_print(arr, indexes, dim);
    index_print(indexes,dim);

    stable(arr,indexes,dim);
    printf("\nafter stable:\n");
    my_print(arr, indexes, dim);
    index_print(indexes,dim);

}

