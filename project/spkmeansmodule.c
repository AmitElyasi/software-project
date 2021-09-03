#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>


/*utils functions*/
int pyMat_to_C_array(PyObject* pyMat, float* mat, int dim){
    int i,j,m,n;
    PyObject* pyVec;
    PyObject* pyItem;
    /* Is it a list? */
    if (!PyList_Check(pyMat))
        return 0;
    /* Get the size of it and build the output list */
    n = PyList_Size(pyMat);  /*  Same as in Python len(_list)  */
    /* Go over each item of the list and reduce it */
    for (i = 0; i < n; i++) {
        pyVec = PyList_GetItem(pyMat, i);
        if (!PyList_Check(pyVec)){  /* We only print lists */
            return 0;
        }
        m = PyList_Size(pyVec);
        for (j = 0; j < m; j++) {
            pyItem = PyList_GetItem(pyVec, j);
            set(mat, i, j , dim, PyFloat_AsDouble(pyItem));

            if (get(mat, i, j, dim) == -1 && PyErr_Occurred()){
                /* float too big to fit in a C double, bail out */
                puts("Error parsing a list to C matrix");
                return 0;
            }
        }
    }
    return 1;
}
/*return pyList ocject in the shape of (n,m) */
PyObject* c_array_to_pyMat(float* mat, int n, int m){
    int i, j;
    PyObject *pyItem, *pyVec, *pyMat;
    pyMat = PyList_New(0);
    for (i=0; i < n; i++){
        pyVec = PyList_New(0);
        for (j = 0; j < m; j++){
            pyItem = Py_BuildValue("d", get(mat, i, j, m));
            PyList_Append(pyVec, pyItem);
        }
        PyList_Append(pyMat, pyVec);
    }
    return pyMat;
}


static PyObject *calc_transformation_matrix(int k, char *goal, PyObject *pyData_points, int dim, int n){
    float x,*data_points, *weighted_adj_mat, *diagonal_mat, *normalized_laplacian, *eigonvalues, *V, *T;
    int i,j, index,*indexes;
    PyObject *item, *pyvec, *pymat;

    data_points = malloc(sizeof(float) * n * dim);
    weighted_adj_mat = malloc(sizeof(float) * n * n);
    
    /*convert python mat to c array*/
    pyMat_to_C_array(pyData_points, weighted_adj_mat, n);
    if(!strcmp(goal, "jacobi")){
        form_weighted_adj_mat(weighted_adj_mat, data_points, dim, n);
        if(strcmp(goal, "wam")){
            return c_array_to_pyMat(weighted_adj_mat, n, n);
        }
        diagonal_mat = malloc(sizeof(float) * n);
        form_diagonal_mat(diagonal_mat, weighted_adj_mat, n);
        if(strcmp(goal, "ddg")){
            return c_array_to_pyMat(diagonal_mat, n, n);
        }
        normalized_laplacian = malloc(sizeof(float) * n * n);
        calc_normalized_laplacian(normalized_laplacian, diagonal_mat, weighted_adj_mat, dim);
        if(strcmp(goal, "lnorm")){
            return c_array_to_pyMat(normalized_laplacian, n, n);
        }
        free(weighted_adj_mat);
        free(diagonal_mat);
    }
    jacobi_algorithm_for_eigenvalues(normalized_laplacian, V, n);
    if(strcmp(goal,"jacobi")){
        pymat = PyList_New(0);
        pyvec = PyList_New(0);
        for(i=0;i<n;i++){
            item = Py_BuildValue("d", get(normalized_laplacian, i, i, n));
            PyList_Append(pyvec, item);
        }
        PyList_Append(pymat,pyvec);
        for (i=0;i<n; i++){
            pyvec = PyList_New(0);
            for (j=0;j<n;j++){
                item = Py_BuildValue("d", get(V, i, j, n));
                PyList_Append(pyvec, item);
            }
            PyList_Append(pymat, pyvec);
        }
        return pymat;
    }
    eigonvalues = malloc(sizeof(float) *n);
    indexes = malloc(sizeof(int) * n);
    diag_to_array(normalized_laplacian, eigonvalues);
    quickSort_indexes(eigonvalues, indexes, n);
    if(k == 0){
        k = calc_eiganvalue_gap;
    }
       T = malloc(sizeof(float) * n * k);
    free(normalized_laplacian);
    free(V);
    normalized_k_eigonvectors(indexes, n, k, T);
    pymat = c_array_to_pyMat(T, n, k);
    free(eigonvalues);
    free(indexes);
    free(T);
    
    return pymat; 
}


static PyObject *calc_transformation_matrix_capi(PyObject *self, PyObject* args){
    PyObject *pyData_points, *pyCentroid;
    int k, dim, n, max_iter;
    char *goal;
    if(!PyArg_ParseTuple(args, "isOii", &k,
                                        &goal,
                                        &pyData_points,
                                        &dim,
                                        &n)){
        return NULL;
    }
    return calc_tranformation_matrix(k, goal, pyData_points, dim, n);
    
}

static PyObject * fit_c(int k, PyObject *pyData_points, PyObject *pyCentroid, int max_iter, int n){
    float *data_points, *centroid, *utl,x;
    Py_ssize_t index;
    int i,j;
    PyObject *item, *pyMat;

    data_points = malloc(sizeof(float) *((k+1) *n));
    centroid = malloc(sizeof(float) * (k * k));
    utl = malloc(sizeof(float) * (k * (k+1)));

    /*convert python lists : pyData_points and pyCentroid to c arrays*/
    pyMat_to_C_array(pyData_points, data_points, k);
    pyMat_to_C_array(pyCentroid, centroid, k);

    kmeans(k,data_points, centroid, utl, max_iter, k, n);

    pyMat = c_array_to_pyMat(centroid, k, k);

    free(data_points);
    free(centroid);
    free(utl);

    return pyMat;
}


static PyObject *fit_capi(PyObject* self, PyObject* args){
    PyObject *pyData_points, *pyCentroid;
    int k, dim, n, max_iter;
    if(!PyArg_ParseTuple(args, "iOOii", &k,
                                        &pyData_points,
                                        &pyCentroid,
                                        &max_iter,
                                        &n)){
        return NULL;
    }
    return fit_c(k, pyData_points, pyCentroid, max_iter, dim, n);

}

/*
 * This array tells Python what methods this module has.
 * We will use it in the next structure
 */
static PyMethodDef spkmeansMethods[] = {
    {"calc_tranformation_matrix",                   /* the Python method name that will be used */
      (PyCFunction) calc_tranformation_matrix_capi, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters accepted for this function */
      PyDoc_STR("calcu")}, /*  The docstring for the function */
    {"fit",                   /* the Python method name that will be used */
      (PyCFunction) fit_capi, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters accepted for this function */
      PyDoc_STR("compute the new matrix under the first k eigenvector base")}, /*  The docstring for the function */
    {NULL, NULL, 0, NULL}     /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};


/* This initiates the module using the above definitions. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeans", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    spkmeansMethods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC
PyInit_spkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}
