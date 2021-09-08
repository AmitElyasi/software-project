#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int pyMat_to_C_array(PyObject*, double*, int);
static PyObject* c_array_to_pyMat(double*, int, int);


/*utils functions*/
static int pyMat_to_C_array(PyObject* pyMat, double* mat, int dim){
    int i,j,m,n;
    PyObject* pyVec = PyList_GetItem(pyMat, 0);
    PyObject* pyItem = PyList_GetItem(pyVec, 0);
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
            set(mat, i, j, dim, PyFloat_AsDouble(pyItem));
    
            if (get(mat, i, j, dim) == -1 && PyErr_Occurred()){
                /* float too big to fit in a C double, bail out */
                puts("An Error Has Occured");
                return 0;
            }
        }
    }
    return 1;
}
/*return pyList ocject in the shape of (n,m) */
static PyObject* c_array_to_pyMat(double* mat, int n, int m){
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
    double *data_points,*target_matrix;
    int rows ,cols;
    PyObject *pymat;


    data_points = malloc(sizeof(double) * n * dim);
    
    /*convert python mat to c array*/
    pyMat_to_C_array(pyData_points, data_points, dim);
    printf("%s \n", goal);
    print_matrix(data_points, n, dim);

    if(strcmp(goal, "wam")){
        target_matrix = wam(data_points, n, dim);
        rows = n;
        cols = n;
    }
    else if(strcmp(goal, "ddg")){
        target_matrix = ddg(data_points, n, dim);
        rows = 1;
        cols = n;
    }
    else if(strcmp(goal, "lnorm")){
        target_matrix = lnorm(data_points, n, dim);
        rows = n;
        cols = n;
    }
    else if(strcmp(goal, "jacobi")){
        target_matrix = jacobi(data_points, n);
        rows = (n+1);
        cols = n;
    }else{
        target_matrix = spk(data_points, n, dim, &k);
        rows = k;
        cols = k;
    }



    pymat = c_array_to_pyMat(target_matrix, rows, cols);
    free(target_matrix);
    
    return pymat; 
}


static PyObject *calc_transformation_matrix_capi(PyObject *self, PyObject* args){
    PyObject *pyData_points;
    int k, dim, n;
    char *goal;
    if(!PyArg_ParseTuple(args, "isOii", &k,
                                        &goal,
                                        &pyData_points,
                                        &dim,
                                        &n)){
        return NULL;
    }
    return calc_transformation_matrix(k, goal, pyData_points, dim, n);
}

static PyObject *fit_c(int k, PyObject *pyData_points, PyObject *pyCentroid, int max_iter, int n){
    double *data_points, *centroid, *utl;
    PyObject *pyMat;

    data_points = malloc(sizeof(double) *((k+1) *n));
    centroid = malloc(sizeof(double) * (k * k));
    utl = malloc(sizeof(double) * (k * (k+1)));

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
    int k, n, max_iter;
    if(!PyArg_ParseTuple(args, "iOOii", &k,
                                        &pyData_points,
                                        &pyCentroid,
                                        &max_iter,
                                        &n)){
        return NULL;
    }
    return fit_c(k, pyData_points, pyCentroid, max_iter, n);

}

/*
 * This array tells Python what methods this module has.
 * We will use it in the next structure
 */
static PyMethodDef _spkmeansMethods[] = {
    {"calc_transformation_matrix",                   /* the Python method name that will be used */
      (PyCFunction) calc_transformation_matrix_capi, /* the C-function that implements the Python function and returns static PyObject*  */
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
    _spkmeansMethods /* the PyMethodDef array from before containing the methods of the extension */
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
