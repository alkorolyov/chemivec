//
// Created by ergot on 09/03/2023.
//

#ifndef CHEMIVEC_CORE_H
#define CHEMIVEC_CORE_H
#endif //CHEMIVEC_CORE_H

#include "Python.h"
#include "indigo.h"
#include "omp.h"

#define PY_ARRAY_UNIQUE_SYMBOL CHEMIVEC_ARRAY_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"


typedef struct {
    qword sid;
    long n_jobs;
} ChemivecOptions;

typedef struct {
    char** pinput;
    npy_bool* poutput;
    int size;
    qword sid;
    int threadid;
} Batch;

PyArrayObject* cstr2numpy(char** strings, int size);

char** numpy2cstr(PyArrayObject * np_array);

int checkReactionSmarts(char* smarts, qword sid);

void reactionMatchBatch(Batch* batch, int query, const char *mode);

void reactionMatchLin(char **in_data, npy_bool *out_data, int size, char *querySmarts, const char *mode);

void reactionMatchVec(char **in_data, npy_bool *out_data, int size, char *querySmarts, const char *mode, int n_jobs);

PyObject* reactionMatchNumPy(PyObject *np_input, char *querySmarts, char *aamMode, int n_jobs);

void structureMatchLin(char **in_data, npy_bool *out_data, int size, char *querySmarts, char *mode);

void structureMatchBatch(Batch* batch, int query, char* mode);

void structureMatchVec(char **in_data, npy_bool *out_data, int size, char *querySmarts, char *mode, int n_jobs);

inline static void finishSearch(int rxn, int matcher, int match) {
    indigoFree(rxn);
    indigoFree(matcher);
    indigoFree(match);
}

#define CREATE_NPY_BOOL_ARRAY(np_input, size)\
    int size = PyArray_SIZE((PyArrayObject*)np_input);\
    npy_intp dims[] = {size};\
    char ** in_data = numpy2cstr((PyArrayObject*)np_input);\
    PyArrayObject* np_output = (PyArrayObject*)PyArray_EMPTY(1, dims, NPY_BOOL, NPY_ARRAY_C_CONTIGUOUS);\
    npy_bool* out_data = (npy_bool*) PyArray_DATA(np_output); // output boolean array

#define RETURN_NPY_BOOL(in_data, np_output)\
    PyMem_Free(in_data);\
    PyArray_XDECREF(np_output);\
    return (PyObject*)np_output;

