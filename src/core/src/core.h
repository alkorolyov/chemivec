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

PyObject* structureMatchNumpy(PyObject *np_input, char* querySmarts, char *mode, int numCores);