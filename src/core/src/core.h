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

int check_reaction_smarts(char* smarts, qword sid);

void reaction_match_batch(Batch* batch, int query, const char *mode);

void reaction_match_lin(char **in_data, npy_bool *out_data, int size, char *querySmarts, const char *mode);

void reaction_match_vec(char **in_data, npy_bool *out_data, int size, char *querySmarts, const char *mode, int n_jobs);

PyObject* reaction_match_numpy(PyObject *np_input, char *querySmarts, char *aamMode, int n_jobs);

int check_structure_smarts(char* smarts, qword sid);

void structure_match_lin(char **in_data, npy_bool *out_data, int size, char *querySmarts, char *mode);

void structure_match_batch(Batch* batch, int query, char* mode);

void structure_match_vec(char **in_data, npy_bool *out_data, int size, char *querySmarts, char *mode, int n_jobs);

PyObject* structure_match_numpy(PyObject *np_input, char* querySmarts, char *mode, int n_jobs);