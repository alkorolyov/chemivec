//
// Created by Alexander Korolyov on 09/03/2023.
//

#define NO_IMPORT_ARRAY // NumPy C-API is already imported
#include "core.h"

/* Macro to create a NumPy array from an input array of strings */
#define CREATE_NUMPY_VEC(np_input, dtype) \
    char **in_data = numpy2cstr((PyArrayObject*)np_input); \
    int size = PyArray_SIZE((PyArrayObject*)np_input); \
    npy_intp dims[] = {size}; \
    PyArrayObject* np_output = (PyArrayObject*)PyArray_EMPTY(1, dims, dtype, NPY_ARRAY_C_CONTIGUOUS); \
    npy_bool* out_data = (npy_bool*) PyArray_DATA(np_output); /* output boolean array */

/* Macro to return a NumPy array and clean up resources */
#define RETURN_NUMPY_VEC \
    free(in_data); \
    PyArray_XDECREF(np_output); \
    return (PyObject*)np_output;

/* Macro to create a batch for threading */
#define CREATE_THREAD_BATCH(in_data, out_data, size) \
    Batch* batch = malloc(sizeof(Batch)); \
    if (!batch) { \
        perror("Failed to allocate memory for Batch"); \
        exit(EXIT_FAILURE); \
    } \
    batch->sid = indigoAllocSessionId(); \
    batch->threadid = omp_get_thread_num(); \
    int n_jobs = omp_get_num_threads(); \
    int batch_size = size / n_jobs; \
    int start_idx = batch->threadid * batch_size; \
    int end_idx = start_idx + batch_size; \
    if (batch->threadid == n_jobs - 1) end_idx = size; \
    batch->pinput = in_data + start_idx; \
    batch->poutput = out_data + start_idx; \
    batch->size = end_idx - start_idx;

/* Macro to free a batch */
#define FREE_BATCH(batch) \
    indigoReleaseSessionId(batch->sid); \
    free(batch);

/* Macro to create a query from SMARTS */
#define CREATE_QUERY_SMARTS(indigoSmartsFuncName, querySmarts) \
    int query = indigoSmartsFuncName(querySmarts); \
    if (query == -1) { \
        fprintf(stderr, "Invalid SMARTS: %s\n", querySmarts); \
        exit(EXIT_FAILURE); \
    } \
    indigoOptimize(query, NULL);

/* Macro to free a query */
#define FREE_QUERY(query) indigoFree(query);

/* Macro to perform matching in a batch for Structures and Reactions */
#define INDIGO_MATCH_BATCH(batch, indigoLoadFunc, indigoMatchFunc, query, mode) \
    indigoSetSessionId(batch->sid); \
    for (int i = 0; i < batch->size; i++) { \
        int obj = indigoLoadFunc(batch->pinput[i]); \
        if (obj == -1) { \
            fprintf(stderr, "Invalid SMILES: %s\n", batch->pinput[i]); \
            batch->poutput[i] = NPY_FALSE; \
            continue; \
        } \
        int matcher = indigoMatchFunc(obj, mode); \
        int match = indigoMatch(matcher, query); \
        batch->poutput[i] = (match != 0) ? NPY_TRUE : NPY_FALSE; \
        indigoFree(obj); \
        indigoFree(matcher); \
        indigoFree(match); \
    }

/* Function to convert C strings to NumPy array */
PyArrayObject *cstr2numpy(char **strings, int size) {
    npy_intp dims[] = {size};
    PyArrayObject* arr = (PyArrayObject*)PyArray_EMPTY(1, dims, NPY_OBJECT, NPY_ARRAY_C_CONTIGUOUS);

    for (npy_intp i = 0; i < size; i++) {
        PyArray_SETITEM(arr, PyArray_GETPTR1(arr, i), PyUnicode_FromString(strings[i]));
    }
    return arr;
}

/* Function to convert NumPy array of Python strings to C strings */
char** numpy2cstr(PyArrayObject* np_array) {
    PyObject** pystr = PyArray_DATA(np_array);
    npy_intp size = PyArray_SIZE(np_array);
    char** cstr = malloc(size * sizeof(char*));
    if (!cstr) {
        perror("Failed to allocate memory for C strings");
        exit(EXIT_FAILURE);
    }
    for (npy_intp i = 0; i < size; i++) {
        cstr[i] = (char*)PyUnicode_AsUTF8(pystr[i]);
    }
    return cstr;
}

/* Function to check reaction SMARTS */
int check_reaction_smarts(char* smarts, qword sid) {
    indigoSetSessionId(sid);
    int query = indigoLoadReactionSmartsFromString(smarts);
    if (query == -1) {
        return -1;
    }
    indigoFree(query);
    return 0;
}

/* Function to check structure SMARTS */
int check_structure_smarts(char* smarts, qword sid) {
    indigoSetSessionId(sid);
    int query = indigoLoadSmartsFromString(smarts);
    if (query == -1) {
        return -1;
    }
    indigoFree(query);
    return 0;
}

/* Function for batch reaction matching */
void reaction_match_batch(Batch* batch, int query, const char *mode) {
    INDIGO_MATCH_BATCH(batch, indigoLoadReactionFromString, indigoSubstructureMatcher, query, mode)
}

/* Linear reaction matching function */
void reaction_match_lin(char **in_data, npy_bool *out_data, int size, char *querySmarts, const char *mode) {
    Batch* batch = malloc(sizeof(Batch));
    if (!batch) {
        perror("Failed to allocate memory for Batch");
        exit(EXIT_FAILURE);
    }
    batch->sid = indigoAllocSessionId();
    batch->threadid = 0;
    batch->pinput = in_data;
    batch->poutput = out_data;
    batch->size = size;

    CREATE_QUERY_SMARTS(indigoLoadReactionSmartsFromString, querySmarts)
    reaction_match_batch(batch, query, mode);
    FREE_QUERY(query)
    FREE_BATCH(batch)
}

/* Vectorized reaction matching function */
void reaction_match_vec(char **in_data, npy_bool *out_data, int size, char *querySmarts, const char *mode, int n_jobs) {
    #pragma omp parallel num_threads(n_jobs)
    {
        CREATE_THREAD_BATCH(in_data, out_data, size)
        CREATE_QUERY_SMARTS(indigoLoadReactionSmartsFromString, querySmarts)
        reaction_match_batch(batch, query, mode);
        FREE_QUERY(query)
        FREE_BATCH(batch)
    }
}

/* Main function for reaction matching with NumPy input */
PyObject* reaction_match_numpy(PyObject *np_input, char *querySmarts, char *aamMode, int n_jobs) {
    CREATE_NUMPY_VEC(np_input, NPY_BOOL)
    reaction_match_vec(in_data, out_data, size, querySmarts, aamMode, n_jobs);
    RETURN_NUMPY_VEC
}

/* Batch structure matching function */
void structure_match_batch(Batch* batch, int query, char* mode) {
    INDIGO_MATCH_BATCH(batch, indigoLoadMoleculeFromString, indigoSubstructureMatcher, query, mode)
}

/* Linear structure matching function */
void structure_match_lin(char **in_data, npy_bool *out_data, int size, char *querySmarts, char *mode) {
    Batch* batch = malloc(sizeof(Batch));
    if (!batch) {
        perror("Failed to allocate memory for Batch");
        exit(EXIT_FAILURE);
    }
    batch->sid = indigoAllocSessionId();
    batch->threadid = 0;
    batch->pinput = in_data;
    batch->poutput = out_data;
    batch->size = size;

    CREATE_QUERY_SMARTS(indigoLoadSmartsFromString, querySmarts)
    structure_match_batch(batch, query, mode);
    FREE_QUERY(query)
    FREE_BATCH(batch)
}

/* Vectorized structure matching function */
void structure_match_vec(char **in_data, npy_bool *out_data, int size, char *querySmarts, char *mode, int n_jobs) {
    #pragma omp parallel num_threads(n_jobs)
    {
        CREATE_THREAD_BATCH(in_data, out_data, size)
        CREATE_QUERY_SMARTS(indigoLoadSmartsFromString, querySmarts)
        structure_match_batch(batch, query, mode);
        FREE_QUERY(query)
        FREE_BATCH(batch)
    }
}

/* Main function for structure matching with NumPy input */
PyObject* structure_match_numpy(PyObject *np_input, char* querySmarts, char *mode, int n_jobs) {
    CREATE_NUMPY_VEC(np_input, NPY_BOOL)
    structure_match_vec(in_data, out_data, size, querySmarts, mode, n_jobs);
    RETURN_NUMPY_VEC
}
