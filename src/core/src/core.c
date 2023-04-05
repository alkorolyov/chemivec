//
// Created by ergot on 09/03/2023.
//

#define NO_IMPORT_ARRAY // NumPy C-API is already imported
#include "core.h"

#define CREATE_NUMPY_VEC(np_input, dtype) \
    char ** in_data = numpy2cstr((PyArrayObject*)np_input); \
    int size = PyArray_SIZE((PyArrayObject*)np_input);      \
    npy_intp dims[] = {size};            \
    PyArrayObject* np_output = (PyArrayObject*)PyArray_EMPTY(1, dims, dtype, NPY_ARRAY_C_CONTIGUOUS); \
    npy_bool* out_data = (npy_bool*) PyArray_DATA(np_output); // output boolean array

#define RETURN_NUMPY_VEC \
    PyMem_Free(in_data); \
    PyArray_XDECREF(np_output); \
    return (PyObject*)np_output;

#define CREATE_THREAD_BATCH(size, n_jobs) \
        Batch* batch = malloc(sizeof(Batch)); \
        batch->sid = indigoAllocSessionId();  \
        batch->threadid = omp_get_thread_num(); \
        int batch_size = size / n_jobs; \
        int start_idx = batch->threadid * batch_size; \
        int end_idx = start_idx + batch_size; \
        if (batch->threadid == n_jobs - 1) end_idx = size; \
        batch->pinput = in_data + start_idx;  \
        batch->poutput = out_data + start_idx;\
        batch->size = end_idx - start_idx;


PyArrayObject *cstr2numpy(char **strings, int size) {
    npy_intp dims[] = {size};

    // create empty 1D numpy array of python objects
    PyArrayObject* arr = (PyArrayObject*)PyArray_EMPTY(1, dims, NPY_OBJECT, NPY_ARRAY_C_CONTIGUOUS);

    // copy strings to numpy array
    for (npy_intp i = 0; i < size; i++) {
        PyArray_SETITEM(arr, PyArray_GETPTR1(arr, i), PyUnicode_FromString(strings[i]));
    }
    return arr;
}

/***
 * Creates array of C strings from Numpy array of python strings
 * @param np_array numpy array of python unicode strings
 * @return pointer array of UTF8 strings
 */
char** numpy2cstr(PyArrayObject* np_array) {
    PyObject** pystr = PyArray_DATA(np_array);
    npy_intp size = PyArray_SIZE(np_array);
    char** cstr = PyMem_Malloc(size * sizeof(char*));
    for (npy_intp i = 0; i < size; i++) {
        cstr[i] = (char*)PyUnicode_AsUTF8(pystr[i]);
    }
    return cstr;
}

int checkReactionSmarts(char* smarts, qword sid){
    indigoSetSessionId(sid);
    int query = indigoLoadReactionSmartsFromString(smarts);
    if (query == -1) {
        return -1;
    }
    indigoFree(query);
    return 0;
}

/***
 * Reaction substructure search for single batch.
 * @param batch pointer to ReactionBatch object
 * @param query handle of indigo Query object
 * @param mode "DAYLIGHT-AAM" or ignored
 */
void reactionMatchBatch(Batch* batch, int query, const char *mode) {
    indigoSetSessionId(batch->sid);
    for (int i = 0; i < batch->size; i++) {
        int rxn = indigoLoadReactionFromString(batch->pinput[i]);
        if (rxn == -1) {
            printf("Invalid reaction SMILES: %s\n", batch->pinput[i]);
            batch->poutput[i] = NPY_FALSE;
            continue;
        }
        int matcher = indigoSubstructureMatcher(rxn, mode);
        int match = indigoMatch(matcher, query);
        if (match != 0)
            batch->poutput[i] = NPY_TRUE;
        else
            batch->poutput[i] = NPY_FALSE;
        //  printf("[%i %i]:\n in =  %s\n out = %i\n", batch->threadid, i, batch->pinput[i], batch->poutput[i]);
        indigoFree(rxn);
        indigoFree(matcher);
        indigoFree(match);
    }
}

void reactionMatchLin(char **in_data, npy_bool *out_data, int size, char *querySmarts, const char *mode) {
    // Single Thread
    
    Batch* batch = PyMem_Malloc(sizeof(Batch));
    batch->sid = indigoAllocSessionId();
    batch->threadid = 0;
    batch->pinput = in_data;
    batch->poutput = out_data;
    batch->size = size;

    int query = indigoLoadReactionSmartsFromString(querySmarts);
    if (query == -1) {
        printf("Invalid query SMARTS: %s", querySmarts);
        exit(EXIT_FAILURE);
    }
    indigoOptimize(query, NULL);

    reactionMatchBatch(batch, query, mode);

    indigoFree(query);
    indigoReleaseSessionId(batch->sid);
    PyMem_Free(batch);
}

/**
 * Vectorized version of reaction match. Creates new
 * boolean NumPy array of the same shape as an output.
 * @param in_data C array of reaction smiles
 * @param querySmarts string of query smarts
 * @param mode "DAYLIGHT-AAM" or NULL
 * @return
 */
void reactionMatchVec(char **in_data, npy_bool *out_data, int size, char *querySmarts, const char *mode, int n_jobs) {

    // Multi Thread

    // NO PYTHON FUNCTIONS HERE
    #pragma omp parallel num_threads(n_jobs)
    {
        // Create batch per each thread
        CREATE_THREAD_BATCH(size, n_jobs)
//        Batch* batch = malloc(sizeof(Batch));
//        batch->sid = indigoAllocSessionId();
//        batch->threadid = omp_get_thread_num();
//        int batch_size = size / n_jobs;
//        int start_idx = batch->threadid * batch_size;
//        int end_idx = start_idx + batch_size;
//        if (batch->threadid == n_jobs - 1) {
//            end_idx = size; // last batch
//        }
//        batch->pinput = in_data + start_idx;
//        batch->poutput = out_data + start_idx;
//        batch->size = end_idx - start_idx;

        // Create query object
        int query = indigoLoadReactionSmartsFromString(querySmarts);
        if (query == -1) {
            printf("Invalid reaction SMARTS %s\n", querySmarts);
            exit(EXIT_FAILURE);
        }
        indigoOptimize(query, NULL);

        reactionMatchBatch(batch, query, mode);

        indigoFree(query);
        indigoReleaseSessionId(batch->sid);
        free(batch);
    }

    return;
}

PyObject* reactionMatchNumPy(PyObject *np_input, char *querySmarts, char *aamMode, int n_jobs) {
    CREATE_NUMPY_VEC(np_input, NPY_BOOL)
//    char ** in_data = numpy2cstr((PyArrayObject*)np_input);

//    int size = PyArray_SIZE((PyArrayObject*)np_input);
//    npy_intp dims[] = {size};
//    PyArrayObject* np_output = (PyArrayObject*)PyArray_EMPTY(1, dims, NPY_BOOL, NPY_ARRAY_C_CONTIGUOUS);
//    npy_bool* out_data = (npy_bool*) PyArray_DATA(np_output); // output boolean array

    reactionMatchVec(in_data, out_data, size, querySmarts, aamMode, n_jobs);
//    reactionMatchLin(in_data, out_data, size, querySmarts, mode);

    RETURN_NUMPY_VEC
//    PyMem_Free(in_data);
//    PyArray_XDECREF(np_output);
//    return (PyObject*)np_output;
}

void structureMatchBatch(Batch* batch, int query, char* mode) {
    indigoSetSessionId(batch->sid);
    for (int i = 0; i < batch->size; i++) {
        int mol = indigoLoadMoleculeFromString(batch->pinput[i]);
        if (mol == -1) {
            printf("Invalid SMILES: %s\n", batch->pinput[i]);
            batch->poutput[i] = NPY_FALSE;
            continue;
        }
        int matcher = indigoSubstructureMatcher(mol, mode);
        int match = indigoMatch(matcher, query);
        if (match != 0)
            batch->poutput[i] = NPY_TRUE;
        else
            batch->poutput[i] = NPY_FALSE;
        //  printf("[%i %i]:\n in =  %s\n out = %i\n", batch->threadid, i, batch->pinput[i], batch->poutput[i]);
        indigoFree(mol);
        indigoFree(matcher);
        indigoFree(match);
    }
}

void structureMatchLin(char **in_data, npy_bool *out_data, int size, char *querySmarts, char *mode) {
    Batch* batch = PyMem_Malloc(sizeof(Batch));
    batch->sid = indigoAllocSessionId();
    batch->threadid = 0;
    batch->pinput = in_data;
    batch->poutput = out_data;
    batch->size = size;

    int query = indigoLoadSmartsFromString(querySmarts);
    if (query == -1) {
        printf("Invalid SMARTS %s\n", querySmarts);
        exit(EXIT_FAILURE);
    }

    indigoOptimize(query, NULL);

    structureMatchBatch(batch, query, mode);

    indigoFree(query);
    indigoReleaseSessionId(batch->sid);
    PyMem_Free(batch);
}

void structureMatchVec(char **in_data, npy_bool *out_data, int size, char *querySmarts, char *mode, int n_jobs) {
    // Multi Thread


    // ONLY THREAD SAFE FUNCTIONS HERE
    #pragma omp parallel num_threads(n_jobs)
    {
        // Create batch per each thread
        CREATE_THREAD_BATCH(size, n_jobs)
//        Batch* batch = malloc(sizeof(Batch));
//        batch->sid = indigoAllocSessionId();
//        batch->threadid = omp_get_thread_num();
//        int start_idx = batch->threadid * batch_size;
//        int end_idx = start_idx + batch_size;
//        if (batch->threadid == n_jobs - 1) {
//            end_idx = size; // last batch
//        }
//        batch->pinput = in_data + start_idx;
//        batch->poutput = out_data + start_idx;
//        batch->size = end_idx - start_idx;
        // Create query object
        int query = indigoLoadSmartsFromString(querySmarts);
        if (query == -1) {
            printf("Invalid SMARTS %s\n", querySmarts);
            exit(EXIT_FAILURE);
        }
        indigoOptimize(query, NULL);

        structureMatchBatch(batch, query, mode);

        indigoFree(query);
        indigoReleaseSessionId(batch->sid);
        free(batch);
    }

    return;

}

PyObject* structureMatchNumpy(PyObject *np_input, char* querySmarts, char *mode, int numCores) {
    CREATE_NUMPY_VEC(np_input, NPY_BOOL)
//    int size = PyArray_SIZE((PyArrayObject*)np_input);
//    npy_intp dims[] = {size};
//    char ** in_data = numpy2cstr((PyArrayObject*)np_input);
//
//    PyArrayObject* np_output = (PyArrayObject*)PyArray_EMPTY(1, dims, NPY_BOOL, NPY_ARRAY_C_CONTIGUOUS);
//    npy_bool* out_data = (npy_bool*) PyArray_DATA(np_output); // output boolean array

    structureMatchLin(in_data, out_data, size, querySmarts, mode);

    RETURN_NUMPY_VEC
//    PyMem_Free(in_data);
//    PyArray_XDECREF(np_output);
//    return (PyObject*)np_output;
}




