//
// Created by ergot on 09/03/2023.
//


#define NO_IMPORT_ARRAY // NumPy C-API is already imported
#include "vec.h"

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


int checkQuerySmarts(char* querySmarts, qword sid){
    indigoSetSessionId(sid);
    qword query = indigoLoadReactionSmartsFromString(querySmarts);
    char message[1024];
    if (query == -1) {
        snprintf(message, strlen(querySmarts) * sizeof(char) + 40, "Invalid SMARTS '%s' for reaction query\n", querySmarts);
        PyErr_SetString(PyExc_ValueError, message);
        return -1;
    }
    indigoFree(query);
//    indigoReleaseSessionId(sid);
    return 0;
}

/***
 * Reaction substructure search for single batch.
 * @param input C string array of rxn smiles
 * @param output boolean array of match results
 * @param size size of batch to process
 * @param query handle of indigo Query object
 * @param mode "DAYLIGHT-AAM" or NULL
 * @param id Indigo session id
 */
void reactionMatchBatch(struct ReactionBatch* batch, int query, const char *mode) {
    indigoSetSessionId(batch->sid);
    for (int i = 0; i < batch->size; i++) {
        int rxn = indigoLoadReactionFromString(batch->pinput[i]);
        if (rxn == -1) {
            printf("[%d %d] Invalid SMILES: %s\n", batch->threadid, i, batch->pinput[i]);
            batch->poutput[i] = NPY_FALSE;
            continue;
        }
        int matcher = indigoSubstructureMatcher(rxn, mode);
        int match = indigoMatch(matcher, query);
        if (match != 0)
            batch->poutput[i] = NPY_TRUE;
        else
            batch->poutput[i] = NPY_FALSE;
        // debug info
//        printf("[%i %i]:\n in =  %s\n out = %i\n", batch->threadid, i, batch->pinput[i], batch->poutput[i]);
        indigoFree(rxn);
        indigoFree(matcher);
        indigoFree(match);
    }
}


void reactionMatchLin(char **in_data, npy_bool *out_data, int size, char *querySmarts, const char *mode) {
    // Single Thread
    
    struct ReactionBatch* batch = PyMem_Malloc(sizeof(struct ReactionBatch));
    batch->sid = indigoAllocSessionId();
    batch->threadid = 0;
    batch->pinput = in_data;
    batch->poutput = out_data;
    batch->size = size;

    qword query = indigoLoadReactionSmartsFromString(querySmarts);
    if (query == -1) {
        printf("Invalid SMARTS: %s", querySmarts);
        return;
    }
    indigoOptimize(query, NULL);

    reactionMatchBatch(batch, query, mode);

    indigoFree(query);
    indigoReleaseSessionId(batch->sid);
    PyMem_Free(batch);
    return;
}

/**
 * Vectorized version of reaction match. Creates new
 * boolean NumPy array of the same shape as an output.
 * @param in_data C array of reaction smiles
 * @param querySmarts string of query smarts
 * @param mode "DAYLIGHT-AAM" or NULL
 * @return
 */
void reactionMatchVec(char **in_data, npy_bool *out_data, int size, char *querySmarts, const char *mode) {

//    if (checkQuerySmarts(querySmarts) == -1) {
//        return NULL;
//    };

    // Multi Thread
    int num_threads = omp_get_max_threads();
    int batch_size = size / num_threads;

    // NO PYTHON FUNCTIONS HERE
    #pragma omp parallel num_threads(num_threads)
    {
        // Create batch per each thread
        struct ReactionBatch* batch = malloc(sizeof(struct ReactionBatch));
        batch->sid = indigoAllocSessionId();
        batch->threadid = omp_get_thread_num();
        int start_idx = batch->threadid * batch_size;
        int end_idx = start_idx + batch_size;
        if (batch->threadid == num_threads - 1) {
            end_idx = size; // last batch
        }
        batch->pinput = in_data + start_idx;
        batch->poutput = out_data + start_idx;
        batch->size = end_idx - start_idx;

        // Create query object
        int query = indigoLoadReactionSmartsFromString(querySmarts);
        if (query == -1) {
            printf("Invalid SMARTS %s\n", querySmarts);
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


PyArrayObject *reactionMatchNumPy(PyArrayObject *np_input, char *querySmarts, char *aam_mode, qword moduleIndigoId) {
//    if (checkQuerySmarts(querySmarts) == -1) {
//        return NULL;
//    }

    // check query SMARTS
    if (checkQuerySmarts(querySmarts, moduleIndigoId) == -1) {
        return NULL;
    };


    int size = PyArray_SIZE(np_input);
    npy_intp dims[] = {size};
    char ** in_data = numpy2cstr(np_input);

    PyArrayObject* np_output = (PyArrayObject*)PyArray_EMPTY(1, dims, NPY_BOOL, NPY_ARRAY_C_CONTIGUOUS);
    npy_bool* out_data = (npy_bool*) PyArray_DATA(np_output); // output boolean array

    reactionMatchVec(in_data, out_data, size, querySmarts, aam_mode);
//    reactionMatchLin(in_data, out_data, size, querySmarts, aam_mode);

    PyMem_Free(in_data);
    PyArray_XDECREF(np_output);
    return np_output;
}


