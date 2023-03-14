//
// Created by ergot on 09/03/2023.
//

#define NO_IMPORT_ARRAY // NumPy C-API is already imported
#include "vec.h"

PyArrayObject *str2numpy(char **strings, int size) {
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
char** numpy2str(PyArrayObject* np_array) {
    PyObject** pystr = PyArray_DATA(np_array);
    int size = PyArray_SIZE(np_array);
    char** cstr = malloc(size * sizeof(char*));
    for (int i = 0; i < PyArray_SIZE(np_array); i++) {
        cstr[i] = (char*)PyUnicode_AsUTF8(pystr[i]);
    }
    return cstr;
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
        int rxn = indigoLoadReactionFromString(batch->pin[i]);
        int matcher = indigoSubstructureMatcher(rxn, mode);
        int match = indigoMatch(matcher, query);
        if (match != 0)
            batch->pout[i] = NPY_TRUE;
        else
            batch->pout[i] = NPY_FALSE;
//        printf("thr[%i]:\n in = %s\n out[%i] = %i\n", batch->threadNum, batch->pin[i], i, batch->pout[i]);
        finishSearch(rxn, matcher, match);
    }
}

/**
 * Vectorized version of reaction match. Creates new
 * boolean NumPy array of the same shape as an output.
 * @param np_input numpy array of reaction smiles
 * @param querySmarts string of query smarts
 * @param mode "DAYLIGHT-AAM" or NULL
 * @return
 */
PyArrayObject* reactionMatchVec(PyArrayObject* np_input, char* querySmarts, const char* mode) {
    // Create output array of the same shape as np_input

    int size = PyArray_SIZE(np_input);
    PyArrayObject* np_output = (PyArrayObject*)PyArray_EMPTY(1, PyArray_DIMS(np_input), NPY_BOOL, NPY_ARRAY_C_CONTIGUOUS);

    // Get pointers to input and output data
    char** in_data = numpy2str(np_input); // C strings of rxn_smiles
    npy_bool* out_data = (npy_bool*) PyArray_DATA(np_output); // output boolean array

    // Single Thread
/*
    qword query = indigoLoadReactionSmartsFromString(querySmarts);
    indigoOptimize(query, NULL);
    struct ReactionBatch batch;
    batch.pin = in_data;
    batch.pout = PyArray_DATA(np_output);
    batch.size = size;
    batch.threadNum = 0;
    batch.sid = indigoAllocSessionId();

    reactionMatchBatch(&batch, query, mode);

    indigoFree(query);
    indigoReleaseSessionId(batch.sid);
*/


    // Multi Thread
    int num_threads = omp_get_max_threads();
    int batch_size = size / num_threads;

    #pragma omp parallel num_threads(num_threads)
    {
        // Create batch per each thread
        struct ReactionBatch* batch = malloc(sizeof(struct ReactionBatch));
        batch->sid = indigoAllocSessionId();
        batch->threadNum = omp_get_thread_num();

        int start_idx = batch->threadNum * batch_size;
        int end_idx = start_idx + batch_size;
        batch->pin = in_data + start_idx;
        batch->pout = out_data + start_idx;
        // last batch
        if (batch->threadNum == num_threads - 1) {
            end_idx = size;
        }
        batch->size = end_idx - start_idx;

        // Create query object
        int query = indigoLoadReactionSmartsFromString(querySmarts);
        indigoOptimize(query, NULL);

        reactionMatchBatch(batch, query, mode);

        indigoFree(query);
        indigoReleaseSessionId(batch->sid);
        free(batch);
    }

//    printf("%s", indigoGetLastError());
    free(in_data);

    return np_output;
}



