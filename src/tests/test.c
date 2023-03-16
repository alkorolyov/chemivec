//
// Created by ergot on 11/03/2023.
//

#include <stdio.h>

#define PY_ARRAY_UNIQUE_SYMBOL CHEMIVEC_ARRAY_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include "vec.h"
#include <assert.h>

void test_numpy_c_string() {
    char* in_strings[] = {"string1", "", "string2", "string3"};
    int size = 4;
    PyArrayObject* np_strings = cstr2numpy(in_strings, size);
    char** out_strings = numpy2cstr(np_strings);

    for (int i = 0; i < size; i++) {
//        printf("%s %s\n", in_strings[i], out_strings[i]);
        assert(in_strings[i] == out_strings[i]);
    }
}

void test_reaction_batch(qword id) {
    char* rxn_smi[] = {"[C:1](=O)C>>[C:1](O)C",
                         "C(=C)C>>C(O)C",
                         "[C:2]=O>>[C:2]O",
                         "C=O>>CO"};
    int size = 4;
    npy_intp dims[] = {size};
    PyArrayObject* np_output = (PyArrayObject*)PyArray_EMPTY(1, dims, NPY_BOOL, NPY_ARRAY_C_CONTIGUOUS);


    char* querySmarts = "[C:1]=[O]>>[C:1]-[OX2]";
    qword query = indigoLoadReactionSmartsFromString(querySmarts);
    indigoOptimize(query, NULL);
    struct ReactionBatch batch;
    batch.pin = rxn_smi;
    batch.pout = PyArray_DATA(np_output);
    batch.size = size;
    batch.sid = id;
    batch.threadNum = 0;

    npy_bool correct_result[] = {1, 0, 1, 0};
    reactionMatchBatch(&batch, query, "DAYLIGHT-AAM");
    for (int i = 0; i < size; i++) {
//        printf("%s - [%i]\n", batch.pin[i], batch.pout[i]);
        assert(batch.pout[i] == correct_result[i]);
    }

}

void test_reaction_vec() {
    char* rxn_smi[] = {"[C:1](=O)C>>[C:1](O)C",
                       "C(=C)C>>C(O)C",
                       "[C:2]=O>>[C:2]O",
                       "C=O>>CO"};
    int size = 4;
    npy_intp dims[] = {size};
    char* querySmarts = "[C:1]=[O]>>[C:1]-[OX2]";
    PyArrayObject* np_result = reactionMatchVec(rxn_smi, size, querySmarts, "DAYLIGHT-AAM");
    npy_bool* result = PyArray_DATA(np_result);
    npy_bool correct_result[] = {1, 0, 1, 0};
    for (int i = 0; i < size; i++) {
//        printf("%s - [%i]\n", batch.pin[i], batch.pout[i]);
        assert(batch.pout[i] == correct_result[i]);
    }
}

int main() {
    Py_Initialize();
    import_array()

    qword id = indigoAllocSessionId();
    if (id == -1) {
        printf("indigoAllocSessionId failed");
        exit(EXIT_FAILURE);
    }

    test_numpy_c_string();
    test_reaction_batch(id);
    test_reaction_vec();

    if (indigoCountReferences() > 0) {
        indigoFreeAllObjects();
    }
    indigoReleaseSessionId(id);
    Py_Finalize();
    return EXIT_SUCCESS;
}
