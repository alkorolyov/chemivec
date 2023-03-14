//
// Created by ergot on 09/03/2023.
//

#ifndef CHEMIVEC_VEC_H
#define CHEMIVEC_VEC_H
#endif //CHEMIVEC_VEC_H
#include "Python.h"
#include "indigo.h"
#include <omp.h>

#define PY_ARRAY_UNIQUE_SYMBOL CHEMIVEC_ARRAY_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

struct ReactionBatch{
    char** pin;
    npy_bool* pout;
    int size;
    qword sid;
    int threadNum;
};

PyArrayObject* str2numpy(char** strings, int size);

char** numpy2str(PyArrayObject * np_array);

void reactionMatchBatch(struct ReactionBatch* batch, int query, const char *mode);

PyArrayObject* reactionMatchVec(PyArrayObject* np_input, char* querySmarts, const char* mode);

inline static void finishSearch(int rxn, int matcher, int match) {
    indigoFree(rxn);
    indigoFree(matcher);
    indigoFree(match);
}


