//
// Created by ergot on 11/03/2023.
//

#include <stdio.h>

#define PY_ARRAY_UNIQUE_SYMBOL CHEMIVEC_ARRAY_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include "core.h"
#include "unity.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */
qword sid; // indigo id for test session

void test_cstr_to_numpy() {
    char* c_strings[] = {"a", "str1", ""};
    int size = 3;
    PyArrayObject* np_strings = cstr2numpy(c_strings, size);
    PyObject**  data = PyArray_DATA(np_strings);
    TEST_ASSERT_EQUAL_STRING(c_strings[0], PyUnicode_AsUTF8(data[0]));
    TEST_ASSERT_EQUAL_STRING(c_strings[1], PyUnicode_AsUTF8(data[1]));
    TEST_ASSERT_EQUAL_STRING(c_strings[2], PyUnicode_AsUTF8(data[2]));

    PyArray_XDECREF(np_strings);
}

void test_numpy_to_cstr() {
    // create numpy str array
    npy_intp dims[1] = {2};
    PyArrayObject* np_array = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_OBJECT);
    PyArray_SETITEM(np_array, PyArray_GETPTR1(np_array, 0), PyUnicode_FromString("foo"));
    PyArray_SETITEM(np_array, PyArray_GETPTR1(np_array, 1), PyUnicode_FromString("bar"));

    char** c_strings = numpy2cstr(np_array);
    TEST_ASSERT_EQUAL_STRING(c_strings[0], "foo");
    TEST_ASSERT_EQUAL_STRING(c_strings[1], "bar");

    free(c_strings);
    PyArray_XDECREF(np_array);
}

void test_reaction_batch() {
    char* input[] = {"[C:1]=O>>[C:1]O",
                     "C=O>>CO"};
    npy_bool output[2];
    int size = 2;

    Batch batch;
    batch.sid = sid;
    batch.pinput = input;
    batch.poutput = output;
    batch.size = size;
    batch.threadid = 0;

    indigoSetSessionId(sid);
    const char* querySmarts = "[C:1]=O>>[C:1]O";
    int query = indigoLoadReactionSmartsFromString(querySmarts);
    indigoOptimize(query, NULL);

    reaction_match_batch(&batch, query, "DAYLIGHT-AAM");
    indigoFree(query);

    TEST_ASSERT_EQUAL(output[0], 1);
    TEST_ASSERT_EQUAL(output[1], 0);
}

void test_reaction_lin() {
    char* input[] = {"[C:1]=O>>[C:1]O",
                       "C=O>>CO"
                    };
    npy_bool output[2];
    int size = 2;

    char* querySmarts = "[C:1]=O>>[C:1]O";
    reaction_match_lin(input, output, size, querySmarts, "DAYLIGHT-AAM");
    TEST_ASSERT_EQUAL(output[0], 1);
    TEST_ASSERT_EQUAL(output[1], 0);
}

void test_reaction_vec() {
    char* input[] = {"[C:1]=O>>[C:1]O",
                       "C=O>>CO"
                    };
    npy_bool output[2];
    int size = 2;

    char* querySmarts = "[C:1]=O>>[C:1]O";
    int max_cores = omp_get_max_threads();
    reaction_match_vec(input, output, size, querySmarts, "DAYLIGHT-AAM", max_cores);
    TEST_ASSERT_EQUAL(output[0], 1);
    TEST_ASSERT_EQUAL(output[1], 0);
}



void test_incorrect_smi_batch() {
    char* input[] = {"[C:1=O>>[C:1]O"};
    npy_bool output[] = {1};
    int size = 1;

    Batch batch;
    batch.sid = sid;
    batch.pinput = input;
    batch.poutput = output;
    batch.size = size;
    batch.threadid = 0;

    indigoSetSessionId(sid);
    const char* querySmarts = "CO>>";
    int query = indigoLoadReactionSmartsFromString(querySmarts);
    indigoOptimize(query, NULL);

    reaction_match_batch(&batch, query, "DAYLIGHT-AAM");
    indigoFree(query);
    TEST_ASSERT_EQUAL(output[0], 0);

}

void test_incorrect_smi_vec() {
    char* input[] = {"[C>>]"};
    npy_bool output[] = {1};
    int size = 1;

    char* querySmarts = "CO>>";
    int max_cores = omp_get_max_threads();
    reaction_match_vec(input, output, size, querySmarts, "DAYLIGHT-AAM", max_cores);
    TEST_ASSERT_EQUAL(output[0], 0);
}

void test_reaction_smarts() {
    indigoSetSessionId(sid);
    TEST_ASSERT_EQUAL(check_reaction_smarts("CO>", sid), -1);
    TEST_ASSERT_EQUAL(check_reaction_smarts("CO>>", sid), 0);
    TEST_ASSERT_EQUAL(check_reaction_smarts("[CH3]>>", sid), 0);
}


void test_substructure_batch() {
    char* input[] = {"CC", "CO"};
    npy_bool output[2];
    int size = 2;

    Batch batch;
    batch.sid = sid;
    batch.pinput = input;
    batch.poutput = output;
    batch.size = size;
    batch.threadid = 0;

    indigoSetSessionId(sid);
    const char* querySmarts = "CC";
    int query = indigoLoadSmartsFromString(querySmarts);
    indigoOptimize(query, NULL);

    structure_match_batch(&batch, query, NULL);
    indigoFree(query);
    TEST_ASSERT_EQUAL(output[0], 1);
    TEST_ASSERT_EQUAL(output[1], 0);
}

void test_substructure_lin() {
    char* input[] = {"CC", "CO"};
    npy_bool output[2];
    int size = 2;

    char* querySmarts = "CC";
    structure_match_lin(input, output, size, querySmarts, NULL);
    TEST_ASSERT_EQUAL(output[0], 1);
    TEST_ASSERT_EQUAL(output[1], 0);
}

void test_substructure_vec() {
    char* input[] = {"CC", "CO"};
    npy_bool output[2];
    int size = 2;

    char* querySmarts = "CC";
    int max_cores = omp_get_max_threads();
    structure_match_vec(input, output, size, querySmarts, NULL, max_cores);
    TEST_ASSERT_EQUAL(output[0], 1);
    TEST_ASSERT_EQUAL(output[1], 0);
}


int main(void) {
    UNITY_BEGIN();

    Py_Initialize();
    import_array()
    
    sid = indigoAllocSessionId();

    RUN_TEST(test_cstr_to_numpy);
    RUN_TEST(test_numpy_to_cstr);
    RUN_TEST(test_reaction_batch);
    RUN_TEST(test_reaction_lin);
    RUN_TEST(test_reaction_vec);
    RUN_TEST(test_incorrect_smi_batch);
    RUN_TEST(test_incorrect_smi_vec);
    RUN_TEST(test_reaction_smarts);
    RUN_TEST(test_substructure_batch);
    RUN_TEST(test_substructure_lin);
    RUN_TEST(test_substructure_vec);

    if (indigoCountReferences() > 0) {
        indigoFreeAllObjects();
    }

    indigoReleaseSessionId(sid);
    Py_Finalize();
    return UNITY_END();
}
