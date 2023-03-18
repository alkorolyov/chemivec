//
// Created by Alex on 09/03/2023.
//

#define PY_SSIZE_T_CLEAN
#include "Python.h"

#define PY_ARRAY_UNIQUE_SYMBOL CHEMIVEC_ARRAY_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#include "core.h"


static ChemivecOptions* options;

static int initOptions(qword SessionId) {
    options = (ChemivecOptions*) PyMem_Malloc(sizeof(options));
    if (options == NULL) {
        return -1;
    }
    options->sid = SessionId;
    options->num_cores = omp_get_max_threads();
    return 0;
}

// Define the module deallocation function
static void freeModule() {
    // Indigo
    indigoSetSessionId(options->sid);
    if (indigoCountReferences() > 0) {
        indigoFreeAllObjects();
    }
    indigoReleaseSessionId(options->sid);

    // Module Options
    PyMem_Free(options);
}


/*


PyArrayObject* createNumpyArrFromString(int size, const char** strings) {
    npy_intp dims[] = {size};

    // create empty 1D numpy array of python objects
    PyArrayObject* arr = (PyArrayObject*)PyArray_EMPTY(1, dims, NPY_OBJECT, NPY_ARRAY_C_CONTIGUOUS);

    // copy strings to numpy array
    for (npy_intp i = 0; i < size; i++) {
        PyArray_SETITEM(arr, PyArray_GETPTR1(arr, i), PyUnicode_FromString(strings[i]));
    }
    return arr;
}

qword initIndigo() {
    printf("Starting indigo session ... ");

    // start session
    qword id = indigoAllocSessionId();

    if (id == -1) {
        printf("FAILED\n");
        exit(1);
    }

    printf("OK \nversion:   %s\n", indigoVersion());
    return id;
}


int main() {

    // Initialize Python interpreter
    Py_Initialize();

    // Initialize NumPy
    initNumpy();

    PyObject* pyobj = PyUnicode_FromString("Hello World!");

    char* buf = calloc(SMILE_BUF_LEN, sizeof(char));
    int len = unicodeAsUTF8(pyobj, buf);
    if (len > 0) {
        printf("length:%i string: %.5s", len, buf);
    }
    free(buf);



    // Initialize Indigo
    qword sessionId = initIndigo();


    // Define C array of string in_data
//    const char* strings[] = {"[C:1](=O)C>>[C:1](O)C",
//                             "C(=C)C>>C(O)C",
//                             "[C:2]=O>>[C:2]O",
//                             "C=O>>CO"};
    const char* strings[] = {"[C:1](=O)C>>[C:1](O)C",
                             "C=O>>CO",
                             "[C:2]=O>>[C:2]O",
                             "[C:1](=O)C>>C(O)[C:1]",
                             "[C:2]=O>>[C:2]O",
                             "[C:1](=O)C>>C(O)[C:1]",
                             };

    // Determine number of elements in the string array
    int size = sizeof(strings) / sizeof(char*);
    npy_intp dims[] = {size};

    // Create input and output numpy arrays
    PyArrayObject* np_input = createNumpyArrFromString(size, strings);
    PyArrayObject* np_output = (PyArrayObject*)PyArray_EMPTY(1, dims, NPY_BOOL, NPY_ARRAY_C_CONTIGUOUS);

    size = PyArray_SIZE(np_input);

    PyObject* repr = PyObject_Repr((PyObject*)np_input);
    printf("Input array:\n%s\n", PyUnicode_AsUTF8(repr));

    // Access the in_data pointer and shape of the array
    PyObject** in_data = (PyObject**) PyArray_DATA(np_input);
    npy_bool* out_data = (npy_bool*) PyArray_DATA(np_output);

    // Create Indigo query object
    const char* querySmarts = "[C:1]=[O]>>[C:1]-[OX2]";
    int query = indigoLoadReactionSmartsFromString(querySmarts);
    indigoOptimize(query, NULL);
    printf("Query: %s\n", querySmarts);


//    int num_threads = omp_get_max_threads();
    int num_threads = 2;

//    Py_BEGIN_ALLOW_THREADS;

    // Divide input in_data into batches
    int batch_size = size / num_threads;

//    #pragma omp parallel num_threads(num_threads)
//    {
//        int thread_num = omp_get_thread_num();
//        int start_idx = thread_num * batch_size;
//        int end_idx = start_idx + batch_size;
//        // last batch
//        if (thread_num == num_threads - 1) {
//            end_idx = size;
//        }
//        int current_size = end_idx - start_idx;
//        reactionSubstructureSearchBatch_py(in_data + start_idx,
//                                        out_data + start_idx,
//                                        current_size,
//                                        query,
//                                        strbuf,
//                                        "DAYLIGHT-AAM");
//    }

//    Py_END_ALLOW_THREADS;

    // Iterate over the in_data pointer to access the Python strings
//    reactionSubstructureSearchBatch_py(in_data, out_data, PyArray_SIZE(np_input), query, NULL);

    // Convert Python strings to C strings
    char ** smiles_array = numpyAsUTF8(np_input);
    reactionMatchBatch_py(smiles_array, out_data, size, query, "DAYLIGHT-AAM");
    free(smiles_array);

    PyArrayObject* np_output_vec = reactionMatchVec(np_input, "[C:1]=[O]>>[C:1]-[OX2]", "DAYLIGHT-AAM");
//    PyArrayObject* np_output_vec = NULL;


    // Print output array

    if (np_output != NULL) {
        repr = PyObject_Repr((PyObject*)np_output);
        printf("Output array:\n%s\n", PyUnicode_AsUTF8(repr));
        Py_DECREF(np_output);
    }

    if (np_output_vec != NULL) {
        repr = PyObject_Repr((PyObject*)np_output_vec);
        printf("Output array vec:\n%s\n", PyUnicode_AsUTF8(repr));
        Py_DECREF(np_output_vec);
    }


    // Decrement the reference count of the NumPy array object
    Py_DECREF(np_input);


    // Shut down Python interpreter
    Py_Finalize();

    // End Indigo session and free memory
    indigoFreeAllObjects();
    indigoReleaseSessionId(sessionId);

    const char* rxnSmarts = "[F:1][C:2]([F:31])([F:30])[C:3]1[CH:4]=[C:5]([C@H:13]2[O:17][C:16](=[O:18])[N:15]([CH2:19][C:20]3[C:25](Br)=[CH:24][N:23]=[C:22]([S:27][CH3:28])[N:21]=3)[C@H:14]2[CH3:29])[CH:6]=[C:7]([C:9]([F:12])([F:11])[F:10])[CH:8]=1.[CH:32]([C:35]1[CH:36]=[C:37](B(O)O)[C:38]([O:41][CH3:42])=[N:39][CH:40]=1)([CH3:34])[CH3:33].C([O-])([O-])=O.[K+].[K+].[NH4+].[Cl-]>C(OCC)(=O)C.[Pd](Cl)Cl.C(P(C(C)(C)C)[C-]1C=CC=C1)(C)(C)C.[C-]1(P(C(C)(C)C)C(C)(C)C)C=CC=C1.[Fe+2]>[F:1][C:2]([F:31])([F:30])[C:3]1[CH:4]=[C:5]([C@H:13]2[O:17][C:16](=[O:18])[N:15]([CH2:19][C:20]3[C:25]([C:37]4[C:38]([O:41][CH3:42])=[N:39][CH:40]=[C:35]([CH:32]([CH3:34])[CH3:33])[CH:36]=4)=[CH:24][N:23]=[C:22]([S:27][CH3:28])[N:21]=3)[C@H:14]2[CH3:29])[CH:6]=[C:7]([C:9]([F:12])([F:11])[F:10])[CH:8]=1";
    const char* querySmarts = "[B;X3,4]-[C,c:1].[C,c:2]-[Cl,Br,I,$([O]-S)]>>[C,c:1]-[C,c:2]";

    int rxn = indigoLoadReactionSmartsFromString(rxnSmarts);
    int query = indigoLoadReactionSmartsFromString(querySmarts);
    int matcher = indigoSubstructureMatcher(rxn, NULL);
    int match = indigoMatch(matcher, query);

    if (reactionSubstructureMatch((char*)rxnSmarts, query))
        printf("%s - SMARTS query matched!\n", querySmarts);
    else
        printf("No matches found\n");



    return 0;
}


*/

// Methods definition

/**
 * Set option by name
 * Option names and values always passed as Unicode Python strings
 * type checking is done on the python side
 * @param self
 * @param args option name, option value as python strings
 * @return None if successful
 */
PyObject* _set_option(PyObject* self, PyObject* args) {
    char* option_name;
    char* option_value;
    if (!PyArg_ParseTuple(args, "ss", &option_name, &option_value)) {
        return NULL;
    }

    if (strcmp(option_name, "num_cores") == 0) {
        options->num_cores = atoi(option_value);
//        printf("set num_cores: %d\n", options->num_cores);
        return Py_None;
    } else {
        printf("Option %s not allowed\n", option_name);
        return NULL;
    }
}

/**
 * Get option by name
 * @param self
 * @param args option name
 * @return option value, always python string
 */
PyObject* _get_option(PyObject* self, PyObject* args) {
    char* option_name;
    if (!PyArg_ParseTuple(args, "s", &option_name)) {
        return NULL;
    }

    if (strcmp(option_name, "num_cores") == 0) {
//        printf("get num_cores: %d\n", options->num_cores);
        PyObject* value = PyUnicode_FromFormat("%d", options->num_cores);
//        Py_DecRef(value);
        return value;
    } else {
        printf("Option %s not found\n");
        return NULL;
    }

}

PyObject* _rxn_match(PyObject* self, PyObject* args, PyObject* kwargs) {
    static char* keywords[] = {"np_input", "query_smarts", "aam_mode", "num_cores", NULL};

    PyArrayObject* np_input;
    char* querySmarts;
    char* aamMode;
    int numCores;

    // Parse the arguments using PyArg_ParseTuple
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Ossi", keywords, &np_input, &querySmarts, &aamMode, &numCores)) {
        return NULL;
    }

    return (PyObject*) reactionMatchNumPy(np_input, querySmarts, aamMode, numCores, options);
}



// Define the module methods
static PyMethodDef methods[] = {
        {"_rxn_match", (PyCFunction) _rxn_match, METH_VARARGS | METH_KEYWORDS, "C-API vecorized reaction match"},
        {"_set_option", (PyCFunction) _set_option, METH_VARARGS, "Set option"},
        {"_get_option", (PyCFunction) _get_option, METH_VARARGS, "Get option"},
        {NULL, NULL, 0, NULL}   // Sentinel value to indicate end of list
};



// Define the module structure
static PyModuleDef module_def = {
        PyModuleDef_HEAD_INIT,
        "Internal \"_chemivec\" module", // Also don't forget to change name here
        "Vectorized cheminformatics module, based on Indigo C-API",
        -1,
        methods,
        NULL, // Optional slot definitions
        NULL, // Optional traversal function
        NULL, // Optional clear function
        freeModule  // Optional module deallocation function
};


// Define module name here by PyInit_<your_modul_ename>
PyMODINIT_FUNC PyInit__chemivec(void) {
    import_array();
    qword SessionId = indigoAllocSessionId();
    initOptions(SessionId);
    PyObject* module = PyModule_Create(&module_def);
    return module;
}