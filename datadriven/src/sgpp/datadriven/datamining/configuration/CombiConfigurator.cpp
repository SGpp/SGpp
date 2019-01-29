/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * CombiConfigurator.cpp
 *
 *  Created on: Jan 29, 2019
 *      Author: nico
 */
#include <sgpp/base/exception/tool_exception.hpp>
#include <sgpp/datadriven/datamining/configuration/CombiConfigurator.hpp>
#include <vector>

namespace sgpp {
namespace datadriven {

const int CombiConfigurator::getStandardCombi(vector<combiConfig> &vec, size_t d, size_t l) {
  PyObject *pName, *pModule, *pFunc;
  PyObject *pArgs, *pValue;
  vec.clear();
  long int dim = static_cast<long int>(d);
  long int level = static_cast<long int>(l);
  Py_Initialize();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.append(\"../spatially-adaptive-combi\")");
  pName = PyUnicode_DecodeFSDefault("SGDEAdapter");
  /* Error checking of pName left out */

  pModule = PyImport_Import(pName);
  Py_DECREF(pName);

  if (pModule != NULL) {
    pFunc = PyObject_GetAttrString(pModule, "getstandardcombi");
    /* pFunc is a new reference */

    if (pFunc && PyCallable_Check(pFunc)) {
      pArgs = PyTuple_New(2);
      pValue = PyLong_FromLong(dim);
      PyTuple_SetItem(pArgs, 0, pValue);
      pValue = PyLong_FromLong(level);
      PyTuple_SetItem(pArgs, 1, pValue);

      pValue = PyObject_CallObject(pFunc, pArgs);
      Py_DECREF(pArgs);
      if (pValue != NULL) {
        for (int j = 0; j < PyList_Size(pValue); j++) {
          combiConfig pair;
          pair.coef = 0.0;
          pair.levels = std::vector<size_t>();
          vec.push_back(pair);
          vec[j].coef = PyFloat_AsDouble(PyList_GetItem((PyList_GetItem(pValue, j)), 0));
          for (int c = 1; c < PyList_Size(PyList_GetItem(pValue, j)); c++) {
            vec[j].levels.push_back(
                static_cast<size_t>(PyLong_AsLong(PyList_GetItem(PyList_GetItem(pValue, j), c))));
          }
        }
        std::cout << "Results of the call" << std::endl;
        Py_DECREF(pValue);
      } else {
        Py_DECREF(pFunc);
        Py_DECREF(pModule);
        PyErr_Print();
        fprintf(stderr, "Call failed\n");
        return 1;
      }
    } else {
      if (PyErr_Occurred()) PyErr_Print();
      fprintf(stderr, "Cannot find function getstandardcombi(dim, level) in SGDEAdapter.py");
    }
    Py_XDECREF(pFunc);
    Py_DECREF(pModule);
  } else {
    PyErr_Print();
    fprintf(stderr, "Failed to load SGDEAdapter.py");
    return 1;
  }
  if (Py_FinalizeEx() < 0) {
    return 120;
  }
  return 0;
  throw base::tool_exception("To make this work, compile with USE_PYTHON_EMBEDDING=1");
}

} /* namespace datadriven */
} /* namespace sgpp */
