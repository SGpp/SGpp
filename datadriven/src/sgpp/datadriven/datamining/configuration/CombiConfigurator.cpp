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

#include <iostream>
#include <sgpp/base/exception/tool_exception.hpp>
#include <sgpp/datadriven/datamining/configuration/CombiConfigurator.hpp>
#include <vector>
using std::cout;

namespace sgpp {
namespace datadriven {

CombiConfigurator::CombiConfigurator() { Initialize(); }

void CombiConfigurator::Initialize() {
  if (initialized) {
    return;
  }
  initialized = true;
  Py_Initialize();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.append(\"../spatially-adaptive-combi\")");

  pName = PyUnicode_DecodeFSDefault("SGDEAdapter");
  pModule = PyImport_Import(pName);

  pFunc = PyObject_GetAttrString(pModule, "getstandardcombi");
  return;
}

void CombiConfigurator::Finalize() {
  if (!initialized) {
    return;
  }
  initialized = false;
  Py_DECREF(pName);
  Py_DECREF(pModule);
  Py_XDECREF(pFunc);
  cout << "Calling Py_Finalize: \n";
  Py_FinalizeEx();
  cout << "Py_Finalize done \n";
  return;
}

void CombiConfigurator::getStandardCombi(vector<combiConfig> &vec, size_t dim, size_t level) {
  PyObject *pArgs, *pValue;
  pArgs = PyTuple_New(2);

  pValue = PyLong_FromLong(dim);
  PyTuple_SetItem(pArgs, 0, pValue);
  pValue = PyLong_FromLong(level);
  PyTuple_SetItem(pArgs, 1, pValue);

  pValue = PyObject_CallObject(pFunc, pArgs);

  for (int j = 0; j < PyList_Size(pValue); j++) {
    combiConfig pair;
    pair.levels = std::vector<size_t>();
    pair.coef = PyFloat_AsDouble(PyList_GetItem((PyList_GetItem(pValue, j)), 0));
    vec.push_back(pair);
    for (int c = 1; c < PyList_Size(PyList_GetItem(pValue, j)); c++) {
      // PyLong_AsSize_t is returning garbage, thats why PyLong_AsLong is used
      vec.at(j).levels.push_back(PyLong_AsLong(PyList_GetItem(PyList_GetItem(pValue, j), c)));
    }
  }

  Py_DECREF(pArgs);
  Py_DECREF(pValue);
  return;

  // throw base::tool_exception("To make this work, compile with USE_PYTHON_EMBEDDING=1");
}

} /* namespace datadriven */
} /* namespace sgpp */
