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

int CombiConfigurator::getStandardCombi(vector<combiConfig> &vec, size_t dim, size_t level) {
  cout << "getStandardCombi";
  PyObject *pName, *pModule, *pFunc;
  PyObject *pArgs, *pValue;

  Py_Initialize();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.append(\"../spatially-adaptive-combi\")");

  pName = PyUnicode_DecodeFSDefault("SGDEAdapter");
  pModule = PyImport_Import(pName);
  Py_DECREF(pName);

  pFunc = PyObject_GetAttrString(pModule, "getstandardcombi");
  pArgs = PyTuple_New(2);

  pValue = PyLong_FromLong(dim);
  PyTuple_SetItem(pArgs, 0, pValue);

  pValue = PyLong_FromLong(level);
  PyTuple_SetItem(pArgs, 1, pValue);

  pValue = PyObject_CallObject(pFunc, pArgs);

  for (int j = 0; j < PyList_Size(pValue); j++) {
    combiConfig pair;
    pair.coef = 0.0;
    pair.levels = std::vector<size_t>();
    vec.push_back(pair);
    cout << "PYBREAK 1";
    vec.at(j).coef = PyFloat_AsDouble(PyList_GetItem((PyList_GetItem(pValue, j)), 0));
    cout << "PYBREAK 2";
    for (int c = 1; c < PyList_Size(PyList_GetItem(pValue, j)); c++) {
      cout << "PYBREAK 3";
      vec.at(j).levels.push_back(PyLong_AsSize_t(PyList_GetItem(PyList_GetItem(pValue, j), c)));
    }
    cout << "PYBREAK 4";
  }
  cout << "PYBREAK 5";

  Py_DECREF(pModule);
  Py_DECREF(pArgs);
  Py_DECREF(pValue);
  Py_XDECREF(pFunc);
  cout << "calling Py_Finalize: ";
  Py_FinalizeEx();
  cout << "Py_Finalize done";
  return 0;

  // throw base::tool_exception("To make this work, compile with USE_PYTHON_EMBEDDING=1");
}

} /* namespace datadriven */
} /* namespace sgpp */
