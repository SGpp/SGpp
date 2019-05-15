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
using sgpp::base::tool_exception;

namespace sgpp {
namespace datadriven {
#ifdef USE_PYTHON_EMBEDDING
namespace spaceconfig {

/**
 * Flag that indicates if the python interpreter was initialized
 */
static bool initialized = false;

/**
 * Stores the path to the SGPP SpACE Adapter
 */
static PyObject *pModule;
}
#endif

CombiConfigurator::CombiConfigurator() { initializePython(); }

void CombiConfigurator::initAdaptiveScheme(size_t dim, size_t level) {
#ifdef USE_PYTHON_EMBEDDING
  PyObject *pFunc = PyObject_GetAttrString(spaceconfig::pModule, "initadaptivescheme");
  PyObject *pArgs, *pValue;
  pArgs = PyTuple_New(2);
  pValue = PyLong_FromLong(dim);
  PyTuple_SetItem(pArgs, 0, pValue);
  pValue = PyLong_FromLong(level);
  PyTuple_SetItem(pArgs, 1, pValue);

  combischeme = PyObject_CallObject(pFunc, pArgs);

  Py_DECREF(pFunc);
  Py_DECREF(pArgs);
  Py_DECREF(pValue);
#endif
#ifndef USE_PYTHON_EMBEDDING
  throw tool_exception("Compiled without USE_PYTHON_EMBEDDING");
#endif
  return;
}

void CombiConfigurator::getCombiScheme(vector<combiConfig> &vec) {
#ifdef USE_PYTHON_EMBEDDING
  vec.clear();
  PyObject *pFunc = PyObject_GetAttrString(spaceconfig::pModule, "getcombischeme");
  PyObject *pArgs, *pValue;
  pArgs = PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, combischeme);

  pValue = PyObject_CallObject(pFunc, pArgs);

  for (int j = 0; j < PyList_Size(pValue); j++) {
    combiConfig pair;
    pair.levels = std::vector<size_t>();
    pair.coef =
        static_cast<ssize_t>(PyFloat_AsDouble(PyList_GetItem((PyList_GetItem(pValue, j)), 0)));
    vec.push_back(pair);
    for (int c = 1; c < PyList_Size(PyList_GetItem(pValue, j)); c++) {
      // PyLong_AsSize_t is returning garbage, thats why PyLong_AsLong is used
      vec.at(j).levels.push_back(PyLong_AsLong(PyList_GetItem(PyList_GetItem(pValue, j), c)));
    }
  }

  Py_DECREF(pFunc);
  /* will this lead to problems? if i deref pArgs i destroy my combischeme*/
  // Py_DECREF(pArgs);
  Py_DECREF(pValue);
#endif
#ifndef USE_PYTHON_EMBEDDING
  throw tool_exception("Compiled without USE_PYTHON_EMBEDDING");
#endif
  return;
}

bool CombiConfigurator::isRefinable(combiConfig levelvec) {
  bool a = false;
#ifdef USE_PYTHON_EMBEDDING
  PyObject *pFunc = PyObject_GetAttrString(spaceconfig::pModule, "isrefinable");
  PyObject *pArgs, *pValue;
  pArgs = PyTuple_New(2);
  PyTuple_SetItem(pArgs, 0, combischeme);
  PyTuple_SetItem(pArgs, 1, combiConfAsPyObj(levelvec));

  pValue = PyObject_CallObject(pFunc, pArgs);

  a = PyObject_IsTrue(pValue);
  Py_DECREF(pFunc);
  // Py_DECREF(pArgs);
  Py_DECREF(pValue);
#endif
#ifndef USE_PYTHON_EMBEDDING
  throw tool_exception("Compiled without USE_PYTHON_EMBEDDING");
#endif
  return a;
}

void CombiConfigurator::refineComponent(combiConfig levelvec) {
#ifdef USE_PYTHON_EMBEDDING
  PyObject *pFunc = PyObject_GetAttrString(spaceconfig::pModule, "refineblock");
  PyObject *pArgs;
  pArgs = PyTuple_New(2);
  PyTuple_SetItem(pArgs, 0, combischeme);
  PyTuple_SetItem(pArgs, 1, combiConfAsPyObj(levelvec));

  combischeme = PyObject_CallObject(pFunc, pArgs);

  Py_DECREF(pFunc);
#endif
#ifndef USE_PYTHON_EMBEDDING
  throw tool_exception("Compiled without USE_PYTHON_EMBEDDING");
#endif
  return;
}

void CombiConfigurator::initializePython() {
#ifdef USE_PYTHON_EMBEDDING
  if (spaceconfig::initialized) {
    return;
  }
  spaceconfig::initialized = true;
  Py_Initialize();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.append(\"../spatially-adaptive-combi\")");

  PyObject *pName;
  pName = PyUnicode_DecodeFSDefault("SGDEAdapter");
  spaceconfig::pModule = PyImport_Import(pName);
  Py_DECREF(pName);
#endif
#ifndef USE_PYTHON_EMBEDDING
  throw tool_exception("Compiled without USE_PYTHON_EMBEDDING");
#endif
  return;
}

void CombiConfigurator::finalizePython() {
#ifdef USE_PYTHON_EMBEDDING
  if (!spaceconfig::initialized) {
    return;
  }
  spaceconfig::initialized = false;
  Py_DECREF(spaceconfig::pModule);
  cout << "Calling Py_Finalize: \n";
  Py_FinalizeEx();
  cout << "Py_Finalize done \n";
#endif
#ifndef USE_PYTHON_EMBEDDING
  throw tool_exception("Compiled without USE_PYTHON_EMBEDDING");
#endif
  return;
}

#ifdef USE_PYTHON_EMBEDDING
inline combiConfig CombiConfigurator::combiConfFromPyObj(PyObject *pValue) {
  combiConfig pair;
  pair.levels = std::vector<size_t>();
  pair.coef = static_cast<ssize_t>(PyFloat_AsDouble(PyList_GetItem(pValue, 0)));
  for (int c = 1; c < PyList_Size(pValue); c++) {
    // PyLong_AsSize_t is returning garbage, thats why PyLong_AsLong is used
    pair.levels.push_back(PyLong_AsLong(PyList_GetItem(pValue, c)));
  }
  return pair;
}
#endif

#ifdef USE_PYTHON_EMBEDDING
inline PyObject *CombiConfigurator::combiConfAsPyObj(combiConfig pair) {
  PyObject *pValue = PyList_New(pair.levels.size() + 1);
  PyList_SetItem(pValue, 0, PyLong_FromLong(pair.coef));
  for (size_t i = 0; i < pair.levels.size(); i++) {
    PyList_SetItem(pValue, i + 1, PyLong_FromLong(pair.levels.at(i)));
  }
  return pValue;
}
#endif

} /* namespace datadriven */
} /* namespace sgpp */
