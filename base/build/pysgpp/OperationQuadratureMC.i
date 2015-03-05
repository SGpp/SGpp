// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%{

namespace SGPP
{
namespace base
{

/* This function has to match the type and parameters of the
 * C++ callback functions that are used. However, the
 * clientdata pointer is used for holding a reference to a 
 * Python callable object. 
 */
  static float_t PythonCallBackFunc(int len, float_t* a, void *clientdata)
{
  PyObject *func, *lst, *arglist;
   PyObject *result;
   float_t    dres = 0;

   // get Python function
   func = (PyObject *) clientdata;
   // build argument list (only Python list, convert float_t* to Python list)
   lst = PyTuple_New(len);        // alternatively: PyList_New(len)
   if (!lst) {
     PyErr_SetString(PyExc_TypeError, "No data provided!");
     return NULL;
   }
   for (int i=0; i<len; i++) {
     // create new Python float_t
     PyObject *num = PyFloat_FromDouble(a[i]);
     if (!num) {
       PyErr_SetString(PyExc_TypeError, "No data in list!");
       Py_DECREF(lst);
       return NULL;
     }
     // steals reference to num:
     PyTuple_SetItem(lst, i, num); // alternatively: PyList_SET_ITEM()
   }
   // build list of one Python object
   arglist = Py_BuildValue("(O)", lst);
   // call Python
   result = PyEval_CallObject(func,arglist);
   // trash arglist and lst
   Py_DECREF(arglist);
   Py_DECREF(lst);
   if (result) {
     dres = PyFloat_AsDouble(result);
   }
   Py_XDECREF(result);
   return dres;
}
}}
%}

// grab a python function object as a Python object
%typemap(in) PyObject *pyfunc {
  if (!PyCallable_Check($input)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      return NULL;
  }
  $1 = $input;
}

%include "base/src/sgpp/base/tools/OperationQuadratureMC.hpp"

// attach a new method to OperationQuadratureMC
%define QUADRATURE_CALLBACK_EXTEND(class)
%extend class {
   // set a Python function object as a callback function
   // overloads functions
   // note : PyObject *pyfunc is remapped with a typemap
   float_t doQuadratureFunc(PyObject *pyfunc) {
     float_t d;
     d = self->doQuadratureFunc(SGPP::base::PythonCallBackFunc, (void *) pyfunc);
     return d;
   }
   float_t doQuadratureL2Error(PyObject *pyfunc, SGPP::base::DataVector& alpha) {
     float_t d;
     d = self->doQuadratureL2Error(SGPP::base::PythonCallBackFunc, (void *) pyfunc, alpha);
     return d;
   }
}
%enddef

%include "base/src/sgpp/base/tools/OperationQuadratureMC.hpp"
QUADRATURE_CALLBACK_EXTEND(SGPP::base::OperationQuadratureMC)


