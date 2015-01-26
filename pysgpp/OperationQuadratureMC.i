/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

%{

namespace sg
{
namespace base
{

/* This function has to match the type and parameters of the
 * C++ callback functions that are used. However, the
 * clientdata pointer is used for holding a reference to a 
 * Python callable object. 
 */
  static double PythonCallBackFunc(int len, double* a, void *clientdata)
{
  PyObject *func, *lst, *arglist;
   PyObject *result;
   double    dres = 0;

   // get Python function
   func = (PyObject *) clientdata;
   // build argument list (only Python list, convert double* to Python list)
   lst = PyTuple_New(len);        // alternatively: PyList_New(len)
   if (!lst) {
     PyErr_SetString(PyExc_TypeError, "No data provided!");
     return NULL;
   }
   for (int i=0; i<len; i++) {
     // create new Python double
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

%include "src/sgpp/base/tools/OperationQuadratureMC.hpp"

// attach a new method to OperationQuadratureMC
%define QUADRATURE_CALLBACK_EXTEND(class)
%extend class {
   // set a Python function object as a callback function
   // overloads functions
   // note : PyObject *pyfunc is remapped with a typemap
   double doQuadratureFunc(PyObject *pyfunc) {
     double d;
     d = self->doQuadratureFunc(sg::base::PythonCallBackFunc, (void *) pyfunc);
     return d;
   }
   double doQuadratureL2Error(PyObject *pyfunc, sg::base::DataVector& alpha) {
     double d;
     d = self->doQuadratureL2Error(sg::base::PythonCallBackFunc, (void *) pyfunc, alpha);
     return d;
   }
}
%enddef

%include "src/sgpp/base/tools/OperationQuadratureMC.hpp"
QUADRATURE_CALLBACK_EXTEND(sg::base::OperationQuadratureMC)

#ifdef SG_MCM
%apply (long long int* IN_ARRAY1, int DIM1) {(long long int* n, int dim)};
%include "src/sgpp/mcm/tools/OperationQuadratureMCAdvanced.hpp"
QUADRATURE_CALLBACK_EXTEND(sg::mcm::OperationQuadratureMCAdvanced)
#endif

