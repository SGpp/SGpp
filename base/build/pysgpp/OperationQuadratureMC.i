// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%{

namespace sgpp
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
   
   // ensure thread state
   PyGILState_STATE d_gstate;
   d_gstate = PyGILState_Ensure();

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
   Py_DECREF(lst);
   
   // call Python
   result = PyObject_CallObject(func,arglist);
   // trash arglist and lst
   Py_DECREF(arglist);
  
   if (result) {
     dres = PyFloat_AsDouble(result);
   }
   Py_XDECREF(result);
   
   // release thread state
   PyGILState_Release(d_gstate);
   
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
   double doQuadratureFunc(PyObject *pyfunc) {
     double d;
     d = self->doQuadratureFunc(sgpp::base::PythonCallBackFunc, (void *) pyfunc);
     return d;
   }
   double doQuadratureL2Error(PyObject *pyfunc, sgpp::base::DataVector& alpha) {
     double d;
     d = self->doQuadratureL2Error(sgpp::base::PythonCallBackFunc, (void *) pyfunc, alpha);
     return d;
   }
}
%enddef

%include "base/src/sgpp/base/tools/OperationQuadratureMC.hpp"
QUADRATURE_CALLBACK_EXTEND(sgpp::base::OperationQuadratureMC)


