%inline %{
/*
* Memory managment function used in SGPP::base::DataVector::__array()
* it simply decrements the number of references to the PyObject datavector 
* every time a referencing ndarray is deleted.
* After reference counter reaches 0, the memory of DataVector will be deallocated 
* automatically. 
*/
void free_array(void* ptr, void* dv){
  // double* vec = (double *) ptr;
      PyObject* datavector = (PyObject*)dv;
      Py_DECREF(datavector);
    }
%}

namespace SGPP
{
namespace base
{

%apply double *OUTPUT { double* min, double* max };
%apply std::string *OUTPUT { std::string& text };
%rename(__str__) DataVector::toString;
%rename(__getitem__) DataVector::get(size_t i) const;
%rename(__setitem__) DataVector::set(size_t i, double value);
%rename(assign) DataVector::operator=;
%rename(__len__) DataVector::getSize;
%ignore DataVector::operator[];
%ignore DataVector::DataVector(std::vector<int> input);
%ignore DataVector::toString(std::string &);


  %extend DataVector{
    // Create a ndarray view from the DataVector data
    // an alternative approach using ARGOUTVIEW will fail since it does not allow to do a proper memory management
    PyObject* __array(PyObject* datavector){
        //Get the data and number of entries
      double *vec = $self->getPointer();
      int n = $self->getSize();

      npy_intp dims[1] = {n};
      
      // Create a ndarray with data from vec
      PyObject* arr = PyArray_SimpleNewFromData(1,dims, NPY_DOUBLE, vec);
      
      // Let the array base point on the original data, free_array is a additional destructor for our ndarray, 
      // since we need to DECREF DataVector object
      PyObject* base = PyCObject_FromVoidPtrAndDesc((void*)vec, (void*)datavector, free_array);
      PyArray_BASE(arr) = base;
      
      // Increase the number of references to PyObject DataVector, after the object the variable is reinitialized or deleted the object
      // will still be on the heap, if the reference counter is positive.
      Py_INCREF(datavector);
      
      return arr;
    }
     %pythoncode
     {
        def array(self):   
          return self.__array(self)
     }
  };
  

}
}
%include "base/src/sgpp/base/datatypes/DataVector.hpp"