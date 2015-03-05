// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%include "base/src/sgpp/globaldef.hpp"

%inline %{
/*
* Memory managment function used in sg::base::DataVector::__array()
* it simply decrements the number of references to the PyObject datavector 
* every time a referencing ndarray is deleted.
* After reference counter reaches 0, the memory of DataVector will be deallocated 
* automatically. 
*/
void free_array(void* ptr, void* dv){
  // SGPP::float_t* vec = (SGPP::float_t *) ptr;
      PyObject* datavector = (PyObject*)dv;
      Py_DECREF(datavector);
    }
%}
namespace SGPP
{
namespace base
{

%apply SGPP::float_t *OUTPUT { SGPP::float_t* min, SGPP::float_t* max };
%apply std::string *OUTPUT { std::string& text };

%rename(__str__) DataVector::toString;

%rename(__getitem__) DataVector::get(size_t i) const;
%rename(__setitem__) DataVector::set(size_t i, SGPP::float_t value);
%rename(assign) DataVector::operator=;
%rename(__len__) DataVector::getSize;


    
class DataVector
{
public:

// typemap allowing to pass sequence of numbers to constructor
/*%typemap(in) (SGPP::float_t *input, int size)
{
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expected a sequence");
    return NULL;
  }
  $2 = PySequence_Size($input);
  $1 = (SGPP::float_t *) malloc(sizeof(SGPP::float_t)*$2);
  for (int i = 0; i < $2; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    if (PyNumber_Check(o)) {
      $1[i] = (SGPP::float_t) PyFloat_AsDouble(o);
    } else {
      PyErr_SetString(PyExc_ValueError,"Sequence elements must be numbers");      
      return NULL;
    }
  }
}
%typemap(freearg) (SGPP::float_t *input, int size)
{
  if ($1) free($1);
}
%typecheck(SWIG_TYPECHECK_FLOAT) (SGPP::float_t *input, int size)
{
$1 = PySequence_Check($input) ? 1 : 0;
}
*/
        // Constructors
  DataVector(size_t size);
  DataVector(DataVector& vec);
  DataVector(SGPP::float_t *input, size_t size);
  DataVector(std::vector<SGPP::float_t> input);
  DataVector(SGPP::base::DataVectorDefinition& DataVectorDef);

  void resize(size_t size);
  void resizeZero(size_t size);
  void addSize(size_t add);
  size_t append();
  size_t append(SGPP::float_t value);
  
  void setAll(SGPP::float_t value);
  
  void copyFrom(const DataVector& vec);
//  void copySmall(const DataVector& vec);
  DataVector& operator=(const DataVector& vec); 
  
  SGPP::float_t get(size_t i) const;
  void set(size_t i, SGPP::float_t value);

  void add(DataVector& vec);
  void sub(DataVector& vec);
  void componentwise_mult(DataVector& vec);
  void componentwise_div(DataVector& vec);
  void mult(SGPP::float_t scalar);
  void sqr();
  void sqrt();
  void abs();
  SGPP::float_t sum();
  SGPP::float_t min();
  SGPP::float_t max();
  void minmax(SGPP::float_t* min, SGPP::float_t* max);
  SGPP::float_t maxNorm();
  SGPP::float_t RMSNorm();
  SGPP::float_t l2Norm();
  SGPP::float_t dotProduct(DataVector& vec);
  
  void axpy(SGPP::float_t alpha, DataVector& x);
  
  size_t getSize();
  size_t getUnused();
  size_t getInc();
  void setInc(size_t inc_elems);
  size_t getNumberNonZero();  
    
  void partitionClasses(SGPP::float_t border);
  void normalize();
  void normalize(SGPP::float_t border);
  
  std::string toString();
  
  %extend {
    // Create a ndarray view from the DataVector data
    // an alternative approach using ARGOUTVIEW will fail since it does not allow to do a proper memory management
    PyObject* __array(PyObject* datavector){
        //Get the data and number of entries
      SGPP::float_t *vec = $self->getPointer();
      int n = $self->getSize();

      npy_intp dims[1] = {n};
      
      // Create a ndarray with data from vec
#if USE_DOUBLE_PRECISION == 1
      PyObject* arr = PyArray_SimpleNewFromData(1,dims, NPY_DOUBLE, vec);
#else
      PyObject* arr = PyArray_SimpleNewFromData(1,dims, NPY_FLOAT, vec);
#endif /* USE_DOUBLE_PRECISION */
      
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
  }
  
};



}
}




//%include "base/src/sgpp/globaldef.hpp"
//
//
//%inline %{
///*
//* Memory managment function used in SGPP::base::DataVector::__array()
//* it simply decrements the number of references to the PyObject datavector 
//* every time a referencing ndarray is deleted.
//* After reference counter reaches 0, the memory of DataVector will be deallocated 
//* automatically. 
//*/
//void free_array(void* ptr, void* dv){
//  // SGPP::float_t* vec = (SGPP::float_t *) ptr;
//      PyObject* datavector = (PyObject*)dv;
//      Py_DECREF(datavector);
//    }
//%}
//
//namespace SGPP
//{
//namespace base
//{
//
//%apply SGPP::float_t *OUTPUT { SGPP::float_t* min, SGPP::float_t* max };
//%apply std::string *OUTPUT { std::string& text };
//%rename(__str__) DataVector::toString;
//%rename(__getitem__) DataVector::get(size_t i) const;
//%rename(__setitem__) DataVector::set(size_t i, SGPP::float_t value);
//%rename(assign) DataVector::operator=;
//%rename(__len__) DataVector::getSize;
//%ignore DataVector::operator[];
//%ignore DataVector::DataVector(std::vector<int> input);
//%ignore DataVector::toString(std::string &);
//
//
//%extend DataVector{
//    // Create a ndarray view from the DataVector data
//    // an alternative approach using ARGOUTVIEW will fail since it does not allow to do a proper memory management
//    PyObject* __array(PyObject* datavector){
//        //Get the data and number of entries
//      SGPP::float_t *vec = $self->getPointer();
//      int n = $self->getSize();
//
//      npy_intp dims[1] = {n};
//      
//      // Create a ndarray with data from vec
//      PyObject* arr = PyArray_SimpleNewFromData(1,dims, NPY_DOUBLE, vec);
//      
//      // Let the array base point on the original data, free_array is a additional destructor for our ndarray, 
//      // since we need to DECREF DataVector object
//      PyObject* base = PyCObject_FromVoidPtrAndDesc((void*)vec, (void*)datavector, free_array);
//      PyArray_BASE(arr) = base;
//      
//      // Increase the number of references to PyObject DataVector, after the object the variable is reinitialized or deleted the object
//      // will still be on the heap, if the reference counter is positive.
//      Py_INCREF(datavector);
//      
//      return arr;
//    }
//     %pythoncode
//     {
//        def array(self):   
//          return self.__array(self)
//     }
//  };
//  
//
//}
//}
//
//%include "base/src/sgpp/base/datatypes/DataVector.hpp"