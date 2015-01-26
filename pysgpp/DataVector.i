/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (Dirk.Pflueger@in.tum.de)
%inline %{
/*
* Memory managment function used in sg::base::DataVector::__array()
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

		
class DataVector
{
public:

// typemap allowing to pass sequence of numbers to constructor
/*%typemap(in) (double *input, int size)
{
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expected a sequence");
    return NULL;
  }
  $2 = PySequence_Size($input);
  $1 = (double *) malloc(sizeof(double)*$2);
  for (int i = 0; i < $2; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    if (PyNumber_Check(o)) {
      $1[i] = (double) PyFloat_AsDouble(o);
    } else {
      PyErr_SetString(PyExc_ValueError,"Sequence elements must be numbers");      
      return NULL;
    }
  }
}
%typemap(freearg) (double *input, int size)
{
  if ($1) free($1);
}
%typecheck(SWIG_TYPECHECK_FLOAT) (double *input, int size)
{
$1 = PySequence_Check($input) ? 1 : 0;
}
*/
        // Constructors
	DataVector(size_t size);
	DataVector(DataVector& vec);
	DataVector(double *input, int size);

	void resize(size_t size);
	void resizeZero(size_t size);
	void addSize(size_t add);
	size_t append();
	size_t append(double value);
	
	void setAll(double value);
	
	void copyFrom(const DataVector& vec);
//	void copySmall(const DataVector& vec);
	DataVector& operator=(const DataVector& vec);	
	
	double get(size_t i) const;
	void set(size_t i, double value);

	void add(DataVector& vec);
	void sub(DataVector& vec);
	void componentwise_mult(DataVector& vec);
	void componentwise_div(DataVector& vec);
	void mult(double scalar);
	void sqr();
	void sqrt();
	void abs();
	double sum();
	double min();
	double max();
	void minmax(double* min, double* max);
	double maxNorm();
	double RMSNorm();
	double l2Norm();
	double dotProduct(DataVector& vec);
	
	void axpy(double alpha, DataVector& x);
	
	size_t getSize();
	size_t getUnused();
	size_t getInc();
	void setInc(size_t inc_elems);
	size_t getNumberNonZero();	
		
	void partitionClasses(double border);
	void normalize();
	void normalize(double border);
	
	std::string toString();
	
	%extend {
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
	}
	
};



}
}
