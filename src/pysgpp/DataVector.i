/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (Dirk.Pflueger@in.tum.de)

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
%typemap(in) (double *input, int size)
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
        // Construcors
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

};


