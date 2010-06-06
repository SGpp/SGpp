/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (Dirk.Pflueger@in.tum.de)

%apply double *OUTPUT { double* min, double* max };
%apply std::string *OUTPUT { std::string& text };
%rename(__str__) DataVector::toString;

%rename(__getitem__) DataVector::get(int i) const;
%rename(__setitem__) DataVector::set(int i, double value);
%rename(assign) DataVector::operator=;
%rename(__len__) DataVector::getTotalSize;

class DataVector
{
public:
	DataVector(size_t size);
//	DataVector(size_t size, size_t dim);
	DataVector(DataVector& vec);
	DataVector(double* input, size_t size);
	
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
//	double get(int row, int col) const;
	void set(size_t i, double value);
//	void set(int row, int col, double value);

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
	double dotProduct(DataVector& vec);
	
	void axpy(double alpha, DataVector& x);
	
//	void getRow(size_t row, DataVector& vec);
//	void setRow(size_t row, DataVector& vec);
//	void getColumn(size_t col, DataVector& vec);
//	void setColumn(size_t col, DataVector& vec);
	
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
