/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

%apply double *OUTPUT { double* min, double* max };
%apply std::string *OUTPUT { std::string& text };
/* %rename(__str__) DataVector::toString; */

/* %rename(__getitem__) DataVector::get(int i) const; */
/* %rename(__setitem__) DataVector::set(int i, double value); */
/* %rename(assign) DataVector::operator=; */
/* %rename(__len__) DataVector::getTotalSize; */

class DataVector
{
public:
	DataVector(int size);
	DataVector(int size, int dim);
	DataVector(DataVector& vec);
	DataVector(double* input, int size, int dim);
	
	void resize(int size);
	void addSize(int add);
	int addValue();
	
	void setAll(double value);
	
	void copyFrom(const DataVector& vec);
	void copySmall(const DataVector& vec);
	DataVector& operator=(const DataVector& vec);	
	
	double get(int i) const;
//	double get(int row, int col) const;
	void set(int i, double value);
//	void set(int row, int col, double value);

	void add(DataVector& vec);
	void sub(DataVector& vec);
	void componentwise_mult(DataVector& vec);
	void componentwise_div(DataVector& vec);
	void mult(double scalar);
	
	void sqr();
	void abs();
	double sum();
	
	void axpy(double alpha, DataVector& x);
	
	void getRow(int row, DataVector& vec);
	void setRow(int row, DataVector& vec);
	void getColumn(int col, DataVector& vec);
	void setColumn(int col, DataVector& vec);
	
	double dotProduct(DataVector& vec);
	
	int getSize();
	int getDim();
	int getTotalSize();	
	inline int getUnused();
		
	void partitionClasses(double border);
	void normalizeDimension(int d);
	void normalizeDimension(int d, double border);
	
	double min(int d);
	double min();
	double max(int d);
	double max();
	void minmax(int d, double* min, double* max);
	
	void toString(std::string& text);

};
