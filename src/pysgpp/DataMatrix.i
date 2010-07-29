/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (Dirk.Pflueger@in.tum.de)

%apply double *OUTPUT { double* min, double* max };
%apply std::string *OUTPUT { std::string& text };
%rename(__str__) DataMatrix::toString;

//%rename(__getitem__) DataMatrix::get(size_t row, size_t col) const;
//%rename(__setitem__) DataMatrix::set(int i, double value);
//%rename(assign) DataMatrix::operator=;
//%rename(__len__) DataMatrix::getSize;

class DataMatrix
{
public:
    // Constructors
    DataMatrix(size_t nrows, size_t ncols);
    DataMatrix(const DataMatrix& matr);
    //@todo (pflueged) Convert python list to double*
    DataMatrix(double* INPUT, size_t nrows, size_t ncols);

    void resize(size_t size);
    void resizeZero(size_t nrows);

    void addSize(size_t inc_nrows);
    size_t appendRow();
    void setAll(double value);
    void copyFrom(const DataMatrix& matr);

    double get(size_t row, size_t col) const;
    void set(size_t row, size_t col, double value);

//    %extend {
//    	    double getxy(int* INPUT) {
//	        	return INPUT[0]*INPUT[1];
//	    };
//    }	    

    void getRow(size_t row, DataVector& vec);
    void setRow(size_t row, DataVector& vec);
    void getColumn(size_t col, DataVector& vec);
    void setColumn(size_t col, DataVector& vec);
    double* getPointer();

    void add(DataMatrix& matr);
    void sub(DataMatrix& matr);
    void componentwise_mult(DataMatrix& matr);
    void componentwise_div(DataMatrix& matr);
    void mult(double scalar);
    void sqr();
    void sqrt();
    void abs();
    double sum();

    size_t getSize();
    size_t getUnused();
    size_t getNumberNonZero();
    size_t getNrows();
    size_t getNcols();
    size_t getInc();
    void setInc(size_t inc_rows);

    void normalizeDimension(size_t d);
    void normalizeDimension(size_t d, double border);

    double min(size_t col);
    double min();
    double max(size_t col);
    double max();
    void minmax(size_t col, double* min, double* max);
    void minmax(double* min, double* max);

//    void toString(std::string& text); // otherwise overloaded duplicate
    std::string toString();

};

