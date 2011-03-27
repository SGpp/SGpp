/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Dirk Pflueger (Dirk.Pflueger@in.tum.de)

#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"
#include "exception/data_exception.hpp"

#include <sstream>
#include <cmath>
#include <algorithm>
#include <cstring>

#include <iostream>

//#include "common/AlignedMemory.hpp"

DataMatrix::DataMatrix(size_t nrows, size_t ncols) :
	nrows(nrows), ncols(ncols), unused(0), inc_rows(100) {
	// create new vector
	this->data = new double[nrows * ncols];
}

DataMatrix::DataMatrix(const DataMatrix &matr) :
	unused(0), inc_rows(100) {
	this->nrows = matr.nrows;
	this->ncols = matr.ncols;
	// create new vector
	this->data = new double[nrows * ncols];
	// copy data
	memcpy(this->data, matr.data, nrows * ncols * sizeof(double));
}

DataMatrix::DataMatrix(double* input, size_t nrows, size_t ncols) :
	nrows(nrows), ncols(ncols), unused(0), inc_rows(100) {
	// create new vector
	this->data = new double[nrows * ncols];
	// copy data
	memcpy(this->data, input, nrows * ncols * sizeof(double));
}

/**
 DataMatrix::DataMatrix(DataMatrixDefinition& DataMatrixDef) {
 this->nrows = DataMatrixDef.nrows;
 this->ncols = DataMatrixDef.ncols;
 this->unused = DataMatrixDef.unused;
 this->data = DataMatrixDef.pointerToData;
 }

 void DataMatrix::getDataMatrixDefinition(DataMatrixDefinition& DataMatrixDef) {
 DataMatrixDef.nrows = this->nrows;
 DataMatrixDef.ncols = this->ncols;
 DataMatrixDef.unused = this->unused;
 DataMatrixDef.pointerToData = this->data;
 }
 **/

void DataMatrix::resize(size_t nrows) {
	// don't do anyhing, if vector already has the correct size
	if (nrows == this->nrows) {
		return;
	}

	// create new vector
	double* newdata = new double[nrows * this->ncols];
	// copy entries of old vector
	memcpy(newdata, this->data, std::min(this->nrows, nrows) * this->ncols
			* sizeof(double));
	delete[] this->data;

	this->data = newdata;
	this->nrows = nrows;
}

void DataMatrix::resizeZero(size_t nrows) {
	// don't do anyhing, if vector already has the correct size
	if (nrows == this->nrows) {
		return;
	}

	// create new vector
	double* newdata = new double[nrows * this->ncols];
	// copy entries of old vector
	memcpy(newdata, this->data, std::min(this->nrows, nrows) * this->ncols
			* sizeof(double));

	// set new elements to zero
	for (size_t i = std::min(this->nrows, nrows) * this->ncols; i < nrows
			* this->ncols; i++) {
		newdata[i] = 0.0;
	}

	delete[] this->data;

	this->data = newdata;
	this->nrows = nrows;
}

void DataMatrix::addSize(size_t inc_rows) {
	// create new vector
	double* newdata = new double[(this->nrows + inc_rows) * this->ncols];
	// copy entries of old vector
	memcpy(newdata, this->data, this->nrows * this->ncols * sizeof(double));

	delete[] this->data;

	this->data = newdata;
	this->unused = inc_rows;
}

size_t DataMatrix::appendRow() {
	// enlarge, if necessary
	if (unused == 0) {
		addSize(this->inc_rows);
	}

	size_t x = nrows;
	nrows++;
	unused--;

	return x;
}

void DataMatrix::transpose()
{
	double* newData = new double[nrows * ncols];

	for (size_t i = 0; i < nrows; i++)
	{
#ifdef __ICC
		#pragma ivdep
#endif
		for (size_t j = 0; j < ncols; j++)
		{
			newData[(j*nrows)+i] = data[(i*ncols)+j];
		}
	}

	delete[] data;
	data = newData;
	size_t tmpRows = nrows;
	nrows = ncols;
	ncols = tmpRows;
}

size_t DataMatrix::appendRow(DataVector& vec) {
	if (vec.getSize() != this->ncols) {
		throw new sg::data_exception(
				"DataMatrix::appendRow : Dimensions do not match");
	}
	size_t x = appendRow();
	// copy data
	memcpy(&this->data[x * this->ncols], vec.getPointer(), this->ncols
			* sizeof(double));
	return x;
}

void DataMatrix::setAll(double value) {
	size_t n = nrows * ncols;
#ifdef __ICC
	#pragma ivdep
	#pragma vector aligned
#endif
	for (size_t i = 0; i < n; i++) {
		data[i] = value;
	}
}

void DataMatrix::getRow(size_t row, DataVector& vec) {
	if (vec.getSize() != this->ncols) {
		throw new sg::data_exception(
				"DataMatrix::getRow : Dimensions do not match");
	}
	for (size_t i = 0; i < this->ncols; i++) {
		vec[i] = this->data[row * ncols + i];
	}
}

void DataMatrix::getRow(size_t row, std::vector<double>& vec) {
    vec.clear();

    for (size_t i = 0; i < this->ncols; i++) {
        vec.push_back(data[row * ncols + i]);
    }
}

void DataMatrix::setRow(size_t row, DataVector& vec) {
	if (vec.getSize() != this->ncols) {
		throw new sg::data_exception(
				"DataMatrix::setRow : Dimensions do not match");
	}
#ifdef __ICC
	#pragma ivdep
#endif
	for (size_t i = 0; i < this->ncols; i++) {
		this->data[row * ncols + i] = vec[i];
	}
}

void DataMatrix::getColumn(size_t col, DataVector& vec) {
	if (vec.getSize() != this->nrows) {
		throw new sg::data_exception(
				"DataMatrix::getColumn : Dimensions do not match");
	}
#ifdef __ICC
	#pragma ivdep
#endif
	for (size_t j = 0; j < this->nrows; j++) {
		vec[j] = data[j * ncols + col];
	}
}

void DataMatrix::setColumn(size_t col, DataVector& vec) {
	if (vec.getSize() != this->nrows) {
		throw new sg::data_exception(
				"DataMatrix::setColumn : Dimensions do not match");
	}
#ifdef __ICC
	#pragma ivdep
#endif
	for (size_t j = 0; j < this->nrows; j++) {
		data[j * ncols + col] = vec[j];
	}
}

void DataMatrix::copyFrom(const DataMatrix& matr) {
	// don't copy from yourself
	if (this == &matr) {
		return;
	}
	/*
	 if (nrows != vec.nrows || ncols != vec.ncols) {
	 delete[] data;
	 nrows = vec.nrows;
	 ncols = vec.ncols;
	 this->data = new double[nrows * ncols];
	 }
	 */
	memcpy(this->data, matr.data, std::min(this->nrows*this->ncols, matr.nrows*matr.ncols) * sizeof(double));
}

/*
 void DataMatrix::copySmall(const DataMatrix& vec) {
 if (this == &vec) {
 return;
 }

 if (vec.ncols != 1 || ncols != 1 || nrows < vec.nrows) {
 return;
 }
 memcpy(this->data, vec.data, vec.nrows * sizeof(double));
 }

 DataMatrix& DataMatrix::operator=(const DataMatrix &vec) {
 if (this == &vec) {
 return *this;
 }

 if (nrows != vec.nrows || ncols != vec.ncols) {
 delete[] data;
 nrows = vec.nrows;
 ncols = vec.ncols;
 this->data = new double[nrows * ncols];
 }
 memcpy(this->data, vec.data, nrows * ncols * sizeof(double));
 return *this;
 }
 */
DataMatrix& DataMatrix::operator=(const DataMatrix &matr) {
    if (this == &matr) {
        return *this;
    }

    if (nrows*ncols != matr.ncols*matr.nrows) {
		throw new sg::data_exception(
				"DataMatrix::= : Dimensions do not match");
    }
    memcpy(this->data, matr.data, nrows*ncols * sizeof(double));
    return *this;
}


void DataMatrix::add(DataMatrix &matr) {
	if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
		throw new sg::data_exception(
				"DataMatrix::add : Dimensions do not match");
	}
	size_t n = nrows * ncols;
#ifdef __ICC
	#pragma ivdep
	#pragma vector aligned
#endif
	for (size_t i = 0; i < n; i++) {
		data[i] += matr.data[i];
	}
}

void DataMatrix::sub(DataMatrix& matr) {
	if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
		throw new sg::data_exception(
				"DataMatrix::sub : Dimensions do not match");
	}
	size_t n = nrows * ncols;

#ifdef __ICC
	#pragma ivdep
	#pragma vector aligned
#endif
	for (size_t i = 0; i < n; i++) {
		data[i] -= matr.data[i];
	}
}

void DataMatrix::componentwise_mult(DataMatrix& matr) {
	if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
		throw new sg::data_exception(
				"DataMatrix::componentwise_mult : Dimensions do not match");
	}
	size_t n = nrows * ncols;

#ifdef __ICC
	#pragma ivdep
	#pragma vector aligned
#endif
	for (size_t i = 0; i < n; i++) {
		data[i] *= matr.data[i];
	}
}

void DataMatrix::componentwise_div(DataMatrix& matr) {
	if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
		throw new sg::data_exception(
				"DataMatrix::componentwise_div : Dimensions do not match");
	}
	size_t n = nrows * ncols;

#ifdef __ICC
	#pragma ivdep
	#pragma vector aligned
#endif
	for (size_t i = 0; i < n; i++) {
		data[i] /= matr.data[i];
	}
}

/*
 void DataMatrix::getLine(int row, DataMatrix& vec) {
 for (int i = 0; i < this->ncols; i++) {
 vec[i] = data[row * ncols + i];
 }
 }

 void DataMatrix::getLine(int row, std::vector<double>& vec) {
 vec.clear();

 for (int i = 0; i < this->ncols; i++) {
 vec.push_back(data[row * ncols + i]);
 }
 }
 */

/*
 double DataMatrix::dotProduct(DataMatrix &vec) {
 double sum = 0.0;

 for (int i = 0; i < nrows; i++) {
 sum += data[i] * vec.data[i];
 }
 return sum;
 }
 */

void DataMatrix::mult(double scalar) {
	size_t n = nrows * ncols;

#ifdef __ICC
	#pragma ivdep
	#pragma vector aligned
#endif
	for (size_t i = 0; i < n; i++) {
		data[i] *= scalar;
	}
}

void DataMatrix::sqr() {
	size_t n = nrows * ncols;

#ifdef __ICC
	#pragma vector aligned
#endif
	for (size_t i = 0; i < n; i++) {
		data[i] = data[i] * data[i];
	}
}

void DataMatrix::sqrt() {
	size_t n = nrows * ncols;

#ifdef __ICC
	#pragma vector aligned
#endif
	for (size_t i = 0; i < n; i++) {
		data[i] = std::sqrt(data[i]);
	}
}

void DataMatrix::abs() {
	size_t n = nrows * ncols;

#ifdef __ICC
	#pragma vector aligned
#endif
	for (size_t i = 0; i < n; i++) {
		data[i] = std::abs(data[i]);
	}
}

double DataMatrix::sum() {
	size_t n = nrows * ncols;
	double result = 0.0;

#ifdef __ICC
	#pragma vector aligned
#endif
	for (size_t i = 0; i < n; i++) {
		result += data[i];
	}
	return result;
}
/*
 double DataMatrix::maxNorm() {
 int n = nrows * ncols;
 double max = 0.0;
 for (int i = 0; i < n; i++) {
 if (max < fabs(data[i]))
 {
 max = fabs(data[i]);
 }
 }
 return max;
 }

 void DataMatrix::partitionClasses(double border) {
 int n = nrows * ncols;
 for (int i = 0; i < n; i++) {
 data[i] = data[i] > border ? 1.0 : -1.0;
 }
 }

 void DataMatrix::axpy(double alpha, DataMatrix& x) {
 if (nrows != x.nrows || ncols != x.ncols) {
 return;
 }
 int n = nrows * ncols;
 double* p_x = x.data;
 double* p_d = data;

 for (int i = 0; i < n; i++) {
 p_d[i] += alpha * p_x[i];
 }
 }
 */

void DataMatrix::normalizeDimension(size_t d) {
	normalizeDimension(d, 0.0);
}

void DataMatrix::normalizeDimension(size_t d, double border) {
	size_t n = nrows * ncols;
	if (ncols <= d) {
		throw new sg::data_exception(
				"DataMatrix::normalizeDimension : Not enough columns in DataMatrix");
	}
	// determine min and max
	double xmin, xmax;
	minmax(d, &xmin, &xmax);

	double delta = (xmax - xmin)/(1-2*border);
	if (delta == 0.0) {
		for (size_t i = d; i < n; i += ncols) {
			data[i] = 0.5;
		}
	}
	else {
		for (size_t i = d; i < n; i += ncols) {
			data[i] = (data[i] - xmin) / delta + border;
		}
	}
}

void DataMatrix::toString(std::string& text) {
	std::stringstream str;
	str << "[";
	for (size_t i = 0; i < nrows; i++) {
		str << "[";
		for (size_t j = 0; j < ncols; j++) {
			if (j != 0) {
				str << ",";
			}
			str << " " << data[i*ncols+j];
		}
		if (i == nrows-1) {
			str << "]";
		}
		else {
			str << "]," << std::endl;
		}
	}
	str << "]";
	text = str.str();
}

std::string DataMatrix::toString() {
	std::string str;
	toString(str);
	return str;
}

double DataMatrix::min(size_t d) {
	size_t n = nrows * ncols;
	double min = data[d];
	for (size_t i = d; i < n; i += ncols) {
		if (min > data[i]) {
			min = data[i];
		}
	}
	return min;
}

double DataMatrix::min() {
	size_t n = nrows * ncols;
	double min = data[0];
	for (size_t i = 1; i < n; i += 1) {
		if (min > data[i]) {
			min = data[i];
		}
	}
	return min;
}

double DataMatrix::max(size_t d) {
	size_t n = nrows * ncols;
	double max = data[d];
	for (size_t i = d; i < n; i += ncols) {
		if (max < data[i]) {
			max = data[i];
		}
	}
	return max;
}

double DataMatrix::max() {
	size_t n = nrows * ncols;
	double max = data[0];
	for (size_t i = 1; i < n; i += 1) {
		if (max < data[i]) {
			max = data[i];
		}
	}
	return max;
}

void DataMatrix::minmax(size_t col, double* min, double* max) {
	size_t n = nrows * ncols;
	if (ncols <= col) {
		throw new sg::data_exception(
				"DataMatrix::minmax : Not enough entries in DataMatrix");
	}
	// find min and max of column col
	double min_t = data[col];
	double max_t = data[col];
	for (size_t i = col; i < n; i += ncols) {
		if (min_t > data[i]) {
			min_t = data[i];
		}
		if (max_t < data[i]) {
			max_t = data[i];
		}
	}
	(*min) = min_t;
	(*max) = max_t;
}

void DataMatrix::minmax(double* min, double* max) {
	size_t n = nrows * ncols;
	if (n == 0) {
		throw new sg::data_exception(
				"DataMatrix::minmax : Empty DataMatrix");
	}
	double min_t = data[0];
	double max_t = data[0];
	for (size_t i = 1; i < n; i += 1) {
		if (min_t > data[i]) {
			min_t = data[i];
		}
		if (max_t < data[i]) {
			max_t = data[i];
		}
	}
	(*min) = min_t;
	(*max) = max_t;
}

double* DataMatrix::getPointer() {
	return data;
}

DataMatrix::~DataMatrix() {
	delete[] data;
}

size_t DataMatrix::getNumberNonZero() {
	size_t n = nrows * ncols;
	size_t nonZero = 0;

	for (size_t i = 0; i < n; i++) {
		if (fabs(data[i]) > 0.0) {
			nonZero++;
		}
	}
	return nonZero;
}
