/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2007 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/


#include "data/DataVector.h"
#include <string.h>

#include <sstream>

#ifdef USEOMP
#include <omp.h>
#endif

DataVector::DataVector()
{

}

/**
 * Default Constructor.
 * @param size the size
 */
DataVector::DataVector(int size) : size(size), dim(1), unused(0)
{
	this->data = new double[size*dim];
}


/**
 * Default Constructor.
 * @param size the size
 * @param dim the dimension
 */
DataVector::DataVector(int size, int dim) : size(size), dim(dim), unused(0)
{
	this->data = new double[size*dim];
}

/**
 * Copy constructor.
 * @param vec the vector to be copied
 * @return A copy of the vector vec
 */
DataVector::DataVector(DataVector &vec) : unused(0)
{
	this->size = vec.size;
	this->dim = vec.dim;
	this->data = new double[size*dim];

	memcpy(this->data, vec.data, size*dim*sizeof(double));
}

/**
 * Default Constructor.
 * @param input
 * @param size
 * @param dim
 */
DataVector::DataVector(double * input, int size, int dim) : size(size), dim(dim), unused(0)
{
	this->data = new double[size*dim];
	memcpy(this->data, input, size*dim*sizeof(double));
}

void DataVector::resize(int size)
{
	if(size < this->size)
	{
		return;
	}

	double* newdata = new double[size * this->dim];

	memcpy(newdata, this->data, this->size*this->dim*sizeof(double));

	for(int i = this->size; i < size; i++)
	{
		newdata[i] = 0.0;
	}

	delete [] this->data;

	this->data = newdata;
	this->size = size;
}


void DataVector::addSize(int add)
{
	double* newdata = new double[(size+add) * this->dim];
	memcpy(newdata, this->data, this->size*this->dim*sizeof(double));

	delete [] this->data;

	this->data = newdata;
	this->unused = add;
}


int DataVector::addValue()
{
	if (unused == 0)
	{
		addSize(size);
	}

	int x = size;

	size++;
	unused--;

	return x;
}

/**
 * Sets the vector to value.
 * @param value
 */
void DataVector::setAll(double value)
{
	int n = size*dim;
#ifdef USEOMP
	#pragma omp parallel for schedule(static)
	for(int i = 0; i < n; i++)
	{
		data[i] = value;
	}
#else
	for(int i = 0; i < n; i++)
	{
		data[i] = value;
	}
#endif
}

/**
 * Returns the i-th value.
 * @param i
 * @return
 */
double DataVector::get(int i) const
{
	return data[i];
}

/**
 * Sets the i-th value.
 * @param i
 * @param value
 */
void DataVector::set(int i, double value)
{
	data[i] = value;
}

/**
 * Get the specified row
 * @param row the row to get
 * @param vec result vector
 */
void DataVector::getRow(int row, DataVector& vec)
{
/*#ifdef USEOMP
	#pragma omp parallel for shared(vec) schedule(static)
	for(int i = 0; i < this->dim; i++)
	{
			vec[i] = data[row*dim+i];
	}
#else*/
	for(int i = 0; i < this->dim; i++)
	{
			vec[i] = data[row*dim+i];
	}
//#endif
}

/**
 * Set the specified row
 * @param row the row to set
 * @param vec vector of values
 */
void DataVector::setRow(int row, DataVector& vec)
{
/*#ifdef USEOMP
	#pragma omp parallel for shared(vec) schedule(static)
	for(int i = 0; i < this->dim; i++)
	{
			data[row*dim+i] = vec[i];
	}
#else*/
	for(int i = 0; i < this->dim; i++)
	{
			data[row*dim+i] = vec[i];
	}
//#endif
}

/**
 * Get the specified column
 * @param col the column to get
 * @param vec result vector
 */
void DataVector::getColumn(int col, DataVector& vec)
{
/*#ifdef USEOMP
	#pragma omp parallel for shared(vec) schedule(static)
	for(int j = 0; j < this->size; j++)
	{
		vec[j] = data[j*dim+col];
	}
#else*/
	for(int j = 0; j < this->size; j++)
	{
		vec[j] = data[j*dim+col];
	}
//#endif
}

/**
 * Set the specified column
 * @param col the column to set
 * @param vec vector of values
 */
void DataVector::setColumn(int col, DataVector& vec)
{
/*#ifdef USEOMP
	#pragma omp parallel for shared(vec) schedule(static)
	for(int j = 0; j < this->size; j++)
	{
		data[j*dim+col] = vec[j];
	}
#else*/
	for(int j = 0; j < this->size; j++)
	{
		data[j*dim+col] = vec[j];
	}
//#endif
}

/**
 * Copies a DataVector.
 * @param vec
 * @return
 */
void DataVector::copyFrom(const DataVector& vec)
{
	if(this == &vec)
	{
		return;
	}

	if(size != vec.size || dim != vec.dim)
	{
		delete [] data;
		size = vec.size;
		dim = vec.dim;
		this->data = new double[size*dim];
	}
	memcpy(this->data, vec.data, size*dim*sizeof(double));
}

void DataVector::copySmall(const DataVector& vec)
{
	if(this == &vec)
	{
		return;
	}

	if(vec.dim != 1 || dim != 1 || size < vec.size)
	{
		return;
	}
	memcpy(this->data, vec.data, vec.size*sizeof(double));
}

/**
 * Overloaded assign operator.
 * @param vec
 * @return
 */

DataVector& DataVector::operator=(const DataVector &vec)
{
	if(this == &vec)
	{
		return *this;
	}

	if(size != vec.size || dim != vec.dim)
	{
		delete [] data;
		size = vec.size;
		dim = vec.dim;
		this->data = new double[size*dim];
	}
	memcpy(this->data, vec.data, size*dim*sizeof(double));
	return *this;
}


/**
 * Adds all values of vec to current vector.
 * @param vec
 */
void DataVector::add(DataVector &vec)
{
	if(size != vec.size || dim != vec.dim)
	{
		return;
	}
	int n = size*dim;

/*#ifdef USEOMP
	#pragma omp parallel for shared(vec) schedule(static)
	for(int i = 0; i < n; i++)
	{
		data[i] += vec.data[i];
	}
#else*/
	for(int i = 0; i < n; i++)
	{
		data[i] += vec.data[i];
	}
//#endif
}

/**
 * Substracts all values of vec from current vector.
 * @param vec
 */
void DataVector::sub(DataVector &vec)
{
	if(size != vec.size || dim != vec.dim)
	{
		return;
	}
	int n = size*dim;

/*#ifdef USEOMP
	#pragma omp parallel for shared(vec) schedule(static)
	for(int i = 0; i < n; i++)
	{
		data[i] -= vec.data[i];
	}
#else*/
	for(int i = 0; i < n; i++)
	{
		data[i] -= vec.data[i];
	}
//#endif
}

/**
 * Extract the specified row
 * @param row
 * @param vec
 */
void DataVector::getLine(int row, DataVector& vec)
{
/*#ifdef USEOMP
	#pragma omp parallel for shared(vec) schedule(static)
	for(int i = 0; i < this->dim; i++)
	{
		vec[i] = data[row*dim+i];
	}
#else*/
	for(int i = 0; i < this->dim; i++)
	{
		vec[i] = data[row*dim+i];
	}
//#endif
}

/**
 * Extract the specified row
 * @param row
 * @param vec
 */
void DataVector::getLine(int row, std::vector<double>& vec)
{
	vec.clear();
/*#ifdef USEOMP
	#pragma omp parallel for shared(vec) schedule(static)
	for(int i = 0; i < this->dim; i++)
	{
		vec.push_back(data[row*dim+i]);
	}
#else*/
	for(int i = 0; i < this->dim; i++)
	{
		vec.push_back(data[row*dim+i]);
	}
//#endif
}


/**
 * Returns the length of the vector.
 * @return
 */
int DataVector::getSize()
{
	return size;
}

 /**
  * Returns the number of dimensions.
  * @return
  */
int DataVector::getDim()
{
	return dim;
}

/**
 * Returns the number of allocated values
 * @return
 */
int DataVector::getTotalSize()
{
	return dim*size;
}

/**
 * Returns the dot product of the two vectors. Only defined for 1 dimensional vectors.
 * @param vec
 * @return
 */
double DataVector::dotProduct(DataVector &vec)
{
	double sum = 0.0;

	for(int i = 0; i < size; i++)
	{
		sum += data[i]*vec.data[i];
	}
	return sum;
}

/**
 * Multiplies vector with scalar.
 * @param scalar
 */
void DataVector::mult(double scalar)
{
	int n = size*dim;
/*#ifdef USEOMP
	#pragma omp parallel for schedule(static)
	for(int i = 0; i < n; i++)
	{
		data[i] *= scalar;
	}
#else*/
	for(int i = 0; i < n; i++)
	{
		data[i] *= scalar;
	}
//#endif
}

/**
 * Squares each component
 */
 void DataVector::sqr()
 {
	int n = size*dim;
/*#ifdef USEOMP
	#pragma omp parallel for schedule(static)
	for(int i = 0; i < n; i++)
	{
		data[i] = data[i] * data[i];
	}
#else*/
	for(int i = 0; i < n; i++)
	{
		data[i] = data[i] * data[i];
	}
//#endif
 }

/**
 * Calculates the sum of the vector
 */
 double DataVector::sum()
 {
	int n = size*dim;
	double result;
	for(int i = 0; i < n; i++)
	{
		result += data[i];
	}
	return result;
 }

/**
 * Partitions vector into two classes using a choosen border.
 * @param border
 */
void DataVector::partitionClasses(double border)
{
	int n = size*dim;
	for(int i = 0; i < n; i++)
	{
		data[i] = data[i] > border ? 1.0 : -1.0;
	}
}

/**
 * Adds alpha*x to current vector.
 * BLAS Level 1 (elementary vector operations) operation: axpy.
 * @param alpha
 * @param x
 */
void DataVector::axpy(double alpha, DataVector& x)
{
	if(size != x.size || dim != x.dim)
	{
		return;
	}
	int n = size*dim;
	double* p_x = x.data;
	double* p_d = data;
/*#ifdef USEOMP
	#pragma omp parallel for shared(p_d, p_x) schedule(static)
	for(int i = 0; i < n; i++)
	{
		p_d[i] += alpha*p_x[i];
	}
#else*/
	for(int i = 0; i < n; i++)
	{
		p_d[i] += alpha*p_x[i];
	}
//#endif

}

/**
 * Normalizes d-th dimension with border 0.0
 * @param d
 */
void DataVector::normalizeDimension(int d)
{
	normalizeDimension(d, 0.0);
}


/**
 * Normalizes d-th dimension with border
 * @param d
 * @param border
 */
void DataVector::normalizeDimension(int d, double border)
{
	int n = size*dim;
	double min = data[d];
	double max = data[d];
	for(int i = d; i < n; i += dim)
	{
		if(min > data[i])
		{
			min = data[i];
		}
		if(max < data[i])
		{
			max = data[i];
		}
	}
	min -= border;
	max += border;
	double delta = max - min;
	for(int i = d; i < n; i += dim)
	{
		data[i] = (data[i] - min) / delta;
	}
}

/**
 * Outputs DataVector as string
 * @param text
 */
void DataVector::toString(std::string& text)
{
	std::stringstream str;
	int n = size*dim;

	str << "[";

	for(int i = 0; i < n; i++)
	{
		if(i != 0)
		{
			str << ",";
		}
		str << " " << data[i];
	}
	str << " ]";
	text = str.str();
}

/**
 * Returns the minimum.
 * @param d
 * @return minimum
 */
double DataVector::min(int d)
{
	int n = size*dim;
	double min = data[d];
	for(int i = d; i < n; i += dim)
	{
		if(min > data[i])
		{
			min = data[i];
		}
	}
	return min;
}

/**
 * Returns the maximum.
 * @param d
 * @return maximum
 */
double DataVector::max(int d)
{
	int n = size*dim;
	double max = data[d];
	for(int i = d; i < n; i += dim)
	{
		if(max < data[i])
		{
			max = data[i];
		}
	}
	return max;
}

/**
 * Returns maximum and minimum.
 * @param d
 * @param min
 * @param max
 */
void DataVector::minmax(int d, double* min, double* max)
{
	int n = size*dim;
	double min_t = data[d];
	double max_t = data[d];
	for(int i = d; i < n; i += dim)
	{
		if(min_t > data[i])
		{
			min_t = data[i];
		}
		if(max_t < data[i])
		{
			max_t = data[i];
		}
	}
	(*min) = min_t;
	(*max) = max_t;
}


/**
 * Returns a borrowed Pointer to the data.
 * @return
 */
double* DataVector::getPointer()
{
	return data;
}


DataVector::~DataVector()
{
	delete [] data;
}
