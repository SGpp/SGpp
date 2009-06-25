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


#ifndef DATAVECTOR_H_
#define DATAVECTOR_H_

#include <string>
#include <vector>
#include "data/DataVectorDefinition.hpp"

/**
 * a data holding class of base function's coefficients
 */
class DataVector
{
public:
	/**
	 * Constructor
	 *
	 * @param size number of elements
	 */
	DataVector(size_t size);

	/**
	 * Constructor
	 *
	 * @param size number of elements per dimension
	 * @param dim dimension of Vector
	 */
	DataVector(size_t size, size_t dim);

	/**
	 * Copy Constructor
	 *
	 * @param vec reference to another instance of DataVector
	 */
	DataVector(DataVector& vec);

	/**
	 * Constructor that construct a DataVector from a double array
	 *
	 * @param input double array that contains the data
	 * @param size number of elements per dimension
	 * @param dim number of dimensions
	 */
	DataVector(double* input, size_t size, size_t dim);

	/**
	 * Constructor that constructs a DataVector from a DataVectorDefinition structure
	 *
	 * @param DataVectorDef reference to a DataVectorDefinition structure
	 */
	DataVector(DataVectorDefinition& DataVectorDef);

	/**
	 * resizes the DataVector
	 *
	 * @param size new size of the DataVector
	 */
	void resize(size_t size);

	/**
	 * add elements to the DataVector
	 *
	 * @param add number of elements that should be added to the data vector
	 */
	void addSize(int add);

	/**
	 * adds one element to data vector
	 */
	int addValue();

	/**
	 * Sets all values of DataVector
	 *
	 * @param value value that is set to all elements
	 */
	void setAll(double value);

	/**
	 * copies the data from another DataVector
	 *
	 * @param vec the DataVector containing the data
	 */
	void copyFrom(const DataVector& vec);

	/**
	 * copies the data from another DataVector
	 *
	 * @param vec the DataVector containing the data
	 */
	void copySmall(const DataVector& vec);

	/**
	 * copies the data from another DataVector
	 *
	 * @param vec the DataVector containing the data
	 */
	DataVector& operator=(const DataVector& vec);

	/**
	 * operator to access an element directly
	 *
	 * @param i position of the element
	 */
	inline double& operator[](int i)
	{
		return data[i];
	};

	/**
	 * operator to get an element
	 *
	 * @param i position of the element
	 */
	double get(int i) const;

	/**
	 * operator to set an element
	 *
	 * @param i position of the element
	 * @param value value that should be set
	 */
	void set(int i, double value);

	/**
	 * gets a row of the DataVector
	 *
	 * @param row the row that should be read
	 * @param vec DataVector into which the data is written
	 */
 	void getRow(int row, DataVector& vec);

	/**
	 * sets a row of the DataVector
	 *
	 * @param row the row that should be written
	 * @param vec DataVector from which the data is read
	 */
 	void setRow(int row, DataVector& vec);

	/**
	 * gets a col of the DataVector
	 *
	 * @param col the col that should be read
	 * @param vec DataVector into which the data is written
	 */
 	void getColumn(int col, DataVector& vec);

	/**
	 * sets a row of the DataVector
	 *
	 * @param col the row that should be written
	 * @param vec DataVector from which the data is read
	 */
 	void setColumn(int col, DataVector& vec);

 	/**
 	 * adds the values from another DataVector
 	 *
 	 * @param vec the DataVector which Data is added
 	 */
	void add(DataVector& vec);

	/**
 	 * subs the values of another DataVector
 	 *
 	 * @param vec the DataVector which Data is subtracted
 	 */
	void sub(DataVector& vec);

	/**
	 * multiplies all elements by a constant factor
	 *
	 * @param scalar the constant
	 */
	void mult(double scalar);

	/**
	 * squares all elements of the DataVector
	 */
	void sqr();

	/**
	 * sums all elements up
	 *
	 * @return the sum of all elements
	 */
	double sum();

	/**
	 * Adds alpha*x to current vector.
	 * BLAS Level 1 (elementary vector operations) operation: axpy.
	 *
	 * @param alpha the constant
	 * @param x reference the the DataVector
	 */
	void axpy(double alpha, DataVector& x);

	/**
	 * gets a line of the DataVector
	 *
	 * @param row the line that should be read
	 * @param vec DataVector into which the data is written
	 */
	void getLine(int row, DataVector& vec);

	/**
	 * gets a line of the DataVector
	 *
	 * @param row the line that should be read
	 * @param vec std vector into which the data is written
	 */
	void getLine(int row, std::vector<double>& vec);

	/**
	 * Returns the dot product of the two vectors. Only defined for 1 dimensional vectors.
	 *
	 * @param vec reference to another vector
	 *
	 * @return the dotProduct
	 */
	double dotProduct(DataVector& vec);

	/**
	 * gets the elements stored in the vector
	 *
	 * @return elements stored in the vector
	 */
	size_t getSize();

	/**
	 * get the dimension of the DataVector
	 *
	 * @return dimension of the DataVector
	 */
	size_t getDim();

	/**
	 * gets number of elements in all dimensions
	 *
	 * @return number of elements in all dimensions
	 */
	size_t getTotalSize();

	/**
	 * gets the unsed entries
	 *
	 * @return unsed entries
	 */
	inline int getUnused() { return unused; };

	/**
	 * Partitions vector into two classes using a choosen border.
	 *
	 * @param border value of the border
	 */
	void partitionClasses(double border);

	/**
	 * Normalizes d-th dimension with border 0.0
	 *
	 * @param d the dimension that should be normalized
	 */
	void normalizeDimension(int d);

	/**
	 * Normalizes d-th dimension with border
	 *
	 * @param d the dimension that should be normalized
	 * @param border value ot the border
	 */
	void normalizeDimension(int d, double border);

	/**
	 * Returns the minimum.
	 *
	 * @param d the dimension in which the minimum should be determined
	 *
	 * @return the minimum
	 */
	double min(int d);

	/**
	 * Returns the maximum.
	 *
	 * @param d the dimension in which the maximum should be determined
	 *
	 * @return the maximum
	 */
	double max(int d);

	/**
	 * gets the minimum and the maximum.
	 *
	 * @param d the dimension in which the minimum & maximum should be determined
	 * @param min reference variable for the minimum
	 * @param max reference variable for the maximum
	 */
	void minmax(int d, double* min, double* max);

	/**
	 * Writes the data stored in the Vector into a string
	 *
	 * @param text string to which the data is written
	 */
	void toString(std::string& text);

	/**
	 * gets a pointer to the data array
	 *
	 * @return pointer to the data array
	 */
	double* getPointer();

	/**
	 * Destructor
	 */
	virtual ~DataVector();

private:
	/**
	 * Standard Constructor
	 */
	DataVector();

	/// array to store the data
	double* data;
	/// number of elements in the data vector
	int size;
	/// number of dimensions of one element in this vector
	int dim;
	/// unused slots in this data vector
	int unused;
};

#endif /*DATAVECTOR_H_*/
