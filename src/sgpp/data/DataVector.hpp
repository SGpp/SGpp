/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (Dirk.Pflueger@in.tum.de)


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
	 * Constructor: one-dimensional DataVector with <it>size</it> rows.
	 *
	 * @param size number of elements
	 */
	DataVector(size_t size);

	/**
	 * Constructor: multi-dimensional DataVector with <it>size</it> rows and
	 * <it>dim</it> columns.
	 *
	 * @param size number of elements per dimension
	 * @param dim dimension of Vector
	 */
	DataVector(size_t size, size_t dim);

	/**
	 * Copy Constructor.
	 *
	 * @param vec reference to another instance of DataVector
	 */
	DataVector(DataVector& vec);

	/**
	 * Copy Constructor.
	 *
	 * @param vec reference to another instance of DataVector
	 */
	DataVector(const DataVector& vec);

	/**
	 * Constructor that construct a DataVector from a double array.
	 * The double array contains the entries rowwise: x00,x01,...,x0dim-1,x10,x11,...
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
	 * Fill DataVectorDefintion with data from the current DataVector
	 *
	 * @param DataVectorDef reference to a DataVectorDefinition struct
	 */
	void getDataVectorDefinition(DataVectorDefinition& DataVectorDef);

	/**
	 * Resizes the DataVector to size*dim elements; sets all new entries to zero.
	 * Does nothing if size smaller than current size.
	 *
	 * @param size new size of the DataVector
	 */
	void resize(size_t size);

	/**
	 * Resizes the DataVector by removing entries. Throws an exception
	 * if boundaries a violated.
	 *
	 * @param remainingIndex vector that contains the remaining indices of the DataVector
	 */
	void restructure(std::vector<size_t>& remainingIndex);

	/**
	 * Add add potentially new elements to the DataVector. The size remains unchanged.
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
	 * Copies the data from another DataVector vec.
	 * Afterwards, the current vector is an exact copy of vec.
	 * If the current DataVector has not the sime dimensions
	 * (size, dim) than vec, it has the same effect than the copy
	 * constructor: The old memory is discarded, new memory
	 * allocated, and the data copied.
	 *
	 * @param vec the DataVector containing the data
	 */
	void copyFrom(const DataVector& vec);

	/**
	 * Copies the data from another, smaller DataVector vec.
	 * Has no effect, if vec has more elements than the current
	 * vector. If it has the same or smaller size, the first
	 * vec.size() elements of the current vector are overwritten.
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
     * Multiplies the current DataVector component-wise with another DataVector.
     * Performs
     * @code
     * for i from 1 to this.getTotalSize()
     *   this[i] *= vec[i]
     * @endcode
     *
     * @param vec the DataVector which is multiplied to current DataVector
     */
    void componentwise_mult(DataVector& vec);

    /**
     * Divides the current DataVector component-wise by another DataVector.
     * Performs
     * @code
     * for i from 1 to this.getTotalSize()
     *   this[i] /= vec[i]
     * @endcode
     * Note: <b>No check for division by zero!</b>
     *
     * @param vec the DataVector which the current DataVector is divided by
     */
    void componentwise_div(DataVector& vec);

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
	 * Sets all elements to their absolute value.
	 *
	 */
	void abs();

	/**
	 * sums all elements up
	 *
	 * @return the sum of all elements
	 */
	double sum();

	/**
	 * calculates the vector's max norm
	 *
	 * @return the vector's max norm
	 */
	double maxNorm();

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

#ifndef LARRABEE
	/**
	 * Returns the minimum.
	 *
	 * @param d the dimension in which the minimum should be determined
	 *
	 * @return the minimum
	 */
	double min(int d);

	/**
	 * Returns the minimum over all entries.
	 *
	 * @return global minimum
	 */
	double min();

	/**
	 * Returns the maximum.
	 *
	 * @param d the dimension in which the maximum should be determined
	 *
	 * @return the maximum
	 */
	double max(int d);

	/**
	 * Returns the maximum over all entries.
	 *
	 * @return global maximum
	 */
	double max();

	/**
	 * gets the minimum and the maximum.
	 *
	 * @param d the dimension in which the minimum & maximum should be determined
	 * @param min reference variable for the minimum
	 * @param max reference variable for the maximum
	 */
	void minmax(int d, double* min, double* max);
#endif

	/**
	 * Writes the data stored in the DataVector into a string
	 *
	 * @param text string to which the data is written
	 */
	void toString(std::string& text);

	/**
	 * Returns a description of the DataVector as a string.
	 *
	 * @returns string of the DataVector
	 */
  std::string toString();

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

	/**
	 * gets the number of none zero elements in the vector
	 *
	 * @return the number of none zero elements
	 */
	size_t getNumberNonZero();

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
