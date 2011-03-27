/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Dirk Pflueger (Dirk.Pflueger@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)


#ifndef DATAVECTOR_H_
#define DATAVECTOR_H_

#include <string>
#include <vector>
#include "data/DataVectorDefinition.hpp"

/**
 * A class to store one-dimensional data.
 * Typically, an object of type DataVector will contain an array
 * of (hierarchical) coefficients (or surplusses), or the coordinates
 * of a data point at which a sparse grid function should be
 * evaluated.
 */
class DataVector
{
public:
	/**
	 * Create a DataVector with @em size elements.
	 *
	 * @param size Number of elements
	 */
	DataVector(size_t size);

	/**
	 * Constructor: multi-dimensional DataVector with <it>size</it> rows and
	 * <it>dim</it> columns.
	 *
	 * @param size number of elements per dimension
	 * @param dim dimension of Vector
	 */
//	DataVector(size_t size, size_t dim);

	/**
	 * Create a new DataVector that is a copy of vec.
	 *
	 * @param vec Reference to another instance of DataMatrix
	 */
	DataVector(DataVector& vec);

	/**
	 * Create a new DataVector that is a copy of vec.
	 *
	 * @param vec Reference to another instance of DataMatrix
	 */
	DataVector(const DataVector& vec);

	/**
	 * Create a new DataVector from a double array with size elements.
	 *
	 * @param input double array that contains the data
	 * @param size number of elements
	 */
	DataVector(double* input, size_t size);

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
	 * Resizes the DataVector to size elements.
	 * All new additional entries are uninitialized.
	 * If nrows is smaller than the current number of rows,
	 * all superfluous entries are removed.
	 *
	 * @param size New number of elements of the DataVector
	 * @todo (pflueged) Check that no wrong usage of reize is left which assumes that new entries are zero.
	 * @todo (pflueged) Optimize implementation, consider unused elements!
	 */
	void resize(size_t size);

    /**
	 * Resizes the DataVector to size elements.
     * All new additional entries are set to zero.
     * If nrows is smaller than the current number of rows,
     * all superfluous entries are removed.
     *
     * @param nrows New number of rows of the DataMatrix
     */
    void resizeZero(size_t nrows);

	/**
	 * Resizes the DataVector by removing entries. Throws an exception
	 * if boundaries a violated.
	 *
	 * @param remainingIndex vector that contains the remaining indices of the DataVector
	 */
	void restructure(std::vector<size_t>& remainingIndex);

	/**
	 * Add add potentially new elements to the DataVector. The size remains unchanged
	 * Reserves memory for potentially inc_elems new elements;
	 * the actual number of elements remains unchanged.
	 * Corresponds to a resize to size+inc_elems new elements while leaving
	 * the current vector's size unchanged.
	 *
	 * @param inc_nrows Number of additional elements for which storage is to be reserved.
	 */
	void addSize(size_t inc_elems);

	/**
	 * Appends a new element and returns index of it.
	 * If the new element does not fit into the reserved memory,
	 * reserves memory for getInc() additional elements.
	 *
	 * @return Index of new element
	 */
	size_t append();

	/**
     * Appends a new element and returns index of new element.
     * If the new element does not fit into the reserved memory,
     * reserves memory for getInc() additional elements.
     *
	 * @param value Value of new element
	 * @return Index of new element
	 */
    size_t append(double value);


	/**
	 * Sets all values of DataVector to value
	 *
	 * @param value New value for all entries
	 */
	void setAll(double value);

	/**
	 * Copies the data from another DataVector vec.
	 * Disregards the number of entries set for the two vectors,
	 * i.e., just copies the data entry by entry.
	 * If the size matches, the current DataVector is an
	 * exact copy of vec. If not, as many elements as possible are
	 * copied, and everything else is left untouched.
	 *
	 * @param vec The source DataVector containing the data
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
//	void copySmall(const DataVector& vec);

	/**
	 * Copies the data from another DataVector.
	 * Dimensions have to match.
	 *
	 * @param vec the DataVector containing the data
	 */
	DataVector& operator=(const DataVector& vec);

	/**
	 * Returns the i-th element.
	 *
	 * @param i position of the element
	 */
	inline double& operator[](size_t i)
	{
		return data[i];
	};

	/**
	 * Returns the i-th element.
	 *
	 * @param i position of the element
	 */
	inline double get(size_t i) const {
		return data[i];
	}

	/**
	 * Sets the element at index i to value.
	 *
	 * @param i Index
	 * @param value New value for element
	 */
	void set(size_t i, double value);

	/**
	 * gets a row of the DataVector
	 *
	 * @param row the row that should be read
	 * @param vec DataVector into which the data is written
	 */
// 	void getRow(size_t row, DataVector& vec) const;

	/**
	 * sets a row of the DataVector
	 *
	 * @param row the row that should be written
	 * @param vec DataVector from which the data is read
	 */
// 	void setRow(int row, DataVector& vec);

	/**
	 * gets a col of the DataVector
	 *
	 * @param col the col that should be read
	 * @param vec DataVector into which the data is written
	 */
// 	void getColumn(int col, DataVector& vec) const;

	/**
	 * sets a row of the DataVector
	 *
	 * @param col the row that should be written
	 * @param vec DataVector from which the data is read
	 */
// 	void setColumn(int col, DataVector& vec);


 	/**
 	 * Adds the values from another DataVector to the current values.
 	 * Modifies the current values.
 	 *
 	 * @param vec The DataVector which is added to the current values
 	 */
	void add(DataVector& vec);

	/**
     * Subtracts the values from another DataMatrix of the current values.
     * Modifies the current values.
 	 *
 	 * @param vec The DataMatrix which is subtracted from the current values
 	 */
	void sub(DataVector& vec);

    /**
     * Multiplies the current DataVector component-wise with another DataVector.
     * Modifies the current values.
     * Performs
     * @code
     * for i from 1 to this.getSize()
     *   this[i] *= vec[i]
     * @endcode
     *
     * @param vec the DataVector which is multiplied to current DataVector
     */
    void componentwise_mult(DataVector& vec);

    /**
     * Divides the current DataVector component-wise by another DataVector.
     * Modifies the current values.
     * Performs
     * @code
     * for i from 1 to this.getSize()
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
	 * Squares all elements of the DataVector
	 */
	void sqr();

	/**
	 * Takes the square root of all elements of the DataVector
	 */
	void sqrt();

	/**
	 * Sets all elements to their absolute value.
	 *
	 */
	void abs();

	/**
	 * Returns the sum of all elements
	 *
	 * @return The sum of all elements
	 */
	double sum();

	/**
	 * calculates the vector's max norm
	 *
	 * @return the vector's max norm
	 */
	double maxNorm();

	/**
	 * calculates the vectors L two norm
	 * function based two norm
	 *
	 * @return the vector's L two norm
	 */
	double RMSNorm();

	/**
	 * calculates the vectors two norm
	 *
	 * @return the vector's two norm
	 */
	double l2Norm();

#ifndef LARRABEE
	/**
	 * Returns the minimum over all entries.
	 *
	 * @return global minimum
	 */
	double min();

	/**
	 * Returns the maximum over all entries.
	 *
	 * @return global maximum
	 */
	double max();

	/**
     * Determines minimum and maximum over all entries.
     *
     * @param min Reference variable for the minimum
     * @param max Reference variable for the maximum
	 */
	void minmax(double* min, double* max);
#endif

	/**
	 * Adds a*x to current vector.
	 * BLAS Level 1 (elementary vector operations) operation: axpy.
	 *
	 * @param a A scalar
	 * @param x Reference to the DataVector
	 */
	void axpy(double a, DataVector& x);

	/**
	 * gets a line of the DataVector
	 *
	 * @param row the line that should be read
	 * @param vec DataVector into which the data is written
	 */
//	void getLine(int row, DataVector& vec);

	/**
	 * gets a line of the DataVector
	 *
	 * @param row the line that should be read
	 * @param vec std vector into which the data is written
	 */
//	void getLine(int row, std::vector<double>& vec);

	/**
	 * Returns the dot product of the two vectors.
	 *
	 * @param vec Reference to another vector
	 *
	 * @return The dot-product
	 */
	double dotProduct(DataVector& vec);


	/**
	 * gets a pointer to the data array
	 *
	 * @return pointer to the data array
	 */
	double* getPointer();

	/**
	 * gets the elements stored in the vector
	 *
	 * @return elements stored in the vector
	 */
	inline size_t getSize() const {
		return size;
	};

	/**
	 * Returns the number of unused elements.
	 *
	 * @return number of unused elements
	 */
	inline size_t getUnused() {
		return unused;
	};

	/**
	 * Determines the number of non-zero elements in the vector.
	 *
	 * @return The number of non-zero elements
	 */
	size_t getNumberNonZero();

    /**
     * Get the current number of elements by which the DataVector is extended,
     * if append() is called and no unused rows are left
     *
     * @return Increment
     */
    inline size_t getInc() const
    {
        return inc_elems;
    }

    /**
     * Sets the current number of elements by which the DataVector is extended,
     * if append() is called and no unused elements are left.
     * Defaults to 100.
     *
     * @param inc_elems Increment
     */
    void setInc(size_t inc_elems)
    {
        this->inc_elems = inc_elems;
    }


	/**
	 * get the dimension of the DataVector
	 *
	 * @return dimension of the DataVector
	 */
//	size_t getDim() const;

	/**
	 * gets number of elements in all dimensions
	 *
	 * @return number of elements in all dimensions
	 */
//	size_t getTotalSize() const;

	/**
	 * Partitions vector into two classes using a choosen border.
	 *
	 * @param threshold value of the border
	 */
	void partitionClasses(double threshold);

	/**
	 * Normalizes d-th dimension with border 0.0
	 *
	 * @param d the dimension that should be normalized
	 */
	void normalize();

	/**
	 * Normalizes d-th dimension with border
	 *
	 * @param d the dimension that should be normalized
	 * @param border value ot the border
	 */
	void normalize(double border);


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
	 * Destructor
	 */
	virtual ~DataVector();



private:
	/**
	 * Standard Constructor
	 */
//	DataVector();

	/// Array to store the data
	double* data;
	/// Number of elements of the data vector
	size_t size;
	/// number of dimensions of one element in this vector
//	int dim;
	/// Number of additional rows for which memory has already been reserved
	size_t unused;
    /// Number of elements by which the reserved memory is increased, if adding a row would exceed the storage reserved so far.
    size_t inc_elems;
};

#endif /*DATAVECTOR_H_*/
