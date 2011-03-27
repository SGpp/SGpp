/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (Dirk.Pflueger@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)


#ifndef DATAMATRIX_H_
#define DATAMATRIX_H_

#include <string>
#include <vector>
#include "data/DataVector.hpp"

/**
 * A class to store two-dimensional data.
 * Typically, DataMatrix would contain a set of (d-dimensional) data or evaluation points, i.e.,
 * the DataMatrix consists of d columns, and each row is one of the points.
 * Thus, typical functionality like obtaining the maximum for a certain dimension (or attribute),
 * or normalizing all data points to the unit interval for a certain dimension are
 * provided.
 */
class DataMatrix
{
public:
	/**
	 * Create a two-dimensional DataMatrix with @em nrows rows and
	 * @em ncols columns.
	 *
	 * @param nrows Number of rows
	 * @param ncols Number of columns
	 */
	DataMatrix(size_t nrows, size_t ncols);

	/**
	 * Create a new DataMatrix that is a copy of matr.
	 *
	 * @param matr Reference to another instance of DataMatrix
	 */
	DataMatrix(const DataMatrix& matr);

	/**
	 * Create a new DataMatrix from a double array.
	 * The double array contains the entries row-wise:
	 * x0_0,x0_1,...,x0_ncol-1,
	 * x1_0,x1_1,...
	 * ...
	 * xnrow_0, xnrow_1,...,xnrow_ncol-1
	 *
	 * @param input double array that contains the data
	 * @param nrows number of rows
	 * @param ncols number of columns
	 */
	DataMatrix(double* input, size_t nrows, size_t ncols);


	/**
	 * Resizes the DataVector to size elements.
	 * All new additional entries are uninitialized.
	 * If nrows is smaller than the current number of rows,
	 * all superfluous entries are removed.
	 *
	 * @param nrows New number of rows of the DataMatrix
	 */
	void resize(size_t nrows);

    /**
     * Resizes the DataMatrix to nrows rows.
     * All new additional entries are set to zero.
     * If nrows is smaller than the current number of rows,
     * all superfluous entries are removed.
     *
     * @param nrows New number of rows of the DataMatrix
     */
    void resizeZero(size_t nrows);

	/**
	 * Reserves memory for potentially inc_nrows new rows;
	 * the actual number of rows remains unchanged.
	 * Corresponds to a resize to nrows+inc_nrows new rows while leaving
	 * the current matrix' size unchanged.
	 *
	 * @param inc_nrows Number of additional rows for which storage is to be reserved.
	 */
	void addSize(size_t inc_nrows);

	/**
	 * Appends a new row and returns index of it.
	 * If the new row does not fit into the reserved memory,
	 * reserves memory for getIncRows() additional rows.
	 *
	 * @return Index of new row
	 */
	size_t appendRow();

	/**
     * Appends a new row with data contained in DataVector vec
     * and returns index of new row.
     * If the new row does not fit into the reserved memory,
     * reserves memory for getIncRows() additional rows.
     *
	 * @param vec DataVector (length has to match getNcols()) with data
	 * @return Index of new row
	 */
    size_t appendRow(DataVector& vec);


    /**
	 * Sets all entries of DataMatrix to value.
	 *
	 * @param value New value for all entries
	 */
	void setAll(double value);

	/**
	 * Copies the data from another DataMatrix matr.
	 * Disregards the number of rows and columns set for the two matrices,
	 * i.e., just copies the data entry by entry (and row by row).
	 * If the dimensions match (nrows, ncols), the current DataMatrix is an
	 * exact copy of matr. If not, as many elements as possible are
	 * copied, and everything else is left untouched.
	 *
	 * @param matr The source DataMatrix containing the data
	 */
	void copyFrom(const DataMatrix& matr);

	/**
	 * Transposes this DataMatrix
	 */
	void transpose();

	/**
	 * Copies the data from another DataMatrix.
	 * Dimensions have to match.
	 *
	 * @param matr the DataMatrix containing the data
	 */
	DataMatrix& operator=(const DataMatrix& matr);

	/**
	 * Returns the i-th element.
	 * For the 5th element in the third row, i would be 2*getNcols()+4.
	 *
	 * @param i position of the element
	 */
	inline double& operator[](size_t i)
	{
		return data[i];
	};

	/**
	 * Returns the value of the element at position [row,col]
	 *
	 * @param row Row
	 * @param col Column
	 * @return Value of the element
	 */
	inline double get(size_t row, size_t col) const {
	    return data[row*ncols+col];
	};

	/**
	 * Sets the element at position [row,col] to value.
	 *
	 * @param row Row
	 * @param col Column
	 * @param value New value for element
	 */
	inline void set(size_t row, size_t col, double value) {
	    data[row*ncols+col] = value;
	};

	/**
	 * Copies the values of a row to the DataVector vec.
	 *
	 * @param row The row
	 * @param vec DataVector into which the data is written
	 */
 	void getRow(size_t row, DataVector& vec);

	/**
	 * Copies the values of a row to the std::vector vec.
	 *
	 * @param row The row
	 * @param vec std::vector into which the data is written
	 */
 	void getRow(size_t row, std::vector<double>& vec);

	/**
	 * Sets a row of the DataMatrix to the values of a DataVector vec.
	 *
	 * @param row The row which is to be overwritten
	 * @param vec DataVector containing the data of the row
	 */
 	void setRow(size_t row, DataVector& vec);

	/**
     * Copies the values of a column to the DataVector vec.
	 *
	 * @param col The column
	 * @param vec DataVector into which the data is written
	 */
 	void getColumn(size_t col, DataVector& vec);

	/**
     * Sets a column of the DataMatrix to the values of a DataVector vec.
	 *
	 * @param col The column which is to be overwritten
	 * @param vec DataVector containing the data of the column
	 */
 	void setColumn(size_t col, DataVector& vec);


 	/**
 	 * Adds the values from another DataMatrix to the current values.
 	 * Modifies the current values.
 	 *
 	 * @param matr The DataMatrix which is added to the current values
 	 */
	void add(DataMatrix& matr);

	/**
     * Subtracts the values from another DataMatrix of the current values.
     * Modifies the current values.
 	 *
 	 * @param matr The DataMatrix which is subtracted from the current values
 	 */
	void sub(DataMatrix& matr);

    /**
     * Multiplies the current DataMatrix component-wise with another DataMatrix.
     * Modifies the current values.
     * Performs
     * @code
     * for i from 1 to this.getSize()
     *   this[i] *= matr[i]
     * @endcode
     *
     * @param matr the DataMatrix which is multiplied to current DataMatrix
     */
    void componentwise_mult(DataMatrix& matr);

    /**
     * Divides the current DataMatrix component-wise by another DataMatrix.
     * Modifies the current values.
     * Performs
     * @code
     * for i from 1 to this.getTotalSize()
     *   this[i] /= matr[i]
     * @endcode
     * Note: <b>No check for division by zero!</b>
     *
     * @param matr the DataMatrix which the current DataMatrix is divided by
     */
    void componentwise_div(DataMatrix& matr);

	/**
	 * Multiplies all elements by a constant factor
	 *
	 * @param scalar the constant
	 */
	void mult(double scalar);

	/**
	 * Squares all elements of the DataMatrix
	 */
	void sqr();

	/**
	 * Takes the square root of all elements of the DataMatrix
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
	 * Returns the minimum value of column col.
	 *
	 * @param col Number of the column
	 *
	 * @return Minimum value
	 */
    double min(size_t col);

    /**
	 * Returns the minimum over all entries.
	 *
	 * @return Minimal value of all entries
	 */
    double min();

    /**
	 * Returns the maximum value of column col.
	 *
	 * @param col Number of the column
	 *
	 * @return Maximum value
	 */
    double max(size_t col);

    /**
	 * Returns the maximum over all entries.
	 *
	 * @return Maximal value of all entries
	 */
    double max();

    /**
	 * Determines minimum and maximum of column col.
	 *
	 * @param col Number of the column
	 * @param min Reference variable for the minimum
	 * @param max Reference variable for the maximum
	 */
    void minmax(size_t col, double *min, double *max);

    /**
     * Determines minimum and maximum over all entries.
     *
     * @param min Reference variable for the minimum
     * @param max Reference variable for the maximum
     */
    void minmax(double *min, double *max);


    /**
 	 * Returns pointer to double array containing underlying data.
 	 *
 	 * @return Pointer to data
 	 */
 	double* getPointer();

	/**
	 * Returns the total number of (used) elements, i.e., getNrows()*getNCols()
	 *
	 * @return Number of elements stored in the matrix
	 */
	inline size_t getSize() {
		return ncols * nrows;
	};

	/**
	 * Returns the number of unused rows.
	 *
	 * @return number of unused rows
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
     * Returns the number of rows of the DataMatrix.
     *
     * @return Number of rows
     */
    inline size_t getNrows() const
    {
        return nrows;
    }

    /**
     * Returns the number of columns of the DataMatrix.
     *
     * @return Number of columns
     */
    inline size_t getNcols() const
    {
        return ncols;
    }

    /**
     * Get the current number of rows by which the DataMatrix is extended,
     * if appendRow() is called and no unused rows are left
     *
     * @return Row increment
     */
    inline size_t getInc() const
    {
        return inc_rows;
    }

    /**
     * Sets the current number of rows by which the DataMatrix is extended,
     * if appendRow() is called and no unused rows are left.
     * Defaults to 100.
     *
     * @param inc_rows Row increment
     */
    void setInc(size_t inc_rows)
    {
        this->inc_rows = inc_rows;
    }


	/**
	 * Normalizes the d-th dimension (entries in the d-th column) to @f$[0,1]@f$.
	 * Considers contents of DataMatrix as a d-dimensional dataset, one
	 * data point per row.
	 *
	 * @param d The dimension (column) that should be normalized (starting with 0)
	 */
	void normalizeDimension(size_t d);

	/**
     * Normalizes the d-th dimension (entries in the d-th column) to @f$[border,1-border]@f$.
     * Considers contents of DataMatrix as a d-dimensional dataset, one
     * data point per row.
     *
     * @param d The dimension (column) that should be normalized (starting with 0)
	 * @param border Width of the border
	 */
	void normalizeDimension(size_t d, double border);


    /**
	 * Writes the data stored in the DataMatrix into a string
	 *
	 * @param text String to which the data is written
	 */
    void toString(std::string& text);

    /**
	 * Returns a description of the DataMatrix as a string.
	 *
	 * @returns string of the DataMatrix
	 */
    std::string toString();

    /**
	 * Destructor
	 */
    virtual ~DataMatrix();

private:
	/// Pointer to the data
	double* data;
	/// Number of rows of the data matrix
	size_t nrows;
	/// Number of columns of the data matrix
	size_t ncols;
	/// Number of additional rows for which memory has already been reserved
	size_t unused;
    /// Number of rows by which the reserved memory is increased, if adding a row would exceed the storage reserved so far.
    size_t inc_rows;
};

#endif /*DATAMATRIX_H_*/
