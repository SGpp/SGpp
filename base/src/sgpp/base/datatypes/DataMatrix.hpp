// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DATAMATRIX_H_
#define DATAMATRIX_H_

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

/**
 * A class to store two-dimensional data.
 * Typically, DataMatrix would contain a set of (d-dimensional) data or evaluation points, i.e.,
 * the DataMatrix consists of d columns, and each row is one of the points.
 * Thus, typical functionality like obtaining the maximum for a certain dimension (or attribute),
 * or normalizing all data points to the unit interval for a certain dimension are
 * provided.
 */
class DataMatrix : public std::vector<double> {
 public:
  /**
   * Creates an empty two-dimensional DataMatrix.
   */
  DataMatrix();

  /**
   * Copy constructor.
   */
  DataMatrix(const DataMatrix&) = default;

  // /**
  //  * Move constructor
  //  */
  // DataMatrix(DataMatrix&&) = default;

  /**
   * Copy assignment operator
   */
  DataMatrix& operator=(const DataMatrix&) = default;

  // /**
  //  * Move assignment operator
  //  */
  // DataMatrix& operator=(DataMatrix&&) = default;

  /**
   * Destructor
   */
  ~DataMatrix() = default;

  /**
   * Create a two-dimensional DataMatrix with @em nrows rows and
   * @em ncols columns (uninitialized values).
   *
   * @param nrows Number of rows
   * @param ncols Number of columns
   */
  DataMatrix(size_t nrows, size_t ncols);

  /**
   * Create a two-dimensional DataMatrix with @em nrows rows and
   * @em ncols columns and initializes all elements with the same value.
   *
   * @param nrows Number of rows
   * @param ncols Number of columns
   * @param value Value for all entries
   */
  DataMatrix(size_t nrows, size_t ncols, double value);

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
  DataMatrix(const double* input, size_t nrows, size_t ncols);

  /**
   * Create a new DataMatrix from a std::vector<double>.
   *
   * @param input std::vector<double> that contains the data
   */
  explicit DataMatrix(std::vector<double> input, size_t nrows);

  /**
   * Create a new DataMatrix from a std::initializer_list<double>.
   *
   * @param input std::initializer_list<double> that contains the data
   */
  explicit DataMatrix(std::initializer_list<double> input, size_t nrows);

  static DataMatrix fromFile(const std::string& fileName);

  static DataMatrix fromString(const std::string& serializedVector);

  /**
   * Resizes the DataMatrix to nrows rows.
   * All new additional entries are uninitialized.
   * If nrows is smaller than the current number of rows,
   * all superfluous entries are removed.
   * \deprecated use resizeRows
   *
   * @param nrows New number of rows of the DataMatrix
   */
  void resize(size_t nrows);

  /**
   * Resizes the DataMatrix to nrows rows.
   * All new additional entries are uninitialized.
   * If nrows is smaller than the current number of rows,
   * all superfluous entries are removed.
   *
   * @param nrows New number of rows of the DataMatrix
   */
  void resizeRows(size_t nrows);

  /**
   * Resizes the DataMatrix to nrows rows and ncols columns.
   * All new additional entries are uninitialized.
   * If nrows*ncols is smaller than the current number of entries,
   * all superfluous entries are removed.
   * \deprecated use resizeRowsCols
   *
   * @param nrows New number of rows of the DataMatrix
   * @param ncols New number of columns of the DataMatrix
   */
  void resize(size_t nrows, size_t ncols);

  /**
   * Resizes the DataMatrix to nrows rows and ncols columns.
   * All new additional entries are uninitialized.
   * If nrows*ncols is smaller than the current number of entries,
   * all superfluous entries are removed.
   *
   * @param nrows New number of rows of the DataMatrix
   * @param ncols New number of columns of the DataMatrix
   */
  void resizeRowsCols(size_t nrows, size_t ncols);

  /**
   * Resizes the quadratic DataMatrix to size rows and size columns.
   * All new additional entries are uninitialized.
   * If size is smaller than the current size of the quadratic DataMatrix,
   * all superfluous entries are removed.
   *
   * @param size New dimension of quadratic data DataMatrix
   */
  void resizeQuadratic(size_t size);

  /**
   * Resizes the DataMatrix to nrows rows.
   * All new additional entries are set to zero.
   * If nrows is smaller than the current number of rows,
   * all superfluous entries are removed.
   * \deprecated use resizeRows
   *
   * @param nrows New number of rows of the DataMatrix
   */
  void resizeZero(size_t nrows);

  /**
   * Resizes the DataMatrix to nrows rows and ncols columns.
   * All new additional entries are set to zero.
   * If nrows*ncols is smaller than the current number of entries,
   * all superfluous entries are removed.
   * \deprecated use resizeRowsCols
   *
   * @param nrows New number of rows of the DataMatrix
   * @param ncols New number of columns of the DataMatrix
   */
  void resizeZero(size_t nrows, size_t ncols);

  /**
   * Resize current matrix to the submatrix Mat[row_1:row_2, col_1:col_2].
   *
   * @param row_1, col_1 corresponding to left upper index of desired submatrix
   * @param row_2, col_2 corresponding to right lower index of desired submatrix
   */
  void resizeToSubMatrix(size_t row_1, size_t col_1, size_t row_2, size_t col_2);

  /**
   * Reserves memory for potentially inc_nrows new rows;
   * the actual number of rows remains unchanged.
   * Corresponds to a resize to nrows+inc_nrows new rows while leaving
   * the current matrix' size unchanged.
   *
   * @param inc_nrows Number of additional rows for which storage is to be reserved.
   */
  void reserveAdditionalRows(size_t inc_nrows);

  /**
   * Appends a new row and returns index of it.
   *
   * @return Index of new row
   */
  size_t appendRow();

  /**
   * Appends a new row with data contained in DataVector vec
   * and returns index of new row.
   *
   * @param vec DataVector (length has to match getNcols()) with data
   * @return Index of new row
   */
  size_t appendRow(const DataVector& vec);

  /**
   * Appends a new Col with data contained in DataVector vec
   * and returns index of new col.
   *
   * @param vec DataVector (length has to match getNcols()) with data
   * @return Index of new col
   */
  size_t appendCol(const DataVector& vec);

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
   * Returns the value of the element at position [row,col]
   *
   * @param row Row
   * @param col Column
   * @return reference to the element
   */
  inline double& operator()(size_t row, size_t col) { return (*this)[row * ncols + col]; }

  /**
   * Returns the value of the element at position [row,col]
   *
   * @param row Row
   * @param col Column
   * @return constant reference to the element
   */
  inline const double& operator()(size_t row, size_t col) const {
    return (*this)[row * ncols + col];
  }

  /**
   * Returns the value of the element at position [row,col]
   *
   * @param row Row
   * @param col Column
   * @return Value of the element
   */
  inline double get(size_t row, size_t col) const { return (*this)[row * ncols + col]; }

  /**
   * Sets the element at position [row,col] to value.
   *
   * @param row Row
   * @param col Column
   * @param value New value for element
   */
  inline void set(size_t row, size_t col, double value) { (*this)[row * ncols + col] = value; }

  /**
   * Copies the values of a row to the DataVector vec.
   *
   * @param row The row
   * @param vec DataVector into which the data is written
   */
  void getRow(size_t row, DataVector& vec) const;

  /**
   * Copies the values of a row to the std::vector vec.
   *
   * @param row The row
   * @param vec std::vector into which the data is written
   */
  void getRow(size_t row, std::vector<double>& vec) const;

  /**
   * Sets a row of the DataMatrix to the values of a DataVector vec.
   *
   * @param row The row which is to be overwritten
   * @param vec DataVector containing the data of the row
   */
  void setRow(size_t row, const DataVector& vec);

  /**
   * Copies the values of a column to the DataVector vec.
   *
   * @param col The column
   * @param vec DataVector into which the data is written
   */
  void getColumn(size_t col, DataVector& vec) const;

  /**
   * Sets a column of the DataMatrix to the values of a DataVector vec.
   *
   * @param col The column which is to be overwritten
   * @param vec DataVector containing the data of the column
   */
  void setColumn(size_t col, const DataVector& vec);

  /**
   * Adds the values from another DataMatrix to the current values.
   * Modifies the current values.
   *
   * @param matr The DataMatrix which is added to the current values
   */
  void add(const DataMatrix& matr);

  /**
   * Subtracts the values from another DataMatrix of the current values.
   * Modifies the current values.
   *
   * @param matr The DataMatrix which is subtracted from the current values
   */
  void sub(const DataMatrix& matr);

  /**
   * Reduce the DataMatrix along the
   * columns by adding all entries in one row.
   *
   * @param reduction DataVector into which the reduce columns are stored
   */
  void addReduce(DataVector& reduction);

  /**
   * Reduce the DataMatrix along the
   * columns by adding all entries in one row.
   *
   * @param reduction DataVector to which the reduce columns are added
   * @param beta vector with length of number of columns beta[i] is multiplied to each element
   *row[j][i]
   * @param start_beta where to start using the beta coefficients
   */
  void addReduce(DataVector& reduction, DataVector& beta, size_t start_beta);

  /**
   * expands a given DataVector into a
   * DataMatrix.
   *
   * @param expand DataVector that should be expanded
   */
  void expand(const DataVector& expand);

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
  void componentwise_mult(const DataMatrix& matr);

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
  void componentwise_div(const DataMatrix& matr);

  /**
   * Multiplies all elements by a constant factor
   *
   * @param scalar the constant
   */
  void mult(double scalar);

  /**
   * Multiplies the matrix with a vector x and stores the result
   * in another vector y.
   *
   * @param[in] x vector to be multiplied
   * @param[out] y vector in which the result should be stored
   */
  void mult(const DataVector& x, DataVector& y);

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
  double sum() const;

  /**
   * Returns the minimum value of column col.
   *
   * @param col Number of the column
   *
   * @return Minimum value
   */
  double min(size_t col) const;

  /**
   * Returns the minimum over all entries.
   *
   * @return Minimal value of all entries
   */
  double min() const;

  /**
   * Returns the maximum value of column col.
   *
   * @param col Number of the column
   *
   * @return Maximum value
   */
  double max(size_t col) const;

  /**
   * Returns the maximum over all entries.
   *
   * @return Maximal value of all entries
   */
  double max() const;

  /**
   * Determines minimum and maximum of column col.
   *
   * @param col Number of the column
   * @param min Reference variable for the minimum
   * @param max Reference variable for the maximum
   */
  void minmax(size_t col, double* min, double* max) const;

  /**
   * Determines minimum and maximum over all entries.
   *
   * @param min Reference variable for the minimum
   * @param max Reference variable for the maximum
   */
  void minmax(double* min, double* max) const;

  /**
   * Returns pointer to double array containing underlying data.
   *
   * @return Pointer to data
   */
  double* getPointer();

  /**
   * Returns const pointer to double array containing underlying data.
   *
   * @return Const pointer to data
   */
  const double* getPointer() const;

  /**
   * Returns the total number of (used) elements, i.e., getNrows()*getNCols()
   *
   * @return Number of elements stored in the matrix
   */
  inline size_t getSize() const { return this->ncols * this->nrows; }

  /**
   * Returns the number of unused rows.
   *
   * @return number of unused rows
   */
  inline size_t getAdditionallyReservedRows() const {
    return this->size() / this->ncols - this->nrows;
  }

  /**
   * Determines the number of non-zero elements in the vector.
   *
   * @return The number of non-zero elements
   */
  size_t getNumberNonZero() const;

  /**
   * Returns the number of rows of the DataMatrix.
   *
   * @return Number of rows
   */
  inline size_t getNrows() const { return this->nrows; }

  /**
   * Returns the number of columns of the DataMatrix.
   *
   * @return Number of columns
   */
  inline size_t getNcols() const { return this->ncols; }

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
  void toString(std::string& text) const;

  void toFile(const std::string& fileName) const;

  /**
   * Returns a description of the DataMatrix as a string.
   *
   * @returns string of the DataMatrix
   */
  std::string toString() const;

 private:
  inline iterator row_begin(size_t row) { return this->begin() + row * this->ncols; }

  inline const_iterator row_begin(size_t row) const { return this->row_cbegin(row); }

  inline const_iterator row_cbegin(size_t row) const { return this->cbegin() + row * this->ncols; }

  inline iterator row_end(size_t row) { return this->row_begin(row + 1); }

  inline const_iterator row_end(size_t row) const { return this->row_cend(row); }

  inline const_iterator row_cend(size_t row) const { return this->row_cbegin(row + 1); }

  /// Number of rows of the data matrix
  size_t nrows;
  /// Number of columns of the data matrix
  size_t ncols;
};

}  // namespace base
}  // namespace sgpp

#endif /*DATAMATRIX_H_*/
