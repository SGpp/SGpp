// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DATAMATRIXSP_HPP
#define DATAMATRIXSP_HPP

#include <sgpp/base/datatypes/DataVectorSP.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

/**
 * A class to store two-dimensional data.
 * Typically, DataMatrixSP would contain a set of (d-dimensional) data or evaluation points, i.e.,
 * the DataMatrixSP consists of d columns, and each row is one of the points.
 * Thus, typical functionality like obtaining the maximum for a certain dimension (or attribute),
 * or normalizing all data points to the unit interval for a certain dimension are
 * provided.
 *
 * This is an re-implementation of the standard DataMatrix
 * for single precision floating point numbers in order to
 * increase support for GPUs.
 */
class DataMatrixSP : public std::vector<float> {
 public:
  /**
   * Creates an empty two-dimensional DataMatrixSP.
   */
  DataMatrixSP();

  /**
   * Copy constructor.
   */
  DataMatrixSP(const DataMatrixSP&) = default;

  /**
   * Move constructor
   */
  DataMatrixSP(DataMatrixSP&&) = default;

  /**
   * Copy assignment operator
   */
  DataMatrixSP& operator=(const DataMatrixSP&) = default;

  /**
   * Move assignment operator
   */
  DataMatrixSP& operator=(DataMatrixSP&&) = default;

  /**
   * Destructor
   */
  ~DataMatrixSP() = default;

  /**
   * Create a two-dimensional DataMatrixSP with @em nrows rows and
   * @em ncols columns (uninitialized values).
   *
   * @param nrows Number of rows
   * @param ncols Number of columns
   */
  DataMatrixSP(size_t nrows, size_t ncols);

  /**
   * Create a two-dimensional DataMatrixSP with @em nrows rows and
   * @em ncols columns and initializes all elements with the same value.
   *
   * @param nrows Number of rows
   * @param ncols Number of columns
   * @param value Value for all entries
   */
  DataMatrixSP(size_t nrows, size_t ncols, float value);

  /**
   * Create a new DataMatrixSP from a float array.
   * The float array contains the entries row-wise:
   * x0_0,x0_1,...,x0_ncol-1,
   * x1_0,x1_1,...
   * ...
   * xnrow_0, xnrow_1,...,xnrow_ncol-1
   *
   * @param input float array that contains the data
   * @param nrows number of rows
   * @param ncols number of columns
   */
  DataMatrixSP(const float* input, size_t nrows, size_t ncols);

  /**
   * Create a new DataMatrixSP from a std::vector<float>.
   *
   * @param input std::vector<float> that contains the data
   * @param nrows number of rows
   */
  explicit DataMatrixSP(std::vector<float> input, size_t nrows);

  /**
   * Create a new DataMatrixSP from a std::initializer_list<float>.
   *
   * @param input std::initializer_list<float> that contains the data
   * @param nrows number of rows
   */
  explicit DataMatrixSP(std::initializer_list<float> input, size_t nrows);

  static DataMatrixSP fromFile(const std::string& fileName);

  static DataMatrixSP fromString(const std::string& serializedVector);

  /**
   * Resizes the DataMatrixSP to nrows rows.
   * All new additional entries are uninitialized.
   * If nrows is smaller than the current number of rows,
   * all superfluous entries are removed.
   * \deprecated use resizeRows
   *
   * @param nrows New number of rows of the DataMatrixSP
   */
  void resize(size_t nrows);

  /**
   * Resizes the DataMatrixSP to nrows rows.
   * All new additional entries are uninitialized.
   * If nrows is smaller than the current number of rows,
   * all superfluous entries are removed.
   *
   * @param nrows New number of rows of the DataMatrixSP
   */
  void resizeRows(size_t nrows);

  /**
   * Resizes the DataMatrixSP to nrows rows and ncols columns.
   * All new additional entries are uninitialized.
   * If nrows*ncols is smaller than the current number of entries,
   * all superfluous entries are removed.
   * \deprecated use resizeRowsCols
   *
   * @param nrows New number of rows of the DataMatrixSP
   * @param ncols New number of columns of the DataMatrixSP
   */
  void resize(size_t nrows, size_t ncols);

  /**
   * Resizes the DataMatrixSP to nrows rows and ncols columns.
   * All new additional entries are uninitialized.
   * If nrows*ncols is smaller than the current number of entries,
   * all superfluous entries are removed.
   *
   * @param nrows New number of rows of the DataMatrixSP
   * @param ncols New number of columns of the DataMatrixSP
   */
  void resizeRowsCols(size_t nrows, size_t ncols);

  /**
   * Resizes the quadratic DataMatrixSP to size rows and size columns.
   * All new additional entries are uninitialized.
   * If size is smaller than the current size of the quadratic DataMatrixSP,
   * all superfluous entries are removed.
   *
   * @param size New dimension of quadratic data DataMatrixSP
   */
  void resizeQuadratic(size_t size);

  /**
   * Resizes the DataMatrixSP to nrows rows.
   * All new additional entries are set to zero.
   * If nrows is smaller than the current number of rows,
   * all superfluous entries are removed.
   * \deprecated use resizeRows
   *
   * @param nrows New number of rows of the DataMatrixSP
   */
  void resizeZero(size_t nrows);

  /**
   * Resizes the DataMatrixSP to nrows rows and ncols columns.
   * All new additional entries are set to zero.
   * If nrows*ncols is smaller than the current number of entries,
   * all superfluous entries are removed.
   * \deprecated use resizeRowsCols
   *
   * @param nrows New number of rows of the DataMatrixSP
   * @param ncols New number of columns of the DataMatrixSP
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
   * Appends a new row with data contained in DataVectorSP vec
   * and returns index of new row.
   *
   * @param vec DataVectorSP (length has to match getNcols()) with data
   * @return Index of new row
   */
  size_t appendRow(const DataVectorSP& vec);

  /**
   * Appends a new Col with data contained in DataVectorSP vec
   * and returns index of new col.
   *
   * @param vec DataVectorSP (length has to match getNcols()) with data
   * @return Index of new col
   */
  size_t appendCol(const DataVectorSP& vec);

  /**
   * Sets all entries of DataMatrixSP to value.
   *
   * @param value New value for all entries
   */
  void setAll(float value);

  /**
   * Copies the data from another DataMatrixSP matr.
   * Disregards the number of rows and columns set for the two matrices,
   * i.e., just copies the data entry by entry (and row by row).
   * If the dimensions match (nrows, ncols), the current DataMatrixSP is an
   * exact copy of matr. If not, as many elements as possible are
   * copied, and everything else is left untouched.
   *
   * @param matr The source DataMatrixSP containing the data
   */
  void copyFrom(const DataMatrixSP& matr);

  /**
   * Transposes this DataMatrixSP
   */
  void transpose();

  /**
   * Returns the value of the element at position [row,col]
   *
   * @param row Row
   * @param col Column
   * @return reference to the element
   */
  inline float& operator()(size_t row, size_t col) { return (*this)[row * ncols + col]; }

  /**
   * Returns the value of the element at position [row,col]
   *
   * @param row Row
   * @param col Column
   * @return constant reference to the element
   */
  inline const float& operator()(size_t row, size_t col) const {
    return (*this)[row * ncols + col];
  }

  /**
   * Returns the value of the element at position [row,col]
   *
   * @param row Row
   * @param col Column
   * @return Value of the element
   */
  inline float get(size_t row, size_t col) const { return (*this)[row * ncols + col]; }

  /**
   * Sets the element at position [row,col] to value.
   *
   * @param row Row
   * @param col Column
   * @param value New value for element
   */
  inline void set(size_t row, size_t col, float value) { (*this)[row * ncols + col] = value; }

  /**
   * Copies the values of a row to the DataVectorSP vec.
   *
   * @param row The row
   * @param vec DataVectorSP into which the data is written
   */
  void getRow(size_t row, DataVectorSP& vec) const;

  /**
   * Copies the values of a row to the std::vector vec.
   *
   * @param row The row
   * @param vec std::vector into which the data is written
   */
  void getRow(size_t row, std::vector<float>& vec) const;

  /**
   * Sets a row of the DataMatrixSP to the values of a DataVectorSP vec.
   *
   * @param row The row which is to be overwritten
   * @param vec DataVectorSP containing the data of the row
   */
  void setRow(size_t row, const DataVectorSP& vec);

  /**
   * Copies the values of a column to the DataVectorSP vec.
   *
   * @param col The column
   * @param vec DataVectorSP into which the data is written
   */
  void getColumn(size_t col, DataVectorSP& vec) const;

  /**
   * Sets a column of the DataMatrixSP to the values of a DataVectorSP vec.
   *
   * @param col The column which is to be overwritten
   * @param vec DataVectorSP containing the data of the column
   */
  void setColumn(size_t col, const DataVectorSP& vec);

  /**
   * Adds the values from another DataMatrixSP to the current values.
   * Modifies the current values.
   *
   * @param matr The DataMatrixSP which is added to the current values
   */
  void add(const DataMatrixSP& matr);

  /**
   * Subtracts the values from another DataMatrixSP of the current values.
   * Modifies the current values.
   *
   * @param matr The DataMatrixSP which is subtracted from the current values
   */
  void sub(const DataMatrixSP& matr);

  /**
   * Reduce the DataMatrixSP along the
   * columns by adding all entries in one row.
   *
   * @param reduction DataVectorSP into which the reduce columns are stored
   */
  void addReduce(DataVectorSP& reduction);

  /**
   * Reduce the DataMatrixSP along the
   * columns by adding all entries in one row.
   *
   * @param reduction DataVectorSP to which the reduce columns are added
   * @param beta vector with length of number of columns beta[i] is multiplied to each element
   *             row[j][i]
   * @param start_beta where to start using the beta coefficients
   */
  void addReduce(DataVectorSP& reduction, DataVectorSP& beta, size_t start_beta);

  /**
   * expands a given DataVectorSP into a
   * DataMatrixSP.
   *
   * @param expand DataVectorSP that should be expanded
   */
  void expand(const DataVectorSP& expand);

  /**
   * Multiplies the current DataMatrixSP component-wise with another DataMatrixSP.
   * Modifies the current values.
   * Performs
   * @code
   * for i from 1 to this.getSize()
   *   this[i] *= matr[i]
   * @endcode
   *
   * @param matr the DataMatrixSP which is multiplied to current DataMatrixSP
   */
  void componentwise_mult(const DataMatrixSP& matr);

  /**
   * Divides the current DataMatrixSP component-wise by another DataMatrixSP.
   * Modifies the current values.
   * Performs
   * @code
   * for i from 1 to this.getTotalSize()
   *   this[i] /= matr[i]
   * @endcode
   * Note: <b>No check for division by zero!</b>
   *
   * @param matr the DataMatrixSP which the current DataMatrixSP is divided by
   */
  void componentwise_div(const DataMatrixSP& matr);

  /**
   * Multiplies all elements by a constant factor
   *
   * @param scalar the constant
   */
  void mult(float scalar);

  /**
   * Multiplies the matrix with a vector x and stores the result
   * in another vector y.
   *
   * @param[in] x vector to be multiplied
   * @param[out] y vector in which the result should be stored
   */
  void mult(const DataVectorSP& x, DataVectorSP& y) const;

  /**
   * Squares all elements of the DataMatrixSP
   */
  void sqr();

  /**
   * Takes the square root of all elements of the DataMatrixSP
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
  float sum() const;

  /**
   * Returns the minimum value of column col.
   *
   * @param col Number of the column
   *
   * @return Minimum value
   */
  float min(size_t col) const;

  /**
   * Returns the minimum over all entries.
   *
   * @return Minimal value of all entries
   */
  float min() const;

  /**
   * Returns the maximum value of column col.
   *
   * @param col Number of the column
   *
   * @return Maximum value
   */
  float max(size_t col) const;

  /**
   * Returns the maximum over all entries.
   *
   * @return Maximal value of all entries
   */
  float max() const;

  /**
   * Determines minimum and maximum of column col.
   *
   * @param col Number of the column
   * @param min Reference variable for the minimum
   * @param max Reference variable for the maximum
   */
  void minmax(size_t col, float* min, float* max) const;

  /**
   * Determines minimum and maximum over all entries.
   *
   * @param min Reference variable for the minimum
   * @param max Reference variable for the maximum
   */
  void minmax(float* min, float* max) const;

  /**
   * Returns pointer to float array containing underlying data.
   *
   * @return Pointer to data
   */
  float* getPointer();

  /**
   * Returns const pointer to float array containing underlying data.
   *
   * @return Const pointer to data
   */
  const float* getPointer() const;

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
   * Returns the number of rows of the DataMatrixSP.
   *
   * @return Number of rows
   */
  inline size_t getNrows() const { return this->nrows; }

  /**
   * Returns the number of columns of the DataMatrixSP.
   *
   * @return Number of columns
   */
  inline size_t getNcols() const { return this->ncols; }

  /**
   * Normalizes the d-th dimension (entries in the d-th column) to @f$[0,1]@f$.
   * Considers contents of DataMatrixSP as a d-dimensional dataset, one
   * data point per row.
   *
   * @param d The dimension (column) that should be normalized (starting with 0)
   */
  void normalizeDimension(size_t d);

  /**
   * Normalizes the d-th dimension (entries in the d-th column) to @f$[border,1-border]@f$.
   * Considers contents of DataMatrixSP as a d-dimensional dataset, one
   * data point per row.
   *
   * @param d The dimension (column) that should be normalized (starting with 0)
   * @param border Width of the border
   */
  void normalizeDimension(size_t d, float border);

  /**
   * Writes the data stored in the DataMatrixSP into a string
   *
   * @param text String to which the data is written
   */
  void toString(std::string& text) const;

  void toFile(const std::string& fileName) const;

  /**
   * Returns a description of the DataMatrixSP as a string.
   *
   * @returns string of the DataMatrixSP
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

#endif /* DATAMATRIXSP_HPP */
