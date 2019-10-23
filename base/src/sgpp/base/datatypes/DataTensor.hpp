// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * A class to store three-dimensional data.
 * It can be thought of as a vector of matrices
 */
class DataTensor : public std::vector<double> {
 public:
  /**
   * Creates an empty two-dimensional DataTensor.
   */
  DataTensor();

  //   /**
  //    * Copy constructor.
  //    */
  //   DataTensor(const DataTensor&) = default;

  //   // /**
  //   //  * Move constructor
  //   //  */
  //   // DataTensor(DataTensor&&) = default;

  //   /**
  //    * Copy assignment operator
  //    */
  //   DataTensor& operator=(const DataTensor&) = default;

  //   // /**
  //   //  * Move assignment operator
  //   //  */
  //   // DataTensor& operator=(DataTensor&&) = default;

  /**
   * Destructor
   */
  ~DataTensor() = default;

  /**
   * Create a thre-dimensional DataTensor of depth @em ndepth with
   * @em nrows rows and @em ncols columns (uninitialized values).
   *
   * @param ndepth  depth of the tensor
   * @param nrows   Number of rows of each entry matrix
   * @param ncols   Number of columns of each entry matrix
   */
  DataTensor(size_t ndepth, size_t nrows, size_t ncols);

  /**
   * Create a three-dimensional DataTensor of depth @em ndepth with
   * @em nrows rows and @em ncols columns and initializes all elements
   * with the same value.
   *
   * @param ndepth  depth of the tensor
   * @param nrows   Number of rows
   * @param ncols   Number of columns
   * @param value   Value for all entries
   */
  DataTensor(size_t ndepth, size_t nrows, size_t ncols, double value);

  /**
   * Resizes the DataTensor to depth ndepth, nrows rows and ncols columns.
   * All new additional entries are set to zero.
   * If nrows*ncols is smaller than the current number of entries,
   * all superfluous entries are removed.
   * \deprecated use resizeRowsCols
   *
   * @param ndepth New depth of the DataTensor
   * @param nrows New number of rows of the DataTensor
   * @param ncols New number of columns of the DataTensor
   */
  void resize(size_t ndepth, size_t nrows, size_t ncols);

  /**
   * Returns the value of the element at position [depth,row,col]
   *
   * @param depth   Depth
   * @param row Row
   * @param col Column
   * @return Value of the element
   */
  double get(size_t depth, size_t row, size_t col);

  /**
   * Sets the element at position [depth,row,col] to value.
   *
   * @param depth Depth
   * @param row Row
   * @param col Column
   * @param value New value for element
   */
  void set(size_t depth, size_t row, size_t col, double value);

  /**
   * Returns the depth of the DataTensor.
   *
   * @return depth
   */
  inline size_t getNdepth() const { return this->ndepth; }

  /**
   * Returns the number of rows of the DataTensor.
   *
   * @return Number of rows
   */
  inline size_t getNrows() const { return this->nrows; }

  /**
   * Returns the number of columns of the DataTensor.
   *
   * @return Number of columns
   */
  inline size_t getNcols() const { return this->ncols; }

  /**
   * Retruns one of the matrix entries of the tensor
   *
   * @param depth   The depth
   * @param matrix  DataMatrix in which the data is written
   */
  void getMatrix(size_t depth, sgpp::base::DataMatrix& matrix);

  /**
   * Returns one row of one of the matrix entries of the tensor
   *
   * @param depth     Depth of the matrix entry
   * @param row       row of the matrix entry
   * @param vec       DataVector in which the data is written
   */
  void getRow(size_t depth, size_t row, sgpp::base::DataVector& vec);

  /**
   * Returns one column of one of the matrix entries of the tensor
   *
   * @param depth     Depth of the matrix entry
   * @param col       column of the matrix entry
   * @param vec       DataVector in which the data is written
   */
  void getColumn(size_t depth, size_t col, sgpp::base::DataVector& vec);

  /**
   * Writes the data stored in the DataTensor into a string
   *
   * @param text String to which the data is written
   */
  void toString(std::string& text) const;

  /**
   * Returns a description of the DataTensor as a string.
   *
   * @returns string of the DataTensor
   */
  std::string toString() const;

 private:
  /// Depth of the tensor
  size_t ndepth;
  /// Number of rows of each entry matrix
  size_t nrows;
  /// Number of columns of each entry matrix
  size_t ncols;
  /// number of entries of each entry matrix
  size_t matrixsize;
  /// the data
  std::vector<DataMatrix> data;
};

}  // namespace base
}  // namespace sgpp
