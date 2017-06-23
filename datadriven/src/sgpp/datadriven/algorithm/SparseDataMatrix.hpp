/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SparseDataMatrix.hpp
 *
 *  Created on: Feb 4, 2017
 *      Author: Michael Lettrich
 */
#ifndef _WIN32

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataMatrix;

/**
 * Simple, rudimentary implementation of a sparse matrix structure based on the CRS (compressed row
 * storage format)
 */
class SparseDataMatrix {
 public:
  /**
   * Default constructor. Empty matrix.
   */
  SparseDataMatrix();

  /**
   * Construct sparse matrix from the contents of a dense matrix
   * @param mat Dense matrix  to convert to sparse matrix
   */
  explicit SparseDataMatrix(DataMatrix& mat);

  /**
   * Construct an empty sparse matrix with a certain size
   * @param nrows number of rows
   * @param ncols number of columns
   */
  SparseDataMatrix(size_t nrows, size_t ncols);

  /**
   * Construct an sparse matrix with a certain size and data
   * @param nrows number of rows
   * @param ncols number of columns
   * @param dataVector vector that holds the actual data values
   * @param colIndexVector vector holding the index of the last non zero entry in the row
   * @param rowPtrVector vector holding the row indices of the non zero entries
   */
  SparseDataMatrix(size_t nrows, size_t ncols, const std::vector<double>& dataVector,
                   const std::vector<size_t>& colIndexVector,
                   const std::vector<size_t>& rowPtrVector);

  /**
   * Construct an sparse matrix with a certain size and use move assignments to avoid copying
   * vectors.
   * @param nrows number of rows
   * @param ncols number of columns
   * @param dataVector vector that holds the actual data values
   * @param colIndexVector vector holding the index of the last non zero entry in the row
   * @param rowPtrVector vector holding the row indices of the non zero entries
   */
  SparseDataMatrix(size_t nrows, size_t ncols,
                   std::vector<double>&& dataVector,      // NOLINT(build/c++11)
                   std::vector<size_t>&& colIndexVector,  // NOLINT(build/c++11)
                   std::vector<size_t>&& rowPtrVector);   // NOLINT(build/c++11)

  /**
   * return the number of rows in the matrix
   * @return number of rows
   */
  size_t getNrows() const;

  /**
   * return the number of columns in the matrix
   * @return number of columns
   */
  size_t getNcols() const;

  /**
   * resize the matrix to new dimensions. Data is cropped accordingly if the new matrix is smaller,
   * scaling up is very cheap.
   * Cropping is very expensive, so use it carefully.
   * @param nrows new amount of rows
   * @param ncols new amount of columns
   */
  void resize(size_t nrows, size_t ncols);

  /**
   * Read only reference to the vector holding the index of the last non zero entry in each row
   * @return vector holding the index of the last non zero entry in each row
   */
  const std::vector<size_t>& getColIndexVector() const;

  /**
   *  Read only reference to the vector that holds the actual non zero values
   *  @return vector that holds the actual non zero values
   */
  const std::vector<double>& getDataVector() const;

  /**
   * Read only reference to the vector holding the row indices of the non zero entries
   * @return vector holding the row indices of the non zero entries
   */
  const std::vector<size_t>& getRowPtrVector() const;

  /**
   * reference to the vector holding the index of the last non zero entry in each row
   * @return vector holding the index of the last non zero entry in each row
   */
  std::vector<size_t>& getColIndexVector();

  /**
   *  Read only reference to the vector that holds the actual non zero values
   *  @return vector that holds the actual non zero values
   */
  std::vector<double>& getDataVector();

  /**
   * Read only reference to the vector holding the row indices of the non zero entries
   * @return vector holding the row indices of the non zero entries
   */
  std::vector<size_t>& getRowPtrVector();

  /**
   * static member function to convert a data matrix to a sparse data matrix.
   * @param in read only reference to the data matrix that should be converted
   * @param out the sparse data matrix that is filled with the non zero entries. Will be resized if
   * necessary.
   * @param threshold use a threshold
   */
  static void fromDataMatrix(const DataMatrix& in, SparseDataMatrix& out, double threshold = 0.0);

  /**
   * static member function to convert a data matrix to a sparse data matrix. only the lower
   * triangular part is read in.
   * @param in read only reference to the data matrix that should be converted
   * @param out the sparse data matrix that is filled with the non zero entries. Will be resized if
   * necessary.
   * @param threshold use a threshold
   */
  static void fromDataMatrixTriangular(const DataMatrix& in, SparseDataMatrix& out,
                                       double threshold = 0.0);

  /**
   * static member function to convert a sparse data matrix to a data matrix.
   * @param in read only reference to the sparse data matrix that should be converted
   * @param out the data matrix that is filled. Will be resized if necessary.
   */
  static void toDataMatrix(const SparseDataMatrix& in, DataMatrix& out);

 private:
  /*
   * number of rows in the matrix
   */
  size_t nrows;
  /**
   * number of columns in the matrix
   */
  size_t ncols;

  /**
   * vector holding the row indices of the non zero entries
   */
  std::vector<double> data;

  /**
   * vector holding the index of the last non zero entry in each row
   */
  std::vector<size_t> colIndex;

  /**
   * vector holding the row indices of the non zero entries
   */
  std::vector<size_t> rowPtr;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* _WIN_32 */
