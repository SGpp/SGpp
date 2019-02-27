/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DataMatrixDistributed.hpp
 *
 * Created on: Jan 14, 2019
 *     Author: Jan Schopohl
 */
#ifdef USE_SCALAPACK

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>
#include <sgpp/datadriven/scalapack/scalapack.hpp>

#include <memory>
#include <vector>

using sgpp::base::DataMatrix;

namespace sgpp {
namespace datadriven {

class DataVectorDistributed;

class DataMatrixDistributed {
 public:
  /**
   * Enum for the datatypes of scalapack matrices, see http://netlib.org/scalapack/slug/node73.html
   */
  enum class DTYPE : int {
    DENSE = 1,
    TRIDIAG_COEFFICIENT = 501,
    TRIDIAG_RHS = 502,
    OUT_OF_CORE = 602
  };

  /**
   * Creates an two-dimensional DataMatrixDistributed filled with a value.
   */
  DataMatrixDistributed(std::shared_ptr<BlacsProcessGrid> grid, int globalRows, int globalColumns,
                        DataMatrixDistributed::DTYPE dtype, int columnBlockSize, int rowBlockSize,
                        double value = 0.0);

  /**
   * Creates an two-dimensional DataMatrixDistributed with the specified input
   */
  DataMatrixDistributed(double* input, std::shared_ptr<BlacsProcessGrid> grid, int globalRows,
                        int globalColumns, DataMatrixDistributed::DTYPE dtype, int columnBlockSize,
                        int rowBlockSize);

  /**
   * Append nrows and ncols filled with zero
   * @param nrows
   * @param ncols
   */
  void appendZero(size_t nrows, size_t ncols);

  /**
   * Appends a new row with data contained in DataVector vec
   * and returns index of new row.
   *
   * @param vec DataVectorDistributed (length has to match getNcols()) with data
   * @return Index of new row
   */
  size_t appendRow(const DataVectorDistributed& vec);

  /**
   * Appends a new Col with data contained in DataVector vec
   * and returns index of new col.
   *
   * @param vec DataVectorDistributed (length has to match getNcols()) with data
   * @return Index of new col
   */
  size_t appendCol(const DataVectorDistributed& vec);

  /**
   * Transposes this matrix
   */
  void transpose();

  /**
   * Returns the value of the element at position [row,col].
   * Only a const version exists as setting an element is not necessarily done in local memory.
   *
   * @param row Row
   * @param col Column
   * @return constant reference to the element
   */
  const double& operator()(size_t row, size_t col) const;

  /**
   * Returns the value of the element at position [row,col]
   *
   * @param row Row
   * @param col Column
   * @return Value of the element
   */
  double get(size_t row, size_t col) const;

  /**
   * Sets the element at global position [row,col] to value.
   *
   * @param row Global row
   * @param col Global column
   * @param value New value for element
   */
  void set(size_t row, size_t col, double value);

  /**
   * Adds the values from another DataMatrix to the current values.
   * Modifies the current values.
   *
   * @param matr The DataMatrix which is added to the current values
   */
  void add(const DataMatrixDistributed& matr);

  /**
   * Subtracts the values from another DataMatrix of the current values.
   * Modifies the current values.
   *
   * @param matr The DataMatrix which is subtracted from the current values
   */
  void sub(const DataMatrixDistributed& matr);

  /**
   * expands a given DataVector into a
   * DataMatrix.
   *
   * @param expand DataVector that should be expanded
   */
  void expand(const DataVectorDistributed& expand);

  /**
   * Multiplies the matrix with a vector x and stores the result
   * in another vector y.
   * sub(y) := alpha*sub(A)'*sub(x) + beta*sub(y)
   *
   * @param[in] x vector to be multiplied
   * @param[out] y vector in which the result should be stored
   * @param[in] alpha factor alpha, default 1.0
   * @param[in] beta factor beta, default 0.0
   */
  void mult(const DataVectorDistributed& x, DataVectorDistributed& y, double alpha = 1.0,
            double beta = 0.0) const;

  /**
   * Uses BLAS to multiply this matrix (called A) with matrix B:
   *
   * sub(C) := alpha*op(sub(A))*op(sub(B)) + beta*sub(C)
   *
   * @param[in] transposeA Whether or not to transpose this matrix
   * @param[in] transposeB Whether or not to transpose matrix B
   * @param[in] b matrix B
   * @param[in, out] c used as output matrix C, can also be added to the product using factor beta
   * @param[in] alpha factor alpha, default 1.0
   * @param[in] beta factor beta, default 0.0
   */
  void mult(bool transposeA, bool transposeB, const DataMatrixDistributed& b,
            DataMatrixDistributed& c, double alpha = 1.0, double beta = 0.0) const;

  /**
   * @return Pointer to the local data of this process
   */
  double* getLocalPointer();

  /**
   * @return const pointer to the local data of this process
   */
  const double* getLocalPointer() const;

  /**
   * @returns the gathered DataMatrix as a normal, not distributed, DataMatrix.
   */
  DataMatrix toLocalDataMatrix() const;

  /**
   * @returns The ScaLAPACK matrix descriptor
   */
  int* getDescriptor();

  /**
   * @returns Const pointer to ScaLAPACK matrix descriptor
   */
  const int* getDescriptor() const;

  /**
   * Returns the number of rows of the DataMatrix.
   *
   * @return Number of rows
   */
  inline size_t getGlobalRows() const { return this->globalRows; }

  /**
   * Returns the number of columns of the DataMatrix.
   *
   * @return Number of columns
   */
  inline size_t getGlobalCols() const { return this->globalColumns; }

  /**
   * @returns number of rows assigned to the current process
   */
  int getLocalRows() const;

  /**
   * @returns number of columns assigned to the current process
   */
  int getLocalColumns() const;

  /**
   * @param processRow row of the process which should print the matrix
   * @param processColumn column of the process which should print the matrix
   */
  void printMatrix(size_t processRow = 0, size_t processCol = 0) const;

 private:
  size_t globalToLocalIndex(size_t globalIndex, size_t numberOfProcesses, size_t blockSize) const;

  size_t globalToProcessIndex(size_t globalIndex, size_t numberOfProcesses, size_t blockSize,
                              size_t processOffset) const;

  // vector to store the local data
  std::vector<double> localData;

  // The blacs process grid this matrix is distributed on
  std::shared_ptr<BlacsProcessGrid> grid;

  // total number of rows of the global matrix
  int globalRows;

  // total number of columns of the global matrix
  int globalColumns;

  // size of blocks in row dimension
  size_t rowBlockSize;

  // size of blocks in column dimension
  size_t columnBlockSize;

  // number of rows in this process
  int localRows;

  // number of columns in this process
  int localColumns;

  // ScaLAPACK matrix descriptor
  int descriptor[dlen_];
};

}  // namespace datadriven
}  // namespace sgpp

#endif  // USE_SCALAPACK