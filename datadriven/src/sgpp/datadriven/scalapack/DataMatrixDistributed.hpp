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

/**
 * Class to represent a DataMatrix which is distributed on a process grid.
 * The class provides a wrapper for ScaLAPACK methods on the matrix.
 */
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
                        int columnBlockSize, int rowBlockSize, double value = 0.0,
                        DTYPE dtype = DTYPE::DENSE);

  /**
   * Creates an two-dimensional DataMatrixDistributed with the specified input.
   * Always call this method from *all* processes, as a BLACS grid is created inside. The init
   * method of a BLACS grid has to be called from all processes, otherwise a deadlock will occur.
   */
  DataMatrixDistributed(double* input, std::shared_ptr<BlacsProcessGrid> grid, int globalRows,
                        int globalColumns, int columnBlockSize, int rowBlockSize,
                        DTYPE dtype = DTYPE::DENSE);

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
  double operator()(size_t row, size_t col) const;

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
   * @param a The DataMatrix which is added to the current values
   */
  void add(const DataMatrixDistributed& a);

  /**
   * Calculates sub(C):=beta*sub(C) - alpha*op(sub(A))
   *
   * @param c matrix C
   * @param a matrix A
   * @param transposeA transpose matrix A if true, default false
   * @param beta scalar factor for matrix C, default 1.0
   * @param alpha scalar factor for matrix A, default 1.0
   */
  static void add(DataMatrixDistributed& c, const DataMatrixDistributed& a, bool transposeA = false,
                  double beta = 1.0, double alpha = 1.0);

  /**
   * Subtracts the values from another DataMatrix of the current values.
   * Modifies the current values.
   *
   * @param a The DataMatrix which is subtracted from the current values
   */
  void sub(const DataMatrixDistributed& a);

  /**
   * Calculates sub(C):=beta*sub(C) - alpha*op(sub(A))
   *
   * @param c matrix C
   * @param a matrix A
   * @param transposeA transpose matrix A if true, default false
   * @param beta scalar factor for matrix C, default 1.0
   * @param alpha scalar factor for matrix A, default 1.0
   */
  static void sub(DataMatrixDistributed& c, const DataMatrixDistributed& a, bool transposeA = false,
                  double beta = 1.0, double alpha = 1.0);

  /**
   * Multiplies the matrix with a vector x and stores the result
   * in another vector y.
   * sub(y) := alpha*sub(this)'*sub(x) + beta*sub(y)
   *
   * @param[in] x vector to be multiplied
   * @param[in, out] y vector in which the result should be stored
   * @param[in] transpose transpose if true, default false
   * @param[in] alpha factor alpha, default 1.0
   * @param[in] beta factor beta, default 0.0
   */
  void mult(const DataVectorDistributed& x, DataVectorDistributed& y, bool transpose = false,
            double alpha = 1.0, double beta = 0.0) const;

  /**
   * Multiplies matrix A with a vector x and stores the result
   * in another vector y.
   * sub(y) := alpha*sub(A)'*sub(x) + beta*sub(y)
   *
   * @param[in] a matrix to be multiplied
   * @param[in] x vector to be multiplied
   * @param[in, out] y vector in which the result should be stored
   * @param[in] transpose transpose if true, default false
   * @param[in] alpha factor alpha, default 1.0
   * @param[in] beta factor beta, default 0.0
   */
  static void mult(const DataMatrixDistributed& a, const DataVectorDistributed& x,
                   DataVectorDistributed& y, bool transpose = false, double alpha = 1.0,
                   double beta = 0.0);

  /**
   * Uses PBLAS to multiply this matrix with matrix B:
   *
   * sub(C) := alpha*op(sub(this))*op(sub(B)) + beta*sub(C)
   *
   * @param[in] b matrix B
   * @param[in, out] c used as output matrix C, can also be added to the product using factor
   * beta
   * @param[in] transposeA Whether or not to transpose this matrix
   * @param[in] transposeB Whether or not to transpose matrix B
   * @param[in] alpha factor alpha, default 1.0
   * @param[in] beta factor beta, default 0.0
   */
  void mult(const DataMatrixDistributed& b, DataMatrixDistributed& c, bool transposeA = false,
            bool transposeB = false, double alpha = 1.0, double beta = 0.0) const;

  /**
   * Uses BLAS to multiply matrix A with matrix B:
   *
   * sub(C) := alpha*op(sub(A))*op(sub(B)) + beta*sub(C)
   *
   * @param[in] a matrix A
   * @param[in] b matrix B
   * @param[in, out] c used as output matrix C, can also be added to the product using factor
   * beta
   * @param[in] transposeA Whether or not to transpose this matrix
   * @param[in] transposeB Whether or not to transpose matrix B
   * @param[in] alpha factor alpha, default 1.0
   * @param[in] beta factor beta, default 0.0
   */
  static void mult(const DataMatrixDistributed& a, const DataMatrixDistributed& b,
                   DataMatrixDistributed& c, bool transposeA = false, bool transposeB = false,
                   double alpha = 1.0, double beta = 0.0);

  /**
   * Convert this matrix to a vector. This operation is only possible if either globalRows or
   * globalColumns equals one.
   *
   * @returns A DataVectorDistributed object with the same size and data
   */
  DataVectorDistributed toVector();

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
  size_t getGlobalRows() const;

  /**
   * Returns the number of columns of the DataMatrix.
   *
   * @return Number of columns
   */
  size_t getGlobalCols() const;

  /**
   * @returns number of rows assigned to the current process
   */
  int getLocalRows() const;

  /**
   * @returns number of columns assigned to the current process
   */
  int getLocalColumns() const;

  /**
   * Prints the matrix on stdout on process 0.
   */
  void printMatrix() const;

  /**
   * @returns true if part of the matrix is mapped to the current process, false otherwise
   */
  bool isProcessMapped() const;

 private:
  /**
   * Distribute the matrix on the process grid according to the 2d block cyclic scheme used in
   * scalapack. More information: http://www.netlib.org/scalapack/slug/node76.html
   *
   * @param matrix Pointer to the local matrix, only relevant for the master process
   * @param masterRow row coordinate of the master process, default 0
   * @param masterCol col coordinate of the master process, default 0
   */
  void distribute(const double* matrix, int masterRow = 0, int masterCol = 0);

  /**
   * Gather the distributed matrix into one local matrix.
   *
   * @param masterRow row coordinate of the master process, default 0
   * @param masterCol col coordinate of the master process, default 0
   */
  DataMatrix gather(int masterRow = 0, int masterCol = 0) const;

  size_t globalToLocalIndex(size_t globalIndex, size_t numberOfProcesses, size_t blockSize) const;

  int globalToProcessIndex(size_t globalIndex, size_t numberOfProcesses, size_t blockSize,
                           int processOffset) const;

  // vector to store the local data
  std::vector<double> localData;

  // The blacs process grid this matrix is distributed on
  std::shared_ptr<BlacsProcessGrid> grid;

  // total number of rows of the global matrix
  int globalRows;

  // total number of columns of the global matrix
  int globalColumns;

  // size of blocks in row dimension
  int rowBlockSize;

  // size of blocks in column dimension
  int columnBlockSize;

  // number of rows in this process
  int localRows;

  // number of columns in this process
  int localColumns;

  // leading dimension
  int leadingDimension;

  // ScaLAPACK matrix descriptor
  int descriptor[dlen_];
};

}  // namespace datadriven
}  // namespace sgpp

#endif  // USE_SCALAPACK