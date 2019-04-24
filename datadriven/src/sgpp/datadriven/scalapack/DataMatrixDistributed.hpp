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
 * See http://netlib.org/scalapack/slug/node76.html for information regarding the distribution
 * scheme.
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
   * Enum that specifies the lower or upper triangular part of a matrix.
   */
  enum class TRIANGULAR { LOWER, UPPER };

  /**
   * Creates an empty DataMatrixDistributed object. Warning: This object cannot be used.
   */
  DataMatrixDistributed();

  /**
   * Creates an two-dimensional DataMatrixDistributed filled with a value.
   * @param grid BLACS process grid this matrix will be distributed on
   * @param globalRows number of global rows of the matrix
   * @param globalColumns number of global columns of the matrix
   * @param rowBlockSize block size in the row dimension
   * @param columnBlockSize block size in the column dimension
   * @param value initial value for all elements of the matrix, default 0
   * @param dtype datatype of the matrix, default DENSE
   */
  DataMatrixDistributed(std::shared_ptr<BlacsProcessGrid> grid, size_t globalRows,
                        size_t globalColumns, size_t rowBlockSize, size_t columnBlockSize,
                        double value = 0.0, DTYPE dtype = DTYPE::DENSE);

  /**
   * Creates an two-dimensional DataMatrixDistributed with the specified input. Call this matrix on
   * all processes in the grid to ensure proper initialization.
   * @param input pointer to data that will be distributed
   * @param grid BLACS process grid this matrix will be distributed on
   * @param globalRows number of global rows of the matrix
   * @param globalColumns number of global columns of the matrix
   * @param rowBlockSize block size in the row dimension
   * @param columnBlockSize block size in the column dimension
   * @param dtype datatype of the matrix, default DENSE
   */
  DataMatrixDistributed(double* input, std::shared_ptr<BlacsProcessGrid> grid, size_t globalRows,
                        size_t globalColumns, size_t rowBlockSize, size_t columnBlockSize,
                        DTYPE dtype = DTYPE::DENSE);

  /**
   * Creates a distributed data matrix from data which is already shared (mirrored) on each process.
   * Avoids network transfers.
   * @param input pointer to data that will be distributed
   * @param grid BLACS process grid this matrix will be distributed on
   * @param globalRows number of global rows of the matrix
   * @param globalColumns number of global columns of the matrix
   * @param rowBlockSize block size in the row dimension
   * @param columnBlockSize block size in the column dimension
   * @param dtype datatype of the matrix, default DENSE
   */
  static DataMatrixDistributed fromSharedData(const double* input,
                                              std::shared_ptr<BlacsProcessGrid> grid,
                                              size_t globalRows, size_t globalColumns,
                                              size_t rowBlockSize, size_t columnBlockSize,
                                              DTYPE dtype = DTYPE::DENSE);

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
   * Set all entries of the matrix to one value
   * @param value
   */
  void setAll(double value);

  /**
   * Copies all values of another distributed data matrix to this object, resizes this object to the
   * size of the other matrix.
   */
  void copyFrom(const DataMatrixDistributed& other);

  /**
   * Transposes this matrix.
   * @returns the transposed version of this matrix
   */
  DataMatrixDistributed transpose();

  /**
   * Transposes matrix A and stores the result in C:
   * sub(C):=beta*sub(C) + alpha*sub(A)'
   *
   * @param[in] a matrix A to transpose
   * @param[out] c result matrix C
   * @param[in] alpha factor for matrix A, default 1.0
   * @param[in] beta factor for matrix C, default 0.0
   */
  static void transpose(const DataMatrixDistributed& a, DataMatrixDistributed& c,
                        double alpha = 1.0, double beta = 0.0);

  /**
   * Adds the values from another DataMatrix to the current values.
   * Modifies the current values.
   *
   * @param a The DataMatrix which is added to the current values
   */
  void add(const DataMatrixDistributed& a);

  /**
   * Calculates sub(C):=beta*sub(C) + alpha*op(sub(A))
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
   * Solves a linear system of equations Ax=b using a previously computed Cholesky decomposition
   * A=LL^T
   *
   * @param[in] l lower triangular matrix L of the Cholesky decomposition
   * @param[in, out] input vector b of the linear system, is overwritten with solution x
   */
  static void solveCholesky(const DataMatrixDistributed& l, DataVectorDistributed& b,
                            TRIANGULAR uplo = TRIANGULAR::LOWER);

  /**
   * Resizes the matrix to rows and cols, data is discarded.
   * @param rows
   * @param cols
   */
  void resize(size_t rows, size_t cols);

  /**
   * Convert this matrix to a vector. This operation is only possible if either globalRows or
   * globalColumns equals one.
   *
   * @returns A DataVectorDistributed object with the same size and data
   */
  DataVectorDistributed toVector() const;

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
   * @returns the whole DataMatrix on broadcasted to all processes in the grid
   */
  DataMatrix toLocalDataMatrixBroadcast() const;

  /**
   * @param[out] the gathered DataMatrix as a normal, not distributed, DataMatrix. Result can only
   * be used on the master process.
   */
  void toLocalDataMatrix(DataMatrix& localMatrix) const;

  /**
   * @param[out] the whole DataMatrix is broadcasted to all processes in the grid
   */
  void toLocalDataMatrixBroadcast(DataMatrix& localMatrix) const;

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
  size_t getLocalRows() const;

  /**
   * @returns number of columns assigned to the current process
   */
  size_t getLocalColumns() const;

  /**
   * @returns the row block size.
   */
  size_t getRowBlockSize() const;

  /**
   * @returns the column block size.
   */
  size_t getColumnBlockSize() const;

  /**
   * @returns the process grid of this matrix.
   */
  std::shared_ptr<BlacsProcessGrid> getProcessGrid() const;

  /**
   * Prints the matrix on stdout on process 0.
   */
  void printMatrix() const;

  /**
   * @returns true if part of the matrix is mapped to the current process, false otherwise
   */
  bool isProcessMapped() const;

  /**
   * Calculates the local row index from the globalRowIndex.
   * @param globalRowIndex
   */
  size_t globalToLocalRowIndex(size_t globalRowIndex) const;

  /**
   * Calculates the local column index from the globalColumnIndex.
   * @param globalColumn
   */
  size_t globalToLocalColumnIndex(size_t globalColumnIndex) const;

  /**
   * Calculates the global row index from the local row index.
   * @param localRowIndex
   */
  size_t localToGlobalRowIndex(size_t localRowIndex) const;

  /**
   * Calculates the global column index from the local column index.
   * @param localColumnIndex
   */
  size_t localToGlobalColumnIndex(size_t localColumnIndex) const;

  /**
   * Calculates the row process index from the global row index.
   * @param globalRowIndex
   */
  size_t globalToRowProcessIndex(size_t globalRowIndex) const;

  /**
   * Calculates the column process index from the global column index.
   * @param globalColumnIndex
   */
  size_t globalToColumnProcessIndex(size_t globalColumnIndex) const;

  /**
   * Distribute the matrix from the master on the process grid according to the 2d block cyclic
   * scheme used in scalapack. More information: http://www.netlib.org/scalapack/slug/node76.html
   *
   * @param matrix Pointer to the local matrix, only relevant for the master process
   * @param masterRow row coordinate of the master process, default 0
   * @param masterCol col coordinate of the master process, default 0
   */
  void distribute(const double* matrix, int masterRow = 0, int masterCol = 0);

 private:
  /**
   * Gather the distributed matrix into one local matrix.
   *
   * @param masterRow row coordinate of the master process, default 0
   * @param masterCol col coordinate of the master process, default 0
   */
  DataMatrix gather(int masterRow = 0, int masterCol = 0) const;

  /**
   * Gather the distributed matrix into one local matrix.
   *
   * @param[out] localMatrix result
   * @param masterRow row coordinate of the master process, default 0
   * @param masterCol col coordinate of the master process, default 0
   */
  void gather(DataMatrix& localMatrix, int masterRow = 0, int masterCol = 0) const;

  /**
   * Broadcasts the whole matrix to all processes in the grid.
   */
  DataMatrix broadcast() const;

  /**
   * Broadcasts the whole matrix to all processes in the grid.
   *
   * @param[out] localMatrix result
   */
  void broadcast(DataMatrix& localMatrix) const;

  /**
   * Calculates the row or column index of the element in the local process.
   * @param globalIndex global index of the element
   * @param numberOfProcesses number of processes in the grid in a row or column
   * @param blockSize size of one block
   */
  size_t globalToLocalIndex(size_t globalIndex, size_t numberOfProcesses, size_t blockSize) const;

  /**
   * Calculates the global row or column index of the element.
   * @param localIndex local index of the element
   * @param process process of the element
   * @param numberOfProcesses number of processes in the grid in a row or column
   * @param blockSize size of one block
   */
  size_t localToGlobalIndex(size_t localIndex, size_t process, size_t numberOfProcesses,
                            size_t blockSize) const;

  /**
   * Calculates the row or column index of the process of an element from its global index
   * @param globalIndex global index of the element
   * @param numberOfProcesses number of processes in the grid in a row or column
   * @param blockSize size of one block
   * @param processOffset process offset in the blacs grid, default 0
   */
  int globalToProcessIndex(size_t globalIndex, size_t numberOfProcesses, size_t blockSize,
                           int processOffset = 0) const;

  // vector to store the local data
  std::vector<double> localData;

  // The blacs process grid this matrix is distributed on
  std::shared_ptr<BlacsProcessGrid> grid;

  // total number of rows of the global matrix
  size_t globalRows;

  // total number of columns of the global matrix
  size_t globalColumns;

  // size of blocks in row dimension
  size_t rowBlockSize;

  // size of blocks in column dimension
  size_t columnBlockSize;

  // number of rows in this process
  size_t localRows;

  // number of columns in this process
  size_t localColumns;

  // leading dimension
  size_t leadingDimension;

  // ScaLAPACK matrix descriptor
  int descriptor[dlen_];
};

}  // namespace datadriven
}  // namespace sgpp
