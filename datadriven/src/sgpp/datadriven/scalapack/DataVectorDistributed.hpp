// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataVector;

class DataVectorDistributed {
 public:
  /**
   * Creates a distributed data vector of specified size and initializes the elements to value.
   * @param grid blacs grid for distribution
   * @param globalSize global size (rows) of this vector
   * @param blockSize size for each block (one process might receive multiple blocks)
   * @param value Initial value of all elements, default 0.0
   */
  DataVectorDistributed(std::shared_ptr<BlacsProcessGrid> grid, size_t globalSize, size_t blockSize,
                        double value = 0.0);

  /**
   * Creates a distributed data vector with specified input data.
   * @param input pointer to input values for this vector
   * @param grid blacs grid for distribution
   * @param globalSize global size (rows) of this vector
   * @param blockSize size for each block (one process might receive multiple blocks)
   */
  DataVectorDistributed(double* input, std::shared_ptr<BlacsProcessGrid> grid, size_t globalSize,
                        size_t blockSize);

  /**
   * Returns the value of the element at position [row].
   * Only a const version exists as setting an element is not necessarily done in local memory.
   *
   * @param row Row
   * @return constant reference to the element
   */
  double operator()(size_t row) const;

  /**
   * Returns the value of the element at position [row]
   *
   * @param row Row
   * @return Value of the element
   */
  double get(size_t row) const;

  /**
   * Sets the element at global position [row] to value.
   *
   * @param row Global row
   * @param value New value for element
   */
  void set(size_t row, double value);

  /**
   * Set all entries of the vector to one value
   * @param value
   */
  void setAll(double value);

  /**
   * Copies all values of another distributed data matrix to this object, resizes this object to the
   * size of the other matrix.
   */
  void copyFrom(const DataVectorDistributed& other);

  /**
   * Adds another vector to this vector, modifies this vector.
   * @param x vector that is added to this vector
   */
  void add(const DataVectorDistributed& x);

  /**
   * Performs the following operation
   * sub(y) := sub(y) + a*sub(x)
   * @param y
   * @param x
   * @param a factor for x, default 1.0
   */
  static void add(DataVectorDistributed& y, const DataVectorDistributed& x, double a = 1.0);

  /**
   * @param y
   * @returns the dot product of this vector transposed and y
   */
  double dot(const DataVectorDistributed& y) const;

  /**
   * @param x
   * @param y
   * @returns dot = sub(x)'*sub(y); vector x is transposed
   */
  static double dot(const DataVectorDistributed& x, const DataVectorDistributed& y);

  /**
   * Scales this vector by factor a.
   * @param a
   */
  void scale(double a);

  /**
   * Resizes the vector to rows, discards the data.
   * @param rows
   */
  void resize(size_t rows);

  /**
   * @return pointer to the local data of this process
   */
  double* getLocalPointer();

  /**
   * @return const pointer to the local data of this process
   */
  const double* getLocalPointer() const;

  /**
   * @returns process grid used by the vector.
   */
  std::shared_ptr<BlacsProcessGrid> getProcessGrid() const;

  /**
   * @returns the gathered DataVector on process 0
   */
  DataVector toLocalDataVector() const;

  /**
   * @returns the broadcasted DataVector on all processes
   */
  DataVector toLocalDataVectorBroadcast() const;

  /**
   * @param[out] localVector the gathered DataVector on process (0, 0)
   */
  void toLocalDataVector(DataVector& localVector) const;

  /**
   * @param[out] localVector the broadcaster DataVector on all processes
   */
  void toLocalDataVectorBroadcast(DataVector& localVector) const;

  /**
   * @return pointer to the descriptor array for the underlying data matrix
   */
  int* getDescriptor();

  /**
   * @return const pointer to the descriptor array for the underlying data matrix
   */
  const int* getDescriptor() const;

  /**
   * Returns the number of rows of the DataMatrix.
   *
   * @return Number of rows
   */
  size_t getGlobalRows() const;

  /**
   * @returns number of rows assigned to the current process
   */
  size_t getLocalRows() const;

  /**
   * @returns the block size
   */
  size_t getBlockSize() const;

  /**
   * Prints the vector on stdout on process 0.
   */
  void printVector() const;

  /**
   * @returns true if part of the vector is mapped to the current process, false otherwise
   */
  bool isProcessMapped() const;

  /**
   * @returns reference to the underlying DataMatrixDistributed object.
   */
  DataMatrixDistributed& getMatrix();

  /**
   * @returns const ref to the underlying DataMatrixDistributed object.
   */
  const DataMatrixDistributed& getMatrix() const;

  /**
   * Distribute the input data to the process grid. Overwrites the current data.
   * For more information, see DataMatrixDistributed::distribute
   *
   * @param input input data
   * @param masterRow row of the process that distributes the data
   * @param masterCol col of the process that distributes the data
   */
  void distribute(double* input, int masterRow = 0, int masterCol = 0);

  /**
   * Calculates the local row index from the globalRowIndex.
   * @param globalRowIndex
   */
  size_t globalToLocalRowIndex(size_t globalRowIndex) const;

  /**
   * Calculates the global row index from the local row index.
   * @param localRowIndex
   */
  size_t localToGlobalRowIndex(size_t localRowIndex) const;

  /**
   * Calculates the row process index from the global row index.
   * @param globalRowIndex
   */
  size_t globalToRowProcessIndex(size_t globalRowIndex) const;

 private:
  // vector is mapped to a matrix with 1 column
  DataMatrixDistributed data;

  // shared pointer to the process grid used by the data matrix
  std::shared_ptr<BlacsProcessGrid> grid;
};

}  // namespace datadriven
}  // namespace sgpp
