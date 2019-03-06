/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DataVectorDistributed.hpp
 *
 * Created on: Feb 19, 2019
 *     Author: Jan Schopohl
 */
#ifdef USE_SCALAPACK

#pragma once

#include <memory>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>

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
   * @return pointer to the local data of this process
   */
  double* getLocalPointer();

  /**
   * @return const pointer to the local data of this process
   */
  const double* getLocalPointer() const;

  /**
   * @returns the gathered DataVector on process 0
   */
  DataVector toLocalDataVector() const;

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
   * Prints the vector on stdout on process 0.
   */
  void printVector() const;

  /**
   * @returns true if part of the vector is mapped to the current process, false otherwise
   */
  bool isProcessMapped() const;

 private:
  // vector is mapped to a matrix with 1 column
  DataMatrixDistributed data;

  // shared pointer to the process grid used by the data matrix
  std::shared_ptr<BlacsProcessGrid> grid;
};

}  // namespace datadriven
}  // namespace sgpp

#endif