/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DistributedDataMatrix.hpp
 *
 * Created on: Jan 14, 2019
 *     Author: Jan Schopohl
 */
#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>
#include <sgpp/datadriven/scalapack/scalapack.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

class DistributedDataMatrix : sgpp::base::DataMatrix {
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
   * Creates an two-dimensional DistributedDataMatrix filled with a value.
   */
  DistributedDataMatrix(std::shared_ptr<BlacsProcessGrid> grid, int globalRows, int globalColumns,
                        DistributedDataMatrix::DTYPE dtype, int columnBlockSize, int rowBlockSize,
                        double value);

  /**
   * Creates an two-dimensional DistributedDataMatrix with the specified input
   */
  DistributedDataMatrix(double* input, std::shared_ptr<BlacsProcessGrid> grid, int globalRows,
                        int globalColumns, DistributedDataMatrix::DTYPE dtype, int columnBlockSize,
                        int rowBlockSize);

  void resize(int rows, int columns);

  /**
   * @returns The ScaLAPACK matrix descriptor
   */
  int* getDescriptor();

  /**
   * @returns number of rows assigned to the current process
   */
  int getLocalRows() const;

  /**
   * @returns number of columns assigned to the current process
   */
  int getLocalColumns() const;

 private:
  // The blacs process grid this matrix is distributed on
  std::shared_ptr<BlacsProcessGrid> grid;

  // total number of rows of the global matrix
  int globalRows;

  // total number of columns of the global matrix
  int globalColumns;

  // number of rows in this process
  int localRows;

  // number of columns in this process
  int localColumns;

  // ScaLAPACK matrix descriptor
  int descriptor[dlen_];
};

}  // namespace datadriven
}  // namespace sgpp