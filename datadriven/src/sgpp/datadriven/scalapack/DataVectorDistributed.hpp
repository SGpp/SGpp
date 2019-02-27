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

#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>

namespace sgpp {
namespace datadriven {

class DataVectorDistributed {
 public:
  /**
   * Creates an distributed data vector of specified size and initializes the elements to value.
   * @param size Size of the vector
   * @param value Initial value of all elements, default 0.0
   */
  DataVectorDistributed(std::shared_ptr<BlacsProcessGrid> grid, size_t globalSize, size_t blockSize,
                        double value = 0.0);

  /**
   * @return pointer to the local data of this process
   */
  double* getLocalPointer();

  /**
   * @return const pointer to the local data of this process
   */
  const double* getLocalPointer() const;

  /**
   * @return pointer to the descriptor array for the underlying data matrix
   */
  int* getDescriptor();

  /**
   * @return const pointer to the descriptor array for the underlying data matrix
   */
  const int* getDescriptor() const;

 private:
  DataMatrixDistributed data;
};

}  // namespace datadriven
}  // namespace sgpp

#endif