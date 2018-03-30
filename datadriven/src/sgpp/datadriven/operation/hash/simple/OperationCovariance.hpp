// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * This class provides the covariance matrix a sparse grid function
 */
class OperationCovariance {
 public:
  /**
   * Constructor
   * @param grid grid
   */
  explicit OperationCovariance(base::Grid& grid) : grid(grid) {}

  /**
   * Destructor
   */
  virtual ~OperationCovariance() {}

  /**
   * Integrate the sparse grid function
   *
   * @param alpha the function's values in the nodal basis
   * @param cov where the covariance matrix will be stored
   * @param bounds describes the boundaries of the hypercube of the original function
   */
  virtual void doQuadrature(base::DataVector& alpha, base::DataMatrix& cov,
                            base::DataMatrix* bounds = nullptr);

 private:
  base::DataMatrix* loadBounds(size_t numDims, base::DataMatrix* bounds, size_t idim,
                               size_t jdim = 0);
  double mean(base::Grid& grid, base::DataVector& alpha, base::DataMatrix* bounds = nullptr);
  double variance(base::Grid& grid, base::DataVector& alpha, base::DataMatrix* bounds = nullptr);

 protected:
  base::Grid& grid;
};

}  // namespace datadriven
}  // namespace sgpp
