// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * Implementation for BsplineBoundary functions of Laplace Operation, linear grids without boundaries
 */
class OperationLaplaceBsplineBoundary : public sgpp::base::OperationMatrix {
 public:
  /**
   * Constructor that creates an own matrix
   * i.e. matrix is destroyed by the destructor of OperationLaplaceBsplineBoundary
   *
   * @param grid the sparse grid
   */
  explicit OperationLaplaceBsplineBoundary(sgpp::base::Grid* grid);

  /**
   * Destructor
   */
  virtual ~OperationLaplaceBsplineBoundary();

  /**
   * Implementation of standard matrix multiplication
   *
   * @param alpha DataVector that is multiplied to the matrix
   * @param result DataVector into which the result of multiplication is stored
   */
  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

 private:
  sgpp::base::Grid* grid;
};

}  // namespace pde
}  // namespace sgpp
