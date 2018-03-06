// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * Implementation for B-spline functions of Laplace Operation, linear grids without boundaries
 */
class OperationLaplaceExplicitLinear : public UpDownOneOpDim {
 public:
  /**
   * Constructor that uses a external matrix pointer to construct the matrix,
   * i.e. matrix is NOT destroyed by the destructor of OperationLaplaceExplicitLinear
   *
   * @param m pointer to datamatrix of size (number of grid point) x (number of grid points)
   * @param storage pointer to the sparse grid storage
   */
  OperationLaplaceExplicitLinear(sgpp::base::DataMatrix* m, sgpp::base::GridStorage* storage);
  /**
   * Constructor that creates an own matrix
   * i.e. matrix is destroyed by the destructor of OperationLaplaceExplicitLinear
   *
   * @param storage pointer to the sparse grid storage
   */
  explicit OperationLaplaceExplicitLinear(sgpp::base::GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~OperationLaplaceExplicitLinear();

  /**
   * Implementation of standard matrix multiplication
   *
   * @param alpha DataVector that is multiplied to the matrix
   * @param result DataVector into which the result of multiplication is stored
   */
  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);
  virtual void specialOP(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim,
                         size_t gradient_dim);

  virtual void up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void downOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void upOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

 private:
  /**
   * This method is used by both constructors to build the matrix
   */
  void buildMatrix(sgpp::base::GridStorage* storage);

  sgpp::base::DataMatrix* m_;
  bool ownsMatrix_;
};

}  // namespace pde
}  // namespace sgpp
