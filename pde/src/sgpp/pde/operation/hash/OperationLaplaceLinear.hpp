// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLAPLACELINEAR_HPP
#define OPERATIONLAPLACELINEAR_HPP

#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * Implementation for linear functions of Laplace Operation, linear grids without boundaries
 *
 */
class OperationLaplaceLinear : public UpDownOneOpDim {
 public:
  /**
   * Constructor of OperationLaplaceLinear
   *
   * @param storage Pointer to the grid's gridstorage obejct
   */
  explicit OperationLaplaceLinear(sgpp::base::GridStorage* storage);

  /**
   * Constructor of OperationLaplaceLinear
   *
   * @param storage Pointer to the grid's gridstorage obejct
   * @param coef reference to a sgpp::base::DataVector object that contains the bilinear form's
   * constant coefficients; one per dimension
   */
  OperationLaplaceLinear(sgpp::base::GridStorage* storage, sgpp::base::DataVector& coef);

  /**
   * Destructor
   */
  virtual ~OperationLaplaceLinear();

  virtual void specialOP(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim,
                         size_t gradient_dim);

  virtual void up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void downOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void upOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);
};
}  // namespace pde
}  // namespace sgpp

#endif /* OPERATIONLAPLACELINEAR_HPP */
