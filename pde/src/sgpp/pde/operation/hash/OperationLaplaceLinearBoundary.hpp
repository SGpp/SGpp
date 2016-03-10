// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLAPLACELINEARBOUNDARY_HPP
#define OPERATIONLAPLACELINEARBOUNDARY_HPP

#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * Implementation of Laplace for linear functions with boundaries
 *
 */
class OperationLaplaceLinearBoundary : public UpDownOneOpDim {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's sgpp::base::GridStorage object
   */
  explicit OperationLaplaceLinearBoundary(sgpp::base::GridStorage* storage);

  /**
   * Constructor of OperationLaplaceLinear
   *
   * @param storage Pointer to the grid's gridstorage obejct
   * @param coef reference to a sgpp::base::DataVector object that contains the bilinear form's
   * constant coefficients; one per dimension
   */
  OperationLaplaceLinearBoundary(sgpp::base::GridStorage* storage, sgpp::base::DataVector& coef);

  /**
   * Destructor
   */
  virtual ~OperationLaplaceLinearBoundary();

 protected:
  virtual void up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void downOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void upOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);
};
}  // namespace pde
}  // namespace sgpp

#endif /* OPERATIONLAPLACELINEARBOUNDARY_HPP */
