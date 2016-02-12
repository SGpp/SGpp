// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLAPLACELINEAR_HPP
#define OPERATIONLAPLACELINEAR_HPP

#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
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
  explicit OperationLaplaceLinear(SGPP::base::GridStorage* storage);

  /**
   * Constructor of OperationLaplaceLinear
   *
   * @param storage Pointer to the grid's gridstorage obejct
   * @param coef reference to a SGPP::base::DataVector object that contains the bilinear form's
   * constant coefficients; one per dimension
   */
  OperationLaplaceLinear(SGPP::base::GridStorage* storage, SGPP::base::DataVector& coef);

  /**
   * Destructor
   */
  virtual ~OperationLaplaceLinear();

  virtual void specialOP(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim,
                         size_t gradient_dim);

  virtual void up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

  virtual void down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

  virtual void downOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

  virtual void upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);
};
}  // namespace pde
}  // namespace SGPP

#endif /* OPERATIONLAPLACELINEAR_HPP */
