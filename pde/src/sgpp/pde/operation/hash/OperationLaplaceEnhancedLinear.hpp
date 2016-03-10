// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLAPLACEENHANCEDLINEAR_HPP
#define OPERATIONLAPLACEENHANCEDLINEAR_HPP

#include <sgpp/pde/algorithm/UpDownOneOpDimEnhanced.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * Implements the Laplace operator based on
 * the UpDownOneOpDimEnhanced method.
 *
 */
class OperationLaplaceEnhancedLinear : public UpDownOneOpDimEnhanced {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's sgpp::base::GridStorage object
   */
  explicit OperationLaplaceEnhancedLinear(sgpp::base::GridStorage* storage);

  /**
   * Constructor of OperationLaplaceLinear
   *
   * @param storage Pointer to the grid's gridstorage obejct
   * @param coef reference to a sgpp::base::DataVector object that contains the bilinear form's
   * constant coefficients; one per dimension
   */
  OperationLaplaceEnhancedLinear(sgpp::base::GridStorage* storage, sgpp::base::DataVector& coef);

  /**
   * Destructor
   */
  virtual ~OperationLaplaceEnhancedLinear();

 protected:
  /**
   * Up-step
   *
   * @param dim dimension in which to apply the up-part
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   */
  virtual void up(sgpp::base::DataMatrix& alpha, sgpp::base::DataMatrix& result, size_t dim);

  /**
   * Down-step
   *
   * @param dim dimension in which to apply the down-part
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   */
  virtual void down(sgpp::base::DataMatrix& alpha, sgpp::base::DataMatrix& result, size_t dim);
};
}  // namespace pde
}  // namespace sgpp

#endif /* OPERATIONLAPLACEENHANCEDLINEAR_HPP */
