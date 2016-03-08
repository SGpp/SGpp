// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLAPLACELINEARSTRETCHEDBOUNDARY_HPP
#define OPERATIONLAPLACELINEARSTRETCHEDBOUNDARY_HPP

#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * Implementation of Laplace for linear functions with boundaries
 *
 */
class OperationLaplaceLinearStretchedBoundary : public UpDownOneOpDim {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's sgpp::base::GridStorage object
   */
  explicit OperationLaplaceLinearStretchedBoundary(sgpp::base::GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~OperationLaplaceLinearStretchedBoundary();

 protected:
  virtual void up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void downOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void upOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);
};
}  // namespace pde
}  // namespace sgpp

#endif /* OPERATIONLAPLACELINEARBOUNDARYSTRETCHED_HPP */
