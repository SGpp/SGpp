// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLAPLACEMODLINEAR_HPP
#define OPERATIONLAPLACEMODLINEAR_HPP

#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * Implementation of Laplace for mod linear functions
 *
 */
class OperationLaplaceModLinear : public UpDownOneOpDim {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's sgpp::base::GridStorage object
   */
  explicit OperationLaplaceModLinear(sgpp::base::GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~OperationLaplaceModLinear();

 protected:
  virtual void up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void downOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void upOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);
};
}  // namespace pde
}  // namespace sgpp

#endif /* OPERATIONLAPLACEMODLINEAR_HPP */
