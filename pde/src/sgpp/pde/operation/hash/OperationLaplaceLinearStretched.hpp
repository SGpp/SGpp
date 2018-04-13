// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLAPLACELINEARSTRETCHED_HPP
#define OPERATIONLAPLACELINEARSTRETCHED_HPP

#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * Implementation for linear functions of Laplace Operation, linear grids without boundaries
 *
 */
class OperationLaplaceLinearStretched : public UpDownOneOpDim {
 public:
  /**
   * Constructor of OperationLaplaceLinearStretched
   *
   * @param storage Pointer to the grid's gridstorage obejct
   */
  explicit OperationLaplaceLinearStretched(sgpp::base::GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~OperationLaplaceLinearStretched();

 protected:
  virtual void specialOP(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim,
                         size_t gradient_dim);

  virtual void up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void downOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  virtual void upOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);
};
}  // namespace pde
}  // namespace sgpp

#endif /* OPERATIONLAPLACELINEARSTRETCHED_HPP */
