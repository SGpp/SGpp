// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLTWODOTLAPLACEVECTORIZEDLINEARBOUNDARYOCL_HPP
#define OPERATIONLTWODOTLAPLACEVECTORIZEDLINEARBOUNDARYOCL_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/parallel/pde/basis/common/OCLPDEKernels.hpp>
#include <sgpp/parallel/pde/operation/OperationParabolicPDEMatrixCombined.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {

/**
 * Implementation for linear functions of LTwoDotLaplace Operation, linear grids with boundaries
 *
 */
class OperationLTwoDotLaplaceVectorizedLinearBoundaryOCL
    : public OperationParabolicPDEMatrixCombined {
 private:
  sgpp::base::GridStorage* storage;
  sgpp::base::DataMatrix* level_;
  sgpp::base::DataMatrix* level_int_;
  sgpp::base::DataMatrix* index_;
  double* lcl_q;
  double* lcl_q_inv;
  sgpp::base::DataVector* lambda;
  OCLPDEKernels OCLPDEKernelsHandle;

  void mult_dirichlet(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

 public:
  /**
   * Construtor of OperationLTwoDotLaplaceLinear
   *
   * @param storage Pointer to the grid's gridstorage obejct
   * @param lambda the lambda parameter which is needed in some cases (Black-Scholes) to modify the
   * dimensional local values
   */
  OperationLTwoDotLaplaceVectorizedLinearBoundaryOCL(sgpp::base::GridStorage* storage,
                                                     sgpp::base::DataVector& lambda);

  /**
   * Construtor of OperationLTwoDotLaplaceLinear
   *
   * @param storage Pointer to the grid's gridstorage obejct
    */
  explicit OperationLTwoDotLaplaceVectorizedLinearBoundaryOCL(sgpp::base::GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~OperationLTwoDotLaplaceVectorizedLinearBoundaryOCL();

  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);
};
}  // namespace parallel
}  // namespace sgpp

#endif /* OPERATIONLTWODOTLAPLACEVECTORIZEDOCLLINEARBOUNDARY_HPP */
