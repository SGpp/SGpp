// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLTWODOTLAPLACEVECTORIZEDLINEAROCL_HPP
#define OPERATIONLTWODOTLAPLACEVECTORIZEDLINEAROCL_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/parallel/pde/basis/common/OCLPDEKernels.hpp>
#include <sgpp/parallel/pde/operation/OperationParabolicPDEMatrixCombined.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {

/**
 * Implementation for linear functions of LTwoDotLaplace Operation, linear grids without boundaries
 * using OpenCL
 *
 */
class OperationLTwoDotLaplaceVectorizedLinearOCL : public OperationParabolicPDEMatrixCombined {
 private:
  sgpp::base::GridStorage* storage;
  sgpp::base::DataMatrix* level_;
  sgpp::base::DataMatrix* level_int_;
  sgpp::base::DataMatrix* index_;
  double* lcl_q;
  double* lcl_q_inv;
  sgpp::base::DataVector* lambda;

  OCLPDEKernels OCLPDEKernelsHandle;
  size_t padding_size;
  size_t sizepad;
  double* subresult;

 public:
  /**
   * Construtor of OperationLTwoDotLaplaceVectorizedLinearOCL
   *
   * @param storage Pointer to the grid's gridstorage obejct
   * @param lambda Vector which contains pre-factors for every dimension of the operator
   */
  OperationLTwoDotLaplaceVectorizedLinearOCL(sgpp::base::GridStorage* storage,
                                             sgpp::base::DataVector& lambda);

  /**
   * Construtor of OperationLTwoDotLaplaceVectorizedLinearOCL
   *
   * @param storage Pointer to the grid's gridstorage obejct
   */
  explicit OperationLTwoDotLaplaceVectorizedLinearOCL(sgpp::base::GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~OperationLTwoDotLaplaceVectorizedLinearOCL();

  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);
};
}  // namespace parallel
}  // namespace sgpp

#endif /* OPERATIONLTWODOTLAPLACEVECTORIZEDLINEAROCL_HPP */
