// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLTWODOTPRODUCTVECTORIZEDLINEAROCL_HPP
#define OPERATIONLTWODOTPRODUCTVECTORIZEDLINEAROCL_HPP

#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/parallel/pde/basis/common/OCLPDEKernels.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {

/**
 * Implements the standard L 2 scalar product on linear grids (no boundaries)
 *
 */
class OperationLTwoDotProductVectorizedLinearOCL : public sgpp::base::OperationMatrix {
 private:
  sgpp::base::GridStorage* storage;
  sgpp::base::DataMatrix* level_;
  sgpp::base::DataMatrix* level_int_;
  sgpp::base::DataMatrix* index_;
  double* lcl_q;
  OCLPDEKernels OCLPDEKernelsHandle;

 public:
  /**
   * Constructor
   *
   * @param storage the grid's sgpp::base::GridStorage object
   */
  explicit OperationLTwoDotProductVectorizedLinearOCL(sgpp::base::GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~OperationLTwoDotProductVectorizedLinearOCL();

 protected:
  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);
};
}  // namespace parallel
}  // namespace sgpp

#endif /* OPERATIONLTWODOTPRODUCTVECTORIZEDLINEAROCL_HPP */
