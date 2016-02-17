// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONINVERSEROSENBLATTTRANSFORMATION_HPP
#define OPERATIONINVERSEROSENBLATTTRANSFORMATION_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

/**
 * Sampling on all dimensions
 */

class OperationInverseRosenblattTransformation {
 public:
  OperationInverseRosenblattTransformation() {}
  virtual ~OperationInverseRosenblattTransformation() {}

  /**
   * Rosenblatt Transformation with mixed starting dimension
   *
   * @param alpha Coefficient vector for current grid
   * @param pointscdf Input DataMatrix (rows: # of samples, columns: # of dims)
   * @param points Output DataMatrix (rows: # of samples, columns: # of dims)
   */
  virtual void doTransformation(base::DataVector* alpha, base::DataMatrix* pointscdf,
                                base::DataMatrix* points) = 0;

  /**
   * Rosenblatt Transformation with fixed starting dimension
   *
   * @param alpha Coefficient vector for current grid
   * @param pointscdf Input DataMatrix (rows: # of samples, columns: # of dims)
   * @param points Output DataMatrix (rows: # of samples, columns: # of dims)
   * @param dim_start starting dimension
   */
  virtual void doTransformation(base::DataVector* alpha, base::DataMatrix* pointscdf,
                                base::DataMatrix* points, size_t dim_start) = 0;
};
}  // namespace datadriven
}  // namespace SGPP
#endif /* OPERATIONINVERSEROSENBLATTTRANSFORMATION_HPP */
