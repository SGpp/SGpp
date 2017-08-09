// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONINVERSEROSENBLATTTRANSFORMATIONPOLY_HPP
#define OPERATIONINVERSEROSENBLATTTRANSFORMATIONPOLY_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * keep applying marginalize to function until it's reduced to only 1 dimension
 */

class OperationInverseRosenblattTransformationPoly
    : public OperationInverseRosenblattTransformation {
 public:
  explicit OperationInverseRosenblattTransformationPoly(base::Grid* grid) : grid(grid) {}
  virtual ~OperationInverseRosenblattTransformationPoly() {}
  /**
   * Transformation with mixed starting dimensions
   *
   * @param alpha Coefficient vector for current grid
     * @param pointscdf Input Matrix
     * @param points Output Matrix
   */
  void doTransformation(base::DataVector* alpha, base::DataMatrix* pointscdf,
                        base::DataMatrix* points);

  /**
     * Transformation with specified starting dimension
     *
     * @param alpha Coefficient vector for current grid
     * @param pointscdf Input Matrix
     * @param points Output Matrix
     * @param dim_start Starting dimension
     */
  void doTransformation(base::DataVector* alpha, base::DataMatrix* pointscdf,
                        base::DataMatrix* points, size_t dim_start);

 protected:
  base::Grid* grid;
  void doTransformation_start_dimX(base::Grid* g_in, base::DataVector* a_in, size_t dim_start,
                                   base::DataVector* cdfs1d, base::DataVector* coords1d);
  void doTransformation_in_next_dim(base::Grid* g_in, base::DataVector* a_in, size_t op_dim,
                                    base::DataVector* cdfs1d, base::DataVector* coords1d,
                                    size_t& curr_dim);
  double doTransformation1D(base::Grid* grid1d, base::DataVector* alpha1d, double coord1d);
};
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONINVERSEROSENBLATTTRANSFORMATIONPOLY_HPP */
