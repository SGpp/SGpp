// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/NakBsplineBoundaryCombigridGrid.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>

namespace sgpp {
namespace combigrid {

/**
 * Implements the standard L 2 scalar product for B splines on an uniform boundary grid created by
 * conversion from an expUniformBoundaryGrid from the combigrid module
 *
 */
class OperationMatrixLTwoDotNakBsplineBoundaryCombigrid : public sgpp::base::OperationMatrix {
 public:
  /**
   * Constructor
   *
   * @param grid sparse grid created by converting a expUniformBoundaryGrid to a sgpp::base::Grid
   */
  explicit OperationMatrixLTwoDotNakBsplineBoundaryCombigrid(sgpp::base::Grid* grid);

  explicit OperationMatrixLTwoDotNakBsplineBoundaryCombigrid(
      sgpp::base::Grid* grid, sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection);

  explicit OperationMatrixLTwoDotNakBsplineBoundaryCombigrid(
      sgpp::base::Grid* grid, sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection,
      sgpp::base::DataVector bounds);

  /**
   * Destructor
   */
  virtual ~OperationMatrixLTwoDotNakBsplineBoundaryCombigrid();

  /**
  * Implementation of standard matrix multiplication
  *
  * @param alpha DataVector that is multiplied to the matrix
  * @param result DataVector into which the result of multiplication is stored
  */
  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

 protected:
  sgpp::base::Grid* grid;
  sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection;
  sgpp::base::DataVector bounds;
};
}  // namespace combigrid
}  // namespace sgpp
