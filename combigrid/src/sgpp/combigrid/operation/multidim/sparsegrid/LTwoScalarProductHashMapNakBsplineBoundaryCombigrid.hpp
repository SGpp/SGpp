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
#include <map>
#include <vector>

namespace sgpp {
namespace combigrid {

class LTwoScalarProductHashMapNakBsplineBoundaryCombigrid {
 public:
  LTwoScalarProductHashMapNakBsplineBoundaryCombigrid() {
    degree = 3;
    grid = nullptr;
    numAdditionalPoints = 0;
    isCustomWeightFunction = false;
  };
  /**
   * Constructor
   *
   * @param grid sparse grid created by converting a expUniformBoundaryGrid to a sgpp::base::Grid
   */
  LTwoScalarProductHashMapNakBsplineBoundaryCombigrid(sgpp::base::Grid* grid);

  LTwoScalarProductHashMapNakBsplineBoundaryCombigrid(
      sgpp::base::Grid* grid, sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection);

  LTwoScalarProductHashMapNakBsplineBoundaryCombigrid(
      sgpp::base::Grid* grid, sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection,
      sgpp::base::DataVector bounds);

  LTwoScalarProductHashMapNakBsplineBoundaryCombigrid(
      sgpp::base::Grid* grid, sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection,
      sgpp::base::DataVector bounds, size_t numAdditionalPoints);

  /**
   * Destructor
   */
  virtual ~LTwoScalarProductHashMapNakBsplineBoundaryCombigrid();

  /**
  * Implementation of standard matrix multiplication
  *
  * @param alpha DataVector that is multiplied to the matrix
  * @param result DataVector into which the result of multiplication is stored
  */
  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  void setWeightFunction(sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection) {
    this->isCustomWeightFunction = true;
    this->weightFunctionsCollection = weightFunctionsCollection;
  };

  void setBounds(sgpp::base::DataVector bounds) { this->bounds = bounds; };

  void updateGrid(sgpp::base::Grid* grid);

  /**
   * Creates hash key from the two level-index pairs of two 1D Bsplines
   */
  void hashLevelIndex(base::level_t li, base::index_t ii, base::level_t lj, base::index_t ij,
                      size_t d, MultiIndex& hashMI);

  /**
   * subroutine to calculate the scalar product of the B splines i and j with level index pairs
   * (li,ii) and (lj,ij)
   */
  double calculateScalarProduct(base::level_t lid, base::index_t iid, base::level_t ljd,
                                base::index_t ijd, size_t quadOrder, size_t d, double offseti_left,
                                double offseti_right, double offsetj_left, double offsetj_right,
                                sgpp::base::index_t hInvik, sgpp::base::index_t hInvjk, double hik,
                                double hjk, size_t pp1h);

 protected:
  sgpp::base::Grid* grid;
  size_t degree;
  sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection;
  bool isCustomWeightFunction;
  sgpp::base::DataVector bounds;
  // if a custom weightFunctionsCollection is used the standard quadrature order p+1 might not
  // suffice. It can be increased with numAdditionalPoints
  size_t numAdditionalPoints;
  // HashMap containing the scalar products of the 1D B spline basis functions. Access via combined
  // MultiIndex of the level-index pair of the two splines
  std::map<MultiIndex, double> innerProducts;
};
}  // namespace combigrid
}  // namespace sgpp
