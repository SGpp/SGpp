// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


/**
 *
 * The following example shows how to integrate a function using Leja quadrature via the
 * combination technique. We first interpolate the function and afterwards we integrate 
 * the interpolant. The entire functionality is implemented in the function quadrature()
 * 
 *
 * We use the two dimensional test function
 * \f[
 *   f\colon [0, 1]^2 \to \mathbb{R},\quad
 *   f(x_0, x_1) := 4 x_0 **2 x_1 (1 - x_1)
 * \f]
 *
 * For instructions on how to compile and run the example, please see \ref installation.
 *
 *
 * This example can be found in the file quadrature.cpp
 */

// include all combigrid headers
#include <sgpp_combigrid.hpp>

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

using sgpp::base::DataVector;
using sgpp::combigrid::MultiFunction;
using sgpp::combigrid::CombigridOperation;
using sgpp::combigrid::AbstractPointHierarchy;
using sgpp::combigrid::NestedPointHierarchy;
using sgpp::combigrid::LejaPointDistribution;
using sgpp::combigrid::IdentityPointOrdering;
using sgpp::combigrid::LinearGrowthStrategy;
using sgpp::combigrid::AbstractLinearEvaluator;
using sgpp::combigrid::FloatArrayVector;
using sgpp::combigrid::FloatScalarVector;
using sgpp::combigrid::ArrayEvaluator;
using sgpp::combigrid::QuadratureEvaluator;
using sgpp::combigrid::CombigridTreeStorage;
using sgpp::combigrid::FullGridTensorEvaluator;
using sgpp::combigrid::FunctionLookupTable;
using sgpp::combigrid::CombigridEvaluator;

/**
 * The function we want to integrate
 */
double f_2D(DataVector v) { return 4.0 * v[0] * v[0] * (v[1] - v[1] * v[1]); }

/**
*  Function that implements Leja quadrature. First, it interpolates the function
*  on Leja points using Lagrange polynomials and afterwards in integrates the interpolant
*/
void quadrature() {
  // dimension of the integration problem
  size_t numDimensions = 2;

  std::vector<DataVector> gridPoints;

  auto dummyFunc = [&](DataVector const &param) -> double {
    gridPoints.push_back(param);
    return 0.0;
  };

  // create grid points hierarchies for all dimensions
  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies(numDimensions);

  // in this case, the grid points grow linearly with the level
  size_t growthFactor = 2; 

  // create a point hierarchy of nested points (Leja points), which are already ordered in the
  // nesting order
  pointHierarchies[0] = std::make_shared<NestedPointHierarchy>(
      std::make_shared<LejaPointDistribution>(),
      std::make_shared<IdentityPointOrdering>(std::make_shared<LinearGrowthStrategy>(growthFactor),
                                              false));  // not yet sorted
  pointHierarchies[1] = pointHierarchies[0];

  // evaluators, i.e. quadrature operators
  std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluators(numDimensions);
  evaluators[0] = std::make_shared<ArrayEvaluator<QuadratureEvaluator>>(
      true);  // true means that the operator needs a parameter
  evaluators[1] = evaluators[0];

  auto storage = std::make_shared<CombigridTreeStorage>(pointHierarchies, MultiFunction(dummyFunc));

  auto fullGridEval = std::make_shared<FullGridTensorEvaluator<FloatArrayVector>>(
      storage, evaluators, pointHierarchies);

  auto combiGridEval =
      std::make_shared<CombigridEvaluator<FloatArrayVector>>(numDimensions, fullGridEval);

  size_t maxLevelSum = 2;

  std::vector<FloatArrayVector> parameters(2);
  parameters[0] = FloatArrayVector(
      std::vector<FloatScalarVector>{FloatScalarVector(0.5), FloatScalarVector(0.3)});
  parameters[1] = FloatArrayVector(
      std::vector<FloatScalarVector>{FloatScalarVector(0.7), FloatScalarVector(0.27)});

  fullGridEval->setParameters(parameters);

  // here, the callback function inserts the grid points to gridPoints
  combiGridEval->addRegularLevels(maxLevelSum);  

  // ----- EVAL FUNCTION
  FunctionLookupTable funcLookup((MultiFunction(f_2D)));

  std::cout << "Grid points:\n";

  for (auto &gridPoint : gridPoints) {
    std::cout << gridPoint[0] << ", " << gridPoint[1] << "\n";
    funcLookup(gridPoint);
  }

  auto realStorage = std::make_shared<CombigridTreeStorage>(
      pointHierarchies, MultiFunction([&](DataVector const &param) { return funcLookup(param); }));

  auto realFullGridEval = std::make_shared<FullGridTensorEvaluator<FloatArrayVector>>(
      realStorage, evaluators, pointHierarchies);

  auto realCombiGridEval =
      std::make_shared<CombigridEvaluator<FloatArrayVector>>(numDimensions, realFullGridEval);

  realFullGridEval->setParameters(parameters);
  realCombiGridEval->addRegularLevels(maxLevelSum);

  FloatArrayVector result = realCombiGridEval->getValue();

  std::cout << "Results:\n";

  for (size_t i = 0; i < result.size(); ++i) {
    std::cout << result.get(i).getValue() << "\n";
  }
}

int main() {
  quadrature();

  return 0;
}