// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/ProbabilityDensityFunction1D.hpp>
#include <sgpp/combigrid/operation/multidim/RegularLevelManager.hpp>
#include <sgpp/combigrid/pce/BsplineStochasticCollocation.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>

/*
 * In this example the mean and variance of an objective function are calculated with the B spline
 * Stochastic Collocation. The objective function is f(x,y) = x^5+y^5 and the moments are calculated
 * assuming a probability density function is used as a weight function w(x) in the integration.
 * This means
 * 			E(f) = \int \int f(x,y) w(x) w(y) dx dy
 * 			V(f) = E(f^2) - E(f)^2
 *
 */

// w.l.o.g. the objective function must be defined on the unit cube
double objectiveFunction(sgpp::base::DataVector const& v) {
  return 4.0 * v[0] - 1.0;  // (std::pow(v[0], 5) + std::pow(v[1], 5));
}

int main() {
  // set the number of dimensions
  size_t numDimensions = 1;
  // set the degree of the B splien basis functions
  size_t degree = 5;
  // set the pointHierarchies (= creation pattern of the grid points) to exponentially growing
  // uniformly spaced grid points including the boudnary
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());

  // set up the grid function. This function takes full grids as its argument and returns the
  // corresponding coefficients for the B spline interpolation
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(
      sgpp::combigrid::MultiFunction(objectiveFunction), pointHierarchies, degree);

  // set up the weight function collection as normally distributed probability density functions
  sgpp::combigrid::ProbabilityDensityFunction1DConfiguration pdf_config;
  pdf_config.pdfParameters.type_ = sgpp::combigrid::ProbabilityDensityFunctionType::NORMAL;
  pdf_config.pdfParameters.mean_ = 1.0;    // => we should obtain mean = 1.0
  pdf_config.pdfParameters.stddev_ = 0.1;  // => we should obtain variance = 0.01
  pdf_config.pdfParameters.lowerBound_ = -1;
  pdf_config.pdfParameters.upperBound_ = 3;
  auto probabilityDensityFunction =
      std::make_shared<sgpp::combigrid::ProbabilityDensityFunction1D>(pdf_config);
  sgpp::combigrid::SingleFunction oneDimensionsalWeightFunction =
      probabilityDensityFunction->getWeightFunction();
  sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection;
  sgpp::base::DataVector bounds;
  for (size_t d = 0; d < numDimensions; d++) {
    bounds.push_back(-1);
    bounds.push_back(3);
    weightFunctionsCollection.push_back(oneDimensionsalWeightFunction);
  }

  // set up the configuration for the B spline Stochastic Collocation
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> storage;
  sgpp::combigrid::CombigridSurrogateModelConfiguration config;
  config.type = sgpp::combigrid::CombigridSurrogateModelsType::BSPLINE_STOCHASTIC_COLLOCATION;
  config.pointHierarchies = pointHierarchies;
  config.levelManager = std::make_shared<sgpp::combigrid::AveragingLevelManager>();
  config.degree = degree;
  config.coefficientStorage = storage;
  config.weightFunctions = weightFunctionsCollection;
  config.bounds = bounds;
  // create the B spline Stochastic Collocation
  sgpp::combigrid::BsplineStochasticCollocation bsc(config);

  //  createVarianceLevelStructure

  sgpp::combigrid::EvaluatorConfiguration EvalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineInterpolation, config.degree);
  sgpp::combigrid::CombiEvaluators::MultiCollection Evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(EvalConfig));
  std::shared_ptr<sgpp::combigrid::LevelManager> varianceLevelManager(
      new sgpp::combigrid::RegularLevelManager());
  sgpp::combigrid::FullGridSummationStrategyType auxiliarySummationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  bool exploitNesting = false;
  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, Evaluators, varianceLevelManager, gf, exploitNesting,
      auxiliarySummationStrategyType);

  // define an upper bound for the number of grid points
  //  Operation->getLevelManager()->addLevelsAdaptive(maximumNumPoints);
  Operation->getLevelManager()->addRegularLevels(1);

  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager = Operation->getLevelManager();
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> coefficientStorage =
      Operation->getStorage();

  //  update config
  config.levelStructure = levelManager->getLevelStructure();
  config.coefficientStorage = coefficientStorage;
  bsc.updateConfig(config);

  double variance = bsc.variance();
  double ev = bsc.mean();
  std::cout << "number of grid points: " << levelManager->numGridPoints() << std::endl;
  std::cout << "mean: " << ev << " variance: " << variance << std::endl;

  return 0;
}
