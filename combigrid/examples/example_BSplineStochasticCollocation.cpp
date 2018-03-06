// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/ProbabilityDensityFunction1D.hpp>
#include <sgpp/combigrid/operation/multidim/RegularLevelManager.hpp>
#include <sgpp/combigrid/pce/BsplineStochasticCollocation.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>

/**
 * \page example_example_BSplineStochasticCollocation_cpp Stochastic Collocation with
 * B-Spline Combigrids
 *
 * In this example the mean and variance of an objective function are calculated with B-Spline
 * Stochastic Collocation. The objective function is f(x) = x transformed to [-1,3] and the moments
 * are calculated assuming a normal probability density function is used as a weight function w(x)
 * in the
 * integration. This means
 *
 * 			\f$E(f) = \int  f(x) w(x) dx \\
 * 			V(f) = E(f^2) - E(f)^2 \f$
 *
 */

/**
 * First we define the objective function,
 */

double objectiveFunction(sgpp::base::DataVector const& v) { return 4.0 * v[0] - 1.0; }

/**
 * set up the B-Splines
 */

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

  /**
   * and intialize the probability density functions for our input.
   */

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
    bounds.push_back(pdf_config.pdfParameters.lowerBound_);
    bounds.push_back(pdf_config.pdfParameters.upperBound_);
    weightFunctionsCollection.push_back(oneDimensionsalWeightFunction);
  }

  /**
   * After that we create a B-Spline stochastic collocation surrogate
   */

  // set up the configuration for the B spline Stochastic Collocation
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> storage;
  sgpp::combigrid::CombigridSurrogateModelConfiguration config;
  config.type = sgpp::combigrid::CombigridSurrogateModelsType::BSPLINE_STOCHASTIC_COLLOCATION;
  config.pointHierarchies = pointHierarchies;
  config.levelManager = std::make_shared<sgpp::combigrid::RegularLevelManager>();
  config.degree = degree;
  config.coefficientStorage = storage;
  config.weightFunctions = weightFunctionsCollection;
  config.bounds = bounds;
  // create the B spline Stochastic Collocation
  sgpp::combigrid::BsplineStochasticCollocation bsc(config);

  /**
   * and finally finish constructing the combigrid operation.
   */

  // create an operation with a regular level manager to create the level structure
  sgpp::combigrid::EvaluatorConfiguration EvalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineInterpolation, config.degree);
  sgpp::combigrid::CombiEvaluators::MultiCollection Evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(EvalConfig));
  std::shared_ptr<sgpp::combigrid::LevelManager> regularLevelManager(
      new sgpp::combigrid::RegularLevelManager());
  sgpp::combigrid::FullGridSummationStrategyType auxiliarySummationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  bool exploitNesting = false;
  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, Evaluators, regularLevelManager, gf, exploitNesting,
      auxiliarySummationStrategyType);

  /**
   * Now we can add some levels to the combigrid,
   */

  // create a regular level structure of level 1. Because the regularLevelManager is part of the
  // above interpolation operation the B spline interpolation coefficients are calculated during the
  // level structure creation. These coefficients can then be used for the quadratures that must be
  // done for the mean and variance calculations
  Operation->getLevelManager()->addRegularLevels(1);

  /**
   * update the surrogate accordingly
   */

  //  update the B Spline Stochastic Collocation configuration with the level stucture and the
  //  interpolation coefficients from the refinement operation
  config.levelStructure = Operation->getLevelManager()->getLevelStructure();
  config.coefficientStorage = Operation->getStorage();
  // update the B Spline Stochastic Collocation with the configuration containing the new level
  // structure and coefficients
  bsc.updateConfig(config);

  /**
   * and calculate mean and variance of the objective function.
   */

  double variance = bsc.variance();
  double mean = bsc.mean();
  std::cout << "number of grid points: " << Operation->getLevelManager()->numGridPoints()
            << std::endl;
  std::cout << "mean: " << mean << " variance: " << variance << std::endl;

  return 0;
}
