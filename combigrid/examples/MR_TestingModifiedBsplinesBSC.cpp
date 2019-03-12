//// Copyright (C) 2008-today The SG++ project
//// This file is part of the SG++ project. For conditions of distribution and
//// use, please see the copyright notice provided with SG++ or at
//// sgpp.sparsegrids.org
//
//#include <sgpp/combigrid/functions/ProbabilityDensityFunction1D.hpp>
//#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
//#include <sgpp/combigrid/operation/CombigridOperation.hpp>
//#include <sgpp/combigrid/operation/multidim/RegularLevelManager.hpp>
//#include <sgpp/combigrid/pce/BsplineStochasticCollocation.hpp>
//#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
//#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
//#include <sgpp/combigrid/utils/BSplineRoutines.hpp>
//#include <sgpp/combigrid/utils/CombigridModifiedNakBSplineBasis.hpp>
//#include "../../base/src/sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp"
//
// double objectiveFunction(sgpp::base::DataVector const& v) { return v[0]; }
//
//// BsplineRoutines has fix BSplineInterpolation or ModifiedBSplineInterpolation
//// BsplienRoutines A is wrong!
//
// int main() {
//  size_t numDimensions = 1;
//  size_t degree = 3;
//  sgpp::combigrid::MultiFunction func(objectiveFunction);
//
//  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
//      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
//  //  sgpp::combigrid::CombiEvaluators::Collection evaluators(
//  //      numDimensions, sgpp::combigrid::CombiEvaluators::ModifiedBSplineInterpolation(degree));
//  //  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager =
//  //      std::make_shared<sgpp::combigrid::AveragingLevelManager>();
//  //  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(func, pointHierarchies,
//  //  degree); bool exploitNesting = false; auto interpolationOperation =
//  //  std::make_shared<sgpp::combigrid::CombigridOperation>(
//  //      pointHierarchies, evaluators, levelManager, gf, exploitNesting);
//  //
//  //  sgpp::base::DataVector v(numDimensions, 0.77);
//  //  //  double res2 = interpolationOperation->evaluate(1, v);
//  //  interpolationOperation->setParameters(v);
//  //  interpolationOperation->getLevelManager()->addRegularLevels(4);
//  //  double res2 = interpolationOperation->getResult();
//  //  std::cout << func(v) << " " << res2 << std::endl;
//
//  //  std::cout << "xxxxxxxxxxxxxxx" << std::endl;
//  //  std::vector<double> points(7);
//  //  points[0] = 0.125;
//  //  points[0] = 0.25;
//  //  points[0] = 0.375;
//  //  points[0] = 0.5;
//  //  points[0] = 0.625;
//  //  points[1] = 0.75;
//  //  points[2] = 0.875;
//  //  for (double i = 0; i < 101; i++) {
//  //    double x = i / 100;
//  //    std::cout << x << " " << sgpp::combigrid::expUniformModifiedNakBspline(x, 3, 0, points) <<
//  "
//  //    "
//  //              << sgpp::combigrid::expUniformModifiedNakBspline(x, 3, 1, points) << "  "
//  //              << sgpp::combigrid::expUniformModifiedNakBspline(x, 3, 2, points) << "  "
//  //              << sgpp::combigrid::expUniformModifiedNakBspline(x, 3, 3, points) << "  "
//  //              << sgpp::combigrid::expUniformModifiedNakBspline(x, 3, 4, points) << "  "
//  //              << sgpp::combigrid::expUniformModifiedNakBspline(x, 3, 5, points) << "  "
//  //              << sgpp::combigrid::expUniformModifiedNakBspline(x, 3, 6, points) << std::endl;
//  //  }
//
//  // set up the weight function collection as normally distributed probability density functions
//  //  sgpp::combigrid::ProbabilityDensityFunction1DConfiguration pdf_config;
//  //  pdf_config.pdfParameters.type_ = sgpp::combigrid::ProbabilityDensityFunctionType::NORMAL;
//  //  pdf_config.pdfParameters.mean_ = 1.0;    // => we should obtain mean = 1.0
//  //  pdf_config.pdfParameters.stddev_ = 0.1;  // => we should obtain variance = 0.01
//  //  pdf_config.pdfParameters.lowerBound_ = -1;
//  //  pdf_config.pdfParameters.upperBound_ = 3;
//  //  auto probabilityDensityFunction =
//  //      std::make_shared<sgpp::combigrid::ProbabilityDensityFunction1D>(pdf_config);
//  //  sgpp::combigrid::SingleFunction oneDimensionsalWeightFunction =
//  //      probabilityDensityFunction->getWeightFunction();
//  //
//  //  // prepare the weight functions and the left and right point of their definition interval
//  //  // here we use the same weight function in every dimension
//  //  sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection;
//  //  sgpp::base::DataVector bounds;
//  //  for (size_t d = 0; d < numDimensions; d++) {
//  //    bounds.push_back(pdf_config.pdfParameters.lowerBound_);
//  //    bounds.push_back(pdf_config.pdfParameters.upperBound_);
//  //    weightFunctionsCollection.push_back(oneDimensionsalWeightFunction);
//  //  }
//
//  sgpp::base::DataVector bounds;
//  for (size_t d = 0; d < numDimensions; d++) {
//    bounds.push_back(0);
//    bounds.push_back(1);
//  }
//
//  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> storage;
//
//  sgpp::combigrid::CombigridSurrogateModelConfiguration config;
//  config.type = sgpp::combigrid::CombigridSurrogateModelsType::BSPLINE_STOCHASTIC_COLLOCATION;
//  config.pointHierarchies = pointHierarchies;
//  config.levelManager = std::make_shared<sgpp::combigrid::RegularLevelManager>();
//  config.degree = degree;
//  config.coefficientStorage = storage;
//  //  config.weightFunctions = weightFunctionsCollection;
//  config.bounds = bounds;
//
//  sgpp::combigrid::BsplineStochasticCollocation bsc(config);
//
//  //  sgpp::combigrid::EvaluatorConfiguration EvalConfig(
//  //      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineInterpolation, config.degree);
//  //  sgpp::combigrid::CombiEvaluators::MultiCollection Evaluators(
//  //      numDimensions,
//  //    sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(EvalConfig));
//  //  std::shared_ptr<sgpp::combigrid::LevelManager> regularLevelManager(
//  //      new sgpp::combigrid::RegularLevelManager());
//  //  sgpp::combigrid::FullGridSummationStrategyType auxiliarySummationStrategyType =
//  //      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
//  //  bool exploitNesting = false;
//  //  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
//  //      pointHierarchies, Evaluators, regularLevelManager, gf, exploitNesting,
//  //      auxiliarySummationStrategyType);
//
//  std::shared_ptr<sgpp::combigrid::CombigridOperation> Operation =
//      sgpp::combigrid::CombigridOperation::createExpUniformBoundaryBsplineInterpolation(
//          numDimensions, sgpp::combigrid::MultiFunction(objectiveFunction), degree);
//
//  Operation->getLevelManager()->addRegularLevels(1);
//
//  config.levelStructure = Operation->getLevelManager()->getLevelStructure();
//  config.coefficientStorage = Operation->getStorage();
//  bsc.updateConfig(config);
//  sgpp::base::DataVector v(numDimensions, 0.337);
//  std::cout << objectiveFunction(v) << " " << bsc.eval(v) << std::endl;
//
//  //  double variance = bsc.variance();
//  //  double mean = bsc.mean();
//  //  std::cout << "# grid points: " << Operation->getLevelManager()->numGridPoints() <<
//  //    std::endl;
//  //  std::cout << "mean: " << mean << " variance: " << variance << std::endl;
//
//  return 0;
//}
