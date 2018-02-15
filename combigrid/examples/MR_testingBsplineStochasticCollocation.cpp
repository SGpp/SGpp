// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/RegularLevelManager.hpp>
#include <sgpp/combigrid/pce/BsplineStochasticCollocation.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>

#include <sgpp/combigrid/utils/Stopwatch.hpp>

void createVarianceLevelStructure(
    size_t numPoints, size_t degree,
    sgpp::combigrid::CombiHierarchies::Collection const& pointHierarchies,
    sgpp::combigrid::GridFunction gf, bool exploitNesting, size_t numthreads,
    std::shared_ptr<sgpp::combigrid::LevelManager>& levelManager,
    std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage>& coefficientStorage,
    size_t numDimensions) {
  sgpp::combigrid::EvaluatorConfiguration EvalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineScalarProduct, degree);
  sgpp::combigrid::CombiEvaluators::MultiCollection Evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(EvalConfig));
  std::shared_ptr<sgpp::combigrid::LevelManager> varianceLevelManager(
      new sgpp::combigrid::AveragingLevelManager());
  sgpp::combigrid::FullGridSummationStrategyType auxiliarySummationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::VARIANCE;

  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, Evaluators, varianceLevelManager, gf, exploitNesting,
      auxiliarySummationStrategyType);

  Operation->getLevelManager()->addRegularLevels(1);
  Operation->getLevelManager()->addLevelsAdaptive(numPoints);

  levelManager = Operation->getLevelManager();
  coefficientStorage = Operation->getStorage();
}

void createRegularLevelStructure(
    size_t const& numlevels, size_t const& degree,
    sgpp::combigrid::CombiHierarchies::Collection const& pointHierarchies,
    sgpp::combigrid::GridFunction gf, bool exploitNesting, size_t numthreads,
    std::shared_ptr<sgpp::combigrid::LevelManager>& levelManager,
    std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage>& coefficientStorage,
    size_t numDimensions) {
  sgpp::combigrid::EvaluatorConfiguration EvalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineInterpolation, degree);
  sgpp::combigrid::CombiEvaluators::MultiCollection Evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(EvalConfig));
  sgpp::combigrid::FullGridSummationStrategyType auxiliarySummationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;

  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, Evaluators, levelManager, gf, exploitNesting,
      auxiliarySummationStrategyType);

  Operation->getLevelManager()->addRegularLevels(numlevels);
  levelManager = Operation->getLevelManager();
  coefficientStorage = Operation->getStorage();
}

double f(sgpp::base::DataVector const& v) { return (std::pow(v[0], 5) + std::pow(v[1], 5)); }
double wcos(double x) { return sin(x); }

int main() {
  size_t numDims = 2;
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDims, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> storage;
  //  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
  //      new sgpp::combigrid::AveragingLevelManager());
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
      new sgpp::combigrid::RegularLevelManager());
  sgpp::combigrid::SingleFunction weightfunction(wcos);
  sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection(0);
  sgpp::base::DataVector bounds(0);
  for (size_t d = 0; d < numDims; d++) {
    bounds.push_back(0);
    bounds.push_back(2);
    weightFunctionsCollection.push_back(weightfunction);
  }

  sgpp::combigrid::CombigridSurrogateModelConfiguration config;
  config.type = sgpp::combigrid::CombigridSurrogateModelsType::BSPLINE_STOCHASTIC_COLLOCATION;
  config.pointHierarchies = pointHierarchies;
  config.storage = storage;
  config.levelManager = levelManager;
  config.degree = 3;
  config.coefficientStorage = storage;
  config.weightFunctions = weightFunctionsCollection;

  config.bounds = bounds;
  sgpp::combigrid::BsplineStochasticCollocation BSC(config);

  std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> levelStructure;
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> coefficientStorage;
  sgpp::combigrid::MultiFunction func(f);
  sgpp::combigrid::GridFunction gf =
      BSplineCoefficientGridFunction(func, pointHierarchies, config.degree);

  //  std::vector<size_t> loopPoints{10, 100, 1000, 10000};
  std::vector<size_t> loopPoints{6, 7, 8};
  for (auto& np : loopPoints) {
    sgpp::combigrid::Stopwatch watch;
    watch.start();
    //    createVarianceLevelStructure(np, config.degree, pointHierarchies, gf, false, 4,
    //    levelManager,
    //                                 coefficientStorage, numDims);
    createRegularLevelStructure(np, config.degree, pointHierarchies, gf, false, 4, levelManager,
                                coefficientStorage, numDims);
    //  update config
    config.levelStructure = levelManager->getLevelStructure();
    config.coefficientStorage = coefficientStorage;
    BSC.updateConfig(config);

    // x^5+y^5  with weight function sin(x) on [0,2]^2
    double realEv = 0.906028496608237;
    double realVar = 0.700571273115382;

    // x^5+y^5 with weight function sin(x) on  [0,1]^2
    //    double realEv = 0.114999004731632;
    //    double realVar = 0.073775657714228;

    sgpp::combigrid::Stopwatch variance_watch;
    variance_watch.start();
    double var = BSC.variance();
    std::cout << levelManager->numGridPoints() << " " << variance_watch.elapsedSeconds()
              << "s  of ";
    double ev = BSC.mean();
    std::cout << watch.elapsedSeconds() << "s | ";
    //    std::cout << "mean: " << ev << " variance: " << var << std::endl;
    std::cout << fabs(ev - realEv) << " " << fabs(fabs(var) - fabs(realVar)) << std::endl;
  }

  return 0;
}
