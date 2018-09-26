#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
#include <sgpp/combigrid/operation/hierarchical/OperationWeightedQuadratureNakBsplineBoundary.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/RegularLevelManager.hpp>
#include <sgpp/combigrid/pce/BsplineStochasticCollocation.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <math.h>
#include <iostream>

#include "../src/sgpp/combigrid/operation/multidim/sparsegrid/LTwoScalarProductNakBsplineBoundary.hpp"
#include "../src/sgpp/combigrid/pce/HierarchicalStochasticCollocation.hpp"

double l2Error(std::shared_ptr<sgpp::base::Grid> surrogateGrid, sgpp::base::DataVector alpha,
               double (*objFunc)(sgpp::base::DataVector), size_t dim) {
  sgpp::optimization::InterpolantScalarFunction sparseGridSurrogate(*surrogateGrid, alpha);
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  size_t level = 10;
  grid->getGenerator().full(level);
  sgpp::base::DataVector diffSquare(gridStorage.getSize(), 0.0);
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
      p[j] = gp.getStandardCoordinate(j);
    }
    diffSquare[i] = std::pow(fabs(objFunc(p) - sparseGridSurrogate.eval(p)), 2);
  }
  std::cout << "l2 error = " << sqrt(diffSquare.sum()) << " calculated on " << gridStorage.getSize()
            << " points." << std::endl;
  return sqrt(diffSquare.sum());
}

double combiRegularBSC(size_t dim, size_t level, size_t degree,
                       double (*objFunc)(sgpp::base::DataVector),
                       sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection,
                       sgpp::base::DataVector bounds, size_t& numGridPoints) {
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      dim, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(
      sgpp::combigrid::MultiFunction(objFunc), pointHierarchies, degree);
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> storage;
  sgpp::combigrid::CombigridSurrogateModelConfiguration config;
  config.type = sgpp::combigrid::CombigridSurrogateModelsType::BSPLINE_STOCHASTIC_COLLOCATION;
  config.pointHierarchies = pointHierarchies;
  config.levelManager = std::make_shared<sgpp::combigrid::RegularLevelManager>();  // dummy
  config.degree = degree;
  config.coefficientStorage = storage;
  config.weightFunctions = weightFunctionsCollection;
  config.bounds = bounds;
  sgpp::combigrid::BsplineStochasticCollocation bsc(config);

  sgpp::combigrid::EvaluatorConfiguration EvalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineInterpolation, config.degree);
  sgpp::combigrid::CombiEvaluators::MultiCollection Evaluators(
      dim, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(EvalConfig));
  std::shared_ptr<sgpp::combigrid::LevelManager> regularLevelManager(
      new sgpp::combigrid::RegularLevelManager());
  sgpp::combigrid::FullGridSummationStrategyType auxiliarySummationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  bool exploitNesting = false;
  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, Evaluators, regularLevelManager, gf, exploitNesting,
      auxiliarySummationStrategyType);

  Operation->getLevelManager()->addRegularLevels(level);

  config.levelStructure = Operation->getLevelManager()->getLevelStructure();
  config.coefficientStorage = Operation->getStorage();
  bsc.updateConfig(config);

  double variance = bsc.variance();
  numGridPoints = bsc.numGridPoints();
  return variance;
}

double combiAdaptiveBSC(size_t dim, size_t maxnumgp, size_t degree,
                        double (*objFunc)(sgpp::base::DataVector),
                        sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection,
                        sgpp::base::DataVector bounds, size_t& numGridPoints) {
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      dim, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(
      sgpp::combigrid::MultiFunction(objFunc), pointHierarchies, degree);
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> storage;
  sgpp::combigrid::CombigridSurrogateModelConfiguration config;
  config.type = sgpp::combigrid::CombigridSurrogateModelsType::BSPLINE_STOCHASTIC_COLLOCATION;
  config.pointHierarchies = pointHierarchies;
  config.levelManager = std::make_shared<sgpp::combigrid::RegularLevelManager>();  // dummy
  config.degree = degree;
  config.coefficientStorage = storage;
  config.weightFunctions = weightFunctionsCollection;
  config.bounds = bounds;
  sgpp::combigrid::BsplineStochasticCollocation bsc(config);

  sgpp::combigrid::EvaluatorConfiguration EvalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineInterpolation, config.degree);
  sgpp::combigrid::CombiEvaluators::MultiCollection Evaluators(
      dim, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(EvalConfig));
  std::shared_ptr<sgpp::combigrid::LevelManager> adaptiveLevelManager(
      new sgpp::combigrid::AveragingLevelManager());
  sgpp::combigrid::FullGridSummationStrategyType auxiliarySummationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  bool exploitNesting = false;
  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, Evaluators, adaptiveLevelManager, gf, exploitNesting,
      auxiliarySummationStrategyType);

  // one initial level
  Operation->getLevelManager()->addRegularLevels(1);
  Operation->getLevelManager()->addLevelsAdaptive(maxnumgp);

  config.levelStructure = Operation->getLevelManager()->getLevelStructure();
  //  std::cout << maxnumgp << ": " << std::endl;
  //  Operation->getLevelManager()->printLevelStructure(
  //      Operation->getLevelManager()->getLevelStructure());
  config.coefficientStorage = Operation->getStorage();
  bsc.updateConfig(config);

  double variance = bsc.variance();
  numGridPoints = bsc.numGridPoints();
  return variance;
}

double func(sgpp::base::DataVector v) { return 4.0 * v[0] - 1.0 + v[1]; }
double genzFunction(sgpp::base::DataVector v) {
  double sum = 0;
  size_t d = v.getSize();
  for (size_t k = 0; k < d; k++) {
    sum += 4.5 * (static_cast<double>(k) + 0.5) / static_cast<double>(d) * v[k];
  }
  return cos(sum);
}
double weight(double x) { return sin(x); }

int main() {
  size_t dim = 2;
  size_t degree = 3;
  double (*objFunc)(sgpp::base::DataVector) = &genzFunction;
  //  sgpp::combigrid::SingleFunction oneDimensionsalWeightFunction(weight);
  sgpp::combigrid::ProbabilityDensityFunction1DConfiguration pdf_config;
  pdf_config.pdfParameters.type_ = sgpp::combigrid::ProbabilityDensityFunctionType::UNIFORM;
  pdf_config.pdfParameters.mean_ = 1.0;
  pdf_config.pdfParameters.stddev_ = 0.1;
  pdf_config.pdfParameters.lowerBound_ = -1;
  pdf_config.pdfParameters.upperBound_ = 3;
  auto probabilityDensityFunction =
      std::make_shared<sgpp::combigrid::ProbabilityDensityFunction1D>(pdf_config);
  sgpp::combigrid::SingleFunction oneDimensionsalWeightFunction =
      probabilityDensityFunction->getWeightFunction();
  sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection(dim);
  sgpp::base::DataVector bounds(2 * dim);
  for (size_t d = 0; d < dim; d++) {
    bounds[2 * d] = pdf_config.pdfParameters.lowerBound_;
    bounds[2 * d + 1] = pdf_config.pdfParameters.upperBound_;
    weightFunctionsCollection[d] = oneDimensionsalWeightFunction;
  }
  // statistics

  std::vector<double> hierError, hierTime, combiError, combiTime;
  std::vector<size_t> numHierGP, numCombiGP;
  size_t numGridPoints = 0;

  // genz 2D, uniform on [-1,3] with mean  =1, stddev  =0.1 (calculated on level 11)
  double realVar = 0.382923024983218;

  //  std::shared_ptr<sgpp::base::Grid> grid(
  //      sgpp::base::Grid::createNakBsplineModifiedGrid(dim, degree));
  //  sgpp::combigrid::HierarchicalStochasticCollocation hBSC(
  //      grid, sgpp::combigrid::MultiFunction(genzFunction), weightFunctionsCollection, bounds);
  //    sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineBoundary;
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineModified;
  sgpp::combigrid::HierarchicalStochasticCollocation hBSC(
      gridType, dim, sgpp::combigrid::MultiFunction(genzFunction), weightFunctionsCollection,
      bounds, degree);
  hBSC.refineRegular(1);
  for (size_t numRefine = 1; numRefine < 20; numRefine++) {
    sgpp::combigrid::Stopwatch watch;
    watch.start();
    hBSC.refineSurplusAdaptive(1);
    double hierVariance = hBSC.variance();
    hierTime.push_back(watch.elapsedSeconds());
    hierError.push_back(fabs(hierVariance - realVar));
    numHierGP.push_back(hBSC.numGridPoints());

    watch.start();
    //    double combiVariance = combiRegularBSC(dim, level, degree, objFunc,
    //    weightFunctionsCollection,
    //                                           bounds, numGridPoints);
    double combiVariance = combiAdaptiveBSC(dim, numHierGP.back(), degree, objFunc,
                                            weightFunctionsCollection, bounds, numGridPoints);
    combiTime.push_back(watch.elapsedSeconds());
    combiError.push_back(fabs(combiVariance - realVar));
    numCombiGP.push_back(numGridPoints);
  }
  for (size_t i = 0; i < hierError.size(); i++) {
    std::cout << i << ": " << hierError[i] << " " << combiError[i] << " | " << hierTime[i] << " "
              << combiTime[i] << " | " << numHierGP[i] << " " << numCombiGP[i] << std::endl;
  }
  return 0;
}
