// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/combigrid/functions/ProbabilityDensityFunction1D.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/pce/BsplineStochasticCollocation.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>

#include <sgpp/combigrid/utils/Stopwatch.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace sgpp {
namespace combigrid {

BsplineStochasticCollocation::BsplineStochasticCollocation(
    sgpp::combigrid::CombigridSurrogateModelConfiguration& config)
    : CombigridSurrogateModel(config),
      weightFunctions(config.weightFunctions),
      currentNumGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0),
      coefficientStorage(config.coefficientStorage),
      scalarProducts() {
  initializeOperations(config.pointHierarchies, coefficientStorage, config.levelManager);
}

BsplineStochasticCollocation::~BsplineStochasticCollocation() {}

// combigridMultiOperation is a BSplineInterpolation
// combigridOperation is a BSplineQuadrature
void BsplineStochasticCollocation::initializeOperations(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::shared_ptr<AbstractCombigridStorage> coefficientStorage,
    std::shared_ptr<LevelManager> levelManager) {
  customWeightFunction = true;
  numDims = weightFunctions.size();
  // If no custom weight function is used the numDimensions MUST be set in the config!
  if (numDims == 0) {
    numDims = config.numDimensions;
    customWeightFunction = false;
  }
  // initialize interpolation operation
  sgpp::combigrid::EvaluatorConfiguration evalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineInterpolation, this->config.degree);
  sgpp::combigrid::CombiEvaluators::MultiCollection evaluators(
      numDims, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(evalConfig));
  sgpp::combigrid::FullGridSummationStrategyType summationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  auto interpolationOperation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, evaluators, levelManager, coefficientStorage, summationStrategyType);

  combigridMultiOperation = interpolationOperation;

  size_t numAdditionalPoints = 0;
  bool normalizeWeights = false;
  sgpp::combigrid::CombiEvaluators::Collection quadEvaluators(0);
  if (customWeightFunction) {
    for (size_t d = 0; d < numDims; d++) {
      quadEvaluators.push_back(sgpp::combigrid::CombiEvaluators::BSplineQuadrature(
          config.degree, weightFunctions[d], numAdditionalPoints, config.bounds[2 * d],
          config.bounds[2 * d + 1], normalizeWeights));
    }

    scalarProducts.setWeightFunction(weightFunctions);
    scalarProducts.setBounds(config.bounds);
  } else {
    for (size_t d = 0; d < numDims; d++) {
      quadEvaluators.push_back(sgpp::combigrid::CombiEvaluators::BSplineQuadrature(config.degree));
    }
  }
  auto quadratureOperation = std::make_shared<sgpp::combigrid::CombigridOperation>(
      pointHierarchies, quadEvaluators, levelManager, coefficientStorage, summationStrategyType);
  combigridOperation = quadratureOperation;

  currentNumGridPoints = 0;
}

// ToDo (rehmemk) tried to use addLevelsFromStructureParallel here, does not work

void BsplineStochasticCollocation::updateConfig(
    sgpp::combigrid::CombigridSurrogateModelConfiguration newConfig) {
  this->config.coefficientStorage = newConfig.coefficientStorage;
  this->config.levelStructure = newConfig.levelStructure;

  combigridMultiOperation = CombigridMultiOperation::createBsplineLinearCoefficientOperation(
      newConfig.degree, numDims, newConfig.coefficientStorage);
  combigridMultiOperation->getLevelManager()->addLevelsFromStructure(newConfig.levelStructure);

  size_t numAdditionalPoints = 0;
  bool normalizeWeights = false;
  sgpp::combigrid::FullGridSummationStrategyType summationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  sgpp::combigrid::CombiEvaluators::Collection quadEvaluators(0);

  if (customWeightFunction) {
    for (size_t d = 0; d < numDims; d++) {
      quadEvaluators.push_back(sgpp::combigrid::CombiEvaluators::BSplineQuadrature(
          config.degree, weightFunctions[d], numAdditionalPoints, config.bounds[2 * d],
          config.bounds[2 * d + 1], normalizeWeights));
    }
  } else {
    for (size_t d = 0; d < numDims; d++) {
      quadEvaluators.push_back(sgpp::combigrid::CombiEvaluators::BSplineQuadrature(config.degree));
    }
  }
  combigridOperation = std::make_shared<sgpp::combigrid::CombigridOperation>(
      config.pointHierarchies, quadEvaluators, config.levelManager, newConfig.coefficientStorage,
      summationStrategyType);
  combigridOperation->getLevelManager()->addLevelsFromStructure(newConfig.levelStructure);

  computedMeanFlag = false;
  computedVarianceFlag = false;
  hierarchicalTransformationFlag = false;
}

bool BsplineStochasticCollocation::updateStatus() {
  if (currentNumGridPoints < combigridMultiOperation->numGridPoints()) {
    coefficientStorage = combigridMultiOperation->getStorage();
    currentNumGridPoints = combigridMultiOperation->numGridPoints();
    computedMeanFlag = false;
    computedVarianceFlag = false;
    computedDiscreteMeanFlag = false;
    computedDiscreteVarianceFlag = false;
    return true;
  } else {
    return false;
  }
}

void BsplineStochasticCollocation::eval(sgpp::base::DataMatrix& xs, sgpp::base::DataVector& res) {
  combigridMultiOperation->setParameters(xs);
  combigridMultiOperation->getLevelManager()->addLevelsFromStructure(config.levelStructure);
  res = combigridMultiOperation->getResult();
}

double BsplineStochasticCollocation::eval(sgpp::base::DataVector& x) {
  sgpp::base::DataMatrix xs(x.size(), 0);
  xs.appendCol(x);
  sgpp::base::DataVector res;
  this->eval(xs, res);
  return res[0];
}

double BsplineStochasticCollocation::computeMean() {
  double mean = combigridOperation->getResult();
  return mean;
}

double BsplineStochasticCollocation::computeDiscreteMean(sgpp::base::DataMatrix discretePoints) {
  if (discretePoints.getNrows() != numDims) {
    throw sgpp::base::data_exception(
        "BsplineStochasticCollocation::computeDiscreteMean : Dimensions do not match");
  }
  sgpp::base::DataVector evaluations(discretePoints.getNcols());
  eval(discretePoints, evaluations);

  return evaluations.sum() / static_cast<double>(discretePoints.getNcols());
}

double BsplineStochasticCollocation::mean() {
  updateStatus();
  if (!computedMeanFlag) {
    ev = computeMean();
    computedMeanFlag = true;
  }
  return ev;
}

double BsplineStochasticCollocation::discreteMean(sgpp::base::DataMatrix discretePoints) {
  updateStatus();
  if (!computedDiscreteMeanFlag) {
    discreteEV = computeDiscreteMean(discretePoints);
    computedDiscreteMeanFlag = true;
  }
  return discreteEV;
}

double BsplineStochasticCollocation::computeVariance() {
  if (!computedMeanFlag) {
    mean();
  }
  // compute hierarchical sparse grid interpolant
  if (hierarchicalTransformationFlag == false) {
    transformToHierarchical();
  }

  sgpp::base::Grid* gridptr = hierarchicalGrid.get();
  sgpp::base::DataVector product(hierarchicalCoefficients.size());
  scalarProducts.updateGrid(gridptr);
  // scalarProducts.setWeightFunction(weightFunctions);
  // scalarProducts.setBounds(config.bounds);
  double variance = 0;
  if (config.degree == 1) {
    // calculate V(u) = E(u^2) - E(u)^2
    // this works for all B spline degrees but may be instable
    //(e.g. ev*ev might be larger than meanSquare => negative variance)
    scalarProducts.mult(hierarchicalCoefficients, product);
    double meanSquare = product.dotProduct(hierarchicalCoefficients);
    variance = meanSquare - ev * ev;
  } else {
    // calculate V(u) = E((u-E(u))^2)
    // this is done by subtracting E(u) from the coefficient of the constant function
    // it does not work for B spline degree 1 because there is no constant function in the basis
    // (We could add a constant basis function on level 0)
    sgpp::base::DataVector copyOfHierarchicalCoefficients(hierarchicalCoefficients);
    copyOfHierarchicalCoefficients[0] -= ev;
    scalarProducts.mult(copyOfHierarchicalCoefficients, product);

    variance = product.dotProduct(copyOfHierarchicalCoefficients);
  }

  return variance;
}

double BsplineStochasticCollocation::computeDiscreteVariance(
    sgpp::base::DataMatrix discretePoints) {
  if (discretePoints.getNrows() != numDims) {
    throw sgpp::base::data_exception(
        "BsplineStochasticCollocation::computeDiscreteVariance : Dimensions do not match");
  }
  if (!computedDiscreteMeanFlag) {
    discreteMean(discretePoints);
  }
  sgpp::base::DataVector evaluations(discretePoints.getNcols());
  eval(discretePoints, evaluations);

  discreteVar = 0.0;
  for (size_t i = 0; i < evaluations.getSize(); i++) {
    discreteVar += std::pow(evaluations[i] - ev, 2);
  }
  discreteVar /= static_cast<double>(discretePoints.getNcols() - 1);

  return discreteVar;
}

double BsplineStochasticCollocation::variance() {
  updateStatus();
  if (!computedVarianceFlag) {
    var = computeVariance();
    computedVarianceFlag = true;
  }
  return var;
}

double BsplineStochasticCollocation::discreteVariance(sgpp::base::DataMatrix discretePoints) {
  updateStatus();
  if (!computedDiscreteVarianceFlag) {
    discreteVar = computeDiscreteVariance(discretePoints);
    computedDiscreteVarianceFlag = true;
  }
  return discreteVar;
}

double BsplineStochasticCollocation::computeVarianceWithCombiTechnique() {
  if (!computedMeanFlag) {
    mean();
  }

  auto pointHierarchies = combigridOperation->getPointHierarchies();
  auto levelManager = combigridOperation->getLevelManager();
  sgpp::combigrid::FullGridSummationStrategyType summationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::QUADRATIC;

  // Scalar products support weight functions. First test some easy case without them!
  sgpp::combigrid::CombiEvaluators::MultiCollection scalarProductEvaluators(0);
  if (customWeightFunction) {
    //    std::cerr
    //        << "the combination technique scalar products currently do not support weight
    //        functions!"
    //        << std::endl;
    size_t numAdditionalPoints = 0;
    bool normalizeWeights = false;
    for (size_t d = 0; d < numDims; d++) {
      scalarProductEvaluators.push_back(sgpp::combigrid::CombiEvaluators::BSplineScalarProduct(
          config.degree, weightFunctions[d], numAdditionalPoints, config.bounds[2 * d],
          config.bounds[2 * d + 1], normalizeWeights));
    }
  } else {
    for (size_t d = 0; d < numDims; d++) {
      scalarProductEvaluators.push_back(
          sgpp::combigrid::CombiEvaluators::BSplineScalarProduct(config.degree));
    }
  }

  auto varianceOperation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, scalarProductEvaluators, levelManager, this->config.coefficientStorage,
      summationStrategyType);

  varianceOperation->getLevelManager()->addLevelsFromStructure(this->config.levelStructure);

  sgpp::base::DataVector meanSquareDataVector = varianceOperation->getResult();
  double meanSquare = meanSquareDataVector[0];

  return meanSquare - ev * ev;
}

// todo (rehmemk) when HierarchicalBsplineStochasticCollocation is finished remove the coarsening
// routines here and refer to HierarchicalBsplineStochasticCollocation

void BsplineStochasticCollocation::coarsen(size_t removements_num, double threshold,
                                           bool recalculateCoefficients) {
  sgpp::base::HashCoarsening coarsen;
  sgpp::base::GridStorage& gridStorage = hierarchicalGrid->getStorage();
  sgpp::base::SurplusCoarseningFunctor functor(hierarchicalCoefficients, removements_num,
                                               threshold);
  if (recalculateCoefficients) {
    std::shared_ptr<sgpp::base::Grid> oldGrid(hierarchicalGrid);
    sgpp::base::DataVector oldCoefficients(hierarchicalCoefficients);
    coarsen.free_coarsen(gridStorage, functor, hierarchicalCoefficients);
    hierarchicalCoefficients = recalculateInterpolationCoefficientsAfterCoarsening(
        oldGrid, hierarchicalGrid, oldCoefficients);
  } else {
    coarsen.free_coarsen(gridStorage, functor, hierarchicalCoefficients);
  }
  computedMeanFlag = false;
  computedVarianceFlag = false;
}

void BsplineStochasticCollocation::coarsenIteratively(double maxVarDiff, double threshold,
                                                      double removements_percentage,
                                                      bool recalculateCoefficients) {
  if ((removements_percentage < 0) || (removements_percentage > 100)) {
    removements_percentage = 10;
  }

  size_t removements_num = static_cast<size_t>(
      floor(static_cast<double>(hierarchicalGrid->getSize()) * removements_percentage / 100));
  double oldVariance = var;

  // coarsen
  sgpp::base::HashCoarsening coarsen;
  sgpp::base::GridStorage& gridStorage = hierarchicalGrid->getStorage();
  sgpp::base::SurplusCoarseningFunctor functor(hierarchicalCoefficients, removements_num,
                                               threshold);
  if (recalculateCoefficients) {
    std::shared_ptr<sgpp::base::Grid> oldGrid(hierarchicalGrid);
    sgpp::base::DataVector oldCoefficients(hierarchicalCoefficients);
    coarsen.free_coarsen(gridStorage, functor, hierarchicalCoefficients);
    hierarchicalCoefficients = recalculateInterpolationCoefficientsAfterCoarsening(
        oldGrid, hierarchicalGrid, oldCoefficients);
  } else {
    coarsen.free_coarsen(gridStorage, functor, hierarchicalCoefficients);
  }
  computedMeanFlag = false;
  computedVarianceFlag = false;

  // check for abort criterion and coarsen further if it's not met
  var = computeVariance();
  double varDiff = fabs(var - oldVariance);
  size_t numDeletedPoints = coarsen.getDeletedPoints().size();
  std::cout << "BSC coarsening, removed " << numDeletedPoints
            << " points, variance difference: " << varDiff << std::endl;

  if ((varDiff <= maxVarDiff) && (numDeletedPoints != 0)) {
    coarsenIteratively(maxVarDiff, threshold, removements_percentage, recalculateCoefficients);
  }
}

void BsplineStochasticCollocation::transformToHierarchical() {
  hierarchicalGrid.reset(
      sgpp::base::Grid::createNakBsplineBoundaryCombigridGrid(numDims, config.degree));
  sgpp::base::GridStorage& gridStorage = hierarchicalGrid->getStorage();
  auto levelStructure = this->config.levelStructure;
  convertexpUniformBoundaryCombigridToHierarchicalSparseGrid(levelStructure, gridStorage);

  hierarchicalCoefficients = calculateInterpolationCoefficientsForConvertedCombigird(
      hierarchicalGrid, gridStorage, combigridMultiOperation, levelStructure);
  hierarchicalTransformationFlag = true;
}

void BsplineStochasticCollocation::getComponentSobolIndices(
    sgpp::base::DataVector& componentSsobolIndices, bool normalized) {
  throw sgpp::base::application_exception(
      "BsplineStochasticCollocation::getComponentSobolIndices - not implemented.");
}
void BsplineStochasticCollocation::getTotalSobolIndices(sgpp::base::DataVector& totalSobolIndices,
                                                        bool normalized) {
  throw sgpp::base::application_exception(
      "PolynomialStochasticCollocation::getTotalSobolIndices - not implemented.");
}

size_t BsplineStochasticCollocation::numGridPoints() { return currentNumGridPoints; }

size_t BsplineStochasticCollocation::numHierarchicalGridPoints() {
  return hierarchicalGrid->getSize();
}

std::shared_ptr<LevelInfos> BsplineStochasticCollocation::getInfoOnAddedLevels() {
  return combigridOperation->getLevelManager()->getInfoOnAddedLevels();
}

void BsplineStochasticCollocation::differenceCTSG(sgpp::base::DataMatrix& xs,
                                                  sgpp::base::DataVector& res) {
  // compute hierarchical sparse grid interpolant
  if (hierarchicalTransformationFlag == false) {
    transformToHierarchical();
  }
  sgpp::optimization::InterpolantScalarFunction sparseGridSurrogate(*hierarchicalGrid,
                                                                    hierarchicalCoefficients);

  // evaluate hierarchical sparse grid surrogate in xs
  size_t numPoints = xs.getNcols();
  size_t dim = xs.getNrows();
  sgpp::base::DataVector point(dim, 0.0);
  sgpp::base::DataVector valuesSG(numPoints, 0.0);
  for (size_t i = 0; i < numPoints; i++) {
    xs.getColumn(i, point);
    valuesSG[i] = sparseGridSurrogate.eval(point);
  }
  // evaluate combi technique surrogate in xs
  sgpp::base::DataVector valuesCT(numPoints, 0.0);
  eval(xs, valuesCT);

  // difference |u_{CT}(p) - u_{SG}(p)|
  valuesCT.sub(valuesSG);
  res.resize(numPoints);
  for (size_t i = 0; i < numPoints; i++) {
    res[i] = fabs(valuesCT[i]);
  }
}

void BsplineStochasticCollocation::sgCoefficientCharacteristics(sgpp::base::DataVector& min,
                                                                sgpp::base::DataVector& max,
                                                                sgpp::base::DataVector& l2norm,
                                                                size_t maxLevel) {
  min.resize(0);
  max.resize(0);
  l2norm.resize(0);

  // compute hierarchical sparse grid interpolant
  if (hierarchicalTransformationFlag == false) {
    transformToHierarchical();
  }
  sgpp::base::GridStorage& gridStorage = hierarchicalGrid->getStorage();

  // sort coefficients by their levelsum
  // level enumeration differs a little for hierarchical Sparse Grids and combigrid module grids
  // becasue of the level 0 and level 1 interchange. Therefore the levelsum here is not equal to the
  // actual combigrid level.
  std::map<size_t, sgpp::base::DataVector> coeffMap_levelsum;
  std::map<std::vector<size_t>, sgpp::base::DataVector> coeffMap_level;
  std::map<size_t, size_t> numPointsPerLevel;
  size_t numGP = gridStorage.getSize();
  size_t numDim = gridStorage.getDimension();
  for (size_t p = 0; p < numGP; p++) {
    size_t levelsum = gridStorage[p].getLevelSum();
    // this is a hack! the level 0/1 mismatch between hierarchical SG and combination technique
    // leads to high levelsums in hierarchical SG where level 0 is sometimes called level 1
    // ToDo(rehmemk) recognize maxLevel automatically / fix mismatch in levelsums
    if (levelsum > maxLevel) {
      levelsum = 0;
      for (size_t d = 0; d < numDim; d++) {
        size_t temp = gridStorage[p].getLevel(d);
        if (temp != 1) {
          levelsum += temp;
        }
      }
    }
    coeffMap_levelsum[levelsum].push_back(hierarchicalCoefficients[p]);
    numPointsPerLevel[levelsum]++;

    std::vector<size_t> level(0);
    for (size_t d = 0; d < numDim; d++) {
      level.push_back(gridStorage[p].getLevel(d));
    }
    coeffMap_level[level].push_back(hierarchicalCoefficients[p]);
  }
  for (auto const& it : coeffMap_levelsum) {
    // print number of points per level
    // std::cout << it.first << ": " << numPointsPerLevel[it.first] << std::endl;
    min.push_back(it.second.min());
    max.push_back(it.second.max());
    l2norm.push_back(it.second.l2Norm());
  }
  // print level structure
  //  std::cout << "BSC: Levels:" << std::endl;
  //  for (auto const& it : coeffMap_level) {
  //    for (size_t d = 0; d < numDim; d++) {
  //      std::cout << it.first[d] << " ";
  //    }
  //    std::cout << "\n";
  //  }
}

sgpp::base::DataMatrix BsplineStochasticCollocation::getHierarchicalGridPoints() {
  sgpp::base::GridStorage& gridStorage = hierarchicalGrid->getStorage();
  sgpp::base::DataMatrix points(gridStorage.getDimension(), gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
      p[j] = gp.getStandardCoordinate(j);
    }
    points.setColumn(i, p);
  }
  return points;
}

double BsplineStochasticCollocation::evalHierarchical(sgpp::base::DataVector& x) {
  sgpp::optimization::InterpolantScalarFunction sparseGridSurrogate(*hierarchicalGrid,
                                                                    hierarchicalCoefficients);
  return sparseGridSurrogate.eval(x);
}

}  // namespace combigrid
}  // namespace sgpp
