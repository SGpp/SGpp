// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "HierarchicalStochasticCollocation.hpp"

namespace sgpp {
namespace combigrid {

HierarchicalStochasticCollocation::HierarchicalStochasticCollocation(
    std::shared_ptr<sgpp::base::Grid> grid, sgpp::combigrid::MultiFunction objectiveFunction,
    WeightFunctionsCollection weightFunctions, sgpp::base::DataVector bounds)
    : grid(grid),
      gridType(grid->getType()),
      objectiveFunction(objectiveFunction),
      weightFunctions(weightFunctions),
      bounds(bounds),
      currentNumGridPoints(grid->getSize()) {
  calculateCoefficients();
}

HierarchicalStochasticCollocation::HierarchicalStochasticCollocation(
    sgpp::base::GridType gridType, size_t dim, sgpp::combigrid::MultiFunction objectiveFunction,
    WeightFunctionsCollection weightFunctions, sgpp::base::DataVector bounds, size_t degree)
    : gridType(gridType),
      objectiveFunction(objectiveFunction),
      weightFunctions(weightFunctions),
      bounds(bounds),
      currentNumGridPoints(0) {
  if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
    grid = std::make_shared<sgpp::base::NakBsplineBoundaryGrid>(dim, degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineBoundaryCombigrid) {
    grid = std::make_shared<sgpp::base::NakBsplineBoundaryCombigridGrid>(dim, degree);
  } else if (gridType == sgpp::base::GridType::NotAKnotBsplineModified) {
    grid = std::make_shared<sgpp::base::NotAKnotBsplineModifiedGrid>(dim, degree);
  } else {
    std::cerr << "HierarchicalStochasticCollocation: grid type currently not supported"
              << std ::endl;
  }
  calculateCoefficients();
}

HierarchicalStochasticCollocation::HierarchicalStochasticCollocation(
    std::shared_ptr<sgpp::base::Grid> grid, sgpp::base::DataVector coefficients,
    WeightFunctionsCollection weightFunctions, sgpp::base::DataVector bounds)
    : grid(grid),
      gridType(grid->getType()),
      coefficients(coefficients),
      weightFunctions(weightFunctions),
      bounds(bounds),
      currentNumGridPoints(grid->getSize()) {}

HierarchicalStochasticCollocation::HierarchicalStochasticCollocation(
    const HierarchicalStochasticCollocation& other)
    : grid(other.grid),
      gridType(other.grid->getType()),
      coefficients(other.coefficients),
      weightFunctions(other.weightFunctions),
      bounds(other.bounds),
      currentNumGridPoints(other.grid->getSize()) {}

HierarchicalStochasticCollocation::~HierarchicalStochasticCollocation() {}

void HierarchicalStochasticCollocation::refineFull(size_t level) {
  grid->getGenerator().full(level);
  updateStatus();
}

void HierarchicalStochasticCollocation::refineRegular(size_t level) {
  grid->getGenerator().regular(level);
  updateStatus();
}

void HierarchicalStochasticCollocation::refineSurplusAdaptive(size_t refinementsNum) {
  sgpp::base::SurplusRefinementFunctor functor(coefficients, refinementsNum);
  grid->getGenerator().refine(functor);
  updateStatus();
}

void HierarchicalStochasticCollocation::refineSurplusAdaptiveByNumGridPoints(
    size_t maxNumGridPoints, size_t refinementsNum) {
  while (currentNumGridPoints < maxNumGridPoints) {
    refineSurplusAdaptive(refinementsNum);
  }
}

void HierarchicalStochasticCollocation::refineSurplusVolumeAdaptive(size_t refinementsNum) {
  sgpp::base::SurplusVolumeRefinementFunctor functor(coefficients, refinementsNum);
  grid->getGenerator().refine(functor);
  updateStatus();
}

void HierarchicalStochasticCollocation::refineSurplusVolumeAdaptiveByNumGridPoints(
    size_t maxNumGridPoints, size_t refinementsNum) {
  while (currentNumGridPoints < maxNumGridPoints) {
    refineSurplusVolumeAdaptive(refinementsNum);
  }
}

bool HierarchicalStochasticCollocation::updateStatus() {
  if (currentNumGridPoints < grid->getSize()) {
    currentNumGridPoints = grid->getSize();
    computedMeanFlag = false;
    computedVarianceFlag = false;
    computedDiscreteMeanFlag = false;
    computedDiscreteVarianceFlag = false;
    calculateCoefficients();
    return true;
  } else {
    return false;
  }
}

void HierarchicalStochasticCollocation::calculateCoefficients() {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  sgpp::base::DataVector f_values(gridStorage.getSize(), 0.0);
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
      p[j] = gp.getStandardCoordinate(j);
      //      std::cout << p[j] << " ";
    }
    //    std::cout << " <-p, f->";
    f_values[i] = objectiveFunction(p);
  }
  //  std::cout << f_values.toString() << std::endl;

  sgpp::optimization::sle_solver::Auto sleSolver;
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  sgpp::base::DataVector alpha(grid->getSize());
  if (!sleSolver.solve(hierSLE, f_values, alpha)) {
    std::cout << "Solving failed!" << std::endl;
  }
  coefficients = alpha;
}

void HierarchicalStochasticCollocation::createGridFromCombiLevelStructure(
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> levelStructure,
    std::shared_ptr<AbstractCombigridStorage> coefficientStorage) {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  size_t degree = 3;
  if (gridType == sgpp::base::GridType::NakBsplineBoundaryCombigrid) {
    degree = dynamic_cast<base::NakBsplineBoundaryCombigridGrid*>(grid.get())->getDegree();
    convertexpUniformBoundaryCombigridToHierarchicalSparseGrid(levelStructure, gridStorage);
  } else if (gridType == sgpp::base::GridType::NotAKnotBsplineModified) {
    degree = dynamic_cast<base::NotAKnotBsplineModifiedGrid*>(grid.get())->getDegree();
    convertCombigridToHierarchicalSparseGrid(levelStructure, gridStorage);
  } else {
    std::cerr << "HierarchicalStochasticCollocation: grid type currently not supported"
              << std::endl;
  }
  std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation =
      CombigridMultiOperation::createBsplineLinearCoefficientOperation(degree, grid->getDimension(),
                                                                       coefficientStorage);

  coefficients = calculateInterpolationCoefficientsForConvertedCombigird(
      grid, gridStorage, combigridMultiOperation, levelStructure);

  currentNumGridPoints = grid->getSize();
}

void HierarchicalStochasticCollocation::eval(sgpp::base::DataMatrix& xs,
                                             sgpp::base::DataVector& res) {
  res.resizeZero(xs.getNcols());
  sgpp::optimization::InterpolantScalarFunction sparseGridSurrogate(*grid, coefficients);
  for (size_t i = 0; i < xs.getNcols(); i++) {
    sgpp::base::DataVector x(grid->getDimension());
    xs.getColumn(i, x);
    res[i] = sparseGridSurrogate.eval(x);
  }
}

double HierarchicalStochasticCollocation::eval(sgpp::base::DataVector& x) {
  sgpp::optimization::InterpolantScalarFunction sparseGridSurrogate(*grid, coefficients);
  return sparseGridSurrogate.eval(x);
}

double HierarchicalStochasticCollocation::computeMean() {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  double mean = 0.0;
  if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
    size_t degree = dynamic_cast<base::NakBsplineBoundaryGrid*>(grid.get())->getDegree();
    OperationWeightedQuadratureNakBsplineBoundary wopQ(gridStorage, degree, weightFunctions,
                                                       bounds);
    mean = wopQ.doQuadrature(coefficients);
  } else if (gridType == sgpp::base::GridType::NakBsplineBoundaryCombigrid) {
    size_t degree = dynamic_cast<base::NakBsplineBoundaryCombigridGrid*>(grid.get())->getDegree();
    OperationWeightedQuadratureNakBsplineBoundaryCombigrid wopQ(gridStorage, degree,
                                                                weightFunctions, bounds);
    mean = wopQ.doQuadrature(coefficients);
  } else if (gridType == sgpp::base::GridType::NotAKnotBsplineModified) {
    size_t degree = dynamic_cast<base::NotAKnotBsplineModifiedGrid*>(grid.get())->getDegree();
    OperationWeightedQuadratureNotAKnotBsplineModified wopQ(gridStorage, degree, weightFunctions,
                                                            bounds);
    mean = wopQ.doQuadrature(coefficients);
  } else {
    std::cerr << "HierarchicalStochasticCollocation: grid type currently not supported"
              << std::endl;
  }

  return mean;
}

double HierarchicalStochasticCollocation::computeDiscreteMean(
    sgpp::base::DataMatrix discretePoints) {
  //  if (discretePoints.getNrows() != numDims) {
  //    throw sgpp::base::data_exception(
  //        "BsplineStochasticCollocation::computeDiscreteMean : Dimensions do not match");
  //  }
  sgpp::base::DataVector evaluations(discretePoints.getNcols());
  eval(discretePoints, evaluations);

  return evaluations.sum() / static_cast<double>(discretePoints.getNcols());
}

double HierarchicalStochasticCollocation::mean() {
  updateStatus();
  if (!computedMeanFlag) {
    ev = computeMean();
    computedMeanFlag = true;
  }
  return ev;
}

double HierarchicalStochasticCollocation::discreteMean(sgpp::base::DataMatrix discretePoints) {
  updateStatus();
  if (!computedDiscreteMeanFlag) {
    discreteEV = computeDiscreteMean(discretePoints);
    computedDiscreteMeanFlag = true;
  }
  return discreteEV;
}

double HierarchicalStochasticCollocation::computeVariance() {
  if (!computedMeanFlag) {
    mean();
  }

  sgpp::base::Grid* gridptr = grid.get();
  sgpp::base::DataVector product(coefficients.getSize());

  double variance = 0;
  if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
    LTwoScalarProductNakBsplineBoundary scalarProducts(gridptr, weightFunctions, bounds);
    size_t degree = dynamic_cast<base::NakBsplineBoundaryGrid*>(grid.get())->getDegree();
    if (degree == 1) {
      // calculate V(u) = E(u^2) - E(u)^2
      // this works for all B spline degrees but may be instable
      //(e.g. ev*ev might be larger than meanSquare => negative variance)
      // It also does not work if the weight function is not a probability density function
      // because then the algebraic formula for the variance does not hold
      scalarProducts.mult(coefficients, product);
      double meanSquare = product.dotProduct(coefficients);
      variance = meanSquare - ev * ev;
    } else {
      // calculate V(u) = E((u-E(u))^2)
      // this is done by subtracting E(u) from the coefficient of the constant function
      // it does not work for B spline degree 1 because there is no constant function in the basis
      // (We could add a constant basis function on level 0)
      sgpp::base::DataVector copyOfCoefficients(coefficients);
      copyOfCoefficients[0] -= ev;
      scalarProducts.mult(copyOfCoefficients, product);
      variance = product.dotProduct(copyOfCoefficients);
    }
  }
  if (gridType == sgpp::base::GridType::NakBsplineBoundaryCombigrid) {
    LTwoScalarProductNakBsplineBoundaryCombigrid scalarProducts(gridptr, weightFunctions, bounds);
    size_t degree = dynamic_cast<base::NakBsplineBoundaryCombigridGrid*>(grid.get())->getDegree();
    if (degree == 1) {
      scalarProducts.mult(coefficients, product);
      double meanSquare = product.dotProduct(coefficients);
      variance = meanSquare - ev * ev;
    } else {
      sgpp::base::DataVector copyOfCoefficients(coefficients);
      copyOfCoefficients[0] -= ev;
      scalarProducts.mult(copyOfCoefficients, product);
      variance = product.dotProduct(copyOfCoefficients);
    }
  } else if (gridType == sgpp::base::GridType::NotAKnotBsplineModified) {
    LTwoScalarProductNotAKnotBsplineModified scalarProducts(gridptr, weightFunctions, bounds);
    size_t degree = dynamic_cast<base::NotAKnotBsplineModifiedGrid*>(grid.get())->getDegree();
    if (degree == 1) {
      scalarProducts.mult(coefficients, product);
      double meanSquare = product.dotProduct(coefficients);
      variance = meanSquare - ev * ev;
    } else {
      sgpp::base::DataVector copyOfCoefficients(coefficients);
      copyOfCoefficients[0] -= ev;
      scalarProducts.mult(copyOfCoefficients, product);
      variance = product.dotProduct(copyOfCoefficients);
    }
  }
  return variance;
}

double HierarchicalStochasticCollocation::computeDiscreteVariance(
    sgpp::base::DataMatrix discretePoints) {
  //  if (discretePoints.getNrows() != numDims) {
  //    throw sgpp::base::data_exception(
  //        "BsplineStochasticCollocation::computeDiscreteVariance : Dimensions do not match");
  //  }
  if (!computedDiscreteMeanFlag) {
    discreteMean(discretePoints);
  }
  sgpp::base::DataVector evaluations(discretePoints.getNcols());
  eval(discretePoints, evaluations);

  discreteVar = 0.0;
  for (size_t i = 0; i < evaluations.getSize(); i++) {
    discreteVar += std::pow(evaluations[i] - discreteEV, 2);
  }
  discreteVar /= static_cast<double>(discretePoints.getNcols() - 1);

  return discreteVar;
}

double HierarchicalStochasticCollocation::variance() {
  updateStatus();
  if (!computedVarianceFlag) {
    var = computeVariance();
    computedVarianceFlag = true;
  }
  return var;
}

double HierarchicalStochasticCollocation::discreteVariance(sgpp::base::DataMatrix discretePoints) {
  updateStatus();
  if (!computedDiscreteVarianceFlag) {
    discreteVar = computeDiscreteVariance(discretePoints);
    computedDiscreteVarianceFlag = true;
  }
  return discreteVar;
}

void HierarchicalStochasticCollocation::coarsen(size_t removements_num, double threshold,
                                                bool recalculateCoefficients) {
  sgpp::base::HashCoarsening coarsen;
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  sgpp::base::SurplusCoarseningFunctor functor(coefficients, removements_num, threshold);
  if (recalculateCoefficients) {
    std::shared_ptr<sgpp::base::Grid> oldGrid(grid);
    sgpp::base::DataVector oldCoefficients(coefficients);
    coarsen.free_coarsen(gridStorage, functor, coefficients);
    coefficients =
        recalculateInterpolationCoefficientsAfterCoarsening(oldGrid, grid, oldCoefficients);
  } else {
    coarsen.free_coarsen(gridStorage, functor, coefficients);
  }
  computedMeanFlag = false;
  computedVarianceFlag = false;
}

void HierarchicalStochasticCollocation::coarsenIteratively(double maxVarDiff, double threshold,
                                                           double removements_percentage,
                                                           bool recalculateCoefficients) {
  if ((removements_percentage < 0) || (removements_percentage > 100)) {
    removements_percentage = 10;
  }

  size_t removements_num = static_cast<size_t>(
      floor(static_cast<double>(grid->getSize()) * removements_percentage / 100));
  double oldVariance = var;

  // coarsen
  sgpp::base::HashCoarsening coarsen;
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  sgpp::base::SurplusCoarseningFunctor functor(coefficients, removements_num, threshold);
  if (recalculateCoefficients) {
    std::shared_ptr<sgpp::base::Grid> oldGrid(grid);
    sgpp::base::DataVector oldCoefficients(coefficients);
    coarsen.free_coarsen(gridStorage, functor, coefficients);
    coefficients =
        recalculateInterpolationCoefficientsAfterCoarsening(oldGrid, grid, oldCoefficients);
  } else {
    coarsen.free_coarsen(gridStorage, functor, coefficients);
  }
  computedMeanFlag = false;
  computedVarianceFlag = false;

  // check for abort criterion and coarsen further if it's not met
  var = computeVariance();
  double varDiff = fabs(var - oldVariance);
  size_t numDeletedPoints = coarsen.getDeletedPoints().size();

  if ((varDiff <= maxVarDiff) && (numDeletedPoints != 0)) {
    coarsenIteratively(maxVarDiff, threshold, removements_percentage, recalculateCoefficients);
  }
}

size_t HierarchicalStochasticCollocation::numGridPoints() { return currentNumGridPoints; }

sgpp::base::DataMatrix HierarchicalStochasticCollocation::getHierarchicalGridPoints() {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  sgpp::base::DataMatrix points(gridStorage.getSize(), gridStorage.getDimension());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
      p[j] = gp.getStandardCoordinate(j);
    }
    points.setRow(i, p);
  }
  return points;
}

}  // namespace combigrid
}  // namespace sgpp
