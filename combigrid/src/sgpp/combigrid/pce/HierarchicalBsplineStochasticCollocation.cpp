// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/pce/HierarchicalBsplineStochasticCollocation.hpp>

namespace sgpp {
namespace combigrid {

HierarchicalBsplineStochasticCollocation::HierarchicalBsplineStochasticCollocation(
    std::shared_ptr<sgpp::base::Grid> grid, size_t degree, sgpp::base::DataVector coefficients,
    WeightFunctionsCollection weightFunctions, sgpp::base::DataVector bounds)
    : grid(grid),
      degree(degree),
      coefficients(coefficients),
      weightFunctions(weightFunctions),
      bounds(bounds),
      currentNumGridPoints(grid->getSize()) {
  initialize(weightFunctions, bounds);
}

HierarchicalBsplineStochasticCollocation::~HierarchicalBsplineStochasticCollocation() {}

void HierarchicalBsplineStochasticCollocation::initialize(WeightFunctionsCollection weightFunctions,
                                                          sgpp::base::DataVector bounds) {
  scalarProducts.setWeightFunction(weightFunctions);
  scalarProducts.setBounds(bounds);
}

bool HierarchicalBsplineStochasticCollocation::updateStatus() {
  if (currentNumGridPoints < grid->getSize()) {
    currentNumGridPoints = grid->getSize();
    computedMeanFlag = false;
    computedVarianceFlag = false;
    return true;
  } else {
    return false;
  }
}

void HierarchicalBsplineStochasticCollocation::eval(sgpp::base::DataMatrix& xs,
                                                    sgpp::base::DataVector& res) {
  res.resizeZero(xs.getNrows());
  sgpp::optimization::InterpolantScalarFunction sparseGridSurrogate(*grid, coefficients);
  for (size_t i = 0; i < xs.getNrows(); i++) {
    sgpp::base::DataVector x;
    xs.getRow(i, x);
    res[i] = sparseGridSurrogate.eval(x);
  }
}

double HierarchicalBsplineStochasticCollocation::eval(sgpp::base::DataVector& x) {
  sgpp::optimization::InterpolantScalarFunction sparseGridSurrogate(*grid, coefficients);
  return sparseGridSurrogate.eval(x);
}

double HierarchicalBsplineStochasticCollocation::computeMean() {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  OperationWeightedQuadratureNakBsplineBoundary wopQ(gridStorage, degree, weightFunctions, bounds);
  double mean = wopQ.doQuadrature(coefficients);
  return mean;
}

double HierarchicalBsplineStochasticCollocation::mean() {
  updateStatus();
  if (!computedMeanFlag) {
    ev = computeMean();
    computedMeanFlag = true;
  }
  return ev;
}

double HierarchicalBsplineStochasticCollocation::computeVariance() {
  if (!computedMeanFlag) {
    mean();
  }

  sgpp::base::Grid* gridptr = grid.get();
  sgpp::base::DataVector product(coefficients.getSize());

  scalarProducts.updateGrid(gridptr);

  double variance = 0;
  if (degree == 1) {
    // calculate V(u) = E(u^2) - E(u)^2
    // this works for all B spline degrees but may be instable
    //(e.g. ev*ev might be larger than meanSquare => negative variance)
    // It also does not work if the weight function is not a probability density function because
    // then the algebraic formula for the variance does not hold
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

  return variance;
}

double HierarchicalBsplineStochasticCollocation::variance() {
  updateStatus();
  if (!computedVarianceFlag) {
    var = computeVariance();
    computedVarianceFlag = true;
  }
  return var;
}

void HierarchicalBsplineStochasticCollocation::coarsen(size_t removements_num, double threshold,
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

void HierarchicalBsplineStochasticCollocation::coarsenIteratively(double maxVarDiff,
                                                                  double threshold,
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
  std::cout << "BSC coarsening, removed " << numDeletedPoints
            << " points, variance difference: " << varDiff << std::endl;

  if ((varDiff <= maxVarDiff) && (numDeletedPoints != 0)) {
    coarsenIteratively(maxVarDiff, threshold, removements_percentage, recalculateCoefficients);
  }
}

size_t HierarchicalBsplineStochasticCollocation::numGridPoints() { return currentNumGridPoints; }

// void HierarchicalBsplineStochasticCollocation::sgCoefficientCharacteristics(
//    sgpp::base::DataVector& min, sgpp::base::DataVector& max, sgpp::base::DataVector& l2norm,
//    size_t maxLevel) {
//  min.resizeZero(0);
//  max.resizeZero(0);
//  l2norm.resizeZero(0);
//  sgpp::base::GridStorage& gridStorage = grid->getStorage();
//
//  // sort coefficients by their levelsum
//  // level enumeration differs a little for hierarchical Sparse Grids and combigrid module grids
//  // becasue of the level 0 and level 1 interchange. Therefore the levelsum here is not equal to
//  the
//  // actual combigrid level.
//  std::map<size_t, sgpp::base::DataVector> coeffMap_levelsum;
//  std::map<std::vector<size_t>, sgpp::base::DataVector> coeffMap_level;
//  std::map<size_t, size_t> numPointsPerLevel;
//  size_t numGP = gridStorage.getSize();
//  size_t numDim = gridStorage.getDimension();
//  for (size_t p = 0; p < numGP; p++) {
//    size_t levelsum = gridStorage[p].getLevelSum();
//    // this is a hack! the level 0/1 mismatch between hierarchical SG and combination technique
//    // leads to high levelsums in hierarchical SG where level 0 is sometimes called level 1
//    // ToDo(rehmemk) recognize maxLevel automatically / fix mismatch in levelsums
//    if (levelsum > maxLevel) {
//      levelsum = 0;
//      for (size_t d = 0; d < numDim; d++) {
//        size_t temp = gridStorage[p].getLevel(d);
//        if (temp != 1) {
//          levelsum += temp;
//        }
//      }
//    }
//    coeffMap_levelsum[levelsum].push_back(coefficients[p]);
//    numPointsPerLevel[levelsum]++;
//
//    std::vector<size_t> level(0);
//    for (size_t d = 0; d < numDim; d++) {
//      level.push_back(gridStorage[p].getLevel(d));
//    }
//    coeffMap_level[level].push_back(coefficients[p]);
//  }
//  for (auto const& it : coeffMap_levelsum) {
//    // print number of points per level
//    // std::cout << it.first << ": " << numPointsPerLevel[it.first] << std::endl;
//    min.push_back(it.second.min());
//    max.push_back(it.second.max());
//    l2norm.push_back(it.second.l2Norm());
//  }
//  // print level structure
//  //  std::cout << "BSC: Levels:" << std::endl;
//  //  for (auto const& it : coeffMap_level) {
//  //    for (size_t d = 0; d < numDim; d++) {
//  //      std::cout << it.first[d] << " ";
//  //    }
//  //    std::cout << "\n";
//  //  }
//}

sgpp::base::DataMatrix HierarchicalBsplineStochasticCollocation::getHierarchicalGridPoints() {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
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

}  // namespace combigrid
}  // namespace sgpp
