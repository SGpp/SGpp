// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/combigrid/functions/ProbabilityDensityFunction1D.hpp>
#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
#include <sgpp/combigrid/operation/hierarchical/OperationWeightedQuadratureNakBsplineBoundary.hpp>
#include <sgpp/combigrid/operation/multidim/sparsegrid/LTwoScalarProductHashMapNakBsplineBoundary.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace sgpp {
namespace combigrid {
/**
 * This class holds a hierarchical B-spline surrogate stored as interpolation coefficients
 * and grid.
 * It can be used to evaluate the surrogate and to calculate mean and variance.
 * The hierarchical interpolant can be coarsened to reduce the number of grid points
 */
class HierarchicalBsplineStochasticCollocation {
 public:
  explicit HierarchicalBsplineStochasticCollocation(std::shared_ptr<sgpp::base::Grid> grid,
                                                    size_t degree,
                                                    sgpp::base::DataVector coefficients,
                                                    WeightFunctionsCollection weightFunctions,
                                                    sgpp::base::DataVector bounds);
  virtual ~HierarchicalBsplineStochasticCollocation();

  void eval(sgpp::base::DataMatrix& xs, sgpp::base::DataVector& res);
  double eval(sgpp::base::DataVector& x);

  double mean();
  double variance();

  /**
   * coarsen the hierarchical sparse grid surrogate
   * @param removements_num Number of grid points which should be removed (if possible - there could
   * be less removable grid points)
   * @param threshold The absolute value of the entries have to be greater or equal than the
   * threshold
   * @param recalculateCoefficients calculate new coefficients or use old (and after coarsening
   * not exact) coefficients)
   */
  void coarsen(size_t removements_num = 1, double threshold = 0.0,
               bool recalculateCoefficients = false);

  /**
   * Applies the coarsen routine and recalcualtes the variance until the difference between old and
   * new variance is larger than maxvarDiff or no more points are removed
   * @param maxVarDiff maximal difference of variance of the original and the coarsened surrogate
   * @param removements_percentage Maximal percentage of total grid points that shall be removed in
   * one coarsening step (if <0 or >100 default value 10 is used)
   * @param threshold The absolute value of the entries have to be greater or equal than the
   * threshold
   * @param recalculateCoefficients calculate new coefficients or use old (and after coarsening
   * not exact) coefficients)
   */
  void coarsenIteratively(double maxVarDiff, double threshold = 0.0,
                          double removements_percentage = 10, bool recalculateCoefficients = false);

  size_t numGridPoints();

  /**
   * returns characteristic values of the sparse grid coefficients, namely minimum, maximum and L2
   * norm per level. These values can be used to estimate the quality of the results if no
   * comparative solution is available
   * @param min vector of the minimal coefficient of each level
   * @param max vector of the maximum coefficient of each level
   * @param l2norm vector of the l2 norm of the coefficients of each level
   * @param maxLevel to get levelsums matching the combination technique levelsums the maximum level
   * must be forwarded
   */
  void sgCoefficientCharacteristics(sgpp::base::DataVector& min, sgpp::base::DataVector& max,
                                    sgpp::base::DataVector& l2norm, size_t maxLevel = 1000000);

  sgpp::base::DataMatrix getHierarchicalGridPoints();

 private:
  void initialize(WeightFunctionsCollection weightFunctions, sgpp::base::DataVector bounds);

  bool updateStatus();
  double computeMean();
  double computeVariance();

  std::shared_ptr<sgpp::base::Grid> grid;
  size_t degree;
  LTwoScalarProductHashMapNakBsplineBoundary scalarProducts;
  sgpp::base::DataVector coefficients;

  // pdf values
  sgpp::combigrid::WeightFunctionsCollection weightFunctions;
  sgpp::base::DataVector bounds;

  size_t currentNumGridPoints;

  // mean and variance storage
  bool computedMeanFlag = false;
  double ev = 0.0;
  bool computedVarianceFlag = false;
  double var = 0.0;
};

} /* namespace combigrid */
}  // namespace sgpp
