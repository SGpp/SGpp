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
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/grid/type/ModPolyGrid.hpp>
#include <sgpp/base/grid/type/PolyBoundaryGrid.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModPoly.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/combigrid/common/GridConversion.hpp>
#include <sgpp/combigrid/functions/ProbabilityDensityFunction1D.hpp>
#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/hierarchical/OperationWeightedQuadratureNakBsplineBoundary.hpp>
#include <sgpp/combigrid/operation/hierarchical/OperationWeightedQuadratureNakBsplineBoundaryCombigrid.hpp>
#include <sgpp/combigrid/operation/multidim/sparsegrid/LTwoScalarProductNakBsplineBoundary.hpp>
#include <sgpp/combigrid/operation/multidim/sparsegrid/LTwoScalarProductNakBsplineBoundaryCombigrid.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <algorithm>
#include <cmath>
#include <vector>
#include "../../../../../base/src/sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp"
#include "../operation/hierarchical/OperationWeightedQuadratureNakBsplineModified.hpp"
#include "../operation/multidim/sparsegrid/LTwoScalarProductNakBsplineModified.hpp"

namespace sgpp {
namespace combigrid {
/**
 * This class holds a hierarchical B-spline surrogate stored as interpolation coefficients
 * and grid.
 * It can be used to evaluate the surrogate and to calculate mean and variance.
 * The hierarchical interpolant can be coarsened to reduce the number of grid points
 */
class HierarchicalStochasticCollocation {
 public:
  explicit HierarchicalStochasticCollocation(std::shared_ptr<sgpp::base::Grid> grid,
                                             sgpp::combigrid::MultiFunction objectiveFunction,
                                             WeightFunctionsCollection weightFunctions,
                                             sgpp::base::DataVector bounds);

  explicit HierarchicalStochasticCollocation(sgpp::base::GridType gridType, size_t dim,
                                             sgpp::combigrid::MultiFunction objectiveFunction,
                                             WeightFunctionsCollection weightFunctions,
                                             sgpp::base::DataVector bounds, size_t degree = 3);

  explicit HierarchicalStochasticCollocation(std::shared_ptr<sgpp::base::Grid> grid,
                                             sgpp::base::DataVector coefficients,
                                             WeightFunctionsCollection weightFunctions,
                                             sgpp::base::DataVector bounds);

  // copy constructor
  // needs additionally either an objectiveFunction or coefficients
  HierarchicalStochasticCollocation(const HierarchicalStochasticCollocation& other);

  virtual ~HierarchicalStochasticCollocation();

  void refineFull(size_t level);
  /*@param level 		 level for inner points
   * @param boundaryLevel 1 + how much levels the boundary is coarser than
   *                      the main axes, 0 means one level finer,
   *                      1 means same level,
   *                      2 means one level coarser, etc.
   */
  void refineRegular(size_t level);
  void refineSurplusAdaptive(size_t refinements_num);
  /**
   * refines the grid depending on the surpluses.
   * @param maxNumGridPoints maximal size of the refined grid
   * @param refinementsNum   number of points refined in each refinement step
   */
  void refineSurplusAdaptiveByNumGridPoints(size_t maxNumGridPoints, size_t refinementsNum = 1);
  void refineSurplusVolumeAdaptive(size_t refinementsNum);
  void refineSurplusVolumeAdaptiveByNumGridPoints(size_t maxNumGridPoints, size_t refinementsNum);

  void createGridFromCombiLevelStructure(
      std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> levelStructure,
      std::shared_ptr<AbstractCombigridStorage> coefficientStorage);

  void eval(sgpp::base::DataMatrix& xs, sgpp::base::DataVector& res);
  double eval(sgpp::base::DataVector& x);

  double mean();
  /**
   * Calculate mean discrete using a set of points. Usually Monte Carlo realisations of the
   * probability densities
   * @param discretePoints the set of discrete points
   * @return mean (expecation value)
   */
  double discreteMean(sgpp::base::DataMatrix discretePoints);
  double variance();
  double discreteVariance(sgpp::base::DataMatrix discretePoints);

  void calculateCoefficients();

  sgpp::base::DataVector leastSquares(sgpp::base::DataMatrix points,
                                      sgpp::base::DataVector functionValues, size_t degree);
  /**
   * coarsen the hierarchical sparse grid surrogate
   * @param removements_num Number of grid points which should be removed (if possible - there
   * could be less removable grid points)
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

  void setObjectivefunction(sgpp::combigrid::MultiFunction objectiveFunction) {
    this->objectiveFunction = objectiveFunction;
  }

  sgpp::base::DataMatrix getHierarchicalGridPoints();

  sgpp::base::DataVector getCoefficients() { return coefficients; }

 private:
  void initialize(WeightFunctionsCollection weightFunctions, sgpp::base::DataVector bounds);

  bool updateStatus();
  double computeMean();
  double computeDiscreteMean(sgpp::base::DataMatrix discretePoints);
  double computeVariance();
  double computeDiscreteVariance(sgpp::base::DataMatrix discretePoints);

  std::shared_ptr<sgpp::base::Grid> grid;
  sgpp::base::GridType gridType;
  //  LTwoScalarProductNakBsplineBoundary scalarProducts;
  sgpp::base::DataVector coefficients;
  sgpp::combigrid::MultiFunction objectiveFunction;

  // pdf weightFunctions and the definition domain bounds.
  // 8the objective function itself is w.l.o.g. assumed to be defined on the unit hypercube)
  sgpp::combigrid::WeightFunctionsCollection weightFunctions;
  sgpp::base::DataVector bounds;

  size_t currentNumGridPoints;

  // mean and variance storage
  bool computedMeanFlag = false;
  double ev = 0.0;
  bool computedVarianceFlag = false;
  double var = 0.0;

  // discrete mean and variance
  bool computedDiscreteMeanFlag = false;
  double discreteEV = 0.0;
  double computedDiscreteVarianceFlag = false;
  double discreteVar = 0.0;
};

} /* namespace combigrid */
}  // namespace sgpp
