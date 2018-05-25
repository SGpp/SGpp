// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/operation/multidim/sparsegrid/LTwoScalarProductHashMapNakBsplineBoundaryCombigrid.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {
/**
 * This class holds a B-spline combination technique surrogate stored as interpolation coefficients
 * and level structure (in the CombigridSurrogateModelConfiguration).
 * It can be used to evaluate the surrogate and to calculate mean and variance.
 * The mean is calcualted as a sparse grid quadrature.
 * The variance is calculated by transforming the level structure to a hierarchical sparse grid and
 * calculation of the neccessary scalar products there
 *
 * Both interpolants: the original one from the combination technique and the hierarchical one can
 * be accessed through this class.
 * The hierarchical one can be coarsened to reduce the number of grid points
 * After coarsening the two interpolants are no longer equivalent!
 */
class BsplineStochasticCollocation : public CombigridSurrogateModel {
 public:
  explicit BsplineStochasticCollocation(
      sgpp::combigrid::CombigridSurrogateModelConfiguration& config);
  virtual ~BsplineStochasticCollocation();

  void eval(sgpp::base::DataMatrix& xs, sgpp::base::DataVector& res) override;
  double eval(sgpp::base::DataVector& x) override;

  double mean() override;
  double variance() override;

  /**
   * Calculates the variance using the combination technique BSplineScalarProductEvaluator
   * This is used to compare runtimes between the hierarchical grid (which should be faster)
   * and the combination technique
   */
  double computeVarianceWithCombiTechnique();

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

  void getComponentSobolIndices(sgpp::base::DataVector& componentSsobolIndices,
                                bool normalized = true) override;
  void getTotalSobolIndices(sgpp::base::DataVector& totalSobolIndices,
                            bool normalized = true) override;

  void updateConfig(sgpp::combigrid::CombigridSurrogateModelConfiguration config) override;

  size_t numGridPoints() override;
  size_t numHierarchicalGridPoints();
  std::shared_ptr<LevelInfos> getInfoOnAddedLevels() override;

  /** calculate the difference between the combination technique surrogate u_{CT} and the
   * corresponding hierarchical sparse grid surrogate u_{SG}. This is close to zero because the
   * surrogates lie in the same space.
   * @param xs matrix of points p of which |u_{CT}(p) - u_{SG}(p)| is calculated
   * @param res returns the absolute differences
   * **/
  void differenceCTSG(sgpp::base::DataMatrix& xs, sgpp::base::DataVector& res);

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
  double evalHierarchical(sgpp::base::DataVector& x);

 private:
  void initializeOperations(std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
                            std::shared_ptr<AbstractCombigridStorage> storage,
                            std::shared_ptr<LevelManager> levelManager);

  bool updateStatus();
  double computeMean();
  double computeVariance();

  /**
   * interpolates the combination technique surrogate with hierarchical B-splines
   * The grid and coefficients are stored in ,hierarchicalGrid and hierarchicalCoefficients
   * respectively
   */
  void transformToHierarchical();

  void countPolynomialTerms();

  double quad(sgpp::combigrid::MultiIndex i);
  double quad(sgpp::combigrid::MultiIndex i, sgpp::combigrid::MultiIndex j);

  // pdf values
  sgpp::combigrid::WeightFunctionsCollection weightFunctions;
  bool customWeightFunction;

  size_t currentNumGridPoints;

  // mean and variance storage
  bool computedMeanFlag;
  double ev;
  bool computedVarianceFlag;
  double var;
  // basis coefficients for Bspline interpolation
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> coefficientStorage;
  LTwoScalarProductHashMapNakBsplineBoundaryCombigrid scalarProducts;

  // operation
  std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation;
  std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation;

  // hierarchical sparse grid interpolant
  bool hierarchicalTransformationFlag = false;
  std::shared_ptr<sgpp::base::Grid> hierarchicalGrid;
  sgpp::base::DataVector hierarchicalCoefficients;
};

} /* namespace combigrid */
}  // namespace sgpp
