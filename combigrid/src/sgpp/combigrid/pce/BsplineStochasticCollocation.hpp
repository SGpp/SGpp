// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

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

  void getComponentSobolIndices(sgpp::base::DataVector& componentSsobolIndices,
                                bool normalized = true) override;
  void getTotalSobolIndices(sgpp::base::DataVector& totalSobolIndices,
                            bool normalized = true) override;

  void updateConfig(sgpp::combigrid::CombigridSurrogateModelConfiguration config) override;

  size_t numGridPoints() override;
  std::shared_ptr<LevelInfos> getInfoOnAddedLevels() override;

  /** calculate the difference between the combination technique surrogate u_{CT} and the
  * corresponding hierarchical sparse grid surrogate u_{SG}. This is close to zero because the
  * surrogates lie in the same space.
  * @param xs matrix of points p of which |u_{CT}(p) - u_{SG}(p)| is calculated
  * @param res returns the absolute differences
  * **/
  void differenceCTSG(sgpp::base::DataMatrix& xs, sgpp::base::DataVector& res);

 private:
  void initializeOperations(std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
                            std::shared_ptr<AbstractCombigridStorage> storage,
                            std::shared_ptr<LevelManager> levelManager);

  bool updateStatus();
  double computeMean();
  double computeVariance();

  void countPolynomialTerms();
  size_t additionalQuadraturePoints(OrthogonalPolynomialBasisType polyType);

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
  //  ToDo(rehmemk)
  size_t numthreads = 4;
  LTwoScalarProductHashMapNakBsplineBoundaryCombigrid scalarProducts;

  // operation
  std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation;
  std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation;
};

} /* namespace combigrid */
} /* namespace sgpp */
