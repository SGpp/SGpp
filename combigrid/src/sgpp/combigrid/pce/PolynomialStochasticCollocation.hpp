// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class PolynomialStochasticCollocation {
 public:
  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation,
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis,
      sgpp::base::DataVector const& bounds = sgpp::base::DataVector(0));

  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation,
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis,
      sgpp::base::DataVector const& bounds = sgpp::base::DataVector(0));

  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation,
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis,
      sgpp::base::DataVector const& bounds = sgpp::base::DataVector(0));

  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation,
      std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases,
      sgpp::base::DataVector const& bounds = sgpp::base::DataVector(0));

  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation,
      std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases,
      sgpp::base::DataVector const& bounds = sgpp::base::DataVector(0));

  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation,
      std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases,
      sgpp::base::DataVector const& bounds = sgpp::base::DataVector(0));

  virtual ~PolynomialStochasticCollocation();

  double mean();
  double variance();

  std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> getCombigridTensorOperation();
  void updateOperation(std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation);
  void updateOperation(
      std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridOperation);
  void updateOperation(
      std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridOperation);

 private:
  void initializeTensorOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::shared_ptr<LevelManager> levelManager);

  void initializeBounds();
  void initializeWeightFunctions();

  bool updateStatus();
  double computeMean();
  double computeVariance();

  void countPolynomialTerms();
  size_t additionalQuadraturePoints(OrthogonalPolynomialBasisType polyType);

  double quad(sgpp::combigrid::MultiIndex i);
  double quad(sgpp::combigrid::MultiIndex i, sgpp::combigrid::MultiIndex j);

  // number of dimensions
  size_t numDims;

  // store the operations for combigrid evaluation
  std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation;
  std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation;
  std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation;

  // global polynomial basis
  std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> legendreBasis;
  // orthogonal basis for pdf values
  std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>> functionBases;
  std::vector<sgpp::combigrid::SingleFunction> weightFunctions;

  sgpp::base::DataVector bounds;
  size_t numGridPoints;
  sgpp::combigrid::FloatTensorVector expansionCoefficients;

  // mean and variance storage
  bool computedMeanFlag;
  double ev;
  bool computedVarianceFlag;
  double var;
};

} /* namespace combigrid */
} /* namespace sgpp */
