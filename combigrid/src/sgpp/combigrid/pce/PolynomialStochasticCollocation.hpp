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

#include <map>
#include <vector>

namespace sgpp {
namespace combigrid {

class LinearTransformation {
 public:
  LinearTransformation() {}

  LinearTransformation(sgpp::base::DataVector& bounds) { initialize(bounds); }

  virtual ~LinearTransformation(){};

  void initialize(sgpp::base::DataVector& bounds) {
    if (bounds.getSize() % 2 != 0) {
      throw sgpp::base::algorithm_exception(
          "LinearTransformation::initialize: bounds have not the correct size");
    }
    size_t numDims = bounds.getSize() >> 1;
    widths.resize(numDims);
    xlower.resize(numDims);

    for (size_t i = 0; i < numDims; i++) {
      widths[i] = bounds[2 * i + 1] - bounds[2 * i];
      xlower[i] = bounds[2 * i];
    }
  }

  void unitToProbabilistic(sgpp::base::DataVector& x) {
    // linear transformation
    for (size_t i = 0; i < x.getSize(); i++) {
      x[i] = unitToProbabilistic(x[i], i);
    }
  }

  double unitToProbabilistic(double x, size_t idim) { return widths[idim] * x + xlower[idim]; }

  double vol() {
    double ans = 1.0;
    for (size_t i = 0; i < widths.getSize(); i++) {
      ans *= widths[i];
    }
    return ans;
  }

  double vol(size_t i) { return widths[i]; }

  double getLowerBound(size_t i) { return xlower[i]; }
  double getUpperBound(size_t i) { return xlower[i] + widths[i]; }

 private:
  sgpp::base::DataVector widths;
  sgpp::base::DataVector xlower;
};

class PolynomialStochasticCollocation {
 public:
  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation,
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis,
      sgpp::base::DataVector& bounds);

  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation,
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis,
      sgpp::base::DataVector& bounds);

  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation,
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis,
      sgpp::base::DataVector& bounds);

  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation,
      std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases,
      sgpp::base::DataVector& bounds);

  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation,
      std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases,
      sgpp::base::DataVector& bounds);

  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation,
      std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases,
      sgpp::base::DataVector& bounds);

#ifdef USE_DAKOTA
  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation,
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis);

  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation,
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis);

  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation,
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis);

  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation,
      std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases);

  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation,
      std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases);

  PolynomialStochasticCollocation(
      std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation,
      std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases);
#endif

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

  void initializeLinearTransformation();

  void countPolynomialTerms();
  size_t additionalQuadraturePoints(OrthogonalPolynomialBasisType polyType);

  double quad(sgpp::combigrid::MultiIndex i);
  double quad(sgpp::combigrid::MultiIndex i, sgpp::combigrid::MultiIndex j);

  void joinMultiIndices(MultiIndex& ix, MultiIndex& jx, MultiIndex& kx);

  size_t numDims;
  std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation;
  std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation;
  std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation;

  std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> legendreBasis;
  std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>> functionBases;

  LinearTransformation trans;
  size_t numGridPoints;
  sgpp::combigrid::FloatTensorVector expansionCoefficients;

  // lookup table for inner products
  std::map<MultiIndex, double> innerProducts;
};

} /* namespace combigrid */
} /* namespace sgpp */
