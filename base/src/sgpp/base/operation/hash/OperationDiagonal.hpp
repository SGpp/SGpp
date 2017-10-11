// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONDIAGONAL_HPP
#define OPERATIONDIAGONAL_HPP

#include <sgpp/base/operation/hash/OperationMatrix.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace base {

/**
 * Implementation of diagonal Operation for all kinds of grids.
 * Corresponds to a Tikhonov matrix of \f$\boldsymbol{\Gamma}_{i,i} = c^{\vert \mathbf{l} \vert_1 - d}\f$
 * which results in the Gaussian prior
 * \f$ \boldsymbol{\alpha} \sim \mathcal{N} (0, \boldsymbol{\Gamma}^{-1}) \f$
 * on the weights \f$ \boldsymbol{\alpha} \f$.
 * The resulting prior for a 2d-grid with level four looks like
 * \image html diagonalRegularization.svg
 */
class OperationDiagonal : public OperationMatrix {
 public:
   /**
   * @brief Constructor of OperationDiagonal
   * @param gridStorage the grid storage
   * @param priorBase is the exponent base of the variance.
   * Corresponds to the reciprocal of \f$ c \f$ above.
   */
  explicit OperationDiagonal(sgpp::base::GridStorage* gridStorage, double priorBase = 0.25)
      : gridStorage(gridStorage), exponentBase(1.0/priorBase) {
    multiplicators = DataVector(0);
  }

  /**
   * Destructor
   */
  ~OperationDiagonal() override {}

  void mult(DataVector& alpha, DataVector& result) override {
    const auto size = alpha.getSize();
    // Reuse multiplicators if grid size hasn't changed!
    if (size != multiplicators.getSize()) {
      calculateMultiplicators(alpha);
    }
    result.resize(size);
#pragma omp parallel for
    for (size_t i = 0; i < size; ++i) {
      result[i] = multiplicators[i] * alpha[i];
    }
  }

 private:
  GridStorage* gridStorage;
  const double exponentBase;
  DataVector multiplicators;

  void calculateMultiplicators(DataVector& alpha) {
    const auto size = alpha.getSize();
    multiplicators.resize(size);
    size_t dimensions = gridStorage->getDimension();
    for (size_t i = 0; i < size; ++i) {
      auto gridIndex = gridStorage->getPoint(i);
      auto levelSum = gridIndex.getLevelSum();
      const double exponent = (static_cast<double>(levelSum) - static_cast<double>(dimensions));
      const double multiplicator = std::pow(exponentBase, exponent);
      multiplicators[i] = multiplicator;
    }
  }
};

}  // namespace base
}  // namespace sgpp
#endif /* OPERATIONDIAGONAL_HPP */
