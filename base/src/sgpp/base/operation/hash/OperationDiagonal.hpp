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
 * Implementation of diagonal Operation for all kinds of grids
 */
class OperationDiagonal : public OperationMatrix {
 public:
  /**
   * Constructor of OperationDiagonal
   */
  explicit OperationDiagonal(sgpp::base::GridStorage* gridStorage,
                             double multiplicationFactor = 0.25)
      : gridStorage(gridStorage), multiplicationFactor(multiplicationFactor) {}

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
    result = DataVector(size);
    for (size_t i = 0; i < size; ++i) {
      result[i] = multiplicators[i] * alpha[i];
    }
  }

 private:
  GridStorage* gridStorage;
  const double multiplicationFactor;
  DataVector multiplicators;

  void calculateMultiplicators(DataVector& alpha) {
    const auto size = alpha.getSize();
    multiplicators = DataVector(size);
    size_t dimensions = gridStorage->getDimension();
    for (size_t i = 0; i < size; ++i) {
      sgpp::base::GridStorage::index_pointer gridIndex = gridStorage->get(i);
      sgpp::base::GridIndex::level_type levelSum = gridIndex->getLevelSum();
      const double exponent = (static_cast<double>(levelSum) - static_cast<double>(dimensions));
      const double multiplicator = std::pow(multiplicationFactor, exponent);
      multiplicators[i] = multiplicator;
    }
  }
};

}  // namespace base
}  // namespace sgpp
#endif /* OPERATIONDIAGONAL_HPP */
