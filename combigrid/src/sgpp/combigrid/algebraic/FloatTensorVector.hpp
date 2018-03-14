// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>

// template <typename T>
// class TreeStorage;

namespace sgpp {
namespace combigrid {

class FloatTensorVector {
  // if d == 0, then values->getNumDimensions() == 1 (to store the only value at index 0)
  size_t d;
  std::shared_ptr<TreeStorage<FloatScalarVector>> values;

  void ensureDim(size_t dim);

 public:
  explicit FloatTensorVector(std::shared_ptr<TreeStorage<FloatScalarVector>> const &values);

  explicit FloatTensorVector(size_t d = 0);

  explicit FloatTensorVector(FloatScalarVector scalar);

  FloatTensorVector(FloatTensorVector const &other);

  FloatTensorVector &operator=(FloatTensorVector const &other);

  std::shared_ptr<TreeStorage<FloatScalarVector>> const &getValues() const { return values; }

  /**
   * This function can be called from python because it does not return a reference.
   */
  FloatScalarVector get(MultiIndex i);

  /**
   * Returns a reference to the i-th FloatScalarVector stored. If there are fewer elements stored,
   * the number of elements is extended via ensureMinimumSize().
   */
  FloatScalarVector &at(MultiIndex i);

  /**
   * Returns a reference to the i-th FloatScalarVector stored. If there are fewer elements stored,
   * the number of elements is extended via ensureMinimumSize().
   */
  FloatScalarVector &operator[](MultiIndex i);

  /**
   * This is an auxiliary function needed in FullGridQuadraticSummationStrategy.
   */
  FloatScalarVector operator[](size_t i);

  // size_t size() const { return values.size(); }

  void add(FloatTensorVector const &other);

  void sub(FloatTensorVector const &other);

  void componentwiseMult(FloatTensorVector const &other);

  void scalarMult(double const &factor);

  double norm() const;

  static FloatTensorVector zero() { return FloatTensorVector(FloatScalarVector(0.0)); }

  static FloatTensorVector one() { return FloatTensorVector(FloatScalarVector(1.0)); }
};

} /* namespace combigrid */
} /* namespace sgpp */
