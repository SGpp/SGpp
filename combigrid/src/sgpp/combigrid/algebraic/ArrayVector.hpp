// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_ARRAYVECTOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_ARRAYVECTOR_HPP_

#include <sgpp/globaldef.hpp>
#include <sgpp/combigrid/algebraic/ScalarVector.hpp>

#include <vector>
#include <algorithm>
#include <cmath>

namespace sgpp {
namespace combigrid {

class FloatArrayVector {
  std::vector<FloatScalarVector> values;

  void ensureMinimumSize(size_t minSize) {
    if (values.size() == 0) {
      values.push_back(FloatScalarVector::zero());
    }
    while (values.size() < minSize) {
      values.push_back(values.back());
    }
  }

 public:
  explicit FloatArrayVector(std::vector<FloatScalarVector> const &values) : values(values) {}

  explicit FloatArrayVector(FloatScalarVector value) : values(1, value) {}

  FloatArrayVector() : values() {}

  FloatArrayVector(FloatArrayVector const &other) : values(other.values) {}

  FloatArrayVector &operator=(FloatArrayVector const &other) {
    values = other.values;
    return *this;
  }

  std::vector<FloatScalarVector> const &getValues() const { return values; }

  FloatScalarVector get(size_t i) const { return (*this)[i]; }

  FloatScalarVector &at(size_t i) {
    ensureMinimumSize(i + 1);
    return values[i];
  }

  FloatScalarVector &operator[](size_t i) {
    ensureMinimumSize(i + 1);
    return values[i];
  }

  /**
   * This operator is unsafe because it does no range checking
   */
  FloatScalarVector const &operator[](size_t i) const { return values[i]; }

  size_t size() const { return values.size(); }

  void add(FloatArrayVector const &other) {
    ensureMinimumSize(other.size());

    size_t otherSize = other.size();
    for (size_t i = 0; i < otherSize; ++i) {
      values[i].add(other.values[i]);
    }

    for (size_t i = otherSize; i < size(); ++i) {
      values[i].add(other.values.back());
    }
  }

  void sub(FloatArrayVector const &other) {
    ensureMinimumSize(other.size());

    size_t otherSize = other.size();
    for (size_t i = 0; i < otherSize; ++i) {
      values[i].sub(other.values[i]);
    }

    for (size_t i = otherSize; i < size(); ++i) {
      values[i].sub(other.values.back());
    }
  }

  void componentwiseMult(FloatArrayVector const &other) {
    ensureMinimumSize(other.size());

    size_t otherSize = other.size();
    for (size_t i = 0; i < otherSize; ++i) {
      values[i].componentwiseMult(other.values[i]);
    }

    for (size_t i = otherSize; i < size(); ++i) {
      values[i].componentwiseMult(other.values.back());
    }
  }

  void scalarMult(double const &factor) {
    for (size_t i = 0; i < values.size(); ++i) {
      values[i].scalarMult(factor);
    }
  }

  double norm() const {
    double result = 0.0;

    for (auto &val : values) {
      double n = val.norm();
      result += n * n;
    }

    return std::sqrt(result);
  }

  static FloatArrayVector zero() { return FloatArrayVector(FloatScalarVector(0.0)); }

  static FloatArrayVector one() { return FloatArrayVector(FloatScalarVector(1.0)); }
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_ARRAYVECTOR_HPP_ */
