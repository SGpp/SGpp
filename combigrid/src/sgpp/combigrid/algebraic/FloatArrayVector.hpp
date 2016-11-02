// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_FLOATARRAYVECTOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_FLOATARRAYVECTOR_HPP_

#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/globaldef.hpp>

#include <cmath>

#include <algorithm>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * This class models a vector of scalars supporting operations such as addition, scalar
 * multiplication and componentwise muliplication. It is intended to have the same interface as
 * FloatScalarVector. Thus, one can use the same code for single-evaluation and multi-evaluation,
 * where the mode is controlled by passing FloatArrayVector or FloatScalarVector as a template
 * parameter.
 * It contains multiple values resulting for example from Interpolation with multiple evaluation
 * points.
 * If an operation on two FloatArrayVector elements with different number of elements is executed,
 * the last element of the vector with less elements is repeated until both of them have the same
 * size. This allows e.g. for a zero vector (with one element) which can adjust to vectors of
 * arbitrary size.
 */
class FloatArrayVector {
  std::vector<FloatScalarVector> values;

  /**
   * Ensures that the current vector has at least minSize elements.
   * If it has less elements, the last element is repeated until minSize is reached.
   */
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

  /**
   * This function can be called from python because it does not return a reference.
   * It does not do range checking
   */
  FloatScalarVector get(size_t i) const { return (*this)[i]; }

  /**
   * Returns a reference to the i-th FloatScalarVector stored. If there are fewer elements stored,
   * the number of elements is extended via ensureMinimumSize().
   */
  FloatScalarVector &at(size_t i) {
    ensureMinimumSize(i + 1);
    return values[i];
  }

  /**
   * Returns a reference to the i-th FloatScalarVector stored. If there are fewer elements stored,
   * the number of elements is extended via ensureMinimumSize().
   */
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

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_FLOATARRAYVECTOR_HPP_ */
