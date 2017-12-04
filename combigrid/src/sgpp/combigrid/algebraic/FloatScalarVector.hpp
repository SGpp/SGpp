// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_FLOATSCALARVECTOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_FLOATSCALARVECTOR_HPP_

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <iostream>

namespace sgpp {
namespace combigrid {

/**
 * This class models a single scalar supporting operations such as addition, scalar
 * multiplication and componentwise muliplication. It is intended to have the same interface as
 * FloatArrayVector. Thus, one can use the same code for single-evaluation and multi-evaluation,
 * where the mode is controlled by passing FloatArrayVector or FloatScalarVector as a template
 * parameter.
 */
class FloatScalarVector {
  double val;

 public:
  FloatScalarVector() : val(double()) {}

  explicit FloatScalarVector(double const &value) : val(value) {}

  FloatScalarVector(FloatScalarVector const &other) : val(other.val) {}

  FloatScalarVector &operator=(FloatScalarVector const &other) {
    val = other.val;
    return *this;
  }

  FloatScalarVector &operator=(double const &value) {
    val = value;
    return *this;
  }

  double const &value() const { return val; }

  double &value() { return val; }

  double getValue() const { return val; }

  void add(FloatScalarVector const &other) { val += other.val; }

  void sub(FloatScalarVector const &other) { val -= other.val; }

  void componentwiseMult(FloatScalarVector const &other) { val *= other.val; }

  void scalarMult(double const &factor) { val *= factor; }

  double norm() const { return std::fabs(val); }

  static FloatScalarVector zero() { return FloatScalarVector(0.0); }

  static FloatScalarVector one() { return FloatScalarVector(1.0); }

  // This is an auxiliary function needed in FullGridQuadraticSummationStrategy.
  FloatScalarVector const operator[](size_t i) const {
    if (i != 0) {
      std::cerr << "FloatScalarVector: operator[] is only defined for argument 0";
    }
    return *this;
  }
};

std::ostream &operator<<(std::ostream &stream, FloatScalarVector v);

std::istream &operator>>(std::istream &stream, FloatScalarVector &v);

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_FLOATSCALARVECTOR_HPP_ */
