// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_SCALARVECTOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_SCALARVECTOR_HPP_

#include <sgpp/globaldef.hpp>
#include <cmath>

namespace sgpp {
namespace combigrid {

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

  double norm() const { return std::abs(val); }

  static FloatScalarVector zero() { return FloatScalarVector(0.0); }

  static FloatScalarVector one() { return FloatScalarVector(1.0); }
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_SCALARVECTOR_HPP_ */
