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

template <typename Scalar>
class ScalarVector {
  Scalar val;

 public:
  ScalarVector() : val(Scalar()) {}

  ScalarVector(Scalar const &value) : val(value) {}

  ScalarVector(ScalarVector const &other) : val(other.val) {}

  ScalarVector<Scalar> &operator=(ScalarVector<Scalar> const &other) {
    val = other.val;
    return *this;
  }

  Scalar const &value() const { return val; }

  Scalar &value() { return val; }

  Scalar getValue() const { return val; }

  void add(ScalarVector<Scalar> const &other) { val += other.val; }

  void sub(ScalarVector<Scalar> const &other) { val -= other.val; }

  void componentwiseMult(ScalarVector<Scalar> const &other) { val *= other.val; }

  void scalarMult(Scalar const &factor) { val *= factor; }

  Scalar norm() const { return std::abs(val); }

  static ScalarVector<Scalar> zero() { return ScalarVector(Scalar(0)); }

  static ScalarVector<Scalar> one() { return ScalarVector(Scalar(1)); }
};

typedef ScalarVector<double> FloatScalarVector;

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_SCALARVECTOR_HPP_ */
