// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
namespace sgpp {
namespace base {

/**
 * One component of a vector-valued function.
 *
 * Given a vector-valued \f$f\colon [0, 1]^d \to \mathbb{R}^m\f$
 * and indices \f$1 \le i_1 < \dotsb < i_n \le d\f$,
 * this class represents a new function
 * \f$g\colon [0, 1]^n \to \mathbb{R}\f$ with
 * \f$g(y_{i_1}, \dotsc, y_{i_n}) := f_k(y_1, \dotsc, y_d)\f$,
 * where \f$y_t\f$ is constant for \f$t \in \{i_1, \dotsc, i_n\}\f$.
 *
 * The resulting function \f$g\f$ is similar to a "slice plot" of
 * the component \f$f_k\f$ of \f$f\f$.
 */
class ComponentScalarFunction : public ScalarFunction {
 public:
  /**
   * Constructor.
   *
   * Use it like this:
   * ComponentScalarFunction g(f, {NAN, NAN, 0.42});
   * where f is a scalar-valued function with 3 parameters.
   * This selects the first two parameters of f, while constantly
   * using 0.42 for the third parameter.
   *
   * @param f             scalar-valued function
   * @param defaultValues Vector of constant default values.
   *                      It can be either empty (the default) or
   *                      a vector of exactly m doubles,
   *                      each of which can be finite or NAN.
   *                      If the vector is empty, it will be initialized
   *                      as m NANs (i.e., no restriction of the
   *                      parameter domain).
   *                      Each NAN represents a free parameter \f$x_t\f$,
   *                      while the finite entries denote the constant
   *                      values for the corresponding parameter.
   */
  ComponentScalarFunction(ScalarFunction& f,
                          std::vector<double> defaultValues = std::vector<double>())
      :

        ScalarFunction((defaultValues.size() > 0)
                           ? std::count(defaultValues.begin(), defaultValues.end(), NAN)
                           : f.getNumberOfParameters()),
        fScalar(&f),
        fVector(nullptr),
        dF(f.getNumberOfParameters()),
        k(0),
        defaultValues((defaultValues.size() > 0) ? defaultValues : std::vector<double>(dF, NAN)),
        tmpVec1(dF),
        tmpVec2(0) {
    initialize();
  }

  /**
   * Constructor.
   *
   * Use it like this:
   * ComponentScalarFunction g(f, 3, {NAN, NAN, 0.42});
   * where f is a vector-valued function with 5 components and
   * 3 parameters.
   * This selects the first two parameters and the fourth component
   * of f, while constantly using 0.42 for the third parameter.
   *
   * @param f             vector-valued function (m components)
   * @param k             index of component \f$f_k\f$ to select
   *                      (between 0 and m - 1)
   * @param defaultValues see other constructor
   */
  ComponentScalarFunction(VectorFunction& f, size_t k,
                          std::vector<double> defaultValues = std::vector<double>())
      :

        ScalarFunction((defaultValues.size() > 0)
                           ? std::count(defaultValues.begin(), defaultValues.end(), NAN)
                           : f.getNumberOfParameters()),
        fScalar(nullptr),
        fVector(&f),
        dF(f.getNumberOfParameters()),
        k(k),
        defaultValues((defaultValues.size() > 0) ? defaultValues : std::vector<double>(dF, NAN)),
        tmpVec1(dF),
        tmpVec2(f.getNumberOfComponents()) {
    initialize();
  }

  /**
   * Destructor.
   */
  ~ComponentScalarFunction() override {}

  /**
   * @param x     evaluation point \f$\vec{x} \in [0, 1]^n\f$
   * @return      \f$g(\vec{x}) := f_k(y_1, \dotsc, y_d)\f$
   *              where \f$(x_1, \dotsc, x_n) =
   *              (y_{i_1}, \dotsc, y_{i_n})\f$
   */
  inline double eval(const base::DataVector& x) override {
    size_t t2 = 0;

    // select entries of x which correspond to NAN entries in
    // defaultValues
    for (size_t t = 0; t < dF; t++) {
      if (std::isnan(defaultValues[t])) {
        tmpVec1[t] = x[t2];
        t2++;
      }
    }

    if (fScalar == nullptr) {
      // evaluate and select component
      fVector->eval(tmpVec1, tmpVec2);
      return tmpVec2[k];
    } else {
      // evaluate
      return fScalar->eval(tmpVec1);
    }
  }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<ScalarFunction>& clone) const override {
    clone = std::unique_ptr<ScalarFunction>(new ComponentScalarFunction(*this));
  }

 protected:
  /// scalar-valued function
  ScalarFunction* fScalar;
  /// vector-valued function
  VectorFunction* fVector;
  /// dimension of underlying function
  size_t dF;
  /// index of component
  size_t k;
  /// vector of default values, indicating free variables with NAN
  std::vector<double> defaultValues;
  /// temporary vector 1
  base::DataVector tmpVec1;
  /// temporary vector 2
  base::DataVector tmpVec2;

  void initialize() {
    // make sure defaultValues has the correct size
    if (defaultValues.size() != dF) {
      throw std::runtime_error("ComponentScalarFunction::initialize(): Invalid defaultValues.");
    }

    // initialize constant non-NAN entries
    for (size_t t = 0; t < dF; t++) {
      if (!std::isnan(defaultValues[t])) {
        tmpVec1[t] = defaultValues[t];
      }
    }
  }
};
}  // namespace base
}  // namespace sgpp
