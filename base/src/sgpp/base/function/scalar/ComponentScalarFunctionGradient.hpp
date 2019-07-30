// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/base/function/vector/VectorFunctionGradient.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
namespace sgpp {
namespace base {

/**
 * One component of a vector-valued function gradient.
 *
 * @see ComponentScalarFunction
 */
class ComponentScalarFunctionGradient : public ScalarFunctionGradient {
 public:
  /**
   * Constructor.
   *
   * Use it like this:
   * ComponentScalarFunctionGradient gGradient(fGradient, {NAN, NAN, 0.42});
   * where fGradient is a scalar-valued function gradient with
   * 3 parameters.
   * This selects the first two parameters of fGradient, while
   * constantly using 0.42 for the third parameter.
   *
   * @param fGradient     scalar-valued function gradient
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
  ComponentScalarFunctionGradient(ScalarFunctionGradient& fGradient,
                                  std::vector<double> defaultValues = std::vector<double>())
      :

        ScalarFunctionGradient((defaultValues.size() > 0)
                                   ? std::count(defaultValues.begin(), defaultValues.end(), NAN)
                                   : fGradient.getNumberOfParameters()),
        fGradientScalar(&fGradient),
        fGradientVector(nullptr),
        dF(fGradient.getNumberOfParameters()),
        k(0),
        defaultValues((defaultValues.size() > 0) ? defaultValues : std::vector<double>(dF, NAN)),
        tmpVec1(dF),
        tmpVec2(dF),
        tmpMat(0, 0) {
    initialize();
  }

  /**
   * Constructor.
   *
   * Use it like this:
   * ComponentScalarFunctionGradient gGradient(fGradient, 3, {NAN, NAN, 0.42});
   * where fGradient is a vector-valued function gradient with
   * 5 components and 3 parameters.
   * This selects the first two parameters and the fourth component
   * of fGradient, while constantly using 0.42 for the third parameter.
   *
   * @param fGradient     vector-valued function gradient (m components)
   * @param k             index of component \f$f_k\f$ to select
   *                      (between 0 and m - 1)
   * @param defaultValues see other constructor
   */
  ComponentScalarFunctionGradient(VectorFunctionGradient& fGradient, size_t k,
                                  std::vector<double> defaultValues = std::vector<double>())
      :

        ScalarFunctionGradient((defaultValues.size() > 0)
                                   ? std::count(defaultValues.begin(), defaultValues.end(), NAN)
                                   : fGradient.getNumberOfParameters()),
        fGradientScalar(nullptr),
        fGradientVector(&fGradient),
        dF(fGradient.getNumberOfParameters()),
        k(k),
        defaultValues((defaultValues.size() > 0) ? defaultValues : std::vector<double>(dF, NAN)),
        tmpVec1(dF),
        tmpVec2(fGradient.getNumberOfComponents()),
        tmpMat(fGradient.getNumberOfComponents(), dF) {
    initialize();
  }

  /**
   * Destructor.
   */
  ~ComponentScalarFunctionGradient() override {}

  /**
   * @param[in] x evaluation point \f$\vec{x} \in [0, 1]^n\f$
   * @param[out] gradient \f$\nabla_{\vec{x}} g(\vec{x})\f$
   * @return      \f$g(\vec{x}) := f_k(y_1, \dotsc, y_d)\f$
   *              where \f$(x_1, \dotsc, x_n) =
   *              (y_{i_1}, \dotsc, y_{i_n})\f$
   */
  inline double eval(const base::DataVector& x, base::DataVector& gradient) override {
    size_t t2 = 0;

    // select entries of x which correspond to NAN entries in
    // defaultValues
    for (size_t t = 0; t < dF; t++) {
      if (std::isnan(defaultValues[t])) {
        tmpVec1[t] = x[t2];
        t2++;
      }
    }

    if (fGradientScalar == nullptr) {
      // evaluate and select component
      fGradientVector->eval(tmpVec1, tmpVec2, tmpMat);
      t2 = 0;

      for (size_t t = 0; t < dF; t++) {
        if (std::isnan(defaultValues[t])) {
          gradient[t2] = tmpMat(k, t);
          t2++;
        }
      }

      return tmpVec2[k];
    } else {
      // evaluate
      const double fx = fGradientScalar->eval(tmpVec1, tmpVec2);
      t2 = 0;

      for (size_t t = 0; t < dF; t++) {
        if (std::isnan(defaultValues[t])) {
          gradient[t2] = tmpVec2[t];
          t2++;
        }
      }

      return fx;
    }
  }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const override {
    clone = std::unique_ptr<ScalarFunctionGradient>(new ComponentScalarFunctionGradient(*this));
  }

 protected:
  /// scalar-valued function gradient
  ScalarFunctionGradient* fGradientScalar;
  /// vector-valued function gradient
  VectorFunctionGradient* fGradientVector;
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
  /// temporary matrix
  base::DataMatrix tmpMat;

  void initialize() {
    // make sure defaultValues has the correct size
    if (defaultValues.size() != dF) {
      throw std::runtime_error(
          "ComponentScalarFunctionGradient::initialize(): Invalid defaultValues.");
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
