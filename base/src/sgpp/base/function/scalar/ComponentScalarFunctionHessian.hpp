// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/function/scalar/ScalarFunctionHessian.hpp>
#include <sgpp/base/function/vector/VectorFunctionHessian.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
namespace sgpp {
namespace base {

/**
 * One component of a vector-valued function Hessian.
 *
 * @see ComponentScalarFunction
 */
class ComponentScalarFunctionHessian : public ScalarFunctionHessian {
 public:
  /**
   * Constructor.
   *
   * Use it like this:
   * ComponentScalarFunctionHessian gHessian(fHessian, {NAN, NAN, 0.42});
   * where fHessian is a scalar-valued function Hessian with
   * 3 parameters.
   * This selects the first two parameters of fHessian, while
   * constantly using 0.42 for the third parameter.
   *
   * @param fHessian      scalar-valued function Hessian
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
  ComponentScalarFunctionHessian(ScalarFunctionHessian& fHessian,
                                 std::vector<double> defaultValues = std::vector<double>())
      :

        ScalarFunctionHessian((defaultValues.size() > 0)
                                  ? std::count(defaultValues.begin(), defaultValues.end(), NAN)
                                  : fHessian.getNumberOfParameters()),
        fHessianScalar(&fHessian),
        fHessianVector(nullptr),
        dF(fHessian.getNumberOfParameters()),
        k(0),
        defaultValues((defaultValues.size() > 0) ? defaultValues : std::vector<double>(dF, NAN)),
        tmpVec1(dF),
        tmpVec2(dF),
        tmpMat(dF, dF),
        tmpVecMat() {
    initialize();
  }

  /**
   * Constructor.
   *
   * Use it like this:
   * ComponentScalarFunctionHessian gHessian(fHessian, 3, {NAN, NAN, 0.42});
   * where fHessian is a vector-valued function Hessian with
   * 5 components and 3 parameters.
   * This selects the first two parameters and the fourth component
   * of fHessian, while constantly using 0.42 for the third parameter.
   *
   * @param fHessian      vector-valued function Hessian (m components)
   * @param k             index of component \f$f_k\f$ to select
   *                      (between 0 and m - 1)
   * @param defaultValues see other constructor
   */
  ComponentScalarFunctionHessian(VectorFunctionHessian& fHessian, size_t k,
                                 std::vector<double> defaultValues = std::vector<double>())
      :

        ScalarFunctionHessian((defaultValues.size() > 0)
                                  ? std::count(defaultValues.begin(), defaultValues.end(), NAN)
                                  : fHessian.getNumberOfParameters()),
        fHessianScalar(nullptr),
        fHessianVector(&fHessian),
        dF(fHessian.getNumberOfParameters()),
        k(k),
        defaultValues((defaultValues.size() > 0) ? defaultValues : std::vector<double>(dF, NAN)),
        tmpVec1(dF),
        tmpVec2(fHessian.getNumberOfComponents()),
        tmpMat(fHessian.getNumberOfComponents(), dF),
        tmpVecMat(std::vector<base::DataMatrix>(fHessian.getNumberOfComponents(),
                                                base::DataMatrix(dF, dF))) {
    initialize();
  }

  /**
   * Destructor.
   */
  ~ComponentScalarFunctionHessian() override {}

  /**
   * @param[in] x evaluation point \f$\vec{x} \in [0, 1]^n\f$
   * @param[out] gradient \f$\nabla_{\vec{x}} g(\vec{x})\f$
   * @param[out] hessian \f$\nabla_{\vec{x}}^2 g(\vec{x})\f$
   * @return      \f$g(\vec{x}) := f_k(y_1, \dotsc, y_d)\f$
   *              where \f$(x_1, \dotsc, x_n) =
   *              (y_{i_1}, \dotsc, y_{i_n})\f$
   */
  inline double eval(const base::DataVector& x, base::DataVector& gradient,
                     base::DataMatrix& hessian) override {
    size_t t2 = 0, t4;

    // select entries of x which correspond to NAN entries in
    // defaultValues
    for (size_t t = 0; t < dF; t++) {
      if (std::isnan(defaultValues[t])) {
        tmpVec1[t] = x[t2];
        t2++;
      }
    }

    if (fHessianScalar == nullptr) {
      // evaluate and select component
      fHessianVector->eval(tmpVec1, tmpVec2, tmpMat, tmpVecMat);
      t2 = 0;

      for (size_t t = 0; t < dF; t++) {
        if (std::isnan(defaultValues[t])) {
          gradient[t2] = tmpMat(k, t);
          t4 = 0;

          for (size_t t3 = 0; t3 < dF; t3++) {
            if (std::isnan(defaultValues[t3])) {
              hessian(t2, t4) = tmpVecMat[k](t, t3);
              t4++;
            }
          }

          t2++;
        }
      }

      return tmpVec2[k];
    } else {
      // evaluate
      const double fx = fHessianScalar->eval(tmpVec1, tmpVec2, tmpMat);
      t2 = 0;

      for (size_t t = 0; t < dF; t++) {
        if (std::isnan(defaultValues[t])) {
          gradient[t2] = tmpVec2[t];
          t4 = 0;

          for (size_t t3 = 0; t3 < dF; t3++) {
            if (std::isnan(defaultValues[t3])) {
              hessian(t2, t4) = tmpMat(t, t3);
              t4++;
            }
          }

          t2++;
        }
      }

      return fx;
    }
  }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<ScalarFunctionHessian>& clone) const override {
    clone = std::unique_ptr<ScalarFunctionHessian>(new ComponentScalarFunctionHessian(*this));
  }

 protected:
  /// scalar-valued function Hessian
  ScalarFunctionHessian* fHessianScalar;
  /// vector-valued function Hessian
  VectorFunctionHessian* fHessianVector;
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
  /// temporary vector of matrices
  std::vector<base::DataMatrix> tmpVecMat;

  void initialize() {
    // make sure defaultValues has the correct size
    if (defaultValues.size() != dF) {
      throw std::runtime_error(
          "ComponentScalarFunctionHessian::initialize(): Invalid defaultValues.");
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
