// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_LINESEARCHARMIJO_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_LINESEARCHARMIJO_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>

#include <cmath>
#include <cstddef>
namespace sgpp {
namespace optimization {
namespace optimizer {

/**
 * Line search (1D optimization on a line) with Armijo's rule used in
 * gradient-based optimization.
 *
 * Armijo's rule calculates \f$\sigma = \beta^k\f$ for
 * \f$k = 0, 1, \dotsc\f$ for a fixed \f$\beta \in (0, 1)\f$ and checks
 * if \f$\vec{y} = \vec{x} + \sigma\vec{s}\f$
 * lies in \f$[0, 1]^d\f$ and whether
 * the objective function value improvement meets the condition
 * \f$f(\vec{x}) - f(\vec{y}) \ge
 * \gamma\sigma (-\nabla f(\vec{x}) \cdot \vec{s})\f$
 * for \f$\gamma \in (0, 1)\f$ fixed.
 *
 * The return value states whether the relative improvement
 * (depending on two tolerances)
 * is big enough to continue the optimization algorithm.
 *
 * @param       f             objective function
 * @param       beta          \f$\beta \in (0, 1)\f$
 * @param       gamma         \f$\gamma \in (0, 1)\f$
 * @param       tol           tolerance 1 (positive)
 * @param       eps           tolerance 2 (positive)
 * @param       x             point to start the line search in
 * @param       fx            objective function value in x
 * @param       gradFx        objective function gradient in x
 * @param       s             search direction (should be normalized)
 * @param[out]  y             new point, must have the same size as x
 *                            before calling this function
 * @param[in,out] evalCounter this variable will be increased by the
 *                            number of evaluations of f while searching
 * @return                    whether the new point reaches an acceptable
 *                            improvement
 */
inline bool lineSearchArmijo(base::ScalarFunction& f, double beta, double gamma, double tol,
                             double eps, const base::DataVector& x, double fx,
                             base::DataVector& gradFx, const base::DataVector& s,
                             base::DataVector& y, size_t& evalCounter) {
  const size_t d = x.getSize();
  double sigma = 1.0;
  const double ip = gradFx.dotProduct(s);
  double fy = fx;

  for (size_t k = 0; k < 100; k++) {
    bool y_in_domain = true;

    // calculate new point
    for (size_t t = 0; t < d; t++) {
      y[t] = x[t] + sigma * s[t];

      if ((y[t] < 0.0) || (y[t] > 1.0)) {
        y_in_domain = false;
        break;
      }
    }

    // check if y lies in [0, 1]^d
    if (y_in_domain) {
      fy = f.eval(y);
      evalCounter++;

      double improvement = fx - fy;
      double rhs = gamma * sigma * (-ip);

      // check if the absolute improvement is big enough
      if (improvement >= rhs) {
        return (std::abs(fx - fy) >= tol * (std::abs(fx) + std::abs(fy) + eps));
      }
    }

    // next sigma
    sigma *= beta;
  }

  return false;
}
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_LINESEARCHARMIJO_HPP */
