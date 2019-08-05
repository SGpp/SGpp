// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>
#include <sgpp/optimization/optimizer/unconstrained/NelderMead.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp>

#include <cstddef>
#include <memory>

namespace sgpp {
namespace optimization {
namespace optimizer {

/**
 * Abstract class for solving constrained optimization problems.
 */
class ConstrainedOptimizer : public UnconstrainedOptimizer {
 public:
  /**
   * Constructor.
   * By default, Nelder-Mead is used as optimization algorithm
   * (gradient-free).
   * The starting point is set to
   * \f$(0.5, \dotsc, 0.5)^{\mathrm{T}}\f$.
   * Depending on the implementation $g$ and/or $h$ may be ignored
   * (if only equality or inequality constraints can be handled
   * by the underlying algorithm).
   *
   * @param f     function to optimize
   * @param g     inequality constraint function
   *              (\f$g(\vec{x}) \le 0\f$)
   * @param h     equality constraint function
   *              (\f$h(\vec{x}) = 0\f$)
   * @param N     maximal number of iterations or
   *              objective function evaluations
   *              (depending on the implementation)
   */
  ConstrainedOptimizer(const base::ScalarFunction& f, const base::VectorFunction& g,
                       const base::VectorFunction& h, size_t N = DEFAULT_N)
      : UnconstrainedOptimizer(f, nullptr, nullptr, N),
        unconstrainedOptimizer(new NelderMead(f)) {
    g.clone(this->g);
    h.clone(this->h);
  }

  /**
   * Constructor.
   * By default, adaptive gradient descent is used as optimization algorithm
   * (gradient-based).
   * The starting point is set to
   * \f$(0.5, \dotsc, 0.5)^{\mathrm{T}}\f$.
   * Depending on the implementation $g$ and/or $h$ may be ignored
   * (if only equality or inequality constraints can be handled
   * by the underlying algorithm).
   *
   * @param f           function to optimize
   * @param fGradient   gradient of f
   * @param g           inequality constraint function
   *                    (\f$g(\vec{x}) \le 0\f$)
   * @param gGradient   gradient of g
   * @param h           equality constraint function
   *                    (\f$h(\vec{x}) = 0\f$)
   * @param hGradient   gradient of h
   * @param N           maximal number of iterations or
   *                    objective function evaluations
   *                    (depending on the implementation)
   */
  ConstrainedOptimizer(const base::ScalarFunction& f,
                       const base::ScalarFunctionGradient& fGradient,
                       const base::VectorFunction& g,
                       const base::VectorFunctionGradient& gGradient,
                       const base::VectorFunction& h,
                       const base::VectorFunctionGradient& hGradient,
                       size_t N = DEFAULT_N)
      : UnconstrainedOptimizer(f, &fGradient, nullptr, N),
        unconstrainedOptimizer(new AdaptiveGradientDescent(f, fGradient)) {
    g.clone(this->g);
    gGradient.clone(this->gGradient);
    h.clone(this->h);
    hGradient.clone(this->hGradient);
  }

  /**
   * Constructor with custom unconstrained optimization algorithm
   * (gradient-free or gradient-based).
   * The starting point is set to
   * \f$(0.5, \dotsc, 0.5)^{\mathrm{T}}\f$.
   * Depending on the implementation $g$ and/or $h$ may be ignored
   * (if only equality or inequality constraints can be handled
   * by the underlying algorithm).
   *
   * @param unconstrainedOptimizer  unconstrained optimizer
   * @param g                       inequality constraint function
   *                                (\f$g(\vec{x}) \le 0\f$)
   * @param gGradient               gradient of g (nullptr to omit)
   * @param h                       equality constraint function
   *                                (\f$h(\vec{x}) = 0\f$)
   * @param hGradient               gradient of h (nullptr to omit)
   * @param N                       maximal number of iterations or
   *                                objective function evaluations
   *                                (depending on the implementation)
   */
  ConstrainedOptimizer(const UnconstrainedOptimizer& unconstrainedOptimizer,
                       const base::VectorFunction& g,
                       const base::VectorFunctionGradient* gGradient,
                       const base::VectorFunction& h,
                       const base::VectorFunctionGradient* hGradient,
                       size_t N = DEFAULT_N)
      : UnconstrainedOptimizer(unconstrainedOptimizer.getObjectiveFunction(),
                               unconstrainedOptimizer.getObjectiveGradient(),
                               unconstrainedOptimizer.getObjectiveHessian(),
                               N) {
    unconstrainedOptimizer.clone(this->unconstrainedOptimizer);
  }

  /**
   * Copy constructor.
   *
   * @param other optimizer to be copied
   */
  ConstrainedOptimizer(const ConstrainedOptimizer& other)
      : ConstrainedOptimizer(*unconstrainedOptimizer,
                             *other.g, other.gGradient.get(),
                             *other.h, other.hGradient.get(), N) {
  }

  /**
   * Destructor.
   */
  ~ConstrainedOptimizer() override {}

  /**
   * @return inequality constraint function
   */
  base::VectorFunction& getInequalityConstraintFunction() const { return *g; }

  /**
   * @param g  inequality constraint function
   */
  void setInequalityConstraintFunction(const base::VectorFunction& g) {
    g.clone(this->g);
  }

  /**
   * @return inequality constraint function gradient
   */
  base::VectorFunctionGradient* getInequalityConstraintGradient() const { return gGradient.get(); }

  /**
   * @param gGradient  inequality constraint function gradient
   */
  void setInequalityConstraintGradient(const base::VectorFunctionGradient* gGradient) {
    if (gGradient != nullptr) {
      gGradient->clone(this->gGradient);
    } else {
      this->gGradient = nullptr;
    }
  }

  /**
   * @return equality constraint function
   */
  base::VectorFunction& getEqualityConstraintFunction() const { return *h; }

  /**
   * @param h  equality constraint function
   */
  void setEqualityConstraintFunction(const base::VectorFunction& h) { h.clone(this->h); }

  /**
   * @return equality constraint function gradient
   */
  base::VectorFunctionGradient* getEqualityConstraintGradient() const { return hGradient.get(); }

  /**
   * @param hGradient  equality constraint function gradient
   */
  void setEqualityConstraintGradient(const base::VectorFunctionGradient* hGradient) {
    if (hGradient != nullptr) {
      hGradient->clone(this->hGradient);
    } else {
      this->hGradient = nullptr;
    }
  }

 protected:
  /// unconstrained optimization algorithm
  std::unique_ptr<UnconstrainedOptimizer> unconstrainedOptimizer;
  /// inequality constraint function
  std::unique_ptr<base::VectorFunction> g;
  /// inequality constraint function gradient
  std::unique_ptr<base::VectorFunctionGradient> gGradient;
  /// equality constraint function
  std::unique_ptr<base::VectorFunction> h;
  /// equality constraint function gradient
  std::unique_ptr<base::VectorFunctionGradient> hGradient;
};
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp
