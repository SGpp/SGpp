// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_CONSTRAINEDOPTIMIZER_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_CONSTRAINEDOPTIMIZER_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

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
  ConstrainedOptimizer(const ScalarFunction& f, const VectorFunction& g, const VectorFunction& h,
                       size_t N = DEFAULT_N)
      : UnconstrainedOptimizer(f, N) {
    g.clone(this->g);
    h.clone(this->h);
  }

  /**
   * Copy constructor.
   *
   * @param other optimizer to be copied
   */
  ConstrainedOptimizer(const ConstrainedOptimizer& other)
      : UnconstrainedOptimizer(other) {
    other.g->clone(g);
    other.h->clone(h);
  }

  /**
   * Destructor.
   */
  ~ConstrainedOptimizer() override {}

  /**
   * @return inequality constraint function
   */
  VectorFunction& getInequalityConstraintFunction() const { return *g; }

  /**
   * @param g  inequality constraint function
   */
  void setInequalityConstraintFunction(const VectorFunction& g) {
    g.clone(this->g);
  }

  /**
   * @return equality constraint function
   */
  VectorFunction& getEqualityConstraintFunction() const { return *h; }

  /**
   * @param h  equality constraint function
   */
  void setEqualityConstraintFunction(const VectorFunction& h) { h.clone(this->h); }

 protected:
  /// inequality constraint function
  std::unique_ptr<VectorFunction> g;
  /// equality constraint function
  std::unique_ptr<VectorFunction> h;
};
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_CONSTRAINEDOPTIMIZER_HPP */
