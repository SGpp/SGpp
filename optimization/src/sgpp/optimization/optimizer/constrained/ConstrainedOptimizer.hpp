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
  ConstrainedOptimizer(const base::ScalarFunction& f, const base::VectorFunction& g,
                       const base::VectorFunction& h, size_t N = DEFAULT_N)
      : UnconstrainedOptimizer(f, N) {
    g.clone(this->g);
    h.clone(this->h);
  }

  /**
   * Copy constructor.
   *
   * @param other optimizer to be copied
   */
  ConstrainedOptimizer(const ConstrainedOptimizer& other) : UnconstrainedOptimizer(other) {
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
  base::VectorFunction& getInequalityConstraintFunction() const { return *g; }

  /**
   * @return equality constraint function
   */
  base::VectorFunction& getEqualityConstraintFunction() const { return *h; }

 protected:
  /// inequality constraint function
  std::unique_ptr<base::VectorFunction> g;
  /// equality constraint function
  std::unique_ptr<base::VectorFunction> h;
};
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp
