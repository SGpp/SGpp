// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_HARTMAN6_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_HARTMAN6_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Hartman6 objective function.
 *
 * Definition:
 * \f{gather*}{
 * \bar{f}(\bar{\vec{x}}) := -\sum_{i=1}^4 a_i
 * \exp\!\left(-\sum_{t=1}^6 b_{i,t} (\bar{x}_t - c_{i,t})^2\right),\\
 * \vec{a} = \begin{pmatrix}1\\1.2\\3\\3.2\end{pmatrix},\qquad
 * B :=
 *      \begin{pmatrix}
 *          10 & 3 & 17 & 3.5 & 1.7 & 8\\
 *          0.05 & 10 & 17 & 0.1 & 8 & 14\\
 *          3 & 3.5 & 1.7 & 10 & 17 & 8\\
 *          17 & 8 & 0.05 & 10 & 0.1 & 14
 *      \end{pmatrix},\\
 * C :=
 *      \begin{pmatrix}
 *          0.1312 & 0.1696 & 0.5569 & 0.0124 & 0.8283 & 0.5886\\
 *          0.2329 & 0.4135 & 0.8307 & 0.3736 & 0.1004 & 0.9991\\
 *          0.2348 & 0.1451 & 0.3522 & 0.2883 & 0.3047 & 0.6650\\
 *          0.4047 & 0.8828 & 0.8732 & 0.5743 & 0.1091 & 0.0381
 *      \end{pmatrix}\f}
 */
class Hartman6Objective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  Hartman6Objective();

  /**
   * Destructor.
   */
  ~Hartman6Objective() override;

  /**
   * @param x     point \f$\vec{x} \in [0, 1]^d\f$
   * @return      \f$f(\vec{x})\f$
   */
  double evalUndisplaced(const base::DataVector& x) override;

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<base::ScalarFunction>& clone) const override;
};

/**
 * Hartman6 unconstrained test problem.
 *
 * * Number of parameters: 6
 * * Domain: \f$\bar{\vec{x}} \in \f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (0.2016895, 0.1500107, 0.4768740,
 *   0.2753324, 0.3116516, 0.6573005)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   -3.322368\f$
 */
class Hartman6 : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  Hartman6();

  /**
   * Destructor.
   */
  ~Hartman6() override;

  /**
   * @return  objective function of the test problem
   */
  TestScalarFunction& getObjectiveFunction() override;

  /**
   * @param[out] x minimal point
   *               \f$\vec{x}_\opt \in [0, 1]^d\f$
   * @return       minimal function value
   *               \f$f(\vec{x}_\opt)\f$
   */
  double getOptimalPointUndisplaced(base::DataVector& x) override;

 protected:
  /// objective function
  Hartman6Objective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_HARTMAN6_HPP */
