// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_HARTMAN3_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_HARTMAN3_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Hartman3 objective function.
 *
 * Definition:
 * \f{gather*}{
 * \bar{f}(\bar{\vec{x}}) := -\sum_{i=1}^4 a_i
 * \exp\!\left(-\sum_{t=1}^3 b_{i,t} (\bar{x}_t - c_{i,t})^2\right),\\
 * \vec{a} = \begin{pmatrix}1\\1.2\\3\\3.2\end{pmatrix},\qquad
 * B :=
 *      \begin{pmatrix}
 *          3 & 10 & 30\\
 *          0.1 & 10 & 35\\
 *          3 & 10 & 30\\
 *          0.1 & 10 & 35
 *      \end{pmatrix},\qquad
 * C :=
 *      \begin{pmatrix}
 *          0.3689 & 0.1170 & 0.2673\\
 *          0.4699 & 0.4387 & 0.7470\\
 *          0.1091 & 0.8732 & 0.5547\\
 *          0.0382 & 0.5743 & 0.8828
 *      \end{pmatrix}\f}
 */
class Hartman3Objective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  Hartman3Objective();

  /**
   * Destructor.
   */
  ~Hartman3Objective() override;

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
 * Hartman3 unconstrained test problem.
 *
 * * Number of parameters: 3
 * * Domain: \f$\bar{\vec{x}} \in [0, 1]^3\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (0.1146398, 0.5556488, 0.8525470)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   -3.862785\f$
 */
class Hartman3 : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  Hartman3();

  /**
   * Destructor.
   */
  ~Hartman3() override;

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
  Hartman3Objective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_HARTMAN3_HPP */
