// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G08_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G08_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp>

namespace SGPP {
  namespace optimization {
    namespace test_problems {

      /**
       * G08 objective function.
       *
       * Definition:
       * \f[\bar{f}(\bar{\vec{x}}) :=
       * -\frac{\sin^3(2\pi \bar{x}_1) \sin(2\pi \bar{x}_2)}{\bar{x}_1^3
       * (\bar{x}_1 + \bar{x}_2)}\f]
       */
      class G08Objective : public TestScalarFunction {
        public:
          /**
           * Constructor.
           */
          G08Objective();

          /**
           * Destructor.
           */
          virtual ~G08Objective() override;

          /**
           * @param x     point \f$\vec{x} \in [0, 1]^d\f$
           * @return      \f$f(\vec{x})\f$
           */
          virtual float_t evalUndisplaced(const base::DataVector& x)
          override;

          /**
           * @param[out] clone pointer to cloned object
           */
          virtual void clone(std::unique_ptr<ScalarFunction>& clone)
          const override;
      };

      /**
       * G08 inequality constraint function.
       *
       * Definition:
       * \f[\bar{\vec{g}}(\bar{\vec{x}}) :=
       * \begin{pmatrix}
       *     \bar{x}_1^2 - \bar{x}_2 + 1\\
       *     1 - \bar{x}_1 + (\bar{x}_2 - 4)^2
       * \end{pmatrix}\f]
       */
      class G08InequalityConstraint :
        public TestVectorFunction {
        public:
          /**
           * Constructor.
           */
          G08InequalityConstraint();

          /**
           * Destructor.
           */
          virtual ~G08InequalityConstraint() override;

          /**
           * @param       x       point \f$\vec{x} \in \mathbb{R}^d\f$
           * @param[out]  value   \f$\vec{f}(\vec{x})\f$
           */
          virtual void evalUndisplaced(const base::DataVector& x,
                                       base::DataVector& value) override;

          /**
           * @param[out] clone pointer to cloned object
           */
          virtual void clone(std::unique_ptr<VectorFunction>& clone)
          const override;
      };

      /**
       * G08 equality constraint function.
       *
       * Definition: empty, i.e., no constraint
       */
      class G08EqualityConstraint :
        public TestVectorFunction {
        public:
          /**
           * Constructor.
           */
          G08EqualityConstraint();

          /**
           * Destructor.
           */
          virtual ~G08EqualityConstraint() override;

          /**
           * @param       x       point \f$\vec{x} \in \mathbb{R}^d\f$
           * @param[out]  value   \f$\vec{f}(\vec{x})\f$
           */
          virtual void evalUndisplaced(const base::DataVector& x,
                                       base::DataVector& value) override;

          /**
           * @param[out] clone pointer to cloned object
           */
          virtual void clone(std::unique_ptr<VectorFunction>& clone)
          const override;
      };

      /**
       * G08 constrained test problem.
       *
       * * Number of parameters: 2
       * * Number of inequality constraints: 2
       * * Number of equality constraints: 0
       * * Domain: \f$\bar{\vec{x}} \in [0, 10]^2\f$
       * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
       *   (1.227971, 4.245373)\f$
       * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
       *   -0.09582504\f$
       */
      class G08 : public ConstrainedTestProblem {
        public:
          /**
           * Constructor.
           */
          G08();

          /**
           * Destructor.
           */
          virtual ~G08() override;

          /**
           * @return  objective function of the test problem
           */
          virtual TestScalarFunction& getObjectiveFunction() override;

          /**
           * @return  inequality function of the test problem
           */
          virtual TestVectorFunction& getInequalityConstraintFunction() override;

          /**
           * @return  equality constraint function of the test problem
           */
          virtual TestVectorFunction& getEqualityConstraintFunction() override;

          /**
           * @param[out] x minimal point
           *               \f$\vec{x}_\opt \in [0, 1]^d\f$
           * @return       minimal function value
           *               \f$f(\vec{x}_\opt)\f$
           */
          virtual float_t getOptimalPointUndisplaced(base::DataVector& x)
          override;

        protected:
          /// objective function
          G08Objective f;
          /// inequality constraint function
          G08InequalityConstraint g;
          /// equality constraint function
          G08EqualityConstraint h;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G08_HPP */
