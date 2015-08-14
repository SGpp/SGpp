// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_TESTFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_TESTFUNCTION_HPP

#include <cstddef>

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/ObjectiveFunction.hpp>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      /**
       * Base class for analytical objective function examples
       * ("test functions").
       * The only difference to ObjectiveFunction is the possibility for
       * pseudorandom displacements of the function and the specification of
       * the (or an) minimal point.
       *
       * Taking the average of results of multiple runs with different
       * displacements makes results more robust and significant.
       * The displaced function is \f$\vec{x} \mapsto f(\vec{x} + \vec{d})\f$
       * for a vector \f$\vec{d}\f$ ("displacement") with each component
       * \f$d_t\f$ distributed normally with mean 0 and a specific standard
       * deviation.
       */
      class TestFunction : public ObjectiveFunction {
        public:
          /// default standard deviation
          static constexpr float_t DEFAULT_STANDARD_DEVIATION = 0.01;

          /**
           * Constructor.
           * The displacement is set to all zeros, so to displace the function
           * call generateDisplacement() afterwards.
           *
           * @param d     dimension of the domain
           */
          TestFunction(size_t d);

          /**
           * Virtual destructor.
           */
          virtual ~TestFunction();

          /**
           * Evaluate displaced function.
           *
           * @param x point \f$\vec{x} \in \mathbb{R}^d\f$
           * @return  \f$f(\vec{x} + \vec{d})\f$
           *          with displacement \f$\vec{d}\f$
           */
          float_t eval(const base::DataVector& x);

          /**
           * Pure virtual method for evaluating the undisplaced function.
           *
           * @param x     point \f$\vec{x} \in \mathbb{R}^d\f$
           * @return      \f$f(\vec{x})\f$
           */
          virtual float_t evalUndisplaced(const base::DataVector& x) = 0;

          /**
           * Returns the minimal point of the displaced function.
           *
           * @param[out] x reverse displaced minimal point
           *               \f$\vec{x}_{\text{opt}} - \vec{d}\f$
           * @return       minimal function value
           *               \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
           */
          float_t getOptimalPoint(base::DataVector& x);

          /**
           * Pure virtual method returning the minimal point
           * (or one of the minimal points, if there are multiple of them)
           * of the test function.
           *
           * @param[out] x    minimal point \f$\vec{x}_{\text{opt}}\f$
           * @return          minimal function value
           *                  \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
           */
          virtual float_t getOptimalPointUndisplaced(base::DataVector& x) = 0;

          /**
           * Generate normally distributed pseudorandom displacement with
           * default standard deviation.
           */
          virtual void generateDisplacement();

          /**
           * Generate normally distributed pseudorandom displacement.
           * This function can be overridden if the minimal points of the
           * test function lie near or on the boundary or if there is the
           * chance of a pole getting in the domain by displacing.
           *
           * @param stdDev standard deviation of the displacement coordinates
           */
          virtual void generateDisplacement(float_t stdDev);

          /**
           * Add the displacement to a vector.
           *
           * @param[in,out] x     vector to be displaced
           */
          void displaceVector(base::DataVector& x) const;

          /**
           * Subtract the displacement from a vector.
           *
           * @param[in,out] x     vector to be reverse displaced
           */
          void reverseDisplaceVector(base::DataVector& x) const;

          /**
           * @return standard deviation of the displacement
           */
          float_t getStandardDeviation() const;

          /**
           * @param[out] displacement     currently used displacement
           */
          void getDisplacement(base::DataVector& displacement) const;

        protected:
          /// standard deviation
          float_t stdDev;
          /// vector displacement
          base::DataVector displacement;
          /// temporary vector for displacing
          base::DataVector xTmp;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_TEST_TESTFUNCTION_HPP */
