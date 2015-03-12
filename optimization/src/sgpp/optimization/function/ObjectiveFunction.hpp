// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_OBJECTIVEFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_OBJECTIVEFUNCTION_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <cstddef>
#include <memory>

namespace SGPP {
  namespace optimization {

    /**
     * Abstract base class for objective functions
     * \f$f\colon [0, 1]^d \to \mathbb{R}\f$.
     * They're used in adaptive grid generation and optimization.
     */
    class ObjectiveFunction {
      public:
        /**
         * Constructor.
         *
         * @param d     dimension of the domain
         */
        ObjectiveFunction(size_t d) : d(d) {
        }

        /**
         * Virtual destructor.
         */
        virtual ~ObjectiveFunction() {
        }

        /**
         * Pure virtual method for calculating \f$f(\vec{x})\f$.
         *
         * @param x     point \f$\vec{x} \in \mathbb{R}^d\f$
         * @return      \f$f(\vec{x})\f$
         */
        virtual float_t eval(const base::DataVector& x) = 0;

        /**
         * @return dimension \f$d\f$ of the domain
         */
        size_t getDimension() const {
          return d;
        }

        /**
         * Pure virtual method for cloning the objective function.
         * It should generate a pointer to the cloned object and
         * it's used for parallel computations
         * (the eval() method might not be thread-safe).
         *
         * @param[out] clone pointer to cloned object
         */
        virtual void clone(
          std::unique_ptr<ObjectiveFunction>& clone) const = 0;

      protected:
        /// dimension of the domain
        size_t d;
    };

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_OBJECTIVEFUNCTION_HPP */
