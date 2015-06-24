// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_CONSTRAINTFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_CONSTRAINTFUNCTION_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <cstddef>

namespace SGPP {
  namespace optimization {

    /**
     * Abstract base class for inequality constraints \f$g(\vec{x}) \le 0\f$
     * defined by \f$g\colon [0, 1]^d \to \mathbb{R}^m\f$.
     */
    class ConstraintFunction {
      public:
        /**
         * Constructor.
         *
         * @param d     dimension of the domain
         * @param m     number of constraints
         */
        ConstraintFunction(size_t d, size_t m) : d(d), m(m) {
        }

        /**
         * Virtual destructor.
         */
        virtual ~ConstraintFunction() {
        }

        /**
         * Pure virtual method for calculating \f$g(\vec{x})\f$.
         *
         * @param[in]  x      evaluation point \f$\vec{x} \in [0, 1]^d\f$
         * @param[out] value  \f$g(\vec{x})\f$
         */
        virtual void eval(const base::DataVector& x,
                          base::DataVector& value) = 0;

        /**
         * @return dimension \f$d\f$ of the domain
         */
        size_t getDimension() const {
          return d;
        }

        /**
         * @return number \f$m\f$ of constraints
         */
        size_t getNumberOfConstraints() const {
          return m;
        }

      protected:
        /// dimension of the domain
        size_t d;
        /// number of constraints
        size_t m;
    };

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_CONSTRAINTFUNCTION_HPP */
