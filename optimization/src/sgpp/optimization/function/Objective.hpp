// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_OBJECTIVE_HPP
#define SGPP_OPTIMIZATION_FUNCTION_OBJECTIVE_HPP

#include <sgpp/globaldef.hpp>

#include <vector>
#include <cstddef>
#include <memory>

namespace SGPP {
  namespace optimization {
    namespace function {

      /**
       * Abstract base class for objective functions \f$f\colon [0, 1]^d \to \mathbb{R}\f$.
       * They're used in adaptive grid generation and optimization.
       */
      class Objective {
        public:
          /**
           * Constructor.
           *
           * @param d     dimension of the domain
           */
          Objective(size_t d) : d(d) {
          }

          /**
           * Virtual destructor.
           */
          virtual ~Objective() {
          }

          /**
           * Pure virtual method for calculating \f$f(\vec{x})\f$.
           *
           * @param x     point \f$\vec{x} \in \mathbb{R}^d\f$
           * @return      \f$f(\vec{x})\f$
           */
          virtual float_t eval(const std::vector<float_t>& x) = 0;

          /**
           * @return dimension \f$d\f$ of the domain
           */
          size_t getDimension() const {
            return d;
          }

          /**
           * Pure virtual method for cloning the objective function.
           * It should return a pointer to the cloned object and it's used for parallel computations
           * (the eval() method might not be thread-safe).
           *
           * @return pointer to cloned object
           */
          virtual Objective* clone() const = 0;

        protected:
          /// dimension of the domain
          size_t d;
      };

    }
  }
}

#endif
