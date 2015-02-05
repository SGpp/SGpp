// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_OBJECTIVEHESSIAN_HPP
#define SGPP_OPTIMIZATION_FUNCTION_OBJECTIVEHESSIAN_HPP

#include <sgpp/globaldef.hpp>

#include <vector>
#include <cstddef>
#include <memory>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

namespace SGPP {
  namespace optimization {
    namespace function {

      /**
       * Abstract base class for objective functions \f$f\colon [0, 1]^d \to \mathbb{R}\f$
       * together with their gradients \f$\nabla f\colon [0, 1]^d \to \mathbb{R}^d\f$
       * and Hessians \f$H_f\colon [0, 1]^d \to \mathbb{R}^{d \times d}\f$.
       * They're used in optimization.
       */
      class ObjectiveHessian {
        public:
          /**
           * Constructor.
           *
           * @param d     dimension of the domain
           */
          ObjectiveHessian(size_t d) : d(d) {
          }

          /**
           * Virtual destructor.
           */
          virtual ~ObjectiveHessian() {
          }

          /**
           * Pure virtual method for calculating \f$f(\vec{x})\f$ together with \f$\nabla f(\vec{x})\f$
           * and \f$H_f(\vec{x}) \in \mathbb{R}^{d \times d}\f$.
           *
           * @param      x            point \f$\vec{x} \in \mathbb{R}^d\f$
           * @param[out] gradient     gradient \f$\nabla f(\vec{x}) \in \mathbb{R}^d\f$
           * @param[out] hessian      Hessian matrix \f$H_f(\vec{x}) \in \mathbb{R}^{d \times d}\f$
           * @return                  \f$f(\vec{x})\f$
           */
          virtual float_t evalHessian(const std::vector<float_t>& x,
                                     base::DataVector& gradient, base::DataMatrix& hessian) = 0;

          /**
           * @return dimension \f$d\f$ of the domain
           */
          size_t getDimension() const {
            return d;
          }

          /**
           * Pure virtual method for cloning the objective Hessian.
           * It should return a pointer to the cloned object and
           * it's used for parallel computations
           * (the evalHessian() method might not be thread-safe).
           *
           * @return pointer to cloned object
           */
          virtual ObjectiveHessian* clone() const = 0;

        protected:
          /// dimension of the domain
          size_t d;
      };

    }
  }
}

#endif
