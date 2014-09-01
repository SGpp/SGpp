/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_OBJECTIVEGRADIENT_HPP
#define SGPP_OPT_FUNCTION_OBJECTIVEGRADIENT_HPP

#include <vector>
#include <cstddef>

#include "base/datatypes/DataVector.hpp"
#include "opt/tools/SmartPointer.hpp"

namespace sg {
  namespace opt {
    namespace function {

      /**
       * Abstract base class for objective functions \f$f\colon [0, 1]^d \to \mathbb{R}\f$
       * together with their gradients \f$\nabla f\colon [0, 1]^d \to \mathbb{R}^d\f$.
       * They're used in optimization.
       */
      class ObjectiveGradient {
        public:
          /**
           * Constructor.
           *
           * @param d     dimension of the domain
           */
          ObjectiveGradient(size_t d) : d(d) {
          }

          /**
           * Virtual destructor.
           */
          virtual ~ObjectiveGradient() {
          }

          /**
           * Pure virtual method for calculating \f$f(\vec{x})\f$ together with \f$\nabla f(\vec{x})\f$.
           *
           * @param      x            point \f$\vec{x} \in \mathbb{R}^d\f$
           * @param[out] gradient     gradient \f$\nabla f(\vec{x}) \in \mathbb{R}^d\f$
           * @return                  \f$f(\vec{x})\f$
           */
          virtual double evalGradient(const std::vector<double>& x,
                                      base::DataVector& gradient) = 0;

          /**
           * @return dimension \f$d\f$ of the domain
           */
          size_t getDimension() const {
            return d;
          }

          /**
           * Pure virtual method for cloning the objective gradient.
           * It should return a pointer to the cloned object and
           * it's used for parallel computations
           * (the evalGradient() method might not be thread-safe).
           *
           * @return smart pointer to cloned object
           */
          virtual tools::SmartPointer<ObjectiveGradient> clone() = 0;

        protected:
          /// dimension of the domain
          size_t d;
      };

    }
  }
}

#endif
