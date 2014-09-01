/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_INTERPOLANTHESSIAN_HPP
#define SGPP_OPT_FUNCTION_INTERPOLANTHESSIAN_HPP

#include "opt/function/ObjectiveHessian.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/operation/OperationNaiveEvalHessian.hpp"

#include <vector>

namespace sg {
  namespace opt {
    namespace function {


      /**
       * Sparse grid interpolant Hessian as an objective Hessian.
       *
       * @see Interpolant
       */
      class InterpolantHessian : public ObjectiveHessian {
        public:
          /**
           * Constructor.
           * Do not destruct the grid before the InterpolantHessian object!
           *
           * @param d     dimension of the domain
           * @param grid  sparse grid
           * @param alpha coefficient vector
           */
          InterpolantHessian(size_t d, base::Grid& grid, base::DataVector& alpha) :
            ObjectiveHessian(d), grid(grid),
            op_eval_hessian(op_factory::createOperationNaiveEvalHessian(grid)),
            alpha(alpha) {
          }

          /**
           * Evaluation of the function, its gradient and its Hessian.
           *
           * @param      x            point \f$\vec{x} \in \mathbb{R}^d\f$
           * @param[out] gradient     gradient \f$\nabla f(\vec{x}) \in \mathbb{R}^d\f$
           * @param[out] hessian      Hessian matrix \f$H_f(\vec{x}) \in \mathbb{R}^{d \times d}\f$
           * @return                  \f$f(\vec{x})\f$
           */
          inline double evalHessian(const std::vector<double>& x,
                                    base::DataVector& gradient, base::DataMatrix& hessian) {
            // copy x, necessary due to non-existing const correctness in sg::base
            std::vector<double> y = x;
            return op_eval_hessian->evalHessian(alpha, y, gradient, hessian);
          }

          /**
           * @return clone of the object
           */
          virtual tools::SmartPointer<ObjectiveHessian> clone() {
            return tools::SmartPointer<ObjectiveHessian>(new InterpolantHessian(d, grid, alpha));
          }

        protected:
          base::Grid& grid;
          /// pointer to evaluation operation
          tools::SmartPointer<base::OperationNaiveEvalHessian> op_eval_hessian;
          /// coefficient vector
          base::DataVector alpha;
      };

    }
  }
}

#endif
