// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_INTERPOLANTHESSIAN_HPP
#define SGPP_OPTIMIZATION_FUNCTION_INTERPOLANTHESSIAN_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/ObjectiveHessian.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalHessian.hpp>

namespace SGPP {
  namespace optimization {


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
         * @param grid  sparse grid
         * @param alpha coefficient vector
         */
        InterpolantHessian(base::Grid& grid, const base::DataVector& alpha) :
          ObjectiveHessian(grid.getStorage()->dim()),
          grid(grid),
          opEvalHessian(op_factory::createOperationNaiveEvalHessian(grid)),
          alpha(alpha) {
        }

        /**
         * Evaluation of the function, its gradient and its Hessian.
         *
         * @param      x        point \f$\vec{x} \in \mathbb{R}^d\f$
         * @param[out] gradient gradient
         *                      \f$\nabla f(\vec{x}) \in \mathbb{R}^d\f$
         * @param[out] hessian  Hessian matrix
         *                      \f$H_f(\vec{x}) \in \mathbb{R}^{d \times d}\f$
         * @return              \f$f(\vec{x})\f$
         */
        inline float_t evalHessian(const base::DataVector& x,
                                   base::DataVector& gradient,
                                   base::DataMatrix& hessian) {
          // copy x, necessary due to non-existing const correctness
          // in SGPP::base
          base::DataVector y(x);
          return opEvalHessian->evalHessian(alpha, y, gradient, hessian);
        }

        /**
         * @param[out] clone pointer to cloned object
         */
        virtual void clone(std::unique_ptr<ObjectiveHessian>& clone) const {
          clone = std::unique_ptr<ObjectiveHessian>(
                    new InterpolantHessian(grid, alpha));
        }

      protected:
        base::Grid& grid;
        /// pointer to evaluation operation
        std::unique_ptr<base::OperationNaiveEvalHessian> opEvalHessian;
        /// coefficient vector
        base::DataVector alpha;
    };

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_INTERPOLANTHESSIAN_HPP */
