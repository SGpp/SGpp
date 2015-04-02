// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_INTERPOLANTFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_INTERPOLANTFUNCTION_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/ObjectiveFunction.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEval.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <cstring>
#include <memory>

namespace SGPP {
  namespace optimization {

    /**
     * Sparse grid interpolant as an objective function.
     *
     * More generally, the function can be any linear combination
     * \f$f\colon [0, 1]^d \to \mathbb{R}\f$,
     * \f$f(\vec{x}) = \sum_{k=1}^N \alpha_k \varphi_k(\vec{x})\f$
     * of the basis functions
     * \f$\varphi_k = \varphi_{\vec{\ell}_k,\vec{i}_k}\f$
     * of a sparse grid with grid points
     * \f$\vec{x}_k = \vec{x}_{\vec{\ell}_k,\vec{i}_k}\f$.
     * But most often, the function (e.g. its coefficients) is constructed
     * as an interpolant at the grid points for some function values.
     */
    class InterpolantFunction : public ObjectiveFunction {
      public:
        /**
         * Constructor.
         * Do not destruct the grid before the InterpolantFunction object!
         *
         * @param grid  sparse grid
         * @param alpha coefficient vector
         */
        InterpolantFunction(base::Grid& grid, const base::DataVector& alpha) :
          ObjectiveFunction(grid.getStorage()->dim()),
          grid(grid),
          opEval(op_factory::createOperationNaiveEval(grid)),
          alpha(alpha) {
        }

        /**
         * Evaluation of the function.
         *
         * @param x     point \f$\vec{x} \in \mathbb{R}^d\f$
         * @return      \f$f(\vec{x})\f$
         */
        inline float_t eval(const base::DataVector& x) {
          // copy x, necessary due to non-existing const correctness
          // in SGPP::base
          base::DataVector y(x);
          return opEval->eval(alpha, y);
        }

        /**
         * @param[out] clone pointer to cloned object
         */
        virtual void clone(std::unique_ptr<ObjectiveFunction>& clone) const {
          clone = std::unique_ptr<ObjectiveFunction>(
                    new InterpolantFunction(grid, alpha));
        }

      protected:
        /// sparse grid
        base::Grid& grid;
        /// pointer to evaluation operation
        std::unique_ptr<base::OperationNaiveEval> opEval;
        /// coefficient vector
        base::DataVector alpha;
    };

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_INTERPOLANTFUNCTION_HPP */
