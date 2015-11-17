// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_SCALAR_INTERPOLANTSCALARFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_SCALAR_INTERPOLANTSCALARFUNCTION_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>

#include <cstring>
#include <memory>

namespace SGPP {
  namespace optimization {

    /**
     * Sparse grid interpolant of a scalar-valued function.
     *
     * More generally, the function can be any linear combination
     * \f$f\colon [0, 1]^d \to \mathbb{R}\f$,
     * \f$f(\vec{x}) = \sum_{k=1}^N \alpha_k \varphi_k(\vec{x})\f$
     * of the basis functions
     * \f$\varphi_k = \varphi_{\vec{\ell}_k,\vec{i}_k}\f$
     * of a sparse grid with grid points
     * \f$\vec{x}_k = \vec{x}_{\vec{\ell}_k,\vec{i}_k}\f$.
     * But most often, the function (e.g., its coefficients) is constructed
     * as an interpolant at the grid points for some function values.
     */
    class InterpolantScalarFunction : public ScalarFunction {
      public:
        /**
         * Constructor.
         * Do not destruct the grid before the
         * InterpolantScalarFunction object!
         *
         * @param grid  sparse grid
         * @param alpha coefficient vector
         */
        InterpolantScalarFunction(base::Grid& grid,
                                  const base::DataVector& alpha) :
          ScalarFunction(grid.getStorage()->dim()),
          grid(grid),
          opEval(op_factory::createOperationNaiveEval(grid)),
          alpha(alpha) {
        }

        /**
         * Destructor.
         */
        virtual ~InterpolantScalarFunction() override {
        }

        /**
         * Evaluation of the function.
         *
         * @param x     evaluation point \f$\vec{x} \in [0, 1]^d\f$
         * @return      \f$f(\vec{x})\f$
         */
        inline virtual float_t eval(const base::DataVector& x) override {
          for (size_t t = 0; t < d; t++) {
            if ((x[t] < 0.0) || (x[t] > 1.0)) {
              return INFINITY;
            }
          }

          return opEval->eval(alpha, x);
        }

        /**
         * @param[out] clone pointer to cloned object
         */
        virtual void clone(std::unique_ptr<ScalarFunction>& clone) const override {
          clone = std::unique_ptr<ScalarFunction>(
                    new InterpolantScalarFunction(grid, alpha));
        }

        /**
         * @return coefficient vector
         */
        const base::DataVector& getAlpha() const {
          return alpha;
        }

        /**
         * @param alpha coefficient vector
         */
        void setAlpha(const base::DataVector& alpha) {
          this->alpha = alpha;
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

#endif /* SGPP_OPTIMIZATION_FUNCTION_SCALAR_INTERPOLANTSCALARFUNCTION_HPP */
