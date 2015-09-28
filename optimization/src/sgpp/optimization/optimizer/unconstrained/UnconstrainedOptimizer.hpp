// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_UNCONSTRAINEDOPTIMIZER_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_UNCONSTRAINEDOPTIMIZER_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>

#include <cstddef>
#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      /**
       * Abstract class for optimizing objective functions.
       */
      class UnconstrainedOptimizer {
        public:
          /// default maximal number of iterations or function evaluations
          static const size_t DEFAULT_N = 1000;

          /**
           * Constructor.
           * The starting point is set to
           * \f$(0.5, \dotsc, 0.5)^{\mathrm{T}}\f$.
           *
           * @param f     function to optimize
           * @param N     maximal number of iterations or function evaluations
           *              (depending on the implementation)
           */
          UnconstrainedOptimizer(ScalarFunction& f, size_t N = DEFAULT_N) :
            f(f),
            N(N),
            x0(f.getDimension(), 0.5),
            xOpt(0),
            fOpt(NAN),
            xHist(0, 0),
            fHist(0) {
          }

          /**
           * Virtual destructor.
           */
          virtual ~UnconstrainedOptimizer() {
          }

          /**
           * Pure virtual method for optimization of the objective function.
           * The result of the optimization process can be obtained by
           * member functions, e.g., getOptimalPoint() and getOptimalValue().
           */
          virtual void optimize() = 0;

          /**
           * @return objective function
           */
          ScalarFunction& getObjectiveFunction() const {
            return f;
          }

          /**
           * @return  maximal number of iterations or function evaluations
           */
          size_t getN() const {
            return N;
          }

          /**
           * @param N maximal number of iterations or function evaluations
           */
          void setN(size_t N) {
            this->N = N;
          }

          /**
           * @return                  starting point
           */
          const base::DataVector& getStartingPoint() const {
            return x0;
          }

          /**
           * @param startingPoint     starting point
           */
          void setStartingPoint(const base::DataVector& startingPoint) {
            this->x0 = startingPoint;
          }

          /**
           * @return result of optimization (location of optimum),
           *         empty vector on error
           */
          const base::DataVector& getOptimalPoint() const {
            return xOpt;
          }

          /**
           * @return result of optimization (optimal function value),
           *         NAN on error
           */
          float_t getOptimalValue() const {
            return fOpt;
          }

          /**
           * @return tall matrix (d columns) in which the k-th row indicates
           *         the best point after k iterations of the optimization,
           *         empty matrix on error or if not supported
           */
          const base::DataMatrix& getHistoryOfOptimalPoints() const {
            return xHist;
          }

          /**
           * @return vector in which the k-th entry indicates the best
           *         function value after k iterations of the optimization,
           *         empty vector on error or if not supported
           */
          const base::DataVector& getHistoryOfOptimalValues() const {
            return fHist;
          }

        protected:
          /// objective function
          ScalarFunction& f;
          /// maximal number of iterations or function evaluations
          size_t N;
          /// starting point
          base::DataVector x0;
          /// result of optimization (location of optimum)
          base::DataVector xOpt;
          /// result of optimization (optimal function value)
          float_t fOpt;
          /// search history matrix (optimal points)
          base::DataMatrix xHist;
          /// search history vector (optimal values)
          base::DataVector fHist;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_UNCONSTRAINEDOPTIMIZER_HPP */
