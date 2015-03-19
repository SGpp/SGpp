// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_OPTIMIZER_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_OPTIMIZER_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/ObjectiveFunction.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <cstddef>
#include <memory>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      /**
       * Abstract class for optimizing objective functions.
       */
      class Optimizer {
        public:
          /// default maximal number of iterations or function evaluations
          static const size_t DEFAULT_N = 200;

          /**
           * Constructor.
           * The starting point is set to
           * \f$(0.5, \dotsc, 0.5)^{\mathrm{T}}\f$.
           *
           * @param f     function to optimize
           * @param N     maximal number of iterations or function evaluations
           *              (depending on the implementation)
           */
          Optimizer(ObjectiveFunction& f, size_t N = DEFAULT_N) :
            f(f), N(N), x0(f.getDimension(), 0.5) {
          }

          /**
           * Virtual destructor.
           */
          virtual ~Optimizer() {
          }

          /**
           * Pure virtual method for optimization of the objective function.
           *
           * @param[out] xOpt optimal point
           * @return          optimal objective function value
           */
          virtual float_t optimize(base::DataVector& xOpt) = 0;

          /**
           * @return objective function
           */
          ObjectiveFunction& getObjectiveFunction() const {
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

        protected:
          /// objective function
          ObjectiveFunction& f;
          /// maximal number of iterations or function evaluations
          size_t N;
          /// starting point
          base::DataVector x0;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_OPTIMIZER_HPP */
