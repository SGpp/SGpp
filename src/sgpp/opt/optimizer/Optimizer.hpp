/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_OPTIMIZER_OPTIMIZER_HPP
#define SGPP_OPT_OPTIMIZER_OPTIMIZER_HPP

#include "opt/function/Objective.hpp"

#include <vector>
#include <cstddef>

#include <iostream>

namespace sg {
  namespace opt {
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
           * The starting point is set to \f$(0.5, \dotsc, 0.5)^{\mathrm{T}}\f$.
           *
           * @param f     function to optimize
           * @param N     maximal number of iterations or function evaluations
           *              (depending on the implementation)
           */
          Optimizer(function::Objective& f, size_t N = DEFAULT_N) :
            f(f.clone()), N(N), x0(std::vector<double>(f.getDimension(), 0.5)) {
          }

          /**
           * Virtual destructor.
           */
          virtual ~Optimizer() {
          }

          /**
           * Pure virtual method for optimization of the objective function.
           *
           * @param[out] xopt optimal point
           * @return          optimal objective function value
           */
          virtual double optimize(std::vector<double>& xopt) = 0;

          /**
           * Pure virtual method for cloning the optimizer.
           * It should return a pointer to the cloned object and it's used for parallel computations.
           *
           * @return smart pointer to cloned object
           */
          virtual tools::SmartPointer<Optimizer> clone() = 0;

          /**
           * @return smart pointer to objective function
           */
          const tools::SmartPointer<function::Objective>& getObjectiveFunction() const {
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
          const std::vector<double>& getStartingPoint() const {
            return x0;
          }

          /**
           * @param starting_point    starting point
           */
          void setStartingPoint(const std::vector<double>& x0) {
            this->x0 = x0;
          }

        protected:
          /// objective function
          tools::SmartPointer<function::Objective> f;
          /// maximal number of iterations or function evaluations
          size_t N;
          /// starting point
          std::vector<double> x0;
      };

    }
  }
}

#endif
