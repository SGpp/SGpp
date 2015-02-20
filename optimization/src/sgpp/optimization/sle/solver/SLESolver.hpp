// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_SLE_SOLVER_SLESOLVER_HPP
#define SGPP_OPTIMIZATION_SLE_SOLVER_SLESOLVER_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/system/SLE.hpp>

#include <vector>

namespace SGPP {
  namespace optimization {
    namespace sle_solver {

      /**
       * Abstract class for solving systems of linear equations.
       */
      class SLESolver {
        public:
          /**
           * Constructor.
           */
          SLESolver() {
          }

          /**
           * Virtual destructor.
           */
          virtual ~SLESolver() {
          }

          /**
           * Pure virtual method for a solving linear system.
           *
           * @param       system  system to be solved
           * @param       b       right-hand side
           * @param[out]  x       solution to the system
           * @return              whether all went well
           *                      (false if errors occurred)
           */
          virtual bool solve(SLE& system, const std::vector<float_t>& b,
                             std::vector<float_t>& x) const = 0;

          /**
           * Virtual method for solving multiple linear systems with
           * different right-hand sides.
           * Defaults to calling the solve() method for a single
           * right-hand side multiple times.
           *
           * @param       system  system to be solved
           * @param       B       vector of right-hand sides
           * @param[out]  X       vector of solutions to the systems
           * @return              whether all went well
           *                      (false if errors occurred)
           */
          virtual bool solve(SLE& system,
                             const std::vector<std::vector<float_t>>& B,
                             std::vector<std::vector<float_t>>& X) const {
            std::vector<float_t> x;
            X.clear();

            for (size_t i = 0; i < B.size(); i++) {
              const std::vector<float_t>& b = B[i];

              if (solve(system, b, x)) {
                X.push_back(x);
              } else {
                X.clear();
                return false;
              }
            }

            return true;
          }
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_SLE_SOLVER_SLESOLVER_HPP */
