// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_SLE_SOLVER_EIGEN_HPP
#define SGPP_OPTIMIZATION_SLE_SOLVER_EIGEN_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/solver/Solver.hpp>

#include <vector>

namespace SGPP {
  namespace optimization {
    namespace sle {
      namespace solver {

        /**
         * Linear system solver using Eigen (direct full solver).
         */
        class Eigen : public Solver {
          public:
            /**
             * @param       system  system to be solved
             * @param       b       right-hand side
             * @param[out]  x       solution to the system
             * @return              whether all went well (false if errors occurred)
             */
            bool solve(system::System& system, const std::vector<float_t>& b, std::vector<float_t>& x) const;

            /**
             * @param       system  system to be solved
             * @param       B       vector of right-hand sides
             * @param[out]  X       vector of solutions to the systems
             * @return              whether all went well (false if errors occurred)
             */
            bool solve(system::System& system, const std::vector<std::vector<float_t> >& B,
                       std::vector<std::vector<float_t> >& X) const;
        };

      }
    }
  }
}

#endif
