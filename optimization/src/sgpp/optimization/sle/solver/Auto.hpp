// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_SLE_SOLVER_AUTO_HPP
#define SGPP_OPTIMIZATION_SLE_SOLVER_AUTO_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/solver/Solver.hpp>

#include <vector>

namespace SGPP {
  namespace optimization {
    namespace sle {
      namespace solver {

        /**
         * Automatic choice of external linear solver.
         */
        class Auto : public Solver {
          public:
            /// maximal matrix dimension to allow use of full solvers
            static const size_t MAX_DIM_FOR_FULL = 30000;
            /// maximal ratio of non-zero entries for sparse solvers
            static const float_t MAX_NNZ_RATIO_FOR_SPARSE;
            /// maximal ratio of non-zero entries to prefer Gmm++ over UMFPACK
            static const float_t MAX_NNZ_RATIO_FOR_GMMPP;
            /// ratio of the rows (e.g. 0.1 = 10%) to use for estimation of the sparsity ratio
            static const float_t ESTIMATE_NNZ_ROWS_SAMPLE_SIZE;

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
