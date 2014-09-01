/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_SLE_SOLVER_AUTO_HPP
#define SGPP_OPT_SLE_SOLVER_AUTO_HPP

#include "opt/sle/solver/Solver.hpp"

#include <vector>

namespace sg {
  namespace opt {
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
            static const double MAX_NNZ_RATIO_FOR_SPARSE;
            /// maximal ratio of non-zero entries to prefer Gmm++ over UMFPACK
            static const double MAX_NNZ_RATIO_FOR_GMMPP;
            /// ratio of the rows (e.g. 0.1 = 10%) to use for estimation of the sparsity ratio
            static const double ESTIMATE_NNZ_ROWS_SAMPLE_SIZE;

            /**
             * @param       system  system to be solved
             * @param       b       right-hand side
             * @param[out]  x       solution to the system
             * @return              whether all went well (false if errors occurred)
             */
            bool solve(system::System& system, const std::vector<double>& b, std::vector<double>& x) const;

            /**
             * @param       system  system to be solved
             * @param       B       vector of right-hand sides
             * @param[out]  X       vector of solutions to the systems
             * @return              whether all went well (false if errors occurred)
             */
            bool solve(system::System& system, const std::vector<std::vector<double> >& B,
                       std::vector<std::vector<double> >& X) const;
        };

      }
    }
  }
}

#endif
