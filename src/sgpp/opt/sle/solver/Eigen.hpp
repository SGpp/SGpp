/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_SLE_SOLVER_EIGEN_HPP
#define SGPP_OPT_SLE_SOLVER_EIGEN_HPP

#include "opt/sle/solver/Solver.hpp"

#include <vector>

namespace sg {
  namespace opt {
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
