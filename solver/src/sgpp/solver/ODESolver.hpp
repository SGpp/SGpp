// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ODESOLVER_HPP
#define ODESOLVER_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/solver/operation/hash/OperationParabolicPDESolverSystem.hpp>

#include <sgpp/solver/SGSolver.hpp>
#include <sgpp/solver/SLESolver.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace solver {

class ODESolver : public SGSolver {
 public:
  /**
   * Std-Constructor
   *
   * @param imax number of maximum executed iterations
   * @param timestepSize the size of one timestep
   */
  ODESolver(size_t imax, float_t timestepSize) : SGSolver(imax, timestepSize) {
  }

  /**
   * Std-Destructor
   */
  virtual ~ODESolver() { }

  /**
   * Pure virtual Function that defines a solve method for an ODE solver
   *
   * @param LinearSystemSolver reference to an instance of a linear system solver that is used by this ODE solver
   * @param System reference to an SGPP::base::OperationMatrix Object that implements the matrix vector multiplication
   * @param bIdentifyLastStep set this to true to tell System the last step
   * @param verbose prints information during execution of the solver
   */
  virtual void solve(SLESolver& LinearSystemSolver,
                     SGPP::solver::OperationParabolicPDESolverSystem& System,
                     bool bIdentifyLastStep = false, bool verbose = false) = 0;
};

}  // namespace solver
}  // namespace SGPP

#endif /* ODESOLVER_HPP */
