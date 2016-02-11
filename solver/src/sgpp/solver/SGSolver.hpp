// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGSOLVER_HPP
#define SGSOLVER_HPP

#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>


namespace SGPP {
namespace solver {

/**
 * Abstract class that defines a solver used in Sparse Grids
 * Applications
 */
class SGSolver {
 protected:
  /// Number of Iterations needed for the solve
  size_t nIterations;
  /// Number of maximum iterations for cg
  size_t nMaxIterations;
  /// residuum
  float_t residuum;
  /// epsilon needed in the, e.g. final error in the iterative solver, or a timestep
  float_t myEpsilon;

 public:
  /**
   * Std-Constructor
   *
   * @param nMaximumIterations number of maximum executed iterations
   * @param epsilon the final error in the iterative solver, or the size of one timestep
   */
  SGSolver(size_t nMaximumIterations,
           float_t epsilon) : nMaxIterations(nMaximumIterations), myEpsilon(epsilon) {
    nIterations = 0;
    residuum = 0.0;
  }

  /**
   * Std-Destructor
   */
  virtual ~SGSolver() { }


  /**
   * function that returns the number of needed solve steps
   *
   * @return the number of needed solve steps of the sovler
   */
  size_t getNumberIterations() {
    return nIterations;
  }

  /**
   * function the returns the residuum (current or final), error of the solver
   *
   * @return the residuum
   */
  float_t getResiduum() {
    return residuum;
  }

  /**
   * resets the number of maximum iterations
   *
   * @param nIterations the new number of maximum iterations
   */
  void setMaxIterations(size_t nIterations) {
    nMaxIterations = nIterations;
  }

  /**
   * resets the epsilon, that is used in the SGSolver
   *
   * @param eps the new value of epsilon
   */
  void setEpsilon(float_t eps) {
    myEpsilon = eps;
  }

  /**
   * gets the the epsilon, that is used in the SGSolver
   *
   * @return the epsilon, used in the solver
   */
  float_t getEpsilon() {
    return myEpsilon;
  }
};

}  // namespace solver
}  // namespace SGPP

#endif /* SGSOLVER_HPP */
