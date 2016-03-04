// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PARABOLICPDESOLVER_HPP
#define PARABOLICPDESOLVER_HPP

#include <sgpp/pde/application/PDESolver.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * This class extends the PDESolver with functions that are needed to
 * solve parabolic PDEs
 *
 */
class ParabolicPDESolver : public PDESolver {
 protected:
  /// the size of one timestep
  // double timestepSize;
  /// The number of timesteps that are executed during solving
  // size_t nTimesteps;

 public:
  /**
   * Std-Constructor of the solver
   */
  ParabolicPDESolver();

  /**
   * Std-Destructor of the solver
   */
  virtual ~ParabolicPDESolver();

  /**
   * Call this routine to use an explicit Euler algorithm to solve the parabolic PDE
   *
   * @param numTimesteps the number of timesteps that should be executed
   * @param timestepsize the size of the interval one timestep moves forward
   * @param maxCGIterations the maximum of interation in the CG solver
   * @param epsilonCG the epsilon used in the CG
   * @param alpha the coefficients of the Sparse Gird's basis functions
   * @param verbose enables verbose output during solving
   * @param generateAnimation set this to true, if you want to generate a grid output in every
   * timestep
   * @param numEvalsAnimation specifies the evaluation per dimension when a animation is created
   */
  virtual void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations,
                                  double epsilonCG, sgpp::base::DataVector& alpha,
                                  bool verbose = false, bool generateAnimation = false,
                                  size_t numEvalsAnimation = 20) = 0;

  /**
   * Call this routine to use an explicit Euler algorithm to solve the parabolic PDE
   *
   * @param numTimesteps the number of timesteps that should be executed
   * @param timestepsize the size of the interval one timestep moves forward
   * @param maxCGIterations the maximum of interation in the CG solver
   * @param epsilonCG the epsilon used in the CG
   * @param alpha the coefficients of the Sparse Gird's basis functions
   * @param verbose enables verbose output during solving
   * @param generateAnimation set this to true, if you want to generate a grid output in every
   * timestep
   * @param numEvalsAnimation specifies the evaluation per dimension when a animation is created
   */
  virtual void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations,
                                  double epsilonCG, sgpp::base::DataVector& alpha,
                                  bool verbose = false, bool generateAnimation = false,
                                  size_t numEvalsAnimation = 20) = 0;

  /**
   * Call this routine to use the Crank Nicolson algorithm to solve the parabolic PDE
   *
   * @param numTimesteps the number of timesteps that should be executed
   * @param timestepsize the size of the interval one timestep moves forward
   * @param maxCGIterations the maximum of interation in the CG solver
   * @param epsilonCG the epsilon used in the CG
   * @param alpha the coefficients of the Sparse Gird's basis functions
   * @param NumImEul specifies how many ImEul steps should be executed before CrNic is used, default
   * is 0
   */
  virtual void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations,
                                  double epsilonCG, sgpp::base::DataVector& alpha,
                                  size_t NumImEul = 0) = 0;
};
}  // namespace pde
}  // namespace sgpp

#endif /* PARABOLICPDESOLVER_HPP */
