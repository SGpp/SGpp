// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STEPSIZECONTROL_HPP
#define STEPSIZECONTROL_HPP

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/solver/ODESolver.hpp>
#include <sgpp/solver/operation/hash/OperationParabolicPDESolverSystem.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace solver {

/**
 * This class implements a step size control using Adams-Bashforth and Crank-Nicolson
 * for solving ordinary partial equations
 *
 */
class StepsizeControl : public ODESolver {
 private:
  /// identify of grid coarsening is used
  bool useCoarsen;

 protected:
  /// Pointer to sgpp::base::ScreenOutput object
  sgpp::base::ScreenOutput* myScreen;

  /// temp. Stepsize Control
  double mySC;

  /// epsilon for the step size control
  double myEps;

  virtual void predictor(SLESolver& LinearSystemSolver,
                         sgpp::solver::OperationParabolicPDESolverSystem& System,
                         double tmp_timestepsize, sgpp::base::DataVector& dv,
                         sgpp::base::DataVector& corr, sgpp::base::DataVector* rhs) = 0;
  virtual void corrector(SLESolver& LinearSystemSolver,
                         sgpp::solver::OperationParabolicPDESolverSystem& System,
                         double tmp_timestepsize, sgpp::base::DataVector& dv,
                         sgpp::base::DataVector* rhs) = 0;

  virtual double norm(sgpp::solver::OperationParabolicPDESolverSystem& System,
                       sgpp::base::DataVector& dv1, sgpp::base::DataVector& dv2);

  virtual double nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old, double norm,
                               double epsilon) = 0;

  double twoNorm(sgpp::solver::OperationParabolicPDESolverSystem& System,
                  sgpp::base::DataVector& dv1, sgpp::base::DataVector& dv2);

  double maxNorm(sgpp::solver::OperationParabolicPDESolverSystem& System,
                  sgpp::base::DataVector& dv1, sgpp::base::DataVector& dv2);

  std::string filename;

  /// damping factor
  double _gamma;

 public:
  /**
   * Std-Constructer
   *
   * @param sc step size
   * @param imax number of maximum executed iterations
   * @param timestepSize the size of one timestep
   * @param eps the epsilon for the step size control
   * @param screen possible pointer to a sgpp::base::ScreenOutput object
   * @param gamma damping factor
   */
  StepsizeControl(size_t imax, double timestepSize, double eps, double sc,
                  sgpp::base::ScreenOutput* screen = nullptr, double gamma = 0.5);

  /**
   * Std-Destructor
   */
  virtual ~StepsizeControl();

  void solve(SLESolver& LinearSystemSolver, sgpp::solver::OperationParabolicPDESolverSystem& System,
             bool bIdentifyLastStep = false, bool verbose = false);
};

}  // namespace solver
}  // namespace sgpp

#endif /* STEPSIZECONTROL_H_ */
