// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>
#include <sgpp/solver/ode/VarTimestep.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/tools/GridPrinter.hpp>
#include <sgpp/base/exception/solver_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

namespace sgpp {
namespace solver {

VarTimestep::VarTimestep(std::string pred, std::string corr, size_t imax, double timestepSize,
                         double eps, sgpp::base::ScreenOutput* screen, double gamma)
    : StepsizeControl(imax, timestepSize, eps, 1.0, screen, gamma),
      _predictor(pred),
      _corrector(corr) {
  std::stringstream fnsstream;

  fnsstream << "Time_"
            << "VaTim" << eps << ".gnuplot";

  filename = fnsstream.str();
}

VarTimestep::~VarTimestep() {}

void VarTimestep::predictor(SLESolver& LinearSystemSolver,
                            sgpp::solver::OperationParabolicPDESolverSystem& System,
                            double tmp_timestepsize, sgpp::base::DataVector& dv,
                            sgpp::base::DataVector& corr, sgpp::base::DataVector* rhs) {
  System.setTimestepSize(tmp_timestepsize);

  System.setODESolver("AdBas");

  // generate right hand side
  rhs = System.generateRHS();

  // solve the system of the current timestep
  LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false, -1.0);

  System.finishTimestep();

  dv.resize(System.getGridCoefficients()->getSize());

  System.getGridCoefficientsForSC(dv);

  System.abortTimestep();
}

void VarTimestep::corrector(SLESolver& LinearSystemSolver,
                            sgpp::solver::OperationParabolicPDESolverSystem& System,
                            double tmp_timestepsize, sgpp::base::DataVector& dv,
                            sgpp::base::DataVector* rhs) {
  System.setODESolver("CrNic");

  // generate right hand side
  rhs = System.generateRHS();

  // solve the system of the current timestep
  LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false, -1.0);

  System.finishTimestep();

  System.getGridCoefficientsForSC(dv);
}

double VarTimestep::nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old,
                                  double norm, double epsilon) {
  double deltaY = norm / (3.0 * (1.0 + tmp_timestepsize / tmp_timestepsize_old));

  return tmp_timestepsize *
         std::max(0.67, std::min(1.5, pow(epsilon / deltaY, 1.0 / 3.0)));
}
}  // namespace solver
}  // namespace sgpp
