/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Peter Hoffmann (peter.hoffmann@mytum.de)

#include "base/grid/common/DirichletUpdateVector.hpp"
#include "solver/ode/VarTimestep.hpp"
#include "base/operation/OperationEval.hpp"
#include "base/tools/GridPrinter.hpp"
#include "base/exception/solver_exception.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

namespace sg {
  namespace solver {

    VarTimestep::VarTimestep(std::string pred, std::string corr, size_t imax, double timestepSize, double eps, sg::base::ScreenOutput* screen, double gamma)
      : StepsizeControl(imax, timestepSize, eps, 1.0, screen, gamma), _predictor(pred), _corrector(corr) {

      std::stringstream fnsstream;

      fnsstream << "Time_" << "VaTim" << eps << ".gnuplot";

      filename = fnsstream.str();

    }

    VarTimestep::~VarTimestep() {
    }

    void VarTimestep::predictor(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System,
                                double tmp_timestepsize, sg::base::DataVector& dv, sg::base::DataVector& corr, sg::base::DataVector* rhs) {
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

    void VarTimestep::corrector(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System, double tmp_timestepsize, sg::base::DataVector& dv, sg::base::DataVector* rhs) {
      System.setODESolver("CrNic");

      // generate right hand side
      rhs = System.generateRHS();

      // solve the system of the current timestep
      LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false, -1.0);

      System.finishTimestep();

      System.getGridCoefficientsForSC(dv);
    }


    double VarTimestep::nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old, double norm, double epsilon) {

      double deltaY = norm / (3.0 * (1.0 + tmp_timestepsize / tmp_timestepsize_old));

      return tmp_timestepsize * std::max(0.67, std::min(1.5, pow(epsilon / deltaY, (double)1.0 / (double)3.0)));

    }
  }
}
