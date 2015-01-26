/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Peter Hoffmann (peter.hoffmann@mytum.de)

#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>
#include <sgpp/solver/ode/StepsizeControlEJ.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/OperationEval.hpp>
#include <sgpp/base/tools/GridPrinter.hpp>
#include <sgpp/base/exception/solver_exception.hpp>


#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

namespace sg {
  namespace solver {

    StepsizeControlEJ::StepsizeControlEJ(std::string odesolver, size_t nTimesteps, double timestepSize, double eps, double sc, sg::base::ScreenOutput* screen, double gamma)
      : StepsizeControl(nTimesteps, timestepSize, eps, sc, screen, gamma), _odesolver(odesolver) {
      std::stringstream fnsstream;
      fnsstream << "Time_" << "SCEJ" << this->myEps << "_" << this->mySC << ".gnuplot";
      filename = fnsstream.str();
    }

    StepsizeControlEJ::~StepsizeControlEJ() {
    }

    void StepsizeControlEJ::predictor(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System,
                                      double tmp_timestepsize, sg::base::DataVector& dv, sg::base::DataVector& corr, sg::base::DataVector* rhs) {
      //pred()
      dv.resize(corr.getSize());
      dv.setAll(0.0);
      dv.add(corr);
    }

    void StepsizeControlEJ::corrector(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System, double tmp_timestepsize, sg::base::DataVector& dv, sg::base::DataVector* rhs) {
      //corr()

      System.setODESolver(_odesolver);

      System.setTimestepSize(tmp_timestepsize);

      // generate right hand side
      rhs = System.generateRHS();

      // solve the system of the current timestep
      LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false, -1.0);

      System.finishTimestep();
      dv.resize(System.getGridCoefficients()->getSize());
      System.getGridCoefficientsForSC(dv);

      // end corr()
    }

    double StepsizeControlEJ::norm(sg::pde::OperationParabolicPDESolverSystem& System, sg::base::DataVector& dv1, sg::base::DataVector& dv2) {
      return maxNorm(System, dv1, dv2);
    }

    double StepsizeControlEJ::nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old, double norm, double epsilon) {
      return tmp_timestepsize * epsilon / norm;
    }
  }
}
