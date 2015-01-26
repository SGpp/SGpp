// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

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

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace solver {

    StepsizeControlEJ::StepsizeControlEJ(std::string odesolver, size_t nTimesteps, double timestepSize, double eps, double sc, SGPP::base::ScreenOutput* screen, double gamma)
      : StepsizeControl(nTimesteps, timestepSize, eps, sc, screen, gamma), _odesolver(odesolver) {
      std::stringstream fnsstream;
      fnsstream << "Time_" << "SCEJ" << this->myEps << "_" << this->mySC << ".gnuplot";
      filename = fnsstream.str();
    }

    StepsizeControlEJ::~StepsizeControlEJ() {
    }

    void StepsizeControlEJ::predictor(SLESolver& LinearSystemSolver, SGPP::pde::OperationParabolicPDESolverSystem& System,
                                      double tmp_timestepsize, SGPP::base::DataVector& dv, SGPP::base::DataVector& corr, SGPP::base::DataVector* rhs) {
      //pred()
      dv.resize(corr.getSize());
      dv.setAll(0.0);
      dv.add(corr);
    }

    void StepsizeControlEJ::corrector(SLESolver& LinearSystemSolver, SGPP::pde::OperationParabolicPDESolverSystem& System, double tmp_timestepsize, SGPP::base::DataVector& dv, SGPP::base::DataVector* rhs) {
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

    double StepsizeControlEJ::norm(SGPP::pde::OperationParabolicPDESolverSystem& System, SGPP::base::DataVector& dv1, SGPP::base::DataVector& dv2) {
      return maxNorm(System, dv1, dv2);
    }

    double StepsizeControlEJ::nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old, double norm, double epsilon) {
      return tmp_timestepsize * epsilon / norm;
    }
  }
}