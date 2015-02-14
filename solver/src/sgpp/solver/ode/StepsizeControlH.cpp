// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>
#include <sgpp/solver/ode/StepsizeControlH.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/tools/GridPrinter.hpp>
#include <sgpp/base/exception/solver_exception.hpp>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace solver {

    StepsizeControlH::StepsizeControlH(std::string odesolver, size_t imax, float_t timestepSize, float_t eps, SGPP::base::ScreenOutput* screen, float_t gamma)
      : StepsizeControl(imax, timestepSize, eps, 1.0, screen, gamma), _odesolver(odesolver) {
      this->residuum = 0.0;
      this->myEps = eps;
      std::stringstream fnsstream;
      fnsstream << "Time_" << "SCH" << this->myEps << ".gnuplot";
      filename = fnsstream.str();
    }

    StepsizeControlH::~StepsizeControlH() {
    }

    void StepsizeControlH::predictor(SLESolver& LinearSystemSolver, SGPP::pde::OperationParabolicPDESolverSystem& System,
                                     float_t tmp_timestepsize, SGPP::base::DataVector& dv, SGPP::base::DataVector& corr, SGPP::base::DataVector* rhs) {
      System.setODESolver(_odesolver);
      System.setTimestepSize(tmp_timestepsize);

      // generate right hand side
      rhs = System.generateRHS();

      // solve the system of the current timesteps
      LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false, -1.0);

      System.finishTimestep();

      dv.resize(System.getGridCoefficients()->getSize());
      System.getGridCoefficientsForSC(dv);

      System.abortTimestep();
    }

    void StepsizeControlH::corrector(SLESolver& LinearSystemSolver, SGPP::pde::OperationParabolicPDESolverSystem& System, float_t tmp_timestepsize, SGPP::base::DataVector& dv, SGPP::base::DataVector* rhs) {
      System.setODESolver(_odesolver);
      System.setTimestepSize(tmp_timestepsize / 2.0);

      // generate right hand side
      rhs = System.generateRHS();

      // solve the system of the current timesteps
      LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false, -1.0);
      System.finishTimestep();

      rhs = System.generateRHS();

      LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false, -1.0);
      System.finishTimestep();

      dv.resize(System.getGridCoefficients()->getSize());
      System.getGridCoefficientsForSC(dv);
    }


    float_t StepsizeControlH::nextTimestep(float_t tmp_timestepsize, float_t tmp_timestepsize_old, float_t norm, float_t epsilon) {
      float_t deltaY = 3.0 * norm / 4.0;

      return tmp_timestepsize * pow(epsilon / deltaY, (float_t)1.0 / (float_t)3.0);

    }
  }
}