/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/algorithm/LaserHeatEquationParabolicPDESolverSystemParallelOMP2D.hpp"
#include "pde/application/LaserHeatEquationSolver2D.hpp"
#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "solver/sle/ConjugateGradients.hpp"
#include "base/grid/Grid.hpp"
#include "base/exception/application_exception.hpp"
#include "base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "base/grid/generation/functors/SurplusCoarseningFunctor.hpp"
#include "stdlib.h"
#include "base/operation/BaseOpFactory.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/tools/SGppStopwatch.hpp"
#include <sstream>
#include <fstream>

namespace sg {
  namespace pde {

    LaserHeatEquationSolver2D::LaserHeatEquationSolver2D(double beam_velocity, double heat_sigma, sg::base::GridStorage::index_type::level_type max_level, double refine_threshold, double coarsen_threshold, double heat) : HeatEquationSolver(), beam_velocity_(beam_velocity), heat_sigma_(heat_sigma), max_level_(max_level), heat_(heat), refine_threshold_(refine_threshold), coarsen_threshold_(coarsen_threshold) {
    }

    LaserHeatEquationSolver2D::~LaserHeatEquationSolver2D() {
    }

    void LaserHeatEquationSolver2D::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, base::DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation) {
      throw new base::application_exception("LaserHeatEquationSolver::solveExplicitEuler : explicit sg::solver::Euler is not supported!");
    }

    void LaserHeatEquationSolver2D::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, base::DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation) {
      if (this->bGridConstructed) {
        this->myScreen->writeStartSolve("2D Laser Heat Solver");
        double dNeededTime;
        sg::solver::Euler* myEuler = new sg::solver::Euler("ImEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, this->myScreen);
        sg::solver::ConjugateGradients* myCG = new sg::solver::ConjugateGradients(maxCGIterations, epsilonCG);
        LaserHeatEquationParabolicPDESolverSystemParallelOMP2D* myHESolver = new LaserHeatEquationParabolicPDESolverSystemParallelOMP2D(this->beam_velocity_, this->heat_sigma_, this->max_level_, this->heat_, this->refine_threshold_, this->coarsen_threshold_, *this->myGrid, alpha, this->a, timestepsize, "ImEul");
        base::SGppStopwatch* myStopwatch = new base::SGppStopwatch();

        myStopwatch->start();
        myEuler->solve(*myCG, *myHESolver, verbose);
        dNeededTime = myStopwatch->stop();

        if (this->myScreen != NULL) {
          std::cout << "Time to solve: " << dNeededTime << " seconds" << std::endl;
          this->myScreen->writeEmptyLines(2);
        }

        delete myStopwatch;
        delete myHESolver;
        delete myCG;
        delete myEuler;
      } else {
        throw new base::application_exception("LaserHeatEquationSolver::solveImplicitEuler : A grid wasn't constructed before!");
      }
    }

    void LaserHeatEquationSolver2D::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, base::DataVector& alpha, size_t NumImEul) {
      throw new base::application_exception("LaserHeatEquationSolver::solveCrankNicolson : Crank Nicolson is not supported!");
    }

    void LaserHeatEquationSolver2D::refineInitialGridWithLaserHeat(base::DataVector& alpha, size_t nRefinements) {
      if (this->bGridConstructed) {
        base::StdNormalDistribution myNormDistr;

        std::cout << std::endl;

        for (size_t r = 0; r < nRefinements; r++) {
          double* dblFuncValues = new double[dim];

          for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
            std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*(this->myBoundingBox));
            std::stringstream coordsStream(coords);

            for (size_t j = 0; j < this->dim; j++) {
              coordsStream >> dblFuncValues[j];
            }

            // check if coordinates at starting point of laser
            alpha[i] =  this->heat_ * (myNormDistr.getDensity(dblFuncValues[0], 0.25, this->heat_sigma_) * myNormDistr.getDensity(dblFuncValues[1], 0.5, this->heat_sigma_));

            //boundaries are set to zero
            if (dblFuncValues[0] == 0.0 || dblFuncValues[1] == 0.0) {
              alpha[i] = 0.0;
            }
          }

          delete[] dblFuncValues;

          // do hierarchisation
          base::OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
          myHierarchisation->doHierarchisation(alpha);
          delete myHierarchisation;

          // do refinement
          base::GridGenerator* myGenerator = this->myGrid->createGridGenerator();
          size_t numRefines = myGenerator->getNumberOfRefinablePoints();
          base::SurplusRefinementFunctor* myRefineFunc = new base::SurplusRefinementFunctor(&alpha, numRefines, this->refine_threshold_);
          myGenerator->refineMaxLevel(myRefineFunc, this->max_level_);
          delete myRefineFunc;
          delete myGenerator;

          std::cout << "Refining grid..." << std::endl;
          std::cout << "New inner grid size): " << this->myGrid->getStorage()->getNumInnerPoints() << std::endl;

          alpha.resize(this->myGridStorage->size());
        }

        std::cout << std::endl << std::endl;

        // Init last grid
        double* dblFuncValues = new double[dim];

        for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
          std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*(this->myBoundingBox));
          std::stringstream coordsStream(coords);

          for (size_t j = 0; j < this->dim; j++) {
            coordsStream >> dblFuncValues[j];
          }

          // check if coordinates at starting point of laser
          alpha[i] =  this->heat_ * (myNormDistr.getDensity(dblFuncValues[0], 0.25, this->heat_sigma_) * myNormDistr.getDensity(dblFuncValues[1], 0.5, this->heat_sigma_));

          //boundaries are set to zero
          if (dblFuncValues[0] == 0.0 || dblFuncValues[1] == 0.0) {
            alpha[i] = 0.0;
          }
        }

        delete[] dblFuncValues;

        // do hierarchisation
        base::OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
        myHierarchisation->doHierarchisation(alpha);
        delete myHierarchisation;
      } else {
        throw new base::application_exception("LaserHeatEquationSolver::refineInitialGridWithPayoff : The grid wasn't initialized before!");
      }
    }

    void LaserHeatEquationSolver2D::initScreen() {
      this->myScreen = new base::ScreenOutput();
      this->myScreen->writeTitle("SGpp - Laser Heat Solver, 1.0.0", "Alexander Heinecke, (C) 2011");
    }

  }
}
