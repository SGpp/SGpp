/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp"
#include "pde/algorithm/HeatEquationParabolicPDESolverSystemParallelOMP.hpp"
#include "pde/application/HeatEquationSolver.hpp"
#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "solver/sle/ConjugateGradients.hpp"
#include "base/grid/Grid.hpp"
#include "base/exception/application_exception.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "stdlib.h"
#include <sstream>
#include <fstream>

using namespace sg::solver;
using namespace sg::base;

namespace sg {
  namespace pde {

    HeatEquationSolver::HeatEquationSolver() : ParabolicPDESolver() {
      this->bGridConstructed = false;
      this->myScreen = NULL;
    }

    HeatEquationSolver::~HeatEquationSolver() {
      if (this->myScreen != NULL) {
        delete this->myScreen;
      }
    }

    void HeatEquationSolver::constructGrid(BoundingBox& BoundingBox, int level) {
      this->dim = BoundingBox.getDimensions();
      this->levels = level;

      this->myGrid = new LinearTrapezoidBoundaryGrid(BoundingBox);

      GridGenerator* myGenerator = this->myGrid->createGridGenerator();
      myGenerator->regular(levels);
      delete myGenerator;

      this->myBoundingBox = this->myGrid->getBoundingBox();
      this->myGridStorage = this->myGrid->getStorage();

      this->bGridConstructed = true;
    }

    void HeatEquationSolver::setHeatCoefficient(double a) {
      this->a = a;
    }

    void HeatEquationSolver::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation) {
      if (this->bGridConstructed) {
        this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
        double dNeededTime;
        Euler* myEuler = new Euler("ExEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, this->myScreen);
        ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
#ifdef _OPENMP
        HeatEquationParabolicPDESolverSystemParallelOMP* myHESolver = new HeatEquationParabolicPDESolverSystemParallelOMP(*this->myGrid, alpha, this->a, timestepsize, "ExEul");
#else
        HeatEquationParabolicPDESolverSystem* myHESolver = new HeatEquationParabolicPDESolverSystem(*this->myGrid, alpha, this->a, timestepsize, "ExEul");
#endif
        SGppStopwatch* myStopwatch = new SGppStopwatch();

        myStopwatch->start();
        myEuler->solve(*myCG, *myHESolver, verbose);
        dNeededTime = myStopwatch->stop();

        if (this->myScreen != NULL) {
          std::cout << "Time to solve: " << dNeededTime << " seconds" << std::endl;
          this->myScreen->writeEmptyLines(2);
        }

        delete myStopwatch;
        delete myCG;
        delete myEuler;
      } else {
        throw new application_exception("HeatEquationSolver::solveExplicitEuler : A grid wasn't constructed before!");
      }
    }

    void HeatEquationSolver::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation) {
      if (this->bGridConstructed) {
        this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
        double dNeededTime;
        Euler* myEuler = new Euler("ImEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, this->myScreen);
        ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
#ifdef _OPENMP
        HeatEquationParabolicPDESolverSystemParallelOMP* myHESolver = new HeatEquationParabolicPDESolverSystemParallelOMP(*this->myGrid, alpha, this->a, timestepsize, "ImEul");
#else
        HeatEquationParabolicPDESolverSystem* myHESolver = new HeatEquationParabolicPDESolverSystem(*this->myGrid, alpha, this->a, timestepsize, "ImEul");
#endif
        SGppStopwatch* myStopwatch = new SGppStopwatch();

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
        throw new application_exception("HeatEquationSolver::solveImplicitEuler : A grid wasn't constructed before!");
      }
    }

    void HeatEquationSolver::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul) {
      if (this->bGridConstructed) {
        this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
        double dNeededTime;
        ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
#ifdef _OPENMP
        HeatEquationParabolicPDESolverSystemParallelOMP* myHESolver = new HeatEquationParabolicPDESolverSystemParallelOMP(*this->myGrid, alpha, this->a, timestepsize, "CrNic");
#else
        HeatEquationParabolicPDESolverSystem* myHESolver = new HeatEquationParabolicPDESolverSystem(*this->myGrid, alpha, this->a, timestepsize, "CrNic");
#endif
        SGppStopwatch* myStopwatch = new SGppStopwatch();

        size_t numCNSteps;
        size_t numIESteps;

        numCNSteps = numTimesteps;

        if (numTimesteps > NumImEul) {
          numCNSteps = numTimesteps - NumImEul;
        }

        numIESteps = NumImEul;

        Euler* myEuler = new Euler("ImEul", numIESteps, timestepsize, false, 0, this->myScreen);
        CrankNicolson* myCN = new CrankNicolson(numCNSteps, timestepsize);

        myStopwatch->start();

        if (numIESteps > 0) {
          myEuler->solve(*myCG, *myHESolver, false);
        }

        myCN->solve(*myCG, *myHESolver, false);
        dNeededTime = myStopwatch->stop();

        if (this->myScreen != NULL) {
          std::cout << "Time to solve: " << dNeededTime << " seconds" << std::endl;
          this->myScreen->writeEmptyLines(2);
        }

        delete myStopwatch;
        delete myHESolver;
        delete myCG;
        delete myCN;
        delete myEuler;
      } else {
        throw new application_exception("HeatEquationSolver::solveCrankNicolson : A grid wasn't constructed before!");
      }
    }

    void HeatEquationSolver::initGridWithSmoothHeat(DataVector& alpha, double mu, double sigma, double factor) {
      if (this->bGridConstructed) {
        double tmp;
        double* dblFuncValues = new double[this->dim];

        for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
          std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
          std::stringstream coordsStream(coords);

          for (size_t j = 0; j < this->dim; j++) {
            coordsStream >> tmp;

            dblFuncValues[j] = tmp;
          }

          tmp = 1.0;

          for (size_t j = 0; j < this->dim; j++) {
            tmp *=  factor * factor * ((1.0 / (sigma * 2.0 * 3.145)) * exp((-0.5) * ((dblFuncValues[j] - mu) / sigma) * ((dblFuncValues[j] - mu) / sigma)));
          }

          alpha[i] = tmp;
        }

        delete[] dblFuncValues;

        OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
        myHierarchisation->doHierarchisation(alpha);
        delete myHierarchisation;
      } else {
        throw new application_exception("HeatEquationSolver::initGridWithSmoothHeat : A grid wasn't constructed before!");
      }
    }

    void HeatEquationSolver::initScreen() {
      this->myScreen = new ScreenOutput();
      this->myScreen->writeTitle("SGpp - Heat Equation Solver, 1.0.0", "Alexander Heinecke, (C) 2009-2011");
    }

    void HeatEquationSolver::storeInnerMatrix(DataVector& alpha, std::string tFilename, double timestepsize) {
      if (this->bGridConstructed) {
        HeatEquationParabolicPDESolverSystem* myHESolver = new HeatEquationParabolicPDESolverSystem(*this->myGrid, alpha, this->a, timestepsize, "ImEul");
        SGppStopwatch* myStopwatch = new SGppStopwatch();

        std::string mtx = "";

        myStopwatch->start();
        std::cout << "Generating matrix in MatrixMarket format..." << std::endl;
        myHESolver->getInnerMatrix(mtx);

        std::ofstream outfile(tFilename.c_str());
        outfile << mtx;
        outfile.close();
        std::cout << "Generating matrix in MatrixMarket format... DONE! (" << myStopwatch->stop() << " s)" << std::endl << std::endl << std::endl;

        delete myHESolver;
      } else {
        throw new application_exception("HeatEquationSolver::storeInnerMatrix : A grid wasn't constructed before!");
      }
    }

    void HeatEquationSolver::storeInnerMatrixDiagonal(DataVector& alpha, std::string tFilename, double timestepsize) {
      if (this->bGridConstructed) {
        HeatEquationParabolicPDESolverSystem* myHESolver = new HeatEquationParabolicPDESolverSystem(*this->myGrid, alpha, this->a, timestepsize, "ImEul");
        SGppStopwatch* myStopwatch = new SGppStopwatch();

        std::string mtx = "";

        myStopwatch->start();
        std::cout << "Generating systemmatrix's diagonal in MatrixMarket format..." << std::endl;
        myHESolver->getInnerMatrixDiagonal(mtx);

        std::ofstream outfile(tFilename.c_str());
        outfile << mtx;
        outfile.close();
        std::cout << "Generating systemmatrix's diagonal in MatrixMarket format... DONE! (" << myStopwatch->stop() << " s)" << std::endl << std::endl << std::endl;

        delete myHESolver;
      } else {
        throw new application_exception("HeatEquationSolver::storeInnerMatrix : A grid wasn't constructed before!");
      }
    }

    void HeatEquationSolver::storeInnerMatrixDiagonalRowSum(DataVector& alpha, std::string tFilename, double timestepsize) {
      if (this->bGridConstructed) {
        HeatEquationParabolicPDESolverSystem* myHESolver = new HeatEquationParabolicPDESolverSystem(*this->myGrid, alpha, this->a, timestepsize, "ImEul");
        SGppStopwatch* myStopwatch = new SGppStopwatch();

        std::string mtx = "";

        myStopwatch->start();
        std::cout << "Generating systemmatrix rowsum as diagonal matrix in MatrixMarket format..." << std::endl;
        myHESolver->getInnerMatrixDiagonalRowSum(mtx);

        std::ofstream outfile(tFilename.c_str());
        outfile << mtx;
        outfile.close();
        std::cout << "Generating systemmatrix rowsum as diagonal matrix in MatrixMarket format... DONE! (" << myStopwatch->stop() << " s)" << std::endl << std::endl << std::endl;

        delete myHESolver;
      } else {
        throw new application_exception("HeatEquationSolver::storeInnerMatrix : A grid wasn't constructed before!");
      }
    }

    void HeatEquationSolver::storeInnerRHS(DataVector& alpha, std::string tFilename, double timestepsize) {
      if (this->bGridConstructed) {
        HeatEquationParabolicPDESolverSystem* myHESolver = new HeatEquationParabolicPDESolverSystem(*this->myGrid, alpha, this->a, timestepsize, "ImEul");
        SGppStopwatch* myStopwatch = new SGppStopwatch();

        myStopwatch->start();
        std::cout << "Exporting inner right-hand-side..." << std::endl;
        DataVector* rhs_inner = myHESolver->generateRHS();

        size_t nCoefs = rhs_inner->getSize();
        std::ofstream outfile(tFilename.c_str());

        for (size_t i = 0; i < nCoefs; i++) {
          outfile << std::scientific << rhs_inner->get(i) << std::endl;
        }

        outfile.close();
        std::cout << "Exporting inner right-hand-side... DONE! (" << myStopwatch->stop() << " s)" << std::endl << std::endl << std::endl;

        delete myHESolver;
      } else {
        throw new application_exception("HeatEquationSolver::storeInnerMatrix : A grid wasn't constructed before!");
      }
    }

    void HeatEquationSolver::storeInnerSolution(DataVector& alpha, size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, std::string tFilename) {
      if (this->bGridConstructed) {
        Euler* myEuler = new Euler("ImEul", numTimesteps, timestepsize, false, 0, this->myScreen);
        ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
        HeatEquationParabolicPDESolverSystem* myHESolver = new HeatEquationParabolicPDESolverSystem(*this->myGrid, alpha, this->a, timestepsize, "ImEul");
        SGppStopwatch* myStopwatch = new SGppStopwatch();

        myStopwatch->start();
        std::cout << "Exporting inner solution..." << std::endl;
        myEuler->solve(*myCG, *myHESolver, false);

        DataVector* alpha_solve = myHESolver->getGridCoefficientsForCG();
        size_t nCoefs = alpha_solve->getSize();
        std::ofstream outfile(tFilename.c_str());

        for (size_t i = 0; i < nCoefs; i++) {
          outfile << std::scientific << alpha_solve->get(i) << std::endl;
        }

        outfile.close();

        std::cout << "Exporting inner solution... DONE!" << std::endl;

        delete myHESolver;
        delete myCG;
        delete myEuler;
      } else {
        throw new application_exception("HeatEquationSolver::solveImplicitEuler : A grid wasn't constructed before!");
      }

    }


  }
}
