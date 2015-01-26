#include <sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp>
#include <sgpp/pde/application/HeatEquationSolverWithStretching.hpp>
#include <sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystemParallelOMP.hpp>
#include <sgpp/solver/ode/Euler.hpp>
#include <sgpp/solver/ode/CrankNicolson.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

using namespace sg::base;
using namespace sg::solver;

namespace sg {
  namespace pde {

    HeatEquationSolverWithStretching::HeatEquationSolverWithStretching() : ParabolicPDESolver() {
      this->bGridConstructed = false;
      this->myScreen = NULL;
    }

    HeatEquationSolverWithStretching::~HeatEquationSolverWithStretching() {
      if (this->myScreen != NULL) {
        delete this->myScreen;
      }
    }

    void HeatEquationSolverWithStretching::constructGrid(Stretching& stretching, int level) {
      this->dim = stretching.getDimensions();
      this->levels = level;

      this->myGrid = new LinearStretchedTrapezoidBoundaryGrid(stretching);

      GridGenerator* myGenerator = this->myGrid->createGridGenerator();
      myGenerator->regular(this->levels);
      delete myGenerator;

      this->myStretching = this->myGrid->getStretching();
      this->myGridStorage = this->myGrid->getStorage();

      this->bGridConstructed = true;
    }

    void HeatEquationSolverWithStretching::constructGrid(BoundingBox& BoundingBox, int level) {
      std::cout << "I'm not supposed to be here, me is constructGrid\n";
    }

    void HeatEquationSolverWithStretching::setHeatCoefficient(double a) {
      this->a = a;
    }

    void HeatEquationSolverWithStretching::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation) {
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
        throw new application_exception("HeatEquationSolverWithStretching::solveExplicitEuler : A grid wasn't constructed before!");
      }
    }

    void HeatEquationSolverWithStretching::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation) {
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
        throw new application_exception("HeatEquationSolverWithStretching::solveImplicitEuler : A grid wasn't constructed before!");
      }
    }

    void HeatEquationSolverWithStretching::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul) {
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
        throw new application_exception("HeatEquationSolverWithStretching::solveCrankNicolson : A grid wasn't constructed before!");
      }
    }

    void HeatEquationSolverWithStretching::initGridWithSmoothHeat(DataVector& alpha, double mu, double sigma, double factor) {
      if (this->bGridConstructed) {
        double tmp;
        double* dblFuncValues = new double[this->dim];

        for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
          std::string coords = this->myGridStorage->get(i)->getCoordsStringStretching(*this->myStretching);
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
        throw new application_exception("HeatEquationSolverWithStretching::initGridWithSmoothHeat : A grid wasn't constructed before!");
      }
    }

    void HeatEquationSolverWithStretching::initScreen() {
      this->myScreen = new ScreenOutput();
      this->myScreen->writeTitle("SGpp - Heat Equation Solver With Stretching, 1.0.1", "Alexander Heinecke, Sarpkan Selcuk (C) 2009-2011");
    }

    void HeatEquationSolverWithStretching::printGrid(DataVector& alpha, double PointesPerDimension, std::string tfilename) const {
      GridPrinterForStretching myPrinter(*this->myGrid);
      myPrinter.printGrid(alpha, tfilename, static_cast<size_t>(PointesPerDimension));
    }

    void HeatEquationSolverWithStretching::printGridDomain(DataVector& alpha, double PointesPerDimension, BoundingBox& GridArea, std::string tfilename) const {
      throw new application_exception("HeatEquationSolverWithStretching::printGridDomain : BoundingBox not supported with this solver, use printGridDomainStretching instead ");
    }

    void HeatEquationSolverWithStretching::printGridDomainStretching(DataVector& alpha, double PointesPerDimension, Stretching& GridArea, std::string tfilename) const {
      GridPrinterForStretching myPrinter(*this->myGrid);
      myPrinter.printGridDomainStretching(alpha, tfilename, GridArea, static_cast<size_t>(PointesPerDimension));
    }

    void HeatEquationSolverWithStretching::printSparseGrid(DataVector& alpha, std::string tfilename, bool bSurplus) const {
      GridPrinterForStretching myPrinter(*this->myGrid);
      myPrinter.printSparseGrid(alpha, tfilename, bSurplus);
    }

    void HeatEquationSolverWithStretching::printSparseGridExpTransform(DataVector& alpha, std::string tfilename, bool bSurplus) const {
      GridPrinterForStretching myPrinter(*this->myGrid);
      myPrinter.printSparseGridExpTransform(alpha, tfilename, bSurplus);
    }

  }
}
