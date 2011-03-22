/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

/*#include "algorithm/pde/HeatEquationODESolverSystem.hpp"
#include "application/pde/HeatEquationSolverWithStretching.hpp"
#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "solver/sle/ConjugateGradients.hpp"
#include "grid/Grid.hpp"
#include "exception/application_exception.hpp"
#include "stdlib.h"
#include <sstream>

namespace sg
{

HeatEquationSolverWithStretching::HeatEquationSolverWithStretching() : ParabolicPDESolver()
{
	this->bGridConstructed = false;
	this->myScreen = NULL;
}

HeatEquationSolverWithStretching::~HeatEquationSolverWithStretching()
{
	if (this->bGridConstructed)
	{
		delete this->myGrid;
	}
	if (this->myScreen != NULL)
	{
		delete this->myScreen;
	}
}

void HeatEquationSolverWithStretching::constructGrid(Stretching& stretching, size_t level)
{
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

void HeatEquationSolverWithStretching::constructGrid(BoundingBox& BoundingBox, size_t level){
	std::cout<<"I'm not supposed to be here, me is constructGrid\n";
}

void HeatEquationSolverWithStretching::setHeatCoefficient(double a)
{
	this->a = a;
}

void HeatEquationSolverWithStretching::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	if (this->bGridConstructed)
	{
		this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
		Euler* myEuler = new Euler("ExEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, this->myScreen);
		ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
		HeatEquationODESolverSystem* myHESolver = new HeatEquationODESolverSystem(*this->myGrid, alpha, this->a, timestepsize, "ExEul");

		myEuler->solve(*myCG, *myHESolver, verbose);

		delete myHESolver;
		delete myCG;
		delete myEuler;
	}
	else
	{
		throw new application_exception("HeatEquationSolverWithStretching::solveExplicitEuler : A grid wasn't constructed before!");
	}
}

void HeatEquationSolverWithStretching::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	if (this->bGridConstructed)
	{
		this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
		Euler* myEuler = new Euler("ImEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, this->myScreen);
		ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
		HeatEquationODESolverSystem* myHESolver = new HeatEquationODESolverSystem(*this->myGrid, alpha, this->a, timestepsize, "ImEul");

		myEuler->solve(*myCG, *myHESolver, verbose);

		delete myHESolver;
		delete myCG;
		delete myEuler;
	}
	else
	{
		throw new application_exception("HeatEquationSolverWithStretching::solveImplicitEuler : A grid wasn't constructed before!");
	}
}

void HeatEquationSolverWithStretching::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul)
{
	if (this->bGridConstructed)
	{
		this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
		ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
		HeatEquationODESolverSystem* myHESolver = new HeatEquationODESolverSystem(*(this->myGrid), alpha, this->a, timestepsize, "CrNic");
		std::cout<<"Solver and Solver system generated @CrankNicolson"<<std::endl;
		size_t numCNSteps;
		size_t numIESteps;

		numCNSteps = numTimesteps;
		if (numTimesteps > NumImEul)
		{
			numCNSteps = numTimesteps - NumImEul;
		}
		numIESteps = NumImEul;

		Euler* myEuler = new Euler("ImEul", numIESteps, timestepsize, false, 0, this->myScreen);
		CrankNicolson* myCN = new CrankNicolson(numCNSteps, timestepsize);

		if (numIESteps > 0)
		{
//			myEuler->solve(*myCG, *myHESolver, false);
			myEuler->solve(*myCG, *myHESolver, true);
		}
		std::cout<<"Solve operation starts @CrankNicolson"<<std::endl;
		myCN->solve(*myCG, *myHESolver, false);
//		myCN->solve(*myCG, *myHESolver, true);
//		std::cout<<"Solve operation ends @CrankNicolson"<<std::endl;

//		if(myHESolver != NULL)
//		delete myHESolver;
		delete myCG;
//		delete myCN;
//		delete myEuler;
	}
	else
	{
		throw new application_exception("HeatEquationSolverWithStretching::solveCrankNicolson : A grid wasn't constructed before!");
	}
}

void HeatEquationSolverWithStretching::initGridWithSingleHeat(DataVector& alpha, double heat)
{
	double tmp;
	double tmp2;

	if (this->bGridConstructed)
	{
		if (this->dim == 1)
		{
			for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
			{
				tmp = atof(this->myGridStorage->get(i)->getCoordsStringStretching(*this->myStretching).c_str());

				if (tmp == 0.5)
				{
					alpha[i] = heat;
				}
				else
				{
					alpha[i] = 0.0;
				}
			}

			OperationHierarchisation* myHierarchisation = this->myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
			delete myHierarchisation;
		}
		else if (this->dim == 2)
		{
			for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
			{
					std::string coords = this->myGridStorage->get(i)->getCoordsStringStretching(*this->myStretching);
					std::stringstream coordsStream(coords);

					coordsStream >> tmp;
					coordsStream >> tmp2;

					if (tmp == 0.5 && tmp2 == 0.5)
					{
						alpha[i] = heat;
					}
					else
					{
						alpha[i] = 0.0;
					}
			}

			OperationHierarchisation* myHierarchisation = this->myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
			delete myHierarchisation;
		}
		else
		{
			throw new application_exception("HeatEquationSolverWithStretching::initGridWithSingleHeat : The constructed grid has more than two dimensions!");
		}
	}
	else
	{
		throw new application_exception("HeatEquationSolverWithStretching::initGridWithSingleHeat : A grid wasn't constructed before!");
	}
}

void HeatEquationSolverWithStretching::initGridWithSmoothHeat(DataVector& alpha, double mu, double sigma, double factor)
{
	double tmp;
	double tmp2;

	if (this->bGridConstructed)
	{
		if (this->dim == 1)
		{
			for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
			{
				tmp = atof(this->myGridStorage->get(i)->getCoordsStringStretching(*this->myStretching).c_str());

				alpha[i] = factor*(1.0/(sigma*2.0*3.145))*exp((-0.5)*((tmp-mu)/sigma)*((tmp-mu)/sigma));
			}

			OperationHierarchisation* myHierarchisation = this->myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
			delete myHierarchisation;
		}
		else if (this->dim == 2)
		{
			for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
			{
				std::string coords = this->myGridStorage->get(i)->getCoordsStringStretching(*this->myStretching);
				std::stringstream coordsStream(coords);

				coordsStream >> tmp;
				coordsStream >> tmp2;

				alpha[i] = factor*factor*((1.0/(sigma*2.0*3.145))*exp((-0.5)*((tmp-mu)/sigma)*((tmp-mu)/sigma))) * ((1.0/(sigma*2.0*3.145))*exp((-0.5)*((tmp2-mu)/sigma)*((tmp2-mu)/sigma)));
			}

			OperationHierarchisation* myHierarchisation = this->myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
			delete myHierarchisation;
		}
		else
		{
			throw new application_exception("HeatEquationSolverWithStretching::initGridWithSmoothHeat : The constructed grid has more than two dimensions!");
		}
	}
	else
	{
		throw new application_exception("HeatEquationSolverWithStretching::initGridWithSmoothHeat : A grid wasn't constructed before!");
	}
}

void HeatEquationSolverWithStretching::initGridWithConstantHeat(DataVector& alpha, double constHeat)
{
	double tmp;
	//double tmp2;

	if (this->bGridConstructed)
	{
		if (this->dim == 1)
		{
			for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
			{
				tmp = atof(this->myGridStorage->get(i)->getCoordsStringStretching(*this->myStretching).c_str());

				alpha[i] = constHeat;
			}

			OperationHierarchisation* myHierarchisation = this->myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
			delete myHierarchisation;
		}
		else
		{
			throw new application_exception("HeatEquationSolverWithStretching::initGridWithConstantHeat : The constructed grid has more than one dimension!");
		}
	}
	else
	{
		throw new application_exception("HeatEquationSolverWithStretching::initGridWithConstantHeat : A grid wasn't constructed before!");
	}
}

void HeatEquationSolverWithStretching::initScreen()
{
	this->myScreen = new ScreenOutput();
	this->myScreen->writeTitle("SGpp - Heat Equation Solver with Stretching, 1.0.0", "Alexander Heinecke, Sarpkan Selcuk, (C) 2009-2011");
}*/

#include "algorithm/pde/HeatEquationParabolicPDESolverSystem.hpp"
#include "algorithm/pde/HeatEquationParabolicPDESolverSystemParallelOMP.hpp"
#include "application/pde/HeatEquationSolverWithStretching.hpp"
#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "solver/sle/ConjugateGradients.hpp"
#include "grid/Grid.hpp"
#include "exception/application_exception.hpp"
#include "stdlib.h"
#include <sstream>

namespace sg
{

HeatEquationSolverWithStretching::HeatEquationSolverWithStretching() : ParabolicPDESolver()
{
	this->bGridConstructed = false;
	this->myScreen = NULL;
}

HeatEquationSolverWithStretching::~HeatEquationSolverWithStretching()
{
	if (this->myScreen != NULL)
	{
		delete this->myScreen;
	}
}

void HeatEquationSolverWithStretching::constructGrid(Stretching& stretching, size_t level)
{
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

void HeatEquationSolverWithStretching::constructGrid(BoundingBox& BoundingBox, size_t level){
	std::cout<<"I'm not supposed to be here, me is constructGrid\n";
}

void HeatEquationSolverWithStretching::setHeatCoefficient(double a)
{
	this->a = a;
}

void HeatEquationSolverWithStretching::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	if (this->bGridConstructed)
	{
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

		if (this->myScreen != NULL)
		{
			std::cout << "Time to solve: " << dNeededTime << " seconds" << std::endl;
			this->myScreen->writeEmptyLines(2);
		}

		delete myStopwatch;
		delete myCG;
		delete myEuler;
	}
	else
	{
		throw new application_exception("HeatEquationSolverWithStretching::solveExplicitEuler : A grid wasn't constructed before!");
	}
}

void HeatEquationSolverWithStretching::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	if (this->bGridConstructed)
	{
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

		if (this->myScreen != NULL)
		{
			std::cout << "Time to solve: " << dNeededTime << " seconds" << std::endl;
			this->myScreen->writeEmptyLines(2);
		}

		delete myStopwatch;
		delete myHESolver;
		delete myCG;
		delete myEuler;
	}
	else
	{
		throw new application_exception("HeatEquationSolverWithStretching::solveImplicitEuler : A grid wasn't constructed before!");
	}
}

void HeatEquationSolverWithStretching::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul)
{
	if (this->bGridConstructed)
	{
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
		if (numTimesteps > NumImEul)
		{
			numCNSteps = numTimesteps - NumImEul;
		}
		numIESteps = NumImEul;

		Euler* myEuler = new Euler("ImEul", numIESteps, timestepsize, false, 0, this->myScreen);
		CrankNicolson* myCN = new CrankNicolson(numCNSteps, timestepsize);

		myStopwatch->start();
		if (numIESteps > 0)
		{
			myEuler->solve(*myCG, *myHESolver, false);
		}

		myCN->solve(*myCG, *myHESolver, false);
		dNeededTime = myStopwatch->stop();

		if (this->myScreen != NULL)
		{
			std::cout << "Time to solve: " << dNeededTime << " seconds" << std::endl;
			this->myScreen->writeEmptyLines(2);
		}

		delete myStopwatch;
		delete myHESolver;
		delete myCG;
		delete myCN;
		delete myEuler;
	}
	else
	{
		throw new application_exception("HeatEquationSolverWithStretching::solveCrankNicolson : A grid wasn't constructed before!");
	}
}

void HeatEquationSolverWithStretching::initGridWithSmoothHeat(DataVector& alpha, double mu, double sigma, double factor)
{
	if (this->bGridConstructed)
	{
		double tmp;
		double* dblFuncValues = new double[this->dim];

		for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
		{
			std::string coords = this->myGridStorage->get(i)->getCoordsStringStretching(*this->myStretching);
			std::stringstream coordsStream(coords);

			for (size_t j = 0; j < this->dim; j++)
			{
				coordsStream >> tmp;

				dblFuncValues[j] = tmp;
			}

			tmp = 1.0;
			for (size_t j = 0; j < this->dim; j++)
			{
				tmp *=  factor*factor*((1.0/(sigma*2.0*3.145))*exp((-0.5)*((dblFuncValues[j]-mu)/sigma)*((dblFuncValues[j]-mu)/sigma)));
			}

			alpha[i] = tmp;
		}

		delete[] dblFuncValues;

		OperationHierarchisation* myHierarchisation = this->myGrid->createOperationHierarchisation();
		myHierarchisation->doHierarchisation(alpha);
		delete myHierarchisation;
	}
	else
	{
		throw new application_exception("HeatEquationSolverWithStretching::initGridWithSmoothHeat : A grid wasn't constructed before!");
	}
}

void HeatEquationSolverWithStretching::initScreen()
{
	this->myScreen = new ScreenOutput();
	this->myScreen->writeTitle("SGpp - Heat Equation Solver, 1.0.0", "Alexander Heinecke, (C) 2009-2011");
}

void HeatEquationSolverWithStretching::printGrid(DataVector& alpha, double PointesPerDimension, std::string tfilename)
{
	GridPrinterForStretching myPrinter(*this->myGrid);
	myPrinter.printGrid(alpha, tfilename, PointesPerDimension);
}

void HeatEquationSolverWithStretching::printGridDomain(DataVector& alpha, double PointesPerDimension, Stretching& GridArea, std::string tfilename)
{
	GridPrinterForStretching myPrinter(*this->myGrid);
	myPrinter.printGridDomain(alpha, tfilename, GridArea, PointesPerDimension);
}

void HeatEquationSolverWithStretching::printSparseGrid(DataVector& alpha, std::string tfilename, bool bSurplus)
{
	GridPrinterForStretching myPrinter(*this->myGrid);
	myPrinter.printSparseGrid(alpha, tfilename, bSurplus);
}

void HeatEquationSolverWithStretching::printSparseGridExpTransform(DataVector& alpha, std::string tfilename, bool bSurplus)
{
	GridPrinterForStretching myPrinter(*this->myGrid);
	myPrinter.printSparseGridExpTransform(alpha, tfilename, bSurplus);
}

}
