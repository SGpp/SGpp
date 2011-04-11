/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

/// default number of Implicit Euler steps when using Crank Nicolson
#define CRNIC_IMEUL_STEPS 3
#define GUNPLOT_RESOLUTION 51
#define SOLUTION_FRAMES 100

#define DIV_SIGMA 4.0
#define DISTRI_FACTOR 5.0

//#define EXPORT_MATRIX_FILES

#include <cstdlib>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <iomanip>

#include "sgpp.hpp"

/**
 * Writes a DataMatrix into a file
 *
 * @param data the DataMatrix that should be written into a file
 * @param tFile the file into which the data is written
 *
 * @return error code
 */
int writeDataMatrix(DataMatrix& data, std::string tFile)
{
	std::ofstream file;
	file.open(tFile.c_str());

	if(!file.is_open())
	{
		std::cout << "Error cannot write file: " << tFile << std::endl;
		return -1;
	}

	for (size_t i = 0; i < data.getNrows(); i++)
	{
		for (size_t j = 0; j < data.getNcols(); j++)
		{
			file << std::scientific << std::setprecision( 16 ) << data.get(i,j) << " ";
		}
		file << std::endl;
	}

	file.close();

	return 0;
}


/**
 * Writes a DataVector into a file
 *
 * @param data the DataVector that should be written into a file
 * @param tFile the file into which the data is written
 *
 * @return error code
 */
int writeDataVector(DataVector& data, std::string tFile)
{
	std::ofstream file;
	file.open(tFile.c_str());

	if(!file.is_open())
	{
		std::cout << "Error cannot write file: " << tFile << std::endl;
		return -1;
	}

	for (size_t i = 0; i < data.getSize(); i++)
	{

		file << std::scientific << std::setprecision( 16 ) << data.get(i) << " " << std::endl;
	}

	file.close();

	return 0;
}

/**
 * Calls the writeHelp method in the BlackScholesSolver Object
 * after creating a screen.
 */
void writeHelp()
{
	sg::HeatEquationSolver* myHESolver = new sg::HeatEquationSolver();

	myHESolver->initScreen();

	delete myHESolver;

	std::stringstream mySStream;

	mySStream << "Some instructions for the use of Poisson/ Heat Equation Solver:" << std::endl;
	mySStream << "---------------------------------------------------------------" << std::endl << std::endl;
	mySStream << "Available execution modes are:" << std::endl;
	mySStream << "  HeatEquation        Solves Heat Equation on a quadratic" << std::endl;
	mySStream << "                      d-dimensional domain" << std::endl << std::endl;
	mySStream << "  PoissonEquation     Solves Poisson Equation on a quadratic" << std::endl;
	mySStream << "                      d-dimensional domain" << std::endl << std::endl << std::endl;

	mySStream << "Execution modes descriptions:" << std::endl;
	mySStream << "-----------------------------------------------------" << std::endl;
	mySStream << "HeatEquation" << std::endl << "------" << std::endl;
	mySStream << "the following options must be specified:" << std::endl;
	mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
	mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
	mySStream << "	left_bound: x_i of left boundary" << std::endl;
	mySStream << "	right_bound: x_i of right boundary" << std::endl;
	mySStream << "	a: thermal diffusivity" << std::endl;
	mySStream << "	initHeat: initial heat distribution" << std::endl;
	mySStream << "	T: time to solve" << std::endl;
	mySStream << "	dT: timestep size" << std::endl;
	mySStream << "	Solver: the solver to use: ExEul, ImEul, CrNic" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "HESolver HeatEquation 3 5 0.0 3.0 1.0 smooth 1.0 0.1 ImEul 0.00001 400" << std::endl;
	mySStream << std::endl;
	mySStream << "Remark: This test generates following files (gnuplot):" << std::endl;
	mySStream << "	heatStart.gnuplot: the start condition" << std::endl;
	mySStream << "	heatSolved.gnuplot: the numerical solution" << std::endl;
	mySStream << std::endl << std::endl;

	mySStream << "PoissonEquation" << std::endl << "------" << std::endl;
	mySStream << "the following options must be specified:" << std::endl;
	mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
	mySStream << "	start_level: number of start-levels of the Sparse Grid" << std::endl;
	mySStream << "	end_level: number of max. levels of the Sparse Grid" << std::endl;
	mySStream << "	left_bound: x_i of left boundary" << std::endl;
	mySStream << "	right_bound: x_i of right boundary" << std::endl;
	mySStream << "	initHeat: initial heat distribution" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "HESolver PoissonEquation 3 3 5 0.0 3.0 smooth 0.00001 400" << std::endl;
	mySStream << std::endl;
	mySStream << "Remark: This test generates following files (gnuplot):" << std::endl;
	mySStream << "	poissonStart.gnuplot: the start condition" << std::endl;
	mySStream << "	poissonSolved.gnuplot: the numerical solution" << std::endl;
	mySStream << std::endl << std::endl;

	std::cout << mySStream.str() << std::endl;
}

void testHeatEquation(size_t dim, size_t level, double bound_left, double bound_right, double a,
						std::string initFunc, double T, double dt, std::string ODESolver,
						double cg_eps, size_t cg_its)
{
	size_t timesteps = (size_t)(T/dt);

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];

	// set the bounding box
	for (size_t i = 0; i < dim; i++)
	{
		myBoundaries[i].leftBoundary = bound_left;
		myBoundaries[i].rightBoundary = bound_right;
		myBoundaries[i].bDirichletLeft = true;
		myBoundaries[i].bDirichletRight = true;
	}

	sg::HeatEquationSolver* myHESolver = new sg::HeatEquationSolver();
	sg::BoundingBox* myBoundingBox = new sg::BoundingBox(dim, myBoundaries);
	delete[] myBoundaries;

	// init Screen Object
	myHESolver->initScreen();

	// Construct a grid
	myHESolver->constructGrid(*myBoundingBox, level);

	// init the basis functions' coefficient vector (start solution)
	DataVector* alpha = new DataVector(myHESolver->getNumberGridPoints());
	if (initFunc == "smooth")
	{
		myHESolver->initGridWithSmoothHeat(*alpha, bound_right, bound_right/DIV_SIGMA, DISTRI_FACTOR);
	}
	else
	{
		writeHelp();
	}

	// Print the initial heat function into a gnuplot file
	if (dim < 3)
	{
		myHESolver->printGrid(*alpha, GUNPLOT_RESOLUTION, "heatStart.gnuplot");
	}

	// set heat coefficient
	myHESolver->setHeatCoefficient(a);

	// Start solving the Heat Equation
	if (ODESolver == "ExEul")
	{
		myHESolver->solveExplicitEuler(timesteps, dt, cg_its, cg_eps, *alpha, true, false, std::max(timesteps/SOLUTION_FRAMES,(size_t)1));
	}
	else if (ODESolver == "ImEul")
	{
		myHESolver->solveImplicitEuler(timesteps, dt, cg_its, cg_eps, *alpha, true, false, std::max(timesteps/SOLUTION_FRAMES,(size_t)1));
	}
	else if (ODESolver == "CrNic")
	{
		myHESolver->solveCrankNicolson(timesteps, dt, cg_its, cg_eps, *alpha, CRNIC_IMEUL_STEPS);
	}

	// Print the solved Heat Equation into a gnuplot file
	if (dim < 3)
	{
		myHESolver->printGrid(*alpha, GUNPLOT_RESOLUTION, "heatSolved.gnuplot");
	}

#ifdef EXPORT_MATRIX_FILES
	// print inner matrix
	std::stringstream mtxFile;
	mtxFile.clear();
	alpha->setAll(0.0);
	if (initFunc == "smooth")
	{
		myHESolver->initGridWithSmoothHeat(*alpha, bound_right, bound_right/DIV_SIGMA, DISTRI_FACTOR);
	}
	else
	{
		writeHelp();
	}
	mtxFile << "SG_HeatEquation_InnerMatrix_" << dim << "d_" << level << "l.mtx";
	myHESolver->storeInnerMatrix(*alpha, mtxFile.str(), dt);

	// print inner matrix, diagonal
	std::stringstream mtxFileDiagonal;
	mtxFileDiagonal.clear();
	alpha->setAll(0.0);
	if (initFunc == "smooth")
	{
		myHESolver->initGridWithSmoothHeat(*alpha, bound_right, bound_right/DIV_SIGMA, DISTRI_FACTOR);
	}
	else
	{
		writeHelp();
	}
	mtxFileDiagonal << "SG_HeatEquation_InnerMatrixDiagonal_" << dim << "d_" << level << "l.mtx";
	myHESolver->storeInnerMatrixDiagonal(*alpha, mtxFileDiagonal.str(), dt);

	// print inner matrix, diagonal row sum
	std::stringstream mtxFileDiagonalRowSum;
	mtxFileDiagonalRowSum.clear();
	alpha->setAll(0.0);
	if (initFunc == "smooth")
	{
		myHESolver->initGridWithSmoothHeat(*alpha, bound_right, bound_right/DIV_SIGMA, DISTRI_FACTOR);
	}
	else
	{
		writeHelp();
	}
	mtxFileDiagonalRowSum << "SG_HeatEquation_InnerMatrixDiagonalRowSum_" << dim << "d_" << level << "l.mtx";
	myHESolver->storeInnerMatrixDiagonalRowSum(*alpha, mtxFileDiagonalRowSum.str(), dt);

	// print inner rhs
	std::stringstream rhsFile;
	rhsFile.clear();
	alpha->setAll(0.0);
	if (initFunc == "smooth")
	{
		myHESolver->initGridWithSmoothHeat(*alpha, bound_right, bound_right/DIV_SIGMA, DISTRI_FACTOR);
	}
	else
	{
		writeHelp();
	}
	rhsFile << "SG_HeatEquation_InnerRHS_" << dim << "d_" << level << "l.vec";
	myHESolver->storeInnerRHS(*alpha, rhsFile.str(), dt);

	// print inner solution
	std::stringstream solFile;
	solFile.clear();
	alpha->setAll(0.0);
	if (initFunc == "smooth")
	{
		myHESolver->initGridWithSmoothHeat(*alpha, bound_right, bound_right/DIV_SIGMA, DISTRI_FACTOR);
	}
	else
	{
		writeHelp();
	}
	solFile << "SG_HeatEquation_InnerSolution_" << dim << "d_" << level << "l.vec";
	myHESolver->storeInnerSolution(*alpha, timesteps, dt, cg_its, cg_eps, solFile.str());
	std::cout << std::endl << std::endl;
#endif

	delete myHESolver;
	delete myBoundingBox;
	delete alpha;
}

void testPoissonEquation(size_t dim, size_t start_level, size_t end_level, double bound_left, double bound_right,
						std::string initFunc, double cg_eps, size_t cg_its)
{
	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];
	DataMatrix EvalPoints(1, dim);
	std::string tFileEvalCuboid = "EvalPointsPoisson";
	std::string tFileEvalCuboidValues = "EvalValuesPoisson";
	size_t evalPoints = 21;
	std::vector<DataVector> results;

	// set the bounding box
	for (size_t i = 0; i < dim; i++)
	{
		myBoundaries[i].leftBoundary = bound_left;
		myBoundaries[i].rightBoundary = bound_right;
		myBoundaries[i].bDirichletLeft = true;
		myBoundaries[i].bDirichletRight = true;
	}

	sg::PoissonEquationSolver* myPoisSolver = new sg::PoissonEquationSolver();
	sg::BoundingBox* myBoundingBox = new sg::BoundingBox(dim, myBoundaries);
	delete[] myBoundaries;

	sg::EvalCuboidGenerator* myEvalCuboidGen = new sg::EvalCuboidGenerator();

	// init Screen Object
	myPoisSolver->initScreen();

	for (size_t l = start_level; l <= end_level; l++)
	{
		// Construct a grid
		myPoisSolver->constructGrid(*myBoundingBox, l);

		// in first iteration -> calculate the evaluation points
		if (l == start_level)
		{
			myEvalCuboidGen->getEvaluationCuboid(EvalPoints, *myBoundingBox, evalPoints);

			writeDataMatrix(EvalPoints, tFileEvalCuboid);
		}

		// init the basis functions' coefficient vector (start solution)
		DataVector* alpha = new DataVector(myPoisSolver->getNumberGridPoints());
		if (initFunc == "smooth")
		{
			myPoisSolver->initGridWithSmoothHeat(*alpha, bound_right, bound_right/DIV_SIGMA, DISTRI_FACTOR);
		}
		else
		{
			writeHelp();
		}

		// Print the initial heat function into a gnuplot file
		if (dim < 3)
		{
			myPoisSolver->printGrid(*alpha, GUNPLOT_RESOLUTION, "poissonStart.gnuplot");
		}

		// solve Poisson Equation
		myPoisSolver->solvePDE(*alpha, *alpha, cg_its, cg_eps, true);

		// Print the solved Heat Equation into a gnuplot file
		if (dim < 3)
		{
			myPoisSolver->printGrid(*alpha, GUNPLOT_RESOLUTION, "poissonSolved.gnuplot");
		}

		// Calculate Norms
		// Evaluate Cuboid
		DataVector PoisEvals(EvalPoints.getNrows());
		myPoisSolver->evaluateCuboid(*alpha, PoisEvals, EvalPoints);
		results.push_back(PoisEvals);

		// write solution in a additional file
		std::stringstream level_string;
		level_string << l;
		writeDataVector(PoisEvals, tFileEvalCuboidValues+".level_"+ level_string.str());
		writeDataVector(PoisEvals, tFileEvalCuboidValues);

		if (l > start_level)
		{
			std::cout << "=====================================================================" << std::endl;
			std::cout << "=====================================================================" << std::endl << std::endl;
			std::cout << "Calculating norms of relative errors to a grid" << std::endl;
			std::cout << "with " << l << " levels and testing-coboid" << std::endl;
			std::cout << "with the bounding box:" << std::endl;
			for (size_t j = 0; j < dim; j++)
			{
				std::cout << myBoundingBox->getBoundary(j).leftBoundary << " " << myBoundingBox->getBoundary(j).rightBoundary << std::endl;
			}
			std::cout << std::endl << std::endl;

			double oldMaxNorm = 0.0;
			double oldTwoNorm = 0.0;

			// Calculate relative errors and some norms
			for (size_t j = 0; j < l-start_level; j++)
			{
				DataVector maxLevel(results[l-start_level]);
				DataVector relError(results[j]);
				double maxNorm = 0.0;
				double l2Norm = 0.0;

				// calculate relative error
				relError.sub(maxLevel);
				relError.componentwise_div(maxLevel);

				// calculate max. norm of relative error
				maxNorm = relError.maxNorm();

				// calculate two norm of relative error
				l2Norm = relError.RMSNorm();

				// Printing norms
				std::cout << "Level " << j + start_level << ": max-norm(rel-error)=" << maxNorm << "; two-norm(rel-error)=" << l2Norm << "; rate max-norm: " << log(oldMaxNorm/maxNorm) << "; rate two-norm: " << log(oldTwoNorm/l2Norm) << std::endl;

				oldMaxNorm = maxNorm;
				oldTwoNorm = l2Norm;
			}
		}
		std::cout << std::endl << std::endl;


	#ifdef EXPORT_MATRIX_FILES
		// print inner matrix
		std::stringstream mtxFile;
		mtxFile.clear();
		mtxFile << "SG_Poisson_InnerMatrix_" << dim << "d_" << level << "l.mtx";
		myPoisSolver->storeInnerMatrix(mtxFile.str());

		// print inner rhs
		std::stringstream rhsFile;
		rhsFile.clear();
		alpha->setAll(0.0);
		if (initFunc == "smooth")
		{
			myPoisSolver->initGridWithSmoothHeat(*alpha, bound_right, bound_right/DIV_SIGMA, DISTRI_FACTOR);
		}
		else
		{
			writeHelp();
		}
		rhsFile << "SG_Poisson_InnerRHS_" << dim << "d_" << level << "l.vec";
		myPoisSolver->storeInnerRHS(*alpha, rhsFile.str());

		// print inner matrix, diagonal only
		std::stringstream mtxDiagFile;
		mtxDiagFile.clear();
		mtxDiagFile << "SG_Poisson_InnerMatrixDiagonal_" << dim << "d_" << level << "l.mtx";
		myPoisSolver->storeInnerMatrixDiagonal(mtxDiagFile.str());

		// print inner matrix, diagonal containing row sum
		std::stringstream mtxDiagRowSumFile;
		mtxDiagRowSumFile.clear();
		mtxDiagRowSumFile << "SG_Poisson_InnerMatrixDiagonalRowSum_" << dim << "d_" << level << "l.mtx";
		myPoisSolver->storeInnerMatrixDiagonalRowSum(mtxDiagRowSumFile.str());

		// print inner solution
		std::stringstream solFile;
		solFile.clear();
		alpha->setAll(0.0);
		if (initFunc == "smooth")
		{
			myPoisSolver->initGridWithSmoothHeat(*alpha, bound_right, bound_right/DIV_SIGMA, DISTRI_FACTOR);
		}
		else
		{
			writeHelp();
		}
		solFile << "SG_Poisson_InnerSolution_" << dim << "d_" << level << "l.vec";
		myPoisSolver->storeInnerSolution(*alpha, cg_its, cg_eps, solFile.str());
		std::cout << std::endl << std::endl;
	#endif

		// Iteration cleanup
		myPoisSolver->deleteGrid();
		delete alpha;
	}

	delete myPoisSolver;
	delete myBoundingBox;
}

int main(int argc, char *argv[])
{
	std::string option;

	if (argc == 1)
	{
		writeHelp();
		return 0;
	}

	option.assign(argv[1]);

	if (option == "HeatEquation")
	{
		if (argc != 13)
		{
			writeHelp();
			return 0;
		}

		size_t dim;
		size_t level;
		double bound_left;
		double bound_right;
		double a;
		std::string initFunc;
		double T;
		double dt;
		std::string ODESolver;
		double cg_eps;
		size_t cg_its;

		dim = atoi(argv[2]);
		level = atoi(argv[3]);
		bound_left = atof(argv[4]);
		bound_right = atof(argv[5]);
		a = atof(argv[6]);
		initFunc.assign(argv[7]);
		T = atof(argv[8]);
		dt = atof(argv[9]);
		ODESolver.assign(argv[10]);
		cg_eps = atof(argv[11]);
		cg_its = atoi(argv[12]);

		testHeatEquation(dim, level, bound_left, bound_right, a, initFunc, T, dt, ODESolver, cg_eps, cg_its);
	}
	else if (option == "PoissonEquation")
	{
		if (argc != 10)
		{
			writeHelp();
			return 0;
		}

		size_t dim;
		size_t start_level;
		size_t end_level;

		double bound_left;
		double bound_right;
		std::string initFunc;
		double cg_eps;
		size_t cg_its;

		dim = atoi(argv[2]);
		start_level = atoi(argv[3]);
		end_level = atoi(argv[4]);
		bound_left = atof(argv[5]);
		bound_right = atof(argv[6]);
		initFunc.assign(argv[7]);
		cg_eps = atof(argv[8]);
		cg_its = atoi(argv[9]);
		testPoissonEquation(dim, start_level, end_level, bound_left, bound_right, initFunc, cg_eps, cg_its);
	}
	else
	{
		writeHelp();
	}
}
