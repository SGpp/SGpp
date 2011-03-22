/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#ifndef HEATEQUATIONSOLVERWITHSTRETCHING_HPP
#define HEATEQUATIONSOLVERWITHSTRETCHING_HPP

#include "sgpp.hpp"

#include "application/pde/ParabolicPDESolver.hpp"

//#include "grid/type/LinearTrapezoidBoundaryGrid.hpp"
//#include "grid/type/LinearGrid.hpp"
#include "grid/common/BoundingBox.hpp"

#include "grid/type/LinearStretchedTrapezoidBoundaryGrid.hpp"
#include "grid/type/LinearStretchedGrid.hpp"
#include "grid/common/Stretching.hpp"

#include "tools/common/StdNormalDistribution.hpp"

#include "application/common/ScreenOutput.hpp"
#include "tools/common/GridPrinterForStretching.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

namespace sg
{

/**
 * This class provides a simple-to-use solver of the multi dimensional
 * Heat Equation that uses Sparse Grids.
 *
 * The class's aim is, to hide all complex details of solving the
 * Heat Equation with Sparse Grids!
 *
 * @version $HEAD$
 */
class HeatEquationSolverWithStretching : public ParabolicPDESolver
{
private:
	/// the heat coefficient
	double a;
	/// screen object used in this solver
	ScreenOutput* myScreen;
	////BoundingBox replacement
	Stretching* myStretching;

public:
	/**
	 * Std-Constructor of the solver
	 */
	HeatEquationSolverWithStretching();

	/**
	 * Std-Destructor of the solver
	 */
	virtual ~HeatEquationSolverWithStretching();

	void constructGrid(Stretching& myStretching, size_t level);

	void constructGrid(BoundingBox& myStretching, size_t level);

	void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

	void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

	void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul = 0);

	/**
	 * This method sets the heat coefficient of the regarded material
	 *
	 * @param a the heat coefficient
	 */
	void setHeatCoefficient(double a);

	/**
	 * Inits the grid in the middle of the whole domain with one single heat
	 *
	 * @param alpha reference to the coefficients vector
	 * @param heat the value of the heat in the middle of the domain
	 */
//	void initGridWithSingleHeat(DataVector& alpha, double heat);

	/**
	 * Inits the grid in the middle the domain with an smooth heat distribution that the
	 * normal distribution formula
	 *
	 * @param alpha reference to the coefficients vector
	 * @param mu the exspected value of the normal distribution
	 * @param sigma the sigma of the normal distribution
	 * @param factor a factor that is used to stretch the function values
	 */
	void initGridWithSmoothHeat(DataVector& alpha, double mu, double sigma, double factor);

	/**
	 * Inits the grid with a constant heat
	 *
	 * @param alpha reference to the coefficients vector
	 * @param constHeat the temperature of the constant heat
	 */
//	void initGridWithConstantHeat(DataVector& alpha, double constHeat);

	/**
	 * Inits the screen object
	 */
	void initScreen();

	/**
	 * This is some kind of debug functionality. It writes a file,
	 * that can be used with gnuplot the print the grid.
	 *
	 * Is only implemented for 1D and 2D grids!
	 *
	 * @param alpha the coefficients of the Sparse Gird's basis functions
	 * @param PointesPerDimension the distance between evaluation points
	 * @param tfilename absolute path to file into which the grid's evaluation is written
	 */
	void printGrid(DataVector& alpha, double PointesPerDimension, std::string tfilename);

	/**
	 * This is some kind of debug functionality. It writes a file,
	 * that can be used with gnuplot the print the grid.
	 *
	 * Is only implemented for 2D grids!
	 *
	 * @param alpha the coefficients of the Sparse Gird's basis functions
	 * @param PointesPerDimension the distance between evaluation points
	 * @param GridArea the area in which the function should be plotted
	 * @param tfilename absolute path to file into which the grid's evaluation is written
	 */
	void printGridDomain(DataVector& alpha, double PointesPerDimension, Stretching& GridArea, std::string tfilename);

	/**
	 * Prints the Grid Points of the Sparse Grid either with their node basis value
	 * or their hierarchical surplus
	 *
	 * This function is available for all dimensions
	 *
	 * @param alpha the coefficients of the grid's ansatzfunctions
	 * @param tfilename absoulte path to the file the grid is written into
	 * @param bSurplus specifies whether the surplus (true) or the node basis value (false) is written
	 */
	void printSparseGrid(DataVector& alpha, std::string tfilename, bool bSurplus);

	/**
	 * Prints the Grid Points of the Sparse Grid either with their node basis value
	 * or their hierarchical surplus
	 *
	 * This function is available for all dimensions.
	 *
	 * The coordinates of the grid points are pushed the exp function. So
	 * log transformed grids can be plotted in cartesion coordinates.
	 *
	 * @param alpha the coefficients of the grid's ansatzfunctions
	 * @param tfilename absoulte path to the file the grid is written into
	 * @param bSurplus specifies whether the surplus (true) or the node basis value (false) is written
	 */
	void printSparseGridExpTransform(DataVector& alpha, std::string tfilename, bool bSurplus);

};

}

#endif /* HEATEQUATIONSOLVERWITHSTRETCHING_HPP */
