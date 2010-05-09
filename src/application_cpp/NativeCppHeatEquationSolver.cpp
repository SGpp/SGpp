/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#include "sgpp.hpp"

void testHeatEquation()
{
	size_t dim = 2;
	size_t level = 7;

	size_t timesteps = 100;
	double stepsize = 0.1;
	size_t CGiterations = 8000;
	double CGepsilon = 0.000001;

	double a = 1.0;

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];

	// set the bounding box
	for (size_t i = 0; i < dim; i++)
	{
		myBoundaries[i].leftBoundary = 0.0;
		myBoundaries[i].rightBoundary = 3.0;
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

	// init the basis functions' coefficient vector
	DataVector* alpha = new DataVector(myHESolver->getNumberGridPoints());

	//myHESolver->initGridWithSingleHeat(*alpha, 100.0);
	myHESolver->initGridWithSmoothHeat(*alpha, 3.0, 1.5, 5.0);
	//myHESolver->initGridWithConstantHeat(*alpha, 4.0);

	// Print the initial heat function into a gnuplot file
	myHESolver->printGrid(*alpha, 50, "heatStart.gnuplot");

	// set heat coefficient
	myHESolver->setHeatCoefficient(a);

	// Start solving the Heat Equation
	//myHESolver->solveExplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, true, true, 50);
	myHESolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, true, false, 50);
	//myHESolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha);

	// Print the solved Heat Equation into a gnuplot file
	myHESolver->printGrid(*alpha, 50, "solvedHeat.gnuplot");

	delete myHESolver;
	delete myBoundingBox;
	delete alpha;
}


int main(int argc, char *argv[])
{
	testHeatEquation();
}
