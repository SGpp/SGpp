/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
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
	size_t dim = 1;
	size_t level = 8;

	size_t timesteps = 10000;
	double stepsize = 0.0001;
	size_t CGiterations = 200;
	double CGepsilon = 0.00000001;

	double a = 0.5;

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];

	// set the bounding box
	for (size_t i = 0; i < dim; i++)
	{
		myBoundaries[i].leftBoundary = 0.0;
		myBoundaries[i].rightBoundary = 1.0;
		myBoundaries[i].bDirichletLeft = true;
		myBoundaries[i].bDirichletRight = true;
	}

	sg::HeatEquationSolver* myHESolver = new sg::HeatEquationSolver();
	sg::BoundingBox* myBoundingBox = new sg::BoundingBox(dim, myBoundaries);
	delete[] myBoundaries;

	// Construct a grid
	myHESolver->constructGrid(*myBoundingBox, level);

	// init the basis functions' coefficient vector
	DataVector* alpha = new DataVector(myHESolver->getNumberGridPoints());

	//myHESolver->initGridWithSingleHeat(*alpha, 100.0);
	myHESolver->initGridWithSmoothHeat(*alpha, 0.5, 0.08);
	//myHESolver->initGridWithConstantHeat(*alpha, 4.0);

	// Print the payoff function into a gnuplot file
	myHESolver->printGrid(*alpha, 1000, "heatStart.gnuplot");

	// Start solving the Black Scholes Equation
	myHESolver->solveExplicitEuler(timesteps, stepsize, a, *alpha, true);
	//myHESolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, a, *alpha, true);
	//myHESolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, a, *alpha);

	// Print the solved Black Scholes Equation into a gnuplot file
	myHESolver->printGrid(*alpha, 1000, "solvedHeat.gnuplot");

	delete myHESolver;
	delete myBoundingBox;
	delete alpha;
}


int main(int argc, char *argv[])
{
	testHeatEquation();
}
