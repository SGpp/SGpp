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

void testOneUnderlying()
{
	size_t dim = 1;
	size_t level = 5;
	double* strike = new double[dim];
	strike[0] = 65.0;

	size_t timesteps = 2500;
	double stepsize = 0.0001;
	size_t CGiterations = 150;
	double CGepsilon = 0.001;

	DataVector mu(1);
	DataVector sigma(1);
	DataVector rho(1);

	double r = 0.05;
	mu.set(0, 0.05);
	sigma.set(0, 0.15);
	rho.set(0, 1.0);

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];

	// set the bounding box
	for (size_t i = 0; i < dim; i++)
	{
		myBoundaries[i].leftBoundary = 0.0;
		myBoundaries[i].rightBoundary = 100.0;
		myBoundaries[i].bDirichletLeft = true;
		myBoundaries[i].bDirichletRight = true;
	}

	sg::BlackScholesSolver* myBSSolver = new sg::BlackScholesSolver();
	sg::BoundingBox* myBoundingBox = new sg::BoundingBox(dim, myBoundaries);
	delete[] myBoundaries;

	// Construct a grid
	myBSSolver->constructGrid(*myBoundingBox, level);

	// init the basis functions' coefficient vector
	DataVector* alpha = new DataVector(myBSSolver->getNumberGridPoints());

	// Init the grid with on payoff function
	myBSSolver->initGridWithPayoff(*alpha, strike);

	// Print the payoff function into a gnuplot file
	myBSSolver->printGrid(*alpha, 0.1, "payoff.gnuplot");

	// Set stochastic data
	myBSSolver->setStochasticData(mu, sigma, rho, r);

	// Start solving the Black Scholes Equation
	myBSSolver->solveEuler(timesteps, stepsize, *alpha);
	//myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha);

	// Print the solved Black Scholes Equation into a gnuplot file
	myBSSolver->printGrid(*alpha, 0.1, "solvedBS.gnuplot");

	delete myBSSolver;
	delete myBoundingBox;
	delete alpha;
}

void testTwoUnderlyings()
{
	size_t dim = 2;
	size_t level = 6;
	double* strike = new double[dim];
	strike[0] = 65.0;
	strike[1] = 55.0;

	size_t timesteps = 100;
	double stepsize = 0.0001;
	size_t CGiterations = 150;
	double CGepsilon = 0.001;

	DataVector mu(2);
	DataVector sigma(2);
	DataVector rho(2,2);

	double r = 0.00;
	mu.set(0, 0.00);
	mu.set(1, 0.00);
	sigma.set(0, 0.05);
	sigma.set(1, 0.05);
	rho.set(0, 1.0);
	rho.set(1, 0.1);
	rho.set(2, 0.1);
	rho.set(3, 1.0);

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];

	// set the bounding box
	for (size_t i = 0; i < dim; i++)
	{
		myBoundaries[i].leftBoundary = 0.0;
		myBoundaries[i].rightBoundary = 100.0;
		myBoundaries[i].bDirichletLeft = true;
		myBoundaries[i].bDirichletRight = true;
	}

	sg::BlackScholesSolver* myBSSolver = new sg::BlackScholesSolver();
	sg::BoundingBox* myBoundingBox = new sg::BoundingBox(dim, myBoundaries);
	delete[] myBoundaries;

	// Construct a grid
	myBSSolver->constructGrid(*myBoundingBox, level);

	// init the basis functions' coefficient vector
	DataVector* alpha = new DataVector(myBSSolver->getNumberGridPoints());

	// Init the grid with on payoff function
	myBSSolver->initGridWithPayoff(*alpha, strike);

	// Print the payoff function into a gnuplot file
	myBSSolver->printGrid(*alpha, 2, "payoff.gnuplot");

	// Set stochastic data
	myBSSolver->setStochasticData(mu, sigma, rho, r);

	// Start solving the Black Scholes Equation
	myBSSolver->solveEuler(timesteps, stepsize, *alpha);
	//myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha);

	// Print the solved Black Scholes Equation into a gnuplot file
	myBSSolver->printGrid(*alpha, 2, "solvedBS.gnuplot");

	delete myBSSolver;
	delete myBoundingBox;
	delete alpha;
}

int main(int argc, char *argv[])
{
	testTwoUnderlyings();
}
