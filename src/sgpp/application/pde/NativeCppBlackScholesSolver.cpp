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
#include <iostream>
#include <string>
#include <stdlib.h>

void testOneUnderlying(size_t l, std::string fileStoch, std::string fileBound, double strike1, double riskfree, size_t timeSt,
						double dt, size_t CGIt, double CGeps, bool animation)
{
	size_t dim = 1;
	size_t level = l;
	double* strike = new double[dim];
	strike[0] = strike1;

	size_t timesteps = timeSt;
	double stepsize = dt;
	size_t CGiterations = CGIt;
	double CGepsilon = CGeps;

	DataVector mu(1);
	DataVector sigma(1);
	DataVector rho(1);

	double r = riskfree;

	mu.set(0, 0.00);
	sigma.set(0, 0.40);
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

	// init Screen Object
	myBSSolver->initScreen();

	// Construct a grid
	myBSSolver->constructGrid(*myBoundingBox, level);

	// init the basis functions' coefficient vector
	DataVector* alpha = new DataVector(myBSSolver->getNumberGridPoints());

	// Init the grid with on payoff function
	myBSSolver->initGridWithPayoff(*alpha, strike);

	// Print the payoff function into a gnuplot file
	myBSSolver->printGrid(*alpha, 50, "payoff.gnuplot");

	// Set stochastic data
	myBSSolver->setStochasticData(mu, sigma, rho, r);

	// Start solving the Black Scholes Equation
	//myBSSolver->solveExplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, animation, 50);
	myBSSolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, animation, 50);
	//myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha);

	// Print the solved Black Scholes Equation into a gnuplot file
	myBSSolver->printGrid(*alpha, 50, "solvedBS.gnuplot");

	// Do analytic test
	std::vector< std::pair<double, double> >premium;
	double t = (((double)timesteps)*stepsize);
	myBSSolver->solve1DAnalytic(premium, myBoundaries[0].rightBoundary, myBoundaries[0].rightBoundary/50.0, strike[0], t);
	myBSSolver->print1DAnalytic(premium, "analyticBS.gnuplot");

	delete[] myBoundaries;
	delete myBSSolver;
	delete myBoundingBox;
	delete alpha;
}

void testTwoUnderlyings(size_t l, std::string fileStoch, std::string fileBound, double strike1, double strike2, double riskfree, size_t timeSt,
		double dt, size_t CGIt, double CGeps, bool animation)
{
	size_t dim = 2;
	size_t level = l;
	double* strike = new double[dim];
	strike[0] = strike1;
	strike[1] = strike2;

	size_t timesteps = timeSt;
	double stepsize = dt;
	size_t CGiterations = CGIt;
	double CGepsilon = CGeps;

	DataVector mu(2);
	DataVector sigma(2);
	DataVector rho(2,2);

	double r = riskfree;

	mu.set(0, 0.00);
	mu.set(1, 0.00);
	sigma.set(0, 0.40);
	sigma.set(1, 0.60);
	rho.set(0, 1.0);
	rho.set(1, 0.1);
	rho.set(2, 0.1);
	rho.set(3, 1.0);

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];

	// set the bounding box
	for (size_t i = 0; i < dim; i++)
	{
		myBoundaries[i].leftBoundary = 0.0;
		myBoundaries[i].rightBoundary = 1.0;
		myBoundaries[i].bDirichletLeft = true;
		myBoundaries[i].bDirichletRight = true;
	}

	sg::BlackScholesSolver* myBSSolver = new sg::BlackScholesSolver();
	sg::BoundingBox* myBoundingBox = new sg::BoundingBox(dim, myBoundaries);
	delete[] myBoundaries;

	// init Screen Object
	myBSSolver->initScreen();

	// Construct a grid
	myBSSolver->constructGrid(*myBoundingBox, level);

	// init the basis functions' coefficient vector
	DataVector* alpha = new DataVector(myBSSolver->getNumberGridPoints());

	// Init the grid with on payoff function
	myBSSolver->initGridWithPayoff(*alpha, strike);

	// Print the payoff function into a gnuplot file
	myBSSolver->printGrid(*alpha, 20, "payoff.gnuplot");

	// Set stochastic data
	myBSSolver->setStochasticData(mu, sigma, rho, r);

	// Start solving the Black Scholes Equation
	//myBSSolver->solveExplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, animation, 20);
	myBSSolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, animation, 20);
	//myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha);

	// Print the solved Black Scholes Equation into a gnuplot file
	myBSSolver->printGrid(*alpha, 20, "solvedBS.gnuplot");

	delete myBSSolver;
	delete myBoundingBox;
	delete alpha;
}

void writeHelp()
{
	sg::BlackScholesSolver* myBSSolver = new sg::BlackScholesSolver();

	// init Screen Object
	myBSSolver->initScreen();

	// print Help
	myBSSolver->writeHelp();

	delete myBSSolver;
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

	if (option == "test1D")
	{
		if (argc != 12)
		{
			writeHelp();
		}
		else
		{
			std::string fileStoch;
			std::string fileBound;
			std::string ani;
			bool animation;

			fileStoch.assign(argv[4]);
			fileBound.assign(argv[3]);
			ani.assign(argv[11]);
			if (ani == "animation")
			{
				animation = true;
			}
			else
			{
				animation = false;
			}
			testOneUnderlying(atoi(argv[2]), fileStoch, fileBound, atof(argv[5]), atof(argv[6]), (size_t)(atof(argv[7])/atof(argv[8])), atof(argv[8]), atoi(argv[9]), atof(argv[10]), animation);
		}
	}
	else if (option == "test2D")
	{
		if (argc != 13)
		{
			writeHelp();
		}
		else
		{
			std::string fileStoch;
			std::string fileBound;
			std::string ani;
			bool animation;

			fileStoch.assign(argv[4]);
			fileBound.assign(argv[3]);
			ani.assign(argv[12]);
			if (ani == "animation")
			{
				animation = true;
			}
			else
			{
				animation = false;
			}
			testTwoUnderlyings(atoi(argv[2]), fileStoch, fileBound, atof(argv[5]), atof(argv[6]), atof(argv[7]), (size_t)(atof(argv[8])/atof(argv[9])), atof(argv[9]), atoi(argv[10]), atof(argv[11]), animation);
		}
	}
	else if (option == "solveBonn")
	{

	}
	else
	{
		writeHelp();
	}

	return 0;
}
