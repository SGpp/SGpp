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

#include "grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "grid/common/BoundingBox.hpp"

#include "algorithm/finance/BlackScholesTimestepMatrix.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

namespace sg
{

/**
 * This class provides a simple-to-use solver of the multi dimensional Black
 * Scholes Equation that uses Sparse Grids.
 *
 * The class's aim is, to hide all complex details of solving the Black Scholes
 * Equation with Sparse Grids!
 */
class BlackScholesSolver
{
private:
	/// vector that contains the sparse's grid coefficients
	DataVector* alpha;
	/// vector that contains the expected values
	DataVector* mus;
	/// vector that contains the standard deviations
	DataVector* sigmas;
	/// Matrix that contains the correlations
	DataVector* rhos;
	/// the riskfree rate
	double r;
	/// the size of one timestep
	double timestepSize;
	/// The number of timesteps that are executed during solving
	size_t nTimesteps;
	/// The Sparse Grid needed in this classificator
	Grid* myGrid;
	/// the number of levels used for an regular grid
	size_t levels;
	/// the dimension of the grid
	size_t dim;
	/// stores if the stochastic asset data was passed to the solver
	bool bStochasticDataAlloc;
	/// stores if the grid was created inside the solver
	bool bGridConstructed;

public:
	/**
	 * Std-Constructor of the solver
	 */
	BlackScholeSolver();

	/**
	 * Std-Destructor of the solver
	 */
	~BlackScholesSolver();

	/**
	 * Use this routine the construct a regular grid to solve the multi-dimensional Black Scholes Equation
	 *
	 * @param myBoundingBox reference to a bounding box that describes the grid
	 * @param level number of the regular's grid levels
	 */
	void constructGrid(BoundingBox& myBoundingBox, size_t level);

	/**
	 * Use this routine if you want to solve a problem stored in the format provided by the solving system
	 * released by the University of Bonn, Germany
	 *
	 * @param tfilename absolute path of the file
	 */
	void constructGrid(std::string tfilename);

	/**
	 * In order to solve the multi dimensional Black Scholes Equation you have to provided
	 * some statistical data about the underlying (expected values, standard deviation
	 * and the correlation between them). This function allows you to set this data.
	 *
	 * @param mus a DataVector that contains the underlyings' expected values
	 * @param sigmas a DataVector that contains the underlyings' standard deviations
	 * @param rhos a DataVector that contains the correlations between the underlyings
	 * @param r the riskfree rate used in the market model
	 */
	void setStochasticData(DataVector& mus, DataVector& sigmas, DataVector& rhos, double r);

	/**
	 * Call this routine to use an explicit Euler algorithm to solve the multi dimensional
	 * Black Scholes Equation.
	 *
	 * @param numTimestpes the number of timesteps that should be executed
	 * @param timestepsize the size of the interval one timestep moves forward
	 * @param alpha the coefficients of the Sparse Gird's basis functions
	 */
	void solveEuler(size_t numTimesteps, double timestepsize, DataVector& alpha);

	/**
	 * This is some kind of debug functionality. It writes a file,
	 * that can be used with gnuplot the print the grid.
	 *
	 * Is only implemented for 1D and 2D grids!
	 *
	 * @param alpha the coefficients of the Sparse Gird's basis functions
	 * @param resolution the distance between evalution points
	 * @param filename absolute path to file into which the grid's evaluation is written
	 *
	 * @todo (heinecke) move this into a class that is located in the folder tool, e.g. GridPrinter
	 */
	void printGrid(DataVector& alpha, double resolution, std::string filename);
};

}
