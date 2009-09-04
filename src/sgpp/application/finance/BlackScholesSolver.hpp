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
#include "exception/operation_exception.hpp"
#include "grid/type/LinearGrid.hpp"
#include "grid/type/LinearBoundaryGrid.hpp"
#include "grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "grid/type/ModLinearGrid.hpp"

#include "algorithm/finance/BlackScholesTimestepMatrix.hpp"

#include <iostream>

namespace sg
{

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
	/// the type of grid used
	std::string GridType;
	/// the dimension of the grid
	size_t dim;

public:
	BlackScholeSolver();

	~BlackScholesSolver();

	void constructGrid(size_t dim, size_t level);

	void setStatData(DataVector& mus, DataVector& sigmas, DataVector& rhos, double r);

	void solveEuler(size_t numTimesteps, double timestepsize, DataVector& alpha);

	void printGrid(DataVector& alpha);

	void get1Dpayoff(DataVector& alpha, size_t levels, double strike);

};

}
