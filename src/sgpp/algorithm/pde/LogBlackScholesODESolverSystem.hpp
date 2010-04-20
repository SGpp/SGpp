/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
/*               2010      Stefanie Schraufstetter (schraufs@in.tum.de)      */
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

#ifndef LOGBLACKSCHOLESODESOLVERSYSTEM_HPP
#define LOGBLACKSCHOLESODESOLVERSYSTEM_HPP

#include "grid/Grid.hpp"
#include "data/DataVector.hpp"
#include "algorithm/pde/BlackScholesODESolverSystem.hpp"
#include "grid/common/DirichletUpdateVector.hpp"


namespace sg
{

/**
 * implements the ODESolverSystem for the log-transformed Black-Scholes equation
 */
class LogBlackScholesODESolverSystem : public BlackScholesODESolverSystem
{
private:
	/**
	 * Build the coefficients for the Gamma Operation
	 * (term for second derivative in PDE), which are the
	 * assets' covariance matrix multiplied by 0.5
	 *
	 * this routine handles also the symmetry of the
	 * gamma operation
	 */
	void buildGammaCoefficients();

	/**
	 * Build the coefficients for the combined Delta Operation
	 * (term for first derivative in PDE), which are the drift factors
	 */
	void buildDeltaCoefficients();


public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param alpha the ansatzfunctions' coefficients
	 * @param mu reference to the mus
	 * @param sigma reference to the sigmas
	 * @param rho reference to the rhos
	 * @param r the riskfree interest rate
	 * @param TimestepSize the size of one timestep used in the ODE Solver
	 * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for explicit Euler,
	 *  							ImEul for implicit Euler, CrNic for Crank Nicolson solver
	 */
	LogBlackScholesODESolverSystem(Grid& SparseGrid, DataVector& alpha, DataVector& mu, DataVector& sigma, DataVector& rho, double r, double TimestepSize, std::string OperationMode = "ExEul");
};

}

#endif /* LOGBLACKSCHOLESODESOLVERSYSTEM_HPP */
