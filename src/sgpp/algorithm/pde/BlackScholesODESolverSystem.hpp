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

#ifndef BLACKSCHOLESODESOLVERSYSTEM_HPP
#define BLACKSCHOLESODESOLVERSYSTEM_HPP

#include "grid/Grid.hpp"
#include "data/DataVector.hpp"
#include "operation/common/OperationODESolverSystem.hpp"
#include "grid/common/DirichletUpdateVector.hpp"
#include "grid/common/DirichletGridConverter.hpp"

namespace sg
{

/**
 * This class implements the ODESolverSystem for the BlackScholes
 * Equation.
 */
class BlackScholesODESolverSystem : public OperationODESolverSystem
{
private:
	/// the riskfree interest rate
	double r;
	/// the delta Operation, on boundary grid
	OperationMatrix* OpDeltaBound;
	/// the Gamma Operation, on boundary grid
	OperationMatrix* OpGammaBound;
	/// the LTwoDotProduct Operation (Mass Matrix), on boundary grid
	OperationMatrix* OpLTwoBound;
	/// the delta Operation, on Inner grid
	OperationMatrix* OpDeltaInner;
	/// the Gamma Operation, on Inner grid
	OperationMatrix* OpGammaInner;
	/// the LTwoDotProduct Operation (Mass Matrix), on Inner grid
	OperationMatrix* OpLTwoInner;
	/// Pointer to the alphas (ansatzfunctions' coefficients)
	DataVector* alpha_complete;
	/// Pointer to the alphas (ansatzfunctions' coefficients; inner points only)
	DataVector* alpha_inner;
	/// Pointer to the mus
	DataVector* mus;
	/// Pointer to the sigmas
	DataVector* sigmas;
	/// Pointer to the rhos;
	DataVector* rhos;
	/// Pointer to the coefficients of operation Delta
	DataVector* deltaCoef;
	/// Pointer to the coefficients ot operation Gamma
	DataVector* gammaCoef;
	/// Routine to modify the boundaries of the grid
	DirichletUpdateVector* BoundaryUpdate;
	/// Class that allows a simple conversion between a grid with and a without boundary points
	DirichletGridConverter* GridConverter;
	/// DateVector to store the right hand side
	DataVector* rhs;

	/**
	 *  specifies in which solver this matrix is used, valid values are:
	 *  ExEul for explicit Euler
	 *  ImEul for implicit Euler
	 *  CrNic for Crank Nicolson solver
	 */
	std::string tOperationMode;
	/// the size of one timestep used in the ODE Solver
	double TimestepSize;
	/// Pointer to the boundary grid object
	Grid* BoundGrid;
	/// Pointer to the inner grid object
	Grid* InnerGrid;

	/**
	 * Do Matrix mutlitplication with the Black Scholes Systemmatrix, inner grid points only
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions (inner grid points)
	 * @param return reference to the DataVector into which the result is written (inner grid points)
	 */
	void applyLOperatorInner(DataVector& alpha, DataVector& result);

	/**
	 * Do Matrix mutlitplication with the Black Scholes Systemmatrix, complete boundary grid
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	void applyLOperatorComplete(DataVector& alpha, DataVector& result);

	/**
	 * Do Matrix mutlitplication with left hand side mass matrix, inner grid points only
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions (inner grid points)
	 * @param return reference to the DataVector into which the result is written (inner grid points)
	 */
	void applyMassMatrixInner(DataVector& alpha, DataVector& result);

	/**
	 * Do Matrix mutlitplication with left hand side mass matrix
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	void applyMassMatrixComplete(DataVector& alpha, DataVector& result);

	/**
	 * Build the coefficients for the Gamma Operation, which
	 * are the assets' covariance matrix multiplied by 0.5
	 *
	 * this routine handles also the symmtrie of the
	 * gamma operation
	 */
	void buildGammaCoefficients();

	/**
	 * Build the coefficients for the combined Delta Operation
	 */
	void buildDeltaCoefficients();

	/**
	 * Implements some start jobs of every timestep, e.g.discounting boundaries
	 */
	void startTimestep();

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
	BlackScholesODESolverSystem(Grid& SparseGrid, DataVector& alpha, DataVector& mu, DataVector& sigma, DataVector& rho, double r, double TimestepSize, std::string OperationMode = "ExEul");

	/**
	 * Std-Destructor
	 */
	virtual ~BlackScholesODESolverSystem();

	virtual void mult(DataVector& alpha, DataVector& result);

	virtual DataVector* generateRHS();

	virtual void finishTimestep();

	virtual Grid* getGrid();

	virtual DataVector* getGridCoefficientsForCG();

	virtual DataVector* getGridCoefficients();
};

}

#endif /* BLACKSCHOLESODESOLVERSYSTEM_HPP */
