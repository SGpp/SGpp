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

#ifndef HEATEQUATIONODESOLVERSYSTEM_HPP
#define HEATEQUATIONODESOLVERSYSTEM_HPP

#include "data/DataVector.hpp"
#include "grid/Grid.hpp"
#include "operation/common/OperationODESolverSystem.hpp"
#include "grid/common/DirichletUpdateVector.hpp"
#include "grid/common/DirichletGridConverter.hpp"

#include <string>

namespace sg
{

/**
 * This class implements the ODESolverSystem for the
 * Heat Equation.
 */
class HeatEquationODESolverSystem : public OperationODESolverSystem
{
private:
	/// the heat coefficient
	double a;
	/// Pointer to the alphas (ansatzfunctions' coefficients)
	DataVector* alpha_complete;
	/// Pointer to the alphas (ansatzfunctions' coefficients; inner points only)
	DataVector* alpha_inner;
	/// the Laplace Operation (Stiffness Matrix), on boundary grid
	OperationMatrix* OpLaplaceBound;
	/// the LTwoDotProduct Operation (Mass Matrix), on boundary grid
	OperationMatrix* OpMassBound;
	/// the Laplace Operation (Stiffness Matrix), on inner grid
	OperationMatrix* OpLaplaceInner;
	/// the LTwoDotProduct Operation (Mass Matrix), on inner grid
	OperationMatrix* OpMassInner;
	/**
	 *  specifies in which solver this matrix is used, valid values are:
	 *  ExEul for explicit Euler
	 *  ImEul for implicit Euler
	 *  CrNic for Crank Nicolson solver
	 */
	std::string tOperationMode;
	/// the size of one timestep used in the ODE Solver
	double TimestepSize;
	/// Routine to modify the boundaries/inner points of the grid
	DirichletUpdateVector* BoundaryUpdate;
	/// Class that allows a simple conversion between a grid with and a without boundary points
	DirichletGridConverter* GridConverter;
	/// DateVector to store the right hand side
	DataVector* rhs;
	/// Pointer to the grid object
	Grid* BoundGrid;
	/// Pointer to the inner grid object
	Grid* InnerGrid;


	/**
	 * applies the mass matrix of the Heat Equation, on complete grid - with boundaries
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	void applyMassMatrixComplete(DataVector& alpha, DataVector& result);

	/**
	 * applies the system matrix of the Heat Equation, on complete grid - with boundaries
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	void applyLOperatorComplete(DataVector& alpha, DataVector& result);

	/**
	 * applies the mass matrix of the Heat Equation, on inner grid only
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	void applyMassMatrixInner(DataVector& alpha, DataVector& result);

	/**
	 * applies the system matrix of the Heat Equation, on inner grid only
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	void applyLOperatorInner(DataVector& alpha, DataVector& result);

	/**
	 * Implements some adjustments needed before soling a timestep
	 */
	void startTimestep();

public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param a the heat coefficient
	 * @param TimestepSize the size of one timestep used in the ODE Solver
	 * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for explicit Euler,
	 *  							ImEul for implicit Euler, CrNic for Crank Nicolson solver
	 */
	HeatEquationODESolverSystem(Grid& SparseGrid, DataVector& alpha, double a, double TimestepSize, std::string OperationMode = "ExEul");

	/**
	 * Std-Destructor
	 */
	virtual ~HeatEquationODESolverSystem();

	virtual void mult(DataVector& alpha, DataVector& result);

	virtual DataVector* generateRHS();

	virtual void finishTimestep();

	virtual Grid* getGrid();

	virtual DataVector* getGridCoefficientsForCG();

	virtual DataVector* getGridCoefficients();
};

}

#endif /* HEATEQUATIONODESOLVERSYSTEM_HPP */
