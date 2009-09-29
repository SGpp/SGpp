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

#ifndef HEATEQUATIONTIMESTEPMATRIX_HPP
#define HEATEQUATIONTIMESTEPMATRIX_HPP

#include "data/DataVector.hpp"
#include "grid/Grid.hpp"
#include "operation/common/OperationODESolverMatrix.hpp"

#include <string>

namespace sg
{

/**
 * @todo (heinecke) add description here
 */
class HeatEquationTimestepMatrix : public OperationODESolverMatrix
{
private:
	/// the heat coefficient
	double a;
	/// the Laplace Operation (Stiffness Matrix)
	OperationMatrix* OpLaplace;
	/// the LTwoDotProduct Operation (Mass Matrix)
	OperationMatrix* OpMass;
	/**
	 *  specifies in which solver this matrix is used, valid values are:
	 *  ExEul for explicit Euler
	 *  ImEul for implicit Euler
	 *  CrNic for Crank Nicolson solver
	 */
	std::string tOperationMode;
	/// the size of one timestep used in the ODE Solver
	double TimestepSize;
	/// Pointer to the grid object
	Grid* myGrid;

	/**
	 * @todo (heinecke) add description
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	void applyMassMatrix(DataVector& alpha, DataVector& result);

	/**
	 * @todo (heinecke) add description
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	void applyLOperator(DataVector& alpha, DataVector& result);

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
	HeatEquationTimestepMatrix(Grid& SparseGrid, double a, double TimestepSize, std::string OperationMode = "ExEul");

	/**
	 * Std-Destructor
	 */
	virtual ~HeatEquationTimestepMatrix();

	virtual void mult(DataVector& alpha, DataVector& result);

	virtual void generateRHS(DataVector& data, DataVector& rhs);

	virtual void finishTimestep(DataVector& alpha);

	virtual void startTimestep(DataVector& alpha);

	virtual Grid* getGrid();
};

}

#endif /* HEATEQUATIONTIMESTEPMATRIX_HPP */
