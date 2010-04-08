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

#include "algorithm/pde/HeatEquationODESolverSystem.hpp"
#include "exception/algorithm_exception.hpp"

namespace sg
{

HeatEquationODESolverSystem::HeatEquationODESolverSystem(Grid& SparseGrid, DataVector& alpha, double a, double TimestepSize, std::string OperationMode)
{
	this->OpLaplace = SparseGrid.createOperationLaplace();
	this->OpMass = SparseGrid.createOperationLTwoDotProduct();
	this->a = a;
	this->tOperationMode = OperationMode;
	this->TimestepSize = TimestepSize;
	this->myGrid = &SparseGrid;
	this->alpha_complete = &alpha;
	this->alpha_inner = &alpha;
}

HeatEquationODESolverSystem::~HeatEquationODESolverSystem()
{
	delete this->OpLaplace;
	delete this->OpMass;
}

void HeatEquationODESolverSystem::mult(DataVector& alpha, DataVector& result)
{
	if (this->tOperationMode == "ExEul")
	{
		result.setAll(0.0);

		applyMassMatrix(alpha, result);
	}
	else if (this->tOperationMode == "ImEul")
	{
		result.setAll(0.0);

		DataVector temp(alpha.getSize());

		temp.setAll(0.0);
		applyMassMatrix(alpha, temp);
		result.add(temp);

		temp.setAll(0.0);
		applyLOperator(alpha, temp);
		result.axpy((-1.0)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "CrNic")
	{
		result.setAll(0.0);

		DataVector temp(alpha.getSize());

		temp.setAll(0.0);
		applyMassMatrix(alpha, temp);
		result.add(temp);

		temp.setAll(0.0);
		applyLOperator(alpha, temp);
		result.axpy((-0.5)*this->TimestepSize, temp);
	}
	else
	{
		throw new algorithm_exception(" HeatEquationTimestepMatrix::mult : An unknown operation mode was specified!");
	}
}

void HeatEquationODESolverSystem::generateRHS(DataVector& data, DataVector& rhs)
{
	if (this->tOperationMode == "ExEul")
	{
		DataVector temp(data.getSize());
		temp.setAll(0.0);

		applyMassMatrix(data, temp);
		rhs.add(temp);

		temp.setAll(0.0);
		applyLOperator(data, temp);
		rhs.axpy(this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "ImEul")
	{
		rhs.setAll(0.0);

		applyMassMatrix(data, rhs);
	}
	else if (this->tOperationMode == "CrNic")
	{
		rhs.setAll(0.0);

		DataVector temp(data.getSize());

		temp.setAll(0.0);
		applyMassMatrix(data, temp);
		rhs.add(temp);

		temp.setAll(0.0);
		applyLOperator(data, temp);
		rhs.axpy((0.5)*this->TimestepSize, temp);
	}
	else
	{
		throw new algorithm_exception(" HeatEquationTimestepMatrix::generateRHS : An unknown operation mode was specified!");
	}
}

void HeatEquationODESolverSystem::applyMassMatrix(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	// Apply the mass matrix
	this->OpMass->mult(alpha, temp);

	result.add(temp);
}

void HeatEquationODESolverSystem::applyLOperator(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	// Apply the laplace Operator rate
	this->OpLaplace->mult(alpha, temp);
	result.axpy((-1.0)*this->a,temp);
}

void HeatEquationODESolverSystem::finishTimestep(DataVector& alpha)
{
}

void HeatEquationODESolverSystem::startTimestep(DataVector& alpha)
{
}

Grid* HeatEquationODESolverSystem::getGrid()
{
	return myGrid;
}

DataVector* HeatEquationODESolverSystem::getGridCoefficients()
{
	return this->alpha_complete;
}

DataVector* HeatEquationODESolverSystem::getGridCoefficientsForCG()
{
	return this->alpha_inner;
}

}
