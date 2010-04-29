/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/HeatEquationTimestepMatrix.hpp"
#include "exception/algorithm_exception.hpp"

namespace sg
{

HeatEquationTimestepMatrix::HeatEquationTimestepMatrix(Grid& SparseGrid, double a, double TimestepSize, std::string OperationMode)
{
	this->OpLaplace = SparseGrid.createOperationLaplace();
	this->OpMass = SparseGrid.createOperationLTwoDotProduct();
	this->a = a;
	this->tOperationMode = OperationMode;
	this->TimestepSize = TimestepSize;
	this->myGrid = &SparseGrid;
}

HeatEquationTimestepMatrix::~HeatEquationTimestepMatrix()
{
	delete this->OpLaplace;
	delete this->OpMass;
}

void HeatEquationTimestepMatrix::mult(DataVector& alpha, DataVector& result)
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

void HeatEquationTimestepMatrix::generateRHS(DataVector& data, DataVector& rhs)
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

void HeatEquationTimestepMatrix::applyMassMatrix(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	// Apply the mass matrix
	this->OpMass->mult(alpha, temp);

	result.add(temp);
}

void HeatEquationTimestepMatrix::applyLOperator(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	// Apply the laplace Operator rate
	this->OpLaplace->mult(alpha, temp);
	result.axpy((-1.0)*this->a,temp);
}

void HeatEquationTimestepMatrix::finishTimestep(DataVector& alpha)
{
}

void HeatEquationTimestepMatrix::startTimestep(DataVector& alpha)
{
}

Grid* HeatEquationTimestepMatrix::getGrid()
{
	return myGrid;
}

}
