/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "operation/pde/OperationODESolverSystemNeumann.hpp"
#include "exception/algorithm_exception.hpp"

namespace sg
{

OperationODESolverSystemNeumann::OperationODESolverSystemNeumann()
{
	this->numSumGridpointsInner = 0;
	this->numSumGridpointsComplete = 0;
}

OperationODESolverSystemNeumann::~OperationODESolverSystemNeumann()
{
}

void OperationODESolverSystemNeumann::mult(DataVector& alpha, DataVector& result)
{
	if (this->tOperationMode == "ExEul")
	{
		applyMassMatrix(alpha, result);
	}
	else if (this->tOperationMode == "ImEul")
	{
		result.setAll(0.0);

		DataVector temp(alpha.getSize());

		applyMassMatrix(alpha, temp);
		result.add(temp);

		applyLOperator(alpha, temp);
		result.axpy((-1.0)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "CrNic")
	{
		result.setAll(0.0);

		DataVector temp(alpha.getSize());

		applyMassMatrix(alpha, temp);
		result.add(temp);

		applyLOperator(alpha, temp);
		result.axpy((-0.5)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "AdBas")
	{
		result.setAll(0.0);

		applyMassMatrix(alpha, result);
	}
	else
	{
		throw new algorithm_exception("OperationODESolverSystemNeumann::mult : An unknown operation mode was specified!");
	}
}

DataVector* OperationODESolverSystemNeumann::generateRHS()
{
	DataVector rhs_complete(this->alpha_complete->getSize());

	if (this->tOperationMode == "ExEul")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(this->alpha_complete->getSize());

		applyMassMatrix(*this->alpha_complete, temp);
		rhs_complete.add(temp);

		applyLOperator(*this->alpha_complete, temp);
		rhs_complete.axpy(this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "ImEul")
	{
		rhs_complete.setAll(0.0);

		applyMassMatrix(*this->alpha_complete, rhs_complete);
	}
	else if (this->tOperationMode == "CrNic")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(this->alpha_complete->getSize());

		applyMassMatrix(*this->alpha_complete, temp);
		rhs_complete.add(temp);

		applyLOperator(*this->alpha_complete, temp);
		rhs_complete.axpy((0.5)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "AdBas")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(this->alpha_complete->getSize());

		applyMassMatrix(*this->alpha_complete, temp);
		rhs_complete.add(temp);

		applyLOperator(*this->alpha_complete, temp);

		temp.mult((2.0)+this->TimestepSize/this->TimestepSize_old);

		DataVector temp_old(this->alpha_complete->getSize());
		applyMassMatrix(*this->alpha_complete_old, temp_old);
		applyLOperator(*this->alpha_complete_old, temp_old);
		temp_old.mult(this->TimestepSize/this->TimestepSize_old);
		temp.sub(temp_old);

		rhs_complete.axpy((0.5)*this->TimestepSize, temp);
	}
	else
	{
		throw new algorithm_exception("OperationODESolverSystemNeumann::generateRHS : An unknown operation mode was specified!");
	}

	// Now we have the right hand side, lets apply the riskfree rate for the next timestep
	this->startTimestep();

	if (this->rhs != NULL)
	{
		delete this->rhs;
	}

	this->rhs = new DataVector(rhs_complete);

	return this->rhs;
}

/*
void OperationODESolverSystemNeumann::finishTimestep(bool isLastTimestep)
{
	// Replace the inner coefficients on the boundary grid
	this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);
}

void OperationODESolverSystemNeumann::startTimestep()
{
}
*/

DataVector* OperationODESolverSystemNeumann::getGridCoefficientsForCG()
{
	return this->alpha_complete;
}

}

