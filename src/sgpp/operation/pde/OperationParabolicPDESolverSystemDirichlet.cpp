/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "operation/pde/OperationParabolicPDESolverSystemDirichlet.hpp"
#include "exception/algorithm_exception.hpp"
using namespace sg::base;

namespace sg
{

OperationParabolicPDESolverSystemDirichlet::OperationParabolicPDESolverSystemDirichlet()
{
	this->numSumGridpointsInner = 0;
	this->numSumGridpointsComplete = 0;
}

OperationParabolicPDESolverSystemDirichlet::~OperationParabolicPDESolverSystemDirichlet()
{
}

void OperationParabolicPDESolverSystemDirichlet::mult(DataVector& alpha, DataVector& result)
{
	if (this->tOperationMode == "ExEul")
	{
		applyMassMatrixInner(alpha, result);
	}
	else if (this->tOperationMode == "ImEul")
	{
		result.setAll(0.0);

		DataVector temp(alpha.getSize());

		applyMassMatrixInner(alpha, temp);
		result.add(temp);

		applyLOperatorInner(alpha, temp);
		result.axpy((-1.0)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "CrNic")
	{
		result.setAll(0.0);

		DataVector temp(alpha.getSize());

		applyMassMatrixInner(alpha, temp);
		result.add(temp);

		applyLOperatorInner(alpha, temp);
		result.axpy((-0.5)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "AdBas")
	{
		result.setAll(0.0);

		applyMassMatrixInner(alpha, result);
	}
	else if (this->tOperationMode == "BDF2")
	{
		double tDiff=this->TimestepSize/this->TimestepSize_old;
		double alpha0 = (2.0*tDiff+1.0)/(tDiff+1.0);
		result.setAll(0.0);

		DataVector temp(alpha.getSize());

		applyMassMatrixInner(alpha, temp);

		temp.mult(alpha0);
		result.add(temp);

		applyLOperatorInner(alpha, temp);
		result.axpy((-1.0)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "F23")
	{
		result.setAll(0.0);
		double tDiff=this->TimestepSize/this->TimestepSize_old;
		double alpha0 = 1.0/(1.0+tDiff);

		applyMassMatrixInner(alpha, result);
		result.mult(alpha0);

	}
	else
	{
		throw new algorithm_exception("OperationParabolicPDESolverSystem::mult : An unknown operation mode was specified!");
	}
}

DataVector* OperationParabolicPDESolverSystemDirichlet::generateRHS()
{
	DataVector rhs_complete(this->alpha_complete->getSize());

	if (this->tOperationMode == "ExEul")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(this->alpha_complete->getSize());

		applyMassMatrixComplete(*this->alpha_complete, temp);
		rhs_complete.add(temp);

		applyLOperatorComplete(*this->alpha_complete, temp);
		rhs_complete.axpy(this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "ImEul")
	{
		rhs_complete.setAll(0.0);

		applyMassMatrixComplete(*this->alpha_complete, rhs_complete);
	}
	else if (this->tOperationMode == "CrNic")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(this->alpha_complete->getSize());

		applyMassMatrixComplete(*this->alpha_complete, temp);
		rhs_complete.add(temp);

		applyLOperatorComplete(*this->alpha_complete, temp);
		rhs_complete.axpy((0.5)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "AdBas")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(this->alpha_complete->getSize());

		applyMassMatrixComplete(*this->alpha_complete, temp);
		rhs_complete.add(temp);

		applyLOperatorComplete(*this->alpha_complete, temp);

		temp.mult((2.0)+this->TimestepSize/this->TimestepSize_old);

		DataVector temp_old(this->alpha_complete->getSize());
		applyMassMatrixComplete(*this->alpha_complete_old, temp_old);
		applyLOperatorComplete(*this->alpha_complete_old, temp_old);
		temp_old.mult(this->TimestepSize/this->TimestepSize_old);
		temp.sub(temp_old);

		rhs_complete.axpy((0.5)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "BDF2")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(this->alpha_complete->getSize());

		applyMassMatrixComplete(*this->alpha_complete, temp);

		double tDiff = this->TimestepSize/this->TimestepSize_old;

		double alpha1 = tDiff +1.0;
		temp.mult(alpha1);
		rhs_complete.add(temp);

		DataVector temp_old(this->alpha_complete->getSize());
		applyMassMatrixComplete(*this->alpha_complete_old, temp_old);

		double alpha2 = tDiff*tDiff/(1.0+tDiff);
		temp_old.mult(alpha2);
		rhs_complete.sub(temp_old);
	}
	else if (this->tOperationMode == "F23")
	{
		rhs_complete.setAll(0.0);
		double tDiff = this->TimestepSize/this->TimestepSize_old;
		double alpha0 = (1.0+tDiff);
		double alpha1 = alpha0*(tDiff -1.0);
		double alpha2 = -alpha0*(tDiff*tDiff/(tDiff+1.0));


		DataVector temp(this->alpha_complete->getSize());
		DataVector temp_old(this->alpha_complete->getSize());

		applyMassMatrixComplete(*this->alpha_complete, temp);
		temp.mult(alpha1);
		rhs_complete.sub(temp);

		applyLOperatorComplete(*this->alpha_complete, temp);
		temp.mult(alpha0*this->TimestepSize);

		applyMassMatrixComplete(*this->alpha_complete_old, temp_old);
		temp_old.mult(alpha2);
		rhs_complete.sub(temp_old);

		rhs_complete.add(temp);
	}
	else
	{
		throw new algorithm_exception("OperationParabolicPDESolverSystem::generateRHS : An unknown operation mode was specified!");
	}

	// Now we have the right hand side, lets apply the riskfree rate for the next timestep
	this->startTimestep();

	// Now apply the boundary ansatzfunctions to the inner ansatzfunctions
	DataVector result_complete(this->alpha_complete->getSize());
	DataVector alpha_bound(*this->alpha_complete);

	result_complete.setAll(0.0);

	this->BoundaryUpdate->setInnerPointsToZero(alpha_bound);

	// apply CG Matrix
	if (this->tOperationMode == "ExEul")
	{
		applyMassMatrixComplete(alpha_bound, result_complete);
	}
	else if (this->tOperationMode == "ImEul")
	{
		DataVector temp(alpha_bound.getSize());

		applyMassMatrixComplete(alpha_bound, temp);
		result_complete.add(temp);

		applyLOperatorComplete(alpha_bound, temp);
		result_complete.axpy((-1.0)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "CrNic")
	{
		DataVector temp(alpha_bound.getSize());

		applyMassMatrixComplete(alpha_bound, temp);
		result_complete.add(temp);

		applyLOperatorComplete(alpha_bound, temp);
		result_complete.axpy((-0.5)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "AdBas")
	{
		applyMassMatrixComplete(alpha_bound, result_complete);
	}
	else if (this->tOperationMode == "BDF2")
	{
		double tDiff =this->TimestepSize/this->TimestepSize_old;
		double alpha0 = (2.0*tDiff+1.0)/(tDiff+1.0);
		DataVector temp(alpha_bound.getSize());
		applyMassMatrixComplete(alpha_bound, temp);
		temp.mult(alpha0);
		result_complete.add(temp);

		applyLOperatorComplete(alpha_bound, temp);
		result_complete.axpy((-1.0)*this->TimestepSize, temp);

	}
	else if (this->tOperationMode == "F23")
	{

		double tDiff=this->TimestepSize/this->TimestepSize_old;
		double alpha0 = 1.0/(1.0+tDiff);

		applyMassMatrixComplete(alpha_bound, result_complete);
		result_complete.mult(alpha0);

	}
	else
	{
		throw new algorithm_exception("OperationParabolicPDESolverSystem::generateRHS : An unknown operation mode was specified!");
	}
	rhs_complete.sub(result_complete);

	if (this->rhs != NULL)
	{
		delete this->rhs;
	}

	this->rhs = new DataVector(this->alpha_inner->getSize());
	this->GridConverter->calcInnerCoefs(rhs_complete, *this->rhs);

	return this->rhs;
}

DataVector* OperationParabolicPDESolverSystemDirichlet::getGridCoefficientsForCG()
{
	this->GridConverter->calcInnerCoefs(*this->alpha_complete, *this->alpha_inner);

	return this->alpha_inner;
}

}
