/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "operation/pde/OperationEllipticPDESolverSystemDirichlet.hpp"
#include "exception/algorithm_exception.hpp"

namespace sg
{

OperationEllipticPDESolverSystemDirichlet::OperationEllipticPDESolverSystemDirichlet(Grid& SparseGrid, DataVector& rhs) : OperationEllipticPDESolverSystem(SparseGrid, rhs)
{
	this->BoundaryUpdate = new DirichletUpdateVector(SparseGrid.getStorage());
	this->GridConverter = new DirichletGridConverter();

	this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->rhs, &(this->InnerGrid), &(this->rhs_inner));

	this->numGridpointsInner = this->InnerGrid->getSize();

	this->alpha_inner = NULL;
}

OperationEllipticPDESolverSystemDirichlet::~OperationEllipticPDESolverSystemDirichlet()
{
	delete this->alpha_inner;
	delete this->rhs_inner;
	delete this->InnerGrid;
	delete this->BoundaryUpdate;
	delete this->GridConverter;
}

void OperationEllipticPDESolverSystemDirichlet::mult(DataVector& alpha, DataVector& result)
{
	applyLOperatorInner(alpha, result);
}

DataVector* OperationEllipticPDESolverSystemDirichlet::generateRHS()
{
	if (this->InnerGrid != NULL)
	{
		DataVector alpha_tmp_complete(*(this->rhs));
		DataVector rhs_tmp_complete(*(this->rhs));

		this->BoundaryUpdate->setInnerPointsToZero(alpha_tmp_complete);
		applyLOperatorComplete(alpha_tmp_complete, rhs_tmp_complete);

		this->GridConverter->calcInnerCoefs(rhs_tmp_complete, *(this->rhs_inner));
		this->rhs_inner->mult(-1.0);
	}
	else
	{
		throw new algorithm_exception("OperationEllipticPDESolverSystemDirichlet::generateRHS : No inner grid exists!");
	}

	return this->rhs_inner;
}

DataVector* OperationEllipticPDESolverSystemDirichlet::getGridCoefficientsForCG()
{
	if (this->InnerGrid != NULL)
	{
		if (this->alpha_inner != NULL)
		{
			delete this->alpha_inner;
		}

		this->alpha_inner = new DataVector(this->InnerGrid->getSize());
		this->alpha_inner->setAll(0.0);
	}
	else
	{
		throw new algorithm_exception("OperationEllipticPDESolverSystemDirichlet::getGridCoefficientsForCG : No inner grid exists!");
	}

	return this->alpha_inner;
}

void OperationEllipticPDESolverSystemDirichlet::getSolutionBoundGrid(DataVector& Solution, DataVector& SolutionInner)
{
	Solution = *(this->rhs);
	this->GridConverter->updateBoundaryCoefs(Solution, SolutionInner);
}

}
