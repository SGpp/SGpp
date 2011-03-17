/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "operation/pde/OperationEllipticPDESolverSystemDirichlet.hpp"
#include "exception/algorithm_exception.hpp"

#include <sstream>

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

size_t OperationEllipticPDESolverSystemDirichlet::getMatrix(std::string& mtxString, bool complete)
{
	size_t vector_size;
	if (complete == true)
	{
		vector_size = this->numGridpointsComplete;
	}
	else
	{
		vector_size = this->numGridpointsInner;
	}

	DataVector alpha(vector_size);
	DataVector result(vector_size);
	std::stringstream mtxStream;
	std::stringstream mtxHeaderStream;
	size_t nonZeros = 0;

	mtxHeaderStream.clear();
	mtxStream.clear();

	// loop over the system matrices columns
	for (size_t i = 0; i < vector_size; i++)
	{
		alpha.setAll(0.0);
		result.setAll(0.0);
		alpha.set(i, 1.0);

		// calculate column via Up/Down
		if (complete == true)
		{
			applyLOperatorComplete(alpha, result);
		}
		else
		{
			applyLOperatorInner(alpha, result);
		}

		// serialize result into mtxStream
		for (size_t j = 0; j < vector_size; j++)
		{
			if (result[j] != 0.0)
			{
				mtxStream << j << " " << i << " " << result[j] << std::endl;
				nonZeros++;
			}
		}
	}

	// Generate Header Line
	mtxHeaderStream << vector_size << " " << vector_size << " " << nonZeros << std::endl;

	mtxString = mtxHeaderStream.str() + mtxStream.str();

	return nonZeros;
}

size_t OperationEllipticPDESolverSystemDirichlet::getCompleteMatrix(std::string& mtxString)
{
	return getMatrix(mtxString, true);
}

size_t OperationEllipticPDESolverSystemDirichlet::getInnerMatrix(std::string& mtxString)
{
	return getMatrix(mtxString, false);
}

}
