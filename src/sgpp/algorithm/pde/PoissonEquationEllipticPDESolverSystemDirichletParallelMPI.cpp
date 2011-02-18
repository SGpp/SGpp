/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "tools/MPI/SGppMPITools.hpp"

#include "algorithm/pde/PoissonEquationEllipticPDESolverSystemDirichletParallelMPI.hpp"
#include "exception/algorithm_exception.hpp"

#include "algorithm/pde/StdUpDown.hpp"
#include "algorithm/pde/UpDownOneOpDim.hpp"

namespace sg
{

PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::PoissonEquationEllipticPDESolverSystemDirichletParallelMPI(Grid& SparseGrid, DataVector& rhs) : OperationEllipticPDESolverSystemDirichlet(SparseGrid, rhs)
{
	this->Laplace_Complete = this->BoundGrid->createOperationLaplace();
	this->Laplace_Inner = this->InnerGrid->createOperationLaplace();
}

PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::~PoissonEquationEllipticPDESolverSystemDirichletParallelMPI()
{
	delete this->Laplace_Complete;
	delete this->Laplace_Inner;
}

void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::applyLOperatorInner(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
	size_t nDims = algoDims.size();

	// Apply Laplace, parallel in Dimensions
	for (size_t i = 0; i < nDims; i++)
	{
		if (i % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank())
		{
			DataVector myResult(result.getSize());

			/// @todo (heinecke) discuss methods in order to avoid this cast
			((UpDownOneOpDim*)(this->Laplace_Inner))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

			result.add(myResult);
		}
	}

}

void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::applyLOperatorComplete(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
	size_t nDims = algoDims.size();

	// Apply Laplace, parallel in Dimensions
	for (size_t i = 0; i < nDims; i++)
	{
		if (i % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank())
		{
			DataVector myResult(result.getSize());

			/// @todo (heinecke) discuss methods in order to avoid this cast
			((UpDownOneOpDim*)(this->Laplace_Complete))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

			result.add(myResult);
		}
	}
}

void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::mult(DataVector& alpha, DataVector& result)
{
	// distribute the current grid coefficients
	if (myGlobalMPIComm->getMyRank() == 0)
	{
		myGlobalMPIComm->broadcastGridCoefficients(alpha);
	}
	else
	{
		myGlobalMPIComm->receiveGridCoefficients(alpha);
	}

	this->applyLOperatorInner(alpha, result);

	// aggregate all results
	if (myGlobalMPIComm->getMyRank() == 0)
	{
		myGlobalMPIComm->aggregateGridCoefficients(result);
	}
	else
	{
		myGlobalMPIComm->sendGridCoefficients(result, 0);
	}
}

DataVector* PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::generateRHS()
{
	if (this->InnerGrid != NULL)
	{
		DataVector alpha_tmp_complete(*(this->rhs));
		DataVector rhs_tmp_complete(*(this->rhs));

		// distribute the current grid coefficients
		if (myGlobalMPIComm->getMyRank() == 0)
		{
			this->BoundaryUpdate->setInnerPointsToZero(alpha_tmp_complete);
			myGlobalMPIComm->broadcastGridCoefficients(alpha_tmp_complete);
		}
		else
		{
			myGlobalMPIComm->receiveGridCoefficients(alpha_tmp_complete);
		}

		applyLOperatorComplete(alpha_tmp_complete, rhs_tmp_complete);

		// aggregate all results
		if (myGlobalMPIComm->getMyRank() == 0)
		{
			myGlobalMPIComm->aggregateGridCoefficients(rhs_tmp_complete);
		}
		else
		{
			myGlobalMPIComm->sendGridCoefficients(rhs_tmp_complete, 0);
		}

		this->GridConverter->calcInnerCoefs(rhs_tmp_complete, *(this->rhs_inner));
		this->rhs_inner->mult(-1.0);
	}
	else
	{
		myGlobalMPIComm->Abort();
		throw new algorithm_exception("OperationEllipticPDESolverSystemDirichlet::generateRHS : No inner grid exists!");
	}

	return this->rhs_inner;
}

}
