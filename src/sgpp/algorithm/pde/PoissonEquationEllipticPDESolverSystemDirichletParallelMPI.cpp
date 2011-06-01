/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "tools/MPI/SGppMPITools.hpp"

#include "algorithm/pde/PoissonEquationEllipticPDESolverSystemDirichletParallelMPI.hpp"
#include "exception/algorithm_exception.hpp"
#include "basis/operations_factory.hpp"

#include "algorithm/pde/StdUpDown.hpp"
#include "algorithm/pde/UpDownOneOpDim.hpp"
using namespace sg::pde;
using namespace sg::base;
using namespace sg::GridOperationFactory;

namespace sg
{
namespace parallel
{

PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::PoissonEquationEllipticPDESolverSystemDirichletParallelMPI(Grid& SparseGrid, DataVector& rhs) : OperationEllipticPDESolverSystemDirichlet(SparseGrid, rhs)
{
	this->Laplace_Complete = createOperationLaplace(*this->BoundGrid);
	this->Laplace_Inner = createOperationLaplace(*this->InnerGrid);
}

PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::~PoissonEquationEllipticPDESolverSystemDirichletParallelMPI()
{
	delete this->Laplace_Complete;
	delete this->Laplace_Inner;
}

void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::applyLOperatorInner(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	#pragma omp parallel shared(alpha, result)
	{
		#pragma omp single nowait
		{
			std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
			size_t nDims = algoDims.size();

			// Apply Laplace, parallel in Dimensions
			for (size_t i = 0; i < nDims; i++)
			{
				if (i % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank())
				{
					#pragma omp task firstprivate(i) shared(alpha, result, algoDims)
					{
						DataVector myResult(result.getSize());

						/// @todo (heinecke) discuss methods in order to avoid this cast
						((UpDownOneOpDim*)(this->Laplace_Inner))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

						#pragma omp critical
						{
							result.add(myResult);
						}
					}
				}
			}

			#pragma omp taskwait
		}
	}

}

void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::applyLOperatorComplete(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	#pragma omp parallel shared(alpha, result)
	{
		#pragma omp single nowait
		{
			std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
			size_t nDims = algoDims.size();

			// Apply Laplace, parallel in Dimensions
			for (size_t i = 0; i < nDims; i++)
			{
				if (i % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank())
				{
					#pragma omp task firstprivate(i) shared(alpha, result, algoDims)
					{

						DataVector myResult(result.getSize());

						/// @todo (heinecke) discuss methods in order to avoid this cast
						((UpDownOneOpDim*)(this->Laplace_Complete))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

						#pragma omp critical
						{
							result.add(myResult);
						}
					}
				}
			}

			#pragma omp taskwait
		}
	}
}

void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::mult(DataVector& alpha, DataVector& result)
{
	// distribute the current grid coefficients
	myGlobalMPIComm->broadcastGridCoefficientsFromRank0(alpha);

	this->applyLOperatorInner(alpha, result);

	// aggregate all results
	myGlobalMPIComm->reduceGridCoefficientsOnRank0(result);
}

DataVector* PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::generateRHS()
{
	if (this->InnerGrid != NULL)
	{
		DataVector alpha_tmp_complete(*(this->rhs));
		DataVector rhs_tmp_complete(*(this->rhs));

		this->BoundaryUpdate->setInnerPointsToZero(alpha_tmp_complete);
		// distribute the current grid coefficients
		myGlobalMPIComm->broadcastGridCoefficientsFromRank0(alpha_tmp_complete);

		applyLOperatorComplete(alpha_tmp_complete, rhs_tmp_complete);

		// aggregate all results
		myGlobalMPIComm->reduceGridCoefficientsOnRank0(rhs_tmp_complete);

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
}
