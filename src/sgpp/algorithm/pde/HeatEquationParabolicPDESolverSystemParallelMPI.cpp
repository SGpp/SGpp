/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "tools/MPI/SGppMPITools.hpp"

#include "algorithm/pde/HeatEquationParabolicPDESolverSystemParallelMPI.hpp"
#include "exception/algorithm_exception.hpp"

#include "algorithm/pde/StdUpDown.hpp"
#include "algorithm/pde/UpDownOneOpDim.hpp"

namespace sg
{

HeatEquationParabolicPDESolverSystemParallelMPI::HeatEquationParabolicPDESolverSystemParallelMPI(Grid& SparseGrid, DataVector& alpha, double a, double TimestepSize, std::string OperationMode)
{
	this->a = a;
	this->tOperationMode = OperationMode;
	this->TimestepSize = TimestepSize;
	this->BoundGrid = &SparseGrid;
	this->alpha_complete = &alpha;
	this->InnerGrid = NULL;
	this->alpha_inner = NULL;

	this->BoundaryUpdate = new DirichletUpdateVector(SparseGrid.getStorage());
	this->GridConverter = new DirichletGridConverter();

	this->OpLaplaceBound = SparseGrid.createOperationLaplace();
	this->OpMassBound = SparseGrid.createOperationLTwoDotProduct();

	// create the inner grid
	this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);

	//Create needed operations, on inner grid
	this->OpLaplaceInner = this->InnerGrid->createOperationLaplace();
	this->OpMassInner = this->InnerGrid->createOperationLTwoDotProduct();

	// right hand side if System
	this->rhs = NULL;
}

HeatEquationParabolicPDESolverSystemParallelMPI::~HeatEquationParabolicPDESolverSystemParallelMPI()
{
	delete this->OpLaplaceBound;
	delete this->OpMassBound;
	delete this->OpLaplaceInner;
	delete this->OpMassInner;

	delete this->BoundaryUpdate;
	delete this->GridConverter;
	if (this->InnerGrid != NULL)
	{
		delete this->InnerGrid;
	}
	if (this->alpha_inner != NULL)
	{
		delete this->alpha_inner;
	}
	if (this->rhs != NULL)
	{
		delete this->rhs;
	}
}

void HeatEquationParabolicPDESolverSystemParallelMPI::applyMassMatrixComplete(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);
	size_t nDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions().size();

	if (nDims % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank())
	{
		DataVector temp(alpha.getSize());

		((StdUpDown*)(this->OpMassBound))->multParallelBuildingBlock(alpha, temp);

		result.add(temp);
	}
}

void HeatEquationParabolicPDESolverSystemParallelMPI::applyLOperatorComplete(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	DataVector temp(alpha.getSize());
	temp.setAll(0.0);

	std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
	size_t nDims = algoDims.size();

	// Apply Laplace, parallel in Dimensions
	for (size_t i = 0; i < nDims; i++)
	{
		if (i % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank())
		{
			#pragma omp task firstprivate(i) shared(alpha, temp, result, algoDims)
			{
				DataVector myResult(result.getSize());

				/// @todo (heinecke) discuss methods in order to avoid this cast
				((UpDownOneOpDim*)(this->OpLaplaceBound))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

				// semaphore
				#pragma omp critical
				{
					temp.add(myResult);
				}
			}
		}
	}

	#pragma omp taskwait

	result.axpy((-1.0)*this->a,temp);
}

void HeatEquationParabolicPDESolverSystemParallelMPI::applyMassMatrixInner(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);
	size_t nDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions().size();

	if (nDims % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank())
	{
		DataVector temp(alpha.getSize());

		((StdUpDown*)(this->OpMassInner))->multParallelBuildingBlock(alpha, temp);

		result.add(temp);
	}
}

void HeatEquationParabolicPDESolverSystemParallelMPI::applyLOperatorInner(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	DataVector temp(alpha.getSize());
	temp.setAll(0.0);

	std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
	size_t nDims = algoDims.size();

	// Apply Laplace, parallel in Dimensions
	for (size_t i = 0; i < nDims; i++)
	{
		if (i % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank())
		{
			#pragma omp task firstprivate(i) shared(alpha, temp, result, algoDims)
			{
				DataVector myResult(result.getSize());

				/// @todo (heinecke) discuss methods in order to avoid this cast
				((UpDownOneOpDim*)(this->OpLaplaceInner))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

				// semaphore
				#pragma omp critical
				{
					temp.add(myResult);
				}
			}
		}
	}

	#pragma omp taskwait

	result.axpy((-1.0)*this->a,temp);
}

void HeatEquationParabolicPDESolverSystemParallelMPI::finishTimestep(bool isLastTimestep)
{
	// Replace the inner coefficients on the boundary grid
	this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);
}

void HeatEquationParabolicPDESolverSystemParallelMPI::startTimestep()
{
}

void HeatEquationParabolicPDESolverSystemParallelMPI::mult(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	if (this->tOperationMode == "ExEul")
	{
		applyMassMatrixInner(alpha, result);
	}
	else if (this->tOperationMode == "ImEul")
	{
		DataVector temp(result.getSize());
		DataVector temp2(result.getSize());

		#pragma omp parallel shared(alpha, result, temp, temp2)
		{
			#pragma omp single nowait
			{
				#pragma omp task shared (alpha, temp)
				{
					applyMassMatrixInner(alpha, temp);
				}

				#pragma omp task shared (alpha, temp2)
				{
					applyLOperatorInner(alpha, temp2);
				}

				#pragma omp taskwait
			}
		}

		result.add(temp);
		result.axpy((-1.0)*this->TimestepSize, temp2);
	}
	else if (this->tOperationMode == "CrNic")
	{
		DataVector temp(result.getSize());
		DataVector temp2(result.getSize());

		#pragma omp parallel shared(alpha, result, temp, temp2)
		{
			#pragma omp single nowait
			{
				#pragma omp task shared (alpha, temp)
				{
					applyMassMatrixInner(alpha, temp);
				}

				#pragma omp task shared (alpha, temp2)
				{
					applyLOperatorInner(alpha, temp2);
				}

				#pragma omp taskwait
			}
		}

		result.add(temp);
		result.axpy((-0.5)*this->TimestepSize, temp2);
	}
	else
	{
		throw new algorithm_exception(" HeatEquationParabolicPDESolverSystemParallelOMP::mult : An unknown operation mode was specified!");
	}
}

DataVector* HeatEquationParabolicPDESolverSystemParallelMPI::generateRHS()
{
	DataVector rhs_complete(this->alpha_complete->getSize());

	if (this->tOperationMode == "ExEul")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(rhs_complete.getSize());
		DataVector temp2(rhs_complete.getSize());
		DataVector myAlpha(*this->alpha_complete);

		#pragma omp parallel shared(myAlpha, temp, temp2)
		{
			#pragma omp single nowait
			{
				#pragma omp task shared (myAlpha, temp)
				{
					applyMassMatrixComplete(myAlpha, temp);
				}

				#pragma omp task shared (myAlpha, temp2)
				{
					applyLOperatorComplete(myAlpha, temp2);
				}

				#pragma omp taskwait
			}
		}

		rhs_complete.add(temp);
		rhs_complete.axpy(this->TimestepSize, temp2);
	}
	else if (this->tOperationMode == "ImEul")
	{
		rhs_complete.setAll(0.0);

		applyMassMatrixComplete(*this->alpha_complete, rhs_complete);
	}
	else if (this->tOperationMode == "CrNic")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(rhs_complete.getSize());
		DataVector temp2(rhs_complete.getSize());
		DataVector myAlpha(*this->alpha_complete);

		#pragma omp parallel shared(myAlpha, temp, temp2)
		{
			#pragma omp single nowait
			{
				#pragma omp task shared (myAlpha, temp)
				{
					applyMassMatrixComplete(myAlpha, temp);
				}

				#pragma omp task shared (myAlpha, temp2)
				{
					applyLOperatorComplete(myAlpha, temp2);
				}

				#pragma omp taskwait
			}
		}

		rhs_complete.add(temp);
		rhs_complete.axpy((0.5)*this->TimestepSize, temp2);
	}
	else
	{
		throw new algorithm_exception("HeatEquationParabolicPDESolverSystemParallelOMP::generateRHS : An unknown operation mode was specified!");
	}

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
		DataVector temp2(alpha_bound.getSize());

		#pragma omp parallel shared(alpha_bound, temp, temp2)
		{
			#pragma omp single nowait
			{
				#pragma omp task shared (alpha_bound, temp)
				{
					applyMassMatrixComplete(alpha_bound, temp);
				}

				#pragma omp task shared (alpha_bound, temp2)
				{
					applyLOperatorComplete(alpha_bound, temp2);
				}

				#pragma omp taskwait
			}
		}

		result_complete.add(temp);
		result_complete.axpy((-1.0)*this->TimestepSize, temp2);
	}
	else if (this->tOperationMode == "CrNic")
	{
		DataVector temp(alpha_bound.getSize());
		DataVector temp2(alpha_bound.getSize());

		#pragma omp parallel shared(alpha_bound, temp, temp2)
		{
			#pragma omp single nowait
			{
				#pragma omp task shared (alpha_bound, temp)
				{
					applyMassMatrixComplete(alpha_bound, temp);
				}

				#pragma omp task shared (alpha_bound, temp2)
				{
					applyLOperatorComplete(alpha_bound, temp2);
				}

				#pragma omp taskwait
			}
		}

		result_complete.add(temp);
		result_complete.axpy((-0.5)*this->TimestepSize, temp2);
	}
	else
	{
		throw new algorithm_exception("HeatEquationParabolicPDESolverSystemParallelOMP::generateRHS : An unknown operation mode was specified!");
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

}
