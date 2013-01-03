/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "DMSystemMatrixMPITypeFactory.hpp"

#include "parallel/datadriven/basis/linear/noboundary/operation/impl/SPX86SimdLinearMult.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/SPX86SimdLinearMultTranspose.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinearMult.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinearMultTranspose.hpp"

#include "parallel/datadriven/basis/modlinear/operation/impl/SPX86SimdModLinearMult.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/SPX86SimdModLinearMultTranspose.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/X86SimdModLinearMult.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/X86SimdModLinearMultTranspose.hpp"

#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityMPI.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityAsyncMPI.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityOnesidedMPI.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityAllreduce.hpp"

#include <cstring>
#include "base/exception/factory_exception.hpp"

#include "parallel/tools/MPI/SGppMPITools.hpp"

namespace sg {
namespace parallel {

template<typename MultType, typename MultTransType>
datadriven::DMSystemMatrixBase *DMSystemMatrixMPITypeFactory::createDMSystemMatrixMPIType(base::Grid &grid, base::DataMatrix &trainDataset, double lambda, VectorizationType vecType)
{
#define MPI_TYPE_ALLTOALLV 1
#define MPI_TYPE_ALLREDUCE 2
#define MPI_TYPE_ASYNC 3
#define MPI_TYPE_ONESIDED 4

	int mpi_type = MPI_TYPE_ASYNC;

	std::string parallelizationType;
	datadriven::DMSystemMatrixBase *result = 0;
	if(mpi_type == MPI_TYPE_ALLTOALLV){
		parallelizationType = "Alltoallv";
		result = new sg::parallel::DMSystemMatrixVectorizedIdentityMPI(grid, trainDataset, lambda, vecType);
	} else if(mpi_type == MPI_TYPE_ALLREDUCE) {
		parallelizationType = "Allreduce";
		result = new sg::parallel::DMSystemMatrixVectorizedIdentityAllreduce<MultType, MultTransType>(grid, trainDataset, lambda, vecType);
	} else if(mpi_type == MPI_TYPE_ASYNC) {
		parallelizationType = "Asynchronous Communication";
		result = new sg::parallel::DMSystemMatrixVectorizedIdentityAsyncMPI<MultType, MultTransType>(grid, trainDataset, lambda, vecType);
	} else if(mpi_type == MPI_TYPE_ONESIDED) {
		parallelizationType = "Onesided Communication";
		result = new sg::parallel::DMSystemMatrixVectorizedIdentityOneSidedMPI<MultType, MultTransType>(grid, trainDataset, lambda, vecType);
	} else {
		throw new sg::base::factory_exception("not implemented");
	}

	if(sg::parallel::myGlobalMPIComm->getMyRank() == 0){
		std::cout << "Using MPI Parallelization: " << parallelizationType << " ("<< sg::parallel::myGlobalMPIComm->getNumRanks() <<" Processes)" << std::endl;
		size_t thread_count = 1;
#ifdef _OPENMP
#pragma omp parallel
		{
			thread_count = omp_get_num_threads();
		}
#endif
		std::cout << "OpenMP: " << thread_count << " Threads active" << std::endl;
	}

	return result;
}

datadriven::DMSystemMatrixBase *DMSystemMatrixMPITypeFactory::getDMSystemMatrix(base::Grid &grid, base::DataMatrix &trainDataset, double lambda, VectorizationType vecType)
{
	if(strcmp(grid.getType(), "linear") == 0 || strcmp(grid.getType(), "linearBoundary") == 0
			|| strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
	{
		return createDMSystemMatrixMPIType<sg::parallel::X86SimdLinearMult, sg::parallel::X86SimdLinearMultTranspose>
				(grid, trainDataset, lambda, vecType);
	}
	else if(strcmp(grid.getType(), "modlinear") == 0)
	{
		return createDMSystemMatrixMPIType<sg::parallel::X86SimdModLinearMult, sg::parallel::X86SimdModLinearMultTranspose>
				(grid, trainDataset, lambda, vecType);
	}
	else
	{
		throw base::factory_exception("OperationMultipleEvalVectorizedSP is not implemented for this grid type.");
	}
}

}
}
