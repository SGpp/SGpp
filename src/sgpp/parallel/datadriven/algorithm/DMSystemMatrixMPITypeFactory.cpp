/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "DMSystemMatrixMPITypeFactory.hpp"

#include "parallel/datadriven/basis/linear/noboundary/operation/impl/SPX86SimdLinear.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinear.hpp"

#include "parallel/datadriven/basis/modlinear/operation/impl/SPX86SimdModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/SPX86SimdModLinearMask.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/X86SimdModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/X86SimdModLinearMask.hpp"

#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityMPI.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityAsyncMPI.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityTrueAsyncMPI.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityTrueAsyncMPIAlltoallv.hpp"
//#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityOnesidedMPI.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityAllreduce.hpp"

#include <cstring>
#include "base/exception/factory_exception.hpp"

#include "parallel/tools/MPI/SGppMPITools.hpp"

namespace sg {
namespace parallel {

template<typename KernelImplementation>
datadriven::DMSystemMatrixBase *DMSystemMatrixMPITypeFactory::createDMSystemMatrixMPIType(base::Grid &grid, base::DataMatrix &trainDataset, double lambda, VectorizationType vecType)
{
#define MPI_TYPE_ALLTOALLV 1
#define MPI_TYPE_ALLREDUCE 2
#define MPI_TYPE_ASYNC 3
#define MPI_TYPE_TRUEASYNC 4
#define MPI_TYPE_TRUEASYNC_ALLTOALLV 5
#define MPI_TYPE_ONESIDED 6

	int mpi_type = MPI_TYPE_ASYNC;

	std::string parallelizationType;
	datadriven::DMSystemMatrixBase *result = 0;
	if(mpi_type == MPI_TYPE_ALLTOALLV){
		parallelizationType = "Alltoallv";
		result = new sg::parallel::DMSystemMatrixVectorizedIdentityMPI(grid, trainDataset, lambda, vecType);
	} else if(mpi_type == MPI_TYPE_ALLREDUCE) {
		parallelizationType = "Allreduce";
		result = new sg::parallel::DMSystemMatrixVectorizedIdentityAllreduce<KernelImplementation>(grid, trainDataset, lambda, vecType);
	} else if(mpi_type == MPI_TYPE_ASYNC) {
		parallelizationType = "Asynchronous Communication";
		result = new sg::parallel::DMSystemMatrixVectorizedIdentityAsyncMPI<KernelImplementation>(grid, trainDataset, lambda, vecType);
	} else if(mpi_type == MPI_TYPE_TRUEASYNC) {
		parallelizationType = "True Asynchronous Communication";
		result = new sg::parallel::DMSystemMatrixVectorizedIdentityTrueAsyncMPI<KernelImplementation>(grid, trainDataset, lambda, vecType);
	} else if(mpi_type == MPI_TYPE_TRUEASYNC_ALLTOALLV) {
		parallelizationType = "True Asynchronous Communication with alltoall end";
		result = new sg::parallel::DMSystemMatrixVectorizedIdentityTrueAsyncMPIAlltoallv<KernelImplementation>(grid, trainDataset, lambda, vecType);
	} /*else if(mpi_type == MPI_TYPE_ONESIDED) {
		parallelizationType = "Onesided Communication";
		result = new sg::parallel::DMSystemMatrixVectorizedIdentityOneSidedMPI<KernelImplementation>(grid, trainDataset, lambda, vecType);
	} */ else {
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
		return createDMSystemMatrixMPIType<X86SimdLinear>
				(grid, trainDataset, lambda, vecType);
	}
	else if(strcmp(grid.getType(), "modlinear") == 0)
	{
		return createDMSystemMatrixMPIType<X86SimdModLinear>
				(grid, trainDataset, lambda, vecType);
	}
	else
	{
		throw base::factory_exception("OperationMultipleEvalVectorizedSP is not implemented for this grid type.");
	}
}

}
}
