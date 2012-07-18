/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "parallel/tools/MPI/SGppMPITools.hpp"
#include "parallel/datadriven/operation/OperationMultipleEvalVectorizedSP.hpp"
#include "base/exception/operation_exception.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif

namespace sg
{
namespace parallel
{

void OperationMultipleEvalVectorizedSP::calcOpenMPLoopDistribution(int processStart, int processEnd, int blockSize, size_t *start, size_t *end)
{
	// check for valid input
	size_t chunkSizeProc = processEnd - processStart;
	if (chunkSizeProc % blockSize != 0 )
	{
		std::cout << "chunkSize: " << blockSize << "; chunkSizeProc: " << chunkSizeProc << std::endl;
		throw sg::base::operation_exception("processed vector segment must fit to chunkSize, but it does not!");
	}

	// doing further calculations with complete blocks
	size_t blockCount = chunkSizeProc/blockSize;
	int blockFragmentSize = blockCount;
	int blockFragmentOffset=0;

	int threadsCount = 1;
	int myThreadNum = 0;
#ifdef _OPENMP
	threadsCount = omp_get_num_threads();
	myThreadNum = omp_get_thread_num();
#endif
	//calculate distribution of blocks
	sg::parallel::myGlobalMPIComm->calcDistributionFragment(blockCount, threadsCount, myThreadNum, &blockFragmentSize, &blockFragmentOffset);

	// set return values
	*start = processStart + blockFragmentOffset*blockSize;
	*end = *start+blockFragmentSize*blockSize;
}

void OperationMultipleEvalVectorizedSP::adaptDatasetBoundaries()
{
	debugMPI(sg::parallel::myGlobalMPIComm, "passed the following bounds: grid:" << m_gridFrom << " - " << m_gridTo << "; dataset: " << m_datasetFrom << " - " << m_datasetTo)

	//check for valid sized dataset already here
	if ( this->dataset_->getNcols() % CHUNKDATAPOINTS_SP_X86 != 0 )
	{
		throw sg::base::operation_exception("For iterative mult transpose an even number of instances is required!");
	}

	//round down to previous CHUNKDATAPOINTS_X86 border
	if(m_datasetFrom%CHUNKDATAPOINTS_SP_X86 != 0) {
		int remainder = m_datasetFrom%CHUNKDATAPOINTS_SP_X86;
		m_datasetFrom -= remainder;
	}

	//round up to next CHUNKDATAPOINTS_X86 border
	if(m_datasetTo%CHUNKDATAPOINTS_SP_X86 != 0) {
		int remainder = m_datasetTo%CHUNKDATAPOINTS_SP_X86;
		m_datasetTo += CHUNKDATAPOINTS_SP_X86-remainder;
	}

	debugMPI(sg::parallel::myGlobalMPIComm, "doing calculations with the following dataset bounds: " << m_datasetFrom << " - " << m_datasetTo);

	// now that both m_datasetFrom and m_datasetTo are aligned to multiples of CHUNKDATAPOINTS_X86, also the
	// chunksize for this process is a multiple of CHUNKDATAPOINTS_X86
}

}
}
