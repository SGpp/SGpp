/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "parallel/tools/MPI/SGppMPITools.hpp"
#include "parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp"
#include "base/exception/operation_exception.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif

namespace sg
{
namespace parallel
{

void OperationMultipleEvalVectorized::calcOpenMPLoopDistribution(int processStart, int processEnd, int chunkSize, size_t *start, size_t *end)
{
#ifdef _OPENMP
	size_t chunkSizeProc = processEnd - processStart;
	if (chunkSizeProc % chunkSize != 0 )
	{
		std::cout << "chunkSize: " << chunkSize << "; chunkSizeProc: " << chunkSizeProc << std::endl;
		throw sg::base::operation_exception("processed vector segment must fit to chunkSize, but it does not!");
	}
	//doing further calculations with complete blocks
	size_t blockCount = chunkSizeProc/chunkSize;

	int chunkFragmentSize, chunkFragmentOffset;
	sg::parallel::myGlobalMPIComm->calcDistributionFragment(blockCount, omp_get_num_threads(), omp_get_thread_num(), &chunkFragmentSize, &chunkFragmentOffset);

	*start = processStart + chunkFragmentOffset*chunkSize;
	*end = *start+chunkFragmentSize*chunkSize;
	//std::cout << "[OpenMP Thread " << omp_get_thread_num() << "] [mult] start: " << start << "; end: " << end << std::endl;
#else
	*start = processStart;
	*end = processEnd;
#endif
}

void OperationMultipleEvalVectorized::adaptDatasetBoundaries()
{
	debugMPI(sg::parallel::myGlobalMPIComm, "passed the following bounds: grid:" << m_gridFrom << " - " << m_gridTo << "; dataset: " << m_datasetFrom << " - " << m_datasetTo)

	//check for valid sized dataset already here
	if ( this->dataset_->getNcols() % CHUNKDATAPOINTS_X86 != 0 )
	{
		throw sg::base::operation_exception("For iterative mult transpose an even number of instances is required!");
	}

	//round down to previous CHUNKDATAPOINTS_X86 border
	if(m_datasetFrom%CHUNKDATAPOINTS_X86 != 0) {
		int remainder = m_datasetFrom%CHUNKDATAPOINTS_X86;
		m_datasetFrom -= remainder;
	}

	//round up to next CHUNKDATAPOINTS_X86 border
	if(m_datasetTo%CHUNKDATAPOINTS_X86 != 0) {
		int remainder = m_datasetTo%CHUNKDATAPOINTS_X86;
		m_datasetTo += CHUNKDATAPOINTS_X86-remainder;
	}

	debugMPI(sg::parallel::myGlobalMPIComm, "doing calculations with the following dataset bounds: " << m_datasetFrom << " - " << m_datasetTo);

	// now that both m_datasetFrom and m_datasetTo are aligned to multiples of CHUNKDATAPOINTS_X86, also the
	// chunksize for this process is a multiple of CHUNKDATAPOINTS_X86
}

}
}
