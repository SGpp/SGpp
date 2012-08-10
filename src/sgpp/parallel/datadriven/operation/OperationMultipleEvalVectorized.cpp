/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp"
#include "base/exception/operation_exception.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif

namespace sg
{
namespace parallel
{

void OperationMultipleEvalVectorized::adaptDatasetBoundaries()
{
	//debugMPI(sg::parallel::myGlobalMPIComm, "passed the following bounds: grid:" << m_gridFrom << " - " << m_gridTo << "; dataset: " << m_datasetFrom << " - " << m_datasetTo)

//	//check for valid sized dataset and upper and lower border
//	if ( this->dataset_->getNcols() % CHUNKDATAPOINTS_X86 != 0 )
//	{
//		throw sg::base::operation_exception("For iterative mult an even number of dataset instances is required!");
//	}
//	if(m_datasetFrom%CHUNKDATAPOINTS_X86 != 0) {
//		throw sg::base::operation_exception("For iterative mult an even number of dataset instances is required and the borders have to be aligned to a multiple of the ");
//	}
//	if(m_datasetTo%CHUNKDATAPOINTS_X86 != 0) {
//		throw sg::base::operation_exception("For iterative mult an even number of dataset instances is required!");
//	}



//	round down to previous CHUNKDATAPOINTS_X86 border
//		if(m_datasetFrom%CHUNKDATAPOINTS_X86 != 0) {
//			int remainder = m_datasetFrom%CHUNKDATAPOINTS_X86;
//			m_datasetFrom -= remainder;
//		}

//	//round up to next CHUNKDATAPOINTS_X86 border
//	if(m_datasetTo%CHUNKDATAPOINTS_X86 != 0) {
//		int remainder = m_datasetTo%CHUNKDATAPOINTS_X86;
//		m_datasetTo += CHUNKDATAPOINTS_X86-remainder;
//	}

	//debugMPI(sg::parallel::myGlobalMPIComm, "doing calculations with the following dataset bounds: " << m_datasetFrom << " - " << m_datasetTo);

	// now that both m_datasetFrom and m_datasetTo are aligned to multiples of CHUNKDATAPOINTS_X86, also the
	// chunksize for this process is a multiple of CHUNKDATAPOINTS_X86
}

}
}
