/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "tools/common/TwoPartitionAutoTuning.hpp""
#include <iostream>

namespace sg
{
namespace base
{

TwoPartitionAutoTuning::TwoPartitionAutoTuning(size_t problemSize, size_t partition2Divider)
	: _problemSize(problemSize), _partition2Divider(partition2Divider), _timePartition1(1.0), _timePartition2(0.0), _oldSizePartition1(problemSize)
{
}

TwoPartitionAutoTuning::~TwoPartitionAutoTuning()
{
}

size_t TwoPartitionAutoTuning::getProblemSize()
{
	return _problemSize;
}

size_t TwoPartitionAutoTuning::getPartition1Size()
{
	size_t partition1 = 0;
	size_t partition2 = 0;

	// don't use an accelerate for smallest problems
	if (2*_partition2Divider > _problemSize)
	{
		partition1 = _problemSize;
	}
	else
	{
		if (_timePartition2 == 0.0)
		{
			partition1 = partition2 = _problemSize/2;
			// ensure that whole problem is partitioned
			partition1 = _problemSize - partition2;
		}
		else
		{
			partition1 = _oldSizePartition1;
			partition2 = _problemSize - _oldSizePartition1;

			// very unlikely ... do nothing since optimal
			if (_timePartition1 == _timePartition2)
			{
			}
			else
			{
				double factor;

				if (_timePartition1 < _timePartition2)
				{
					factor = ((double)_timePartition2)/((double)_timePartition1);
				}
				else
				{
					factor = ((double)_timePartition1)/((double)_timePartition2);
				}

				partition1 = (size_t)(((double)partition1)*factor);
				partition2 = _problemSize - partition1;
			}
		}

		// Check if partition2 is dividable
	    size_t partition2_remainder = partition2 % _partition2Divider;
	    partition2 -= partition2_remainder;
	    partition1 = _problemSize - partition2;

	    std::cout <<  "AUTOTUNING-PARTITION-SIZES (" << _problemSize << "): Size1: " << partition1 << ", Size2: " << partition2 << std::endl;
	}

	_oldSizePartition1 = partition1;

	return partition1;
}

void TwoPartitionAutoTuning::setProblemSize(size_t problemSize)
{
	_problemSize = problemSize;
}

void TwoPartitionAutoTuning::setPartition2Divider(size_t partition2Divider)
{
	_partition2Divider = partition2Divider;
}

void TwoPartitionAutoTuning::setExecutionTimes(double timePartition1, double timePartition2)
{
	_timePartition1 = timePartition1;
	_timePartition2 = timePartition2;
}

}

}
