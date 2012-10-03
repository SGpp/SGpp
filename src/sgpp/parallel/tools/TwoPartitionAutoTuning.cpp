/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/tools/TwoPartitionAutoTuning.hpp"
#include <iostream>
#include <algorithm>

#define INITIAL_SPEEDUP_PARTITION_2 10.0;

namespace sg
{
namespace parallel
{

TwoPartitionAutoTuning::TwoPartitionAutoTuning(size_t problemSize, size_t partition2Divider, size_t retune_cycles, double damping, double maxPercent)
	: _problemSize(problemSize), _partition2Divider(partition2Divider), _timePartition1(0.0), _timePartition2(0.0),
	  _oldSizePartition1(problemSize), _isFirstTuning(true),
	  _tuneCounter(0), _retune(retune_cycles), _isStatic(false),
	  _percentPartion1(0.0), _staticOutputCounter(0), _staticOutputFreq(20)
{
}

TwoPartitionAutoTuning::TwoPartitionAutoTuning(size_t problemSize, double percentPartion1, size_t partition2Divider, size_t OutputFreq)
	: _problemSize(problemSize), _partition2Divider(partition2Divider), _timePartition1(0.0), _timePartition2(0.0),
	  _oldSizePartition1(problemSize), _isFirstTuning(false),
	  _tuneCounter(0), _retune(50), _isStatic(true),
	  _percentPartion1(percentPartion1), _staticOutputCounter(0), _staticOutputFreq(OutputFreq)
{
	rescaleAutoTuning(_problemSize);
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
	if (!_isStatic)
	{
		size_t partition1 = 0;
		size_t partition2 = 0;

		// don't use an accelerator for smallest problems
		if (4*_partition2Divider > _problemSize)
		{
			partition1 = _problemSize;
			_oldSizePartition1 = _problemSize;
			partition2 = 0;
		}
		else
		{
			double partition2_speedup;
			double timeSave1 = _timePartition1;
			double timeSave2 = _timePartition2;

			if (((_tuneCounter % _retune) == 0) || _isFirstTuning == true)
			{
				if ( _isFirstTuning == false)
				{
					double partition1_element_time = static_cast<double>(_oldSizePartition1)/static_cast<double>(_timePartition1);
					double partition2_element_time = static_cast<double>(_problemSize-_oldSizePartition1)/static_cast<double>(_timePartition2);

					partition2_speedup = partition2_element_time/partition1_element_time;

					_timePartition1 = 0.0;
					_timePartition2 = 0.0;
					_tuneCounter = 0;
				}
				else
				{
					_isFirstTuning = false;
					partition2_speedup = INITIAL_SPEEDUP_PARTITION_2;

					_timePartition1 = 0.0;
					_timePartition2 = 0.0;
					_tuneCounter = 0;
				}

				double normalized_workingset = static_cast<double>(_problemSize)/(partition2_speedup+0.93);

				partition2 = static_cast<size_t>(normalized_workingset*partition2_speedup);

				size_t partition2_remainder = partition2 % _partition2Divider;
				if (partition2 + (_partition2Divider - partition2_remainder) > _problemSize)
				{
					partition2 -=  partition2_remainder;
				}
				else
				{
					partition2 +=  (_partition2Divider - partition2_remainder);
				}
				partition1 = _problemSize - partition2;

				_oldSizePartition1 = partition1;

				std::cout << "AUTOTUNING-PARTITION-SIZES (" << _problemSize << "): Time1: " << timeSave1 << " Size1: " << _oldSizePartition1 << "(" << 100.0*(double)_oldSizePartition1/(double)_problemSize << "%); Time2: " << timeSave2 << " Size2: " << _problemSize-_oldSizePartition1 << " (" << 100.0*(double)(_problemSize-_oldSizePartition1)/(double)_problemSize << "%)" << std::endl;
			}
		}
	}
	else
	{
		if (_staticOutputCounter % _staticOutputFreq == 0)
			std::cout << "AUTOTUNING-PARTITION-SIZES (" << _problemSize << "): Time1: " << _timePartition1 << " Size1: " << _oldSizePartition1 << "(" << 100.0*(double)_oldSizePartition1/(double)_problemSize << "%); Time2: " << _timePartition2 << " Size2: " << _problemSize-_oldSizePartition1 << " (" << 100.0*(double)(_problemSize-_oldSizePartition1)/(double)_problemSize << "%)" << std::endl;

		_staticOutputCounter++;
	}

	return _oldSizePartition1;
}

void TwoPartitionAutoTuning::setProblemSize(size_t problemSize)
{
	rescaleAutoTuning(problemSize);
	_problemSize = problemSize;
}

void TwoPartitionAutoTuning::setPartition2Divider(size_t partition2Divider)
{
	_partition2Divider = partition2Divider;
	rescaleAutoTuning(_problemSize);
}

void TwoPartitionAutoTuning::setExecutionTimes(double timePartition1, double timePartition2)
{
	_timePartition1 += timePartition1;
	_timePartition2 += timePartition2;
	_tuneCounter++;
}

void TwoPartitionAutoTuning::resetAutoTuning()
{
	_timePartition1 = 0.0;
	_timePartition2 = 0.0;
	_isFirstTuning = true;
	_tuneCounter = 0;
	_staticOutputCounter = 0;
}

void TwoPartitionAutoTuning::softResetAutoTuning()
{
	_timePartition1 = 0.0;
	_timePartition2 = 0.0;
	_tuneCounter = 0;
	_staticOutputCounter = 0;
}

void TwoPartitionAutoTuning::rescaleAutoTuning(size_t newProblemSize)
{
	if (newProblemSize == _problemSize)
	{
		return;
	}

	if (!_isStatic)
	{
		_problemSize = newProblemSize;
		_oldSizePartition1 = _problemSize;
		_isFirstTuning = true;
		_timePartition1 = 0.0;
		_timePartition2 = 0.0;
		_tuneCounter = 0;
	}
	else
	{
		size_t partition1 = 0;
		size_t partition2 = 0;

		partition1 = (size_t)(((double)newProblemSize)*_percentPartion1);

		partition1 = (size_t)std::min<double>(((double)newProblemSize)*_percentPartion1, (double)newProblemSize);
		partition2 = newProblemSize - partition1;

		size_t partition2_remainder = partition2 % _partition2Divider;
		partition2 -=  partition2_remainder;
		partition1 = newProblemSize - partition2;

		_problemSize = newProblemSize;
		_oldSizePartition1 = partition1;
		_staticOutputCounter = 0;
	}
}

}

}
