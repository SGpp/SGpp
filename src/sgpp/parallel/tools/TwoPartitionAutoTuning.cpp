/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/tools/TwoPartitionAutoTuning.hpp"
#include <iostream>
#include <algorithm>

namespace sg
{
namespace parallel
{

TwoPartitionAutoTuning::TwoPartitionAutoTuning(size_t problemSize, size_t partition2Divider, size_t retune_cycles, double damping, double maxPercent)
	: _problemSize(problemSize), _partition2Divider(partition2Divider), _timePartition1(1.0), _timePartition2(0.0),
	  _oldSizePartition1(problemSize), _testPartition1(true), _testPartition2(true), _isFirstTuning(true),
	  _tuneCounter(0), _retune(retune_cycles), _damping(damping), _maxPercent(maxPercent), _isStatic(false), _percentPartion1(0.0)
{
}

TwoPartitionAutoTuning::TwoPartitionAutoTuning(size_t problemSize, double percentPartion1, size_t partition2Divider)
	: _problemSize(problemSize), _partition2Divider(partition2Divider), _timePartition1(0.0), _timePartition2(0.0),
	  _oldSizePartition1(problemSize), _testPartition1(false), _testPartition2(false), _isFirstTuning(false),
	  _tuneCounter(0), _retune(50), _damping(0.0), _maxPercent(0.0), _isStatic(true), _percentPartion1(percentPartion1)
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

		// don't use an accelerate for smallest problems
		if (2*_partition2Divider > _problemSize)
		{
			partition1 = _problemSize;
			_oldSizePartition1 = _problemSize;
			partition2 = 0;
		}
		else
		{
			if (_testPartition1 == true)
			{
				partition1 = _problemSize;
				partition2 = 0;

				_oldSizePartition1 = partition1;
			}
			else if (_testPartition2 == true)
			{
				partition1 = 0;
				partition2 = _problemSize;
				size_t partition2_remainder = partition2 % _partition2Divider;
				partition2 -=  partition2_remainder;
				partition1 = _problemSize - partition2;

				_oldSizePartition1 = partition1;
			}
			else if (_isFirstTuning == true)
			{
				_isFirstTuning = false;

				double maxtime = std::max<double>(_timePartition1, _timePartition2);

				double part1 = maxtime / _timePartition1;
				double part2 = maxtime / _timePartition2;
				double parts = part1 + part2;
				double factor = _damping * part1 / parts;

				partition1 = (size_t)std::min<double>(((double)_problemSize)*factor, (double)_problemSize);
				partition2 = _problemSize - partition1;

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

				std::cout << "AUTOTUNING-PARTITION-SIZES (" << _problemSize << "): Time1: " << _timePartition1 << " Size1: " << _oldSizePartition1 << "(" << 100.0*(double)_oldSizePartition1/(double)_problemSize << "%); Time2: " << _timePartition2 << " Size2: " << _problemSize-_oldSizePartition1 << " (" << 100.0*(double)(_problemSize-_oldSizePartition1)/(double)_problemSize << "%)" << std::endl;

				_timePartition1 = 0.0;
				_timePartition2 = 0.0;
				_tuneCounter = 0;
			}
			else if ((_tuneCounter % _retune) == 0 && _tuneCounter != 0)
			{
				double factor = _damping * _timePartition2/_timePartition1;

				// only allow 3% change
				if (factor < (1.0-(_maxPercent/100.0)))
				{
					factor = 1.0-(_maxPercent/100.0);
				}
				if (factor > (1.0+(_maxPercent/100.0)))
				{
					factor = 1.0+(_maxPercent/100.0);
				}

				partition1 = (size_t)std::min<double>(((double)_oldSizePartition1)*factor, (double)_problemSize);
				partition2 = _problemSize - partition1;

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

				std::cout << "AUTOTUNING-PARTITION-SIZES (" << _problemSize << "): Time1: " << _timePartition1 << " Size1: " << _oldSizePartition1 << "(" << 100.0*(double)_oldSizePartition1/(double)_problemSize << "%); Time2: " << _timePartition2 << " Size2: " << _problemSize-_oldSizePartition1 << " (" << 100.0*(double)(_problemSize-_oldSizePartition1)/(double)_problemSize << "%)" << std::endl;

				_timePartition1 = 0.0;
				_timePartition2 = 0.0;
				_tuneCounter = 0;
			}
		}
	}
	else
	{
		std::cout << "AUTOTUNING-PARTITION-SIZES (" << _problemSize << "): Time1: " << _timePartition1 << " Size1: " << _oldSizePartition1 << "(" << 100.0*(double)_oldSizePartition1/(double)_problemSize << "%); Time2: " << _timePartition2 << " Size2: " << _problemSize-_oldSizePartition1 << " (" << 100.0*(double)(_problemSize-_oldSizePartition1)/(double)_problemSize << "%)" << std::endl;
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
	if (_testPartition1 == true)
	{
		_testPartition1 = false;
		_timePartition1 = timePartition1;
	}
	else if (_testPartition2 == true)
	{
		_testPartition2 = false;
		_timePartition2 = timePartition2;
	}
	else
	{
		_timePartition1 += timePartition1;
		_timePartition2 += timePartition2;
		_tuneCounter++;
	}
}

void TwoPartitionAutoTuning::resetAutoTuning()
{
	_testPartition1 = true;
	_testPartition2 = true;
	_isFirstTuning = true;
	_tuneCounter = 0;
}

void TwoPartitionAutoTuning::softResetAutoTuning()
{
	_timePartition1 = 0.0;
	_timePartition2 = 0.0;
	_tuneCounter = 0;
}

void TwoPartitionAutoTuning::rescaleAutoTuning(size_t newProblemSize)
{
	if (!_isStatic)
	{
		if (_testPartition1 == false && _testPartition2 == false && _isFirstTuning == false)
		{
			double factor = _damping*((double)_oldSizePartition1)/((double)_problemSize);

			size_t partition1 = 0;
			size_t partition2 = 0;

			partition1 = (size_t)(((double)newProblemSize)*factor);

			partition1 = (size_t)std::min<double>(((double)newProblemSize)*factor, (double)newProblemSize);
			partition2 = newProblemSize - partition1;

			size_t partition2_remainder = partition2 % _partition2Divider;
			partition2 -=  partition2_remainder;
			partition1 = newProblemSize - partition2;

			_problemSize = newProblemSize;
			_oldSizePartition1 = partition1;
		}
		else
		{
			_problemSize = newProblemSize;
			_oldSizePartition1 = _problemSize;
			_testPartition1 = true;
			_testPartition2 = true;
			_isFirstTuning = true;
		}
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
	}
}

}

}
