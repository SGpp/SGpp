/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef TWOPARTIONAUTOTUNING_HPP
#define TWOPARTIONAUTOTUNING_HPP

#include <cstddef>

namespace sg
{
namespace parallel
{

/**
 * This class provides functionality to split a workload range
 * into two parts based on timing results from previous executions. The
 * second partition is forced to be dividable by a given divider.
 *
 * For instance this class is used in datadriven application for solving
 * systems of linear equations concurrently on hybrid platforms using standard
 * CPUs and accelerators.
 */
class TwoPartitionAutoTuning
{
public:
	/**
	 * Constructor for dynamic load balancing
	 *
	 * @param problemSize contains the overall size which should be partitioned
	 * @param partition2Divider the second partition divider, partition2's size is a multiple
	 * @param retune_cycles number of iteration after which the problem's separation is re-considered
	 */
	TwoPartitionAutoTuning(size_t problemSize, size_t partition2Divider, size_t retune_cycles);

	/**
	 * Constructor for static load balancing
	 *
	 * @param problemSize contains the overall size which should be partitioned
	 * @param percentPartion1 how big is the first, non accelerated portion?
	 * @param partition2Divider the second partition divider, partition2's size is a multiple
	 * @param OutputFreq how often should we print timings?
	 */
	TwoPartitionAutoTuning(size_t problemSize, double percentPartion1, size_t partition2Divider, size_t OutputFreq);

	/**
	 * Destructor
	 */
	virtual ~TwoPartitionAutoTuning();

	/**
	 * get problem size
	 *
	 * @return problem size to should partitioned
	 */
	size_t getProblemSize();

	/**
	 * gets size of partition 1 based on the currently stored
	 * runtimes for partition 1 and 2
	 *
	 * @return size of partition 1
	 */
	virtual size_t getPartition1Size();

	/**
	 * Update execution times in order to allow
	 * a new calculation of the partition sizes
	 *
	 * @param timePartition1 time needed for partition 1
	 * @param timePartition2 time needed for partition 2
	 */
	void setExecutionTimes(double timePartition1, double timePartition2);

	/**
	 * sets the problem size
	 *
	 * @param problemSize problem size
	 */
	void setProblemSize(size_t problemSize);

	/**
	 * set the possible divider of partition 2
	 *
	 * @param partition2Divider the divider of partition 2
	 */
	void setPartition2Divider(size_t partition2Divider);

	/**
	 * resets all auto tuning parameter
	 */
	void resetAutoTuning();

	/**
	 * resets only temp. auto tuning data
	 */
	void softResetAutoTuning();

protected:
	/// store problemsize
	size_t _problemSize;

	/// store required divider of partition 2
	size_t _partition2Divider;

	/// time needed to execute partition 1
	double _timePartition1;

	/// time needed to execute partition 2
	double _timePartition2;

	/// (old) size of partition1
	size_t _oldSizePartition1;

	/// first run, do initial calibration
	bool _isFirstTuning;

	/// counter for timer updates
	size_t _tuneCounter;

	/// number of updates that cause a tuning update
	size_t _retune;

	/// if static load balancing is enabled
	bool _isStatic;

	/// static, percent threshold of partition 1
	double _percentPartion1;

	/**
	 * re-scale the data and tuning parameter due to workload change
	 *
	 * @param newProblemSize new workload size
	 */
	void rescaleAutoTuning(size_t newProblemSize);
};

}

}

#endif /* TWOPARTIONAUTOTUNING_HPP */
