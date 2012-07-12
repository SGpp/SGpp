/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef DMSYSTEMMATRIXVECTORIZEDIDENTITYMPI_HPP
#define DMSYSTEMMATRIXVECTORIZEDIDENTITYMPI_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/grid/Grid.hpp"

#include "datadriven/algorithm/DMSystemMatrixBase.hpp"

#include "parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp"
#include "parallel/tools/TypesParallel.hpp"

#include <string>

namespace sg
{
namespace parallel
{

/**
 * Class that implements the virtual class sg::base::OperationMatrix for the
 * application of classification for the Systemmatrix
 *
 * The Identity matrix is used as regularization operator.
 *
 * For the Operation B's mult and mutlTransposed functions
 * vectorized formulations are used.
 */
class DMSystemMatrixVectorizedIdentityMPI : public sg::datadriven::DMSystemMatrixBase
{
private:
	/// vectorization mode
	VectorizationType vecMode_;
	/// vector width, class internal variable to enable padding and patching of vectors
	size_t vecWidth_;
	/// Number of orignal training instances
	size_t numTrainingInstances_;
	/// Number of patched and used training instances
	size_t numPatchedTrainingInstances_;
	/// OperationB for calculating the data matrix
	sg::parallel::OperationMultipleEvalVectorized* B_;

public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param trainData reference to sg::base::DataMatrix that contains the training data
	 * @param lambda the lambda, the regression parameter
	 * @param vecMode vectorization mode
	 */
    DMSystemMatrixVectorizedIdentityMPI(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode);

	/**
	 * Std-Destructor
	 */
    virtual ~DMSystemMatrixVectorizedIdentityMPI();

	virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

	virtual void generateb(sg::base::DataVector& classes, sg::base::DataVector& b);

	virtual void rebuildLevelAndIndex();

private:
    /// how to distribute storage array
    int* _mpi_storage_sizes;
    int* _mpi_storage_offsets;

    /// how to distribute grid
    int* _mpi_data_sizes;
    int* _mpi_data_offsets;

//    /// which part of the storage to send to all the other processes, will hold MPI_SIZE identical values
//    int* _mpi_storage_send_sizes;
//    int* _mpi_storage_send_offsets;

//    /// which part of the grid to send to all the other processes, will hold MPI_SIZE identical values
//    int* _mpi_data_send_sizes;
//    int* _mpi_data_send_offsets;

    /**
     * Wrapper function that handles communication after calculation and time measurement
     */
    void multVec(base::DataVector &alpha, base::DataVector &result);

    /**
     * Wrapper function that handles communication after calculation and time measurement
     */
    void multTransposeVec(base::DataVector &source, base::DataVector &result);

    /**
     * calculate size and offset of fragment number rank for a domain with a size of totalSize for a distribution over procCount Fragments
     * distribute domain as equal as possible
     * (the difference between the minimum and maximum number of items is at most 1)
     * otherwise, a worker could have almost twice as much to do as all the others (example: total=127, proccount = 16)
     *
     * @param totalSize size of domain to distribute
     * @param procCount number of fragments (processes) to distribute domain into
     * @param rank specifies the number of the fragment for which to calculate size and offset
     * @param size output variable to put resulting size into
     * @param offset output variable to put resulting offset into
     */
    void calcDistributionFragment(int totalSize, int procCount, int rank, int* size, int *offset);

    /**
     * calculates the distribution for the current MPI setting for a domain of
     * size totalSize and stores the result into the arrays sizes and offsets. These
     * arrays must have a size equal to the number of MPI processes currently running.
     *
     * @param totalSize size of domain to distribute
     * @param sizes output array to store resulting distribution sizes (array size must match the number of MPI processes)
     * @param offsets output array to store resulting distribution offsets (array size must match the number of MPI processes)
     *
     */
    void calcDistribution(int totalSize, int* sizes, int* offsets);
};

}
}

#endif /* DMSYSTEMMATRIXVECTORIZEDIDENTITYMPI_HPP */
