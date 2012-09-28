/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef DMSYSTEMMATRIXVECTORIZEDIDENTITYONESIDEDMPI_HPP
#define DMSYSTEMMATRIXVECTORIZEDIDENTITYONESIDEDMPI_HPP

#include "parallel/tools/MPI/SGppMPITools.hpp"

#include "base/datatypes/DataVector.hpp"
#include "base/grid/Grid.hpp"

#include "datadriven/algorithm/DMSystemMatrixBase.hpp"

#include "parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinearMult.h"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinearMultTranspose.h"
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
class DMSystemMatrixVectorizedIdentityOneSidedMPI : public sg::datadriven::DMSystemMatrixBase
{
private:
	/// vectorization mode
	VectorizationType vecMode_;
	/// Number of orignal training instances
	size_t numTrainingInstances_;
	/// Number of patched and used training instances
	size_t numPatchedTrainingInstances_;

	sg::base::DataMatrix* level_;
	/// Member to store the sparse grid's indices for better vectorization
	sg::base::DataMatrix* index_;


public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param trainData reference to sg::base::DataMatrix that contains the training data
	 * @param lambda the lambda, the regression parameter
	 * @param vecMode vectorization mode
	 */
	DMSystemMatrixVectorizedIdentityOneSidedMPI(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode);

	/**
	 * Std-Destructor
	 */
	virtual ~DMSystemMatrixVectorizedIdentityOneSidedMPI();

	virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

	virtual void generateb(sg::base::DataVector& classes, sg::base::DataVector& b);

	virtual void rebuildLevelAndIndex();

private:
	/// how to distribute storage array across processes
	int * _mpi_grid_sizes;
	int * _mpi_grid_offsets;

	/// reference to grid. needed to get new grid size after it changes
	sg::base::Grid& m_grid;

	/// how to distribute dataset across processes
	int * _mpi_data_sizes;
	int * _mpi_data_offsets;

	/// which chunks belong to which process
	int * _mpi_data_sizes_global;
	int * _mpi_data_offsets_global;

	/// which chunks belong to which process
	int * _mpi_grid_sizes_global;
	int * _mpi_grid_offsets_global;

	/// into how many chunks should data and grid be partitioned
	size_t _chunkCountData;
	size_t _chunkCountGrid;

	/// MPI windows
	sg::base::DataVector* _mpi_grid_window_buffer;
	sg::base::DataVector* _mpi_data_window_buffer;
	MPI_Win* _mpi_grid_window;
	MPI_Win* _mpi_data_window;

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
	void calcDistribution(int totalSize, int numChunks, int* sizes, int* offsets, size_t blocksize);
};

}
}

#endif /* DMSYSTEMMATRIXVECTORIZEDIDENTITYONESIDEDMPI_HPP */
