/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "parallel/tools/MPI/SGppMPITools.hpp"
#include "base/exception/operation_exception.hpp"

#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityMPI.hpp"
#include "parallel/operation/ParallelOpFactory.hpp"
#include "parallel/tools/PartitioningTool.hpp"

namespace sg
{
namespace parallel
{

DMSystemMatrixVectorizedIdentityMPI::DMSystemMatrixVectorizedIdentityMPI(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode)
	: DMSystemMatrixBase(trainData, lambda), vecMode_(vecMode), vecWidth_(0), numTrainingInstances_(0), numPatchedTrainingInstances_(0), m_grid(SparseGrid)
{
	// handle unsupported vector extensions
	// @TODO (heinecke) refactor: better way to set vector width
	if (this->vecMode_ == X86SIMD)
	{
		this->vecWidth_ = 24;
	}
	else if (this->vecMode_ == OpenCL)
	{
		this->vecWidth_ = 128;
	}
	else if (this->vecMode_ == Hybrid_X86SIMD_OpenCL)
	{
		this->vecWidth_ = 128;
	}
	else if (this->vecMode_ == ArBB)
	{
		this->vecWidth_ = 16;
	}
	else if (this->vecMode_ == MIC)
	{
		this->vecWidth_ = 96;
	}
	else if (this->vecMode_ == Hybrid_X86SIMD_MIC)
	{
		this->vecWidth_ = 96;
	}
	else
	{
		throw new sg::base::operation_exception("DMSystemMatrixVectorizedIdentity : un-supported vector extensions!");
	}

	this->dataset_ = new sg::base::DataMatrix(trainData);

	this->numTrainingInstances_ = this->dataset_->getNrows();

    // Assure that data has a "even" number of instances -> padding might be needed
	size_t remainder = this->dataset_->getNrows() % this->vecWidth_;
	size_t loopCount = this->vecWidth_ - remainder;

	if (loopCount != this->vecWidth_)
	{
		sg::base::DataVector lastRow(this->dataset_->getNcols());
		for (size_t i = 0; i < loopCount; i++)
		{
			this->dataset_->getRow(this->dataset_->getNrows()-1, lastRow);
			this->dataset_->resize(this->dataset_->getNrows()+1);
			this->dataset_->setRow(this->dataset_->getNrows()-1, lastRow);
		}
	}

	this->numPatchedTrainingInstances_ = this->dataset_->getNrows();

	if (this->vecMode_ != OpenCL && this->vecMode_ != ArBB && this->vecMode_ != Hybrid_X86SIMD_OpenCL)
	{
		this->dataset_->transpose();
	}

    int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();
    int mpi_rank = sg::parallel::myGlobalMPIComm->getMyRank();

    // arrays for distribution settings
	_mpi_grid_sizes = new int[mpi_size];
	_mpi_grid_offsets = new int[mpi_size];
	_mpi_data_sizes= new int[mpi_size];
	_mpi_data_offsets = new int[mpi_size];

    // calculate distribution
	calcDistribution(m_grid.getStorage()->size(), _mpi_grid_sizes, _mpi_grid_offsets, 1);
	calcDistribution(this->numPatchedTrainingInstances_, _mpi_data_sizes, _mpi_data_offsets, this->vecWidth_);

    debugMPI(sg::parallel::myGlobalMPIComm, "storage: " << _mpi_grid_offsets[mpi_rank] << " -- " << _mpi_grid_offsets[mpi_rank] + _mpi_grid_sizes[mpi_rank] - 1 << "size: " <<  _mpi_grid_sizes[mpi_rank]);
    debugMPI(sg::parallel::myGlobalMPIComm, "data:" << _mpi_data_offsets[mpi_rank]  << " -- " <<_mpi_data_offsets[mpi_rank] + _mpi_data_sizes[mpi_rank] - 1 << "size: " <<  _mpi_data_sizes[mpi_rank]);

    // mult: distribute calculations over storage
    // multTranspose: distribute calculations over dataset

    //std::cout << "gridtype: " << SparseGrid.getType() << std::endl;

	this->B_ = sg::op_factory::createOperationMultipleEvalVectorized(m_grid, this->vecMode_, this->dataset_,
                _mpi_grid_offsets[mpi_rank],
                _mpi_grid_offsets[mpi_rank] + _mpi_grid_sizes[mpi_rank],
                _mpi_data_offsets[mpi_rank],
                _mpi_data_offsets[mpi_rank] + _mpi_data_sizes[mpi_rank]
                );
}

DMSystemMatrixVectorizedIdentityMPI::~DMSystemMatrixVectorizedIdentityMPI()
{
	delete this->B_;
	delete this->dataset_;

    delete[] this->_mpi_grid_sizes;
    delete[] this->_mpi_grid_offsets;
    delete[] this->_mpi_data_sizes;
	delete[] this->_mpi_data_offsets;
}

void DMSystemMatrixVectorizedIdentityMPI::mult(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(this->numPatchedTrainingInstances_);

	// Operation B
	multVec(alpha, temp);

	// patch result -> set additional entries zero
	if (this->numTrainingInstances_ != temp.getSize())
	{
		for (size_t i = 0; i < (temp.getSize()-this->numTrainingInstances_); i++)
		{
			temp.set(temp.getSize()-(i+1), 0.0f);
		}
	}

	multTransposeVec(temp, result);

	result.axpy(static_cast<double>(this->numTrainingInstances_)*this->lambda_, alpha);
}

void DMSystemMatrixVectorizedIdentityMPI::generateb(sg::base::DataVector& classes, sg::base::DataVector& b)
{
	sg::base::DataVector myClasses(classes);

	// Apply padding
	if (this->numPatchedTrainingInstances_ != myClasses.getSize())
	{
		myClasses.resizeZero(this->numPatchedTrainingInstances_);
	}

	multTransposeVec(myClasses, b);

	debugMPI(sg::parallel::myGlobalMPIComm, "end generate b");
}

void DMSystemMatrixVectorizedIdentityMPI::rebuildLevelAndIndex()
{
    this->B_->rebuildLevelAndIndex();

	calcDistribution(m_grid.getStorage()->size(), _mpi_grid_sizes, _mpi_grid_offsets, 1);
	int mpi_rank = sg::parallel::myGlobalMPIComm->getMyRank();

	this->B_->updateGridComputeBoundaries(_mpi_grid_offsets[mpi_rank],
										  _mpi_grid_offsets[mpi_rank] + _mpi_grid_sizes[mpi_rank]);
}

void DMSystemMatrixVectorizedIdentityMPI::multVec(base::DataVector &alpha, base::DataVector &result)
{
	this->myTimer_->start();

	this->computeTimeMult_ += this->B_->multVectorized(alpha, result);

	debugMPI(sg::parallel::myGlobalMPIComm, "_mpi_data_offsets[" << sg::parallel::myGlobalMPIComm->getMyRank() << "] = " << _mpi_data_offsets[sg::parallel::myGlobalMPIComm->getMyRank()] << std::endl);
	debugMPI(sg::parallel::myGlobalMPIComm, "_mpi_data_sizes[" << sg::parallel::myGlobalMPIComm->getMyRank() << "] = " << _mpi_data_sizes[sg::parallel::myGlobalMPIComm->getMyRank()] << std::endl);

	debugMPI(sg::parallel::myGlobalMPIComm, "result.size() = " << result.getSize());

	sg::parallel::myGlobalMPIComm->dataVectorAllToAll(result, _mpi_data_offsets, _mpi_data_sizes);

	this->completeTimeMult_ += this->myTimer_->stop();
}

void DMSystemMatrixVectorizedIdentityMPI::multTransposeVec(base::DataVector &source, base::DataVector &result)
{
	this->myTimer_->start();
	this->computeTimeMultTrans_ += this->B_->multTransposeVectorized(source, result);

	sg::parallel::myGlobalMPIComm->dataVectorAllToAll(result, _mpi_grid_offsets, _mpi_grid_sizes);

	this->completeTimeMultTrans_ += this->myTimer_->stop();
}

void DMSystemMatrixVectorizedIdentityMPI::calcDistribution(int totalSize, int *sizes, int *offsets, size_t blocksize)
{
	for(int rank = 0; rank < sg::parallel::myGlobalMPIComm->getNumRanks(); ++rank){
		size_t size;
		size_t offset;
		sg::parallel::PartitioningTool::getPartitionSegment(totalSize, sg::parallel::myGlobalMPIComm->getNumRanks(), rank, &size, &offset, blocksize);
		sizes[rank] = size;
		offsets[rank] = offset;
	}
}

}
}
