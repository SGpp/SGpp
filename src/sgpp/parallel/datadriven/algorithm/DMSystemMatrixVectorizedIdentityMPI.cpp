/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
#ifdef USE_MPI
#include <mpi.h>
#endif

#include "base/exception/operation_exception.hpp"

#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityMPI.hpp"
#include "parallel/operation/ParallelOpFactory.hpp"
#include "parallel/tools/MPI/SGppMPITools.hpp"

namespace sg
{
namespace parallel
{

DMSystemMatrixVectorizedIdentityMPI::DMSystemMatrixVectorizedIdentityMPI(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode)
	: DMSystemMatrixBase(trainData, lambda), vecMode_(vecMode), vecWidth_(0), numTrainingInstances_(0), numPatchedTrainingInstances_(0)
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

	if (this->vecMode_ != OpenCL && this->vecMode_ != ArBB  && this->vecMode_ != Hybrid_X86SIMD_OpenCL)
	{
		this->dataset_->transpose();
	}

	this->myTimer_ = new sg::base::SGppStopwatch();

    int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();
    int mpi_rank = sg::parallel::myGlobalMPIComm->getMyRank();

    // arrays for distribution settings
    _mpi_storage_sizes = new int[mpi_size];
    _mpi_storage_offsets = new int[mpi_size];
    _mpi_data_sizes= new int[mpi_size];
    _mpi_data_offsets = new int[mpi_size];
    _mpi_storage_send_sizes = new int[mpi_size];
    _mpi_storage_send_offsets = new int[mpi_size];
    _mpi_data_send_sizes = new int[mpi_size];
    _mpi_data_send_offsets = new int[mpi_size];

    // calculate distribution
    calcDistribution(SparseGrid.getStorage()->size(), _mpi_storage_sizes, _mpi_storage_offsets);
    calcDistribution(this->dataset_->getNcols(), _mpi_data_sizes, _mpi_data_offsets);

    // put values into send array
    for(int rank = 0; rank<mpi_size; rank++){
        _mpi_storage_send_sizes[rank] = _mpi_storage_sizes[mpi_rank];
        _mpi_storage_send_offsets[rank] = _mpi_storage_offsets[mpi_rank];
        _mpi_data_send_sizes[rank] = _mpi_data_sizes[mpi_rank];
        _mpi_data_send_offsets[rank] = _mpi_data_offsets[mpi_rank];
    }
    std::cout << "[rank " << mpi_rank << "] storage: " << _mpi_storage_offsets[mpi_rank] << " -- " << _mpi_storage_offsets[mpi_rank] + _mpi_storage_sizes[mpi_rank] - 1 << "size: " <<  _mpi_storage_sizes[mpi_rank] <<std::endl;
    std::cout << "[rank " << mpi_rank << "] data:" << _mpi_data_offsets[mpi_rank]  << " -- " <<_mpi_data_offsets[mpi_rank] + _mpi_data_sizes[mpi_rank] - 1 << "size: " <<  _mpi_data_sizes[mpi_rank]  << std::endl;

    std::cout << "my rank:  #" << mpi_rank << std::endl;

    // mult: distribute calculations over storage
    // multTranspose: distribute calculations over dataset

    //std::cout << "gridtype: " << SparseGrid.getType() << std::endl;

    this->B_ = sg::op_factory::createOperationMultipleEvalVectorized(SparseGrid, this->vecMode_, this->dataset_,
                _mpi_storage_offsets[mpi_rank],
                _mpi_storage_offsets[mpi_rank] + _mpi_storage_sizes[mpi_rank],
                _mpi_data_offsets[mpi_rank],
                _mpi_data_offsets[mpi_rank] + _mpi_data_sizes[mpi_rank]
                );
}

DMSystemMatrixVectorizedIdentityMPI::~DMSystemMatrixVectorizedIdentityMPI()
{
	delete this->B_;
	delete this->dataset_;

    delete this->_mpi_storage_sizes;
    delete this->_mpi_storage_offsets;
    delete this->_mpi_data_sizes;
    delete this->_mpi_data_offsets;
    delete this->_mpi_storage_send_sizes;
    delete this->_mpi_storage_send_offsets;
    delete this->_mpi_data_send_sizes;
    delete this->_mpi_data_send_offsets;
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

    std::cout << "end generate b"<< sg::parallel::myGlobalMPIComm->getMyRank() << std::endl;
}

void DMSystemMatrixVectorizedIdentityMPI::rebuildLevelAndIndex()
{
    this->B_->rebuildLevelAndIndex();
}

void DMSystemMatrixVectorizedIdentityMPI::multVec(base::DataVector &alpha, base::DataVector &result)
{
    this->myTimer_->start();
    this->computeTimeMult_ += this->B_->multVectorized(alpha, result);

    sg::parallel::myGlobalMPIComm->dataVectorAllToAll(result, _mpi_data_offsets, _mpi_data_sizes);

    this->completeTimeMult_ += this->myTimer_->stop();
}

void DMSystemMatrixVectorizedIdentityMPI::multTransposeVec(base::DataVector &source, base::DataVector &result)
{
    this->myTimer_->start();
    this->computeTimeMultTrans_ += this->B_->multTransposeVectorized(source, result);
    //std::cout << "size of result: " << result.getSize() << std::endl;

    sg::parallel::myGlobalMPIComm->dataVectorAllToAll(result, _mpi_storage_offsets, _mpi_storage_sizes);


    this->completeTimeMultTrans_ += this->myTimer_->stop();
}

void DMSystemMatrixVectorizedIdentityMPI::calcDistributionFragment(int totalSize, int procCount, int rank, int *size, int *offset)
{
    int result_size = totalSize / procCount;
    int remainder = totalSize - result_size*procCount;
    int result_offset = 0;
    if(rank < remainder){
        result_size++;
        result_offset = result_size * rank;
    } else {
        result_offset = remainder * (result_size + 1) + (rank - remainder)*result_size;
    }

    *size = result_size;
    *offset = result_offset;
}

void DMSystemMatrixVectorizedIdentityMPI::calcDistribution(int totalSize, int *sizes, int *offsets)
{
    for(int rank = 0; rank < sg::parallel::myGlobalMPIComm->getNumRanks(); ++rank){
        calcDistributionFragment(totalSize, sg::parallel::myGlobalMPIComm->getNumRanks(), rank, &sizes[rank], &offsets[rank]);
    }
}

}
}
