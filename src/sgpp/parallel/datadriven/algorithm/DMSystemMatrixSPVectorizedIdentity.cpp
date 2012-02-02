/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/exception/operation_exception.hpp"

#include "parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentity.hpp"
#include "parallel/operation/ParallelOpFactory.hpp"

namespace sg
{
namespace parallel
{

DMSystemMatrixSPVectorizedIdentity::DMSystemMatrixSPVectorizedIdentity(sg::base::Grid& SparseGrid, sg::base::DataMatrixSP& trainData, float lambda, VectorizationType vecMode)
	:  DMSystemMatrixBaseSP(trainData, lambda), vecMode_(vecMode), vecWidth_(0), numTrainingInstances_(0), numPatchedTrainingInstances_(0)
{
	// handle unsupported vector extensions and set vector width
	// @TODO (heinecke) refactor: better way to set vector width
	if (this->vecMode_ == X86SIMD)
	{
		this->vecWidth_ = 48;
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
		this->vecWidth_ = 192;
	}
	else if (this->vecMode_ == Hybrid_X86SIMD_MIC)
	{
		this->vecWidth_ = 192;
	}
	else
	{
		throw new sg::base::operation_exception("DMSystemMatrixSPVectorizedIdentity : un-supported vector extension!");
	}

	// create the operations needed in ApplyMatrix
	this->dataset_ = new sg::base::DataMatrixSP(trainData);

	this->numTrainingInstances_ = this->dataset_->getNrows();

	// Assure that data has a even number of instances -> padding might be needed
	size_t remainder = this->dataset_->getNrows() % this->vecWidth_;
	size_t loopCount = this->vecWidth_ - remainder;

	if (loopCount != this->vecWidth_)
	{
		sg::base::DataVectorSP lastRow(this->dataset_->getNcols());
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

	this->B_ = sg::op_factory::createOperationMultipleEvalVectorizedSP(SparseGrid, this->vecMode_, this->dataset_);
}

DMSystemMatrixSPVectorizedIdentity::~DMSystemMatrixSPVectorizedIdentity()
{
	delete this->B_;
	delete this->dataset_;
}

void DMSystemMatrixSPVectorizedIdentity::mult(sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& result)
{
	sg::base::DataVectorSP temp(this->numPatchedTrainingInstances_);

    // Operation B
	this->myTimer_->start();
	this->computeTimeMult_ += this->B_->multVectorized(alpha, temp);
    this->completeTimeMult_ += this->myTimer_->stop();

    // patch result -> set additional entries zero
    if (this->numTrainingInstances_ != temp.getSize())
    {
    	for (size_t i = 0; i < (temp.getSize()-this->numTrainingInstances_); i++)
    	{
    		temp.set(temp.getSize()-(i+1), 0.0f);
    	}
    }

    this->myTimer_->start();
    this->computeTimeMultTrans_ += this->B_->multTransposeVectorized(temp, result);
    this->completeTimeMultTrans_ += this->myTimer_->stop();

    result.axpy(static_cast<float>(this->numTrainingInstances_)*this->lambda_, alpha);
}

void DMSystemMatrixSPVectorizedIdentity::generateb(sg::base::DataVectorSP& classes, sg::base::DataVectorSP& b)
{
	sg::base::DataVectorSP myClasses(classes);

	// Apply padding
	if (this->numPatchedTrainingInstances_ != myClasses.getSize())
	{
		myClasses.resizeZero(this->numPatchedTrainingInstances_);
	}

	this->myTimer_->start();
	this->computeTimeMultTrans_ += this->B_->multTransposeVectorized(myClasses, b);
	this->completeTimeMultTrans_ += this->myTimer_->stop();
}

void DMSystemMatrixSPVectorizedIdentity::rebuildLevelAndIndex()
{
	this->B_->rebuildLevelAndIndex();
}

}
}
