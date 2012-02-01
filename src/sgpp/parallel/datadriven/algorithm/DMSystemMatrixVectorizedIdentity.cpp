/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/exception/operation_exception.hpp"

#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp"
#include "parallel/operation/ParallelOpFactory.hpp"

namespace sg
{
namespace parallel
{

DMSystemMatrixVectorizedIdentity::DMSystemMatrixVectorizedIdentity(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode)
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

	// Assure that data has a even number of instances -> padding might be needed
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

	this->B_ = sg::op_factory::createOperationMultipleEvalVectorized(SparseGrid, this->vecMode_, this->dataset_);
}

DMSystemMatrixVectorizedIdentity::~DMSystemMatrixVectorizedIdentity()
{
	delete this->B_;
	delete this->dataset_;
}

void DMSystemMatrixVectorizedIdentity::mult(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(this->numPatchedTrainingInstances_);

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

    result.axpy(static_cast<double>(this->numTrainingInstances_)*this->lambda_, alpha);
}

void DMSystemMatrixVectorizedIdentity::generateb(sg::base::DataVector& classes, sg::base::DataVector& b)
{
	sg::base::DataVector myClasses(classes);

	// Apply padding
	if (this->numPatchedTrainingInstances_ != myClasses.getSize())
	{
		myClasses.resizeZero(this->numPatchedTrainingInstances_);
	}

	this->myTimer_->start();
	this->computeTimeMultTrans_ += this->B_->multTransposeVectorized(myClasses, b);
	this->completeTimeMultTrans_ += this->myTimer_->stop();
}

void DMSystemMatrixVectorizedIdentity::rebuildLevelAndIndex()
{
	this->B_->rebuildLevelAndIndex();
}

}
}
