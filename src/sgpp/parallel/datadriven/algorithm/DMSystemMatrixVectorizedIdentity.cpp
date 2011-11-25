/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp"
#include "base/exception/operation_exception.hpp"
#include "parallel/operation/ParallelOpFactory.hpp"
#include "datadriven/operation/DatadrivenOpFactory.hpp"

namespace sg
{
namespace parallel
{

DMSystemMatrixVectorizedIdentity::DMSystemMatrixVectorizedIdentity(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, std::string vecMode)
{
	// handle unsupported vector extensions
	if (vecMode != "X86SIMD" && vecMode != "OCL" && vecMode != "ArBB" && vecMode != "HYBRID_X86SIMD_OCL")
	{
		throw new sg::base::operation_exception("DMSystemMatrixVectorizedIdentity : Only X86SIMD or OCL or ArBB or HYBRID_X86SIMD_OCL are supported vector extensions!");
	}

	resetTimers();

	// create the operations needed in ApplyMatrix
	this->vecMode = vecMode;
	this->lamb = lambda;
	this->data = new sg::base::DataMatrix(trainData);

	if (this->vecMode == "X86SIMD")
	{
		this->vecWidth = 24;
	}
	else if (this->vecMode == "OCL")
	{
		this->vecWidth = 128;
	}
	else if (this->vecMode == "HYBRID_X86SIMD_OCL")
	{
		this->vecWidth = 128;
	}
	else if (this->vecMode == "ArBB")
	{
		this->vecWidth = 16;
	}
	// should not happen because this exception should have been thrown some lines upwards!
	else
	{
		throw new sg::base::operation_exception("DMSystemMatrixVectorizedIdentity : Only X86SIMD or OCL or ArBB or HYBRID_X86SIMD_OCL are supported vector extensions!");
	}

	numTrainingInstances = data->getNrows();

	// Assure that data has a even number of instances -> padding might be needed
	size_t remainder = data->getNrows() % this->vecWidth;
	size_t loopCount = this->vecWidth - remainder;

	if (loopCount != this->vecWidth)
	{
		sg::base::DataVector lastRow(data->getNcols());
		for (size_t i = 0; i < loopCount; i++)
		{
			data->getRow(data->getNrows()-1, lastRow);
			data->resize(data->getNrows()+1);
			data->setRow(data->getNrows()-1, lastRow);
		}
	}

	numPatchedTrainingInstances = data->getNrows();

	if (this->vecMode != "OCL" && this->vecMode != "ArBB"  && this->vecMode != "HYBRID_X86SIMD_OCL")
	{
		data->transpose();
	}

	this->myTimer = new sg::base::SGppStopwatch();

	this->B = sg::op_factory::createOperationMultipleEvalVectorized(SparseGrid, this->vecMode, this->data);
}

DMSystemMatrixVectorizedIdentity::~DMSystemMatrixVectorizedIdentity()
{
	delete this->B;
	delete this->data;
	delete this->myTimer;
}

void DMSystemMatrixVectorizedIdentity::mult(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(numPatchedTrainingInstances);

    // Operation B
	this->myTimer->start();
	this->computeTimeMult += this->B->multVectorized(alpha, temp);
    this->completeTimeMult += this->myTimer->stop();

    // patch result -> set additional entries zero
    if (numTrainingInstances != temp.getSize())
    {
    	for (size_t i = 0; i < (temp.getSize()-numTrainingInstances); i++)
    	{
    		temp.set(temp.getSize()-(i+1), 0.0f);
    	}
    }

    this->myTimer->start();
    this->computeTimeMultTrans += this->B->multTransposeVectorized(temp, result);
    this->completeTimeMultTrans += this->myTimer->stop();

    result.axpy(numTrainingInstances*this->lamb, alpha);
}

void DMSystemMatrixVectorizedIdentity::generateb(sg::base::DataVector& classes, sg::base::DataVector& b)
{
	sg::base::DataVector myClasses(classes);

	// Apply padding
	if (numPatchedTrainingInstances != myClasses.getSize())
	{
		myClasses.resizeZero(numPatchedTrainingInstances);
	}

	this->myTimer->start();
	this->computeTimeMultTrans += this->B->multTransposeVectorized(myClasses, b);
	this->completeTimeMultTrans += this->myTimer->stop();
}

void DMSystemMatrixVectorizedIdentity::rebuildLevelAndIndex()
{
	this->B->rebuildLevelAndIndex();
}

void DMSystemMatrixVectorizedIdentity::resetTimers()
{
	this->completeTimeMult = 0.0;
	this->computeTimeMult = 0.0;
	this->completeTimeMultTrans = 0.0;
	this->computeTimeMultTrans = 0.0;
}

void DMSystemMatrixVectorizedIdentity::getTimers(double& timeMult, double& computeMult, double& timeMultTrans, double& computeMultTrans)
{
	timeMult = this->completeTimeMult;
	computeMult = this->computeTimeMult;
	timeMultTrans = this->completeTimeMultTrans;
	computeMultTrans = this->computeTimeMultTrans;
}

}
}
