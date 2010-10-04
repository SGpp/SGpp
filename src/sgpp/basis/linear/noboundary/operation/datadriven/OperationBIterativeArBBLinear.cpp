/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationBIterativeArBBLinear.hpp"
#include "exception/operation_exception.hpp"

namespace sg
{

OperationBIterativeArBBLinear::OperationBIterativeArBBLinear(GridStorage* storage) : storage(storage)
{
	Level = new DataMatrix(storage->size(), storage->dim());
	Index = new DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*Level, *Index);

	myTimer = new SGppStopwatch();
	myArBBKernels = new ArBBKernels();
}

OperationBIterativeArBBLinear::~OperationBIterativeArBBLinear()
{
	delete Level;
	delete Index;
	delete myTimer;
	delete myArBBKernels;
}

void OperationBIterativeArBBLinear::rebuildLevelAndIndex()
{
	delete Level;
	delete Index;

	Level = new DataMatrix(storage->size(), storage->dim());
	Index = new DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*Level, *Index);

	myArBBKernels->resetKernels();
}

double OperationBIterativeArBBLinear::multVectorized(DataVector& alpha, DataMatrix& data, DataVector& result)
{
	size_t sourceSize = alpha.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    result.setAll(0.0);

    double* ptrSource = alpha.getPointer();
    double* ptrData = data.getPointer();
    double* ptrLevel = this->Level->getPointer();
    double* ptrIndex = this->Index->getPointer();
    double* ptrGlobalResult = result.getPointer();

    if (data.getNrows() % 16 != 0 || sourceSize != data.getNrows())
    {
    	throw operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

    double time = myArBBKernels->multArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, sourceSize, storageSize, dims);

	return time;
}

double OperationBIterativeArBBLinear::multTransposeVectorized(DataVector& alpha, DataMatrix& data, DataVector& result)
{
	size_t result_size = result.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    result.setAll(0.0);

    double* ptrAlpha = alpha.getPointer();
    double* ptrData = data.getPointer();
    double* ptrResult = result.getPointer();
    double* ptrLevel = this->Level->getPointer();
    double* ptrIndex = this->Index->getPointer();

    if (data.getNrows() % 16 != 0 || result_size != data.getNrows())
    {
    	throw operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

    double time = myArBBKernels->multTransArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);

   	return time;
}

}
