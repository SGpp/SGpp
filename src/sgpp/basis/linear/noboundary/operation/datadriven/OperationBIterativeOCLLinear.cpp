/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationBIterativeOCLLinear.hpp"
#include "exception/operation_exception.hpp"

#ifdef USEOMP
#include "omp.h"
#endif

namespace sg
{

OperationBIterativeOCLLinear::OperationBIterativeOCLLinear(GridStorage* storage) : storage(storage)
{
	Level = new DataMatrix(storage->size(), storage->dim());
	Index = new DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*Level, *Index);

	myTimer = new SGppStopwatch();
	myOCLKernels = new OCLKernels();
}

OperationBIterativeOCLLinear::~OperationBIterativeOCLLinear()
{
	delete Level;
	delete Index;
	delete myTimer;
	delete myOCLKernels;
}

void OperationBIterativeOCLLinear::rebuildLevelAndIndex()
{
	delete Level;
	delete Index;

	Level = new DataMatrix(storage->size(), storage->dim());
	Index = new DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*Level, *Index);
}

double OperationBIterativeOCLLinear::multVectorized(DataVector& alpha, DataMatrix& data, DataVector& result)
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

    if (data.getNcols() % 8 != 0 || sourceSize != data.getNcols())
    {
    	throw operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

    double time = myOCLKernels->multOCL(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, sourceSize, storageSize, dims);

	return time;
}

double OperationBIterativeOCLLinear::multTransposeVectorized(DataVector& alpha, DataMatrix& data, DataVector& result)
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

    if (data.getNcols() % 8 != 0 || result_size != data.getNcols())
    {
    	throw operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

    double time = myOCLKernels->multTransOCL(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);

   	return time;
}

}
