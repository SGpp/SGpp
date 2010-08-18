/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/datadriven/DMSystemMatrixAVXIdentity.hpp"
#include "exception/operation_exception.hpp"


namespace sg
{

DMSystemMatrixAVXIdentity::DMSystemMatrixAVXIdentity(Grid& SparseGrid, DataMatrix& trainData, double lambda)
{
	// create the operations needed in ApplyMatrix
	this->B = SparseGrid.createOperationBVectorized("AVX");
	this->lamb = lambda;
	this->data = new DataMatrix(trainData);

	numTrainingInstances = data->getNrows();

	// Assure that data has a even number of instances -> padding might be needed
	size_t remainder = data->getNrows() % 4;
	size_t loopCount = 4 - remainder;

	if (loopCount != 4)
	{
		DataVector lastRow(data->getNcols());
		for (size_t i = 0; i < loopCount; i++)
		{
			data->getRow(data->getNrows()-1, lastRow);
			data->resize(data->getNrows()+1);
			data->setRow(data->getNrows()-1, lastRow);
		}
	}
	data->transpose();
}

DMSystemMatrixAVXIdentity::~DMSystemMatrixAVXIdentity()
{
	delete this->B;
	delete this->data;
}

void DMSystemMatrixAVXIdentity::mult(DataVector& alpha, DataVector& result)
{
	DataVector temp((*data).getNcols());

    // Operation B
    this->B->multTransposeVectorized(alpha, (*data), temp);
    // patch result -> set additional entries zero
    if (numTrainingInstances != temp.getSize())
    {
    	for (size_t i = 0; i < (temp.getSize()-numTrainingInstances); i++)
    	{
    		temp.set(temp.getSize()-(i+1), 0.0);
    	}
    }
    this->B->multVectorized(temp, (*data), result);

    result.axpy(numTrainingInstances*this->lamb, alpha);
}

void DMSystemMatrixAVXIdentity::generateb(DataVector& classes, DataVector& b)
{
	DataVector myClasses(classes);

	// Apply padding
	size_t numCols = (*data).getNcols();
	if (numCols != myClasses.getSize())
	{
		myClasses.resizeZero(numCols);
	}
	this->B->multVectorized(myClasses, (*data), b);
}

void DMSystemMatrixAVXIdentity::rebuildLevelAndIndex()
{
	this->B->rebuildLevelAndIndex();
}

}
