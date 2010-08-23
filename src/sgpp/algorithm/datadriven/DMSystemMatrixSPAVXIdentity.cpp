/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/datadriven/DMSystemMatrixSPAVXIdentity.hpp"
#include "exception/operation_exception.hpp"


namespace sg
{

DMSystemMatrixSPAVXIdentity::DMSystemMatrixSPAVXIdentity(Grid& SparseGrid, DataMatrixSP& trainData, float lambda)
{
	// create the operations needed in ApplyMatrix
	this->B = SparseGrid.createOperationBVectorizedSP("AVX");
	this->lamb = lambda;
	this->data = new DataMatrixSP(trainData);

	numTrainingInstances = data->getNrows();

	// Assure that data has a even number of instances -> padding might be needed
	size_t remainder = data->getNrows() % 8;
	size_t loopCount = 8 - remainder;

	if (loopCount != 8)
	{
		DataVectorSP lastRow(data->getNcols());
		for (size_t i = 0; i < loopCount; i++)
		{
			data->getRow(data->getNrows()-1, lastRow);
			data->resize(data->getNrows()+1);
			data->setRow(data->getNrows()-1, lastRow);
		}
	}
	data->transpose();
}

DMSystemMatrixSPAVXIdentity::~DMSystemMatrixSPAVXIdentity()
{
	delete this->B;
	delete this->data;
}

void DMSystemMatrixSPAVXIdentity::mult(DataVectorSP& alpha, DataVectorSP& result)
{
	DataVectorSP temp((*data).getNcols());

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

void DMSystemMatrixSPAVXIdentity::generateb(DataVectorSP& classes, DataVectorSP& b)
{
	DataVectorSP myClasses(classes);

	// Apply padding
	size_t numCols = (*data).getNcols();
	if (numCols != myClasses.getSize())
	{
		myClasses.resizeZero(numCols);
	}
	this->B->multVectorized(myClasses, (*data), b);
}

void DMSystemMatrixSPAVXIdentity::rebuildLevelAndIndex()
{
	this->B->rebuildLevelAndIndex();
}

}
