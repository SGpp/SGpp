/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/datadriven/DMSystemMatrixSPSSEIdentity.hpp"
#include "exception/operation_exception.hpp"


namespace sg
{

DMSystemMatrixSPSSEIdentity::DMSystemMatrixSPSSEIdentity(Grid& SparseGrid, DataMatrixSP& trainData, float lambda)
{
	// create the operations needed in ApplyMatrix
	this->B = SparseGrid.createOperationBVectorizedSP("SSE");
	this->lamb = lambda;
	this->data = new DataMatrixSP(trainData);

	numTrainingInstances = data->getNrows();

	// Assure that data has a even number of instances -> padding might be needed
	if (data->getNrows() % 2 != 0)
	{
		DataVectorSP lastRow(data->getNcols());
		data->getRow(data->getNrows()-1, lastRow);
		data->resize(data->getNrows()+1);
		data->setRow(data->getNrows()-1, lastRow);
	}
	data->transpose();
}

DMSystemMatrixSPSSEIdentity::~DMSystemMatrixSPSSEIdentity()
{
	delete this->B;
	delete this->data;
}

void DMSystemMatrixSPSSEIdentity::mult(DataVectorSP& alpha, DataVectorSP& result)
{
	DataVectorSP temp((*data).getNcols());

    // Operation B
    this->B->multTransposeVectorized(alpha, (*data), temp);
    // patch result -> set additional entries zero
    if (numTrainingInstances != temp.getSize())
    {
    	temp.set(temp.getSize()-1, 0.0);
    }
    this->B->multVectorized(temp, (*data), result);

    result.axpy(numTrainingInstances*this->lamb, alpha);
}

void DMSystemMatrixSPSSEIdentity::generateb(DataVectorSP& classes, DataVectorSP& b)
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

void DMSystemMatrixSPSSEIdentity::rebuildLevelAndIndex()
{
	this->B->rebuildLevelAndIndex();
}

}
