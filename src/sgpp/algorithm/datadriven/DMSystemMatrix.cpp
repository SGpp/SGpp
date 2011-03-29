/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/datadriven/DMSystemMatrix.hpp"
#include "exception/operation_exception.hpp"
using namespace sg::base;


namespace sg
{
namespace datadriven
{

DMSystemMatrix::DMSystemMatrix(Grid& SparseGrid, DataMatrix& trainData, OperationMatrix& C, double lambda)
{
	// create the operations needed in ApplyMatrix
	this->C = &C;
	this->lamb = lambda;
	this->data = &trainData;
	this->B = SparseGrid.createOperationMultipleEval(this->data);
}

DMSystemMatrix::~DMSystemMatrix()
{
	delete this->B;
}

void DMSystemMatrix::mult(DataVector& alpha, DataVector& result)
{
  DataVector temp((*data).getNrows());
  size_t M = (*data).getNrows();

    // Operation B
    this->B->mult(alpha, temp);
    this->B->multTranspose(temp, result);

	DataVector temptwo(alpha.getSize());
	this->C->mult(alpha, temptwo);
	result.axpy(M*this->lamb, temptwo);
}

void DMSystemMatrix::generateb(DataVector& classes, DataVector& b)
{
	this->B->multTranspose(classes, b);
}

}
}
