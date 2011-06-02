/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/datadriven/DMSystemMatrix.hpp"
#include "exception/operation_exception.hpp"
#include "basis/operations_factory.hpp"


namespace sg
{
namespace datadriven
{

DMSystemMatrix::DMSystemMatrix(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, sg::base::OperationMatrix& C, double lambda)
{
	// create the operations needed in ApplyMatrix
	this->C = &C;
	this->lamb = lambda;
	this->data = &trainData;
	this->B = sg::GridOperationFactory::createOperationMultipleEval(SparseGrid, this->data);
}

DMSystemMatrix::~DMSystemMatrix()
{
	delete this->B;
}

void DMSystemMatrix::mult(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
  sg::base::DataVector temp((*data).getNrows());
  size_t M = (*data).getNrows();

    // Operation B
    this->B->mult(alpha, temp);
    this->B->multTranspose(temp, result);

	sg::base::DataVector temptwo(alpha.getSize());
	this->C->mult(alpha, temptwo);
	result.axpy(M*this->lamb, temptwo);
}

void DMSystemMatrix::generateb(sg::base::DataVector& classes, sg::base::DataVector& b)
{
	this->B->multTranspose(classes, b);
}

}
}
