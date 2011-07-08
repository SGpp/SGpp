/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/

// @author Zhongwen Song (tjsongzw@gmail.com)
// @author Benjamin (pehersto@in.tum.de)
#include "algorithm/datadriven/DMWeightMatrix.hpp"
#include "exception/operation_exception.hpp"
#include "basis/operations_factory.hpp"


namespace sg
{
namespace datadriven
{

DMWeightMatrix::DMWeightMatrix(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, sg::base::OperationMatrix& C, double lambda, sg::base::DataVector& w)
{
	// create the operations needed in ApplyMatrix
	this->C = &C;
	this->lamb = lambda;
	this->data = &trainData;
	//this->B = SparseGrid.createOperationMultipleEval(this->data);
	this->B = sg::GridOperationFactory::createOperationMultipleEval(SparseGrid, this->data);
    this->weight = &w;
}

DMWeightMatrix::~DMWeightMatrix()
{
	delete this->B;
}


void DMWeightMatrix::mult(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp((*data).getNrows());
	sg::base::DataVector tempwithweight((*data).getNrows());
        //size_t M = (*data).getNrows();
        //// Operation B
	this->B->mult(alpha, temp);
    tempwithweight = temp;
    tempwithweight.componentwise_mult(*weight);

    this->B->multTranspose(tempwithweight, result);

	sg::base::DataVector temptwo(alpha.getSize());
	this->C->mult(alpha, temptwo);
	result.axpy(this->lamb, temptwo);
}

void DMWeightMatrix::generateb(sg::base::DataVector& classes, sg::base::DataVector& b)
{
	sg::base::DataVector classeswithweight(classes);
    classeswithweight.componentwise_mult(*weight);
	this->B->multTranspose(classeswithweight, b);
}

}
}
