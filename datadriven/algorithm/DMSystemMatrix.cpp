/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
#include "datadriven/algorithm/DMSystemMatrix.hpp"
#include "base/exception/operation_exception.hpp"
//#include "datadriven/DatadrivenOpFactory.hpp"
#include "base/operation/BaseOpFactory.hpp"

namespace sg {
namespace datadriven {

//TODO where do I get the kernel from (who constructs the kernel?)

DMSystemMatrix::DMSystemMatrix(sg::base::Grid& grid, sg::base::DataMatrix& trainData, sg::base::OperationMatrix& C,
		double lambda) :
		DMSystemMatrixBase(trainData, lambda), grid(grid) {
	// create the operations needed in ApplyMatrix
	this->C = &C;
	//this->B = sg::op_factory::createOperationMultiEval(grid);
	this->B = sg::op_factory::createOperationMultipleEval(grid, *(this->dataset_));
}

DMSystemMatrix::~DMSystemMatrix() {
	delete this->B;
}

void DMSystemMatrix::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
	//temp variable may be resized during computation
	sg::base::DataVector temp(this->dataset_->getNrows());
	size_t M = this->dataset_->getNrows();

	// Operation B
//      this->B->mult((*this->dataset_), alpha, temp);
//      this->B->multTranspose((*this->dataset_), temp, result);

//this->B->mult(alpha, temp);

	base::OperationMultipleEval *op = sg::op_factory::createOperationMultipleEval(grid, *(this->dataset_));
	sg::base::DataVector temp2(this->dataset_->getNrows());
	//op->mult(*(this->dataset_), alpha, temp2);
	op->mult(alpha, temp);

	/*if (temp.getSize() != temp2.getSize()) {
	 std::cout << "error: sizes don't match" << std::endl;
	 std::cout << "ref: " << temp.getSize() << " mine: " << temp2.getSize() << std::endl;
	 }*/

	/*for (size_t i = 0; i < temp.getSize(); i++) {
	 if (temp[i] != temp2[i]) {
	 if (fabs(temp[i] - temp2[i]) > 0.0000001) {
	 std::cout << "ref: " << temp[i] << " mine: "	<< temp2[i]	<< " i: " << i << std::endl;
	 }
	 }
	 }*/

//      this->B->multTranspose(temp, result);
//      sg::base::DataVector result2(result);
//      op->multTranspose(*(this->dataset_), temp, result2);
//      if (result.getSize() != result2.getSize()) {
//    	  std::cout << "ref: " << result.getSize() << " mine: " << result2.getSize() << std::endl;
//      }
//      for (size_t i = 0; i < result.getSize(); i++) {
//    	  if (result[i] != result2[i]) {
//    		  if (fabs(result[i] - result2[i]) > 0.0000001) {
//    			  std::cout << "ref: " << result[i] << " mine: "	<< result2[i]	<< " i: " << i << std::endl;
//    		  }
//    	  }
//      }
	op->multTranspose(temp, result);

	delete op;

	sg::base::DataVector temptwo(alpha.getSize());
	this->C->mult(alpha, temptwo);
	result.axpy(static_cast<double>(M) * this->lambda_, temptwo);
}

void DMSystemMatrix::generateb(sg::base::DataVector& classes, sg::base::DataVector& b) {
	//this->B->multTranspose((*this->dataset_), classes, b);
	//this->B->multTranspose(classes, b);

	base::OperationMultipleEval *op = sg::op_factory::createOperationMultipleEval(grid, *(this->dataset_));
	op->multTranspose(classes, b);
	delete op;
}

}

}
