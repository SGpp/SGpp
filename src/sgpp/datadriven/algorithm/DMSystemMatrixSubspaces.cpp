#include <iostream>

#include "base/exception/operation_exception.hpp"

#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
#include "DMSystemMatrixSubspaces.hpp"
#include "datadriven/DatadrivenOpFactory.hpp"

using namespace std;

namespace sg {
namespace datadriven {

DMSystemMatrixSubspaces::DMSystemMatrixSubspaces(sg::base::Grid& grid, sg::base::DataMatrix& trainData, double lambda) :
		DMSystemMatrixBase(trainData, lambda), instances(0), paddedInstances(0), grid(grid)  {

	this->dataset_ = new sg::base::DataMatrix(trainData);
	this->instances = this->dataset_->getNrows();
	//this->paddedInstances = PaddingAssistant::padDataset(*(this->dataset_));
	//sg::datadriven::OperationMultipleEvalType type = sg::datadriven::OperationMultipleEvalType::SUBSPACELINEAR;
	this->B = sg::op_factory::createOperationMultipleEval(grid, *(this->dataset_), this->implementationConfiguration);
	// padded during Operator construction, fetch new size
	this->paddedInstances = this->dataset_->getNrows();
}

DMSystemMatrixSubspaces::~DMSystemMatrixSubspaces() {
	delete this->B;
	delete this->dataset_;
}

void DMSystemMatrixSubspaces::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
	sg::base::DataVector temp(this->paddedInstances);

	// Operation B
	this->myTimer_->start();
	this->B->mult(alpha, temp);
	this->completeTimeMult_ += this->myTimer_->stop();

	//TODO don't pad here, because padding should be done within the operations to ensure consistency
	// patch result -> set additional entries zero
//	if (this->instances != temp.getSize()) {
//		for (size_t i = 0; i < (temp.getSize() - this->instances); i++) {
//			temp.set(temp.getSize() - (i + 1), 0.0f);
//		}
//	}

	this->myTimer_->start();
	this->B->multTranspose(temp, result);
	this->completeTimeMultTrans_ += this->myTimer_->stop();

	result.axpy(static_cast<double>(this->instances) * this->lambda_, alpha);

}

void DMSystemMatrixSubspaces::generateb(sg::base::DataVector& classes, sg::base::DataVector& b) {
	sg::base::DataVector myClasses(classes);

	//TODO don't pad here, because padding should be done within the operations to ensure consistency
	// Apply padding
//	if (this->paddedInstances != myClasses.getSize()) {
//		myClasses.resizeZero(this->paddedInstances);
//	}

	this->myTimer_->start();
	this->B->multTranspose(myClasses, b);
	this->completeTimeMultTrans_ += this->myTimer_->stop();
}

void DMSystemMatrixSubspaces::rebuildLevelAndIndex() {
	//TODO call to prepare not modeled
	//this->B->prepareExecution();
}

}
}
