#include <iostream>

#include "base/exception/operation_exception.hpp"

#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
#include "DMSystemMatrixSubspaces.hpp"
#include "datadriven/DatadrivenOpFactory.hpp"

using namespace std;

namespace sg {
namespace datadriven {

DMSystemMatrixSubspaces::DMSystemMatrixSubspaces(sg::base::Grid& grid, sg::base::DataMatrix& trainData, double lambda) :
		DMSystemMatrixBase(trainData, lambda), instances(0), paddedInstances(0) {

	this->dataset_ = new sg::base::DataMatrix(trainData);
	this->instances = this->dataset_->getNrows();
	//this->paddedInstances = PaddingAssistant::padDataset(*(this->dataset_));
	sg::datadriven::OperationMultipleEvalType type = sg::datadriven::OperationMultipleEvalType::SUBSPACELINEAR;
	this->B = sg::op_factory::createOperationMultipleEval(grid, *(this->dataset_), type);
	// padded during Operator construction, fetch new size
	this->paddedInstances = this->dataset_->getNrows();
	//cout << "paddedInstances: " << this->paddedInstances << endl;
}

DMSystemMatrixSubspaces::~DMSystemMatrixSubspaces() {
	delete this->B;
	delete this->dataset_;
}

void DMSystemMatrixSubspaces::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
	sg::base::DataVector temp(this->paddedInstances);

	// Operation B
	this->myTimer_->start();
	this->B->multTranspose(alpha, temp);
	this->computeTimeMultTrans_ += this->B->getLastOperationDuration();
	this->completeTimeMult_ += this->myTimer_->stop();

	/*for (size_t i = 0; i < temp.getSize(); i++) {
	 cout << "endi: " << i << endl;
	 cout << " => " << temp.get(i) << endl;
	 }*/

	// patch result -> set additional entries zero
	if (this->instances != temp.getSize()) {
		for (size_t i = 0; i < (temp.getSize() - this->instances); i++) {
			temp.set(temp.getSize() - (i + 1), 0.0f);
		}
	}

	this->myTimer_->start();
	this->B->mult(temp, result);
	this->computeTimeMult_ += this->B->getLastOperationDuration();

	this->completeTimeMultTrans_ += this->myTimer_->stop();

	/*cout << "result: " << result.toString() << endl;
	 cout << "instances: " << this->instances << endl;
	 cout << "lambda: " << this->lambda_ << endl;*/

	result.axpy(static_cast<double>(this->instances) * this->lambda_, alpha);
	//cout << "result2: " << result.toString() << endl;
}

void DMSystemMatrixSubspaces::generateb(sg::base::DataVector& classes, sg::base::DataVector& b) {
	sg::base::DataVector myClasses(classes);

	// Apply padding
	if (this->paddedInstances != myClasses.getSize()) {
		myClasses.resizeZero(this->paddedInstances);
	}

	this->myTimer_->start();
	this->B->mult(myClasses, b);
	this->computeTimeMultTrans_ += this->B->getLastOperationDuration();
	this->completeTimeMultTrans_ += this->myTimer_->stop();
}

void DMSystemMatrixSubspaces::rebuildLevelAndIndex() {
	//TODO call to prepare not modeled
	//this->B->prepareExecution();
}

}
}
