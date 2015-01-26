#include <iostream>

#include "base/exception/operation_exception.hpp"

#include "datadriven/DatadrivenOpFactory.hpp"
#include "SystemMatrixLeastSquaresIdentity.hpp"

using namespace std;

namespace sg {
namespace datadriven {

SystemMatrixLeastSquaresIdentity::SystemMatrixLeastSquaresIdentity(sg::base::Grid& grid,
        sg::base::DataMatrix& trainData, double lambda) :
    DMSystemMatrixBase(trainData, lambda), instances(0), paddedInstances(0), grid(grid) {

    this->dataset_ = new sg::base::DataMatrix(trainData);
    this->instances = this->dataset_->getNrows();
    //this->paddedInstances = PaddingAssistant::padDataset(*(this->dataset_));
    //sg::datadriven::OperationMultipleEvalType type = sg::datadriven::OperationMultipleEvalType::SUBSPACELINEAR;
    this->B = sg::op_factory::createOperationMultipleEval(grid, *(this->dataset_), this->implementationConfiguration);
    // padded during Operator construction, fetch new size
    this->paddedInstances = this->dataset_->getNrows();
}

SystemMatrixLeastSquaresIdentity::~SystemMatrixLeastSquaresIdentity() {
    delete this->B;
    delete this->dataset_;
}

void SystemMatrixLeastSquaresIdentity::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
    sg::base::DataVector temp(this->paddedInstances);

    // Operation B
    this->myTimer_->start();
    this->B->mult(alpha, temp);
    this->completeTimeMult_ += this->myTimer_->stop();

    this->myTimer_->start();
    this->B->multTranspose(temp, result);
    this->completeTimeMultTrans_ += this->myTimer_->stop();

    result.axpy(static_cast<double>(this->instances) * this->lambda_, alpha);

}

void SystemMatrixLeastSquaresIdentity::generateb(sg::base::DataVector& classes, sg::base::DataVector& b) {
    sg::base::DataVector myClasses(classes);

    this->myTimer_->start();
    this->B->multTranspose(myClasses, b);
    this->completeTimeMultTrans_ += this->myTimer_->stop();
}

void SystemMatrixLeastSquaresIdentity::rebuildLevelAndIndex() {
    //TODO call to prepare not modeled
    //this->B->prepareExecution();
}

}
}
