// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <iostream>

#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include "SystemMatrixLeastSquaresIdentity.hpp"

using namespace std;

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace datadriven {

SystemMatrixLeastSquaresIdentity::SystemMatrixLeastSquaresIdentity(SGPP::base::Grid& grid,
        SGPP::base::DataMatrix& trainData, double lambda) :
    DMSystemMatrixBase(trainData, lambda), instances(0), paddedInstances(0), grid(grid) {

    this->dataset_ = new SGPP::base::DataMatrix(trainData);
    this->instances = this->dataset_->getNrows();
    //this->paddedInstances = PaddingAssistant::padDataset(*(this->dataset_));
    //SGPP::datadriven::OperationMultipleEvalType type = SGPP::datadriven::OperationMultipleEvalType::SUBSPACELINEAR;
    this->B = SGPP::op_factory::createOperationMultipleEval(grid, *(this->dataset_), this->implementationConfiguration);
    // padded during Operator construction, fetch new size
    this->paddedInstances = this->dataset_->getNrows();
}

SystemMatrixLeastSquaresIdentity::~SystemMatrixLeastSquaresIdentity() {
    delete this->B;
    delete this->dataset_;
}

void SystemMatrixLeastSquaresIdentity::mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
    SGPP::base::DataVector temp(this->paddedInstances);

    // Operation B
    this->myTimer_->start();
    this->B->mult(alpha, temp);
    this->completeTimeMult_ += this->myTimer_->stop();

    this->myTimer_->start();
    this->B->multTranspose(temp, result);
    this->completeTimeMultTrans_ += this->myTimer_->stop();

    result.axpy(static_cast<double>(this->instances) * this->lambda_, alpha);

}

void SystemMatrixLeastSquaresIdentity::generateb(SGPP::base::DataVector& classes, SGPP::base::DataVector& b) {
    SGPP::base::DataVector myClasses(classes);

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
