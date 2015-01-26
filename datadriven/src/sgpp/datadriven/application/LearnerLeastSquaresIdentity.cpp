/******************************************************************************
 * Copyright (C) 2012 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
#include <sgpp/datadriven/application/LearnerLeastSquaresIdentity.hpp>

#include "../algorithm/SystemMatrixLeastSquaresIdentity.hpp"
#include <sgpp/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp>
//#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

//TODO how do I construct the operation here (type/subtype)

#include <sgpp/base/exception/factory_exception.hpp>

namespace sg {
namespace datadriven {

LearnerLeastSquaresIdentity::LearnerLeastSquaresIdentity(const bool isRegression, const bool verbose) :
    sg::datadriven::LearnerBase(isRegression, verbose) {
}

LearnerLeastSquaresIdentity::LearnerLeastSquaresIdentity(const std::string tGridFilename,
        const std::string tAlphaFilename, const bool isRegression, const bool verbose) :
    sg::datadriven::LearnerBase(tGridFilename, tAlphaFilename, isRegression, verbose) {
}

LearnerLeastSquaresIdentity::~LearnerLeastSquaresIdentity() {
}

sg::datadriven::DMSystemMatrixBase* LearnerLeastSquaresIdentity::createDMSystem(sg::base::DataMatrix& trainDataset,
        double lambda) {
    if (this->grid_ == NULL)
        return NULL;

    sg::datadriven::SystemMatrixLeastSquaresIdentity *systemMatrix = new sg::datadriven::SystemMatrixLeastSquaresIdentity(*(this->grid_),
            trainDataset, lambda);
    systemMatrix->setImplementation(this->implementationConfiguration);
    return systemMatrix;
}

void LearnerLeastSquaresIdentity::postProcessing(const sg::base::DataMatrix& trainDataset,
        const sg::solver::SLESolverType& solver, const size_t numNeededIterations) {
    LearnerVectorizedPerformance currentPerf = LearnerVectorizedPerformanceCalculator::getGFlopAndGByte(*this->grid_,
            trainDataset.getNrows(), solver, numNeededIterations, sizeof(double));

    this->GFlop_ += currentPerf.GFlop_;
    this->GByte_ += currentPerf.GByte_;

    // Calculate GFLOPS and GBytes/s and write them to console
    if (this->isVerbose_) {
        std::cout << std::endl;
        std::cout << "Current GFlop/s: " << this->GFlop_ / this->execTime_ << std::endl;
        std::cout << "Current GByte/s: " << this->GByte_ / this->execTime_ << std::endl;
        std::cout << std::endl;
    }
}

sg::base::DataVector LearnerLeastSquaresIdentity::predict(sg::base::DataMatrix& testDataset) {
    sg::base::DataVector classesComputed(testDataset.getNrows());

    sg::base::OperationMultipleEval* MultEval = sg::op_factory::createOperationMultipleEval(*(this->grid_), testDataset,
            this->implementationConfiguration);
    MultEval->mult(*alpha_, classesComputed);
    delete MultEval;

    return classesComputed;
}

double LearnerLeastSquaresIdentity::testRegular(const sg::base::RegularGridConfiguration& GridConfig,
        sg::base::DataMatrix& testDataset) {

    InitializeGrid(GridConfig);

    sg::base::OperationMultipleEval* MultEval = sg::op_factory::createOperationMultipleEval(*(this->grid_), testDataset,
            this->implementationConfiguration);

    sg::base::DataVector classesComputed(testDataset.getNrows());

    classesComputed.setAll(0.0);

    execTime_ = 0.0;

    sg::base::SGppStopwatch* myStopwatch = new sg::base::SGppStopwatch();
    myStopwatch->start();

    //TODO could be wrong
    MultEval->mult(*alpha_, classesComputed);
    double stopTime = myStopwatch->stop();
    this->execTime_ += stopTime;
    std::cout << "execution duration: " << this->execTime_ << std::endl;
    delete MultEval;

    return stopTime;
}

std::vector<std::pair<size_t, double> > LearnerLeastSquaresIdentity::getRefinementExecTimes() {
    return this->ExecTimeOnStep;
}

}
}
