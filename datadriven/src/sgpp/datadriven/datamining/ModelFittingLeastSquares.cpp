// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "ModelFittingLeastSquares.hpp"

#include <sgpp/datadriven/algorithm/SystemMatrixLeastSquaresIdentity.hpp>
#include <sgpp/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/exception/application_exception.hpp>

//TODO: use the system matrix with flexible regularization

namespace SGPP {
namespace datadriven {

ModelFittingLeastSquares::ModelFittingLeastSquares(
SGPP::datadriven::SampleProvider& sampleProvider,
SGPP::datadriven::DataMiningConfigurationLeastSquares config) :
		datadriven::ModelFittingBase(sampleProvider), configuration(config) {
}

ModelFittingLeastSquares::~ModelFittingLeastSquares() {
}

datadriven::DMSystemMatrixBase*
ModelFittingLeastSquares::createDMSystem(base::DataMatrix& trainDataset,
		float_t lambda) {
	if (this->grid) {
		//TODO: improve
		throw base::application_exception("grid is null");
	}

	datadriven::SystemMatrixLeastSquaresIdentity* systemMatrix =
			new datadriven::SystemMatrixLeastSquaresIdentity(*(this->grid),
					trainDataset, lambda);
	systemMatrix->setImplementation(this->implementationConfiguration);
	return systemMatrix;
}

base::DataVector ModelFittingLeastSquares::predict(
		base::DataMatrix& testDataset) {
	base::DataVector classesComputed(testDataset.getNrows());

	base::OperationMultipleEval* MultEval =
			op_factory::createOperationMultipleEval(*(this->grid), testDataset,
					this->implementationConfiguration);
	MultEval->mult(*alpha, classesComputed);
	delete MultEval;

	return classesComputed;
}

}
}
