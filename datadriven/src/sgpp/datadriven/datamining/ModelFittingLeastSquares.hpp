// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include "ModelFittingBase.hpp"

#include "DataMiningConfigurationLeastSquares.hpp"

#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>

namespace SGPP {
namespace datadriven {

/**
 * This class implements standard sparse grid regression
 * with an Identity matrix as regularization operator.
 *
 * Furthermore this Learner provides support for several
 * vectorization approaches covering GPUs, CPUs and coprocessors.
 */
class ModelFittingLeastSquares: public ModelFittingBase {
private:
	std::vector<std::pair<size_t, float_t> > ExecTimeOnStep;

	base::OperationMultipleEval* kernel = nullptr;

	datadriven::OperationMultipleEvalConfiguration implementationConfiguration;

	DataMiningConfigurationLeastSquares configuration;
protected:

	virtual datadriven::DMSystemMatrixBase* createDMSystem(
			base::DataMatrix& trainDataset, float_t lambda);


public:
	/**
	 * Constructor
	 *
	 * @param isRegression set to true if a regression task should be executed
	 * @param isVerbose set to true in order to allow console output
	 */
	ModelFittingLeastSquares(SGPP::datadriven::SampleProvider& sampleProvider,
	SGPP::datadriven::DataMiningConfigurationLeastSquares config);

	/**
	 * Destructor
	 */
	virtual ~ModelFittingLeastSquares();

    void fit() override;

	void setImplementation(
			datadriven::OperationMultipleEvalConfiguration operationConfiguration) {
		this->implementationConfiguration = operationConfiguration;
	}
};

}

}

