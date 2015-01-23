#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "test/TestCase.hpp"
#include "base/grid/Grid.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "base/operation/OperationMultipleEval.hpp"
#include "datadriven/application/LearnerLeastSquaresIdentity.hpp"
#include "datadriven/operation/DatadrivenOperationCommon.hpp"
#include "datadriven/tools/ARFFTools.hpp"

namespace sg {
namespace test {

class LearnerLeastSquaresIdentityFilesTestCase: public TestCase {
private:
	std::vector<std::string> testFiles = { "tests/data/chess_5d_2000.arff", "tests/data/friedman_4d_2000.arff" };

	std::vector<double> mseEpisilons = { 0.2, 7000.0 };

	size_t dim;
	size_t instances;

	sg::base::DataVector *classesVector = nullptr;
	sg::base::DataMatrix *trainingData = nullptr;

	void readTrainingData(std::string datasetFileName) {
		sg::datadriven::ARFFTools arffTools;
		this->dim = arffTools.getDimension(datasetFileName);
		this->instances = arffTools.getNumberInstances(datasetFileName);

		if (this->classesVector == nullptr) {
			this->classesVector = new sg::base::DataVector(0);
		}
		classesVector->resize(instances);

		arffTools.readClasses(datasetFileName, *classesVector);

		if (this->trainingData != nullptr) {
			delete this->trainingData;
		}
		this->trainingData = new sg::base::DataMatrix(instances, dim);
		arffTools.readTrainingData(datasetFileName, *this->trainingData);
	}
public:

	void run() override {
		sg::datadriven::OperationMultipleEvalConfiguration configuration;
		configuration.type = sg::datadriven::OperationMultipleEvalType::DEFAULT;
		configuration.subType = sg::datadriven::OperationMultipleEvalSubType::DEFAULT;
		testWithFiles(configuration);
		configuration.type = sg::datadriven::OperationMultipleEvalType::STREAMING;
		configuration.subType = sg::datadriven::OperationMultipleEvalSubType::DEFAULT;
		testWithFiles(configuration);
		configuration.type = sg::datadriven::OperationMultipleEvalType::SUBSPACELINEAR;
		configuration.subType = sg::datadriven::OperationMultipleEvalSubType::COMBINED; // ==COMBINED
		testWithFiles(configuration);
		configuration.type = sg::datadriven::OperationMultipleEvalType::SUBSPACELINEAR;
		configuration.subType = sg::datadriven::OperationMultipleEvalSubType::SIMPLE;
		testWithFiles(configuration);
	}

	void testWithFiles(sg::datadriven::OperationMultipleEvalConfiguration &configuration) {

		sg::solver::SLESolverConfiguration SLESolverConfigRefine;
		sg::solver::SLESolverConfiguration SLESolverConfigFinal;
		sg::base::AdpativityConfiguration adaptivityConfiguration;

		// Set solver during refinement
		SLESolverConfigRefine.eps_ = 0;
		SLESolverConfigRefine.maxIterations_ = 20;
		SLESolverConfigRefine.threshold_ = -1.0;
		SLESolverConfigRefine.type_ = sg::solver::CG;

		// Set solver for final step
		SLESolverConfigFinal.eps_ = 0;
		SLESolverConfigFinal.maxIterations_ = 100;
		SLESolverConfigFinal.threshold_ = -1.0;
		SLESolverConfigFinal.type_ = sg::solver::CG;

		// Set Adaptivity
		adaptivityConfiguration.maxLevelType_ = false;
		adaptivityConfiguration.noPoints_ = 80;
		adaptivityConfiguration.numRefinements_ = 7;
		adaptivityConfiguration.percent_ = 100.0;
		adaptivityConfiguration.threshold_ = 0.0;

		for (size_t i = 0; i < testFiles.size(); i++) {
			std::string datasetFileName = this->testFiles[i];
			double mseEpsilon = this->mseEpisilons[i];
			this->readTrainingData(datasetFileName);

			sg::base::RegularGridConfiguration gridConfig;
			gridConfig.dim_ = this->dim;
			gridConfig.type_ = sg::base::Linear;
			gridConfig.level_ = 2;

			sg::datadriven::LearnerLeastSquaresIdentity learner(false, false);
			learner.setImplementation(configuration);

			learner.train(*this->trainingData, *this->classesVector, gridConfig, SLESolverConfigRefine,
					SLESolverConfigFinal, adaptivityConfiguration, false, 0.000001);

			sg::base::DataVector result = learner.predict(*this->trainingData);

			this->assertMeanSquareError(result, *this->classesVector, mseEpsilon);

		}

	}

	std::string getName() override {
		return "LearnerLeastSquaresIdentityFilesTestCase";
	}

};

}
}
