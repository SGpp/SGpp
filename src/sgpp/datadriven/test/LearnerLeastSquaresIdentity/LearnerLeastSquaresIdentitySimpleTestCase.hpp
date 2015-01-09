#pragma once

#include <iostream>
#include <string>

#include "base/test/TestCase.hpp"
#include "base/grid/Grid.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "base/operation/OperationMultipleEval.hpp"
#include "datadriven/application/LearnerLeastSquaresIdentity.hpp"
#include "datadriven/operation/DatadrivenOperationCommon.hpp"

namespace sg {
namespace test {

class LearnerLeastSquaresIdentitySimpleTestCase: public TestCase {
public:

	void run() override {
		std::cout << "testSubspaceLinearSimple();" << std::endl;
		testSubspaceLinearSimple();
		std::cout << "testSubspaceLinearCombined();" << std::endl;
		testSubspaceLinearCombined();
		std::cout << "testStreaming()" << std::endl;
		testStreaming();
	}

	void testSubspaceLinearCombined() {

		sg::base::RegularGridConfiguration gridConfig;
		gridConfig.dim_ = 2;
		gridConfig.type_ = sg::base::Linear;
		gridConfig.level_ = 4;

		//sg::solver::SLESolverConfiguration SLESolverConfigRefine;
		sg::solver::SLESolverConfiguration SLESolverConfigFinal;

//		// Set solver during refinement
//		SLESolverConfigRefine.eps_ = 0;
//		SLESolverConfigRefine.maxIterations_ = 10;
//		SLESolverConfigRefine.threshold_ = -1.0;
//		SLESolverConfigRefine.type_ = sg::solver::CG;

// Set solver for final step
		SLESolverConfigFinal.eps_ = 0;
		SLESolverConfigFinal.maxIterations_ = 5;
		SLESolverConfigFinal.threshold_ = -1.0;
		SLESolverConfigFinal.type_ = sg::solver::CG;

		sg::datadriven::OperationMultipleEvalConfiguration configuration;
		configuration.type = sg::datadriven::OperationMultipleEvalType::SUBSPACELINEAR;
		configuration.subType = sg::datadriven::OperationMultipleEvalSubType::COMBINED;

		size_t numberDataPoints = 3;

		sg::base::DataVector alpha(numberDataPoints);
		for (size_t i = 0; i < alpha.getSize(); i++) {
			alpha[i] = static_cast<double>(i + 1);
		}

		//points = [[0.5, 0.1], [0.3, 0.4], [0.9, 0.7]]
		double pointsArray[] = { 0.5, 0.1, 0.3, 0.4, 0.9, 0.7 };
		sg::base::DataMatrix evalPoints(pointsArray, numberDataPoints, gridConfig.dim_);
		double classesArray[] = { 1.0, -1.0, 1.0 };
		sg::base::DataVector classes(classesArray, numberDataPoints);

		sg::datadriven::LearnerLeastSquaresIdentity learner(false, false);
		learner.setImplementation(configuration);

		learner.train(evalPoints, classes, gridConfig, SLESolverConfigFinal, 0.000000001);

		sg::base::DataVector result = learner.predict(evalPoints);

		this->assertEqual(result[0], 1.0, 1E-7);
		this->assertEqual(result[1], -1.0, 1E-7);
		this->assertEqual(result[2], 1.0, 1E-7);

	}

	void testSubspaceLinearSimple() {

		sg::base::RegularGridConfiguration gridConfig;
		gridConfig.dim_ = 2;
		gridConfig.type_ = sg::base::Linear;
		gridConfig.level_ = 4;

		//sg::solver::SLESolverConfiguration SLESolverConfigRefine;
		sg::solver::SLESolverConfiguration SLESolverConfigFinal;

//		// Set solver during refinement
//		SLESolverConfigRefine.eps_ = 0;
//		SLESolverConfigRefine.maxIterations_ = 10;
//		SLESolverConfigRefine.threshold_ = -1.0;
//		SLESolverConfigRefine.type_ = sg::solver::CG;

// Set solver for final step
		SLESolverConfigFinal.eps_ = 0;
		SLESolverConfigFinal.maxIterations_ = 5;
		SLESolverConfigFinal.threshold_ = -1.0;
		SLESolverConfigFinal.type_ = sg::solver::CG;

		sg::datadriven::OperationMultipleEvalConfiguration configuration;
		configuration.type = sg::datadriven::OperationMultipleEvalType::SUBSPACELINEAR;
		configuration.subType = sg::datadriven::OperationMultipleEvalSubType::SIMPLE;

		size_t numberDataPoints = 3;

		sg::base::DataVector alpha(numberDataPoints);
		for (size_t i = 0; i < alpha.getSize(); i++) {
			alpha[i] = static_cast<double>(i + 1);
		}

		//points = [[0.5, 0.1], [0.3, 0.4], [0.9, 0.7]]
		double pointsArray[] = { 0.5, 0.1, 0.3, 0.4, 0.9, 0.7 };
		sg::base::DataMatrix evalPoints(pointsArray, numberDataPoints, gridConfig.dim_);
		double classesArray[] = { 1.0, -1.0, 1.0 };
		sg::base::DataVector classes(classesArray, numberDataPoints);

		sg::datadriven::LearnerLeastSquaresIdentity learner(false, false);
		learner.setImplementation(configuration);

		learner.train(evalPoints, classes, gridConfig, SLESolverConfigFinal, 0.000000001);

		sg::base::DataVector result = learner.predict(evalPoints);

		this->assertEqual(result[0], 1.0, 1E-7);
		this->assertEqual(result[1], -1.0, 1E-7);
		this->assertEqual(result[2], 1.0, 1E-7);

	}

	void testStreaming() {

		sg::base::RegularGridConfiguration gridConfig;
		gridConfig.dim_ = 2;
		gridConfig.type_ = sg::base::Linear;
		gridConfig.level_ = 4;

		//sg::solver::SLESolverConfiguration SLESolverConfigRefine;
		sg::solver::SLESolverConfiguration SLESolverConfigFinal;

//		// Set solver during refinement
//		SLESolverConfigRefine.eps_ = 0;
//		SLESolverConfigRefine.maxIterations_ = 10;
//		SLESolverConfigRefine.threshold_ = -1.0;
//		SLESolverConfigRefine.type_ = sg::solver::CG;

// Set solver for final step
		SLESolverConfigFinal.eps_ = 0;
		SLESolverConfigFinal.maxIterations_ = 5;
		SLESolverConfigFinal.threshold_ = -1.0;
		SLESolverConfigFinal.type_ = sg::solver::CG;

		sg::datadriven::OperationMultipleEvalConfiguration configuration;
		configuration.type = sg::datadriven::OperationMultipleEvalType::STREAMING;
		configuration.subType = sg::datadriven::OperationMultipleEvalSubType::DEFAULT;

		size_t numberDataPoints = 3;

		sg::base::DataVector alpha(numberDataPoints);
		for (size_t i = 0; i < alpha.getSize(); i++) {
			alpha[i] = static_cast<double>(i + 1);
		}

		//points = [[0.5, 0.1], [0.3, 0.4], [0.9, 0.7]]
		double pointsArray[] = { 0.5, 0.1, 0.3, 0.4, 0.9, 0.7 };
		sg::base::DataMatrix evalPoints(pointsArray, numberDataPoints, gridConfig.dim_);
		double classesArray[] = { 1.0, -1.0, 1.0 };
		sg::base::DataVector classes(classesArray, numberDataPoints);

		sg::datadriven::LearnerLeastSquaresIdentity learner(false, false);
		learner.setImplementation(configuration);

		learner.train(evalPoints, classes, gridConfig, SLESolverConfigFinal, 0.000000001);

		sg::base::DataVector result = learner.predict(evalPoints);

		this->assertEqual(result[0], 1.0, 1E-7);
		this->assertEqual(result[1], -1.0, 1E-7);
		this->assertEqual(result[2], 1.0, 1E-7);

	}

	std::string getName() override {
		return "LearnerLeastSquaresIdentitySimpleTestCase";
	}

};

}
}
