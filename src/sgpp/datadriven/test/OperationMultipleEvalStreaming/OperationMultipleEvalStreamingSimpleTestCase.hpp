#pragma once

#include <iostream>
#include <string>

#include "base/test/TestCase.hpp"
#include "base/grid/Grid.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "datadriven/DatadrivenOpFactory.hpp"
#include "datadriven/operation/DatadrivenOperationCommon.hpp"
#include "base/operation/OperationMultipleEval.hpp"

namespace sg {
namespace test {

class OperationMultipleEvalStreamingSimpleTestCase: public TestCase {
public:

	void run() override {

		sg::base::RegularGridConfiguration gridConfig;
		gridConfig.dim_ = 2;
		gridConfig.type_ = sg::base::Linear;
		gridConfig.level_ = 2;

		sg::base::Grid *grid = sg::base::Grid::createLinearGrid(gridConfig.dim_);
		sg::base::GridGenerator* myGenerator = grid->createGridGenerator();
		myGenerator->regular(gridConfig.level_);
		delete myGenerator;

		sg::base::GridStorage *gridStorage = grid->getStorage();
		size_t gridSize = gridStorage->size();

		sg::base::DataVector alpha(gridSize);
		for (size_t i = 0; i < alpha.getSize(); i++) {
			alpha[i] = static_cast<double>(i + 1);
		}

		size_t numberDataPoints = 3;

		//points = [[0.5, 0.1], [0.3, 0.4], [0.9, 0.7]]
		double pointsArray[] = { 0.5, 0.1, 0.3, 0.4, 0.9, 0.7 };
		sg::base::DataMatrix evalPoints(pointsArray, numberDataPoints, gridConfig.dim_);

		sg::base::DataVector result(numberDataPoints);

		sg::datadriven::OperationMultipleEvalConfiguration configuration;
		configuration.type = sg::datadriven::OperationMultipleEvalType::STREAMING;
		sg::base::OperationMultipleEval *multiEvalOp = sg::op_factory::createOperationMultipleEval(*grid, evalPoints,
				configuration);

		multiEvalOp->mult(alpha, result);
		delete multiEvalOp;

		this->assertEqual(result[0], 1.8, 1E-10);
		this->assertEqual(result[1], 2.72, 1E-10);
		this->assertEqual(result[2], 1.64, 1E-10);

		delete grid;
	}

	std::string getName() override {
		return "OperationMultipleEvalStreamingSimpleTestCase";
	}

};

}
}
