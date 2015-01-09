/*
 * OperationMultipleEvalTestSuite.hpp
 *
 *  Created on: Dec 19, 2014
 *      Author: pfandedd
 */

#pragma once

#include <string>

#include "datadriven/test/OperationMultipleEvalSubspace/OperationMultipleEvalSubspaceSimpleTestCase.hpp"

namespace sg {
namespace test {

class OperationMultipleEvalSubspaceTestSuite : public TestSuite {
public:
	std::string getName() override {
		return "LearnerLeastSquaresIdentityTestSuite";
	}

	void runAllTests() override {
		OperationMultipleEvalSubspaceSimpleTestCase simple;
		this->run(simple);
	}

};

}
}
