/*
 * OperationMultipleEvalTestSuite.hpp
 *
 *  Created on: Dec 19, 2014
 *      Author: pfandedd
 */

#pragma once

#include <string>

#include "datadriven/test/OperationMultipleEval/OperationMultipleEvalSubspaceSimpleTestCase.hpp"

namespace sg {
namespace test {

class OperationMultipleEvalSubspaceTestSuite : public TestSuite {
public:
	std::string getName() override {
		return "OperationMultipleEvalSubspaceTestSuite";
	}

	void runAllTests() override {
		OperationMultipleEvalSubspaceSimpleTestCase simple;
		this->run(simple);
	}

};

}
}
