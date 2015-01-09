/*
 * OperationMultipleEvalTestSuite.hpp
 *
 *  Created on: Dec 19, 2014
 *      Author: pfandedd
 */

#pragma once

#include <string>

#include "datadriven/test/OperationMultipleEvalSimple/OperationMultipleEvalBaseSimpleTestCase.hpp"

namespace sg {
namespace test {

class OperationMultipleEvalBaseTestSuite : public TestSuite {
public:
	std::string getName() override {
		return "OperationMultipleEvalBaseTestSuite";
	}

	void runAllTests() override {
		OperationMultipleEvalBaseSimpleTestCase simple;
		this->run(simple);
	}

};

}
}
