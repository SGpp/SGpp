#pragma once

#include <string>

#include "datadriven/test/LearnerLeastSquaresIdentity/LearnerLeastSquaresIdentitySimpleTestCase.hpp"
#include "datadriven/test/LearnerLeastSquaresIdentity/LearnerLeastSquaresIdentityFilesTestCase.hpp"

namespace sg {
namespace test {

class LearnerLeastSquaresIdentityTestSuite : public TestSuite {

public:
	std::string getName() override {
		return "LearnerLeastSquaresIdentityTestSuite";
	}

	void runAllTests() override {
		LearnerLeastSquaresIdentitySimpleTestCase simple;
		LearnerLeastSquaresIdentityFilesTestCase files;
		this->run(simple);
		this->run(files);
	}

};

}
}
