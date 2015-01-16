/*
 * TestCase.hpp
 *
 *  Created on: Dec 19, 2014
 *      Author: pfandedd
 */

#ifndef SRC_SGPP_BASE_APPLICATION_TESTCASE_HPP_
#define SRC_SGPP_BASE_APPLICATION_TESTCASE_HPP_

#include <string>
#include <iomanip>
#include <cmath>
#include "TestException.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
namespace test {

class TestCase {
private:
	bool success = 0;
public:
	TestCase() {
		success = true;
	}

	virtual ~TestCase() {

	}

	virtual void setup() {

	}

	//don't call the run method directly, use TestSuite.runTest() instead
	virtual void run() = 0;

	virtual std::string getName() = 0;

	virtual void tearDown() {

	}

	bool getSuccessStatus() {
		return success;
	}

	void assertTrue(bool expression) {
		if (!expression) {
			success = false;
		}
	}

	void assertEqual(double value, double reference, double epsilon = 0.0) {
		double temp = fabs(value - reference);
		if (temp > epsilon || (temp != 0.0 && !std::isnormal(temp))) {
			success = false;
			std::cout << std::setprecision(10);
			std::cout << std::fixed;
			std::cout << "assertEqual: value(" << value << " - " << reference << ") = " << temp << " (epsilon: " << epsilon << ")" << std::endl;
		}
	}

	void assertMeanSquareError(sg::base::DataVector &left, sg::base::DataVector &right, double epsilon = 0.0) {
		if (left.getSize() != right.getSize()) {
			success = false;
			std::cout << "assertMeanSquareError: input sizes don't match" << std::endl;
		}

		double mse = 0.0;
		for (size_t i = 0; i < left.getSize(); i++) {
			double temp = left[i] - right[i];
			temp *= temp;
			mse += temp;
		}
		mse = mse / static_cast<double>(left.getSize());

		if (mse > epsilon || (mse != 0.0 && !std::isnormal(mse))) {
			success = false;
			std::cout << std::setprecision(10);
			std::cout << std::fixed;
			std::cout << "assertMeanSquareError: mse = " << mse << " (epsilon: " << epsilon << ")" << std::endl;
		}
	}
};

}
}

#endif /* SRC_SGPP_BASE_APPLICATION_TESTCASE_HPP_ */
