/*
 * TestCase.hpp
 *
 *  Created on: Dec 19, 2014
 *      Author: pfandedd
 */

#ifndef SRC_SGPP_BASE_APPLICATION_TESTCASE_HPP_
#define SRC_SGPP_BASE_APPLICATION_TESTCASE_HPP_

#include <string>
#include <cmath>

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
		if (fabs(value - reference) > epsilon) {
			success = false;
		}
	}
};

}
}

#endif /* SRC_SGPP_BASE_APPLICATION_TESTCASE_HPP_ */
