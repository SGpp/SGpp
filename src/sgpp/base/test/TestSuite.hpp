/*
 * TestSuite.hpp
 *
 *  Created on: Dec 19, 2014
 *      Author: pfandedd
 */

#ifndef TESTSUITE_HPP_
#define TESTSUITE_HPP_

#include <string>
#include <iostream>

#include "base/test/TestCase.hpp"
#include "base/test/TestException.hpp"

namespace sg {
namespace test {

class TestSuite {
private:
	bool success;
public:
	TestSuite() {
		success = true;
	}

	virtual ~TestSuite() {}

	virtual void runAllTests() = 0;

	virtual std::string getName() = 0;

	bool getSuccessState() {
		return success;
	}

	void run(TestCase &test) {
		try {
			test.setup();
		} catch (TestException &e) {
			std::cout << "test " << test.getName() << "failed during setup(): " << e.getMessage() << std::endl;
			success = false;
			return;
		}

		try {
			test.run();
		} catch (TestException &e) {
			std::cout << "test " << test.getName() << "failed during run(): " << e.getMessage() << std::endl;
			success = false;
			return;
		}

		try {
			test.tearDown();
		} catch (TestException &e) {
			std::cout << "test " << test.getName() << "failed during tearDown() : " << e.getMessage() << std::endl;
			success = false;
			return;
		}
		std::cout << "  test " << test.getName() << ": OK" << std::endl;
	}




};

}
}

#endif /* TESTSUITE_HPP_ */
