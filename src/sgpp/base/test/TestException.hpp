/*
 * TestException.hpp
 *
 *  Created on: Dec 19, 2014
 *      Author: pfandedd
 */

#ifndef TESTEXCEPTION_HPP_
#define TESTEXCEPTION_HPP_

#include <string>
#include <exception>

namespace sg {
namespace test {

class TestException : public std::exception {
private:
	std::string message;
public:
	TestException(std::string message) {
		this->message = message;
	}

	std::string getMessage() {
		return message;
	}
};

}
}

#endif /* TESTEXCEPTION_HPP_ */
