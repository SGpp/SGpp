/*
 * OperationMultipleEvalTestSuite.hpp
 *
 *  Created on: Dec 19, 2014
 *      Author: pfandedd
 */

#pragma once

#include <string>

#include "datadriven/test/OperationMultipleEvalStreaming/OperationMultipleEvalStreamingSimpleTestCase.hpp"

namespace sg {
namespace test {

class OperationMultipleEvalStreamingTestSuite: public TestSuite {
public:
  std::string getName() override {
    return "OperationMultipleEvalStreamingTestSuite";
  }

  void runAllTests() override {
    OperationMultipleEvalStreamingSimpleTestCase simple;
    this->run(simple);
  }

};

}
}
