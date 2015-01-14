#pragma once

#include <string>

#include "datadriven/test/MetaLearner/MetaLearnerLearnReferenceTestCase.hpp"

namespace sg {
namespace test {

class MetaLearnerTestSuite : public TestSuite {
public:
    std::string getName() override {
        return "MetaLearnerTestSuite";
    }

    void runAllTests() override {
        MetaLearnerLearnReferenceTestCase learnReference;
        this->run(learnReference);
    }

};

}
}
