#if USE_OCL == 1

#include <sgpp/globaldef.hpp>
#include <sgpp/datadriven/opencl/OCLConfigurationParameters.hpp>
#include <sgpp/datadriven/application/StaticParameterTuner.hpp>
#include "sgpp/datadriven/application/MetaLearner.hpp"
#include <limits>

namespace SGPP {
namespace datadriven {

StaticParameterTuner::StaticParameterTuner() :
        writeStatistics(false), outFileName("") {

}

//write a file with detailed stats for the optimization
StaticParameterTuner::StaticParameterTuner(std::string outFileName) :
        writeStatistics(true), outFileName(outFileName) {

}

void StaticParameterTuner::addFixedParameter(std::string name, std::string value) {
    this->fixedParameters.set(name, value);
}

void StaticParameterTuner::addParameter(std::string name, std::vector<std::string> valueRange) {
    this->tunableParameters.push_back(TunableParameter(name, valueRange));
}

SGPP::base::OCLConfigurationParameters StaticParameterTuner::tuneParameters(
SGPP::datadriven::LearnerScenario scenario) {

    SGPP::base::OCLConfigurationParameters currentParameters = fixedParameters;

    //create initial parameter combination
    std::vector<size_t> valueIndices(tunableParameters.size());
    for (size_t i = 0; i < valueIndices.size(); i++) {
        valueIndices[i] = 0;
        TunableParameter &parameter = tunableParameters[i];
        currentParameters.set(parameter.getName(), parameter.getValues()[0]);
    }

    // evaluate initial parameter combination
    double shortestDuration = evaluateSetup(scenario, currentParameters);
    SGPP::base::OCLConfigurationParameters bestParameters = currentParameters;

    size_t parameterIndex = 0;
    while (parameterIndex < tunableParameters.size()) {

        // enumerate next parameter combination
        // can the current index be increased?
        TunableParameter &parameter = tunableParameters[parameterIndex];
        if (valueIndices[parameterIndex] + 1 < parameter.getValues().size()) {
            valueIndices[parameterIndex] += 1;
            currentParameters.set(parameter.getName(), parameter.getValues()[valueIndices[parameterIndex]]);
            //reset lower indices
            for (size_t i = 0; i < parameterIndex; i++) {
                valueIndices[i] = 0;
                TunableParameter &parameterForReset = tunableParameters[i];
                currentParameters.set(parameterForReset.getName(), parameterForReset.getValues()[0]);
            }
            parameterIndex = 0;

            // evaluate current parameter combination
            double duration = evaluateSetup(scenario, currentParameters);
            if (duration < shortestDuration) {
                shortestDuration = duration;
                bestParameters = currentParameters;
            }
        } else {
            parameterIndex += 1;
        }
    }

    return bestParameters;
}

double StaticParameterTuner::evaluateSetup(SGPP::datadriven::LearnerScenario scenario,
SGPP::base::OCLConfigurationParameters currentParameters) {
//    std::cout << "-----------------------" << std::endl;
//    std::vector<std::string> keys = currentParameters.getKeys();
//    for (std::string key : keys) {
//        if (key.compare(0, 7, "KERNEL_", 0, 7) == 0) {
//            std::cout << "key: " << key << " value: " << currentParameters.get(key) << std::endl;
//        }
//    }
//    static uint64_t reachedCounter = 0;
//    reachedCounter += 1;
//    std::cout << "reached: " << reachedCounter << std::endl;

    bool verbose = false;
    SGPP::datadriven::MetaLearner learner(scenario.getGridConfig(), scenario.getSolverConfigurationRefine(),
            scenario.getSolverConfigurationFinal(), scenario.getAdaptivityConfiguration(), scenario.getLambda(),
            verbose);

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMP, currentParameters);

    std::string fileName = scenario.getDatasetFileName();

    double duration = std::numeric_limits<double>::max();
    try {
        learner.learn(configuration, fileName);

        LearnerTiming timing = learner.getLearnerTiming();

        duration = timing.timeComplete_;
    } catch (SGPP::base::operation_exception &exception) {
        std::cout << "invalid combination detected" << std::endl;
    }

    std::cout << "duration: " << duration << std::endl;
    return duration;
}

}
}

#endif
